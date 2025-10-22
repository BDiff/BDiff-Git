import os
import re
from collections import OrderedDict
from pprint import pprint
import copy
import numpy as np
from scipy.optimize import linear_sum_assignment
import myers


# Levenshtein_ratio
def levenshtein_ratio(s1, s2):
    if len(s1) == 0 and len(s2) == 0:
        return 1.0

    # 初始化动态规划表
    dp = np.zeros((len(s1) + 1, len(s2) + 1), dtype=float)
    for i in range(len(s1) + 1):
        dp[i, 0] = i * 1.0  # 删除代价
    for j in range(len(s2) + 1):
        dp[0, j] = j * 1.0  # 插入代价

    # 填充DP表
    for i in range(1, len(s1) + 1):
        for j in range(1, len(s2) + 1):
            cost = 0 if s1[i - 1] == s2[j - 1] else 2.0
            dp[i, j] = min(
                dp[i - 1, j] + 1,  # 删除
                dp[i, j - 1] + 1,  # 插入
                dp[i - 1, j - 1] + cost  # 替换
            )

    distance = dp[len(s1), len(s2)]
    total_length = len(s1) + len(s2)
    return (total_length - distance) / total_length


# 多匹配更新时的筛选时，同交点数随意删除可能使得筛掉的进入km算法被筛掉
# move和copy带更新时，计算权重需要加上更新的部分
def W_BESTI_LINE(src_line_no, dest_line_no, src_lines, dest_lines, ctx_length=4, line_sim_weight=0.6,
                 sim_threshold=0.5):
    # sim = line_sim * 0.6 + context_sim * 0.4
    # 先匹配相等
    src_upper_ctx, src_under_ctx = [], []
    if src_lines[src_line_no - 1].strip() == "" and dest_lines[dest_line_no - 1].strip() == "":
        # 如果两行一样，返回False
        if src_lines[src_line_no - 1] == dest_lines[dest_line_no - 1]:
            return False, 0
        line_sim = 1
    elif not src_lines[src_line_no - 1] or not dest_lines[dest_line_no - 1]:
        return False, 0
    else:
        line_sim = levenshtein_ratio(src_lines[src_line_no - 1].strip(), dest_lines[dest_line_no - 1].strip())
    if src_line_no <= ctx_length:
        src_upper_ctx = src_lines[:src_line_no - 1]
    else:
        src_upper_ctx = src_lines[src_line_no - ctx_length - 1: src_line_no - 1]
    if src_line_no + ctx_length > len(src_lines):
        src_under_ctx = src_lines[src_line_no:]
    else:
        src_under_ctx = src_lines[src_line_no: src_line_no + ctx_length]
    if dest_line_no <= ctx_length:
        dest_upper_ctx = dest_lines[:dest_line_no - 1]
    else:
        dest_upper_ctx = dest_lines[dest_line_no - ctx_length - 1: dest_line_no - 1]
    if dest_line_no + ctx_length >= len(dest_lines):
        dest_under_ctx = dest_lines[dest_line_no:]
    else:
        dest_under_ctx = dest_lines[dest_line_no: dest_line_no + ctx_length]
    upper_group = [group[0].strip() == group[1].strip() for group in zip(src_upper_ctx, dest_upper_ctx)]
    under_group = [group[0].strip() == group[1].strip() for group in zip(src_under_ctx, dest_under_ctx)]
    if len(upper_group) + len(under_group) == 0:
        return True if line_sim >= sim_threshold else False, round(line_sim, 3)
    ctx_sim = (upper_group.count(True) + under_group.count(True)) / (len(upper_group) + len(under_group))
    synthetic_sim = line_sim * line_sim_weight + ctx_sim * (1 - line_sim_weight)
    return True if synthetic_sim >= sim_threshold else False, round(synthetic_sim, 3)


def is_line_change(src_line_no, dest_line_no, src_lines, dest_lines, ctx_length=4, line_sim_threshold=0.6,
                   sim_threshold=0.5):
    # sim = line_sim * 0.6 + context_sim * 0.4
    # 先匹配相等，
    src_ctx, dest_ctx = "", ""
    src_index, src_length = src_line_no - 2, len(src_lines)
    i = 0
    while i < ctx_length and src_index >= 0:
        if src_lines[src_index].strip() != '':
            src_ctx = src_lines[src_index] + src_ctx
            i += 1
        src_index -= 1
    i = 0
    src_index = src_line_no
    while i < ctx_length and src_index < src_length:
        if src_lines[src_index].strip() != '':
            src_ctx = src_ctx + src_lines[src_index]
            i += 1
        src_index += 1
    j = 0
    dest_index, dest_length = dest_line_no - 2, len(dest_lines)
    while j < ctx_length and dest_index >= 0:
        if dest_lines[dest_index].strip() != '':
            dest_ctx = dest_lines[dest_index] + dest_ctx
            j += 1
        dest_index -= 1
    j = 0
    dest_index = dest_line_no
    while j < ctx_length and dest_index < dest_length:
        if dest_lines[dest_index].strip() != '':
            dest_ctx = dest_ctx + dest_lines[dest_index]
            j += 1
        dest_index += 1
    ctx_sim = levenshtein_ratio(src_ctx, dest_ctx)
    line_sim = levenshtein_ratio(src_lines[src_line_no - 1], dest_lines[dest_line_no - 1])
    sim = line_sim * line_sim_threshold + ctx_sim * (1 - line_sim_threshold)
    return True if sim >= sim_threshold else False, sim

def construct_line_data(diffs, indent_tabs_size):
    hunks, hunk_dels, hunk_adds = [], [], []
    diff_line_dict_src = OrderedDict()
    diff_line_dict_added = OrderedDict()
    a_line_no = 0
    b_line_no = 0
    diff_scripts = []
    hunk = 0
    counting_hunk = False
    for mode, line in diffs:
        if mode == 'k':
            if counting_hunk:
                counting_hunk = False
                hunks.append([hunk_dels, hunk_adds])
                hunk_dels, hunk_adds = [], []
            a_line_no += 1
            b_line_no += 1
            diff_line_dict_src[a_line_no] = (
                line.lstrip().rstrip('\n'), compute_line_indent(line, indent_tabs_size), 'k')
            diff_scripts.append('k' + str(a_line_no))
        elif mode == 'r':
            if not counting_hunk:
                hunk += 1
                counting_hunk = True
            a_line_no += 1
            hunk_dels.append(a_line_no)
            diff_line_dict_src[a_line_no] = (
                line.lstrip().rstrip('\n'), compute_line_indent(line, indent_tabs_size), 'r', hunk)
            diff_scripts.append('r' + str(a_line_no))
        elif mode == 'i':
            if not counting_hunk:
                hunk += 1
                counting_hunk = True
            b_line_no += 1
            hunk_adds.append(b_line_no)
            diff_line_dict_added[b_line_no] = (
                line.lstrip().rstrip('\n'), compute_line_indent(line, indent_tabs_size), 'i', hunk)
            diff_scripts.append('i' + str(b_line_no))
    return diff_line_dict_src, diff_line_dict_added, diff_scripts, hunks


def _find_same_left(a: str, b: str, /, min_len: int) -> int:
    """查找两个字符串左侧最大相同区域的字符数"""
    low, high = 0, min_len

    while low < high:
        mid = (low + high) >> 1
        if a[low:mid + 1] == b[low:mid + 1]:
            low = mid + 1
        else:
            high = mid

    return low


def _find_diff_area(a: str, b: str, /) -> tuple[int, int]:
    """查找两个字符串不相同的区域，返回左右相同字符数"""
    min_len = min(len(a), len(b))

    left = _find_same_left(a, b, min_len)
    right = _find_same_left(a[::-1], b[::-1], min_len)

    return left, min(right, min_len - left)


def find_diff_area(a: str, b: str, /) -> list[list[list[int]]]:
    """查找两个字符串不相同区域，左闭右开"""
    start, end = _find_diff_area(a, b)

    area_a = [start, len(a) - end - 1]
    area_b = [start, len(b) - end - 1]

    if area_a[0] > area_a[1]:
        area_a = []

    if area_b[0] > area_b[1]:
        area_b = []

    return [[area_a], [area_b]]


def construct_str_diff_data(src_tuple: tuple, dest_tuple: tuple):
    left_range, right_range = find_diff_area(src_tuple[0], dest_tuple[0])
    if not left_range[0] and not right_range[0]:
        left_range[0].append(0)
        left_range[0].append(src_tuple[1][1] + src_tuple[1][2] - 1)
        right_range[0].append(0)
        right_range[0].append(dest_tuple[1][1] + dest_tuple[1][2] - 1)
        return [left_range, right_range]
    if left_range[0]:
        left_range[0][0] = left_range[0][0] + src_tuple[1][1] + src_tuple[1][2]
        left_range[0][1] = left_range[0][1] + src_tuple[1][1] + src_tuple[1][2]
    if right_range[0]:
        right_range[0][0] = right_range[0][0] + dest_tuple[1][1] + dest_tuple[1][2]
        right_range[0][1] = right_range[0][1] + dest_tuple[1][1] + dest_tuple[1][2]
    return [left_range, right_range]


def construct_str_myers_diff_data(diffs, src_indent, dest_indent):
    a_line_no = -1
    b_line_no = -1
    src_indexes = []
    dest_indexes = []
    cur_mode = None
    cur_start = None
    for mode, line in diffs:
        if mode == 'k':
            if cur_mode == 'r':
                src_indexes.append((cur_start + src_indent, a_line_no + dest_indent))
            elif cur_mode == 'i':
                dest_indexes.append((cur_start + src_indent, b_line_no + dest_indent))
            cur_mode = 'k'
            a_line_no += 1
            b_line_no += 1
        elif mode == 'r':
            a_line_no += 1
            if cur_mode == 'k':
                cur_start = a_line_no
                cur_mode = 'r'
        elif mode == 'i':
            b_line_no += 1
            if cur_mode == 'k':
                cur_start = b_line_no
                cur_mode = 'i'
            elif cur_mode == 'r':
                src_indexes.append((cur_start + src_indent, a_line_no + src_indent))
                cur_start = b_line_no
                cur_mode = 'i'
    if cur_mode == 'r':
        src_indexes.append((cur_start + src_indent, a_line_no + src_indent))
    elif cur_mode == 'i':
        dest_indexes.append((cur_start + dest_indent, b_line_no + dest_indent))
    return src_indexes, dest_indexes


def compute_line_indent(diff_line, indent_tabs_size):
    # return: total, n_spaces, n_tabs
    if diff_line.startswith(" ") or diff_line.startswith("\t"):
        if diff_line.lstrip() == "":
            first_chara_index = len(diff_line) + 1
        else:
            first_chara_index = diff_line.find(diff_line.lstrip()[0])
        n_spaces = diff_line[:first_chara_index].count(" ")
        n_tabs = diff_line[:first_chara_index].count("\t")
        return n_spaces + n_tabs * indent_tabs_size, n_spaces, n_tabs
    else:
        return 0, 0, 0


def is_pure_punctuation(s):
    if not s:
        return True
    pattern = r'^[~`!@#$%^&*()-_+={}\[\]|\\:;"\'<,>.?/\n\s]+$'
    return bool(re.match(pattern, s))


def pure_block_len(block_length, src_start, src_lines_list, added_start, added_lines_list, pure_mv_block_contain_punc,
                   pure_cp_block_contain_punc, mode):
    i = 0
    pure_block_length = block_length
    while i < block_length:
        if not src_lines_list[src_start - 1] and not added_lines_list[added_start - 1]:
            pure_block_length -= 1
        elif ((not pure_mv_block_contain_punc and mode == 'r') or (
                not pure_cp_block_contain_punc and mode == 'k')) and is_pure_punctuation(
            src_lines_list[src_start - 1]) and is_pure_punctuation(added_lines_list[added_start - 1]):
            pure_block_length -= 1
        src_start += 1
        added_start += 1
        i += 1
    return pure_block_length


def mapping_block_move(src_lines, added_lines, src_all_lines, dest_all_lines, min_block_length, diff_scripts,
                       pure_mv_block_contain_punc, count_mv_block_update):
    mappings = []
    checked_dict = {}
    # 对added_lines分块
    for added_line in added_lines:
        if added_lines[added_line][0] == "":
            continue
        for src_line in src_lines:
            if src_lines[src_line][0] == "" or (src_line, added_line) in checked_dict or src_lines[src_line][2] == "k":
                continue
            else:
                checked_dict[(src_line, added_line)] = True
            cur_src_line_no = src_line
            cur_added_line_no = added_line
            indent_diff = added_lines[added_line][1][0] - src_lines[src_line][1][0]
            src_mode = "r"
            block_length = 0
            pure_block_length = 0
            edit_actions = 2
            m_updates = []
            while cur_src_line_no in src_lines and cur_added_line_no in added_lines and src_lines[cur_src_line_no][
                2] == src_mode and (src_lines[cur_src_line_no][0] == added_lines[cur_added_line_no][0] or (
                    src_lines[cur_src_line_no][0] != added_lines[cur_added_line_no][0] and count_mv_block_update and
                    levenshtein_ratio(src_lines[cur_src_line_no][0], added_lines[cur_added_line_no][0]) >= 0.6)) and (
                    (added_lines[cur_added_line_no][0] != "" and added_lines[cur_added_line_no][1][0] -
                     src_lines[cur_src_line_no][1][0] == indent_diff) or added_lines[cur_added_line_no][0] == ""):
                if count_mv_block_update and src_lines[cur_src_line_no][0] != added_lines[cur_added_line_no][0]:
                    edit_actions += 1
                    m_updates.append([cur_src_line_no, cur_added_line_no])
                # 非空行数大于等于2
                if src_lines[cur_src_line_no][0] != "" and added_lines[cur_added_line_no][0] != "":
                    if pure_mv_block_contain_punc or not (
                            is_pure_punctuation(src_lines[cur_src_line_no][0]) and is_pure_punctuation(
                        added_lines[cur_added_line_no][0])):
                        pure_block_length += 1
                # 添加已经判断过的行映射，后面无需再判断
                checked_dict[(cur_src_line_no, cur_added_line_no)] = True
                cur_src_line_no += 1
                cur_added_line_no += 1
                block_length += 1
            if pure_block_length >= min_block_length and not is_pure_punctuation(
                    "".join([src_lines[line][0] for line in range(src_line, src_line + block_length)])):
                # 回溯看前面是否有空格
                cur_src_line_no = src_line - 1
                cur_added_line_no = added_line - 1
                while cur_src_line_no >= 1 and cur_added_line_no >= 1:
                    if cur_src_line_no in src_lines and cur_added_line_no in added_lines and src_lines[cur_src_line_no][
                        2] == src_mode and src_lines[cur_src_line_no][0] == "" and added_lines[cur_added_line_no][0] == "":
                        src_line = cur_src_line_no
                        added_line = cur_added_line_no
                        block_length += 1
                        cur_src_line_no -= 1
                        cur_added_line_no -= 1
                    else:
                        break
                ctx_similarity = context_similarity(src_line, added_line, block_length, src_all_lines, dest_all_lines)
                # move_type: h:平移，d:下移，u:上移
                if src_lines[src_line][3] == added_lines[added_line][3]:
                    move_type = "h"
                elif src_lines[src_line][3] < added_lines[added_line][3]:
                    move_type = "d"
                elif src_lines[src_line][3] > added_lines[added_line][3]:
                    move_type = "u"
                # 非平移有缩进时edit_action要+1
                if move_type == 'h' and indent_diff == 0:
                    continue
                if indent_diff != 0 and move_type != "h":
                    edit_actions += 1
                rd = relative_distance(src_line, added_line, block_length, diff_scripts)
                candidate = {"mode": src_mode, "block_length": block_length,
                             "src_start": src_line, "added_start": added_line, "context_similarity": ctx_similarity,
                             "weight": (edit_actions / block_length + (1 - ctx_similarity) / 10 + rd / 100),
                             "move_type": move_type,
                             "updates": m_updates,
                             "indent_diff": indent_diff,
                             "edit_actions": edit_actions,
                             "relative_distance": rd}
                # 如果不是平移只有update，就不算move
                mappings.append(candidate)
    return mappings


def mapping_block_copy(src_lines, added_lines, src_all_lines, dest_all_lines, min_copy_block_length, hunks, diff_scripts,
                       pure_cp_block_contain_punc, count_cp_block_update):
    mappings = []
    # 对added_lines分块
    checked_dict = {}
    for added_line in added_lines:
        if added_lines[added_line][0] == "":
            continue
        candidates = {}  # key: src_line, value: mapping
        for src_line in src_lines:
            if src_lines[src_line][0] == "" or (src_line, added_line) in checked_dict:
                continue
            else:
                checked_dict[(src_line, added_line)] = True
            cur_src_line_no = src_line
            cur_added_line_no = added_line
            indent_diff = added_lines[added_line][1][0] - src_lines[src_line][1][0]
            src_mode = "k"
            block_length = 0
            pure_block_length = 0
            edit_actions = 4
            c_updates = []
            while cur_src_line_no in src_lines and cur_added_line_no in added_lines and \
                    (src_lines[cur_src_line_no][0] == added_lines[cur_added_line_no][0] or (
                            src_lines[cur_src_line_no][0] != added_lines[cur_added_line_no][
                        0] and count_cp_block_update and
                            levenshtein_ratio(src_lines[cur_src_line_no][0],
                                              added_lines[cur_added_line_no][0]) >= 0.6)) and (
                    (added_lines[cur_added_line_no][0] != "" and added_lines[cur_added_line_no][
                        1][0] - src_lines[cur_src_line_no][1][0] == indent_diff) or added_lines[cur_added_line_no][
                        0] == ""):
                # 添加已经判断过的行映射，后面无需再判断
                checked_dict[(cur_src_line_no, cur_added_line_no)] = True
                if count_cp_block_update and src_lines[cur_src_line_no][0] != added_lines[cur_added_line_no][0]:
                    edit_actions += 1
                    c_updates.append([cur_src_line_no, cur_added_line_no])
                if src_lines[cur_src_line_no][0] != "" and added_lines[cur_added_line_no][0] != "":
                    if pure_cp_block_contain_punc or not (
                            is_pure_punctuation(src_lines[cur_src_line_no][0]) and is_pure_punctuation(
                        added_lines[cur_added_line_no][0])):
                        pure_block_length += 1
                cur_src_line_no += 1
                cur_added_line_no += 1
                block_length += 1
            if pure_block_length >= min_copy_block_length and not copy_block_in_hunk(
                    {"mode": src_mode, "block_length": block_length, "src_start": src_line, "added_start": added_line},
                    hunks) and not is_pure_punctuation(
                "".join([src_lines[line][0] for line in range(src_line, src_line + block_length)])):
                # 回溯看前面是否有空格
                cur_src_line_no = src_line - 1
                cur_added_line_no = added_line - 1
                while cur_src_line_no >= 1 and cur_added_line_no >= 1:
                    if cur_src_line_no in src_lines and cur_added_line_no in added_lines and src_lines[cur_src_line_no][
                        0] == "" and added_lines[cur_added_line_no][0] == "":
                        src_line = cur_src_line_no
                        added_line = cur_added_line_no
                        block_length += 1
                        cur_src_line_no -= 1
                        cur_added_line_no -= 1
                    else:
                        break
                if indent_diff != 0:
                    edit_actions += 1
                ctx_similarity = context_similarity(src_line, added_line, block_length, src_all_lines, dest_all_lines)
                rd = relative_distance(src_line, added_line, block_length, diff_scripts)
                weight = edit_actions / block_length + (1 - ctx_similarity) / 10 + rd / 100
                # 判断candidates中是否存在同样块大小且的mapping
                for s in candidates:
                    if candidates[s]['block_length'] == block_length:
                        if candidates[s]['weight'] > weight:
                            # 删除这个键值对，新增另一个
                            del candidates[s]
                            candidate = {"mode": src_mode, "block_length": block_length,
                                         "src_start": src_line, "added_start": added_line,
                                         "context_similarity": ctx_similarity,
                                         "weight": weight,
                                         "updates": c_updates,
                                         "indent_diff": indent_diff,
                                         "edit_actions": edit_actions,
                                         "relative_distance": rd}
                            candidates[src_line] = candidate
                            break
                        else:
                            break
                else:
                    candidate = {"mode": src_mode, "block_length": block_length,
                                 "src_start": src_line, "added_start": added_line, "context_similarity": ctx_similarity,
                                 "weight": weight,
                                 "updates": c_updates,
                                 "indent_diff": indent_diff,
                                 "edit_actions": edit_actions,
                                 "relative_distance": rd}
                    candidates[src_line] = candidate
        for s in candidates:
            mappings.append(candidates[s])
    return mappings


def context_similarity(src_start, dest_start, block, src_lines, dest_lines):
    # 参考LHDiff
    src_context = construct_context(src_start, block, src_lines)
    dest_context = construct_context(dest_start, block, dest_lines)
    return levenshtein_ratio(src_context, dest_context)


def construct_context(start, block_length, lines):
    context = ""
    # 前溯
    i, j = 1, 1
    start_ptr = start - 2
    while i < 5 and start_ptr >= 0:
        if lines[start_ptr].strip() == "":
            start_ptr -= 1
            continue
        else:
            context = lines[start_ptr].strip() + " " + context
            start_ptr -= 1
            i += 1
    # 后溯
    start_ptr = start + block_length - 1
    while j < 5 and start_ptr < len(lines):
        if lines[start_ptr].strip() == "":
            start_ptr += 1
            continue
        else:
            context = context + " " + lines[start_ptr].strip()
            start_ptr += 1
            j += 1
    return context


def judge_overlap_type(assigned_start, assigned_block_length, overlapped_start, overlapped_block_length):
    """
    equal: e
    cover: c
    inner: i
    up: u
    down: d
    """
    if assigned_start == overlapped_start and assigned_block_length == overlapped_block_length:
        return "e"
    elif assigned_start >= overlapped_start and (assigned_start + assigned_block_length) <= (
            overlapped_start + overlapped_block_length):
        return "c"
    elif assigned_start <= overlapped_start and (assigned_start + assigned_block_length) >= (
            overlapped_start + overlapped_block_length):
        return "i"
    # elif overlapped_start <= assigned_start < (overlapped_start + overlapped_block_length - 1):
    elif overlapped_start <= assigned_start and overlapped_start + overlapped_block_length - 1 >= assigned_start:
        return "u"
    elif assigned_start <= overlapped_start and assigned_start + assigned_block_length - 1 >= overlapped_start:
        return "d"
    else:
        return None


def km_compute(mappings, src_all_lines, dest_all_lines, min_move_block_length=2, min_copy_block_length=2,
               pure_mv_block_contain_punc=True, pure_cp_block_contain_punc=True):
    mappings = [x for i, x in enumerate(mappings) if x not in mappings[:i]]
    grouped_mappings_src = []
    grouped_mappings_added = []
    mappings.sort(key=lambda x: x["src_start"])
    km_start = 0
    km_end = 0
    for index, mapping in enumerate(mappings):
        # state: a: assigned, d: deleted, None: waiting for assigned, s: sliced
        mapping['state'] = None
        found_src_mapping = False
        for group_src in grouped_mappings_src:
            for mapping_src in group_src:
                if not (
                        (mapping['src_start'] + mapping['block_length'] - 1) < mapping_src['src_start']
                        or mapping['src_start'] > (mapping_src['src_start'] + mapping_src['block_length'] - 1)) and \
                        mapping['mode'] != 'k' and mapping_src['mode'] != 'k':
                    mapping['km_start'] = mapping_src['km_start']
                    group_src.append(mapping)
                    found_src_mapping = True
                    break
            if found_src_mapping:
                break
        if not found_src_mapping:
            mapping['km_start'] = km_start
            grouped_mappings_src.append([mapping])
            km_start += 1
    mappings.sort(key=lambda x: x["added_start"])
    for index, mapping in enumerate(mappings):
        found_added_mapping = False
        for group_added in grouped_mappings_added:
            for mapping_added in group_added:
                if not (
                        (mapping['added_start'] + mapping['block_length'] - 1) < mapping_added['added_start']
                        or mapping['added_start'] > (mapping_added['added_start'] + mapping_added['block_length'] - 1)):
                    mapping['km_end'] = mapping_added['km_end']
                    group_added.append(mapping)
                    found_added_mapping = True
                    break
            if found_added_mapping:
                break
        if not found_added_mapping:
            mapping['km_end'] = km_end
            km_end += 1
            grouped_mappings_added.append([mapping])
    # 计算km的mapping、以及去除km-mapping后的mapping
    cost_matrix = np.full((km_start, km_end), 1000.0)
    for mapping in mappings:
        x = mapping['km_start']
        y = mapping['km_end']
        if cost_matrix[x][y] != 1000:
            if mapping['weight'] < cost_matrix[x][y]:
                cost_matrix[x][y] = mapping['weight']
        else:
            cost_matrix[x][y] = mapping['weight']
    # 计算最优分配（会忽略掉所有成本为无穷大的分配）
    row_ind, col_ind = linear_sum_assignment(cost_matrix)
    assignments = list(zip(row_ind, col_ind))
    # 查找剩余没有匹配上的
    km_matches = []
    remain_mappings = []
    for assginment in assignments:
        present_assignment = {}
        max_weight = len(src_all_lines) * 2
        for mapping1 in mappings:
            # 已是映射
            if not mapping1['state'] and mapping1['km_start'] == assginment[0] and mapping1['km_end'] == assginment[
                1] and mapping1['weight'] < max_weight:
                present_assignment = mapping1
                max_weight = mapping1['weight']
        if present_assignment:
            present_assignment['state'] = 'a'
            km_matches.append(present_assignment)
        # 没有映射的情况
        else:
            continue
        for mapping2 in mappings:
            # 已是映射
            if mapping2['state'] or (
                    mapping2['km_start'] == assginment[0] and mapping2['km_end'] == assginment[1] and mapping2[
                'mode'] != 'u'):
                continue
            # 左节点重叠
            elif mapping2['km_start'] == assginment[0]:
                # cover出问题的需要去掉
                overlap_type = judge_overlap_type(present_assignment['src_start'], present_assignment['block_length'],
                                                  mapping2['src_start'], mapping2['block_length'])
                if overlap_type == 'e' or overlap_type == 'i':
                    mapping2['state'] = 'd'
                    continue
                elif overlap_type == None:
                    mapping2['state'] = 's'
                    remain_mappings.append(mapping2)
                elif overlap_type == 'c':
                    # 两段, 先过完src，再过一遍dest
                    # 先判断updates是否在片段中
                    # 优化：先判断up_offset是否大于pure_block_length
                    cannot_be_sliced = True
                    up_offset = present_assignment['src_start'] - mapping2['src_start']
                    pure_up_offset = pure_block_len(up_offset, mapping2['src_start'], src_all_lines,
                                                    mapping2['added_start'], dest_all_lines, pure_mv_block_contain_punc,
                                                    pure_cp_block_contain_punc, mapping2['mode'])
                    edit_actions = 1 if mapping2['mode'] == 'u' else 2 if mapping2['mode'] == 'r' else 4
                    if (mapping2['indent_diff'] != 0 and mapping2['mode'] == 'k') or (
                            mapping2['indent_diff'] != 0 and mapping2['mode'] == 'r' and mapping2['move_type'] != 'h'):
                        edit_actions += 1
                    updates = []
                    for ud in mapping2['updates']:
                        if ud[0] in range(mapping2['src_start'], present_assignment['src_start']):
                            updates.append(ud)
                            edit_actions += 1
                    if (mapping2['mode'] == 'r' and pure_up_offset >= min_move_block_length) or (
                            mapping2['mode'] == 'k' and pure_up_offset >= min_copy_block_length):
                        ctx_similarity = context_similarity(mapping2['src_start'], mapping2['added_start'], up_offset,
                                                            src_all_lines,
                                                            dest_all_lines)
                        mapping2['state'] = 's'
                        cannot_be_sliced = False
                        remain_mappings.append({'mode': mapping2['mode'], 'block_length': up_offset,
                                                'context_similarity': ctx_similarity,
                                                'weight': edit_actions / up_offset + (1 - ctx_similarity) / 10 +
                                                          mapping2['relative_distance'] / 100,
                                                'src_start': mapping2['src_start'],
                                                'added_start': mapping2['added_start'],
                                                'km_start': mapping2['km_start'],
                                                'km_end': mapping2['km_end'],
                                                'indent_diff': mapping2['indent_diff'],
                                                'edit_actions': edit_actions,
                                                'updates': updates,
                                                'relative_distance': mapping2['relative_distance'],
                                                'state': None
                                                })
                        if mapping2['mode'] == 'r':
                            remain_mappings[-1]['move_type'] = mapping2['move_type']
                    down_offset = mapping2['src_start'] + mapping2['block_length'] - (
                            present_assignment['src_start'] + present_assignment['block_length'])
                    pure_down_offset = pure_block_len(down_offset, present_assignment['src_start'] + present_assignment[
                        'block_length'], src_all_lines, mapping2['added_start'] + (
                                                              present_assignment['src_start'] + present_assignment[
                                                          'block_length'] - mapping2['src_start']), dest_all_lines,
                                                      pure_mv_block_contain_punc, pure_cp_block_contain_punc,
                                                      mapping2['mode'])
                    edit_actions = 1 if mapping2['mode'] == 'u' else 2 if mapping2['mode'] == 'r' else 4
                    # 有缩进
                    if (mapping2['indent_diff'] != 0 and mapping2['mode'] == 'k') or (
                            mapping2['indent_diff'] != 0 and mapping2['mode'] == 'r' and mapping2['move_type'] != 'h'):
                        edit_actions += 1
                    updates = []
                    for ud in mapping2['updates']:
                        if ud[0] in range(present_assignment['src_start'] + present_assignment['block_length'],
                                          present_assignment['src_start'] + present_assignment[
                                              'block_length'] + down_offset):
                            # 有更新
                            updates.append(ud)
                            edit_actions += 1
                    if (mapping2['mode'] == 'r' and pure_down_offset >= min_move_block_length) or (
                            mapping2['mode'] == 'k' and pure_down_offset >= min_copy_block_length):
                        cannot_be_sliced = False
                        ctx_similarity = context_similarity(present_assignment['src_start'] + present_assignment[
                            'block_length'], mapping2['added_start'] + present_assignment['src_start'] +
                                                            present_assignment[
                                                                'block_length'] - mapping2['src_start'], down_offset,
                                                            src_all_lines, dest_all_lines)
                        mapping2['state'] = 's'
                        remain_mappings.append({'mode': mapping2['mode'], 'block_length': down_offset,
                                                'context_similarity': ctx_similarity,
                                                'weight': edit_actions / down_offset + (1 - ctx_similarity) / 10 +
                                                          mapping2['relative_distance'] / 100,
                                                'src_start': present_assignment['src_start'] + present_assignment[
                                                    'block_length'],
                                                'added_start': mapping2['added_start'] + (
                                                        present_assignment['src_start'] + present_assignment[
                                                    'block_length'] - mapping2['src_start']),
                                                'km_start': mapping2['km_start'],
                                                'km_end': mapping2['km_end'],
                                                'indent_diff': mapping2['indent_diff'],
                                                'edit_actions': edit_actions,
                                                'updates': updates,
                                                "relative_distance": mapping2['relative_distance'],
                                                'state': None
                                                })
                        if mapping2['mode'] == 'r':
                            remain_mappings[-1]['move_type'] = mapping2['move_type']
                    if cannot_be_sliced:
                        mapping2['state'] = 's'
                elif overlap_type == 'u':
                    up_offset = present_assignment['src_start'] - mapping2['src_start']
                    pure_up_offset = pure_block_len(up_offset, mapping2['src_start'], src_all_lines,
                                                    mapping2['added_start'], dest_all_lines, pure_mv_block_contain_punc,
                                                    pure_cp_block_contain_punc, mapping2['mode'])
                    edit_actions = 1 if mapping2['mode'] == 'u' else 2 if mapping2['mode'] == 'r' else 4
                    if (mapping2['indent_diff'] != 0 and mapping2['mode'] == 'k') or (
                            mapping2['indent_diff'] != 0 and mapping2['mode'] == 'r' and mapping2['move_type'] != 'h'):
                        edit_actions += 1
                    updates = []
                    for ud in mapping2['updates']:
                        if ud[0] in range(mapping2['src_start'], present_assignment['src_start']):
                            updates.append(ud)
                            edit_actions += 1
                    ctx_similarity = context_similarity(mapping2['src_start'], mapping2['added_start'], up_offset,
                                                        src_all_lines,
                                                        dest_all_lines)
                    mapping2['state'] = 's'
                    if (mapping2['mode'] == 'r' and pure_up_offset >= min_move_block_length) or (
                            mapping2['mode'] == 'k' and pure_up_offset >= min_copy_block_length):
                        remain_mappings.append({'mode': mapping2['mode'], 'block_length': up_offset,
                                                'context_similarity': ctx_similarity,
                                                'weight': edit_actions / up_offset + (1 - ctx_similarity) / 10 +
                                                          mapping2['relative_distance'] / 100,
                                                'src_start': mapping2['src_start'],
                                                'added_start': mapping2['added_start'],
                                                'km_start': mapping2['km_start'],
                                                'km_end': mapping2['km_end'],
                                                'indent_diff': mapping2['indent_diff'],
                                                'edit_actions': edit_actions,
                                                'updates': updates,
                                                "relative_distance": mapping2['relative_distance'],
                                                'state': None
                                                })
                        if mapping2['mode'] == 'r':
                            remain_mappings[-1]['move_type'] = mapping2['move_type']
                    else:
                        mapping2['state'] = 's'
                elif overlap_type == 'd':
                    down_offset = mapping2['src_start'] + mapping2['block_length'] - (
                            present_assignment['src_start'] + present_assignment['block_length'])
                    pure_down_offset = pure_block_len(down_offset, present_assignment['src_start'] + present_assignment[
                        'block_length'], src_all_lines, mapping2['added_start'] + (
                                                              present_assignment['src_start'] + present_assignment[
                                                          'block_length'] - mapping2['src_start']), dest_all_lines,
                                                      pure_mv_block_contain_punc, pure_cp_block_contain_punc,
                                                      mapping2['mode'])
                    edit_actions = 1 if mapping2['mode'] == 'u' else 2 if mapping2['mode'] == 'r' else 4
                    # 有缩进
                    if (mapping2['indent_diff'] != 0 and mapping2['mode'] == 'k') or (
                            mapping2['indent_diff'] != 0 and mapping2['mode'] == 'r' and mapping2['move_type'] != 'h'):
                        edit_actions += 1
                    updates = []
                    for ud in mapping2['updates']:
                        if ud[0] in range(present_assignment['src_start'] + present_assignment['block_length'],
                                          present_assignment['src_start'] + present_assignment[
                                              'block_length'] + down_offset):
                            # 有更新
                            updates.append(ud)
                            edit_actions += 1
                    if (mapping2['mode'] == 'r' and pure_down_offset >= min_move_block_length) or (
                            mapping2['mode'] == 'k' and pure_down_offset >= min_copy_block_length):
                        ctx_similarity = context_similarity(present_assignment['src_start'] + present_assignment[
                            'block_length'], mapping2['added_start'] + present_assignment['src_start'] +
                                                            present_assignment['block_length'] - mapping2['src_start'],
                                                            down_offset,
                                                            src_all_lines, dest_all_lines)
                        mapping2['state'] = 's'
                        remain_mappings.append({'mode': mapping2['mode'], 'block_length': down_offset,
                                                'context_similarity': ctx_similarity,
                                                'weight': edit_actions / down_offset + (1 - ctx_similarity) / 10 +
                                                          mapping2['relative_distance'] / 100,
                                                'src_start': present_assignment['src_start'] + present_assignment[
                                                    'block_length'],
                                                'added_start': mapping2['added_start'] + (
                                                        present_assignment['src_start'] +
                                                        present_assignment[
                                                            'block_length'] - mapping2[
                                                            'src_start']),
                                                'km_start': mapping2['km_start'],
                                                'km_end': mapping2['km_end'],
                                                'indent_diff': mapping2['indent_diff'],
                                                'edit_actions': edit_actions,
                                                'updates': updates,
                                                "relative_distance": mapping2['relative_distance'],
                                                'state': None
                                                })
                        if mapping2['mode'] == 'r':
                            remain_mappings[-1]['move_type'] = mapping2['move_type']
                    else:
                        mapping2['state'] = 's'
    for assginment2 in assignments:
        for mapping2 in mappings:
            # start不相交、end相交的情形，交给下一步切割
            if not mapping2['state'] and mapping2['km_end'] == assginment2[1] and (
                    (mapping2['mode'] == 'k' and mapping2['block_length'] >= min_copy_block_length) or mapping2[
                'mode'] == 'u' or mapping2['mode'] == 'r'):
                # 判断remain_mappings是否已有mapping2
                if mapping2 not in remain_mappings:
                    remain_mappings.append(mapping2)
    # 再对end扫描分析
    final_remain_mappings = []
    for remain_mapping in remain_mappings:
        for km_match in km_matches:
            if remain_mapping['state'] != 'd' and remain_mapping['km_end'] == km_match['km_end'] and (
                    remain_mapping['km_start'] != km_match['km_start'] or (
                    km_match['mode'] == 'u' and remain_mapping['mode'] == 'u')):
                # 剥离end
                end_overlap_type = judge_overlap_type(km_match['added_start'], km_match['block_length'],
                                                      remain_mapping['added_start'], remain_mapping['block_length'])
                if end_overlap_type == 'e' or end_overlap_type == 'i':
                    continue
                # 无交集时，直接添加到final中
                elif end_overlap_type is None:
                    final_remain_mappings.append(remain_mapping)
                elif end_overlap_type == 'c':
                    up_offset = km_match['added_start'] - remain_mapping['added_start']
                    pure_up_offset = pure_block_len(up_offset, remain_mapping['src_start'], src_all_lines,
                                                    remain_mapping['added_start'], dest_all_lines,
                                                    pure_mv_block_contain_punc, pure_cp_block_contain_punc,
                                                    remain_mapping['mode'])
                    edit_actions = 1 if remain_mapping['mode'] == 'u' else 2 if remain_mapping['mode'] == 'r' else 3
                    if (remain_mapping['indent_diff'] != 0 and remain_mapping['mode'] == 'k') or (
                            remain_mapping['indent_diff'] != 0 and remain_mapping['mode'] == 'r' and remain_mapping[
                        'move_type'] != 'h'):
                        edit_actions += 1
                    updates = []
                    for ud in remain_mapping['updates']:
                        if ud[1] in range(remain_mapping['added_start'], km_match['added_start']):
                            updates.append(ud)
                            edit_actions += 1
                    if (remain_mapping['mode'] == 'r' and pure_up_offset >= min_move_block_length) or (
                            remain_mapping['mode'] == 'k' and pure_up_offset >= min_copy_block_length):
                        ctx_similarity = context_similarity(remain_mapping['src_start'], remain_mapping['added_start'],
                                                            up_offset, src_all_lines,
                                                            dest_all_lines)
                        final_remain_mappings.append({'mode': remain_mapping['mode'], 'block_length': up_offset,
                                                      'context_similarity': ctx_similarity,
                                                      'weight': edit_actions / up_offset + (1 - ctx_similarity) / 10 +
                                                                remain_mapping['relative_distance'] / 100,
                                                      'src_start': remain_mapping['src_start'],
                                                      'added_start': remain_mapping['added_start'],
                                                      'km_start': remain_mapping['km_start'],
                                                      'km_end': remain_mapping['km_end'],
                                                      'indent_diff': remain_mapping['indent_diff'],
                                                      'edit_actions': edit_actions,
                                                      'updates': updates,
                                                      'relative_distance': remain_mapping['relative_distance'],
                                                      'state': None
                                                      })
                        if remain_mapping['mode'] == 'r':
                            final_remain_mappings[-1]['move_type'] = remain_mapping['move_type']
                    down_offset = remain_mapping['added_start'] + remain_mapping['block_length'] - (
                            km_match['added_start'] + km_match['block_length'])
                    pure_down_offset = pure_block_len(down_offset, remain_mapping['src_start'] + km_match[
                        'added_start'] + km_match['block_length'] - remain_mapping['added_start'], src_all_lines,
                                                      km_match['added_start'] + km_match['block_length'],
                                                      dest_all_lines, pure_mv_block_contain_punc,
                                                      pure_cp_block_contain_punc, remain_mapping['mode'])
                    edit_actions = 1 if remain_mapping['mode'] == 'u' else 2 if remain_mapping['mode'] == 'r' else 3
                    # 有缩进
                    if (remain_mapping['indent_diff'] != 0 and remain_mapping['mode'] == 'k') or (
                            remain_mapping['indent_diff'] != 0 and remain_mapping['mode'] == 'r' and remain_mapping[
                        'move_type'] != 'h'):
                        edit_actions += 1
                    updates = []
                    for ud in remain_mapping['updates']:
                        if ud[1] in range(km_match['added_start'] + km_match['block_length'],
                                          km_match['added_start'] + km_match[
                                              'block_length'] + down_offset):
                            # 有更新
                            updates.append(ud)
                            edit_actions += 1
                    if (remain_mapping['mode'] == 'r' and pure_down_offset >= min_move_block_length) or (
                            remain_mapping['mode'] == 'k' and pure_down_offset >= min_copy_block_length):
                        ctx_similarity = context_similarity(
                            remain_mapping['src_start'] + km_match['added_start'] + km_match['block_length'] -
                            remain_mapping['added_start'], km_match['added_start'] + km_match[
                                'block_length'], down_offset, src_all_lines, dest_all_lines)
                        final_remain_mappings.append({'mode': remain_mapping['mode'], 'block_length': down_offset,
                                                      'context_similarity': ctx_similarity,
                                                      'weight': edit_actions / down_offset + (1 - ctx_similarity) / 10 +
                                                                remain_mapping['relative_distance'] / 100,
                                                      'src_start': remain_mapping['src_start'] + km_match[
                                                          'added_start'] + km_match['block_length'] - remain_mapping[
                                                                       'added_start'],
                                                      'added_start': km_match['added_start'] + km_match['block_length'],
                                                      'km_start': remain_mapping['km_start'],
                                                      'km_end': remain_mapping['km_end'],
                                                      'indent_diff': remain_mapping['indent_diff'],
                                                      'edit_actions': edit_actions,
                                                      'updates': updates,
                                                      "relative_distance": remain_mapping['relative_distance'],
                                                      'state': None
                                                      })
                        if remain_mapping['mode'] == 'r':
                            final_remain_mappings[-1]['move_type'] = remain_mapping['move_type']
                elif end_overlap_type == 'u':
                    up_offset = km_match['added_start'] - remain_mapping['added_start']
                    pure_up_offset = pure_block_len(up_offset, remain_mapping['src_start'], src_all_lines,
                                                    remain_mapping['added_start'], dest_all_lines,
                                                    pure_mv_block_contain_punc, pure_cp_block_contain_punc,
                                                    remain_mapping['mode'])
                    edit_actions = 1 if remain_mapping['mode'] == 'u' else 2 if remain_mapping['mode'] == 'r' else 3
                    if (remain_mapping['indent_diff'] != 0 and remain_mapping['mode'] == 'k') or (
                            remain_mapping['indent_diff'] != 0 and remain_mapping['mode'] == 'r' and remain_mapping[
                        'move_type'] != 'h'):
                        edit_actions += 1
                    updates = []
                    for ud in remain_mapping['updates']:
                        if ud[1] in range(remain_mapping['added_start'], km_match['added_start']):
                            updates.append(ud)
                            edit_actions += 1
                    if (remain_mapping['mode'] == 'r' and pure_up_offset >= min_move_block_length) or (
                            remain_mapping['mode'] == 'k' and pure_up_offset >= min_copy_block_length):
                        ctx_similarity = context_similarity(remain_mapping['src_start'], remain_mapping['added_start'],
                                                            up_offset, src_all_lines, dest_all_lines)
                        final_remain_mappings.append({'mode': remain_mapping['mode'], 'block_length': up_offset,
                                                      'context_similarity': ctx_similarity,
                                                      'weight': edit_actions / up_offset + (1 - ctx_similarity) / 10 +
                                                                remain_mapping['relative_distance'] / 100,
                                                      'src_start': remain_mapping['src_start'],
                                                      'added_start': remain_mapping['added_start'],
                                                      'km_start': remain_mapping['km_start'],
                                                      'km_end': remain_mapping['km_end'],
                                                      'indent_diff': remain_mapping['indent_diff'],
                                                      'edit_actions': edit_actions,
                                                      'updates': updates,
                                                      "relative_distance": remain_mapping['relative_distance'],
                                                      'state': None
                                                      })
                        if remain_mapping['mode'] == 'r':
                            final_remain_mappings[-1]['move_type'] = remain_mapping['move_type']
                elif end_overlap_type == 'd':
                    down_offset = remain_mapping['added_start'] + remain_mapping['block_length'] - (
                            km_match['added_start'] + km_match['block_length'])
                    pure_down_offset = pure_block_len(down_offset, remain_mapping['src_start'] + (
                            km_match['added_start'] + km_match['block_length'] -
                            remain_mapping['added_start']), src_all_lines,
                                                      km_match['added_start'] + km_match['block_length'],
                                                      dest_all_lines, pure_mv_block_contain_punc,
                                                      pure_cp_block_contain_punc, remain_mapping['mode'])
                    edit_actions = 1 if remain_mapping['mode'] == 'u' else 2 if remain_mapping['mode'] == 'r' else 3
                    # 有缩进
                    if (remain_mapping['indent_diff'] != 0 and remain_mapping['mode'] == 'k') or (
                            remain_mapping['indent_diff'] != 0 and remain_mapping['mode'] == 'r' and remain_mapping[
                        'move_type'] != 'h'):
                        edit_actions += 1
                    updates = []
                    for ud in remain_mapping['updates']:
                        if ud[1] in range(km_match['added_start'] + km_match['block_length'],
                                          km_match['added_start'] + km_match[
                                              'block_length'] + down_offset):
                            # 有更新
                            updates.append(ud)
                            edit_actions += 1
                    if (remain_mapping['mode'] == 'r' and pure_down_offset >= min_move_block_length) or (
                            remain_mapping['mode'] == 'k' and pure_down_offset >= min_copy_block_length):
                        ctx_similarity = context_similarity(
                            remain_mapping['src_start'] + km_match['added_start'] + km_match[
                                'block_length'] - remain_mapping['added_start'], km_match['added_start'] + km_match[
                                'block_length'], down_offset, src_all_lines, dest_all_lines)
                        final_remain_mappings.append({'mode': remain_mapping['mode'], 'block_length': down_offset,
                                                      'context_similarity': ctx_similarity,
                                                      'weight': edit_actions / down_offset + (1 - ctx_similarity) / 10 +
                                                                remain_mapping['relative_distance'] / 100,
                                                      'src_start': remain_mapping['src_start'] + (
                                                              km_match['added_start'] + km_match['block_length'] -
                                                              remain_mapping['added_start']),
                                                      'added_start': km_match['added_start'] + km_match['block_length'],
                                                      'km_start': remain_mapping['km_start'],
                                                      'km_end': remain_mapping['km_end'],
                                                      'indent_diff': remain_mapping['indent_diff'],
                                                      'edit_actions': edit_actions,
                                                      'updates': updates,
                                                      "relative_distance": remain_mapping['relative_distance'],
                                                      'state': None
                                                      })
                        if remain_mapping['mode'] == 'r':
                            final_remain_mappings[-1]['move_type'] = remain_mapping['move_type']
    return km_matches, final_remain_mappings


def generate_edit_action(mode, *args):
    if mode == 'move':
        if args[3] < 0:
            move_direction = " with moving left " + str(abs(args[3])) + " whitespaces."
        elif args[3] == 0:
            move_direction = ""
        else:
            move_direction = " with moving right " + str(args[3]) + " whitespaces."
        if args[0] == 1:
            return "Move 1 line from line " + str(args[1]) + " to line " + str(args[2]) + move_direction
        else:
            return "Move a " + str(args[0]) + "-line block from line " + str(args[1]) + " to line " + str(
                args[2]) + move_direction
    elif mode == 'copy':
        if args[3] < 0:
            move_direction = " with moving left " + str(abs(args[3])) + " whitespaces."
        elif args[3] == 0:
            move_direction = ""
        else:
            move_direction = " with moving right " + str(args[3]) + " whitespaces."
        return "Copy a " + str(args[0]) + "-line block from line " + str(args[1]) + " to line " + str(
            args[2]) + move_direction
    elif mode == 'm_update':
        return "Update line " + str(args[0]) + " to line " + str(args[1])
    elif mode == 'c_update':
        return "Update line " + str(args[0]) + " to line " + str(args[1])
    elif mode == 'update':
        if args[2] < 0:
            move_direction = " with moving left " + str(abs(args[2])) + " whitespaces."
        elif args[2] == 0:
            move_direction = ""
        else:
            move_direction = " with moving right " + str(args[2]) + " whitespaces."
        return "Update line " + str(args[0]) + " to line " + str(args[1]) + move_direction
    elif mode == 'insert':
        return "Insert line " + str(args[0])
    elif mode == 'delete':
        return "Delete line " + str(args[0])
    elif mode == "split":
        return "Split line " + str(args[0]) + " to lines " + str(args[1][0]) + "-" + str(args[1][-1])
    elif mode == "merge":
        return "Merge lines " + str(args[0][0]) + "-" + str(args[0][-1]) + " to line " + str(args[1])


def generate_edit_scripts_from_match(km_matches, diff_scripts, src_lines, added_lines, splits_merges, hunks, src_len,
                                     dest_len):
    diff_scripts_dict = {}
    k_pairs = {}
    src_line_no, dest_line_no = 1, 1
    for ds in diff_scripts:
        diff_scripts_dict[ds] = None
        if ds[0] == 'k':
            k_pairs[int(ds[1:])] = dest_line_no
            src_line_no += 1
            dest_line_no += 1
        elif ds[0] == 'r':
            src_line_no += 1
        else:
            dest_line_no += 1
    edit_scripts = []
    for split_merge in splits_merges:
        if len(split_merge[0]) == 1:
            edit_action = generate_edit_action("split", split_merge[0][0], split_merge[1])
            edit_scripts.append({"src_line": split_merge[0][0], "block_length": len(split_merge[1]),
                                 "dest_line": split_merge[1][0], "mode": "split",
                                 "edit_action": edit_action})
            diff_scripts_dict['r' + str(split_merge[0][0])] = "split-" + str(split_merge[1][0])
            for d_no in range(split_merge[1][0], split_merge[1][0] + len(split_merge[1])):
                diff_scripts_dict['i' + str(d_no)] = "split-" + str(split_merge[0][0])
        else:
            edit_action = generate_edit_action("merge", split_merge[0], split_merge[1][0])
            edit_scripts.append({"src_line": split_merge[0][0], "block_length": len(split_merge[0]),
                                 "dest_line": split_merge[1][0], "mode": "merge",
                                 "edit_action": edit_action})
            diff_scripts_dict['i' + str(split_merge[1][0])] = "merge-" + str(split_merge[0][0])
            for s_no in range(split_merge[0][0], split_merge[0][0] + len(split_merge[0])):
                diff_scripts_dict['r' + str(s_no)] = "merge-" + str(split_merge[1][0])
    for km_match in km_matches:
        if km_match['mode'] == 'k':
            edit_action = generate_edit_action("copy", km_match['block_length'], km_match['src_start'],
                                               km_match['added_start'], added_lines[km_match['added_start']][1][0] -
                                               src_lines[km_match['src_start']][1][0])
            edit_scripts.append({"src_line": km_match['src_start'], "block_length": km_match['block_length'],
                                 "dest_line": km_match['added_start'], "mode": "copy",
                                 "indent_offset": km_match['indent_diff'],
                                 "edit_action": edit_action,
                                 "updates": km_match['updates']})
            for d_no in range(km_match['added_start'], km_match['added_start'] + km_match['block_length']):
                diff_scripts_dict['i' + str(d_no)] = "copy-" + str(km_match['src_start'])
            # 块拷贝更新
            for update in km_match['updates']:
                copy_update = {}
                edit_action = generate_edit_action("c_update", update[0], update[1])
                copy_update["src_line"] = update[0]
                copy_update["dest_line"] = update[1]
                copy_update["mode"] = "c_update"
                copy_update["edit_action"] = edit_action
                copy_update["str_diff"] = construct_str_diff_data(src_lines[update[0]], added_lines[update[1]])
                edit_scripts.append(copy_update)
        elif km_match['mode'] == 'r':
            edit_action = generate_edit_action("move", km_match['block_length'], km_match['src_start'],
                                               km_match['added_start'], added_lines[km_match['added_start']][1][0] -
                                               src_lines[km_match['src_start']][1][0])
            edit_scripts.append({"src_line": km_match['src_start'], "block_length": km_match['block_length'],
                                 "dest_line": km_match['added_start'], "mode": "move",
                                 "indent_offset": km_match['indent_diff'], "edit_action": edit_action,
                                 "move_type": km_match['move_type'],
                                 "updates": km_match['updates']})
            for bl in range(km_match['block_length']):
                r_line_no = km_match['src_start'] + bl
                i_line_no = km_match['added_start'] + bl
                diff_scripts_dict['r' + str(r_line_no)] = "move-" + str(km_match['added_start'])
                diff_scripts_dict['i' + str(i_line_no)] = "move-" + str(km_match['src_start'])
            # 块移动更新
            for update in km_match['updates']:
                move_update = {}
                edit_action = generate_edit_action("m_update", update[0], update[1])
                move_update["src_line"] = update[0]
                move_update["dest_line"] = update[1]
                move_update["mode"] = "m_update"
                move_update["edit_action"] = edit_action
                move_update["str_diff"] = construct_str_diff_data(src_lines[update[0]], added_lines[update[1]])
                edit_scripts.append(move_update)
        elif km_match['mode'] == 'u':
            edit_action = generate_edit_action("update", km_match['src_start'], km_match['added_start'],
                                               added_lines[km_match['added_start']][1][0] -
                                               src_lines[km_match['src_start']][1][0])
            diff_scripts_dict['r' + str(km_match['src_start'])] = "update-" + str(km_match['added_start'])
            diff_scripts_dict['i' + str(km_match['added_start'])] = "update-" + str(km_match['src_start'])
            edit_scripts.append({"src_line": km_match['src_start'],
                                 "dest_line": km_match['added_start'], "mode": "update", "str_diff":
                                     construct_str_diff_data(src_lines[km_match['src_start']],
                                                             added_lines[km_match['added_start']]),
                                 "indent_offset": added_lines[km_match['added_start']][1][0] -
                                                  src_lines[km_match['src_start']][1][0],
                                 "edit_action": edit_action})
    for hunk in hunks:
        if not hunk[0]:
            hunk_last_index = diff_scripts.index("i" + str(hunk[1][-1]))
            if hunk_last_index == len(diff_scripts) - 1:
                src_line_no = src_len + 1
            else:
                src_line_no = int(diff_scripts[hunk_last_index + 1][1:])
            for i_line_no in hunk[1]:
                if not diff_scripts_dict["i" + str(i_line_no)]:
                    diff_scripts_dict['i' + str(i_line_no)] = "insert"
                    edit_action = generate_edit_action("insert", i_line_no)
                    edit_scripts.append(
                        {"mode": "insert", "dest_line": i_line_no, "src_line": src_line_no, "edit_action": edit_action})
        elif not hunk[1]:
            hunk_last_index = diff_scripts.index("r" + str(hunk[0][-1]))
            if hunk_last_index == len(diff_scripts) - 1:
                dest_line_no = dest_len + 1
            else:
                dest_line_no = k_pairs[int(diff_scripts[hunk_last_index + 1][1:])]
            for r_line_no in hunk[0]:
                if not diff_scripts_dict['r' + str(r_line_no)]:
                    diff_scripts_dict['r' + str(r_line_no)] = "delete"
                    edit_action = generate_edit_action("delete", r_line_no)
                    edit_scripts.append(
                        {"mode": "delete", "dest_line": dest_line_no, "src_line": r_line_no,
                         "edit_action": edit_action})
        else:
            hunk_last_index = diff_scripts.index("i" + str(hunk[1][-1]))
            if hunk_last_index == len(diff_scripts) - 1:
                cur_left_line = src_len + 1
                cur_right_line = dest_len + 1
            else:
                cur_left_line = int(diff_scripts[hunk_last_index + 1][1:])
                cur_right_line = k_pairs[cur_left_line]
            for rs in hunk[1][::-1]:
                cur_i_ds = 'i' + str(rs)
                if not diff_scripts_dict[cur_i_ds]:
                    diff_scripts_dict[cur_i_ds] = "insert"
                    edit_action = generate_edit_action("insert", rs)
                    edit_scripts.append(
                        {"mode": "insert", "dest_line": rs, "src_line": cur_left_line,
                         "edit_action": edit_action})
                    cur_right_line = rs
                else:
                    s_line_no = int(diff_scripts_dict[cur_i_ds].split("-")[1])
                    if diff_scripts_dict[cur_i_ds].startswith("update") or diff_scripts_dict[cur_i_ds].startswith(
                            "split") or diff_scripts_dict[cur_i_ds].startswith("merge") or (
                            diff_scripts_dict[cur_i_ds].startswith("move") and (
                            int(diff_scripts_dict['r' + str(s_line_no)].split("-")[1]) - cur_right_line) == (
                                    s_line_no - cur_left_line)):
                        cur_left_line = s_line_no
                        cur_right_line = int(diff_scripts_dict['r' + str(s_line_no)].split("-")[1])
                    else:
                        cur_right_line = rs
            if hunk_last_index == len(diff_scripts) - 1:
                cur_left_line = src_len + 1
                cur_right_line = dest_len + 1
            else:
                cur_left_line = int(diff_scripts[hunk_last_index + 1][1:])
                cur_right_line = k_pairs[cur_left_line]
            for ls in hunk[0][::-1]:
                cur_r_ds = 'r' + str(ls)
                if not diff_scripts_dict[cur_r_ds]:
                    diff_scripts_dict[cur_r_ds] = "delete"
                    edit_action = generate_edit_action("delete", ls)
                    edit_scripts.append(
                        {"mode": "delete", "dest_line": cur_right_line, "src_line": ls,
                         "edit_action": edit_action})
                    cur_left_line = ls
                else:
                    d_line_no = int(diff_scripts_dict[cur_r_ds].split("-")[1])
                    if diff_scripts_dict[cur_r_ds].startswith("update") or diff_scripts_dict[cur_r_ds].startswith(
                            "split") or diff_scripts_dict[cur_r_ds].startswith("merge") or (
                            diff_scripts_dict[cur_r_ds].startswith("move") and (
                            int(diff_scripts_dict['i' + str(d_line_no)].split("-")[1]) - cur_left_line) == (
                                    d_line_no - cur_right_line)):
                        cur_right_line = d_line_no
                        cur_left_line = int(diff_scripts_dict['i' + str(d_line_no)].split("-")[1])
    edit_scripts.sort(key=lambda x: (x["src_line"], x["dest_line"]))
    # 修复
    # 修复删除新增交叉问题
    for esr in edit_scripts:
        if esr['mode'] == 'delete':
            for esi in edit_scripts:
                if esi['mode'] == "insert" and esr["dest_line"] > esi['dest_line'] and esr["src_line"] < esi[
                    'src_line']:
                    esr['dest_line'] = esi['dest_line']
    return edit_scripts


def generate_edit_scripts_from_diff(diff_scripts):
    # 从diff_scripts生成
    src_line_no = 1
    dest_line_no = 1
    edit_scripts = []
    for diff_script in diff_scripts:
        if diff_script[0] == 'r':
            r_line_no = int(diff_script[1:])
            edit_action = generate_edit_action("delete", r_line_no)
            edit_scripts.append(
                {"mode": "delete", "src_line": r_line_no, "dest_line": dest_line_no, "edit_action": edit_action})
            src_line_no += 1
        elif diff_script[0] == 'i':
            i_line_no = int(diff_script[1:])
            edit_action = generate_edit_action("insert", i_line_no)
            edit_scripts.append({"mode": "insert", "dest_line": int(diff_script[1:]), "src_line": src_line_no,
                                 "edit_action": edit_action})
            dest_line_no += 1
        else:
            src_line_no += 1
            dest_line_no += 1
    return edit_scripts


def exists_hunk_inter(changes):
    for key, value in changes.items():
        if value:
            return key
    return None


def mapping_line_update(src_lines_list, dest_lines_list, hunks, ctx_length, line_sim_weight, sim_threshold):
    change_diffs = []
    for hunk in hunks:
        if hunk[0] and hunk[1]:
            changes = OrderedDict()
            # 找出所有的映射对
            for r_line_no in hunk[0]:
                for i_line_no in hunk[1]:
                    syn_sim = W_BESTI_LINE(r_line_no, i_line_no, src_lines_list, dest_lines_list, ctx_length,
                                           line_sim_weight,
                                           sim_threshold)
                    if syn_sim[0]:
                        changes[(r_line_no, i_line_no, 1 - syn_sim[1])] = []
            for change1 in changes:
                for change2 in changes:
                    if change1 == change2:
                        continue
                    if (change2[1] - change1[1]) * (change2[0] - change1[0]) < 0:
                        changes[change1].append(change2)
            changes = OrderedDict(sorted(changes.items(), key=lambda x: (len(x[1]), x[0][2])))
            # 重复删除交叉点最多的线，第二排序为综合相似度
            while changes:
                if list(changes.items())[-1][1]:
                    last_item = changes.popitem()
                    for change in changes:
                        if last_item[0] in changes[change]:
                            changes[change].remove(last_item[0])
                if not exists_hunk_inter(changes):
                    break
                changes = OrderedDict(sorted(changes.items(), key=lambda x: (len(x[1]), x[0][2])))
            for src_line_no, dest_line_no, syn_sim in sorted(changes.keys(), key=lambda x: x[0]):
                # ctx_sim = context_similarity(src_line_no, dest_line_no, 1,
                #                              src_lines_list,
                #                              dest_lines_list)
                change_diffs.append(
                    {'src_start': src_line_no, 'added_start': dest_line_no, 'context_similarity': None, 'mode': 'u',
                     'block_length': 1, 'weight': 1 + syn_sim / 10})
    return change_diffs


def remove_move_changes_from_diffs(km_matches, diff_scripts):
    diffs = diff_scripts[:]
    for km_match in km_matches:
        for line in range(km_match['block_length']):
            if km_match['mode'] == 'r':
                diffs.remove('r' + str(km_match['src_start'] + line))
                diffs.remove('i' + str(km_match['added_start'] + line))
            if km_match['mode'] == 'k':
                diffs.remove('i' + str(km_match['added_start'] + line))
    return diffs


def identify_splits_per_hunk(hunk, src_lines, added_lines, max_split_lines=8):
    # 每次看是否startswith，是的话，去掉前段，再判断是否startswith，知道判断remain相等，空行不计数max_split_lines
    results = []
    left_lines = hunk[0]
    right_lines = hunk[1]
    traverse_start = right_lines[0]
    for left_line_no in left_lines:
        blank_first_line = True
        left_line = src_lines[left_line_no][0].strip()
        right_line_no_start = traverse_start
        cur_right_line_no = right_line_no_start
        if cur_right_line_no not in added_lines:
            break
        cur_right_line = added_lines[cur_right_line_no][0].strip()
        lines = 1
        if not right_lines:
            break
        while cur_right_line_no <= right_lines[-1]:
            if cur_right_line == "":
                if blank_first_line:
                    right_line_no_start += 1
                cur_right_line_no += 1
                if cur_right_line_no not in right_lines or cur_right_line_no > right_lines[-1]:
                    break
                cur_right_line = added_lines[cur_right_line_no][0].strip()
                continue
            if cur_right_line == left_line and lines > 1:
                results.append([[left_line_no], list(range(right_line_no_start, cur_right_line_no + 1))])
                # 删除已匹配的
                for split_line in range(right_line_no_start, cur_right_line_no + 1):
                    right_lines.remove(split_line)
                    # diff_scripts.remove('i' + str(split_line))
                    del added_lines[split_line]
                del src_lines[left_line_no]
                traverse_start = cur_right_line_no + 1
                # diff_scripts.remove('r' + str(left_line_no))
                break
            elif left_line.startswith(cur_right_line) and lines <= max_split_lines:
                blank_first_line = False
                left_line = left_line[len(cur_right_line):].lstrip()
                cur_right_line_no += 1
                if cur_right_line_no > right_lines[-1]:
                    break
                if cur_right_line_no not in added_lines:
                    while cur_right_line_no not in added_lines and cur_right_line_no <= right_lines[-1]:
                        cur_right_line_no += 1
                    if cur_right_line_no == right_lines[-1]:
                        break
                    else:
                        right_line_no_start = cur_right_line_no
                        cur_right_line = added_lines[right_line_no_start][0].strip()
                        left_line = src_lines[left_line_no][0].strip()
                        lines = 1
                else:
                    cur_right_line = added_lines[cur_right_line_no][0].strip()
                    lines += 1
            else:
                if cur_right_line_no == right_lines[-1]:
                    break
                else:
                    if right_line_no_start == cur_right_line_no:
                        right_line_no_start += 1
                        if right_line_no_start not in added_lines:
                            while right_line_no_start not in added_lines and right_line_no_start <= right_lines[-1]:
                                right_line_no_start += 1
                            if right_line_no_start == right_lines[-1]:
                                break
                        cur_right_line_no = right_line_no_start
                    else:
                        right_line_no_start = cur_right_line_no
                    cur_right_line = added_lines[right_line_no_start][0].strip()
                    left_line = src_lines[left_line_no][0].strip()
                    lines = 1
    # 删除left
    for result in results:
        left_lines.remove(result[0][0])
    return results


def identify_merges_per_hunk(hunk, src_lines, added_lines, max_merge_lines=8):
    # 每次看是否startswith，是的话，去掉前段，再判断是否startswith，知道判断remain相等
    results = []
    left_lines = hunk[0]
    right_lines = hunk[1]
    traverse_start = left_lines[0]
    for right_line_no in right_lines:
        right_line = added_lines[right_line_no][0].strip()
        left_line_no_start = traverse_start
        cur_left_line_no = left_line_no_start
        if cur_left_line_no not in src_lines:
            break
        cur_left_line = src_lines[cur_left_line_no][0].strip()
        lines = 1
        if not left_lines:
            break
        while cur_left_line_no <= left_lines[-1]:
            if cur_left_line == "":
                cur_left_line_no += 1
                if cur_left_line_no not in left_lines or cur_left_line_no > left_lines[-1]:
                    break
                cur_left_line = src_lines[cur_left_line_no][0].strip()
                continue
            if cur_left_line == right_line:
                if lines > 1:
                    results.append([list(range(left_line_no_start, cur_left_line_no + 1)), [right_line_no]])
                    # 删除已匹配的
                    for split_line in range(left_line_no_start, cur_left_line_no + 1):
                        left_lines.remove(split_line)
                        # diff_scripts.remove('r' + str(split_line))
                        del src_lines[split_line]
                    del added_lines[right_line_no]
                    # diff_scripts.remove('i' + str(right_line_no))
                    traverse_start = cur_left_line_no + 1
                    break
                else:
                    if left_line_no_start == cur_left_line_no:
                        left_line_no_start += 1
                        cur_left_line_no = left_line_no_start
                    else:
                        left_line_no_start = cur_left_line_no
                    if cur_left_line_no not in src_lines:
                        break
                    cur_left_line = src_lines[cur_left_line_no][0].strip()
                    right_line = added_lines[right_line_no][0].strip()
                    lines = 1
            elif right_line.startswith(cur_left_line) and lines <= max_merge_lines:
                right_line = right_line[len(cur_left_line):].lstrip()
                cur_left_line_no += 1
                if cur_left_line_no > left_lines[-1]:
                    break
                if cur_left_line_no not in src_lines:
                    while cur_left_line_no not in src_lines and cur_left_line_no <= left_lines[-1]:
                        cur_left_line_no += 1
                    if cur_left_line_no == left_lines[-1]:
                        break
                    else:
                        left_line_no_start = cur_left_line_no
                        cur_left_line = src_lines[left_line_no_start][0].strip()
                        right_line = added_lines[right_line_no][0].strip()
                        lines = 1
                else:
                    cur_left_line = src_lines[cur_left_line_no][0].strip()
                    lines += 1
            else:
                if cur_left_line_no == left_lines[-1]:
                    break
                else:
                    if left_line_no_start == cur_left_line_no:
                        left_line_no_start += 1
                        if left_line_no_start not in src_lines:
                            while left_line_no_start not in src_lines and left_line_no_start <= left_lines[-1]:
                                left_line_no_start += 1
                            if left_line_no_start == left_lines[-1]:
                                break
                        cur_left_line_no = left_line_no_start
                    else:
                        left_line_no_start = cur_left_line_no
                    cur_left_line = src_lines[cur_left_line_no][0].strip()
                    right_line = added_lines[right_line_no][0].strip()
                    lines = 1
    # 删除right
    for result in results:
        right_lines.remove(result[1][0])
    return results


def mapping_merges(hunks, src_lines, added_lines, max_merge_lines):
    merges = []
    for hunk in hunks:
        if hunk[0] and hunk[1]:
            # 设置函数专门判断
            if len(hunk[0]) > 1:
                merges_per_hunk = identify_merges_per_hunk(hunk, src_lines, added_lines, max_merge_lines)
                merges = merges + merges_per_hunk
    return merges


def mapping_splits(hunks, src_lines, added_lines, max_split_lines):
    splits = []
    for hunk in hunks:
        if hunk[0] and hunk[1]:
            # 设置函数专门判断
            if len(hunk[1]) > 1:
                splits_per_hunk = identify_splits_per_hunk(hunk, src_lines, added_lines, max_split_lines)
                splits = splits + splits_per_hunk
    return splits


def copy_block_in_hunk(copy_block, hunks):
    for hunk in hunks:
        src = hunk[0]
        dest = hunk[1]
        if src and dest and copy_block["src_start"] >= src[0] and copy_block["src_start"] + copy_block[
            "block_length"] - 1 <= src[-1] and copy_block["added_start"] >= dest[0] and copy_block["added_start"] + \
                copy_block["block_length"] - 1 <= dest[-1]:
            return True
    return False


def relative_distance(src_line, dest_line, block_length, diff_scripts):
    # relative distance: src_line与dest_line中间有多少‘k’以及增和删的最大值
    src_index = 0
    while src_index < len(diff_scripts):
        if diff_scripts[src_index] == 'k' + str(src_line) or diff_scripts[src_index] == 'r' + str(src_line):
            break
        src_index += 1
    dest_index = diff_scripts.index('i' + str(dest_line))
    k_num, i_num, r_num = 0, 0, 0
    if src_index <= dest_index:
        for i in range(src_index + block_length, dest_index):
            if diff_scripts[i][0] == 'k':
                k_num += 1
            elif diff_scripts[i][0] == 'i':
                i_num += 1
            elif diff_scripts[i][0] == 'r':
                r_num += 1
    else:
        for i in range(dest_index + block_length, src_index):
            if diff_scripts[i][0] == 'k':
                k_num += 1
            elif diff_scripts[i][0] == 'i':
                i_num += 1
            elif diff_scripts[i][0] == 'r':
                r_num += 1
    return k_num + max(r_num, i_num)

def myers_diff(src, dest):
    with open(src, "r", encoding="utf8") as lf:
        left_lines = [line.rstrip() for line in lf]
    with open(dest, "r", encoding="utf8") as rf:
        right_lines = [line.rstrip() for line in rf]
    return myers.diff(left_lines, right_lines)

def BDiff(src, dest, src_lines_list, dest_lines_list, diff_algorithm="Histogram", indent_tabs_size=4,
          min_move_block_length=2, min_copy_block_length=2, ctx_length=4, line_sim_weight=0.6,
          sim_threshold=0.5, max_merge_lines=8, max_split_lines=8, pure_mv_block_contain_punc=False,
          pure_cp_block_contain_punc=False, count_mv_block_update=True, count_cp_block_update=True, identify_move=True,
          identify_copy=True, identify_update=True, identify_split=True, identify_merge=True):
    # edit effort: move、copy后有缩进+1，有n行update就+n
    # update: 1
    # move: 2
    # copy: 4
    # 1）运行Myers算法：Myers.diff，得到编辑脚本以及带行号的增、减行：Myers.diff_line_no
    env = os.environ.copy()
    # os.chdir(r'/app')
    # 将 /usr/bin 添加到 PATH 中，确保 git 命令可用

    env["PATH"] = "/usr/bin:" + env["PATH"]
    diffs = myers_diff(src, dest)
    # 2）构建数据列表：construct_line_data
    src_lines, added_lines, diff_scripts, hunks = construct_line_data(diffs, indent_tabs_size)
    src_lines_copy = src_lines.copy()
    # 3）先生成change的变更
    if added_lines:
        move_mappings, copy_mappings, splits, merges, update_mappings = [], [], [], [], []
        hunks_copy = copy.deepcopy(hunks)
        if identify_split:
            splits = mapping_splits(hunks, src_lines, added_lines, max_split_lines)
        if identify_merge:
            merges = mapping_merges(hunks, src_lines, added_lines, max_merge_lines)
        splits_merges = splits + merges
        if identify_move:
            move_mappings = mapping_block_move(src_lines, added_lines, src_lines_list, dest_lines_list,
                                               min_move_block_length, diff_scripts, pure_mv_block_contain_punc,
                                               count_mv_block_update)
        if identify_copy:
            copy_mappings = mapping_block_copy(src_lines_copy, added_lines, src_lines_list, dest_lines_list,
                                               min_copy_block_length, hunks, diff_scripts, pure_cp_block_contain_punc,
                                               count_cp_block_update)
        if identify_update:
            update_mappings = mapping_line_update(src_lines_list, dest_lines_list, hunks, ctx_length, line_sim_weight,
                                                  sim_threshold)
        # 去除分割合并与更新的交叉的情形，删更新
        update_mappings_copy = update_mappings[:]
        for split_merge in splits_merges:
            for update_change in update_mappings_copy:
                if (split_merge[0][0] - update_change['src_start']) * (
                        split_merge[1][0] - update_change['added_start']) < 0 and update_change in update_mappings:
                    update_mappings.remove(update_change)
        all_mappings = move_mappings[:]
        for copy_mapping in copy_mappings:
            for move_mapping in move_mappings:
                if copy_mapping['src_start'] == move_mapping['src_start'] and copy_mapping['added_start'] == \
                        move_mapping['added_start'] and copy_mapping['block_length'] == move_mapping['block_length']:
                    break
            else:
                all_mappings.append(copy_mapping)
        all_mappings = all_mappings + update_mappings
        km_matches = []
        if all_mappings:
            km_matches, remaining_mappings = km_compute(all_mappings, src_lines_list, dest_lines_list,
                                                        min_move_block_length, min_copy_block_length,
                                                        pure_mv_block_contain_punc, pure_cp_block_contain_punc)
            while remaining_mappings:
                additional_matches, remaining_mappings = km_compute(remaining_mappings, src_lines_list, dest_lines_list,
                                                                    min_move_block_length, min_copy_block_length,
                                                                    pure_mv_block_contain_punc,
                                                                    pure_cp_block_contain_punc)
                km_matches = km_matches + additional_matches
            km_matches.sort(key=lambda x: x['src_start'])
        edit_scripts = generate_edit_scripts_from_match(km_matches, diff_scripts, src_lines_copy, added_lines,
                                                        splits_merges, hunks_copy, len(src_lines_list),
                                                        len(dest_lines_list))
        return edit_scripts
    edit_scripts = generate_edit_scripts_from_diff(diff_scripts)
    return edit_scripts


def BDiffFile(src, dest):
    src_infile = open(src, 'r', encoding='utf-8')
    dest_infile = open(dest, 'r', encoding='utf-8')
    src_lines_list = src_infile.read().splitlines()
    dest_lines_list = dest_infile.read().splitlines()
    src_infile.close()
    dest_infile.close()
    pprint(BDiff(src, dest, src_lines_list, dest_lines_list))
