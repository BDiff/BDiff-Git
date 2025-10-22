# encoding:utf-8
import os
import shutil
import sys
import webbrowser
import time
from output_html import run

# 目录配置
UPLOADS_DIR = "bdiff_files"
SRC_DIR = os.path.join(UPLOADS_DIR, "src")
DEST_DIR = os.path.join(UPLOADS_DIR, "dest")
HTML_DIR = os.path.abspath(os.path.join(UPLOADS_DIR, "html"))
INDEX_FILE = os.path.join(HTML_DIR, "index.html")


def ensure_dirs():
    os.makedirs(UPLOADS_DIR, exist_ok=True)
    os.makedirs(HTML_DIR, exist_ok=True)
    os.makedirs(SRC_DIR, exist_ok=True)
    os.makedirs(DEST_DIR, exist_ok=True)


def update_index(html_file, display_name):
    """更新索引页，将新的 diff 文件添加进去"""
    links = []

    # 如果 index 文件存在，先读取已有内容
    if os.path.exists(INDEX_FILE):
        with open(INDEX_FILE, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if line.startswith("<li><a href="):
                    links.append(line)

    # 添加新文件
    rel_path = os.path.relpath(html_file, HTML_DIR)
    links.append(f'<li><a href="{rel_path}" target="_blank">{display_name}</a></li>')

    # 写入 index.html
    with open(INDEX_FILE, 'w', encoding='utf-8') as f:
        f.write("<html><head><meta charset='utf-8'><title>BDiff Index</title></head><body>\n")
        f.write("<h1>Diff 文件列表</h1>\n<ul>\n")
        for link in links:
            f.write(f"{link}\n")
        f.write("</ul>\n</body></html>")


def main():
    if len(sys.argv) != 3:
        print("Usage: python diff.py <LOCAL> <REMOTE>")
        sys.exit(1)

    ensure_dirs()

    local_file = sys.argv[1]
    remote_file = sys.argv[2]

    # 获取原始文件名，去掉 latest_src_ / latest_dest_ 前缀
    local_filename = os.path.basename(local_file)
    remote_filename = os.path.basename(remote_file)

    # 保存到各自目录
    saved_local_path = os.path.join(SRC_DIR, local_filename)
    saved_remote_path = os.path.join(DEST_DIR, remote_filename)

    if os.path.abspath(local_file) != os.path.abspath(saved_local_path):
        shutil.copy2(local_file, saved_local_path)

    if os.path.abspath(remote_file) != os.path.abspath(saved_remote_path):
        shutil.copy2(remote_file, saved_remote_path)

    # 生成带时间戳的 HTML 文件名
    timestamp = int(time.time())
    base_name = f"{timestamp}_{os.path.splitext(local_filename)[0]}--{os.path.splitext(remote_filename)[0]}.html"
    html_file = os.path.join(HTML_DIR, base_name)

    # 生成单文件 diff HTML
    print("生成单文件 diff HTML=========================")
    html_file_path = run(saved_local_path, saved_remote_path, output_name=html_file)
    print(f"Diff HTML generated: {html_file_path}")

    # 更新索引页
    display_name = f"{os.path.splitext(local_filename)[0]}--{os.path.splitext(remote_filename)[0]}"
    update_index(html_file_path, display_name)
    print(f"Index page updated: {INDEX_FILE}")

    # 自动打开最新 diff
    webbrowser.open(f"file://{os.path.abspath(html_file_path)}")


if __name__ == "__main__":
    main()
