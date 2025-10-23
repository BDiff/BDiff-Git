# encoding:utf-8
import os
from bs4 import BeautifulSoup
from BDiff_Myers import BDiff
import asyncio
import sys


def save_file(filename, content, output_path="html"):
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    filepath = os.path.join(output_path, filename)
    print(f"{'覆写' if os.path.exists(filepath) else '创建'} 文件 {filename}")
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(content)
    return filepath


def get_html_template():
    with open('index.html', 'r', encoding='utf-8') as f:
        return BeautifulSoup(f.read(), 'html.parser')


def build_diff_script(src, dest, src_lines, dest_lines, edit_scripts):
    script = BeautifulSoup("", "html.parser").new_tag("script")
    content1 = repr("\n".join(src_lines))
    content2 = repr("\n".join(dest_lines))
    filename1 = repr(src)
    filename2 = repr(dest)
    script.string = (
        f"content1={content1};"
        f"filename1={filename1};"
        f"content2={content2};"
        f"filename2={filename2};"
        f"diffJson={edit_scripts};"
    )
    return str(script)


def run(src, dest, output_name=''):
    if not output_name:
        output_name = f'{os.path.splitext(os.path.basename(src))[0]}--{os.path.splitext(os.path.basename(dest))[0]}.html'

    with open(src, 'r', encoding='utf-8') as f:
        src_lines = f.read().splitlines()
    with open(dest, 'r', encoding='utf-8') as f:
        dest_lines = f.read().splitlines()
    print("Python version:", sys.version)
    print("Python executable:", sys.executable)
    edit_scripts = BDiff(os.path.abspath(src), os.path.abspath(dest), src_lines, dest_lines)
    print(edit_scripts)
    print(output_name)
    soup = get_html_template()
    soup.body.append(BeautifulSoup(build_diff_script(src, dest, src_lines, dest_lines, edit_scripts), "html.parser"))

    return save_file(output_name, soup.prettify())
