# bdiff git difftool配置
### 环境配置
本项目采用Python 3.12.10 ,如果命令行使用的是python3和pip3 ,请将下列Python与pip指令换成相应的指令
本项目会在git
```
git clone

cd bdiff  #当前项目地址

# 在项目目录创建 venv
python3 -m venv .venv

#激活虚拟环境
source .venv/bin/activate
   
mkdir bdiff_files &&  mkdir bdiff_html

echo "\nbdiff_files" >> .gitignore && echo "\nbdiff_html" >> .gitignore
pip install -r requirements.txt
```

### 配置git difftool 命令
调用 Python 来执行它（推荐方式）。

```
git config --global diff.tool bdiff

注意：/path/to/bdiff/diff.py 替换为脚本的实际路径
/Library/Frameworks/Python.framework/Versions/3.12/bin/python3 替换为您的 python 环境
git config --global difftool.bdiff.cmd  '/Library/Frameworks/Python.framework/Versions/3.12/bin/python3 /path/to/bdiff/diff.py "$LOCAL" "$REMOTE"'
git config --global difftool.prompt false

```

在仓库下执行
```
git difftool 
```
git config --global difftool.bdiff.cmd  '/Library/Frameworks/Python.framework/Versions/3.12/bin/python3 /Users/feng/NUDT/bdiff-git/diff.py "$LOCAL" "$REMOTE"'


