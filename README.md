# bdiff git difftool配置
### 环境配置
本项目采用Python 3.12.10 ,如果命令行使用的是python3和pip3 ,请将下列Python与pip指令换成相应的指令
本项目会在git
```
git clone
cd bdiff-git  #当前项目地址
pip install -r requirements.txt
```

如果环境存在问题则：
```
# 在项目目录创建 venv
python3 -m venv .venv

#激活虚拟环境
source .venv/bin/activate

```

### 配置git difftool 命令
调用 Python 来执行它（推荐方式）。

```
git config --global diff.tool bdiff

注意：注意：注意：
/path/to/bdiff/diff.py 替换为脚本的实际路径
/path/to/python3 替换为python 环境

git config --global difftool.bdiff.cmd  '/path/to/python3 /path/to/bdiff/diff.py "$LOCAL" "$REMOTE"'
git config --global difftool.prompt false

```

### 使用
在某个代码仓库下执行

```
cd xxxxx 
#屏蔽目录(在项目目录下执行，不影响 git 仓库)
echo "bdiff_files/" >> .git/info/exclude  

git difftool 
```


