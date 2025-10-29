# bdiff git difftool Configuration

### Environment Setup
This project uses Python 3.12.10. If your command line uses `python3` and `pip3`, replace the following Python and pip commands with the corresponding ones.

```bash
git clone <repository-url>
cd bdiff-git  # Current project directory
pip install -r requirements.txt
```

If you encounter environment issues:
```
# Create a virtual environment in the project directory
python3 -m venv .venv

# Activate the virtual environment
source .venv/bin/activate
```




# bdiff git difftool配置
### 环境配置
本项目采用Python 3.12.10 ,如果命令行使用的是python3和pip3 ,请将下列Python与pip指令换成相应的指令
本项目会在git
```
git clone
cd bdiff-git  #当前项目地址
pip install -r requirements.txt
```

### Configure git difftool Command
Execute it via Python (recommended method).
```
git config --global diff.tool bdiff

NOTE: NOTE: NOTE:
Replace /path/to/bdiff/diff.py with the actual path to the script
Replace /path/to/python3 with your Python environment path

git config --global difftool.bdiff.cmd  '/path/to/python3 /path/to/bdiff/diff.py "$LOCAL" "$REMOTE"'
git config --global difftool.prompt false
```

### Usage

在某个代码仓库下执行
Run the following commands in a code repository:
```
cd xxxxx 
# Exclude directory (run in the project directory, no impact on the git repository)
echo "bdiff_files/" >> .git/info/exclude  

git difftool 
```
