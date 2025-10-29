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
Run the following commands in a code repository:
```
cd xxxxx 
# Exclude directory (run in the project directory, no impact on the git repository)
echo "bdiff_files/" >> .git/info/exclude  

git difftool 
```
