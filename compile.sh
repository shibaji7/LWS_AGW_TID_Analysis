rm -rf `find -type d -name .ipynb_checkpoints`:
isort -rc -sl .
autoflake --remove-all-unused-imports -i -r .
isort -rc -m 3 .
black .
