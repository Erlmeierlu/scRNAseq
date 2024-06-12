reticulate::install_python(version = '3.11.2:latest')

#creation
reticulate::virtualenv_create("renv/python/virtualenvs/renv-python-3.11.2")
renv::use_python("renv/python/virtualenvs/renv-python-3.11.2")

#usage in script
reticulate::use_virtualenv("renv/python/virtualenvs/renv-python-3.11.2")
