# MacOS

```
brew install python
```

Maybe restart your shell, make sure that `which python` points to homebrew's python in `/usr/local/bin`

```
pip install virtualenv
mkdir -p ~/.virtualenvs
cd ~/.virtualenvs
virtualenv aptamer-runtime
cd aptamer-runtime
source bin/activate
```

Now you should be locked inside your aptamer virtual environment
and you can install dependencies.

```
pip install --upgrade pip
pip install numpy
pip install scipy
pip install biopython
pip install python-Levenshtein
```

In order to run `predict_structures` which calls out to a program calld `rnafold` you must install
the `viennarna` homebrew package. This package also depends on `xquartz`

_you will need sudo permissions to install xquartz_

```
brew cask install xquartz
brew install viennarna
```








