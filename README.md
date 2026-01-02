# fractionalstd2inner

A **high-performance Cython + GMP**, very basic and custom implementation of Python’s built-in `fractions` library.  
It is **~1.5× to 4× faster** than the standard Python implementation and was made to be used internally by  
[`td2inner`](https://pypi.org/project/td2inner/).

> ⚠️ **This package must be compiled locally.**  
> **No prebuilt wheels are provided.**

---

## ⚠️ Important Notice (Read This First)

`fractionalstd2inner` contains a **native Cython extension** that depends on  
**GMP (GNU Multiple Precision Arithmetic Library)**.

Because of this:

- ❌ **No prebuilt wheels**
- ✅ **Development tools are required**
- ✅ **Compilation occurs during `pip install`**

This is **intentional**.

---

## 🧩 Requirements

### General
- Python **≥ 3.10**
- C / C++ compiler
- **Cython ≥ 3.2**
- **GMP development headers**

---

## INSTALLMENT STEPS

### 1 INSTALL DEPENDENCIES ASSUMING YOU ALREADY HAVE PYTHON AND PIP

#### Arch Linux
```bash
sudo pacman -S gmp base-devel cython
```

#### Ubuntu/Debian/Mint
```bash
sudo apt update
sudo apt install -y libgmp-dev build-essential python3-dev cython3
```

#### Fedora
```bash
sudo dnf install -y gcc gcc-c++ make automake gmp-devel python3-devel cython
```

#### CentOS/ RHEL
```bash
sudo yum groupinstall "Development Tools" -y
sudo yum install -y gmp-devel python3-devel cython
```

#### openSUSE
```bash
sudo zypper install -y gcc gcc-c++ make gmp-devel python3-devel cython
```

#### Windows(YOU NEED TO DOWNLOAD MSYS2) THEN
```bash
pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-make mingw-w64-x86_64-gmp mingw-w64-x86_64-python3-cython
```

#### MacOS(YOU NEED TO DOWNLOAD Homebrew(brew)) THEN
```bash
xcode-select --install
brew install gmp cython python
```
### 2 INSTALL THE LIBRARY ITSELF VIA PIP

```bash
python3 -m pip install fractionalstd2inner
```
### 3 VERIFY INSTALLATION VIA

```bash[dist](dist)
python3 -c "import fractionalstd2inner; print(fractionalstd2inner.__version__)"
```