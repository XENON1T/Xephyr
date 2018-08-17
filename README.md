
# Requirements:
- g++ >= 4.8.1
- gcc >= 4.8.1
- cmake >= 3.9
- ROOT >= 6.0.0 but reccomended is >= 6.10.04 to get all latest jupyter notebook features.
- python >=3  (used for the package manager) 
- if you use Anaconda ROOT, you will pay for your sins :fire:.


# Xephyr

```bash
├── Xephyr
    ├── examples    ---> Where all supported examples are, you can install them! have a look at the docs
    ├── notebooks   ---> Where all tutorials are (as notebooks) you can just run them after setup (described below)
    ├── pacman      ---> contains all cool scripts for generating makefiles for you
    └── src	    ---> Source code.... Can't touch this! ;)
```

# Userfriendly documentation, tutorials etc. at: https://xenon1t.github.io/Xephyr/docs/
# Doxigen style code reference, classes and description at: https://xenon1t.github.io/Xephyr/class_reference/

# Super quick setup
Setting Up XEPHYR is quite simple, have a look below. However we reccommend you to really go trhough the Docs!!

```bash
mkdir XEPHYR_PKG
cd XEPHYR_PKG
export XEPHYR_DIR=$(pwd)/
echo "export XEPHYR_DIR=$(pwd)/" >> ~/.bashrc
git clone git@github.com:XENON1T/Xephyr.git .

# now compile and load the library (you need to load this each time)
root -l Xephyr/loadXephyr.C

#That's it, now you can run a script that uses XEPHYR classes.
```
