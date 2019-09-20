# Xephyr
Frequentist inference using ROOT.

## Requirements:
If you want to use Xephyr with your local ROOT, make sure you have all the requirements.
- g++ >= 4.8.1
- gcc >= 4.8.1
- cmake >= 3.9
- ROOT >= 6.0.0 but reccomended is >= 6.10.04 to get all latest jupyter notebook features.
- python >=3  (used for the package manager) 
- if you use Anaconda ROOT, you will pay for your sins :fire:.

Alternatively, none of the above are needed when using Xephyr in a Docker container.

## Xephyr Directory Structure

```bash
├── Xephyr
    ├── examples    ---> Where all supported examples are, you can install them! have a look at the docs
    ├── notebooks   ---> Where all tutorials are (as notebooks) you can just run them after setup (described below)
    ├── pacman      ---> contains all cool scripts for generating makefiles for you
    └── src	    ---> Source code.... Can't touch this! ;)
```

## Documentation
 - Readthedocs style [documentation](https://xenon1t.github.io/Xephyr/docs/), tutorials etc.
 - Doxigen style [code reference](https://xenon1t.github.io/Xephyr/class_reference/) with class hierarchy tree and description.

## Native setup using local ROOT
Setting Up XEPHYR is quite simple, have a look below. But we reccommend you setup while following the Docs.
```bash
mkdir XEPHYR_PKG
cd XEPHYR_PKG
export XEPHYR_DIR=$(pwd)/
echo "export XEPHYR_DIR=$(pwd)/" >> ~/.bashrc
git clone git@github.com:XENON1T/Xephyr.git .

# now compile and load the library (you need to load this each time)
root -l Xephyr/loadXephyr.C

#That's it, now you can run scripts that use XEPHYR classes.
```

## Isolated setup with docker
 Getting ROOT to work on jupyter notebooks can be a hassle since root needs to be compiled with the same compiler as the python being used. Since most people have many python versions and distributions on their pc, getting ROOT to see only the "right" python can be time consuming and very fragile, especially when automated environment-setup scripts such as conda activate are invlolved. For this reason it is recomended to use containers to isolate your ROOT-Xephyr-Python environment from the rest of your system. A pre-built [Docker](https://docs.docker.com/get-started/) image is available at docker hub for your convenience.
### Running the Docker image
 Make sure you have Docker [installed](https://docs.docker.com/install/) and working. Then run

```bash
docker run --rm -p 8080:8080 yossimo/xephyr:latest
```
You should see the jupyter notebook output, copy the address with the key. 


You can mount a local directory to the docker image with the -v command
```bash
docker run --rm -p 8080:8080 -v dir_local_path:/home/xephyrian/dir_name_in_docker_image yossimo/xephyr:latest
```
Here the local folder was mounted to the /home/xephyrian/ directory in the docker image since that is the default work dir. 