# Examples

This section contains real life code examples. This is still under construction and more example will come in the future.

# Content

```bash
├── likelihood1D    ---> Example that shows how to "hack" Xephyr for 1D unbinned likelihood
└── SR1Like	    ---> SR1 enhanced example: this not only show how to build a real life 
		         SR1 likelihood but also allows you to easily plug your own signal model, 
			 keeping the rest of SR1 models unchanged.

```

# Installation

These example (although they work from the dir within Xephyr) can be installed, so that you have your own copy (and can mess up with the code)
in the workdir. It is very easy, just do:

```bash
cd $XEPHYR_DIR
source Xephyr/pacman/installExample.sh  name_of_example
```

Where **name_of_example** parameter can be: **likelihood1D** or **SR1Like**.
