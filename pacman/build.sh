cp Xephyr/pacman/CMakeLists.txt .
python Xephyr/pacman/check.py
mkdir build
cd build
cmake ..
make
cd ..
