cp Xephyr/pacman/CMakeLists.txt .
python Xephyr/pacman/check.py
cd build
cmake ..
make
cd ..
