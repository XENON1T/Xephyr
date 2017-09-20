cp Xephyr/pacman/CMakeLists.txt .
python3.4 Xephyr/pacman/check.py
cd build
cmake ..
make
cd ..
