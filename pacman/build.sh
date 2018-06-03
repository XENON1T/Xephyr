cp Xephyr/pacman/CMakeLists.txt .
python3.4 Xephyr/pacman/check.py
mkdir -p build
cd build
cmake ..
make
cd ..
