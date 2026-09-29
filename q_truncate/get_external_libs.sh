# Download required external libraries

mkdir -p libraries

download() {
    if command -v wget >/dev/null 2>&1; then
        wget -O "$1" "$2"
    elif command -v curl >/dev/null 2>&1; then
        curl -L --fail --retry 2 -o "$1" "$2"
    else
        echo "error: neither wget nor curl is available" >&2
        return 1
    fi
}

# EIGEN (https://eigen.tuxfamily.org/)
download eigen-3.4.0.tar.gz https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
tar -xf eigen-3.4.0.tar.gz
mv eigen-3.4.0/Eigen libraries 
rm -rf eigen-3.4.0 eigen-3.4.0.tar.gz

# PCG (https://www.pcg-random.org/)
mkdir -p libraries/PCG
download pcg-cpp-0.98.zip https://www.pcg-random.org/downloads/pcg-cpp-0.98.zip
unzip -o pcg-cpp-0.98.zip
mv pcg-cpp-0.98/include libraries/PCG/pcg-cpp-0.98
rm -rf pcg-cpp-0.98 pcg-cpp-0.98.zip
