LINK_TO_HRLSPG_REPO = "https://github.com/alex-fricker/Hyper-Reduced-LSPG.git"
echo "export HRLSPG_INSTALL_DIR=$HOME/Codes" >> ~/.bashrc
source ~/.bashrc

cd ${HRLSPG_INSTALL_DIR} && mkdir Libraries && cd Libraries

git clone https://gitlab.com/libeigen/eigen.git
cd eigen && mkdir build && cd build
cmake ..
sudo make install

git clone LINK_TO_HRLSPG_REPO
mkdir build && cd build
cmake ..
make
