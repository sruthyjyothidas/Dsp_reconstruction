#script to install ePIC software and Run D0 and Ds Reconstruction

# Step 1: download the container using the command below.
curl --location https://eicweb.phy.anl.gov/containers/eic_container/-/raw/master/install.sh | bash
cd eic/ and do ./eic-shell

git clone https://github.com/eic/epic.git
cd epic/
cmake -B build -S . -DCMAKE_INSTALL_PREFIX=install
cmake --build build -- install -j8

source install/bin/thisepic.sh
cd ../

# Step 2: compile reconstruction software
git clone https://github.com/eic/EICrecon.git
cd EICrecon
source /opt/detector/epic-main/bin/thisepic.sh
cmake -B build -S . -DCMAKE_INSTALL_PREFIX=install
cmake --build build --target install -- -j8
cd ../
source EICrecon/install/bin/eicrecon-this.sh

# Step 3a: copy D0 Reconstruction software
git clone https://github.com/eic/snippets.git
cd snippets/JetsAndHF/SecondaryVertex_Chi2
make

# Step 3b: copy Ds Reconstruction software
git clone https://github.com/sruthyjyothidas/Dsp_reconstruction.git
cd snippets/JetsAndHF/Dsp_Reconstruction
make

# Input files
./make_list.sh

# Running D0 macro
./analysis

# Running Ds macro
./Dsp_Reconstruction


