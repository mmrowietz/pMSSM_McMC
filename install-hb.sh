cd packages 
git clone https://gitlab.com/higgsbounds/higgsbounds.git
cd higgsbounds && git checkout 5.9.1
rm -rf example_programs && cp -r ../example_programs .
mkdir build && cd build 
cmake .. -DFeynHiggs_ROOT=/pMSSM_McMC/packages/FeynHiggs-2.16.1 -DLEP_CHISQ=ON && make
