cd packages/
git clone https://gitlab.com/higgsbounds/higgssignals.git
cd higgssignals && git checkout 2.6.0
mkdir build && cd build 
cmake .. -DFeynHiggs_ROOT=/pMSSM_McMC/packages/FeynHiggs-2.16.1 && make
