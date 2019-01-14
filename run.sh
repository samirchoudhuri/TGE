path="/data2/samir/Estimator/TGE_github/srcbin"
export path

$path/v2psfitsc ../visdiff.fits GV.fits input.v2psfitsm Mg.fits k2gg.fits k2gg1.fits 1
$path/binv2psfitscm GV.fits input.v2psfitsm powerspec.dat Mg.fits 
rm -rf k2gg*.fits