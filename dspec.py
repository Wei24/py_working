from suncasa.utils import dspec2 as ds 
specdata = ds.get_dspec("/srg/ywei/data/eovsa/ms/slfcal/IDB20170906T190319-195320.ms.corrected.xx.slfcal0", specfile="./IDB20170906T190319-195320.ms.corrected.xx.slfcal0.dspec.npz", domedian=True, verbose=True, savespec=True) 
