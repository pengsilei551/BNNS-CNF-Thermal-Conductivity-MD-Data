import sys
sys.path.append(r"E:\python\pysimm-stable")
from pysimm import system, lmps, forcefield


def run(test=False):
    acetone = system.read_mol("chuyeqi-14.mol")
    #acetone = system.read_mol2("model.mol2")

    f = forcefield.Dreiding()
    
    acetone.apply_forcefield(f, charges='gasteiger')
    #lmps.quick_min(acetone, min_style='fire')
        
   
    acetone.write_lammps('mixture.lmps')


if __name__ == '__main__':
    run()
