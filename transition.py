import sys, os
from Bio import PDB
from Bio.PDB import PDBIO
import datetime
from logging import getLogger, FileHandler, DEBUG, StreamHandler

class Transition():

    def __init__(self, logger, pdb_1, pdb_2):
        self.logger = logger
        self.print_logo()
        self.transition(pdb_1, pdb_2)


    def print_logo(self):
        self.logger.debug("""
    ---------------------------------------------
                       Transition
    ---------------------------------------------
    """)
        self.logger.debug("\n \nTransition started from " + str(dt))

    def transition(self, pdb_1, pdb_2):
        pdb_parser = PDB.PDBParser()
        structure_1 = pdb_parser.get_structure('X', pdb_1)
        structure_2 = pdb_parser.get_structure('Y', pdb_2)
        structure_moved = pdb_parser.get_structure('Z', pdb_1)

        chain_A_model_1 = structure_1[0]["A"]
        chain_A_model_2 = structure_2[0]["A"]

        for i in range(0, 101):
            for residue_1 in chain_A_model_1:
                res_id_1 = residue_1.get_full_id()[3][1]
                for residue_2 in chain_A_model_2.get_list():
                    res_id_2 = residue_2.get_full_id()[3][1]

                    if res_id_2 == res_id_1:
                        for atom_2 in residue_2.get_list():
                            atom_id = atom_2.get_id()
                            atom_2_coord = residue_2[atom_id].get_coord()

                            if atom_id in residue_1:
                                atom_1_coord = residue_1[atom_id].get_coord()
                                atom_moved_coord = (atom_2_coord - atom_1_coord)*(i/100) + atom_1_coord
                                if res_id_1 in structure_moved[0]["A"]:
                                    structure_moved[0]["A"][res_id_1][atom_id].set_coord(atom_moved_coord)

            try:
                os.makedirs("transition")
            except:
                pass
            w = PDBIO()
            w.set_structure(structure_moved)
            w.save("transition/moved_"+str(i)+".pdb")

            self.logger.debug("State_"+str(i)+" created. @"+str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")))

        body = 'for idx in range(0,100):cmd.load("transition/moved_%i.pdb"%idx,"mov")'
        f = open("pymol_transition.py", "w")
        f.write(body)
        f.close()

        self.logger.debug("transition.py finished @"+str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")))

if __name__=="__main__":
    args = sys.argv
    now = datetime.datetime.now()
    dt = now.strftime("%Y%m%d_%H%M%S")
    logfile = "transition_" + str(dt) + ".log"
    logger = getLogger(__name__)
    handler = FileHandler(filename=logfile)
    handler2 = StreamHandler()
    handler.setLevel(DEBUG)
    logger.setLevel(DEBUG)
    logger.addHandler(handler)
    logger.addHandler(handler2)
    logger.propagate = False

    pdb_1 = args[1]
    pdb_2 = args[2]

    Transition(logger, pdb_1, pdb_2)
