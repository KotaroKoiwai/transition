# transition
Create PDB models in transition between two models.

[Description]

transition.py creates 100 intermediate-models between 2 input PDB models. 
When created pymol_transition.py is loaded on pymol, the created models are loaded as single model, "mov", and can be visualized as a short movie.

[Usage]

python3 transition.py pdb_1 pdb_2


[Note]

Ver0.0 calculates intermediate-models of only chain A.
