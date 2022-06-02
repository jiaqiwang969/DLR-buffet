fluent3DMeshToFoam  0x1.msh
1. autogrid5 mesh in numeca
2. igg mesh refine in numeca
3. export into openfoam
4. change the direction `transformPoints -rotate '((1 0 0) (0 0 -1))'
transformPoints "rotate=((1 0 0) (0 0 -1))" in openfoam 9
5. createPatch
6. renumberMesh -overwrite


1. copy 01-classic in sure perodic verctor is correct.
2. checkMesh
3. transformPoints "rotate=((1 0 0) (0 0 -1))" in openfoam 9
4. createPatch -overwrite
5. flattenMesh

norenumber is fast than renumber
