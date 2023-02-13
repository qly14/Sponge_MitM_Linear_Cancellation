from sage.misc.lazy_import import lazy_import
lazy_import('sage.geometry.polyhedron.base', 'Polyhedron')

condb0b4=[[0,1,0,1,0,1,0,1,0,0,0], 
[0,1,1,1,0,1,0,1,0,0,0], 
[0,1,1,1,0,1,1,1,0,1,0], 
[1,0,1,0,1,0,1,0,0,0,0], 
[1,0,1,1,1,0,1,0,0,0,0], 
[1,0,1,1,1,0,1,1,0,1,0], 
[1,1,0,1,0,1,0,1,0,0,0], 
[1,1,0,1,1,1,0,1,0,0,1], 
[1,1,0,1,0,1,1,1,1,0,0], 
[1,1,0,1,1,1,1,1,1,0,1], 
[1,1,1,0,1,0,1,0,0,0,0], 
[1,1,1,0,1,1,1,0,0,0,1], 
[1,1,1,0,1,0,1,1,1,0,0], 
[1,1,1,0,1,1,1,1,1,0,1], 
[1,1,1,1,1,1,1,1,0,0,0]] 

condb3=[[1,1,1,1,0],[0,1,0,1,0],[0,1,1,1,1],[1,0,1,0,0],[1,0,1,1,1]]



file=open('linear_condb3.txt','w')
poly_test = Polyhedron(vertices = condb3)
for v in poly_test.Hrepresentation():
    #print (v)
    file.write(str(v)+'\n')
file.close()
