import dolfin
import sys, pprint
import math as m
from pprint import pprint

initial_vertices = ( ( 0.0, 0.0, 1.0 ),
                     ( 1.0, 0.0, 0.0 ),
                     ( 0.0,-1.0, 0.0 ),
                     (-1.0, 0.0, 0.0 ),
                     ( 0.0, 1.0, 0.0 ),
                     ( 0.0, 0.0,-1.0 ) )

triangles = ( (0, 1, 2 ),
              (0, 2, 3 ),
              (0, 3, 4 ),
              (0, 4, 1 ),
              (5, 1, 4 ),
              (5, 2, 1 ),
              (5, 3, 2 ),
              (5, 4, 3 ) )

edges = ( (0, 1),
          (0, 2),
          (0, 3),
          (0, 4),
          (1, 2),
          (2, 3),
          (3, 4),
          (1, 4),
          (1, 5),
          (2, 5),
          (3, 5),
          (4, 5) )


def ordinary(ref, l) :
    x1, y1, z1 = initial_vertices[ref[0]]
    x2, y2, z2 = initial_vertices[ref[1]]
    x3, y3, z3 = initial_vertices[ref[2]]

    l1, l2, l3 = l

    x = l1*x1 + l2*x2 + l3*x3
    y = l1*y1 + l2*y2 + l3*y3
    z = l1*z1 + l2*z2 + l3*z3

    return x, y, z

def get_edge_point(a, b, f) :
    #print "Get point {} between {} and {}".format(f, a, b)
    a = dolfin.Point(*a)
    b = dolfin.Point(*b)
    e = b-a
    return a + e*f

def get_edge_vertex(a, b, i) :
    v = edge_vertices[(min(a, b), max(a, b), i if a < b else N-i)]
    print "getting edge vertex ({}, {}, {}) = {}".format(a, b, i, v)
    return v

def add_cell(editor, cell_no, v0, v1, v2) :
    print "    Adding cell {}: {}, {}, {}".format(cell_count, v0, v1, v2)
    editor.add_cell(cell_no, v0, v1, v2)

N = int(sys.argv[1])

m = dolfin.Mesh()
editor = dolfin.MeshEditor()
editor.open(m, 2, 3)
num_vertices = len(initial_vertices) + (N-1)*len(edges) + sum([i for i in range(1,N-1)])*len(triangles)
print "Num vertices: ", num_vertices
editor.init_vertices(num_vertices)
editor.init_cells(len(triangles)*N**2)

vertex_count = 0

# Add the corner vertices
for v in initial_vertices :
    print "Adding vertex {} at {}".format(vertex_count, v)
    p = dolfin.Point(*v)
    editor.add_vertex(vertex_count, p/p.norm())
    vertex_count += 1


edge_vertices = {}
# Add the "inner" vertices along the edges
for i in range(1, N) :
    for e in edges :
        edge_vertices[(e[0], e[1], i)] = vertex_count
        v = get_edge_point(initial_vertices[e[0]], initial_vertices[e[1]], float(i)/N) 
        print "Adding vertex {} at {}".format(vertex_count, v.str())
        
        editor.add_vertex(vertex_count, v/v.norm())
        vertex_count += 1

cell_count = 0
for triangle_no, triangle in enumerate(triangles) :
    print "Processing triangle ({}, {}, {})".format(*triangle)
    vertex_start = vertex_count

    for i in range(1, N) :
        print "i = {}".format(i)

        l1 = float(i)/N
        for j in range(1, N-i) :
            print "  j = {}".format(j)

            # Don't edge along initial vertices
            if i+j == N :
                continue

            l2 = float(j)/N
            l3 = 1.0 - l1 - l2
            
            p = ordinary(triangle, (l1, l2, l3))
            pp = dolfin.Point(p[0], p[1], p[2])
            pp /= pp.norm()
            print "Adding vertex {} at {}".format(vertex_count, pp.str())
            print "    ({}, {}, {}): {}".format(l1, l2, l3, pp.str())
            editor.add_vertex(vertex_count, pp)
            #vertices.append(vertex_count)
            vertex_count += 1


    # add the "corner" facets
    print "  Adding corner facets"
    add_cell(editor, cell_count,
             triangle[0],
             edge_vertices[(min(triangle[0], triangle[1]), max(triangle[0], triangle[1]), 1 if triangle[0] < triangle[1] else N-1)],
             edge_vertices[(min(triangle[0], triangle[2]), max(triangle[0], triangle[2]),  1 if triangle[0] < triangle[2] else N-1)])
    cell_count += 1

    add_cell(editor, cell_count,
             triangle[1], 
             edge_vertices[(min(triangle[1], triangle[2]), max(triangle[1], triangle[2]), 1 if triangle[1] < triangle[2] else N-1)],
             edge_vertices[(min(triangle[1], triangle[0]), max(triangle[1], triangle[0]),  1 if triangle[1] < triangle[0] else N-1)])
    cell_count += 1

    add_cell(editor, cell_count,
             triangle[2],
             edge_vertices[(min(triangle[2], triangle[0]), max(triangle[2], triangle[0]), 1 if triangle[2] < triangle[0] else N-1)],
             edge_vertices[(min(triangle[2], triangle[1]), max(triangle[2], triangle[1]), 1 if triangle[2] < triangle[1] else N-1)])
    cell_count += 1


    if N == 2 :
        print "  Refinement level 2: Adding interior cell"
        add_cell(editor, cell_count,
                 get_edge_vertex(triangle[0], triangle[1], 1),
                 get_edge_vertex(triangle[1], triangle[2], 1),
                 get_edge_vertex(triangle[2], triangle[0], 1))
        cell_count += 1

    else :
        print "  Add facets incident to original edges"
        add_cell(editor, cell_count,
                 vertex_start,
                 get_edge_vertex(triangle[2], triangle[1], 1),
                 get_edge_vertex(triangle[2], triangle[0], 1))
        cell_count += 1

        add_cell(editor, cell_count,
                 vertex_start+N-3,
                 get_edge_vertex(triangle[1], triangle[2], 1),
                 get_edge_vertex(triangle[1], triangle[0], 1))
        cell_count += 1 

        add_cell(editor, cell_count,
                 vertex_count - 1,
                 get_edge_vertex(triangle[0], triangle[2], 1),
                 get_edge_vertex(triangle[0], triangle[1], 1))
        cell_count += 1 


        print "  Add the facets along the original edges"
        # along edge(triangle[1] <--> triangle[2])
        add_cell(editor, cell_count,
                 vertex_start,
                 get_edge_vertex(triangle[2], triangle[1], 2),
                 get_edge_vertex(triangle[2], triangle[1], 1))
        cell_count += 1

        for i in range(3, N) :
            add_cell(editor, cell_count,
                     vertex_start+i-3,
                     vertex_start+i-2,
                     get_edge_vertex(triangle[2], triangle[1], i-1))
            cell_count += 1

            add_cell(editor, cell_count,
                     vertex_start+i-2,
                     get_edge_vertex(triangle[2], triangle[1], i),
                     get_edge_vertex(triangle[2], triangle[1], i-1))
            cell_count += 1
            
        # along edge(triangle[2], triangle[0])
        add_cell(editor, cell_count,
                 vertex_start,
                 get_edge_vertex(triangle[2], triangle[0], 1),
                 get_edge_vertex(triangle[2], triangle[0], 2))
        cell_count += 1

        current_inner_vertex = 0
        for i in range(3, N) :
            step = N-2-i+3
            add_cell(editor, cell_count,
                     vertex_start+current_inner_vertex+step,
                     vertex_start+current_inner_vertex,
                     get_edge_vertex(triangle[2], triangle[0], i-1))
            cell_count += 1

            add_cell(editor, cell_count,
                     get_edge_vertex(triangle[2], triangle[0], i-1),
                     get_edge_vertex(triangle[2], triangle[0], i),
                     vertex_start+current_inner_vertex+step)
            cell_count += 1
            
            current_inner_vertex += step

        # along edge(triangle[1], triangle[0])
        add_cell(editor, cell_count,
                 vertex_start+N-3,
                 get_edge_vertex(triangle[1], triangle[0], 2),
                 get_edge_vertex(triangle[1], triangle[0], 1))
        cell_count += 1

        current_inner_vertex = N-3
        for i in range(3, N) :
            step = N-3-i+3
            add_cell(editor, cell_count,
                     vertex_start+current_inner_vertex,
                     vertex_start+current_inner_vertex+step,
                     get_edge_vertex(triangle[1], triangle[0], i-3+2))
            cell_count += 1

            add_cell(editor, cell_count,
                     vertex_start+current_inner_vertex+step,
                     get_edge_vertex(triangle[1], triangle[0], i-3+2),
                     get_edge_vertex(triangle[1], triangle[0], i-3+3))
            cell_count += 1
            
            current_inner_vertex += step
   
    print "  Add the inner facets that don't touch initial edges"
    
    row_offset = 0
    for i in range(1, N-2) :
        row_length = N-1-i
        add_cell(editor, cell_count,
                 vertex_start+row_offset,
                 vertex_start+row_offset+row_length,
                 vertex_start+row_offset+1)
        cell_count += 1

        for j in range(N-3-i) :
            add_cell(editor, cell_count,
                     vertex_start+row_offset+row_length+j,
                     vertex_start+row_offset+row_length+j+1,
                     vertex_start+row_offset+j+1)
            cell_count += 1

            add_cell(editor, cell_count,
                     vertex_start+row_offset+row_length+j+1,
                     vertex_start+row_offset+j+2,
                     vertex_start+row_offset+j+1)
            cell_count += 1


        row_offset += row_length


editor.close()
print "{} {} {}".format(len(initial_vertices), (N-1)*len(edges), sum([i for i in range(1,N-2)])*len(triangles))
print num_vertices
print vertex_count
print cell_count
dolfin.info(m)

# Write off file
with open("sphere.off", "w") as f :
    f.write("OFF\n")
    f.write("{} {} 0\n\n".format(vertex_count, cell_count))
    for v in m.coordinates() :
        f.write("{} {} {}\n".format(v[0], v[1], v[2]))

    for c in m.cells() :
        f.write("3 {} {} {}\n".format(c[0], c[1], c[2]))

dolfin.plot(m)
dolfin.interactive()
