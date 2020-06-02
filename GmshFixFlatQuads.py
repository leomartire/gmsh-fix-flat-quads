import argparse
import numpy as np
import sys

################################################################
################################################################

debug_dofix = True # Debugging switch. Leave as is.
debug_dofix__dont_add_elements = False # Debug flag, leave as is.

threshold = 5e-3
ngnod = 4
xunit = np.array([-1,1,1,-1])
zunit = np.array([-1,-1,1,1])

################################################################
################################################################

def prepareArgumentParser():
  parser = argparse.ArgumentParser(description='Try to fix the flat elements in a GMSH msh file.')
  required=parser.add_argument_group('required arguments')
  required.add_argument("-i", "--input", type=str, help=".msh file input", required=True)
  required.add_argument("-o", "--output", type=str, help=".msh file output", required=True)
  parser.add_argument("-v", "--verbose", type=int, default=0, help="Activate verbosity.", choices=[0, 1])
  return(parser)

def error(msg='Error, stopping script.'):
  sys.exit(msg)

################################################################

def read_msh_file(filename):
  if(verbose):
    print('Loading GMSH file "'+filename+'".')
  fid = open(filename, 'r')
  lines = [line.rstrip('\n') for line in fid]
  reading_nodes = False
  read_number_of_nodes = False
  reading_elements = False
  read_number_of_elements = False
  header = [];
  for l in lines:
    if(l=='$Nodes'):
      reading_nodes = True
      continue
    if(l=='$EndNodes'):
      reading_nodes = False
    if(l=='$Elements'):
      reading_elements = True
      continue
    if(l=='$EndElements'):
      reading_elements = False
    if(reading_nodes):
      if(read_number_of_nodes):
        # read next node
        nodes[c, :] = np.array(l.split(' ')[1:3]).astype(np.float)
        #print(nodes[c,:])
        c = c+1;
        #stop
      else:
        nnodes = int(l)
        nodes = np.zeros((nnodes, 2))
        c = 0
        read_number_of_nodes = True
    elif(reading_elements):
      if(read_number_of_elements):
        # read next element
        curel = np.array(l.split(' ')[1:]).astype(np.int)
        if(curel[0]==1):
          # line, 6 items
          elements[c, 0:6] = curel
        elif(curel[0]==3):
          # quadrilateral, 8 items
          elements[c, 0:8] = curel
        else:
          print(curel)
          error()
        c = c+1;
      else:
        nelements = int(l)
        elements = np.zeros((nelements, 8), dtype=int)
        c = 0
        read_number_of_elements = True
    else:
      header.append(l)
  fid.close()
  header.pop() # pop "$EndNodes"
  header.pop() # pop "$EndElements"
  if(verbose):
    print('Finished loading GMSH file "'+filename+'".')
  return (header, nodes, elements)

################################################################

def write_msh_file(filename, header, nodes, elements):
  nnodes = np.shape(nodes)[0]
  nelements = np.shape(elements)[0]
  if(verbose):
    print('Writing '+str(nnodes)+' nodes and '+str(nelements)+' elements to file "'+filename+'".')
  fid = open(filename, 'w')
  fid.write('\n'.join(header))
  fid.write('\n')
  fid.write('$Nodes\n')
  fid.write('%d\n' % nnodes)
  for i in range(nnodes):
    fid.write('%d %.6f %.6f 0\n' % (i+1, nodes[i, 0], nodes[i, 1]))
    # Note: append a "0" at the end, since we are dealing with 2D.
  fid.write('$EndNodes\n')
  fid.write('$Elements\n')
  fid.write('%d\n' % nelements)
  for i in range(nelements):
    el = elements[i, :]
    if(el[0]==1):
      # line, 6 items
      fid.write('%d %d %d %d %d %d %d\n' % (i+1, el[0], el[1], el[2], el[3], el[4], el[5]))
    elif(el[0]==3):
      # quadrialteral, 8 items
      fid.write('%d %d %d %d %d %d %d %d %d\n' % (i+1, el[0], el[1], el[2], el[3], el[4], el[5], el[6], el[7]))
    else:
      print(el)
      error()
  fid.write('$EndElements\n')
  fid.write('\n')
  fid.close()
  if(verbose):
    print('Finished writing to file "'+filename+'".')
  return(0)

################################################################

def checkCol(p1, p2, p3):
  # Check colinearity of three points.
  # Use shoelace: x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
  diag = np.abs(p1[0]*(p2[1]-p3[1])+p2[0]*(p3[1]-p1[1])+p3[0]*(p1[1]-p2[1]))
  #print(diag)
  return(diag)
def find_colinear_within_el(el, nodes):
  # Find which three points are colinear within a given element.
  if(el[0]==3):
    # quadrilateral
    v = el[-4:]
    cc012 = checkCol(nodes[v[0]-1,:],nodes[v[1]-1,:],nodes[v[2]-1, :])
    cc123 = checkCol(nodes[v[1]-1,:],nodes[v[2]-1,:],nodes[v[3]-1, :])
    cc230 = checkCol(nodes[v[2]-1,:],nodes[v[3]-1,:],nodes[v[0]-1, :])
    cc301 = checkCol(nodes[v[3]-1,:],nodes[v[0]-1,:],nodes[v[1]-1, :])
    #print(cc012); print(cc123); print(cc230); print(cc301)
    #stop
    ccgrpsort = np.argsort(np.array([cc012, cc123, cc230, cc301]))
    #print(ccgrpsort)
    if(ccgrpsort[0]==0):
      colids = [0,1,2]
    elif(ccgrpsort[0]==1):
      colids = [1,2,3]
    elif(ccgrpsort[0]==2):
      colids = [2,3,0]
    elif(ccgrpsort[0]==3):
      colids = [3,0,1]
    else:
      colids = []
    colids = np.array(colids, dtype=int)
    healthyid = np.setdiff1d(np.arange(0,4), colids)
    return(healthyid[0], colids, v, nodes[v-1, :])
  else:
    print(el)
    error()

def find_neighbours(elid, elements, queryid, nodes):
  # Find neighbours of a given element.
  el = elements[elid-1, :]
  if(el[0]==3):
    # quadrilateral
    v = el[-4:]
    #print(v)
    el_cont_v0 = np.setdiff1d(np.ndarray.flatten(np.argwhere(np.any(np.isin(elements[:,-4:], v[0]), axis=1)))+1, elid)
    el_cont_v1 = np.setdiff1d(np.ndarray.flatten(np.argwhere(np.any(np.isin(elements[:,-4:], v[1]), axis=1)))+1, elid)
    el_cont_v2 = np.setdiff1d(np.ndarray.flatten(np.argwhere(np.any(np.isin(elements[:,-4:], v[2]), axis=1)))+1, elid)
    el_cont_v3 = np.setdiff1d(np.ndarray.flatten(np.argwhere(np.any(np.isin(elements[:,-4:], v[3]), axis=1)))+1, elid)
    #print(el_cont_v0)
    #print(el_cont_v1)
    neigh_01 = np.intersect1d(el_cont_v0, el_cont_v1)
    neigh_12 = np.intersect1d(el_cont_v1, el_cont_v2)
    neigh_23 = np.intersect1d(el_cont_v2, el_cont_v3)
    neigh_30 = np.intersect1d(el_cont_v3, el_cont_v0)
    if(np.size(neigh_01)==0):
      print('  No neighbour found for edge 0-1 of element '+str(elid)+'.')
      print(el_cont_v0, el_cont_v1)
      stop
    if(np.size(neigh_12)==0):
      print('  No neighbour found for edge 1-2 of element '+str(elid)+'.')
      print(el_cont_v1, el_cont_v2)
      stop
    if(np.size(neigh_23)==0):
      print('  No neighbour found for edge 2-3 of element '+str(elid)+'.')
      print(el_cont_v2, el_cont_v3)
      stop
    if(np.size(neigh_30)==0):
      print('  No neighbour found for edge 3-0 of element '+str(elid)+'.')
      print(el_cont_v3, el_cont_v4)
      stop
    #print(neigh_01)
    #print(neigh_12)
    #print(neigh_23)
    #print(neigh_30)
    # Save also the element sharing the query id on one vertice only.
    if(queryid==0):
      query_el = np.setdiff1d(el_cont_v0, np.array([neigh_01, neigh_30]))
    elif(queryid==1):
      query_el = np.setdiff1d(el_cont_v1, np.array([neigh_01, neigh_12]))
    elif(queryid==2):
      query_el = np.setdiff1d(el_cont_v2, np.array([neigh_12, neigh_23]))
    elif(queryid==3):
      query_el = np.setdiff1d(el_cont_v3, np.array([neigh_23, neigh_30]))
    else:
      queryid
      error()
    #print(neigh_01, neigh_12, neigh_23, neigh_30, query_el)
    if(np.size(query_el)==0):
      if(verbose):
        print('  Query of elements sharing vertice with local number '+str(queryid)+' (global number '+str(v[queryid])+') has returned nothing. Element has only four neighbours. Returning empty list.')
      #baryc = bary(elid, elements, nodes)
      #print('  Check element %d at roughly (x,z)=(%f, %f), with neighbours %d %d %d %d.' % (elid, baryc[0], baryc[1], neigh_01[0], neigh_12[0], neigh_23[0], neigh_30[0]))
      #error()
      return(neigh_01[0], neigh_12[0], neigh_23[0], neigh_30[0], [])
    elif(np.size(query_el)==1):
      if(verbose):
        print('  Query of elements sharing vertice with local number '+str(queryid)+' (global number '+str(v[queryid])+') has returned one element, returning it.')
      return(neigh_01[0], neigh_12[0], neigh_23[0], neigh_30[0], query_el[0])
    else:
      print('  Query of elements sharing vertice with local number '+str(queryid)+' (global number '+str(v[queryid])+') has returned multiple elements, don''t know what to do next.')
      error()
      #return(neigh_01[0], neigh_12[0], neigh_23[0], neigh_30[0], query_el[0])
  else:
    print(el)
    error()

def bary(elid, elements, nodes):
  # Compute barycentre of a given element.
  #print(elements[elid-1,-4:])
  #print(nodes[elements[elid-1,-4:]-1,:])
  return(np.sum(nodes[elements[elid-1,-4:]-1,:], axis=0)/4)

def modify_vertice(elements, elid, oldnode, newnode):
  # Modify a given vertice of a given element, from oldnode to newnode.
  tomodify_loc = np.argwhere(elements[elid-1, -4:]==oldnode)[0][0]
  #print(elements[elid-1, -4:])
  #print(tomodify_loc)
  elements[elid-1, (-4+tomodify_loc)] = newnode
  #print(elements[elid-1, -4:])
  return(elements)

def find_common_vertex(elements, el1, el2):
  # Find the common vertices between two elements.
  return(np.intersect1d(elements[el1-1, -4:], elements[el2-1, -4:]))

def compute_jacobian_one_pt(xquad, zquad, s, t):
  # Compute the Jacobian at the reference point $(s, t)\in[-1,1]^2$.
  ONE = 1.0
  QUARTER = 0.25
  sp = s + ONE; sm = s - ONE; tp = t + ONE; tm = t - ONE
  dershape2D = np.zeros((2, ngnod))
  dershape2D[0, 0] = QUARTER * tm
  dershape2D[0, 1] = - QUARTER * tm
  dershape2D[0, 2] =  QUARTER * tp
  dershape2D[0, 3] = - QUARTER * tp
  dershape2D[1, 0] = QUARTER * sm
  dershape2D[1, 1] = - QUARTER * sp
  dershape2D[1, 2] =  QUARTER * sp
  dershape2D[1, 3] = - QUARTER * sm
  # sum of derivatives of shape functions should be zero
  if (np.abs(np.sum(dershape2D[0,:])) > 0.):
    print('Error deriv xi shape functions'); error()
  if (abs(sum(dershape2D[1,:])) > 0.):
    print('Error deriv gamma shape functions'); error()
  xxi = 0.; zxi = 0.; xgamma = 0.; zgamma = 0.
  for ia in range(ngnod):
    xelm = xquad[ia]
    zelm = zquad[ia]
    xxi = xxi + dershape2D[0,ia]*xelm
    zxi = zxi + dershape2D[0,ia]*zelm
    xgamma = xgamma + dershape2D[1,ia]*xelm
    zgamma = zgamma + dershape2D[1,ia]*zelm
  return (xxi*zgamma-xgamma*zxi)
def compute_jacobian_one_element(coords):
  # The query points are chosen as corners of the reference element, defined by xunit and zunit as global variables.
  jacs = np.zeros(ngnod)
  for i in range(ngnod):
    jacs[i] = compute_jacobian_one_pt(coords[:,0], coords[:,1], xunit[i], zunit[i])
  return (np.array([np.min(jacs), np.mean(jacs), np.max(jacs)]))

def fix_sick_element(nodes, elements, elid):
  # Main function, aimed at modifying the flat quadrilateral.
  
  el = elements[elid-1, :]
  if(verbose):
    print('Trying to heal sick element '+str(elid)+'.')
    print(el)

  # Identify which three vertices are aligned.
  (healthyid, col_ids, vertices, vertices_coord) = find_colinear_within_el(el, nodes)
  if(verbose):
    print('healthyid, col_ids, vertices, vertices_coord: ', healthyid, col_ids, vertices, vertices_coord)
    print('-----------------')

  # Grab neighbouring elements.
  (neigh_01, neigh_12, neigh_23, neigh_30, neigh_at_healthy) = find_neighbours(elid, elements, healthyid, nodes)
  if(verbose):
    print('neighbours: ', neigh_01, neigh_12, neigh_23, neigh_30, neigh_at_healthy)
    print('-----------------')
  n_neighbours_at_healthy = len(neigh_at_healthy)

  # Select the two SIDE elements to be modified.
  # Those are the neighbours on the remaining two sides, opposed to the flat sides.
  if(healthyid==0):
    neigh1 = neigh_01
    neigh2 = neigh_30
  elif(healthyid==1):
    neigh1 = neigh_01
    neigh2 = neigh_12
  elif(healthyid==2):
    neigh1 = neigh_12
    neigh2 = neigh_23
  elif(healthyid==3):
    neigh1 = neigh_23
    neigh2 = neigh_30
  else:
    error()
  if(verbose):
    print('side neighbours: ', neigh1, neigh2)
    print('-----------------')

  # Create the two new points.
  #print(bary(neigh1, elements, nodes))
  #print(np.ndarray.flatten(vertices_coord[healthyid, :]))
  np1 = (np.ndarray.flatten(vertices_coord[healthyid, :]) + 1.5*bary(neigh1, elements, nodes))/2.5
  np2 = (np.ndarray.flatten(vertices_coord[healthyid, :]) + 1.5*bary(neigh2, elements, nodes))/2.5
  #np1 = bary(neigh1, elements, nodes)
  #np2 = bary(neigh2, elements, nodes)
  nodes = np.append(nodes, np.array([np1]), axis=0)
  np1id = np.shape(nodes)[0]
  nodes = np.append(nodes, np.array([np2]), axis=0)
  np2id = np.shape(nodes)[0]
  if(verbose):
    print('new points: ', np1, np2)
    print('-----------------')

  # Modify each side neighbour according to the new points.
  #print(elements)
  elements = modify_vertice(elements, neigh1, vertices[healthyid], np1id)
  #print(elements)
  elements = modify_vertice(elements, neigh2, vertices[healthyid], np2id)
  if(verbose):
    print('----------------- finished modifying side neighbours')

  # Modify sick element: move healthy vertice to neighbour 1
  save_healthy_node = vertices[healthyid]
  elements = modify_vertice(elements, elid, save_healthy_node, np1id)
  
  # Check whether new elements are to be added to left or right, and thus whether we should flip the numbering
  if(np1[0] > nodes[save_healthy_node-1, 0]):
    # added to right, no need to flip
    flip_new_elements = False
  else:
    # added to left, need to flip
    flip_new_elements = True
  
  # Compute distance between neighbour 1 (towards which we'll move the sick element)
  # and the sick colinear points (among which we need to choose one to leave behind, i.e. the furthest)
  dist = np.zeros(3)
  for i in range(np.size(col_ids)):
    ci = col_ids[i]
    #print(ci, vertices[ci], vertices_coord[ci, :])
    dist[i] = np.linalg.norm(np1-vertices_coord[ci, :])
    #print(' ')
  #print(np1id, np1)
  loc_id_sort = np.argsort(dist)
  if(verbose):
    print("col_ids, loc_id_sort: ", col_ids, loc_id_sort)
  emptied_col_node = vertices[col_ids[loc_id_sort[-1]]]; # chose to leav aside the furthest one
  #next_remaining_col_node = vertices[col_ids[loc_id_sort[-2]]]; # next one should be the one next to the emptied one
  if(loc_id_sort[-1]==0 or loc_id_sort[-1]==2):
    # if first or last, remaining should be second
    next_remaining_col_node = vertices[col_ids[1]]
  else:
    # if middle, remaining should be second furthest
    next_remaining_col_node = vertices[col_ids[loc_id_sort[-2]]];
  # Modify sick element: move opposed sick vertice to healthy position
  elements = modify_vertice(elements, elid, emptied_col_node, save_healthy_node)
  if(verbose):
    print('----------------- finished modifying sick element')
  
  # Modify element at healthy point.
  if(n_neighbours_at_healthy==0):
    # Only four neighbours, do not modify side neighbours, we only need to create two new elements.
    newnodes1 = [save_healthy_node, next_remaining_col_node, emptied_col_node, np2id]
    shared_neigh1_neigh2 = find_common_vertex(elements, neigh1, neigh2)
    newnodes2 = [save_healthy_node, np2id, shared_neigh1_neigh2, np1id];
    # Modify slightly the healthy node, such that the new element no. 2 is not too flat, by sliding it slightly towards the aligned points.
    #print(next_remaining_col_node, save_healthy_node, shared_neigh1_neigh2)
    #print(nodes[next_remaining_col_node-1,:], nodes[save_healthy_node-1,:], nodes[shared_neigh1_neigh2-1,:])
    nodes[save_healthy_node-1,:] = (nodes[next_remaining_col_node-1,:] + nodes[save_healthy_node-1,:])/2.
  elif(n_neighbours_at_healthy==1):
    # Five neighbours.
    # Modify the 'neigh_at_healthy' (the element sharing only the healthy vertice) such that it lies sharing neigh1 and one of the new elements.
    elements = modify_vertice(elements, neigh_at_healthy, save_healthy_node, np1id)
    emptied_neigh2_node = find_common_vertex(elements, neigh_at_healthy, neigh2)
    elements = modify_vertice(elements, neigh_at_healthy, emptied_neigh2_node, save_healthy_node)
    # Get the last unmodified node from 'neigh_at_healthy', and name it neigh_at_healthy_remaining_node. It will be a node of the new element.
    neigh_at_healthy_remaining_node = np.setdiff1d(elements[neigh_at_healthy-1, -4:], np.append(find_common_vertex(elements, neigh_at_healthy, neigh1), save_healthy_node))[0]
    newnodes1 = [save_healthy_node, next_remaining_col_node, emptied_col_node, np2id]
    newnodes2 = [save_healthy_node, np2id, emptied_neigh2_node, neigh_at_healthy_remaining_node];
  else:
    error('Don''t know what to do.')
  if(verbose):
    print('----------------- finished modifying summit neighbour')

  # Create new elements.
  if(not(debug_dofix__dont_add_elements)):
    if(flip_new_elements):
      # if new elements added on the "right", flip them to keep numbering in nice order
      newnodes1 = np.flip(newnodes1)
      newnodes2 = np.flip(newnodes2)
    model = elements[elid-1,:];
    ne1 = np.copy(model)
    ne1[-4:] = newnodes1
    ne2 = np.copy(model)
    ne2[-4:] = newnodes2
    if(verbose):
      print('new element 1', ne1)
      print('new element 2', ne2)
    elements = np.append(elements, np.array([ne1]), axis=0)
    ne1id = np.shape(elements)[0]
    elements = np.append(elements, np.array([ne2]), axis=0)
    ne2id = np.shape(elements)[0]
    
    if(verbose):
      print('----------------- finished creating new elements')
  else:
    print('DEBUG: DID NOT ADD NEW ELEMENTS !!!!!!!!!!!!!')
  
  if(verbose):
    print('Finished healing sick element '+str(elid)+'.')
    print('Created nodes '+str(np1id)+' and '+str(np2id)+'.')
    jne1 = compute_jacobian_one_element(nodes[ne1[-4:]-1, :]);
    jne2 = compute_jacobian_one_element(nodes[ne2[-4:]-1, :]);
    print('Created elements '+str(ne1id)+' and '+str(ne2id)+' (min. Jacobian '+str(jne1[0])+' and '+str(jne2[0])+').')
  return (nodes, elements)

def compute_jacobians(elements, nodes):
  nelements = np.shape(elements)[0]
  if(verbose or nelements>1000):
    print('Computing Jacobians.')
  list_problematic = [];
  for elid in range(nelements):
    if(nelements>1000 and elid%int(nelements/20.)==0):
      print('  '+str(int(100*elid/nelements))+' % done.')
    el = elements[elid, :]
    if(el[0]==3):
      # compute only for quads
      #print(el[-4:])
      #print(nodes[el[-4:]-1, :])
      min_avg_max = compute_jacobian_one_element(nodes[el[-4:]-1, :])
      if(min_avg_max[0]<=threshold):
        list_problematic.append(elid+1)
  if(verbose):
    print('problematic elements found: ', list_problematic)
    print('Finished computing Jacobians.')
  return (list_problematic)

################################################################
################################################################

# Parse arguments.
parser=prepareArgumentParser()
args = parser.parse_args()

filename_inp = args.input
filename_out = args.output
verbose = args.verbose

(header, nodes, elements) = read_msh_file(filename_inp)
if(verbose):
  print(' ')
  print('Loaded nodes: ', np.shape(nodes))
  print('Loaded elements: ', np.shape(elements))
  print(' ')

if(debug_dofix):
  # Compute Jacobians and store element IDs having Jacobian<=0.
  list_problematic = compute_jacobians(elements, nodes)
  
  if(verbose):
    print(' ')
  
  # Treat problematic elements.
  for elid in list_problematic:
    (nodes, elements) = fix_sick_element(nodes, elements, elid)
  
  if(verbose):
    print(' ')

# Re-generate .msh file as output.
write_msh_file(filename_out, header, nodes, elements)
