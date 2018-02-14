import numpy as np

# Copied from redcal.py in heracal (https://github.com/HERA-Team/heracal/hercal/redcal.py)
def get_pos_reds(antpos, precisionFactor=1e6):
    """ Figure out and return list of lists of redundant baseline pairs. Ordered by length.
        All baselines have the same orientation with a preference for positive b_y and,
        when b_y==0, positive b_x where b((i,j)) = pos(i) - pos(j).

        Args:
            antpos: dictionary of antenna positions in the form {ant_index: np.array([x,y,z])}.
            precisionFactor: factor that when multiplied by different baseline vectors and rounded
                to integer values, gives unique integer tuples for unique baselines

        Returns:
            reds: list of lists of redundant tuples of antenna indices (no polarizations)
    """

    keys = antpos.keys()
    reds = {}
    array_is_2D = np.all(np.all(np.array(antpos.values())[:,2]==0))
    for i,ant1 in enumerate(keys):
        for ant2 in keys[i+1:]:
            delta = tuple((precisionFactor*2.0 * (np.array(antpos[ant1]) - np.array(antpos[ant2]))).astype(int))
            # Multiply by 2.0 because rounding errors can mimic changes below the grid spacing
            if delta[0] > 0 or (delta[0]==0 and delta[1] > 0) or (delta[0]==0 and delta[1]==0 and delta[2] > 0):
                bl_pair = (ant1,ant2)
            else:
                delta = tuple([-d for d in delta])
                bl_pair = (ant2,ant1)
            # Check to make sure reds doesn't have the key plus or minus rounding error
            p_or_m = (0,-1,1)
            if array_is_2D:
                epsilons = [[dx,dy,0] for dx in p_or_m for dy in p_or_m]
            else:
                epsilons = [[dx,dy,dz] for dx in p_or_m for dy in p_or_m for dz in p_or_m]
            for epsilon in epsilons:
                newKey = (delta[0]+epsilon[0], delta[1]+epsilon[1], delta[2]+epsilon[2])
                if reds.has_key(newKey):
                    reds[newKey].append(bl_pair)
                    break
            if not reds.has_key(newKey):
                reds[delta] = [bl_pair]
    orderedDeltas = [delta for (length,delta) in sorted(zip([np.linalg.norm(delta) for delta in reds.keys()],reds.keys()))]
    return [reds[delta] for delta in orderedDeltas]

def get_uvdata(files,ant,pol):
    """
    Method to extract data from the HERA GPU Correlator UV files.

    Parameters:
    -----------
    files : list of UV files to process
    ant   : the antennas to extract
    pol   : polarization to extract

    Return:
    -------
    lst   : sidereal radians that data belong to
    data  : 2d array (time and frequency) 
    """
    lst = []
    uvdata = {}

    for uvfile in files:
        print 'Reading', uvfile
        uv = ap.miriad.UV(uvfile)
        ap.scripting.uv_selector(uv, ant, pol)
        # Read data from a single UV file
        for (uvw,t,(i,j)),d in uv.all():
            bl = '%d-%d' % (i,j)
            t = uv['lst']
            if not lst or lst[-1]!=t: lst.append(uv['lst'])
            if not uvdata.has_key(bl): uvdata[bl] = []
            uvdata[bl].append(d)
        del(uv)

    return np.asarray(lst), uvdata
    

def get_corr(vis,ants):
    """ Extracts the cross-correlation of the two antennas given and 
        corrects for the antenna gains, if provided.
        
        Parameters:
        -----------
             vis: dict,    Dictionary containing all the crosscorrelations
            ant1:  int,    Antenna number (from HERA position map)
            ant2:  int,    Antenna number (cannot be same as ant1)
    
        Returns: 
            vis: complex array of cross-correlation amplitudes
    """
    tempvis = {}
    for ant1,ant2 in ants:
        if ant1==ant2:
            print ("Auto correlations not available in this version!")
            continue
        try: tempvis['%d-%d'%(ant1,ant2)] = vis['%dN-%dN'%(ant1,ant2)]#*gain[0]*np.conj(gain[1])
        except(KeyError): tempvis['%d-%d'%(ant1,ant2)] = np.conj(vis['%dN-%dN'%(ant2,ant1)])#*gain[1]*np.conj(gains[0])

    return tempvis
