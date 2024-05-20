import numpy as np



def MFA(B_NEC, B_MEAN_NEC, R_NEC):
    def magn(X):
        """
        Return the row-wise magnitude of elements in 2D array 'X' as a single-column array.
        """
        return np.reshape(np.sqrt(np.sum(X**2, axis=1)), (-1, 1))
    """ """

    # if no positional vector is given just assume the direction (0,0,-1) in NEC
    # coordinates, i.e. radial outwards (the magnitude is not necessary, only
    # the direction matters in order to compute its cross product with the mean
    # field component
    if R_NEC is None:
        R_NEC = np.zeros(B_NEC.shape)
        R_NEC[:, 2] = -1
    #print(R_NEC)
    MFA = np.full(B_NEC.shape, np.NaN)
    # create the unitary vector of the mean field
    B_MEAN = magn(B_MEAN_NEC)
    B_MEAN_UNIT = B_MEAN_NEC / B_MEAN
    # find the field along the mean field direction
    MFA[:, 2] = np.sum(B_NEC * B_MEAN_UNIT, axis=1)

    # find the direction of the azimuthal component
    B_AZIM = np.cross(B_MEAN_UNIT, R_NEC)
    B_AZIM_UNIT = B_AZIM / magn(B_AZIM)
    # find the field along the azimuthal direction
    MFA[:, 1] = np.sum(B_NEC * B_AZIM_UNIT, axis=1)

    # find the direction of the poloidal component
    B_POL = np.cross(B_AZIM_UNIT, B_MEAN_UNIT)
    B_POL_UNIT= B_POL / magn(B_POL)
    # no need to 
    MFA[:, 0] = np.sum(B_NEC * B_POL_UNIT, axis=1)
    # test that magnitude is conserved
    return MFA

