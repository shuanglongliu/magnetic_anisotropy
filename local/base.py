
import math
import numpy as np

def sph2cart(sph, rad=False):
    # Convert cartesian coordinate to spherical coordinate
    if rad:
        pass
    else:
        sph[1:3] = np.deg2rad(sph[1:3])
    z = sph[0] * math.cos( sph[1] )
    x = sph[0] * math.sin( sph[1] ) * math.cos( sph[2] )
    y = sph[0] * math.sin( sph[1] ) * math.sin( sph[2] )
    return [x, y, z]

def cart2sph(xyz, rad=False):
    # Convert cartesian coordinate to spherical coordinate
    sph = []
    ptsnew = appendSpherical_np(np.array([xyz]))
    sph.append(ptsnew[0][3])
    if rad:
        sph.append(ptsnew[0][4])
        sph.append(ptsnew[0][5])
    else:
        sph.append(ptsnew[0][4]/math.pi*180)
        sph.append(ptsnew[0][5]/math.pi*180)
    return sph

def get_opposite_direction_sph(direction):
    theta, phi = direction
    theta = 180 - theta
    phi = phi + 180
    phi = phi % 360
    if theta == 0 or theta == 180:
        phi = 0
    return (theta, phi)

def get_angle(a, b, rad=False):
    cosab = np.dot(a, b)/np.linalg.norm(a)/np.linalg.norm(b)
    if cosab > 1:
        cosab = 1
    elif cosab < -1:
        cosab = -1
    if rad:
        theta = math.acos(cosab)
    else:
        theta = 180*math.acos(cosab)/math.pi
    return theta

def get_angle_sph(a, b, rad=False):
    a = sph2cart(a, rad)
    b = sph2cart(b, rad)
    theta = get_angle(a, b, rad)
    return theta

def change_frame_sph(vi_sph, emati, ematf, rad=False):
    """
    emati: initial basis vectors 
    ematf: final basis vectors 
    vi_sph: the vector in the initial frame 
    vi_sph: the same vector in the final frame 
    """
    
    emati_t = np.transpose(emati)
    ematf_t = np.transpose(ematf)
    ematf_t_inv = np.linalg.inv(ematf_t)
    
    vi_cart = sph2cart(vi_sph, rad=rad)
    vf_cart = np.matmul(emati_t, vi_cart)
    vf_cart = np.matmul(ematf_t_inv, vf_cart)
    
    vf_sph = cart2sph(vf_cart, rad=rad)

    if rad:
        pi = np.pi
        twopi = 2*np.pi
    else:
        pi = 180
        twopi = 360

    if abs(vf_sph[1]) < 1.e-9 or abs(vf_sph[1]-pi) < 1.e-9:
        vf_sph[2] = 0 
    if vf_sph[2] < 0:
        vf_sph[2] = vf_sph[2] + twopi 
    elif vf_sph[2] > twopi:
        vf_sph[2] = vf_sph[2] - twopi 
    
    return vf_sph

def get_emat_local(z_direction=(0,0)):
    # let the local z axis be along the z_direction

    if z_direction == (90, 0):
        # local z axis is along the global x axis
        # swap axes: x -> z, y -> x, z -> y
        emat_local = [[0,0,1], [1,0,0], [0,1,0]]
    else:
        # define local x axis such that the global x axis is in the local xz plane
        ez = sph2cart([1, z_direction[0], z_direction[1]])
        ey = np.cross(ez, [1, 0, 0])
        ey = ey/np.linalg.norm(ey)
        ex = np.cross(ey, ez)
        emat_local = [ex, ey, ez]
    emat_local = np.array(emat_local)

    return emat_local

class magmoms:
    def __init__(self, n_site=1, magmoms=[1], directions=[[0,0]]):
        self.n_site = n_site
        if n_site != len(magmoms):
            magmoms = [1 for i in range(n_site)]
        self.magmoms = magmoms
        if n_site != len(directions):
            directions = [[0, 0] for i in range(n_site)]
        self.directions = directions
        self.get_directions_magmom_to_spin()
        self.get_m_mat_car()
        self.get_mstring()
    
    def get_directions_magmom_to_spin(self, reverse=False):
        if reverse:
            self.directions = []
            for i in range(self.n_site):
                self.directions.append( get_opposite_direction_sph(self.directions_of_spin[i]) )
        else:
            self.directions_of_spin = []
            for i in range(self.n_site):
                self.directions_of_spin.append( get_opposite_direction_sph(self.directions[i]) )

    def get_m_mat_car(self):
        m_mat_car = []
        for i in range(self.n_site):
            m_vec_car = sph2cart([self.magmoms[i], self.directions[i][0], self.directions[i][1]])
            m_mat_car.append(m_vec_car)
        self.m_mat_car = np.array(m_mat_car)

    def get_mstring(self):
        mstring = ""
        for i in range(self.n_site):
            m_vec_car = self.m_mat_car[i]
            m_vec_car = (3*"{:15.10f}").format(*m_vec_car)
            mstring = mstring + m_vec_car + "   \\\n"
        self.mstring = mstring + "   3000*0.0"

    def set_directions(self, directions):
        self.directions = directions
        self.get_directions_magmom_to_spin()
        self.get_m_mat_car()
        self.get_mstring()
        #print(self.directions)

    def set_directions_of_spin(self, directions_of_spin):
        self.directions_of_spin = directions_of_spin
        self.get_directions_magmom_to_spin(reverse=True)
        self.get_m_mat_car()
        self.get_mstring()
        #print(self.directions)

    def set_magmoms(self, magmoms):
        self.magmoms = magmoms
        self.get_m_mat_car()
        self.get_mstring()

class sphere:
    def __init__(self, n_theta=1801, n_phi=3600):
        self.n_theta = n_theta
        self.n_phi = n_phi
        self.dtheta = 180/(self.n_theta - 1)
        self.dphi = 360/self.n_phi
        self.n_direction = self.n_phi*(self.n_theta-2) + 2
        self.i_direction = 0
        self.theta = 0.0
        self.phi = 0.0

    def get_i_direction(self):
        i_theta = int( round(self.theta/self.dtheta) )
        i_phi = int( round(self.phi/self.dphi) )
        if i_phi == self.n_phi:
            i_phi = 0
        if i_theta == 0:
            i_direction = 0
        elif i_theta == self.n_theta-1:
            i_direction = self.n_direction - 1
        else:
            i_direction = self.n_phi*(i_theta-1) + i_phi + 1
        self.i_direction = i_direction

    def get_theta(self):
        if self.i_direction == 0:
            self.theta = 0.0
        elif self.i_direction == self.n_direction - 1:
            self.theta = 180.0
        else:
            i_theta = (self.i_direction - 1)//self.n_phi + 1
            self.theta = i_theta*self.dtheta

    def get_phi(self):
        if self.i_direction == 0:
            self.phi = 0.0
        elif self.i_direction == self.n_direction - 1:
            self.phi = 0.0
        else:
            i_phi = (self.i_direction - 1)%self.n_phi
            self.phi = i_phi*self.dphi

    def set_i_direction(self, i_direction=0):
        if i_direction < 0:
            i_direction = 0
        elif i_direction >= self.n_direction:
            i_direction = self.n_direction - 1
        self.i_direction = i_direction
        self.get_theta()
        self.get_phi()

    def set_angle(self, direction):
        theta, phi = direction
        if theta < 0:
            theta = 0.0
        elif theta > 180:
            theta = 180.0
        if phi < 0:
            phi = 0.0
        elif phi > 360:
            phi = 360.0
        self.theta = theta
        self.phi = phi
        self.get_i_direction()
        self.get_theta()
        self.get_phi()

    def set_n_theta_n_phi(self, n_theta=181, n_phi=360):
        self.n_theta = n_theta
        self.n_phi = n_phi
        self.dtheta = 180/(self.n_theta - 1)
        self.dphi = 360/self.n_phi
        self.n_direction = self.n_phi*(self.n_theta-2) + 2
        # Assume the same self.theta and self.phi
        self.get_i_direction()

