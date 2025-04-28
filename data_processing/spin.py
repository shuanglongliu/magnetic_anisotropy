
import numpy as np
from sympy import Rational
from sympy.physics.quantum import spin, hbar, represent

### Spin and identify operators

Spin_rational = Rational(4, 1)

dim = int( 2*Spin_rational + 1 )

Sx = represent(spin.JxOp(), j=Spin_rational)/hbar
Sx = np.array(Sx, dtype=complex)

Sy = represent(spin.JyOp(), j=Spin_rational)/hbar
Sy = np.array(Sy, dtype=complex)

Sz = represent(spin.JzOp(), j=Spin_rational)/hbar
Sz = np.array(Sz, dtype=complex)

S = [Sx, Sy, Sz]

Splus = Sx + complex(0, 1)*Sy
Sminus = Sx - complex(0, 1)*Sy

ID = np.eye(dim, dtype=complex)

Spin = float( Spin_rational )
ss = Spin*(Spin + 1)


### Powers of spin operators

Sx2 = np.matmul(Sx, Sx)
Sy2 = np.matmul(Sy, Sy)
Sz2 = np.matmul(Sz, Sz)
SxSy = np.matmul(Sx, Sy)
SxSz = np.matmul(Sx, Sz)
SySz = np.matmul(Sy, Sz)
SySx = np.matmul(Sy, Sx)
SzSx = np.matmul(Sz, Sx)
SzSy = np.matmul(Sz, Sy)
Sx4 = np.matmul(Sx2, Sx2)
Sy4 = np.matmul(Sy2, Sy2)
Sz4 = np.matmul(Sz2, Sz2)
Sz6 = np.matmul(Sz4, Sz2)


if __name__ == "__main__":

    #print(ID)
    #print(ss)
    #print(Splus)
    #print(Sminus)

    pass


