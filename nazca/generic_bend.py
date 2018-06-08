from scipy.optimize import fminbound, minimize_scalar
from math import sin,cos,radians,sqrt,log,pi
import nazca as nd
from nazca.util import polyline_length

def dist_perp(P, A, B):
    '''Calculate the distance of point P to a line through A and B
    (perpendicular distance).

    Args:
        P (x,y): target point
        A (x,y): first point on line
        B (x,y): second point on line

    Returns:
        distance squared
    '''
    a = [P[0] - A[0], P[1] - A[1]]
    b = [B[0] - A[0], B[1] - A[1]]
    ab = a[0]*b[0]+a[1]*b[1]
    bb = b[0]*b[0]+b[1]*b[1]
    D = [a[0] - ab / bb * b[0],
         a[1] - ab / bb * b[1]]
    return sqrt(D[0]*D[0]+D[1]*D[1])


# This implements a P-curve type connection. Original code by Daniele
# Melati and Emil Kleijn, adapted for Nazca by Xaveer Leijtens. Based on
# Francois Ladouceur and Pierre Labeye, "A New General Approach to Optical
# Waveguide Path Design", Journal of Lightwave Technology, vol. 13, no. 3,
# march 1995 The algorithm tries to construct a curve that matches the
# input and output radius, and position, while maximizing the radius along
# the curve. Patent expired.

# The radius arguments are sensitive to the sign. A positive R indicates a
# counter-clockwise direction.

X = list()
Z = list()

def getCurvature (A, B, L):
    """Return the local curvature at a given point t of parametrized curve
    P(L) = (A(L),B(L)) along the P-Curve

    Args:
        A (array): coefficients
        B (array): coefficients
        L (float): parameter point L

    Returns:
        local curvature
    """
    dx = ddx = dy = ddy = 0
    for i in range(1,6):
      # Because x=A[i]*t^i -> dx/dL = ... dy/dL works similarly
      dx += i*A[i]*L**(i-1)
      dy += i*B[i]*L**(i-1)
      # and (d^2 x)/(d^2 L) = ... dy/dt works similarly
      if i > 1:
        ddx += i*(i-1)*A[i]*L**(i-2)
        ddy += i*(i-1)*B[i]*L**(i-2)
 
    # Local curvature, see
    # http://en.wikipedia.org/wiki/Curvature#Local_expressions
    return (ddx * dy - ddy * dx) * ((dx**2 + dy**2)**(-1.5))

# InvertMatrix helper function returns the solution matrix as a function of L
def InvertMatrix(L):
    # This is the correct matrix. Please not that there is an error in the
    # original paper (entry [3,6])
    matrix = (
        (       1,       0,       0  ,        0,       0,         0),
        (       0,       1,       0  ,        0,       0,         0),
        (       0,       0,       0.5,        0,       0,         0),
        (-10/L**3, -6/L**2, -3/2/L   ,  10/L**3, -4/L**2,  1/2/L   ),
        ( 15/L**4,  8/L**3,  3/2/L**2, -15/L**4,  7/L**3, -1/L**2  ),
        ( -6/L**5, -3/L**4, -1/2/L**3,   6/L**5, -3/L**4,  1/2/L**3))
    return matrix

def curve_AB(L):
    global X, Z
    matrix = InvertMatrix(L)

    # Compute the coefficient
    A = [0, 0, 0, 0, 0, 0]
    B = [0, 0, 0, 0, 0, 0]
    for i in range(6):
        for t in range(6):
            A[i] = A[i] + matrix[i][t] * X[t]
            B[i] = B[i] + matrix[i][t] * Z[t]
    return (A, B)

def maxcurvature(L):
    # Only this many points for the maximum curvature finding.
    Npoints = 100

    global X, Z
    A, B = curve_AB(L)
    curvature = [getCurvature(A,B,L/(Npoints-1)*i) for i in range(Npoints)]
    maxcur =  abs(max(curvature, key=abs))

    # There is an ambiguity when the maximum curvature occurs at the input
    # or output position. Many paths inbetween then satisfy the curvature
    # requirement. We therefore multiply by log(length) to favor short
    # connections without affecting the result much.
    return maxcur * log(L) #  max curvature

def gb_coefficients(xya, Rin=0, Rout=0):
    global X, Z
    # In the algorithm, a negative curvature indicates a counterclockwise
    # direction, while Nazca uses the sign of the angle instead (negative
    # angle in Nazca = clockwise direction). We therefore invert the
    # arguments here.
    _Rin  = -Rin
    _Rout = -Rout
    # Input curvature (when radius larger than a nanometer)
    if abs(_Rin) > 1e-3:
        cin = 1 / _Rin
    else:
        cin = 0
    # Output curvature
    if abs(_Rout) > 1e-3:
        cout = 1 / _Rout
    else:
        cout = 0

    xin = 0
    zin = 0
    thetain = 90

    xout = xya[0] # final position
    zout = xya[1]
    thetaout = -(xya[2]-90.0) # final slope [deg]

    # Search between cartesian distance (straight line between input and
    # output) Note that this length is related to the physical length, but
    # not equal.
    minFindL = sqrt(zout**2+xout**2)/4
    # and 3*the length of a straight line
    maxFindL = sqrt(zout**2+xout**2)*3

    # Calculate initial known terms
    dxin = sin(radians(thetain))
    dzin = cos(radians(thetain))
    ddxin =  cin * cos(radians(thetain))
    ddzin = -cin * sin(radians(thetain))

    dxout = sin(radians(thetaout))
    dzout = cos(radians(thetaout))
    ddxout =  cout * cos(radians(thetaout))
    ddzout = -cout * sin(radians(thetaout))

    # Vectors of known terms
    X = [xin, dxin, ddxin, xout, dxout, ddxout]
    Z = [zin, dzin, ddzin, zout, dzout, ddzout]

    # Determine the optimum curve length so that the minimum Radius along
    # the curve is maximized
    L = fminbound(maxcurvature, minFindL, maxFindL, disp=1)
    # Get the final matrix and extract the polynomial coefficients from it
    A, B = curve_AB(L)
    # The center of the curve is defined by a polynomial with coefficients
    # A[] for x, and B[] for y and the independent variable runs from 0 to L.
    return A,B,L

def gb_point(t, A, B, L):
    x = y = 0
    for i in range(6):
        x += A[i]*(t*L)**i
        y += B[i]*(t*L)**i
    return (x, y)

def cbend_point(t, distance, offset):
    x = t * distance
    if t > 0:
        y = offset / 2 * (1 - cos(pi * t))
    else:
        y = 0
    return (x, y)

def curve2polyline(fie, xya, acc, args=()):
    '''Generate a polyline from a parameterized curve. Use a sufficient
    number of points to ensure that the deviation of the sampled curve from
    the real curve is less than the specified accuracy.

    Args:
        fie (function): the curve function that takes one parameter that
            runs from 0 to 1 as first argument. The curve starts at the
            origin at 0 angle.
        xya: the position (x,y) and angle (a) in degrees, of the end point
            of the curve.
        acc (float): desired accuracy in micrometer.
        args (tuple): additional arguments to be passed to the curve function.

    Returns:
        a polyline approximation of the curve.
    '''
    # As a starting point find the value for t where the y-coordinate is
    # equal to acc. Since the curve always starts horizontal, this gives a
    # good and easy value.
    def fun(x):
        return abs(fie(x, *args)[1]-acc)
    res = minimize_scalar(fun, bounds=(0, 0.04), method='bounded')
    dt0 = res.x
    t0 = dt0 / 2

    P1 = fie(0, *args)
    p = [P1]
    P1 = fie(t0, *args)
    p.append((P1[0],0)) # First section has to be horizontal.
    while t0 + dt0 < 1:
        P0 = P1
        P1 = fie(t0+dt0, *args)
        P2 = fie(t0+2*dt0, *args)
        d = dist_perp(P1, P0, P2)
        f = sqrt(4*acc/d)
        dt1 = max(0.8, min(f, 1.2)) * dt0
        if t0+dt1 < 1:
            P1 = fie(t0+dt1, *args)
            p.append(P1)
        t0 += dt1
        dt0 = dt1
    if t0+dt0-1 > 0.4*dt0:
        # Make room for the point with the right direction.
        p.pop()
    t0 = 1 - dt0/3
    P1 = fie(t0, *args)
    P2 = fie(1, *args)
    lv = sqrt((P2[0]-P1[0])**2 + (P2[1]-P1[1])**2)
    a = radians(xya[2])
    P1 = (P2[0]-lv*cos(a), P2[1]-lv*sin(a))
    p.append(P1)
    p.append(P2)
    return p


if __name__ == "__main__":
    # From (0,0,0) to (x,y,a)
    xya = (1000, 1000, 20)

    # Curve with different accuracy. Placed at an angle to prevent klayout
    # from removing points which are needed to check the number of points.
    A, B, L = gb_coefficients(xya, Rin=0, Rout=0)
    xy = curve2polyline(gb_point, xya, 0.1, (A, B, L))
    xy = nd.util.polyline2polygon(xy, width=2)
    nd.Polygon(layer=202, points=xy).put(0,0,10)
    xy = curve2polyline(gb_point, xya, 0.001, (A, B, L))
    xy = nd.util.polyline2polygon(xy, width=2)
    nd.Polygon(layer=203, points=xy).put(0,0,10)

    # Curve with different accuracy and curvature at start and finish.
    A, B, L = gb_coefficients(xya, Rin=300, Rout=-100)
    xy = curve2polyline(gb_point, xya, 0.1, (A, B, L))
    nd.Polyline(layer=1, width=2, points=xy).put(0,0,10)
    xy = curve2polyline(gb_point, xya, 0.001, (A, B, L))
    nd.Polyline(layer=2, width=2, points=xy).put(0,0,10)

    # Curves with different "length".
    A, B, L = gb_coefficients(xya, Rin=0, Rout=0) # sets X and Z
    for L in range(300, 2000, 100):
        A, B = curve_AB(L)
        xy = [ gb_point(t/1000, A, B, L) for t in range(1001) ]
        nd.Polyline(layer=2, width=2, points=xy).put(500,0,10)

    A, B, L = gb_coefficients(xya, Rin=300, Rout=-100) # sets X and Z
    for L in range(300, 2000, 100):
        A, B = curve_AB(L)
        xy = [ gb_point(t/1000, A, B, L) for t in range(1001) ]
        nd.Polyline(layer=3, width=2, points=xy).put(500,0,10)

    # Length calculation
    A, B, L = gb_coefficients(xya, Rin=300, Rout=-100) # sets X and Z
    xy = curve2polyline(gb_point, xya, 0.001, (A, B, L))
    length = polyline_length(xy)
    print("Length of curve: {:.3f} µm".format(length))

    # Cosine bend
    xy = curve2polyline(cbend_point, (100,20,0), 0.001, (100,20))
    nd.Polyline(layer=2, width=2, points=xy).put(-300, 0, 0)
    xy = nd.util.polyline2polygon(xy, width=2)
    nd.Polygon(layer=203, points=xy).put(-300,0,0)

    nd.export_gds()

