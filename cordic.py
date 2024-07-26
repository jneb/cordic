"""CORDIC algorithms optimized for BCD arithmetic on hardware
New version that better shows the limited resources it needs
Assumes hardware with very limited functionality:
    Add/subtract (decimally) shifted values
    add/subtract tabled values
    compare with zero
that can then be used for all kinds of operations.

This demonstrates using decimal CORDIC for:
    multiplication
    division
    octal -> decimal conversion
    decimal -> octal conversion
    tan
    sin/cos
    atan
    tanh
    sinh/cosh
    exp
    atanh
    log
"""
#NOTE: the doctests are showing the accuracy, so they all fail if you change DIGITS
#that's not a problem; I should think of doctest output that's clear and still
#does not depend on DIGITS

import math
from itertools import chain
from collections import namedtuple

# CORDIC flags
TRIG, LINEAR, HYP = -1, 0, 1
ROTATE, VECTOR = 0,1

# decimal computation values
DIGITS = 7
MAX = 10 ** DIGITS
SCALE = MAX // 10

# output from cordic routine
CordicResult = namedtuple("CordicResult", "a, x, y")

class Fixed:
    """A number in the hardware
    with limited operations
    Assumed to be fixed point #.###### in decimal
    Range is from -10 to 10,
    with DIGITS - 1 digits after the comma
    """
    def __init__(self, x):
        """
        >>> Fixed(math.pi)
        3.141593
        >>> Fixed(-3)
        -0.000003
        """
        if isinstance(x, int):
            # for simplicity, I allow internal values to overflow beyond 10
            self.value = x
        elif isinstance(x, Fixed):
            self.value = x.value
        elif isinstance(x, float):
            if abs(x) >= 10:
                raise ValueError("Overflow")
            self.value = round(x * SCALE)
        else:
            raise TypeError(f"Can't make Fixed from {type(x)}")

    def __repr__(self):
        return format(self.value / SCALE, '.' + str(DIGITS-1) + 'f')

    def __str__(self):
        return f"Fixed({self.value})"

    def __add__(self, other):
        """
        >>> Fixed(3.) + Fixed(-4.)
        -1.000000
        """
        return Fixed(self.value + other.value)

    def __sub__(self, other):
        return Fixed(self.value - other.value)

    def __rshift__(self, s):
        return Fixed(self.value // 10 ** s)

    def __lshift__(self, s):
        return Fixed(self.value * 10 ** s)

    def __gt__(self, zero):
        """
        >>> Fixed(5) > 0
        True
        """
        assert zero == 0, "Can only compare with 0"
        return self.value > 0

    def __rmul__(self, smallInt):
        """
        >>> -1 * Fixed(0.25)
        -0.250000
        """
        if smallInt == -1:
            return -self
        elif smallInt == 0:
            return F0
        elif smallInt == 1:
            return self
        else:
            raise InternalError("Can only multiply -1, 0, 1 with Fixed")

    def __neg__(self):
        """
        >>> -Fixed(.23)
        -0.230000
        """
        return Fixed(-self.value)

    def accuracy(self, exactValue):
        """Compute accuracy with respect to correct value
        """
        scaled = self.value / SCALE
        if not scaled: return float(DIGITS)
        return round(
            -math.log10(abs((scaled - exactValue) / exactValue)),
            1)


# some convenient values
F0, F1 = Fixed(0), Fixed(SCALE)


def makeCordic(kind, table, start=0, repeat=9):
    """Build a CORDIC routine for doing calculations with a given CORDIC table
    repeat is the number of operations per digit
        determined by the biggest quotient of table values, minus one
        (which means 9 is good most of the time)
    start is the index of the first table element
        normally 0, but must be >0 for hyp because atanh(1) doesn't exist
    kind determines the iteration of the basic matrix:
    kind    operation            table values     stretch factor
    TRIG    rotate vector in 2D  atan(10**-i)   sqrt(1+10**-2i)
    LINEAR  operate on y only    10**-i
    HYP     hyperbolic           atanh(10**-i)  sqrt(1-10**-2i)
    """
    def cordic(mode, *, a=F0, x=F1, y=F0):
        """The routine doing the vector rotation
        mode
            ROTATE: input in a (range -10..10), output in x,y (int)
            VECTOR: input in x,y, output in a
        (x, y) is rotated while a is adjusted according to the table.
        """
        for i,tablei in enumerate(table, start=start):
            for j in range(repeat):
                if (-y if mode else a) > 0:
                    a -= tablei
                    x,y = (x + kind * (y >> i), y + (x >> i))
                else:
                    a += tablei
                    x,y = (x - kind * (y >> i), y - (x >> i))
        return CordicResult(a, x, y)
    return cordic

__test__ = {
'trig': '''
    >>> trig = makeCordic(TRIG, [
    ...    Fixed(math.atan(10 ** -i))
    ...    for i in range(DIGITS - 1)
    ...    ])
    >>> trigConst = math.prod(
    ...    1 + 100 ** -i
    ...    for i in range(DIGITS - 1)
    ...    ) ** (9 / 2)
    >>> trig(ROTATE, a=F1).y.accuracy(math.sin(1) * trigConst)
    5.6
    >>> trig(ROTATE, a=F1).x.accuracy(math.cos(1) * trigConst)
    5.3
    >>> trig(VECTOR, y=Fixed(3.0)).a.accuracy(math.atan(3))
    5.2
    >>> trig(VECTOR, x=Fixed(3.0), y=Fixed(4.0)).x.accuracy(5 * trigConst)
    7.0
    ''',
}


def cordic(s, table, depth=9):
    """cordic routine
    One single routine for doing everything,
    using only addition, subtraction, decimal shift
    and comparison with 0
    s: step direction (-1 for trig, 0 for linear, 1 for hyp)
    Table has values for
    2, 1.1, 1.01, ..., 1+10**-DIGITS
    m: mode
    ROTATE: input in a (range -10..10), output in x,y (int)
    VECTOR: input in x,y, output in a
    """
    def f(m, a=0, x0=1, y0=0):
        x, y = int(x0*SCALE), int(y0*SCALE)
        for i,tablei in enumerate(table):
            if not tablei: continue
            for j in range(depth):
                if (a, -y)[m] > 0:
                    a -= tablei
                    x,y = (
                        x + s * y // 10 ** i,
                        y + x // 10 ** i)
                else:
                    a += tablei
                    x,y = (
                        x - s * y // 10 ** i,
                        y - x // 10 ** i)
        return a, x, y
    return f

def accuracy(approx, value):
    """Compute number of digits of accuracy
    """
    return round(
        -math.log10(abs(
            (approx - value) / value
            or 1 / SCALE
        )),
        1)
        
# Let's start simple: multiplication and division
linTable = [
    10 ** -i
    for i in range(DIGITS)
    ]
linCordic = cordic(0, linTable)

def mul(a, b):
    """
    >>> accuracy(mul(5, 7),35)
    6.7
    """
    a, x, y = linCordic(ROTATE, a=a, x0=b)
    return y / SCALE

def div(a, b):
    """
    >>> accuracy(div(1, 3), 1/3)
    6.0
    """
    a, x, y = linCordic(VECTOR, x0=b, y0=a)
    return a

octTable = [
    8 ** (DIGITS - 1 - i)
    for i in range(DIGITS)
    ]
octCordic = cordic(0, octTable)
def oct2dec(n):
    """
    >>> oct2dec(0o12345)
    ... #doctest:+SKIP
    12345
    """
    a, x, y = octCordic(ROTATE, a=n, x0=.1)
    return x

def dec2oct(n):
    """
    >>> oct(dec2oct(123456))
    ... #doctest:+SKIP
    '0o123456'
    """
    return octCordic(VECTOR, x0=.1, y0=n)

# a bit of trigonometry
trigTable = [
    math.atan(10 ** -i)
    for i in range(DIGITS)
    ]
trigCordic = cordic(-1, trigTable)
trigConst = math.prod(
    (1 + 100 ** -i)
    for i in range(DIGITS)
    ) ** (9 / 2) * SCALE

def sin(x):
    """ accuracy(sin(1), math.sin(1))
    6.1
    """
    a, x, y = trigCordic(ROTATE, a=x)
    print(a, x, y)
    return y / trigConst

def cos(x):
    """ accuracy(cos(1), math.cos(1))
    6.1
    """
    a, x, y = trigCordic(ROTATE, a=x)
    print(a, x, y)
    return x / trigConst

def tan(x):
    """
    >>> accuracy(tan(1), math.tan(1))
    ... #doctest:+SKIP
    6.5
    """
    a, x, y = trigCordic(ROTATE, a=x)
    return y / x


def hypot(a, b):
    """
    >>> accuracy(hypot(.3,.4), .5)
    ... #doctest:+SKIP
    7.3
    """
    a, x, y = trigCordic(VECTOR, x0=a, y0=b)
    return x / trigConst


def atan(x):
    """
    >>> accuracy(atan(math.sqrt(3)), math.pi / 3)
    ... #doctest:+SKIP
    6.3
    """
    a, x, y = trigCordic(VECTOR, y0=x)
    return a

# now do the hyperbolics and logarithm
hypTable = [
    math.atanh(10 ** -i) if i else 0
    for i in range(DIGITS)
    ]
hypCordic = cordic(1, hypTable, 10)
hypConst = math.prod(
    1 - 100 ** -i or 1
    for i in range(DIGITS)
    ) ** (10 / 2) * SCALE
def tanh(x):
    """
    >>> accuracy(tanh(1),  math.tanh(1))
    ... #doctest:+SKIP
    6.3
    """
    a, x, y = hypCordic(ROTATE, a=x)
    return y / x

def sinh(x):
    """
    >>> accuracy(sinh(1), math.sinh(1))
    ... #doctest:+SKIP
    5.8
    """
    a, x, y = hypCordic(ROTATE, a=x)
    return y / hypConst

def cosh(x):
    """
    >>> accuracy(cosh(1), math.cosh(1))
    ... #doctest:+SKIP
    6.0
    """
    a, x, y = hypCordic(ROTATE, a=x)
    return x / hypConst
 
def exp(x):
    """
    >>> accuracy(exp(1), math.exp(1))
    ... #doctest:+SKIP
    6.4
    """
    a, x, y = hypCordic(ROTATE, a=x/2)
    return (x + y) / (x - y)

def atanh(x):
    """
    >>> accuracy(atanh(.5), math.atanh(.5))
    ... #doctest:+SKIP
    7.0
    """
    a, x, y = hypCordic(VECTOR, y0=x)
    return a

def log(x):
    """
    >>> accuracy(log(1/3), math.log(1/3))
    ... #doctest:+SKIP
    7.0
    >>> accuracy(log(math.pi), math.log(math.pi))
    ... #doctest:+SKIP
    6.1
    """
    a, x, y =  hypCordic(VECTOR, x0=1+x, y0=1-x)
    return -2 * a

if __name__ == "__main__":
    import doctest
    doctest.testmod()
