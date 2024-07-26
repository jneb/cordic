"""CORDIC algorithms optimized for BCD arithmetic on hardware
New version that better shows the limited resources it needs
Assumes hardware with very limited functionality:
    Add/subtract (decimally) shifted values
    add/subtract tabled values
    compare with zero
that can then be used for all kinds of operations.
Operattion takes 9 * DIGITS passes through the loop
which is three additions, two shifts, and a compare

This demonstrates using decimal CORDIC for:
    multiplication
    division
    tan *, tandeg *
    sin, sindeg
    cos, cosdeg
    atan, atandeg
    hypot *
    tanh *
    sinh
    cosh
    exp
    atanh
    log
    sqrt
* requires extra multiplication or division
"""

# Ideas for future version that includes asin/acos/asinh/atanh
# make expansion estimator 1+10**(-2i)/2
# do this my only multiplying 1+10**(-2i) if j odd
# (and of course starting x just the right value to make the result work)
# then make cordic version that stops on x or y value,
# expanding the stop value with input expanded by estimator
# this is the decimal version of https://ieeexplore.ieee.org/ielx7/8920/4358609/10082948.pdf

import math, operator, random
from itertools import chain
from collections import namedtuple

# CORDIC flags
TRIG, LINEAR, HYP = -1, 0, 1
ROTATE, VECTOR = 0,1

# decimal computation values
# total digits including the first
DIGITS = 1 + 6
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
        if not s: return self
        divider = 10 ** s
        t, rem = divmod(self.value, divider)
        # round as good as possible
        if rem * 2 > divider: t += 1
        elif rem * 2 < divider: pass
        elif t % 2: t += 1
        return Fixed(t)

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
        """Allow multiplication with -1,0,1
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
        """Compute absolute accuracy with respect to correct value
        Gives a good idea of the performance
        """
        scaled = self.value / SCALE
        if not scaled: return float(DIGITS)
        return round(
            -math.log10(abs(scaled - exactValue)),
            1)

    def ulpErr(self, exactValue):
        """Compute absolute accuracy  expressed in ULP
        >>> sin(4.0).ulpErr(math.sin(4))
        11
        """
        return abs((self - Fixed(exactValue)).value)

    def __float__(self):
        return self.value / SCALE


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
    assert max(
        c2.value / c1.value
        for c1,c2 in zip(table[:-1], table[1:])
        ) <= repeat + 1, "repeat too low"
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


linear = makeCordic(LINEAR, [
    Fixed(SCALE // 10 ** i)
    for i in range(DIGITS)
    ])

def mul(a, b):
    """First argument cannot be larger than 10!
    >>> mul(1.3, 1.5).accuracy(1.95)
    5.5
    """
    return linear(ROTATE, a=Fixed(a), x=Fixed(b)).y

def div(a, b):
    """
    >>> div(.5, .7).accuracy(5 / 7)
    5.9
    """
    return linear(VECTOR, y=Fixed(a), x=Fixed(b)).a


# trig in rad
trig = makeCordic(TRIG, [
    Fixed(math.atan(10 ** -i))
   for i in range(DIGITS - 1)
   ])

# correction for trig
trigConst = Fixed(1 / math.prod(
        (1 + 100 ** -i)
        for i in range(DIGITS - 1)
        ) ** (9 / 2))

def sin(x):
    """
    >>> sin(1.).accuracy(math.sin(1))
    6.0
    """
    return trig(ROTATE, a=Fixed(x), x=trigConst).y

def cos(x):
    """
    >>> cos(1.2).accuracy(math.cos(1.2))
    5.5
    """
    return trig(ROTATE, a=Fixed(x), x=trigConst).x

def tan(x):
    """
    >>> tan(5.3).accuracy(math.tan(5.3))
    4.4
    """
    t = trig(ROTATE, a=Fixed(x))
    return div(t.y, t.x)

def atan(x):
    """
    >>> atan(5.3).accuracy(math.atan(5.3))
    5.8
    """
    return trig(VECTOR, y=Fixed(x)).a

def hypot(a, b):
    """pretty terrible accuracy, alas
    >>> hypot(.333, .444).accuracy(.555)
    5.0
    """
    return mul(
        trigConst,
        trig(VECTOR, x=Fixed(a), y=Fixed(b)).x,
        )

# trig in tens of degrees
deg = makeCordic(TRIG, [
    Fixed(math.atan(10 ** -i) * 18 / math.pi)
    for i in range(DIGITS - 0)
    ])

# correction for trig
degConst = Fixed(1 / math.prod(
        1 + (100 ** -i)
        for i in range(DIGITS - 0)
        ) ** (9 / 2))

def sindeg(x):
    """
    >>> sindeg(1.1).accuracy(math.sin(math.pi / 180 * 11))
    5.7
    """
    return deg(ROTATE, a=Fixed(x), x=degConst).y

def cosdeg(x):
    """
    >>> cosdeg(2.2).accuracy(math.cos(math.pi / 180 * 22))
    5.1
    """
    return deg(ROTATE, a=Fixed(x), x=degConst).x

def tandeg(x):
    """
    >>> tandeg(3.3).accuracy(math.tan(math.pi / 180 * 33))
    6.2
    """
    t = deg(ROTATE, a=Fixed(x))
    return div(t.y, t.x)

def atandeg(x):
    """
    >>> atandeg(5.3).accuracy(math.atan(5.3) * 18 / math.pi)
    5.0
    """
    return deg(VECTOR, y=Fixed(x)).a


# hyperbolic functions
# works only for abs(a) < 1
hyp = makeCordic(HYP, [
   Fixed(math.atanh(10 ** -i))
   for i in range(1, DIGITS - 1)
   ], repeat=9, start=1)

# correction for hyp and sqrt
hypConst = math.prod(
        (1 - 100 ** -i)
        for i in range(1, DIGITS - 1)
        ) ** (9 / 2)

sqrtConst = Fixed(hypConst ** -2 / 4)
hypConst = Fixed(1 / hypConst)

def sinh(x):
    """
    >>> sinh(.4).accuracy(math.sinh(.4))
    6.5
    """
    return hyp(ROTATE, a=Fixed(x), x=hypConst).y

def cosh(x):
    """
    >>> cosh(.2).accuracy(math.cosh(.2))
    5.6
    """
    return hyp(ROTATE, a=Fixed(x), x=hypConst).x

def tanh(x):
    """
    >>> tanh(.31).accuracy(math.tanh(.31))
    7.0
    """
    t = hyp(ROTATE, a=Fixed(x))
    return div(t.y, t.x)

def atanh(x):
    """
    >>> atanh(0.53).accuracy(math.atanh(0.53))
    6.8
    """
    return hyp(VECTOR, y=Fixed(x)).a

def exp(x):
    """
    >>> exp(.31).accuracy(math.exp(.31))
    5.0
    """
    t = hyp(ROTATE, a=Fixed(x), x=hypConst)
    return t.x + t.y

def log(x):
    """
    >>> log(1/3).accuracy(math.log(1/3))
    5.6
    >>> log(math.pi).accuracy(math.log(math.pi))
    6.9
    """
    result = hyp(VECTOR, x=F1+Fixed(x), y=F1-Fixed(x)).a
    return -result - result

def sqrt(x):
    """
    >>> sqrt(.5).accuracy(math.sqrt(.5))
    6.7
    >>> sqrt(1.9).accuracy(math.sqrt(1.9))
    5.7
    """
    return hyp(VECTOR, x=sqrtConst+Fixed(x), y=sqrtConst-Fixed(x)).x

def testRoutine(myfn, sysfn, params=1, repeat=1000):
    """
    Test CORDIC routine with different values
    and compute range and RMS relative error in digits
    >>> testRoutine(sin, math.sin)
    ... #doctest:+SKIP
    sin     -8.1.. 8.1  5.7   78%
    ([-8.083713], [8.076818], [])    """
    errors = []
    sumsq = accurate = 0
    # maximum allowed error in result to be counted
    # 4 digits absolute precision loss accepted for OK
    margin = 10 ** (3 - DIGITS)
    count = 0
    if myfn.__name__.endswith(('h','exp')): low,high = -0.9, 0.9
    elif myfn.__name__ in ('log, sqrt, hypot'): low,high = 0, 2
    else: low,high = -9, 9
    minOK = [low] * params
    maxOK = [high] * params
    worstError = 0
    while count < repeat:
        point = [round(random.uniform(low, high), DIGITS - 1) for j in range(params)]
        if myfn.__name__ == 'div':
            point[1] = abs(point[1])
        try:
            exact = sysfn(*point)
        except (ValueError, OverflowError, ZeroDivisionError):
            continue
        count += 1
        try:
            result = float(myfn(*point))
            if (error := abs((result - exact))) < margin:
                sumsq += ((result - exact) / exact) ** 2
                accurate += 1
            else:
                #if error > worstError:
                #    worstError = error
                #    print(point, error)
                for i,x in enumerate(point):
                    if minOK[i] < x < 0: minOK[i] = x
                    if maxOK[i] > x > 0: maxOK[i] = x
        except Exception:
            errors.append(point)
    if accurate and sumsq:
        precision = -math.log10(math.sqrt(sumsq) / accurate)
    else:
        precision = 0
    #precision = math.sqrt(sumsq) / accurate * SCALE
    print(
        f"{myfn.__name__:<7s}"
        f"{max(minOK):5.1f}..{min(maxOK):4.1f}"
        f"{precision:5.1f}"
        f"  {accurate / repeat:4.0%}")
    return minOK, maxOK, errors

def testPerf():
    print("Routine  range    #digits  OK")
    for fn in (
            "sin cos tan atan "
            "sindeg cosdeg tandeg atandeg "
            "sinh cosh tanh atanh "
            "exp log sqrt ".split()):
        if fn == 'atandeg':
            mfn = lambda x: getattr(math, fn[:-3])(x) / math.pi * 18
        elif fn.endswith('deg'):
            mfn = lambda x: getattr(math, fn[:-3])(x / 18 * math.pi)
        else:
            mfn = getattr(math, fn)
        testRoutine(globals()[fn], mfn, 1)
    testRoutine(mul, operator.mul, 2)
    testRoutine(div, operator.truediv, 2)
    testRoutine(hypot, math.hypot, 2)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
    testPerf()
