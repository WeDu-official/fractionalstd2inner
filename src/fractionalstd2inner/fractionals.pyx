# fractionals.pyx
# distutils: language=c++
# cython: language_level=3

from decimal import Decimal, localcontext
from .fractionals cimport *
cdef int _is_perfect_power(mpz_t x, unsigned long int n, mpz_t root):
        """
        Check if x is a perfect nth power.
        If yes, sets root = x^(1/n) and returns 1, else returns 0.
        """
        return mpz_root(root, x, n)  # mpz_root returns 1 if exact, 0 otherwise
cdef class Fraction:

    def __cinit__(self):
        mpq_init(self.value)

    def __dealloc__(self):
        mpq_clear(self.value)

    # -----------------------
    # Helper: set fraction from Python ints (arbitrary size)
    # -----------------------
    cdef void _set_from_bigints(self, object num, object den):
        cdef mpz_t num_mpz, den_mpz
        mpz_init(num_mpz)
        mpz_init(den_mpz)
        mpz_set_str(num_mpz, str(num).encode(), 10)
        mpz_set_str(den_mpz, str(den).encode(), 10)
        mpq_set_num(self.value, num_mpz)
        mpq_set_den(self.value, den_mpz)
        mpq_canonicalize(self.value)
        mpz_clear(num_mpz)
        mpz_clear(den_mpz)

    # -----------------------
    # Initialize from Python object
    # -----------------------
    cdef void _set_from_object(self, object a):
        cdef Fraction bf

        if isinstance(a, int):
            self._set_from_bigints(a, 1)

        elif isinstance(a, float):
            s = format(a, ".50g")
            mpq_set_str(self.value, s.encode(), 10)
            mpq_canonicalize(self.value)

        elif isinstance(a, Decimal):
            try:
                numerator, denominator = a.as_integer_ratio()
                self._set_from_bigints(numerator, denominator)
            except OverflowError:
                s = format(a, "f")
                if '.' in s:
                    integer_part, frac_part = s.split('.')
                    numerator_str = integer_part + frac_part
                    denominator_str = '1' + '0' * len(frac_part)
                else:
                    numerator_str = s
                    denominator_str = '1'
                self._set_from_bigints(numerator_str, denominator_str)

        elif isinstance(a, str):
            if '/' in a:
                parts = a.split('/')
                if len(parts) != 2:
                    raise ValueError(f"Invalid fraction string: {a}")
                self._set_from_bigints(parts[0].strip(), parts[1].strip())
            else:
                mpq_set_str(self.value, a.encode(), 10)
                mpq_canonicalize(self.value)

        elif isinstance(a, Fraction):
            bf = <Fraction>a
            mpq_set(self.value, bf.value)

        else:
            raise TypeError(f"Cannot convert type {type(a)} to Fraction")

    # -----------------------
    # Constructor
    # -----------------------
    def __init__(self, *args):
        if len(args) == 0:
            mpq_set_si(self.value, 0, 1)
        elif len(args) == 1:
            self._set_from_object(args[0])
        elif len(args) == 2:
            num, den = args
            if den == 0:
                raise ZeroDivisionError("Denominator cannot be zero")
            self._set_from_bigints(num, den)
        else:
            raise TypeError("Fraction() takes 0, 1, or 2 arguments")

    # -----------------------
    # Arithmetic helpers
    # -----------------------
    cdef Fraction _binary_op(self, object other,
                             void (*op)(mpq_t, mpq_srcptr, mpq_srcptr)):
        cdef Fraction result = Fraction()
        cdef Fraction o
        if isinstance(other, Fraction):
            o = <Fraction>other
        else:
            o = Fraction(other)
        op(result.value, self.value, o.value)
        return result

    def __add__(self, other): return self._binary_op(other, mpq_add)
    def __sub__(self, other): return self._binary_op(other, mpq_sub)
    def __mul__(self, other): return self._binary_op(other, mpq_mul)
    def __truediv__(self, other): return self._binary_op(other, mpq_div)
    # -----------------------
    # Powers and roots
    # -----------------------

    def __pow__(self, object exponent):
        """
        Fraction ** exponent
        Supports integer and Fraction exponents exactly.
        Raises ValueError if root is not exact.
        """
        cdef Fraction result
        cdef mpz_t num_root, den_root
        mpz_init(num_root)
        mpz_init(den_root)

        # Integer exponent
        if isinstance(exponent, int):
            result = Fraction()
            if exponent >= 0:
                mpz_pow_ui(mpq_numref(result.value), mpq_numref(self.value), exponent)
                mpz_pow_ui(mpq_denref(result.value), mpq_denref(self.value), exponent)
            else:
                mpz_pow_ui(num_root, mpq_denref(self.value), -exponent)
                mpz_pow_ui(den_root, mpq_numref(self.value), -exponent)
                mpq_set_num(result.value, num_root)
                mpq_set_den(result.value, den_root)
                mpq_canonicalize(result.value)
            mpz_clear(num_root)
            mpz_clear(den_root)
            return result

        # Fraction exponent (a/b)
        elif isinstance(exponent, Fraction):
            n, d = exponent.numerator, exponent.denominator
            # Step 1: take d-th root exactly
            if not _is_perfect_power(mpq_numref(self.value), d, num_root):
                mpz_clear(num_root)
                mpz_clear(den_root)
                raise ValueError(f"Numerator {self.numerator} is not a perfect {d}-th power")
            if not _is_perfect_power(mpq_denref(self.value), d, den_root):
                mpz_clear(num_root)
                mpz_clear(den_root)
                raise ValueError(f"Denominator {self.denominator} is not a perfect {d}-th power")
            # Step 2: raise result to numerator
            result = Fraction()
            mpz_pow_ui(mpq_numref(result.value), num_root, abs(n))
            mpz_pow_ui(mpq_denref(result.value), den_root, abs(n))
            if n < 0:
                # reciprocal
                mpq_inv(result.value, result.value)
            mpz_clear(num_root)
            mpz_clear(den_root)
            return result

        else:
            raise TypeError(f"Exponent must be int or Fraction, got {type(exponent)}")
    def __rpow__(self, object base):
        """
        Computes base ** self exactly.
        base can be int or Fraction.
        Fractional exponents must yield exact roots.
        """
        cdef Fraction base_frac
        if isinstance(base, int):
            base_frac = Fraction(base)
        elif isinstance(base, Fraction):
            base_frac = <Fraction>base
        else:
            raise TypeError(f"Base must be int or Fraction, got {type(base)}")

        # Negative base check for fractional exponents
        if base_frac.numerator < 0 and isinstance(self, Fraction) and self.denominator != 1:
            raise ValueError("Negative base with fractional exponent cannot be represented as exact Fraction")

        # Delegate to __pow__ logic
        return base_frac.__pow__(self)
    def __radd__(self, other): return Fraction(other) + self
    def __rsub__(self, other): return Fraction(other) - self
    def __rmul__(self, other): return Fraction(other) * self
    def __rtruediv__(self, other): return Fraction(other) / self

    # -----------------------
    # Comparisons
    # -----------------------
    cdef int _cmp(self, object other):
        cdef Fraction o
        if isinstance(other, Fraction):
            o = <Fraction>other
        else:
            o = Fraction(other)
        return mpq_cmp(self.value, o.value)

    def __eq__(self, other): return self._cmp(other) == 0
    def __ne__(self, other): return self._cmp(other) != 0
    def __lt__(self, other): return self._cmp(other) < 0
    def __le__(self, other): return self._cmp(other) <= 0
    def __gt__(self, other): return self._cmp(other) > 0
    def __ge__(self, other): return self._cmp(other) >= 0

    # -----------------------
    # String / representation
    # -----------------------
    def __str__(self):
        cdef char* s = mpq_get_str(NULL, 10, self.value)
        py_str = s.decode()
        free(s)
        return py_str

    def __repr__(self):
        return f"Fraction('{str(self)}')"

    # -----------------------
    # Numerator / Denominator
    # -----------------------
    @property
    def numerator(self):
        cdef mpz_ptr num = mpq_numref(self.value)
        cdef char* s = mpz_get_str(NULL, 10, num)
        py_str = s.decode()
        free(s)
        return int(py_str)

    @property
    def denominator(self):
        cdef mpz_ptr den = mpq_denref(self.value)
        cdef char* s = mpz_get_str(NULL, 10, den)
        py_str = s.decode()
        free(s)
        return int(py_str)

    # -----------------------
    # Conversion
    # -----------------------
    def __float__(self):
        return self.numerator / self.denominator

    def __int__(self):
        return self.numerator // self.denominator

    def __round__(self, n=0):
        return round(float(self), n)

    def __trunc__(self):
        return self.numerator // self.denominator

    def as_integer_ratio(self):
        return (self.numerator, self.denominator)

    def __bool__(self):
        return self.numerator != 0

    def __hash__(self):
        return hash((self.numerator, self.denominator))

    def to_decimal(self, precision=None):
        if precision is None:
            return Decimal(self.numerator) / Decimal(self.denominator)
        else:
            with localcontext() as ctx:
                ctx.prec = precision
                return Decimal(self.numerator) / Decimal(self.denominator)
