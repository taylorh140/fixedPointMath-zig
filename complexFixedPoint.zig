const std = @import("std");
const fp = @import("fixedPoint.zig");

/// A complex number consisting of a real an imaginary part. T must be a fixed-point value.
pub fn Complex(comptime T: type) type {
    return struct {
        const Self = @This();

        /// Real part.
        re: T,

        /// Imaginary part.
        im: T,

        /// Create a new Complex number from the given real and imaginary parts.
        pub fn init(re: T, im: T) Self {
            return Self{
                .re = re,
                .im = im,
            };
        }

        /// Returns the sum of two complex numbers.
        pub fn add(self: Self, other: Self) Self {
            return Self{
                .re = T.add(self.re, other.re),
                .im = T.add(self.im, other.im),
            };
        }

        /// Returns the subtraction of two complex numbers.
        pub fn sub(self: Self, other: Self) Self {
            return Self{
                .re = T.sub(self.re, other.re),
                .im = T.sub(self.im, other.im),
            };
        }

        /// Returns the product of two complex numbers.
        pub fn mul(self: Self, other: Self) Self {
            return Self{
                .re = T.sub(T.mul(self.re, other.re), T.mul(self.im, other.im)),
                .im = T.add(T.mul(self.im, other.re), T.mul(self.re, other.im)),
            };
        }

        /// Returns the quotient of two complex numbers.
        pub fn div(self: Self, other: Self) Self {
            const re_num = T.add(T.mul(self.re, other.re), T.mul(self.im, other.im));
            const im_num = T.sub(T.mul(self.im, other.re), T.mul(self.re, other.im));
            const den = T.add(T.mul(other.re, other.re), T.mul(other.im, other.im));

            return Self{
                .re = T.div(re_num, den),
                .im = T.div(im_num, den),
            };
        }

        /// Returns the complex conjugate of a number.
        pub fn conjugate(self: Self) Self {
            return Self{
                .re = self.re,
                .im = T.neg(self.im),
            };
        }

        /// Returns the negation of a complex number.
        pub fn neg(self: Self) Self {
            return Self{
                .re = T.neg(self.re),
                .im = T.neg(self.im),
            };
        }

        /// Returns the product of complex number and i=sqrt(-1)
        pub fn mulbyi(self: Self) Self {
            return Self{
                .re = T.neg(self.im),
                .im = self.re,
            };
        }

        /// Returns the reciprocal of a complex number.
        pub fn reciprocal(self: Self) Self {
            const m = T.add(T.mul(self.re, self.re), T.mul(self.im, self.im));
            return Self{
                .re = T.div(self.re, m),
                .im = T.neg(T.div(self.im, m)),
            };
        }

        /// Returns the magnitude of a complex number.
        pub fn magnitude(self: Self) T {
            return T.sqrt(T.add(T.mul(self.re, self.re), T.mul(self.im, self.im)));
        }

        pub fn squaredMagnitude(self: Self) T {
            return T.add(T.mul(self.re, self.re), T.mul(self.im, self.im));
        }
    };
}

test "multiply" {
    const e_t = fp.i32p10;
    const c_t = Complex(e_t);
    const a = c_t.init(.fromFloat(-1.0), .fromFloat(2.0));
    const b = c_t.init(.fromFloat(-3.0), .fromFloat(7.0));

    const p = a.mul(b);

    //std.debug.print("{}\n", .{p});
    std.debug.assert(p.re.raw == e_t.fromFloat(-11).raw);
    std.debug.assert(p.im.raw == e_t.fromFloat(-13).raw);
}

test "divide" {
    const e_t = fp.i32p10;
    const c_t = Complex(e_t);
    const a = c_t.init(.fromFloat(-18.0), .fromFloat(42.0));
    const b = c_t.init(.fromFloat(-3.0), .fromFloat(7.0));

    const p = a.div(b);

    //std.debug.print("{}\n", .{p});
    std.debug.assert(p.re.raw == e_t.fromFloat(6.0).raw);
    std.debug.assert(p.im.raw == e_t.fromFloat(0.0).raw);
}

test "add" {
    const e_t = fp.i32p10;
    const c_t = Complex(e_t);
    const a = c_t.init(.fromFloat(-18.0), .fromFloat(42.0));
    const b = c_t.init(.fromFloat(-3.0), .fromFloat(7.0));

    const p = a.add(b);

    //std.debug.print("{}\n", .{p});
    std.debug.assert(p.re.raw == e_t.fromFloat(-21.0).raw);
    std.debug.assert(p.im.raw == e_t.fromFloat(49.0).raw);
}

test "sub" {
    const e_t = fp.i32p10;
    const c_t = Complex(e_t);
    const a = c_t.init(.fromFloat(-18.0), .fromFloat(42.0));
    const b = c_t.init(.fromFloat(-3.0), .fromFloat(7.0));

    const p = a.sub(b);

    //std.debug.print("{}\n", .{p});
    std.debug.assert(p.re.raw == e_t.fromFloat(-15.0).raw);
    std.debug.assert(p.im.raw == e_t.fromFloat(35.0).raw);
}

test "conj" {
    const e_t = fp.i32p10;
    const c_t = Complex(e_t);
    const a = c_t.init(.fromFloat(-18.0), .fromFloat(42.0));

    const p = a.conjugate();

    //std.debug.print("{}\n", .{p});
    std.debug.assert(p.re.raw == e_t.fromFloat(-18.0).raw);
    std.debug.assert(p.im.raw == e_t.fromFloat(-42.0).raw);
}

test "sqmag" {
    const e_t = fp.i32p10;
    const c_t = Complex(e_t);
    const a = c_t.init(.fromFloat(-18.0), .fromFloat(42.0));

    const p = a.squaredMagnitude();

    std.debug.print("{}\n", .{p});
    std.debug.assert(p.raw == e_t.fromFloat(2088.0).raw);
}

test "mag" {
    const e_t = fp.i32p16;
    const c_t = Complex(e_t);
    const a = c_t.init(.fromFloat(-18.0), .fromFloat(42.0));

    const p = a.magnitude();

    std.debug.print("{}\n", .{p});
    std.debug.assert(p.raw == e_t.fromFloat(45.69463863518345).raw);
}
