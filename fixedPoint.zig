const std = @import("std");

/// high range, low precision fixed point value between -2.000.000 and +2.000.000, with a precision of roughly 0.001
pub const i32p10 = FixedPoint(32, 1024);

/// medium range, medium precision fixed point value between -32000 and +32000, with a precision of roughly 0.000015
pub const i32p16 = FixedPoint(32, 65536);

/// high precision fixed point with i32 integer range and a precision of roughly 0.00000000025
pub const i64p32 = FixedPoint(64, 1 << 32);

/// Configurable fixed point implementation.
/// - `bits` is the total number of bits in the fixed point value
/// - `scaling` is a
pub fn FixedPoint(comptime bits: comptime_int, comptime scaling: comptime_int) type {
    if (scaling < 1)
        @compileError("scaling must be a positive, non-zero integer!");
    const BaseInt = @Type(.{
        .int = .{ .bits = bits, .signedness = .signed },
    });
    const BaseIntUnsigned = @Type(.{
        .int = .{ .bits = bits - 1, .signedness = .unsigned },
    });
    if (scaling > std.math.maxInt(BaseInt))
        @compileError(std.fmt.comptimePrint("scaling must be less than {}", .{std.math.maxInt(BaseInt)}));

    const scaling_bits_max: comptime_int = std.math.log2_int_ceil(BaseIntUnsigned, scaling);

    const significant_digits = std.math.ceil(std.math.log10(@as(comptime_float, scaling)));

    return struct {
        pub const Int = BaseInt;

        pub const Int2 = @Type(.{
            .int = .{ .bits = 2 * bits, .signedness = .signed },
        });

        pub const IntPart = @Type(.{
            .int = .{ .bits = bits - scaling_bits_max, .signedness = .signed },
        });

        pub const precision = 2.0 / @as(comptime_float, scaling);

        // comptime {
        //     @compileLog(bits, scaling, precision, Int, Int2, IntPart);
        // }

        const F = @This();

        raw: Int,

        // conversion operators

        pub fn fromFloat64(v: f64) F {
            // std.debug.print("fromFloat({}, {d})\n", .{ Int, v });
            return .{ .raw = @as(Int, @intFromFloat(scaling * v)) };
        }

        pub fn fromFloat64Nearest(v: f64) F {
            // std.debug.print("fromFloat({}, {d})\n", .{ Int, v });
            return .{ .raw = @as(Int, @intFromFloat((scaling * v) + @as(f64, if (v > 0) 0.5 else -0.5))) };
        }

        pub fn fromFloat(v: f32) F {
            // std.debug.print("fromFloat({}, {d})\n", .{ Int, v });
            return .{ .raw = @as(Int, @intFromFloat(scaling * v)) };
        }

        pub fn fromFloatNearest(v: f32) F {
            // std.debug.print("fromFloat({}, {d})\n", .{ Int, v });
            return .{ .raw = @as(Int, @intFromFloat((scaling * v) + @as(f32, if (v > 0) 0.5 else -0.5))) };
        }

        pub fn toFloat(v: F, comptime T: type) T {
            // std.debug.print("toFloat({}, {})\n", .{ Int, v.raw });
            _ = @typeInfo(T).float;
            return @as(T, @floatFromInt(v.raw)) / scaling;
        }

        pub fn fromInt(i: IntPart) F {
            return .{ .raw = scaling * @as(Int, i) };
        }

        pub fn toInt(f: F) IntPart {
            // std.debug.print("toInt({}, {})\n", .{ Int, f.raw });
            return @as(IntPart, @intCast(@divTrunc(f.raw, scaling)));
        }

        // arithmetic operators:

        pub fn add(a: F, b: F) F {
            return .{ .raw = a.raw + b.raw };
        }

        pub fn sub(a: F, b: F) F {
            return .{ .raw = a.raw - b.raw };
        }

        pub fn neg(a: F) F {
            return .{ .raw = -a.raw };
        }

        pub fn mul(a: F, b: F) F {
            return .{ .raw = @as(Int, @intCast(scaleDown(@as(Int2, a.raw) * @as(Int2, b.raw)))) };
        }

        pub fn div(a: F, b: F) F {
            return .{ .raw = @as(Int, @intCast(@divTrunc(scaleUp(a.raw), b.raw))) };
        }

        pub fn mod(a: F, b: F) F {
            return .{ .raw = @mod(a.raw, b.raw) };
        }

        pub fn sqrt(a: F) F {
            // Handle non-positive inputs
            if (a.raw <= 0) return .{ .raw = 0 };

            // Initial guess: using 'a' may be too large, so let's try a simple shift
            var x = F{ .raw = @max(1, a.raw >> (precision_bits / 2)) };

            // Number of iterations; 5 should suffice for convergence with good guess
            const max_iterations = 22;

            // Compile-time loop for efficiency
            inline for (0..max_iterations) |_| {
                // Compute a/x with fixed-point division
                const quotient = div(a, x);

                // Use Int2 to avoid overflow and preserve precision
                const sum = @as(Int2, x.raw) + @as(Int2, quotient.raw);

                // Divide by 2 and cast back to Int
                const x_new_raw = @as(Int, @intCast(sum >> 1));

                // Check for convergence
                if (x_new_raw == x.raw) break;
                x.raw = x_new_raw;
            }

            return x;
        }

        // relational operators:

        pub fn lessThan(a: F, b: F) bool {
            return a.raw < b.raw;
        }

        pub fn greaterThan(a: F, b: F) bool {
            return a.raw > b.raw;
        }

        pub fn lessOrEqual(a: F, b: F) bool {
            return a.raw <= b.raw;
        }

        pub fn greaterOrEqual(a: F, b: F) bool {
            return a.raw >= b.raw;
        }

        pub fn eql(a: F, b: F) bool {
            return a.raw == b.raw;
        }

        pub const tau = fromFloat64Nearest(std.math.tau);
        pub const pi = fromFloat64Nearest(std.math.pi);

        // Reduce angle to [-π, π)
        fn reduceAngle(theta: F) F {
            var reduced = theta;

            // Perform modulo 2π
            if (reduced.raw >= 0) {
                reduced.raw = @mod(reduced.raw, tau.raw);
            } else {
                reduced.raw = -@mod(-reduced.raw, tau.raw);
            }

            // Adjust result into the symmetric interval [-π, π)
            if (reduced.raw >= pi.raw) { // Note: using >= to exclude π
                reduced.raw -= tau.raw;
            } else if (reduced.raw < -pi.raw) {
                reduced.raw += tau.raw;
            }

            return reduced;
        }

        // Sine function using Taylor series with fixed-point arithmetic.
        fn sinFixedPoint(theta: F, max_terms: usize) F {
            const x = reduceAngle(theta); // Reduce angle to [-π, π]
            var sin_x = x; // Initial term T₀ = x
            const x2 = mul(x, x); // Pre-compute x²
            var term = x; // T₀ also (this will be updated)

            // Taylor series: sin(x) = x - x^3/3! + x^5/5! - x^7/7! + ...
            for (1..max_terms + 1) |n| {
                // Update term using recurrence: Tₙ = Tₙ₋₁ * (-x²) / ((2n-1) * (2n))
                term = div(mul(term, neg(x2)), fromInt(@intCast((2 * n) * (2 * n + 1))));
                if (term.raw == fromInt(0).raw) break;
                sin_x = add(sin_x, term);
            }
            return sin_x;
        }

        // Sine function
        pub fn sin(theta: F) F {
            return sinFixedPoint(theta, 8);
        }

        // Cosine function using Taylor series with fixed-point arithmetic.
        fn cosFixedPoint(theta: F, max_terms: usize) F {
            const x = reduceAngle(theta); // Reduce angle to [-π, π]
            //std.debug.print("Angle:{}\n", .{x});
            var cos_x = fromInt(1); // Initial term T₀ = 1
            var term = fromInt(1); // T₀ also (this will be updated)
            const x2 = mul(x, x); // Pre-compute x²

            // Compute Taylor series: cos(x) = 1 - x²/2! + x⁴/4! - x⁶/6! + ...
            for (1..max_terms + 1) |n| {
                // Update term using recurrence: Tₙ = Tₙ₋₁ * (-x²) / ((2n-1) * (2n))
                term = div(mul(term, neg(x2)), fromInt(@intCast(((2 * n - 1) * (2 * n)))));
                if (term.raw == fromInt(0).raw) break;
                cos_x = add(cos_x, term);
            }
            return cos_x;
        }

        // Cosine function
        pub fn cos(theta: F) F {
            return cosFixedPoint(theta, 8);
        }

        // format api

        pub fn format(f: F, comptime fmt: []const u8, options: std.fmt.FormatOptions, writer: anytype) !void {
            if (comptime std.mem.eql(u8, fmt, "x")) {
                return std.fmt.formatFloatHexadecimal(f.toFloat(f32), options, writer);
            }

            const mode: std.fmt.format_float.Format = comptime blk: {
                if (fmt.len == 0 or std.mem.eql(u8, fmt, "any") or std.mem.eql(u8, fmt, "d")) {
                    break :blk .decimal;
                } else if (std.mem.eql(u8, fmt, "e")) {
                    break :blk .scientific;
                } else {
                    @compileError(std.fmt.comptimePrint("Invalid fmt for FixedPoint({},{}): {{{s}}}", .{ bits, scaling, fmt }));
                }
            };

            const foptions = std.fmt.format_float.FormatOptions{
                .mode = mode,
                .precision = if (mode == .decimal) blk: {
                    break :blk if (options.precision) |p| p else significant_digits;
                } else options.precision,
            };

            var buf: [std.fmt.format_float.bufferSize(mode, f32)]u8 = undefined;
            const s = try std.fmt.format_float.formatFloat(&buf, f.toFloat(f32), foptions);
            try std.fmt.formatBuf(s, options, writer);
        }

        // implement guaranteed shift semantics for POT scalings

        const is_pot = std.math.isPowerOfTwo(scaling);
        const precision_bits: comptime_int = if (is_pot)
            std.math.log2_int(u64, scaling)
        else
            @compileError("scaling is not a power a power of two.");
        fn scaleUp(in: Int) Int2 {
            return if (is_pot)
                return @as(Int2, in) << precision_bits
            else
                return @as(Int2, in) * scaling;
        }

        fn scaleDown(in: Int2) Int {
            return if (is_pot)
                @as(Int, @intCast(in >> precision_bits))
            else
                @as(Int, @intCast(@divTrunc(in, scaling)));
        }
    };
}

test {
    _ = TestSuite(i32p10);
    _ = TestSuite(i32p16);
    _ = TestSuite(i64p32);
    _ = TestSuite(FixedPoint(32, 1000));
    _ = TestSuite(FixedPoint(64, 1000));
}

fn TestSuite(comptime FP: type) type {
    return struct {
        const float_test_vals = [_]f32{
            0.0,
            1.0,
            std.math.e,
            std.math.pi,
            42.0,
            42.1337,
            21.0,
            13.37,
            100.0,
            -100.0,
            // limit is roughly 100, as we're doing a val*val*2.5 and will get a overflow otherwise for a i32p16
        };
        test "float conversion" {
            for (float_test_vals) |val| {
                const fp = FP.fromFloat(val);
                const f = fp.toFloat(f32);
                try std.testing.expectApproxEqAbs(val, f, FP.precision);
            }
        }

        test "int conversion" {
            const test_vals = [_]FP.IntPart{
                0, 1,  2,  3,  4,  5,  6,  7,  2000,  20_000,  30_000,  std.math.maxInt(i16), std.math.maxInt(FP.IntPart),
                0, -1, -2, -3, -4, -5, -6, -7, -2000, -20_000, -30_000, std.math.minInt(i16), std.math.minInt(FP.IntPart),
            };

            for (test_vals) |val| {
                const fp = FP.fromInt(val);
                const f = fp.toInt();
                try std.testing.expectEqual(val, f);
            }
        }

        test "add arithmetic" {
            for (float_test_vals) |val| {
                const fp = FP.fromFloat(val);
                const f = fp.add(fp).add(FP.fromFloat(10)).toFloat(f32);
                try std.testing.expectApproxEqAbs(2.0 * val + 10, f, FP.precision);
            }
        }

        test "sub arithmetic" {
            for (float_test_vals) |val| {
                const fp = FP.fromFloat(val);
                const f = fp.sub(FP.fromFloat(10)).toFloat(f32);
                try std.testing.expectApproxEqAbs(val - 10, f, FP.precision);
            }
        }

        test "mul arithmetic" {
            for (float_test_vals) |val| {
                const fp = FP.fromFloat(val);
                const f = fp.mul(fp).mul(FP.fromFloat(2.5)).toFloat(f32);
                try std.testing.expectApproxEqRel(val * val * 2.5, f, @sqrt(FP.precision));
            }
        }

        test "div arithmetic" {
            for (float_test_vals) |val| {
                const fp = FP.fromFloat(val);
                const f = fp.div(FP.fromFloat(2.5)).toFloat(f32);
                try std.testing.expectApproxEqRel(val / 2.5, f, @sqrt(FP.precision));
            }
        }

        test "mod arithmetic" {
            const test_vals = [_]f32{
                0.0,
                1.0,
                std.math.e,
                std.math.pi,
            };

            for (test_vals) |val| {
                const fp = FP.fromFloat(val);
                const f = fp.mod(FP.fromFloat(2.5)).toFloat(f32);
                try std.testing.expectApproxEqRel(@mod(val, 2.5), f, @sqrt(FP.precision));
            }
        }
    };
}

test "basic printing" {
    const Meter = FixedPoint(32, 1000); // millimeter precision meter units, using 32 bits

    const position_1 = Meter.fromFloat(10); // 10m
    const position_2 = position_1.add(Meter.fromFloat(0.01)); // add 1cm
    const position_3 = position_2.add(Meter.fromFloat(0.09)); // add 9cm
    const distance = position_3.sub(position_1);

    var buffer: [64]u8 = undefined;

    try std.testing.expectEqualStrings("Distance = 0.100m", try std.fmt.bufPrint(&buffer, "Distance = {}m", .{distance}));
    try std.testing.expectEqualStrings("Distance = 0.10m", try std.fmt.bufPrint(&buffer, "Distance = {d:0.2}m", .{distance}));
    try std.testing.expectEqualStrings("Distance = 0.1m", try std.fmt.bufPrint(&buffer, "Distance = {d:0.1}m", .{distance}));
}

test "cos" {
    const Time = FixedPoint(32, 1000);

    const position_1 = Time.fromFloat(0);
    const position_2 = Time.tau;
    const position_3 = Time.pi;
    const position_4 = Time.div(Time.pi, .fromInt(2));
    const position_5 = Time.div(Time.pi, .fromInt(4));

    std.debug.print("{}\n", .{Time.cos(position_1)});
    std.debug.print("{}\n", .{Time.cos(position_2)});
    std.debug.print("{}\n", .{Time.cos(position_3)});
    std.debug.print("{}\n", .{Time.cos(position_4)});
    std.debug.print("{}\n", .{Time.cos(position_5)});
}

test "sin" {
    const Time = FixedPoint(32, 1000);

    const position_1 = Time.fromFloat(0);
    const position_2 = Time.tau;
    const position_3 = Time.pi;
    const position_4 = Time.div(Time.pi, .fromInt(2));
    const position_5 = Time.div(Time.pi, .fromInt(4));

    std.debug.print("{}\n", .{Time.sin(position_1)});
    std.debug.print("{}\n", .{Time.sin(position_2)});
    std.debug.print("{}\n", .{Time.sin(position_3)});
    std.debug.print("{}\n", .{Time.sin(position_4)});
    std.debug.print("{}\n", .{Time.sin(position_5)});
}
