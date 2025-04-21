const std = @import("std");
const fp = @import("fixedPoint.zig");

/// Generates a trigonometric lookup table type for sine and cosine with the specified fixed-point parameters.
/// - total_bits: Total number of bits for the fixed-point type (e.g., 32).
/// - fractional_bits: Number of fractional bits (e.g., derived from std.math.pow(usize, 2, 11) = 2048, implying 11 fractional bits).
/// - M: Number of points in the lookup table, must be a power of 2 (e.g., 256).
/// - interpolation: Interpolation method, either .nearest or .linear.
pub fn TrigTable(F: type, comptime M: u32, comptime interpolation: enum { nearest, linear }) type {
    const total_bits = @typeInfo(F.Int).int.bits;
    const int_bits = @typeInfo(F.IntPart).int.bits;
    const fractional_bits = total_bits - int_bits;

    // Ensure M is a power of 2
    if ((M & (M - 1)) != 0) @compileError("M must be a power of 2");
    const k = @ctz(M); // log2(M)
    if (total_bits < k) @compileError("total_bits must be at least log2(M)");

    // Generate sine table at comptime
    const sin_table = blk: {
        @setEvalBranchQuota(M * 10);
        var table: [M]F = undefined;
        for (0..M) |i| {
            const theta = std.math.tau * @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(M));
            const sin_val = std.math.sin(theta);
            table[i] = .fromFloat64Nearest(sin_val);
            //@compileLog(table[i].toFloat(f32));
        }

        break :blk table;
    };

    return struct {
        /// Number of points in the lookup table
        pub const num_entries = M;
        /// Computes the sine of the given fixed-point angle
        pub fn sin(angle: F) F {
            if (interpolation == .nearest) {
                const scale: usize = @intCast(F.div(angle, F.tau).raw);
                const preRound = (scale >> (fractional_bits - k - 1)); // Leave one bit for rounding
                const rounded = ((preRound) + 1) >> 1;
                const index = (rounded) & (M - 1);
                return sin_table[index];
            } else { // linear
                const angleV: usize = @intCast(F.div(angle, F.tau).raw);
                const index0 = (angleV >> (fractional_bits - k)) & (M - 1);
                const index1 = (index0 + 1) & (M - 1);
                const mask = (1 << (fractional_bits - k)) - 1;
                const remainder = angleV & mask;
                const frac = if (fractional_bits - k >= fractional_bits)
                    remainder >> (fractional_bits - k - fractional_bits)
                else
                    remainder << (fractional_bits - (fractional_bits - k));
                const a = sin_table[(index0)];
                const b = sin_table[(index1)];
                const diff = b.sub(a);

                const interpolated = a.add(diff.mul(F{ .raw = @intCast(frac) }));
                return interpolated;
            }
        }

        fn cos_idx(idx: usize) usize {
            const offset = M + M / 4;
            return (idx + offset) & (M - 1);
        }

        /// Computes the cosine of the given fixed-point angle
        pub fn cos(angle: F) F {
            if (interpolation == .nearest) {
                const scale: usize = @intCast(F.div(angle, F.tau).raw);
                const preRound = (scale >> (fractional_bits - k - 1)); // Leave one bit for rounding
                const rounded = ((preRound) + 1) >> 1;
                const index = (rounded) & (M - 1);
                return sin_table[cos_idx(index)];
            } else { // linear
                const angleV: usize = @intCast(F.div(angle, F.tau).raw);
                const index0 = (angleV >> (fractional_bits - k)) & (M - 1);
                const index1 = (index0 + 1) & (M - 1);

                const mask = (1 << (fractional_bits - k)) - 1;
                const remainder = angleV & mask;
                const frac = if (fractional_bits - k >= fractional_bits)
                    remainder >> (fractional_bits - k - fractional_bits)
                else
                    remainder << (fractional_bits - (fractional_bits - k));
                const a = sin_table[cos_idx(index0)];
                const b = sin_table[cos_idx(index1)];
                const diff = b.sub(a);
                //std.debug.print(" >>> angle:{b} frac:{} a:{} b:{}  idx0:{} idx1:{} mask{b}\n", .{ angleV, frac, a, b, index0, index1, mask });
                const interpolated = a.add(diff.mul(F{ .raw = @intCast(frac) }));
                return interpolated;
            }
        }

        /// Computes both cosine and sine of the given fixed-point angle
        pub fn cosSin(angle: F) struct { cos: F, sin: F } {
            if (interpolation == .nearest) {
                const scale: usize = @intCast(F.div(angle, F.tau).raw);
                const preRound = (scale >> (fractional_bits - k - 1)); // Leave one bit for rounding
                const rounded = ((preRound) + 1) >> 1;
                const index = (rounded) & (M - 1);
                return .{ sin_table[cos_idx(index)], sin_table[index] };
            } else { // linear
                const angleV: usize = @intCast(F.div(angle, F.tau).raw & ((1 << fractional_bits) - 1));
                const index0 = (angleV >> (fractional_bits - k)) & (M - 1);
                const index1 = (index0 + 1) & (M - 1);

                const mask = (1 << (fractional_bits - k)) - 1;
                const remainder = angleV & mask;
                const frac = if (fractional_bits - k >= fractional_bits)
                    remainder >> (fractional_bits - k - fractional_bits)
                else
                    remainder << (fractional_bits - (fractional_bits - k));
                const cos_val = blk: {
                    const a = sin_table[cos_idx(index0)];
                    const b = sin_table[cos_idx(index1)];
                    const diff = b.sub(a);
                    //std.debug.print(" >>> angle:{b} frac:{} a:{} b:{}  idx0:{} idx1:{} mask{b}\n", .{ angleV, frac, a, b, index0, index1, mask });
                    const interpolated = a.add(diff.mul(F{ .raw = @intCast(frac) }));
                    break :blk interpolated;
                };
                const sin_val = blk: {
                    const a = sin_table[index0];
                    const b = sin_table[index1];
                    const diff = b.sub(a);
                    //std.debug.print(" >>> angle:{b} frac:{} a:{} b:{}  idx0:{} idx1:{} mask{b}\n", .{ angleV, frac, a, b, index0, index1, mask });
                    const interpolated = a.add(diff.mul(F{ .raw = @intCast(frac) }));
                    break :blk interpolated;
                };
                return .{ .cos = cos_val, .sin = sin_val };
            }
        }
    };
}

// Example usage
pub fn main() void {
    const F = fp.i32p16;
    // Define the LUT with 32 total bits, 11 fractional bits (2^11 = 2048), 256 points, and nearest interpolation
    const LUT = TrigTable(F, 256, .linear);
    const extra = 4;

    for (0..LUT.num_entries * extra) |i| {
        const t: f64 = @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(LUT.num_entries * extra));
        const num = t * std.math.tau;
        const frequency = 1; //Hz
        const inner = frequency * num;
        const val = LUT.cosSin(.fromFloat64Nearest(inner));
        // _ = val;
        std.io.getStdOut().writer().print("{any}:{any}:{any}:{any}:{any}:{any}\n", .{ inner, F.fromFloat64Nearest(inner).raw, val.cos.raw, val.sin.raw, F.fromFloat64Nearest(std.math.cos(inner)).raw, F.fromFloat64Nearest(std.math.sin(inner)).raw }) catch {};
    }

    // std.debug.print("{}\n", .{@as(u32, @intFromFloat(0.49999999))});
    // std.debug.print("{}\n", .{@as(u32, @intFromFloat(0.5001))});

    // std.debug.print("{}\n", .{@as(i32, @intFromFloat(0.49999999 + 0.5))});
    // std.debug.print("{}\n", .{@as(i32, @intFromFloat(0.5001 + 0.5))});

    // std.debug.print("{}\n", .{@as(i32, @intFromFloat(-0.49999999 - 0.5))});
    // std.debug.print("{}\n", .{@as(i32, @intFromFloat(-0.5001 - 0.5))});
}
