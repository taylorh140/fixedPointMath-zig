const std = @import("std");
const Complex = @import("complexFixedPoint.zig").Complex;
const fp = @import("fixedPoint.zig");
const lut = @import("trigLut.zig");

pub fn GetFFT(fftBits: comptime_int, FixedType_t: type, TrigLut: ?type) type {
    return struct {
        const fftSize_t = @Type(.{ .int = .{ .bits = fftBits, .signedness = .unsigned } });

        const F = FixedType_t;
        const N: usize = std.math.maxInt(fftSize_t) + 1; // Fixed size for FFT/IFFT

        // FFT implementation for fixed size N=256
        pub fn fft(data: *[N]Complex(F)) void {
            // Bit-reversal permutation using precomputed table
            for (0..N) |i| {
                const j = @bitReverse(@as(fftSize_t, @intCast(i)));
                if (j > i) {
                    std.mem.swap(Complex(F), &data[i], &data[j]);
                }
            }

            // Butterfly operations
            var M: usize = 1;
            while (M < N) : (M <<= 1) {
                const M2 = M << 1;
                for (0..M) |j| {
                    const theta = F.div(F.mul(F.tau.neg(), F.fromInt(@intCast(j))), F.fromInt(@intCast(M2)));
                    // std.debug.print("THETA: {} J:{} M2:{}\n", .{ theta, j, M2 });
                    const twiddle = if (TrigLut) |LUT| blk: {
                        const sinCos = LUT.cosSin(theta);
                        break :blk Complex(F).init(sinCos.cos, sinCos.sin);
                    } else Complex(F).init(F.cos(theta), F.sin(theta));

                    for (0..N / M2) |k| {
                        const idx1 = k * M2 + j;
                        const idx2 = idx1 + M;
                        const t = data[idx2].mul(twiddle);
                        data[idx2] = data[idx1].sub(t);
                        data[idx1] = data[idx1].add(t);
                    }
                }
            }
        }

        // IFFT implementation for fixed size N=256
        pub fn ifft(data: *[N]Complex(F)) void {
            // Conjugate the data
            for (data) |*x| {
                x.im = x.im.neg();
            }

            // Perform FFT
            fft(data);

            // Conjugate again and scale by 1/N
            const scale = F.div(F.fromFloat(1.0), F.fromInt(N));
            for (data) |*x| {
                x.re = F.mul(x.re, scale);
                x.im = F.mul(x.im, scale.neg()); // Conjugate and scale
            }
        }
    };
}

// Example usage
pub fn main() !void {
    const ftype = fp.FixedPoint(32, std.math.pow(usize, 2, 20));
    const cftype = Complex(ftype);
    const LUT = lut.TrigTable(ftype, 256, .linear);
    const fft_t = GetFFT(5, ftype, LUT);

    var data: [fft_t.N]cftype = undefined;

    // Initialize with some data (e.g., real input, imaginary = 0)
    for (&data, 0..) |*x, i| {
        const t = ftype.fromInt(@intCast(i)).div(.fromInt(32));
        const num = t.mul(ftype.fromFloat(std.math.tau));
        const inner = ftype.mul(.fromFloat(1), num);
        std.debug.print("{any}:{any}:{any}\n", .{ inner, ftype.sin(inner), ftype.fromFloat(std.math.sin(inner.toFloat(f32))) });
        x.* = cftype.init(ftype.sin(inner), ftype.fromFloat(0.0));
    }

    std.debug.print("Time Domain:\n", .{});
    for (data) |pt|
        std.debug.print("{any}\n", .{pt.re});

    // Perform FFT
    fft_t.fft(&data);

    std.debug.print("Frequency Domain:\n", .{});
    for (data) |pt|
        std.debug.print("{any}\n", .{pt.magnitude()});

    // Perform IFFT
    fft_t.ifft(&data);

    std.debug.print("Time Domain:\n", .{});
    for (data) |pt|
        std.debug.print("{any}\n", .{pt.re});

    // Results could be printed here for verification
}
