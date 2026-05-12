using System;
using System.Numerics;
using System.Security.Cryptography;
using System.Text;

class Program {
    static void Main(string[] args) {
        var username = "a";
        var hex      = "ED14B132D4BBD428";

        // Curve parameters (custom short Weierstrass over prime field)
        BigInteger p     = H("FE4382C5413A02FF");
        BigInteger a     = H("5BA3091245C856AB");
        BigInteger b     = H("C2AB76EF7FE1D7F2");
        BigInteger order = 18_321_631_499_947_426_219UL;
        BigInteger priv  = 2_826_820_123_527_714_983UL;

        // Recover user public key from token (X coord -> Y via curve equation)
        var ux = H(hex);
        var uy = SqrtMod((BigInteger.ModPow(ux, 3, p) + a * ux + b) % p, p);
        if (uy.IsEven) uy = p - uy;   // tinyec parity convention: odd Y preferred

        Console.WriteLine($"Username: {username}");
        Console.WriteLine($"Token (hex): {hex}");
        Console.WriteLine($"Ux: {ux:x}");
        Console.WriteLine($"Uy: {uy:x}");

        // ECDH: shared secret = priv * userPublicPoint
        var shared = Mul(priv, (ux, uy));

        // Shared secret X as even-length hex string
        var sk = shared.x.ToString("x");
        if (sk.Length % 2 != 0) sk = "0" + sk;

        Console.WriteLine($"Shared Key X: {sk}");

        // Double HMAC-SHA1 password derivation
        string password = CalcPassword(sk, username);
        Console.WriteLine($"Password: {password}");
    }

    static string CalcPassword(string sk, string user) {
        var mid = HMACSHA1.HashData(new byte[20], Convert.FromHexString(sk));
        using var h = new HMACSHA1(mid);
        h.TransformBlock(Encoding.UTF8.GetBytes(user), 0, user.Length, null, 0);
        h.TransformFinalBlock(new byte[] { 0x01 }, 0, 1);
        return Convert.ToHexString(h.Hash!, 0, 8).ToLowerInvariant();
    }

    // Parse hex string to positive BigInteger (leading zero avoids sign-bit issues)
    static BigInteger H(string s) =>
        BigInteger.Parse("0" + s, System.Globalization.NumberStyles.HexNumber);

    // Modular inverse via extended Euclidean algorithm
    static BigInteger Inv(BigInteger v) {
        BigInteger p     = H("FE4382C5413A02FF");
        var (g, x, y, t) = (p, BigInteger.Zero, BigInteger.One, v);
        while (t != 0) { var q = g / t; (g, t) = (t, g - q * t); (x, y) = (y, x - q * y); }
        return ((x % p) + p) % p;
    }

    // Affine point addition; null represents the point at infinity
    static (BigInteger x, BigInteger y)? Add((BigInteger x, BigInteger y)? p1, (BigInteger x, BigInteger y)? p2) {
        BigInteger p     = H("FE4382C5413A02FF");
        BigInteger a     = H("5BA3091245C856AB");
        
        if (p1 is null) return p2;
        if (p2 is null) return p1;
        var (x1, y1) = p1.Value;
        var (x2, y2) = p2.Value;
        BigInteger lam;
        if (x1 == x2) {
            if (y1 != y2) return null;
            lam = (3 * x1 * x1 + a) % p * Inv(2 * y1 % p) % p;  // point doubling
        } else {
            lam = (y2 - y1 + p) % p * Inv((x2 - x1 + p) % p) % p;
        }
        var xr = ((lam * lam - x1 - x2) % p + p) % p;
        var yr = ((lam * (x1 - xr) - y1) % p + p) % p;
        return (xr, yr);
    }

    // Double-and-add scalar multiplication
    static (BigInteger x, BigInteger y) Mul(BigInteger k, (BigInteger x, BigInteger y) pt) {
        BigInteger order = 18_321_631_499_947_426_219UL;
        (BigInteger x, BigInteger y)? result = null;
        var addend = ((BigInteger x, BigInteger y)?)pt;
        k = ((k % order) + order) % order;
        while (k > 0) {
            if (!k.IsEven) result = Add(result, addend);
            addend = Add(addend, addend);
            k >>= 1;
        }
        return result!.Value;
    }

    // Tonelli-Shanks modular square root (r^2 == n mod p)
    static BigInteger SqrtMod(BigInteger n, BigInteger mod) {
        n = ((n % mod) + mod) % mod;
        if (n == 0) return 0;
        if (mod % 4 == 3) return BigInteger.ModPow(n, (mod + 1) / 4, mod);
        var (q, s) = (mod - 1, 0);
        while (q.IsEven) { q >>= 1; s++; }
        var z = new BigInteger(2);
        while (BigInteger.ModPow(z, (mod - 1) / 2, mod) != mod - 1) z++;
        var (m, c, t, r) = ((BigInteger)s, BigInteger.ModPow(z, q, mod),
                            BigInteger.ModPow(n, q, mod), BigInteger.ModPow(n, (q + 1) / 2, mod));
        while (true) {
            if (t == 1) return r;
            if (t == 0) return 0;
            var (tmp, i) = (t, 0);
            while (tmp != 1) { tmp = tmp * tmp % mod; i++; }
            var b = BigInteger.ModPow(c, BigInteger.Pow(2, (int)(m - i - 1)), mod);
            (m, c, t, r) = (i, b * b % mod, t * (b * b % mod) % mod, r * b % mod);
        }
    }
}
