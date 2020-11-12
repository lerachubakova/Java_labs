import java.io.FileWriter;
import java.io.IOException;
import java.math.*;
import java.util.Scanner;

public class Main {

    public static BigInteger MODPOW(BigInteger a, BigInteger e, BigInteger n) {
        BigInteger accum = a;
        if (e.equals(BigInteger.ZERO)) {
            return BigInteger.ONE;
        }
        ;
        int bitptr = e.bitLength() - 1;
        for (bitptr--; bitptr >= 0; bitptr--) {
            accum = accum.multiply(accum).mod(n);
            if (e.testBit(bitptr)) {
                accum = accum.multiply(a).mod(n);
            }
            ;
        }
        ;
        return accum;
    }

    public static BigInteger GCD(BigInteger a, BigInteger b) {
        if (a.equals(BigInteger.ZERO))
            return b;
        if (b.equals(BigInteger.ZERO))
            return a;
        if (a.equals(b))
            return a;
        if (a.compareTo(b) > 0)
            return GCD(a.subtract(b), b);
        return GCD(a, b.subtract(a));
    }

    public static void main(String[] args) throws IOException {

        FileWriter out = new FileWriter("Report.txt");

        BigInteger p = new BigInteger("563036103490583", 10);
        BigInteger q = new BigInteger("1063300642915937", 10);
        BigInteger e = new BigInteger("372585779765210097553647509959", 10);
        BigInteger X1 = new BigInteger("399754188907643924420059310699", 10);
        BigInteger X2;
        BigInteger n = p.multiply(q);
        BigInteger fi = p.subtract(BigInteger.ONE).multiply(q.subtract(BigInteger.ONE));
        BigInteger Y1;
        BigInteger Y2 = new BigInteger("293314580135454643114146935352", 10);
        BigInteger d = e.modInverse(fi);

        if (!GCD(e, fi).equals(BigInteger.ONE)) {
            System.out.println("Wrong e...");
            System.exit(-1);
        }

        System.out.println("Generate...");
        out.write("d: " + d.toString(10) + "\n");

        System.out.println("Encrypt...");
        Y1 = MODPOW(X1, e, n);
        out.write("Y1: " + Y1.toString(10) + "\n");
        out.write(MODPOW(Y1, d, n).equals(X1) ? "Encryption OK\n" : "Encryption failed\n");

        System.out.println("Decrypt...");
        X2 = MODPOW(Y2, d, n);
        out.write("X2: " + X2.toString(10) + "\n");
        out.write(MODPOW(X2, e, n).equals(Y2) ? "Decryption OK\n" : "Decryption failed\n");
        
        System.out.println("Results in Report.txt");
        out.close();
    }
}
