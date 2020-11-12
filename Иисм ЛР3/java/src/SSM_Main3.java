import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.lang.reflect.Array;
import java.util.Arrays;
import static java.lang.Math.*;

public class SSM_Main3 {
    public static double erf(double z) {
        double t = 1.0 / (1.0 + 0.5 * Math.abs(z));

        // use Horner's method
        double ans = 1 - t * Math.exp( -z*z   -   1.26551223 +
                t * ( 1.00002368 +
                        t * ( 0.37409196 +
                                t * ( 0.09678418 +
                                        t * (-0.18628806 +
                                                t * ( 0.27886807 +
                                                        t * (-1.13520398 +
                                                                t * ( 1.48851587 +
                                                                        t * (-0.82215223 +
                                                                                t * ( 0.17087277))))))))));
        if (z >= 0) return  ans;
        else        return -ans;
    }

    public static double G(double x) {
        double[] p = {0.99999999999980993, 676.5203681218851, -1259.1392167224028,
                771.32342877765313, -176.61502916214059, 12.507343278686905,
                -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};
        int g = 7;
        if (x < 0.5) return Math.PI / (Math.sin(Math.PI * x) * G(1 - x));
        x -= 1;
        double a = p[0];
        double t = x + g + 0.5;
        for (int i = 1; i < p.length; i++) {
            a += p[i] / (x + i);
        }
        return Math.sqrt(2 * Math.PI) * Math.pow(t, x + 0.5) * Math.exp(-t) * a;
    }

    public static int fact(int x) {
        if (x == 0) return 1;
        if (x == 1) return 1;
        return x * fact(x - 1);
    }

    public static double g(int m, double x) {
        double res = (double)fact(m - 1);
        res *= exp(-x);
        double sum = 0;
        for (int i = 0; i <= m - 1; i++) {
            sum += pow(x, i) / fact(i);
        }
        res*=sum;
        return 1-res;
    }

    public static double B(double a, double b){
        return G(a)*G(b)/G(a+b);
    }

    public static double I(double x, double a, double b) {
        double res = pow(x, a) * pow((1 - x), b) / (a * B(a, b));
        double sum = 1;
        for (int n = 0; n < 20; n++) {
            sum += B(a + 1, n + 1) * pow(x, n + 1) / B(a + b, n + 1);
        }
        res *= sum;
        return res;
    }

    public static int[] freq(double[] A, int K) {
        Arrays.sort(A);
        double interval = (A[A.length - 1] - A[0]) / K;
        int[] frequencies = new int[K];
        for (int i = 0; i < A.length; i++)
            for (int j = 1; j <= K; j++)
                if (A[i] < A[0] + j * interval)
                    if (A[i] >= A[0] + (j - 1) * interval)
                        frequencies[j - 1]++;
        return frequencies;
    }

    public static double[] mult_cong(long alpha_0, long beta, int n, long M) {
        long buffer[] = new long[n];
        buffer[0] = (long) alpha_0;
        double A[] = new double[n];
        A[0] = (double) buffer[0] / M;
        for (int i = 1; i < n; i++) {
            long tmp = beta * buffer[i - 1];
            buffer[i] = tmp - M * (tmp / M);
            A[i] = (double) buffer[i] / M;
        }
        return A;
    }

    public static double mean(double[] A) {
        long sum = 0;
        for (int i = 0; i < A.length; i++)
            sum += A[i];
        double res = sum / (double) A.length;
        return res;
    }

    public static double dispersion(double[] A) {
        double mean_A = mean(A);
        double res = 0;
        for (int i = 0; i < A.length; i++)
            res += pow((A[i] - mean_A), 2);
        res = res / (double) (A.length - 1);
        return res;
    }

    public static double empiricalFunction(double[] A, double x) {
        double result = 0;
        for (int i = 0; i < A.length; i++) {
            if (A[i] <= x) result++;
        }
        return result / A.length;
    }


    public static double[] normal(int n, int m, int s_2) {
        double s = sqrt((double) s_2);
        int N = 12;
        double result[] = new double[n];
        double sum, eta;
        long alpha = 65643;
        double[] A = mult_cong(alpha, alpha, n * N, 2147483648L);
        for (int i = 0; i < n; i++) {
            sum = 0;
            for (int j = 0; j < N; j++) {
                sum += A[N * i + j];
            }
            eta = sqrt(12. / N) * (sum - N / 2.);
            result[i] = m + s * eta;
        }
        return result;
    }

    public static double Pearson_normal(double[] A) {
        int K = 10;
        double theor_freq, xi_2 = 0;
        int[] frequencies = freq(A, K);
        double interval = (A[A.length - 1] - A[0]) / K;
        double mean_ = mean(A);
        double s_2 = dispersion(A);

        for (int i = 0; i < K; i++) {
            double x = A[0] + ((i + 0.5) * interval);
            double u = pow(x - mean_, 2);
            double p_k = interval / sqrt(2 * PI * s_2) * exp(-u / (2 * s_2));
            theor_freq = A.length * p_k;
            xi_2 += pow((frequencies[i] - theor_freq), 2) / theor_freq;
        }
        return xi_2;
    }

    public static double Kolmogorov_normal(double[] A, int mu, int s_2) {
        double D_n = 0;
        double empiricalFRes, theoreticalFRes;
        for (int i = 0; i < A.length; i++) {
            empiricalFRes = empiricalFunction(A, A[i]);
            theoreticalFRes = theorNormFunction(A, A[i]);
            D_n = Math.max(D_n, Math.abs(empiricalFRes - theoreticalFRes));
        }
        return D_n*sqrt(A.length);
    }

    public static double theorNormFunction(double[] A, double x) {
        return 0.5 * (1. + erf((x - mean(A)) / sqrt(2 *  dispersion(A))));
    }


    public static double[] chi_2(int n, int m, long alpha) {
        double[] result = new double[n];
        int N = 12;
        double[] A = mult_cong(alpha, alpha, n * m * N, 2147483648L);
        double sum, big_sum;
        double eta;
        for (int i = 0; i < n; i++) {
            big_sum = 0;
            for (int j = 0; j < m; j++) {
                sum = 0;
                for (int k = 0; k < N; k++) {
                    sum += A[12 * m * i + 12 * j + k];
                }
                eta = sqrt(12. / N) * (sum - N / 2.);
                big_sum += pow(eta, 2);
            }
            result[i] = big_sum;
        }
        return result;
    }

    public static double Pearson_chi_2(double[] A, int m) {
        int K = 10;
        int s = 10;
        double theor_freq, xi_2 = 0;
        int[] frequencies = freq(A, K);
        double interval = (A[A.length - 1] - A[0]) / K;
        for (int i = 0; i < K; i++) {
            double p_k = 0;
            for (int j = 0; j < s; j++) {
                double x = A[0] + ((i + (double) j / s) * interval);
                double chisl = pow(x, (m - 2.) / 2) * exp(-x / 2);
                double znam = pow(2, m / 2.) * G(m / 2.);
                p_k += interval / s * chisl / znam;
            }
            theor_freq = A.length * p_k;
            xi_2 += pow((frequencies[i] - theor_freq), 2) / theor_freq;
        }
        return xi_2;
    }

    public static double Kolmogorov_chi_2(double[] A, int m) {
        Arrays.sort(A);
        double D_n = 0;
        double empiricalFRes, theoreticalFRes;
        for (int i = 0; i < A.length; i++) {
            empiricalFRes = empiricalFunction(A, A[i]);
            theoreticalFRes = theorChi2Function(A, A[i], m);
            D_n = Math.max(D_n, Math.abs(empiricalFRes - theoreticalFRes));
           // System.out.println(empiricalFRes+" "+theoreticalFRes);
        }
        return  D_n*sqrt(A.length);
    }

    public static double theorChi2Function(double[] A, double x, int m) {
        return g(m/2,x/2)/G(m / 2.);
    }


    public static double[] Fisher(int n, int l, int m) {
        double[] result = new double[n];
        long alpha1 = 262147;
        long alpha2 = 78125;
        double[] chi_l = chi_2(n, l, alpha1);
        double[] chi_m = chi_2(n, m, alpha2);
        for (int i = 0; i < n; i++) {
            result[i] = (chi_l[i] / (double) l) / (chi_m[i] / (double) m);
        }
        return result;
    }

    public static double Pearson_Fisher(double[] A, int l, int m) {
        int K = 10;
        int s = 50;
        int[] frequencies = freq(A, K);
        double interval = (A[A.length - 1] - A[0]) / K;
        double theor_freq, xi_2 = 0;

        double lm2 = (l + m) / 2.;
        double l22 = (l - 2.) / 2.;
        double lm = (double) l / m;
        double l2 = l / 2.;
        double m2 = m / 2.;

        for (int i = 0; i < K; i++) {
            double p_k = 0;
            for (int j = 0; j < s; j++) {
                double x = A[0] + ((i + (double) j / s) * interval);
                double m1x = 1. + l * x / m;
                p_k += interval / s * G(lm2) * pow(x, l22) * pow(lm, l2) /
                        (G(l2) * G(m2) * pow(m1x, lm2));
            }
            theor_freq = A.length * p_k;
            xi_2 += pow((frequencies[i] - theor_freq), 2) / theor_freq;
        }
        return xi_2;
    }

    public static double Kolmogorov_Fisher(double[] A, int l, int m) {
        double D_n = 0;
        double empiricalFRes, theoreticalFRes;
        for (int i = 0; i < A.length; i++) {
            empiricalFRes = empiricalFunction(A, A[i]);
            theoreticalFRes = theorFisherFunction(A, A[i],l,m);
            D_n = Math.max(D_n, Math.abs(empiricalFRes - theoreticalFRes));
        }
        return D_n;
    }

    public static double theorFisherFunction(double[] A, double x, int l, int m) {
        double k = l*x/(l*x+m);
        return I(k,l/2.,m/2.);
    }


    public static void main(String[] args) {
        try (PrintStream res = new PrintStream("result.txt")) {
            PrintStream out = new PrintStream("out.txt");
            int n = 1000;
            int mu = 0;
            int s_2 = 64;
            double[] A = normal(n, mu, s_2);
            double mean_ = mean(A);
            double sigma_2 = dispersion(A);
            for (int i = 0; i < n; i++) out.println(A[i]);
            res.println("Normal distribution: mu = " + mu + ", s_2 = " + s_2);
            res.println("size = " + n);
            res.println("mu = " + mu);
            res.println("mean = " + mean_);
            res.println("s_2 = " + s_2);
            res.println("dispersion = " + sigma_2);
            res.println("Kolmogorov test: " + Kolmogorov_normal(A, mu, s_2) + " < 1.36");
            res.println("Pearson test: " + Pearson_normal(A) + " < 16.9\n\n");

            PrintStream out1 = new PrintStream("out1.txt");
            int n1 = 1000;
            int m = 4;
            long alpha = 65643;
            double[] B = chi_2(n1, m, alpha);
            for (int i = 0; i < n1; i++) out1.println(B[i]);
            res.println("chi^2 distribution: m = " + m);
            res.println("size = " + n1);
            res.println("mu = " + m);
            res.println("mean = " + mean(B));
            res.println("s_2 = " + 2. * m);
            res.println("dispersion = " + dispersion(B));
            res.println("Kolmogorov test: " + Kolmogorov_chi_2(B,m) + " < 1.36");
            res.println("Pearson test: " + Pearson_chi_2(B, m) + " < 16.9\n\n");

            PrintStream out2 = new PrintStream("out2.txt");
            int n2 = 1000;
            int l = 3;
            int m_ = 5;
            double[] C = Fisher(n2, l, m_);
            for (int i = 0; i < n2; i++) out2.println(C[i]);
            res.println("Fisher's distribution: l = " + l + " m = " + m_);
            res.println("size = " + n2);
            res.println("mu = " + m_ / (m_ - 2));
            res.println("mean = " + mean(C));
            double s2 = (2. * m_ * m_ * (l + m_ - 2.)) /
                    (l * (m_ - 2.) * (m_ - 2.) * (m_ - 4.));
            res.println("s_2 = " + s2);
            res.println("dispersion = " + dispersion(C));
            res.println("Kolmogorov test: " + Kolmogorov_Fisher(C,l,m_) + " < 1.36");
            res.println("Pearson test: " + Pearson_Fisher(C, l, m_) + " < 16.9");
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }
}
