import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * ImpliedTrinomialTree.java
 *
 * Java port of:
 *   Impl_Trinomial_Tree.py    -- builds an implied trinomial tree (Derman-Kani-Chriss)
 *                                from market-observed European call/put prices, recovering
 *                                state prices Q[j,i] and arbitrage-free transition
 *                                probabilities (pu, pm, pd).
 *   Amer_Up&Out_put.py        -- prices an American Up-and-Out put on the calibrated tree
 *                                via backward induction with early-exercise check and
 *                                knock-out at the barrier.
 *
 * The original Python loads European prices from 'Euro_prices_F5.3.xlsx'. To keep this
 * Java file self-contained and runnable without external XLSX dependencies, two input
 * modes are supported:
 *
 *   1) CSV mode: pass two CSV files via --call <path> --put <path>. Each file must be
 *      a (2N+1) x (N+1) grid of European call/put prices (same indexing as the
 *      original spreadsheet, where row 1 is the top-most node).
 *
 *   2) Black-Scholes fallback (default): if no CSVs are provided, European call/put
 *      prices at every (j, i) tree node are generated from a flat-vol Black-Scholes
 *      model so the pipeline still runs end-to-end. Useful as a sanity check --
 *      when fed BS prices, the implied tree should reproduce a (near-)constant local
 *      volatility surface.
 *
 * Compile & run:
 *   javac ImpliedTrinomialTree.java
 *   java  ImpliedTrinomialTree                                   # BS fallback
 *   java  ImpliedTrinomialTree --call calls.csv --put puts.csv   # market data mode
 */
public class ImpliedTrinomialTree {
    static final double T  = 1.0;     // maturity (years)
    static final double S  = 100.0;   // spot
    static final double r  = 0.05;    // risk-free rate
    static final int    N  = 4;       // number of time steps
    static final double dx = 0.2524;  // log-spacing (~ sigma * sqrt(3*dt) for sigma=0.2)

    static final double dt   = T / N;
    static final double edx  = Math.exp(dx);
    static final double infl = Math.exp(r * dt);  // one-step capitalisation factor

    static final double SIGMA_FLAT = 0.20;

    static double[][] callGrid;
    static double[][] putGrid;

    static double call(int x, int y) {
        return callGrid[x - 1][y - 1];
    }

    static double put(int x, int y) {
        return putGrid[x - 1][y - 1];
    }

    public static void main(String[] args) throws IOException {
        String callPath = null, putPath = null;
        for (int i = 0; i < args.length - 1; i++) {
            if (args[i].equals("--call")) callPath = args[i + 1];
            else if (args[i].equals("--put")) putPath = args[i + 1];
        }

        // Allocate the (2N+1) x (N+1) state matrices
        int rows = 2 * N + 1;
        int cols = N + 1;
        double[][] St = new double[rows][cols];   // asset price at terminal column N
        double[][] Q  = new double[rows][cols];   // Arrow-Debreu (state) prices
        double[][] pu = new double[rows][cols];   // up-transition probabilities
        double[][] pm = new double[rows][cols];   // middle-transition probabilities
        double[][] pd = new double[rows][cols];   // down-transition probabilities

        // ---- Build the terminal asset-price column St[*, N] ----
        // Python: St[2N, N] = S*exp(-N*dx); then St[j, N] = St[j+1, N] * edx for j = 2N-1..0
        St[2 * N][N] = S * Math.exp(-N * dx);
        for (int j = 2 * N - 1; j >= 0; j--) {
            St[j][N] = St[j + 1][N] * edx;
        }

        // ---- Load European prices ----
        if (callPath != null && putPath != null) {
            callGrid = readCsvGrid(callPath, rows, cols);
            putGrid  = readCsvGrid(putPath, rows, cols);
            System.out.println("Loaded market call/put prices from CSV.");
        } else {
            // Fall back to Black-Scholes prices on the same node grid (rows = strikes
            // taken from St[*, N], columns = maturities i*dt for i=1..N).
            callGrid = new double[rows][cols];
            putGrid  = new double[rows][cols];
            for (int rr = 0; rr < rows; rr++) {
                double K = St[rr][N];
                for (int cc = 0; cc < cols; cc++) {
                    double tau = (cc + 1) * dt; // mirror the (i+1) offset used in Python
                    callGrid[rr][cc] = bsCall(S, K, r, SIGMA_FLAT, tau);
                    putGrid [rr][cc] = bsPut (S, K, r, SIGMA_FLAT, tau);
                }
            }
            System.out.println("No CSVs supplied -- using Black-Scholes fallback prices "
                    + "(sigma = " + SIGMA_FLAT + ").");
        }

        // ---- Calibrate state prices Q[j, i] from European call/put prices ----
        // This faithfully reproduces the two-pass loop in Impl_Trinomial_Tree.py.
        for (int i = 0; i <= N; i++) {

            // Upper half of the tree: derive Q from European calls struck at St[j+1, N]
            for (int j = N - i; j < N + 1; j++) {
                double sum = 0.0;
                for (int k = j - 2; k < N; k++) {
                    if (k >= 0 && k < rows) {
                        sum += Q[k][i] * (St[k][N] - St[j + 1][N]);
                    }
                }
                double D = St[j][N] - St[j + 1][N];
                Q[j][i] = (call(j + 1, i + 1) - sum) / D;
            }

            // Lower half of the tree: derive Q from European puts struck at St[j-1, N]
            for (int j = N + i; j > N; j--) {
                double sum = 0.0;
                for (int k = N + i; k < j + 2; k++) {
                    if (k <= 8 && k < rows) {           // guard kept from the original
                        sum += Q[k][i] * (St[j - 1][N] - St[k][N]);
                    }
                }
                double D = St[j - 1][N] - St[j][N];
                Q[j][i] = (put(j - 1, i + 1) - sum) / D;
            }
        }

        // ---- Imply transition probabilities pu, pm, pd from Q ----
        // Forward equations relating Q at i+1 to Q at i, plus the forward-price
        // condition E[S_{i+1} | F_i] = infl * S_i.
        for (int i = 0; i < N; i++) {

            // Upper-half nodes: walk j from N-i upward
            for (int j = N - i; j < N + 1; j++) {

                if (j == N - i) {
                    pu[N - i][i] = infl * Q[N - i - 1][i + 1] / Q[N - i][i];
                    pm[j][i] = (infl * St[j][N] - St[j + 1][N]
                            - pu[j][i] * (St[j - 1][N] - St[j + 1][N]))
                            / (St[j][N] - St[j + 1][N]);
                    pd[j][i] = 1.0 - pm[j][i] - pu[j][i];

                } else if (j == N - i + 1) {
                    pu[N - i + 1][i] = (infl * Q[N - i][i + 1]
                            - pm[N - i][i] * Q[N - i][i]) / Q[N - i + 1][i];
                    pm[j][i] = (infl * St[j][N] - St[j + 1][N]
                            - pu[j][i] * (St[j - 1][N] - St[j + 1][N]))
                            / (St[j][N] - St[j + 1][N]);
                    pd[j][i] = 1.0 - pm[j][i] - pu[j][i];

                } else {
                    pu[j][i] = (infl * Q[j - 1][i + 1]
                            - pd[j - 2][i] * Q[j - 2][i]
                            - pm[j - 1][i] * Q[j - 1][i]) / Q[j][i];
                    pm[j][i] = (infl * St[j][N] - St[j + 1][N]
                            - pu[j][i] * (St[j - 1][N] - St[j + 1][N]))
                            / (St[j][N] - St[j + 1][N]);
                    pd[j][i] = 1.0 - pm[j][i] - pu[j][i];
                }
            }

            // Lower-half nodes: walk j from N+i downward
            for (int j = N + i; j > N; j--) {

                if (j == N + i) {
                    pd[N + i][i] = infl * Q[N + i + 1][i + 1] / Q[N + i][i];
                    pm[j][i] = (infl * St[j][N] - St[j - 1][N]
                            - pd[j][i] * (St[j + 1][N] - St[j - 1][N]))
                            / (St[j][N] - St[j - 1][N]);
                    pu[j][i] = 1.0 - pm[j][i] - pd[j][i];

                } else if (j == N + i - 1) {
                    pd[N + i - 1][i] = (infl * Q[N + i][i + 1]
                            - pm[N + i][i] * Q[N + i][i]) / Q[N + i - 1][i];
                    pm[j][i] = (infl * St[j][N] - St[j - 1][N]
                            - pd[j][i] * (St[j + 1][N] - St[j - 1][N]))
                            / (St[j][N] - St[j - 1][N]);
                    pu[j][i] = 1.0 - pm[j][i] - pd[j][i];

                } else {
                    pd[j][i] = (infl * Q[j + 1][i + 1]
                            - pu[j + 2][i] * Q[j + 2][i]
                            - pm[j + 1][i] * Q[j + 1][i]) / Q[j][i];
                    pm[j][i] = (infl * St[j][N] - St[j - 1][N]
                            - pd[j][i] * (St[j + 1][N] - St[j - 1][N]))
                            / (St[j][N] - St[j - 1][N]);
                    pu[j][i] = 1.0 - pm[j][i] - pd[j][i];
                }
            }
        }

        // ---- Diagnostic dumps (mirrors the np.savetxt calls in the Python) ----
        printMatrix("Asset prices St[j, N]",   St);
        printMatrix("State prices Q[j, i]",    Q);
        printMatrix("p_up   (pu[j, i])",       pu);
        printMatrix("p_mid  (pm[j, i])",       pm);
        printMatrix("p_down (pd[j, i])",       pd);

        // ---- Build the spot tree S[j, i] for pricing ----
        // The Python code only stores the terminal column, but for path-dependent
        // products we need the spot at every (j, i). Reconstruct it by stepping
        // forward from S using log-spacing dx (consistent with the lattice).
        double[][] spot = new double[rows][cols];
        spot[N][0] = S;                      // root node sits at row N, column 0
        for (int i = 1; i <= N; i++) {
            for (int j = N - i; j <= N + i; j++) {
                spot[j][i] = S * Math.exp((N - j) * dx);
            }
        }

        // ---- Price an American Up-and-Out put on the calibrated tree ----
        double K       = 100.0;   // strike
        double barrier = 110.0;   // up-and-out barrier
        double rebate  = 0.0;     // knock-out rebate
        double price = priceAmericanUpAndOutPut(spot, pu, pm, pd, K, barrier, rebate);

        System.out.println();
        System.out.printf("American Up-and-Out put%n");
        System.out.printf("  S0      = %.4f%n", S);
        System.out.printf("  K       = %.4f%n", K);
        System.out.printf("  Barrier = %.4f (up-and-out)%n", barrier);
        System.out.printf("  T       = %.4f, N = %d, r = %.4f%n", T, N, r);
        System.out.printf("  Price   = %.6f%n", price);
    }

    // ---------------------------------------------------------------------
    // American Up-and-Out put pricing on the implied trinomial tree
    // ---------------------------------------------------------------------
    /**
     * Backward induction:
     *   - At maturity: payoff = max(K - S, 0), but 0 (or rebate) if S >= barrier
     *     was ever crossed -- on a recombining tree the standard treatment is to
     *     knock out whenever the *current* node breaches the barrier, which is
     *     equivalent for monotone barriers when the lattice is dense enough.
     *   - At each earlier step:
     *         continuation = exp(-r*dt) * (pu*V_up + pm*V_mid + pd*V_down)
     *         exercise     = max(K - S, 0)
     *         V            = max(continuation, exercise)         if S < barrier
     *                      = rebate                              otherwise
     */
    static double priceAmericanUpAndOutPut(double[][] spot,
                                           double[][] pu, double[][] pm, double[][] pd,
                                           double K, double barrier, double rebate) {
        int rows = spot.length;
        double[][] V = new double[rows][N + 1];
        double disc = Math.exp(-r * dt);

        // Terminal payoffs
        for (int j = 0; j <= 2 * N; j++) {
            double sT = spot[j][N];
            if (sT >= barrier) {
                V[j][N] = rebate;
            } else {
                V[j][N] = Math.max(K - sT, 0.0);
            }
        }

        // Backward induction
        for (int i = N - 1; i >= 0; i--) {
            for (int j = N - i; j <= N + i; j++) {
                double s = spot[j][i];
                if (s >= barrier) {
                    V[j][i] = rebate;
                    continue;
                }
                // Children at column i+1: up = j-1, mid = j, down = j+1
                double cont = disc * (pu[j][i] * V[j - 1][i + 1]
                                    + pm[j][i] * V[j    ][i + 1]
                                    + pd[j][i] * V[j + 1][i + 1]);
                double exer = Math.max(K - s, 0.0);
                V[j][i] = Math.max(cont, exer);
            }
        }
        return V[N][0];
    }

    // ---------------------------------------------------------------------
    // Black-Scholes helpers (used only when no CSVs are supplied)
    // ---------------------------------------------------------------------
    static double bsCall(double S, double K, double r, double sigma, double tau) {
        if (tau <= 0) return Math.max(S - K, 0.0);
        double d1 = (Math.log(S / K) + (r + 0.5 * sigma * sigma) * tau) / (sigma * Math.sqrt(tau));
        double d2 = d1 - sigma * Math.sqrt(tau);
        return S * normalCdf(d1) - K * Math.exp(-r * tau) * normalCdf(d2);
    }

    static double bsPut(double S, double K, double r, double sigma, double tau) {
        if (tau <= 0) return Math.max(K - S, 0.0);
        double d1 = (Math.log(S / K) + (r + 0.5 * sigma * sigma) * tau) / (sigma * Math.sqrt(tau));
        double d2 = d1 - sigma * Math.sqrt(tau);
        return K * Math.exp(-r * tau) * normalCdf(-d2) - S * normalCdf(-d1);
    }

    /** Abramowitz-Stegun 7.1.26 approximation to the standard normal CDF. */
    static double normalCdf(double x) {
        double a1 =  0.254829592, a2 = -0.284496736, a3 = 1.421413741;
        double a4 = -1.453152027, a5 =  1.061405429, p  = 0.3275911;
        int sign = x < 0 ? -1 : 1;
        x = Math.abs(x) / Math.sqrt(2.0);
        double t = 1.0 / (1.0 + p * x);
        double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
        return 0.5 * (1.0 + sign * y);
    }

    // ---------------------------------------------------------------------
    // CSV reader -- expects rows separated by newlines, columns by commas.
    // ---------------------------------------------------------------------
    static double[][] readCsvGrid(String path, int rows, int cols) throws IOException {
        double[][] g = new double[rows][cols];
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            List<String> lines = new ArrayList<>();
            String line;
            while ((line = br.readLine()) != null) {
                if (!line.trim().isEmpty()) lines.add(line);
            }
            for (int i = 0; i < Math.min(rows, lines.size()); i++) {
                String[] toks = lines.get(i).split(",");
                for (int j = 0; j < Math.min(cols, toks.length); j++) {
                    String t = toks[j].trim();
                    g[i][j] = t.isEmpty() ? 0.0 : Double.parseDouble(t);
                }
            }
        }
        return g;
    }

    // ---------------------------------------------------------------------
    // Pretty-printer for diagnostic matrices
    // ---------------------------------------------------------------------
    static void printMatrix(String label, double[][] m) {
        System.out.println();
        System.out.println(label + ":");
        for (double[] row : m) {
            StringBuilder sb = new StringBuilder("  ");
            for (double v : row) sb.append(String.format("%10.4f ", v));
            System.out.println(sb.toString());
        }
    }
}
