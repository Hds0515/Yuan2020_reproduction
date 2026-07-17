import com.comsol.model.*;
import com.comsol.model.util.*;

import java.util.HashMap;
import java.util.Locale;
import java.util.Map;

/**
 * COMSOL 6.4 implementation of the frozen V13 fixed-60%-PWM thermal network.
 *
 * This is a solver-independent equation implementation, not a reconstructed
 * fan plenum, manifold, or electrochemical CAD model. No parameter is fitted.
 */
public class V15Fixed60NetworkConfirmation {
  private static final int CELLS = 40;
  private static final int ROWS = 3;
  private static final int COLS = 5;
  private static final int N = CELLS * ROWS * COLS;
  private static final double AMBIENT_C = 23.0;

  private static final double ACTIVE_AREA = 0.0051;
  private static final double CHANNEL_LENGTH = 0.0608;
  private static final double STACK_LENGTH = 0.13642;
  private static final double PLATE_THICKNESS = 0.0025;
  private static final double CHANNEL_WIDTH = 0.0010;
  private static final double CHANNEL_DEPTH = 0.0014;
  private static final int CHANNEL_COUNT = 45;

  private static final double K_IN_PLANE = 60.0;
  private static final double K_STACK = 2.8745487024696486;
  private static final double CP_AIR = 1006.0;
  private static final double TOTAL_CAPACITY = 3999.999989994387;
  private static final double SURFACE_FRACTION = 0.010010000000056184;
  private static final double SURFACE_CORE_G = 8.065428145396455;
  private static final double FLOW_57 = 0.021934345212215003 * 2.072349298596786;
  private static final double PWM_EXPONENT = 2.5917841202975773;
  private static final double HT_MULTIPLIER = 1.5434528081728003;
  private static final double EDGE_H = 13.9999972210056;
  private static final double HEAT_NONUNIFORMITY = 0.3721976849771562;
  private static final double HUB_DEFICIT = 0.449999;
  private static final double AXIAL_SKEW = 0.4128642328339601;

  private static int index(int cell, int row, int col) {
    return (cell * ROWS + row) * COLS + col;
  }

  private static void add(Map<Integer, Double> row, int col, double value) {
    row.put(col, row.getOrDefault(col, 0.0) + value);
  }

  private static void connect(Map<Integer, Double>[] matrix, int a, int b, double g) {
    add(matrix[a], a, -g);
    add(matrix[b], b, -g);
    add(matrix[a], b, g);
    add(matrix[b], a, g);
  }

  @SuppressWarnings("unchecked")
  private static Map<Integer, Double>[] surfaceMatrix() {
    Map<Integer, Double>[] matrix = new Map[N];
    for (int i = 0; i < N; i++) matrix[i] = new HashMap<>();

    double dx = (ACTIVE_AREA / CHANNEL_LENGTH) / COLS;
    double dy = CHANNEL_LENGTH / ROWS;
    double dz = STACK_LENGTH / CELLS;
    double gx = K_IN_PLANE * PLATE_THICKNESS * dy / dx;
    double gy = K_IN_PLANE * PLATE_THICKNESS * dx / dy;
    double gz = K_STACK * dx * dy / dz;
    for (int cell = 0; cell < CELLS; cell++) {
      for (int row = 0; row < ROWS; row++) {
        for (int col = 0; col < COLS; col++) {
          int node = index(cell, row, col);
          if (col + 1 < COLS) connect(matrix, node, index(cell, row, col + 1), gx);
          if (row + 1 < ROWS) connect(matrix, node, index(cell, row + 1, col), gy);
          if (cell + 1 < CELLS) connect(matrix, node, index(cell + 1, row, col), gz);
        }
      }
    }

    double normalized = (0.60 - 0.15) / (0.57 - 0.15);
    double totalFlow = FLOW_57 * (0.08 + 0.92 * Math.pow(normalized, PWM_EXPONENT));
    double[][] weights = new double[CELLS][COLS];
    double sumRawWeight = 0.0;
    for (int cell = 0; cell < CELLS; cell++) {
      double z = (cell + 0.5) / CELLS;
      for (int col = 0; col < COLS; col++) {
        double x = (col + 0.5) / COLS;
        double rz = Math.abs(z - 0.5) / 0.5;
        double rx = Math.abs(x - 0.5) / 0.5;
        double radius = Math.sqrt(rz * rz + 0.35 * rx * rx);
        double raw = 1.0 - HUB_DEFICIT * Math.exp(-Math.pow(radius / 0.38, 2.0));
        weights[cell][col] = raw;
        sumRawWeight += raw;
      }
    }
    double rawMean = sumRawWeight / (CELLS * COLS);
    double sumWeight = 0.0;
    for (int cell = 0; cell < CELLS; cell++) {
      for (int col = 0; col < COLS; col++) {
        weights[cell][col] = 1.0 + (weights[cell][col] / rawMean - 1.0) * (1.0 - 0.55);
        sumWeight += weights[cell][col];
      }
    }
    double meanFlow = totalFlow / (CELLS * COLS);
    double wettedArea = (CHANNEL_COUNT / (double) COLS)
        * (CHANNEL_WIDTH + 2.0 * CHANNEL_DEPTH) * dy;
    for (int cell = 0; cell < CELLS; cell++) {
      for (int col = 0; col < COLS; col++) {
        double mdot = totalFlow * weights[cell][col] / sumWeight;
        double h = 85.0 * HT_MULTIPLIER * Math.pow(mdot / meanFlow, 0.33);
        double effectiveness = 1.0 - Math.exp(-h * wettedArea / (mdot * CP_AIR));
        double mcp = mdot * CP_AIR;
        Map<Integer, Double> upstream = new HashMap<>();
        for (int row = 0; row < ROWS; row++) {
          int node = index(cell, row, col);
          double removal = mcp * effectiveness;
          add(matrix[node], node, -removal);
          for (Map.Entry<Integer, Double> entry : upstream.entrySet()) {
            add(matrix[node], entry.getKey(), removal * entry.getValue());
          }
          Map<Integer, Double> next = new HashMap<>();
          for (Map.Entry<Integer, Double> entry : upstream.entrySet()) {
            next.put(entry.getKey(), (1.0 - effectiveness) * entry.getValue());
          }
          next.put(node, next.getOrDefault(node, 0.0) + effectiveness);
          upstream = next;
        }
      }
    }

    double ambientK = AMBIENT_C + 273.15;
    double referenceK = 55.0 + 273.15;
    double hRad = 0.82 * 5.670374419e-8
        * (referenceK * referenceK + ambientK * ambientK) * (referenceK + ambientK);
    double hExternal = EDGE_H + hRad;
    for (int cell = 0; cell < CELLS; cell++) {
      for (int row = 0; row < ROWS; row++) {
        for (int col = 0; col < COLS; col++) {
          double area = 0.0;
          if (col == 0 || col == COLS - 1) area += dz * dy;
          if (row == 0 || row == ROWS - 1) area += dz * dx;
          if (cell == 0 || cell == CELLS - 1) area += dx * dy;
          add(matrix[index(cell, row, col)], index(cell, row, col), -hExternal * area);
        }
      }
    }
    return matrix;
  }

  private static String heatExpression(int node) {
    int row = (node / COLS) % ROWS;
    int col = node % COLS;
    double[] pattern = {-0.65, 0.10, 0.28, 0.52, -0.25};
    double columnFactor = 1.0 + HEAT_NONUNIFORMITY * pattern[col];
    double rowCoordinate = 2.0 * (row + 0.5) / ROWS - 1.0;
    return String.format(Locale.US,
        "Qgen/%d*(%.17g)*(1+(%.17g)*(Iapp/(40[A]))*(%.17g))",
        N, columnFactor, AXIAL_SKEW, rowCoordinate);
  }

  private static String surfaceEquation(int node, Map<Integer, Double> row) {
    StringBuilder rhs = new StringBuilder("(");
    boolean first = true;
    for (Map.Entry<Integer, Double> entry : row.entrySet()) {
      if (!first && entry.getValue() >= 0.0) rhs.append("+");
      rhs.append(String.format(Locale.US, "%.17g[W/K]*s%d", entry.getValue(), entry.getKey()));
      first = false;
    }
    rhs.append(String.format(Locale.US, "+gsc*(c%d-s%d)+", node, node));
    rhs.append(heatExpression(node)).append(")");
    return String.format(Locale.US, "cs*s%dt-%s", node, rhs);
  }

  private static double initialFor(int node) {
    int row = (node / COLS) % ROWS;
    double[] initialC = {27.43911439114391, 27.793357933579337, 28.14760147601476};
    return initialC[row] - AMBIENT_C;
  }

  private static String regionExpression(int row) {
    int[] cells = {3, 11, 19, 27, 35};
    StringBuilder expression = new StringBuilder("(");
    int count = 0;
    for (int cell : cells) {
      for (int col = 0; col < COLS; col++) {
        if (count++ > 0) expression.append("+");
        expression.append("s").append(index(cell, row, col));
      }
    }
    expression.append(")/").append(count).append("+").append(AMBIENT_C)
        .append("[degC]-273.15[K]");
    return expression.toString();
  }

  private static String storageExpression() {
    StringBuilder expression = new StringBuilder("(");
    for (int node = 0; node < N; node++) {
      if (node > 0) expression.append("+");
      expression.append("cs*s").append(node).append("t+cc*c").append(node).append("t");
    }
    return expression.append(")").toString();
  }

  private static String rejectionExpression(Map<Integer, Double>[] matrix) {
    StringBuilder expression = new StringBuilder("-(");
    boolean first = true;
    for (int row = 0; row < N; row++) {
      for (Map.Entry<Integer, Double> entry : matrix[row].entrySet()) {
        if (!first && entry.getValue() >= 0.0) expression.append("+");
        expression.append(String.format(Locale.US, "%.17g[W/K]*s%d",
            entry.getValue(), entry.getKey()));
        first = false;
      }
    }
    return expression.append(")").toString();
  }

  private static Model build() {
    Model model = ModelUtil.create("Model");
    model.label("v15_fixed60_network_confirmation.mph");
    model.param().set("cs", (TOTAL_CAPACITY * SURFACE_FRACTION / N) + "[J/K]");
    model.param().set("cc", (TOTAL_CAPACITY * (1.0 - SURFACE_FRACTION) / N) + "[J/K]");
    model.param().set("gsc", (SURFACE_CORE_G / N) + "[W/K]");
    model.component().create("comp1", true);
    model.component("comp1").variable().create("var1");
    model.component("comp1").variable("var1").set("Iapp",
        "if(t<1000[s],min(4[A]*(floor(t/100[s])+1),40[A]),"
        + "max(36[A]-4[A]*floor((t-1000[s])/100[s]),0[A]))");
    model.component("comp1").variable("var1").set("Vstack",
        "min(max(33.5[V]-0.34[V/A]*Iapp,19.5[V]),33.5[V])");
    model.component("comp1").variable("var1").set("Qgen",
        "max(Iapp*(59.2[V]-Vstack),0[W])");

    model.component("comp1").physics().create("ge", "GlobalEquations");
    Map<Integer, Double>[] matrix = surfaceMatrix();
    for (int node = 0; node < N; node++) {
      model.component("comp1").physics("ge").feature("ge1")
          .setIndex("name", "s" + node, node, 0);
      model.component("comp1").physics("ge").feature("ge1")
          .setIndex("equation", surfaceEquation(node, matrix[node]), node, 0);
      model.component("comp1").physics("ge").feature("ge1")
          .setIndex("initialValueU", initialFor(node) + "[K]", node, 0);
    }
    for (int node = 0; node < N; node++) {
      int core = N + node;
      model.component("comp1").physics("ge").feature("ge1")
          .setIndex("name", "c" + node, core, 0);
      model.component("comp1").physics("ge").feature("ge1").setIndex("equation",
          String.format(Locale.US, "cc*c%dt-gsc*(s%d-c%d)", node, node, node), core, 0);
      model.component("comp1").physics("ge").feature("ge1")
          .setIndex("initialValueU", initialFor(node) + "[K]", core, 0);
    }

    model.study().create("std1");
    model.study("std1").feature().create("time", "Transient");
    model.study("std1").feature("time").set("tlist", "range(50,50,1950)");
    model.study("std1").feature("time").set("rtol", "1e-7");
    return model;
  }

  private static void writeResults(Model model) throws Exception {
    String storage = storageExpression();
    String rejection = rejectionExpression(surfaceMatrix());
    String residual = "abs(Qgen-(" + rejection + ")-(" + storage
        + "))/max(abs(Qgen),1[W])*100";
    String[] expr = {regionExpression(0), regionExpression(1), regionExpression(2),
        "max(" + regionExpression(0) + ",max(" + regionExpression(1) + "," + regionExpression(2) + "))",
        "Qgen", rejection, storage, residual};
    model.result().numerical().create("eval", "EvalGlobal");
    model.result().numerical("eval").set("expr", expr);
    model.result().table().create("tbl", "Table");
    model.result().numerical("eval").set("table", "tbl");
    model.result().numerical("eval").setResult();
    model.result().export().create("tableExport", "Table");
    model.result().export("tableExport").set("table", "tbl");
    model.result().export("tableExport").set(
        "filename", "outputs_v15/comsol_fixed60_regional_timeseries.csv");
    model.result().export("tableExport").run();
  }

  public static void main(String[] args) throws Exception {
    Model model = build();
    System.out.println("V15_COMSOL_SOLVE_BEGIN states=" + (2 * N));
    long start = System.nanoTime();
    model.study("std1").run();
    double runtime = (System.nanoTime() - start) / 1e9;
    System.out.printf(Locale.US, "V15_COMSOL_SOLVE_CONVERGED runtime_s=%.9g%n", runtime);
    writeResults(model);
    model.save("comsol/v15_reference/v15_fixed60_network_confirmation.mph");
    System.out.println("V15_COMSOL_RUN_COMPLETE");
  }
}
