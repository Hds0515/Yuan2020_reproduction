import com.comsol.model.*;
import com.comsol.model.util.*;

/**
 * V7 prescribed-velocity conjugate heat model for COMSOL 6.4.
 *
 * Table 2 material properties and one consistent representative-cell scaling
 * are frozen.  Only hExternal and emissivity are calibrated, using the 4 and
 * 8 m/s Fig.7 regional temperatures.  Speeds 5, 6, 10 and 12 m/s never enter
 * calibration.
 */
public class Stage2PhysicalHeat {
  private static final double[] VELOCITIES = {4, 5, 6, 8, 10, 12};
  private static final double[][] CAL_TARGET = {
    {55.3365085188, 60.2981218807, 62.9111747982},
    {45.7197598698, 49.1994505495, 51.8807692308}
  };
  private static final String[] METRIC_NAMES = {
    "T_inlet_region_C", "T_middle_region_C", "T_outlet_region_C",
    "Tmax_C", "Tmin_C", "DeltaT_C", "hotspot_location_normalized",
    "hotspot_y_m", "hotspot_z_m", "air_inlet_C", "air_outlet_C",
    "air_enthalpy_gain_W", "natural_convection_loss_W", "radiation_loss_W",
    "conduction_loss_W", "total_heat_input_W", "energy_residual_relative",
    "mdot_cell_kg_s", "mdot_stack_kg_s", "mdot_identity_relative_error"
  };

  private static void box(Model model, String tag, int dim,
      String xmin, String xmax, String ymin, String ymax, String zmin, String zmax) {
    model.component("comp1").selection().create(tag, "Box");
    model.component("comp1").selection(tag).set("entitydim", dim);
    model.component("comp1").selection(tag).set("condition", "inside");
    model.component("comp1").selection(tag).set("xmin", xmin);
    model.component("comp1").selection(tag).set("xmax", xmax);
    model.component("comp1").selection(tag).set("ymin", ymin);
    model.component("comp1").selection(tag).set("ymax", ymax);
    model.component("comp1").selection(tag).set("zmin", zmin);
    model.component("comp1").selection(tag).set("zmax", zmax);
  }

  private static void union(Model model, String tag, int dim, String... inputs) {
    model.component("comp1").selection().create(tag, "Union");
    model.component("comp1").selection(tag).set("entitydim", dim);
    model.component("comp1").selection(tag).set("input", inputs);
  }

  private static Model buildModel(
      String meshName, int nLength, int nWidth, int nAir, int nSolid) {
    ModelUtil.clear();
    Model model = ModelUtil.create("Model");
    model.label("v7_physical_heat_" + meshName + ".mph");

    model.param().set("AinStack", "3962[mm^2]");
    model.param().set("nCells", "40");
    model.param().set("AinCell", "AinStack/nCells");
    model.param().set("Hair", "1[mm]");
    model.param().set("W", "AinCell/Hair");
    model.param().set("AheatCell", "0.0247[m^2]");
    model.param().set("L", "AheatCell/W");
    model.param().set("Hs", "10[mm]");
    model.param().set("uin", "4[m/s]");
    model.param().set("Tin", "25[degC]");
    model.param().set("Tamb", "25[degC]");
    model.param().set("qflux", "2425.5[W/m^2]");
    model.param().set("rhoAir", "1.184[kg/m^3]");
    model.param().set("CpAir", "1007[J/(kg*K)]");
    model.param().set("sigmaSB", "5.670374419e-8[W/(m^2*K^4)]");
    model.param().set("hExternal", "10[W/(m^2*K)]");
    model.param().set("emissivity", "0.8");
    model.param().set("tol", "1e-7[m]");

    model.component().create("comp1", true);
    model.component("comp1").geom().create("geom1", 3);
    model.component("comp1").geom("geom1").lengthUnit("m");
    model.component("comp1").geom("geom1").feature().create("air", "Block");
    model.component("comp1").geom("geom1").feature("air")
        .set("size", new String[]{"L", "W", "Hair"});
    model.component("comp1").geom("geom1").feature("air").set("selresult", "on");
    model.component("comp1").geom("geom1").feature().create("solid", "Block");
    model.component("comp1").geom("geom1").feature("solid")
        .set("size", new String[]{"L", "W", "Hs"});
    model.component("comp1").geom("geom1").feature("solid")
        .set("pos", new String[]{"0", "0", "Hair"});
    model.component("comp1").geom("geom1").feature("solid").set("selresult", "on");
    model.component("comp1").geom("geom1").run();

    box(model, "selAirIn", 2, "-tol", "tol", "-tol", "W+tol", "-tol", "Hair+tol");
    box(model, "selAirOut", 2, "L-tol", "L+tol", "-tol", "W+tol", "-tol", "Hair+tol");
    box(model, "selSource", 2, "-tol", "tol", "-tol", "W+tol", "-tol", "Hair+Hs+tol");
    box(model, "selTarget", 2, "L-tol", "L+tol", "-tol", "W+tol", "-tol", "Hair+Hs+tol");
    box(model, "selHeatTop", 2, "-tol", "L+tol", "-tol", "W+tol",
        "Hair+Hs-tol", "Hair+Hs+tol");
    box(model, "selSolidX0", 2, "-tol", "tol", "-tol", "W+tol", "Hair-tol", "Hair+Hs+tol");
    box(model, "selSolidXL", 2, "L-tol", "L+tol", "-tol", "W+tol", "Hair-tol", "Hair+Hs+tol");
    box(model, "selSolidY0", 2, "-tol", "L+tol", "-tol", "tol", "Hair-tol", "Hair+Hs+tol");
    box(model, "selSolidYW", 2, "-tol", "L+tol", "W-tol", "W+tol", "Hair-tol", "Hair+Hs+tol");
    union(model, "selExternalSolid", 2, "selHeatTop", "selSolidX0", "selSolidXL", "selSolidY0", "selSolidYW");
    model.component("comp1").selection().create("selAllDomains", "Explicit");
    model.component("comp1").selection("selAllDomains").geom("geom1", 3);
    model.component("comp1").selection("selAllDomains").all();

    box(model, "selWidthBottom", 1, "-tol", "tol", "-tol", "W+tol", "-tol", "tol");
    box(model, "selWidthInterface", 1, "-tol", "tol", "-tol", "W+tol", "Hair-tol", "Hair+tol");
    box(model, "selWidthTop", 1, "-tol", "tol", "-tol", "W+tol", "Hair+Hs-tol", "Hair+Hs+tol");
    union(model, "selWidthEdges", 1, "selWidthBottom", "selWidthInterface", "selWidthTop");
    box(model, "selAirZ0", 1, "-tol", "tol", "-tol", "tol", "-tol", "Hair+tol");
    box(model, "selAirZW", 1, "-tol", "tol", "W-tol", "W+tol", "-tol", "Hair+tol");
    union(model, "selAirThickness", 1, "selAirZ0", "selAirZW");
    box(model, "selSolidZ0", 1, "-tol", "tol", "-tol", "tol", "Hair-tol", "Hair+Hs+tol");
    box(model, "selSolidZW", 1, "-tol", "tol", "W-tol", "W+tol", "Hair-tol", "Hair+Hs+tol");
    union(model, "selSolidThickness", 1, "selSolidZ0", "selSolidZW");

    if (model.component("comp1").selection("selAirIn").entities().length != 1
        || model.component("comp1").selection("selAirOut").entities().length != 1
        || model.component("comp1").selection("selHeatTop").entities().length != 1) {
      throw new IllegalStateException("V7 thermal selection audit failed");
    }

    model.component("comp1").material().create("matAir", "Common");
    model.component("comp1").material("matAir").selection().named("geom1_air_dom");
    model.component("comp1").material("matAir").propertyGroup("def").set("density", "rhoAir");
    model.component("comp1").material("matAir").propertyGroup("def").set("heatcapacity", "CpAir");
    model.component("comp1").material("matAir").propertyGroup("def")
        .set("thermalconductivity", "0.0251[W/(m*K)]");
    model.component("comp1").material("matAir").propertyGroup("def")
        .set("dynamicviscosity", "1.849e-5[Pa*s]");

    model.component("comp1").material().create("matGraphite", "Common");
    model.component("comp1").material("matGraphite").selection().named("geom1_solid_dom");
    model.component("comp1").material("matGraphite").propertyGroup("def")
        .set("density", "2250[kg/m^3]");
    model.component("comp1").material("matGraphite").propertyGroup("def")
        .set("heatcapacity", "460[J/(kg*K)]");
    model.component("comp1").material("matGraphite").propertyGroup("def")
        .set("thermalconductivity", "24[W/(m*K)]");

    model.component("comp1").physics().create("ht", "HeatTransferInSolidsAndFluids", "geom1");
    model.component("comp1").physics("ht").feature("fluid1").selection().named("geom1_air_dom");
    model.component("comp1").physics("ht").feature("fluid1").set("u_src", "userdef");
    model.component("comp1").physics("ht").feature("fluid1").set("u", new String[]{"uin", "0", "0"});
    model.component("comp1").physics("ht").feature().create("in1", "Inflow", 2);
    model.component("comp1").physics("ht").feature("in1").selection().named("selAirIn");
    model.component("comp1").physics("ht").feature("in1").set("Tustr", "Tin");
    model.component("comp1").physics("ht").feature().create("heat1", "HeatFluxBoundary", 2);
    model.component("comp1").physics("ht").feature("heat1").selection().named("selHeatTop");
    model.component("comp1").physics("ht").feature("heat1").set("q0", "qflux");
    model.component("comp1").physics("ht").feature().create("conv1", "HeatFluxBoundary", 2);
    model.component("comp1").physics("ht").feature("conv1").selection().named("selExternalSolid");
    model.component("comp1").physics("ht").feature("conv1").set("HeatFluxType", "ConvectiveHeatFlux");
    model.component("comp1").physics("ht").feature("conv1").set("h", "hExternal");
    model.component("comp1").physics("ht").feature("conv1").set("Text", "Tamb");
    model.component("comp1").physics("ht").feature().create("rad1", "SurfaceToAmbientRadiation", 2);
    model.component("comp1").physics("ht").feature("rad1").selection().named("selExternalSolid");
    model.component("comp1").physics("ht").feature("rad1").set("epsilon_rad_mat", "userdef");
    model.component("comp1").physics("ht").feature("rad1").set("epsilon_rad", "emissivity");
    model.component("comp1").physics("ht").feature("rad1").set("Tamb", "Tamb");
    model.component("comp1").physics("ht").feature("init1").set("Tinit", "Tin");

    model.component("comp1").mesh().create("mesh1");
    String[] edgeTags = {"edgWidth", "edgAir", "edgSolid"};
    String[] edgeSelections = {"selWidthEdges", "selAirThickness", "selSolidThickness"};
    int[] edgeCounts = {nWidth, nAir, nSolid};
    for (int i = 0; i < edgeTags.length; i++) {
      model.component("comp1").mesh("mesh1").feature().create(edgeTags[i], "Edge");
      model.component("comp1").mesh("mesh1").feature(edgeTags[i]).selection().named(edgeSelections[i]);
      model.component("comp1").mesh("mesh1").feature(edgeTags[i]).create("dis1", "Distribution");
      model.component("comp1").mesh("mesh1").feature(edgeTags[i]).feature("dis1")
          .set("numelem", edgeCounts[i]);
    }
    model.component("comp1").mesh("mesh1").feature().create("map1", "Map");
    model.component("comp1").mesh("mesh1").feature("map1").selection().named("selSource");
    model.component("comp1").mesh("mesh1").feature().create("swe1", "Sweep");
    model.component("comp1").mesh("mesh1").feature("swe1").selection().named("selAllDomains");
    model.component("comp1").mesh("mesh1").feature("swe1").selection("sourceface").named("selSource");
    model.component("comp1").mesh("mesh1").feature("swe1").selection("targetface").named("selTarget");
    model.component("comp1").mesh("mesh1").feature("swe1").create("dis1", "Distribution");
    model.component("comp1").mesh("mesh1").feature("swe1").feature("dis1").set("numelem", nLength);
    model.component("comp1").mesh("mesh1").run();

    model.component("comp1").cpl().create("intSolid", "Integration");
    model.component("comp1").cpl("intSolid").selection().named("geom1_solid_dom");
    model.component("comp1").cpl().create("intHeat", "Integration");
    model.component("comp1").cpl("intHeat").selection().named("selHeatTop");
    model.component("comp1").cpl().create("intExt", "Integration");
    model.component("comp1").cpl("intExt").selection().named("selExternalSolid");
    model.component("comp1").cpl().create("aveIn", "Average");
    model.component("comp1").cpl("aveIn").selection().named("selAirIn");
    model.component("comp1").cpl().create("aveOut", "Average");
    model.component("comp1").cpl("aveOut").selection().named("selAirOut");
    model.component("comp1").cpl().create("maxSolid", "Maximum");
    model.component("comp1").cpl("maxSolid").selection().named("geom1_solid_dom");
    model.component("comp1").cpl().create("minSolid", "Minimum");
    model.component("comp1").cpl("minSolid").selection().named("geom1_solid_dom");

    model.study().create("std1");
    model.study("std1").feature().create("stat", "Stationary");
    model.study("std1").feature().create("param", "Parametric");
    model.study("std1").feature("param").set("pname", new String[]{"uin"});
    model.study("std1").feature("param").set("plistarr", new String[]{"4 8"});
    model.study("std1").feature("param").set("punit", new String[]{"m/s"});
    model.result().numerical().create("gevMetrics", "EvalGlobal");
    model.result().numerical("gevMetrics").set("expr", expressions());
    return model;
  }

  private static String[] expressions() {
    String z1 = "intSolid(if(x<L/3,T-273.15[K],0[K]))/intSolid(if(x<L/3,1,0))";
    String z2 = "intSolid(if(x>=L/3&&x<2*L/3,T-273.15[K],0[K]))/intSolid(if(x>=L/3&&x<2*L/3,1,0))";
    String z3 = "intSolid(if(x>=2*L/3,T-273.15[K],0[K]))/intSolid(if(x>=2*L/3,1,0))";
    String qIn = "intHeat(qflux)";
    String mdot = "rhoAir*uin*AinCell";
    String qAir = mdot + "*CpAir*(aveOut(T)-Tin)";
    String qConv = "intExt(hExternal*(T-Tamb))";
    String qRad = "intExt(emissivity*sigmaSB*(T^4-Tamb^4))";
    String balance = "abs((" + qIn + ")-(" + qAir + ")-(" + qConv + ")-(" + qRad + "))/abs(" + qIn + ")";
    return new String[]{
      z1, z2, z3, "maxSolid(T)-273.15[K]", "minSolid(T)-273.15[K]",
      "maxSolid(T)-minSolid(T)", "maxSolid(T,x)/L", "maxSolid(T,y)", "maxSolid(T,z)",
      "aveIn(T)-273.15[K]", "aveOut(T)-273.15[K]", qAir, qConv, qRad,
      "0[W]", qIn, balance, mdot, "nCells*(" + mdot + ")",
      "abs(nCells*(" + mdot + ")-rhoAir*uin*AinStack)/(rhoAir*uin*AinStack)"
    };
  }

  private static double[][] evaluate(Model model) {
    // Recreate after each study run so the evaluation binds to the newly
    // generated parametric solution instead of an empty pre-solve dataset.
    model.result().numerical().remove("gevMetrics");
    model.result().numerical().create("gevMetrics", "EvalGlobal");
    model.result().numerical("gevMetrics").set("expr", expressions());
    return model.result().numerical("gevMetrics").getReal();
  }

  private static double calibrationRmse(double[][] values) {
    double sum = 0.0;
    int count = 0;
    for (int speedIndex = 0; speedIndex < 2; speedIndex++) {
      for (int region = 0; region < 3; region++) {
        double error = values[region][speedIndex] - CAL_TARGET[speedIndex][region];
        sum += error * error;
        count++;
      }
    }
    return Math.sqrt(sum / count);
  }

  private static double[] calibrate() {
    Model model = buildModel("calibration_coarse", 30, 12, 2, 4);
    double bestH = 0.0, bestEps = 0.0, bestRmse = Double.POSITIVE_INFINITY;
    for (double h = 0.0; h <= 30.0001; h += 5.0) {
      for (double eps = 0.0; eps <= 1.0001; eps += 0.2) {
        model.param().set("hExternal", h + "[W/(m^2*K)]");
        model.param().set("emissivity", Double.toString(eps));
        model.study("std1").run();
        double rmse = calibrationRmse(evaluate(model));
        System.out.printf("V7_CAL_GRID h=%.6g eps=%.6g rmse_C=%.12g%n", h, eps, rmse);
        if (rmse < bestRmse) { bestH = h; bestEps = eps; bestRmse = rmse; }
      }
    }
    double coarseH = bestH, coarseEps = bestEps;
    for (double h = Math.max(0.0, coarseH - 4.0); h <= Math.min(30.0, coarseH + 4.0) + 1e-9; h += 1.0) {
      for (double eps = Math.max(0.0, coarseEps - 0.15); eps <= Math.min(1.0, coarseEps + 0.15) + 1e-9; eps += 0.05) {
        model.param().set("hExternal", h + "[W/(m^2*K)]");
        model.param().set("emissivity", Double.toString(eps));
        model.study("std1").run();
        double rmse = calibrationRmse(evaluate(model));
        System.out.printf("V7_CAL_REFINE h=%.6g eps=%.6g rmse_C=%.12g%n", h, eps, rmse);
        if (rmse < bestRmse) { bestH = h; bestEps = eps; bestRmse = rmse; }
      }
    }
    System.out.printf("V7_CALIBRATION_FROZEN hExternal_W_m2K=%.12g emissivity=%.12g calibration_rmse_C=%.12g%n",
        bestH, bestEps, bestRmse);
    return new double[]{bestH, bestEps, bestRmse};
  }

  private static double[][] solveFinal(String meshName, int nLength, int nWidth,
      int nAir, int nSolid, double h, double eps, boolean images) throws Exception {
    Model model = buildModel(meshName, nLength, nWidth, nAir, nSolid);
    model.param().set("hExternal", h + "[W/(m^2*K)]");
    model.param().set("emissivity", Double.toString(eps));
    model.study("std1").feature("param").set("plistarr", new String[]{"4 5 6 8 10 12"});
    System.out.println("V7_HEAT_SOLVE_BEGIN mesh=" + meshName);
    model.study("std1").run();
    System.out.println("V7_HEAT_SOLVE_CONVERGED mesh=" + meshName);
    double[][] values = evaluate(model);
    for (int j = 0; j < VELOCITIES.length; j++) {
      System.out.printf("V7_HEAT_RESULT mesh=%s velocity_m_s=%.0f", meshName, VELOCITIES[j]);
      for (int i = 0; i < METRIC_NAMES.length; i++) {
        System.out.printf(" %s=%.12g", METRIC_NAMES[i], values[i][j]);
      }
      System.out.println();
    }
    model.result().table().create("tblFinal", "Table");
    model.result().numerical("gevMetrics").set("table", "tblFinal");
    model.result().numerical("gevMetrics").setResult();
    model.result().export().create("tblExport", "Table");
    model.result().export("tblExport").set("table", "tblFinal");
    model.result().export("tblExport").set("filename", "comsol/v7_runtime/physical_heat_metrics_" + meshName + ".csv");
    model.result().export("tblExport").run();
    if (images) {
      try {
        model.result().create("pgTemp", "PlotGroup3D");
        model.result("pgTemp").feature().create("surf1", "Surface");
        model.result("pgTemp").feature("surf1").set("expr", "T-273.15[K]");
        model.result("pgTemp").feature("surf1").set("colortable", "RainbowLight");
        model.result("pgTemp").feature("surf1").set("rangecoloractive", "on");
        model.result("pgTemp").feature("surf1").set("rangecolormin", "25");
        model.result("pgTemp").feature("surf1").set("rangecolormax", "70");
        model.result().export().create("imgTemp", "pgTemp", "Image");
        model.result().export("imgTemp").set("width", 1200);
        model.result().export("imgTemp").set("height", 700);
        model.result().export("imgTemp").set("zoomextents", "on");
        for (int j = 0; j < VELOCITIES.length; j++) {
          model.result("pgTemp").set("looplevel", new int[]{j + 1});
          model.result("pgTemp").run();
          model.result().export("imgTemp").set("pngfilename", String.format(
              "comsol/v7_runtime/physical_fig7_like_%02.0fms.png", VELOCITIES[j]));
          model.result().export("imgTemp").run();
        }
        System.out.println("V7_IMAGE_EXPORT_PASS");
      } catch (Exception error) {
        System.out.println("V7_IMAGE_EXPORT_FAIL=" + error.getMessage());
      }
    }
    model.save("comsol/v7_runtime/physical_heat_" + meshName + ".mph");
    return values;
  }

  public static void main(String[] args) throws Exception {
    double[] fit = calibrate();
    double[][] coarse = solveFinal("coarse", 30, 12, 2, 4, fit[0], fit[1], false);
    double[][] medium = solveFinal("medium", 60, 24, 4, 8, fit[0], fit[1], true);
    double[][] fine = solveFinal("fine", 100, 36, 6, 12, fit[0], fit[1], false);
    boolean all = true;
    for (int j = 0; j < VELOCITIES.length; j++) {
      double tmaxRelative = Math.abs(medium[3][j] - fine[3][j]) / Math.abs(fine[3][j]);
      double regionDifference = Math.max(Math.abs(medium[0][j] - fine[0][j]),
          Math.max(Math.abs(medium[1][j] - fine[1][j]), Math.abs(medium[2][j] - fine[2][j])));
      boolean pass = tmaxRelative < 0.01 && regionDifference < 0.2
          && fine[16][j] < 0.005 && fine[19][j] < 1e-10;
      all &= pass;
      System.out.printf("V7_MESH_AUDIT velocity_m_s=%.0f Tmax_relative=%.12g max_region_difference_C=%.12g energy_error=%.12g mdot_identity_error=%.12g pass=%s%n",
          VELOCITIES[j], tmaxRelative, regionDifference, fine[16][j], fine[19][j], pass);
    }
    System.out.println("V7_NUMERICAL_ACCEPTANCE_PASS=" + all);
    if (!all) throw new IllegalStateException("V7 numerical mesh/energy gate failed");
  }
}
