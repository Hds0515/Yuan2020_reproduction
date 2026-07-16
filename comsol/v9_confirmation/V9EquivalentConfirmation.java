import com.comsol.model.*;
import com.comsol.model.util.*;

import java.util.Locale;

/**
 * Independent COMSOL 6.4 implementation of the frozen V9 equivalent model.
 *
 * <p>This is a conservative 1-D finite-element model, not a reconstruction of
 * the authors' CAD, channel dimensions, or manifold. The four inverse
 * parameters are copied verbatim from
 * outputs_v9/frozen_equivalent_parameters.json and are never fitted here.</p>
 */
public class V9EquivalentConfirmation {
  private static final double[] VELOCITIES = {4, 5, 6, 8, 10, 12};
  private static final int PROFILE_POINTS = 200;
  private static final String OUTPUT_DIR = "comsol/v9_confirmation";

  private static final String[] METRIC_NAMES = {
    "T_inlet_region_C",
    "T_middle_region_C",
    "T_outlet_region_C",
    "Tmax_C",
    "Tmin_C",
    "DeltaT_C",
    "hotspot_location_normalized",
    "air_outlet_temperature_C",
    "air_enthalpy_gain_W",
    "natural_convection_loss_W",
    "radiation_loss_W",
    "total_heat_input_W",
    "energy_residual_relative",
    "mdot_kg_s"
  };

  private static String[] metricExpressions() {
    String z1 = "intop(if(x<L/3,T-273.15[K],0[K]))/intop(if(x<L/3,1,0))";
    String z2 =
        "intop(if(x>=L/3&&x<2*L/3,T-273.15[K],0[K]))"
        + "/intop(if(x>=L/3&&x<2*L/3,1,0))";
    String z3 =
        "intop(if(x>=2*L/3,T-273.15[K],0[K]))"
        + "/intop(if(x>=2*L/3,1,0))";
    String qAir = "rhoAir*uin*AinTotal*CpAir*(outop(Ta)-Tin)";
    String qNatural = "intop(Pext*hNatural*(T-Tamb))";
    String qRadiation = "intop(Pext*emissivity*sigmaSB*(T^4-Tamb^4))";
    String qInput = "intop(qLine)";
    String residual =
        "abs((" + qInput + ")-(" + qAir + ")-(" + qNatural + ")-("
        + qRadiation + "))/abs(" + qInput + ")";
    return new String[] {
      z1,
      z2,
      z3,
      "maxop(T)-273.15[K]",
      "minop(T)-273.15[K]",
      "maxop(T)-minop(T)",
      "maxop(T,x)/L",
      "outop(Ta)-273.15[K]",
      qAir,
      qNatural,
      qRadiation,
      qInput,
      residual,
      "rhoAir*uin*AinTotal"
    };
  }

  private static Model buildModel(String meshName, int elements) {
    ModelUtil.clear();
    Model model = ModelUtil.create("Model");
    model.label("v9_equivalent_confirmation_" + meshName + ".mph");

    // Frozen V8 extensive scales and Table 2 properties.
    model.param().set("Tamb", "25[degC]");
    model.param().set("Tin", "25[degC]");
    model.param().set("L", "0.24936900555275113[m]");
    model.param().set("AinTotal", "3962[mm^2]");
    model.param().set("Aheat", "0.0247[m^2]");
    model.param().set("qTotal", "59.90985[W]");
    model.param().set("width", "Aheat/L");
    model.param().set("solidThickness", "10[mm]");
    model.param().set("Asolid", "width*solidThickness");
    model.param().set("Aexternal", "0.03166838011105502[m^2]");
    model.param().set("Pext", "Aexternal/L");
    model.param().set("kGraphite", "24[W/(m*K)]");
    model.param().set("rhoAir", "1.184[kg/m^3]");
    model.param().set("CpAir", "1007[J/(kg*K)]");
    model.param().set("kAir", "0.0251[W/(m*K)]");
    model.param().set("muAir", "1.849e-5[Pa*s]");
    model.param().set("emissivity", "0.82");
    model.param().set("sigmaSB", "5.670374419e-8[W/(m^2*K^4)]");
    model.param().set("gravity", "9.80665[m/s^2]");
    model.param().set("uin", "4[m/s]");

    // Frozen V9 inverse parameters. These values must not be tuned here.
    model.param().set("h0", "88.75329934335754[W/(m^2*K)]");
    model.param().set("velExponent", "0.6572949080566769");
    model.param().set("coolingBias", "1.3679347347371371");
    model.param().set("heatSkew", "-0.3118665087792232");

    model.component().create("comp1", true);
    model.component("comp1").geom().create("geom1", 1);
    model.component("comp1").geom("geom1").lengthUnit("m");
    model.component("comp1").geom("geom1").feature().create("interval1", "Interval");
    model.component("comp1").geom("geom1").feature("interval1")
        .set("p1", "0");
    model.component("comp1").geom("geom1").feature("interval1")
        .set("p2", "L");
    model.component("comp1").geom("geom1").run();

    model.component("comp1").selection().create("selDomain", "Explicit");
    model.component("comp1").selection("selDomain").geom("geom1", 1);
    model.component("comp1").selection("selDomain").all();
    model.component("comp1").selection().create("selInlet", "Explicit");
    model.component("comp1").selection("selInlet").geom("geom1", 0);
    model.component("comp1").selection("selInlet").set(new int[] {1});
    model.component("comp1").selection().create("selOutlet", "Explicit");
    model.component("comp1").selection("selOutlet").geom("geom1", 0);
    model.component("comp1").selection("selOutlet").set(new int[] {2});

    model.component("comp1").cpl().create("intop", "Integration");
    model.component("comp1").cpl("intop").selection().named("selDomain");
    model.component("comp1").cpl().create("maxop", "Maximum");
    model.component("comp1").cpl("maxop").selection().named("selDomain");
    model.component("comp1").cpl().create("minop", "Minimum");
    model.component("comp1").cpl("minop").selection().named("selDomain");
    model.component("comp1").cpl().create("outop", "Average");
    model.component("comp1").cpl("outop").selection().named("selOutlet");

    model.component("comp1").variable().create("var1");
    model.component("comp1").variable("var1").selection().named("selDomain");
    model.component("comp1").variable("var1").set("xi", "x/L");
    model.component("comp1").variable("var1").set(
        "hForced",
        "h0*(uin/(8[m/s]))^velExponent*exp(-coolingBias*(xi-0.5))");
    model.component("comp1").variable("var1").set("hP", "hForced*width");
    model.component("comp1").variable("var1").set(
        "qRaw", "1+heatSkew*(2*xi-1)");
    model.component("comp1").variable("var1").set(
        "qLine", "qTotal*qRaw/intop(qRaw)");
    model.component("comp1").variable("var1").set(
        "mdot", "rhoAir*uin*AinTotal");
    model.component("comp1").variable("var1").set("mcp", "mdot*CpAir");
    model.component("comp1").variable("var1").set("nuAir", "muAir/rhoAir");
    model.component("comp1").variable("var1").set(
        "alphaAir", "kAir/(rhoAir*CpAir)");
    model.component("comp1").variable("var1").set("PrAir", "nuAir/alphaAir");
    model.component("comp1").variable("var1").set("filmT", "(T+Tamb)/2");
    model.component("comp1").variable("var1").set(
        "RaL",
        "gravity*(1/filmT)*max(T-Tamb,1e-6[K])*L^3/(nuAir*alphaAir)");
    model.component("comp1").variable("var1").set(
        "NuNatural",
        "(0.825+0.387*RaL^(1/6)"
        + "/(1+(0.492/PrAir)^(9/16))^(8/27))^2");
    model.component("comp1").variable("var1").set(
        "hNatural", "NuNatural*kAir/L");
    model.component("comp1").variable("var1").set(
        "qExternalLine",
        "Pext*(hNatural*(T-Tamb)"
        + "+emissivity*sigmaSB*(T^4-Tamb^4))");

    model.component("comp1").physics().create(
        "solid", "GeneralFormPDE", "geom1", new String[] {"T"});
    model.component("comp1").physics("solid").prop("Units")
        .set("DependentVariableQuantity", "temperature");
    model.component("comp1").physics("solid").prop("Units")
        .setIndex("CustomSourceTermUnit", "W/m", 0, 0);
    model.component("comp1").physics("solid").feature("gfeq1")
        .setIndex("Ga", new String[] {"-kGraphite*Asolid*Tx"}, 0);
    model.component("comp1").physics("solid").feature("gfeq1")
        .setIndex("f", "qLine-hP*(T-Ta)-qExternalLine", 0);
    model.component("comp1").physics("solid").feature("init1")
        .set("T", "45[degC]");

    model.component("comp1").physics().create(
        "air", "GeneralFormPDE", "geom1", new String[] {"Ta"});
    model.component("comp1").physics("air").prop("Units")
        .set("DependentVariableQuantity", "temperature");
    model.component("comp1").physics("air").prop("Units")
        .setIndex("CustomSourceTermUnit", "W/m", 0, 0);
    model.component("comp1").physics("air").feature("gfeq1")
        .setIndex("Ga", new String[] {"0[W]"}, 0);
    model.component("comp1").physics("air").feature("gfeq1")
        .setIndex("f", "hP*(T-Ta)-mcp*Tax", 0);
    model.component("comp1").physics("air").feature("init1")
        .set("Ta", "Tin");
    model.component("comp1").physics("air").feature()
        .create("inletTemperature", "DirichletBoundary", 0);
    model.component("comp1").physics("air").feature("inletTemperature")
        .selection().named("selInlet");
    model.component("comp1").physics("air").feature("inletTemperature")
        .setIndex("r", "Tin", 0);

    model.component("comp1").mesh().create("mesh1");
    model.component("comp1").mesh("mesh1").feature().create("edge1", "Edge");
    model.component("comp1").mesh("mesh1").feature("edge1")
        .selection().named("selDomain");
    model.component("comp1").mesh("mesh1").feature("edge1")
        .create("distribution1", "Distribution");
    model.component("comp1").mesh("mesh1").feature("edge1")
        .feature("distribution1").set("numelem", elements);
    model.component("comp1").mesh("mesh1").run();

    model.study().create("std1");
    model.study("std1").feature().create("stat", "Stationary");
    model.study("std1").feature().create("param", "Parametric");
    model.study("std1").feature("param").set("pname", new String[] {"uin"});
    model.study("std1").feature("param").set(
        "plistarr", new String[] {"4 5 6 8 10 12"});
    model.study("std1").feature("param").set("punit", new String[] {"m/s"});
    return model;
  }

  private static double[][] evaluateMetrics(Model model, String meshName) {
    model.result().numerical().create("metrics", "EvalGlobal");
    model.result().numerical("metrics").set("expr", metricExpressions());
    double[][] result = model.result().numerical("metrics").getReal();
    model.result().table().create("metricsTable", "Table");
    model.result().numerical("metrics").set("table", "metricsTable");
    model.result().numerical("metrics").setResult();
    model.result().export().create("metricsExport", "Table");
    model.result().export("metricsExport").set("table", "metricsTable");
    model.result().export("metricsExport").set(
        "filename", OUTPUT_DIR + "/metrics_" + meshName + ".csv");
    model.result().export("metricsExport").run();
    return result;
  }

  private static void exportProfiles(Model model, String meshName) {
    model.result().export().create("profileExport", "Data");
    model.result().export("profileExport").set(
        "expr", new String[] {"T-273.15[K]", "Ta-273.15[K]"});
    model.result().export("profileExport").set("sort", "on");
    model.result().export("profileExport").set(
        "filename", OUTPUT_DIR + "/profiles_" + meshName + ".csv");
    model.result().export("profileExport").run();
  }

  private static void exportProfileImages(Model model, String meshName) {
    try {
      model.result().create("pgProfile", "PlotGroup1D");
      model.result("pgProfile").label("Equivalent solid and air temperatures");
      model.result("pgProfile").feature().create("solidLine", "LineGraph");
      model.result("pgProfile").feature("solidLine")
          .set("expr", "T-273.15[K]");
      model.result("pgProfile").feature("solidLine")
          .selection().named("selDomain");
      model.result("pgProfile").feature("solidLine")
          .set("descr", "Equivalent solid temperature");
      model.result("pgProfile").feature().create("airLine", "LineGraph");
      model.result("pgProfile").feature("airLine")
          .set("expr", "Ta-273.15[K]");
      model.result("pgProfile").feature("airLine")
          .selection().named("selDomain");
      model.result("pgProfile").feature("airLine")
          .set("descr", "Along-flow air temperature");
      model.result().export().create("profileImage", "pgProfile", "Image");
      model.result().export("profileImage").set("width", 1200);
      model.result().export("profileImage").set("height", 700);
      model.result().export("profileImage").set("zoomextents", "on");
      model.result("pgProfile").run();
      model.result().export("profileImage").set(
          "pngfilename",
          OUTPUT_DIR + "/temperature_profiles_" + meshName + "_all_speeds.png");
      model.result().export("profileImage").run();
      System.out.println("V9_COMSOL_IMAGE_EXPORT_PASS mesh=" + meshName);
    } catch (Exception error) {
      System.out.println(
          "V9_COMSOL_IMAGE_EXPORT_FAIL mesh=" + meshName
          + " message=" + error.getMessage());
    }
  }

  private static void solveMesh(
      String meshName, int elements, boolean exportImages) throws Exception {
    Model model = buildModel(meshName, elements);
    System.out.printf(
        Locale.US,
        "V9_COMSOL_SOLVE_BEGIN mesh=%s elements=%d%n",
        meshName,
        elements);
    long start = System.nanoTime();
    model.study("std1").run();
    double runtimeSeconds = (System.nanoTime() - start) / 1.0e9;
    System.out.printf(
        Locale.US,
        "V9_COMSOL_SOLVE_CONVERGED mesh=%s runtime_s=%.9g%n",
        meshName,
        runtimeSeconds);
    double[][] metrics = evaluateMetrics(model, meshName);
    exportProfiles(model, meshName);
    for (int speed = 0; speed < VELOCITIES.length; speed++) {
      System.out.printf(
          Locale.US,
          "V9_COMSOL_RESULT mesh=%s velocity_m_s=%.0f",
          meshName,
          VELOCITIES[speed]);
      for (int metric = 0; metric < METRIC_NAMES.length; metric++) {
        System.out.printf(
            Locale.US, " %s=%.12g", METRIC_NAMES[metric], metrics[metric][speed]);
      }
      System.out.println();
    }
    if (exportImages) {
      exportProfileImages(model, meshName);
    }
    model.save(OUTPUT_DIR + "/v9_equivalent_confirmation_" + meshName + ".mph");
  }

  public static void main(String[] args) throws Exception {
    solveMesh("coarse", 50, false);
    solveMesh("medium", 100, true);
    solveMesh("fine", 200, false);
    System.out.println("V9_COMSOL_CONFIRMATION_RUN_COMPLETE");
  }
}
