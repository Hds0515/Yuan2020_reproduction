import com.comsol.model.*;
import com.comsol.model.util.*;

import java.util.Locale;

/**
 * Independent COMSOL 6.4 implementation of the frozen V10 literature channel.
 *
 * <p>The geometry and six validation inputs are copied from Shahsavari et al.
 * (2012), Tables 1 and 2. No parameter is fitted here. This is a conservative
 * one-dimensional representative channel, not recovered full-cell CAD.</p>
 */
public class V10LiteratureChannel {
  private static final int CASES = 6;
  private static final String OUTPUT_DIR = "comsol/v10_reference";

  private static final String[] METRIC_NAMES = {
    "T_inlet_region_C",
    "T_middle_region_C",
    "T_outlet_region_C",
    "Tmax_C",
    "Tmin_C",
    "DeltaT_C",
    "hotspot_position_normalized",
    "air_outlet_temperature_C",
    "air_enthalpy_gain_channel_W",
    "channel_heat_input_W",
    "energy_residual_fraction",
    "mass_flow_channel_kg_s",
    "reynolds_number",
    "nusselt_number",
    "forced_h_W_m2K",
    "pressure_drop_Pa"
  };

  private static String nestedCaseExpression(String[] values) {
    String expression = values[values.length - 1];
    for (int index = values.length - 2; index >= 0; index--) {
      expression = "if(caseId<" + (index + 1.5) + "," + values[index]
          + "," + expression + ")";
    }
    return expression;
  }

  private static String[] metricExpressions() {
    String z1 = "intop(if(x<L/3,T-273.15[K],0[K]))/intop(if(x<L/3,1,0))";
    String z2 =
        "intop(if(x>=L/3&&x<2*L/3,T-273.15[K],0[K]))"
        + "/intop(if(x>=L/3&&x<2*L/3,1,0))";
    String z3 =
        "intop(if(x>=2*L/3,T-273.15[K],0[K]))"
        + "/intop(if(x>=2*L/3,1,0))";
    String qAir = "rhoAir*uin*ACh*CpAir*(outop(Ta)-Tin)";
    String qInput = "qCell/nChannels";
    String residual = "abs((" + qInput + ")-(" + qAir + "))/abs(" + qInput + ")";
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
      qInput,
      residual,
      "rhoAir*uin*ACh",
      "rhoAir*uin*Dh/muAir",
      "Nu",
      "hForced",
      "(56.91/(rhoAir*uin*Dh/muAir))*(L/Dh)*rhoAir*uin^2/2"
    };
  }

  private static Model buildModel(String meshName, int elements) {
    ModelUtil.clear();
    Model model = ModelUtil.create("Model");
    model.label("v10_literature_channel_" + meshName + ".mph");

    // Shahsavari et al. (2012), Tables 1 and 2.
    model.param().set("caseId", "1");
    model.param().set("Tin", "21[degC]");
    model.param().set("L", "60[mm]");
    model.param().set("cellWidth", "280[mm]");
    model.param().set("nChannels", "80");
    model.param().set("pitch", "cellWidth/nChannels");
    model.param().set("hCh", "2.5[mm]");
    model.param().set("wLong", "2.5[mm]");
    model.param().set("theta", "80[deg]");
    model.param().set("wShort", "wLong-2*hCh/tan(theta)");
    model.param().set("ACh", "(wLong+wShort)*hCh/2");
    model.param().set("PCh", "wLong+wShort+2*hCh/sin(theta)");
    model.param().set("Dh", "4*ACh/PCh");
    model.param().set("tBP", "5.5[mm]");
    model.param().set("tGDLa", "0.2[mm]");
    model.param().set("tGDLc", "0.2[mm]");
    model.param().set("tCCM", "0.05[mm]");
    model.param().set("kBPip", "60[W/(m*K)]");
    model.param().set("kGDLip", "10[W/(m*K)]");
    model.param().set("kCCM", "1.5[W/(m*K)]");
    model.param().set(
        "kAxial", "pitch*(kBPip*tBP+kGDLip*(tGDLa+tGDLc)+kCCM*tCCM)");

    // Fixed dry-air properties near the 21 degC inlet condition.
    model.param().set("rhoAir", "1.204[kg/m^3]");
    model.param().set("CpAir", "1006[J/(kg*K)]");
    model.param().set("kAir", "0.02514[W/(m*K)]");
    model.param().set("muAir", "1.825e-5[Pa*s]");
    model.param().set("Nu", "3.610224");
    model.param().set("hForced", "Nu*kAir/Dh");

    // Six experimental inputs from Table 1; caseId is only an index.
    model.param().set(
        "uin",
        nestedCaseExpression(new String[] {
          "1.76[m/s]", "2.19[m/s]", "1.63[m/s]", "1.45[m/s]", "0.95[m/s]", "0.94[m/s]"
        }));
    model.param().set(
        "qCell",
        nestedCaseExpression(new String[] {
          "55.7[W]", "37.7[W]", "27.6[W]", "29.1[W]", "6.3[W]", "13.7[W]"
        }));

    model.component().create("comp1", true);
    model.component("comp1").geom().create("geom1", 1);
    model.component("comp1").geom("geom1").lengthUnit("m");
    model.component("comp1").geom("geom1").feature().create("interval1", "Interval");
    model.component("comp1").geom("geom1").feature("interval1").set("p1", "0");
    model.component("comp1").geom("geom1").feature("interval1").set("p2", "L");
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
    model.component("comp1").variable("var1").set("qLine", "qCell/(nChannels*L)");
    model.component("comp1").variable("var1").set("hP", "hForced*PCh");
    model.component("comp1").variable("var1").set("mdot", "rhoAir*uin*ACh");
    model.component("comp1").variable("var1").set("mcp", "mdot*CpAir");
    model.component("comp1").variable("var1").set("Re", "rhoAir*uin*Dh/muAir");
    model.component("comp1").variable("var1").set("fDarcy", "56.91/Re");
    model.component("comp1").variable("var1").set(
        "deltaP", "fDarcy*(L/Dh)*rhoAir*uin^2/2");

    model.component("comp1").physics().create(
        "solid", "GeneralFormPDE", "geom1", new String[] {"T"});
    model.component("comp1").physics("solid").prop("Units")
        .set("DependentVariableQuantity", "temperature");
    model.component("comp1").physics("solid").prop("Units")
        .setIndex("CustomSourceTermUnit", "W/m", 0, 0);
    model.component("comp1").physics("solid").feature("gfeq1")
        .setIndex("Ga", new String[] {"-kAxial*Tx"}, 0);
    model.component("comp1").physics("solid").feature("gfeq1")
        .setIndex("f", "qLine-hP*(T-Ta)", 0);
    model.component("comp1").physics("solid").feature("init1").set("T", "50[degC]");

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
    model.component("comp1").physics("air").feature("init1").set("Ta", "Tin");
    model.component("comp1").physics("air").feature()
        .create("inletTemperature", "DirichletBoundary", 0);
    model.component("comp1").physics("air").feature("inletTemperature")
        .selection().named("selInlet");
    model.component("comp1").physics("air").feature("inletTemperature")
        .setIndex("r", "Tin", 0);

    model.component("comp1").mesh().create("mesh1");
    model.component("comp1").mesh("mesh1").feature().create("edge1", "Edge");
    model.component("comp1").mesh("mesh1").feature("edge1").selection().named("selDomain");
    model.component("comp1").mesh("mesh1").feature("edge1")
        .create("distribution1", "Distribution");
    model.component("comp1").mesh("mesh1").feature("edge1")
        .feature("distribution1").set("numelem", elements);
    model.component("comp1").mesh("mesh1").run();

    model.study().create("std1");
    model.study("std1").feature().create("stat", "Stationary");
    model.study("std1").feature().create("param", "Parametric");
    model.study("std1").feature("param").set("pname", new String[] {"caseId"});
    model.study("std1").feature("param").set("plistarr", new String[] {"1 2 3 4 5 6"});
    model.study("std1").feature("param").set("punit", new String[] {""});
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

  private static void solveMesh(String meshName, int elements) throws Exception {
    Model model = buildModel(meshName, elements);
    System.out.printf(Locale.US, "V10_COMSOL_SOLVE_BEGIN mesh=%s elements=%d%n", meshName, elements);
    long start = System.nanoTime();
    model.study("std1").run();
    double runtimeSeconds = (System.nanoTime() - start) / 1.0e9;
    System.out.printf(Locale.US, "V10_COMSOL_SOLVE_CONVERGED mesh=%s runtime_s=%.9g%n", meshName, runtimeSeconds);
    double[][] metrics = evaluateMetrics(model, meshName);
    exportProfiles(model, meshName);
    for (int caseIndex = 0; caseIndex < CASES; caseIndex++) {
      System.out.printf(Locale.US, "V10_COMSOL_RESULT mesh=%s case_id=%d", meshName, caseIndex + 1);
      for (int metric = 0; metric < METRIC_NAMES.length; metric++) {
        System.out.printf(Locale.US, " %s=%.12g", METRIC_NAMES[metric], metrics[metric][caseIndex]);
      }
      System.out.println();
    }
    model.save(OUTPUT_DIR + "/v10_literature_channel_" + meshName + ".mph");
  }

  public static void main(String[] args) throws Exception {
    solveMesh("coarse", 50);
    solveMesh("medium", 100);
    solveMesh("fine", 200);
    System.out.println("V10_COMSOL_REFERENCE_RUN_COMPLETE");
  }
}
