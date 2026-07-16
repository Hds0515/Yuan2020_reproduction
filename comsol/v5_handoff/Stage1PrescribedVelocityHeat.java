import com.comsol.model.*;
import com.comsol.model.util.*;

/**
 * Prescribed-velocity conjugate heat-transfer model for COMSOL 6.4.
 *
 * The Heat Transfer in Solids and Fluids interface is used without Laminar
 * Flow. The Fluid Heat Transfer Model receives the verified user-defined
 * velocity vector {uin,0,0}. Each mesh runs the fixed velocity set
 * 4, 5, 6, 8, 10, and 12 m/s with COMSOL's generated default solver.
 */
public class Stage1PrescribedVelocityHeat {
  private static final double[] VELOCITIES = {4, 5, 6, 8, 10, 12};

  private static void createBoxBoundary(
      Model model, String tag, String xmin, String xmax,
      String ymin, String ymax, String zmin, String zmax) {
    model.component("comp1").selection().create(tag, "Box");
    model.component("comp1").selection(tag).set("entitydim", 2);
    model.component("comp1").selection(tag).set("condition", "inside");
    model.component("comp1").selection(tag).set("xmin", xmin);
    model.component("comp1").selection(tag).set("xmax", xmax);
    model.component("comp1").selection(tag).set("ymin", ymin);
    model.component("comp1").selection(tag).set("ymax", ymax);
    model.component("comp1").selection(tag).set("zmin", zmin);
    model.component("comp1").selection(tag).set("zmax", zmax);
  }

  private static void createBoxEdge(
      Model model, String tag, String xmin, String xmax,
      String ymin, String ymax, String zmin, String zmax) {
    model.component("comp1").selection().create(tag, "Box");
    model.component("comp1").selection(tag).set("entitydim", 1);
    model.component("comp1").selection(tag).set("condition", "inside");
    model.component("comp1").selection(tag).set("xmin", xmin);
    model.component("comp1").selection(tag).set("xmax", xmax);
    model.component("comp1").selection(tag).set("ymin", ymin);
    model.component("comp1").selection(tag).set("ymax", ymax);
    model.component("comp1").selection(tag).set("zmin", zmin);
    model.component("comp1").selection(tag).set("zmax", zmax);
  }

  private static Model buildModel(
      String meshName, int nLength, int nWidth, int nAir, int nSolid) {
    ModelUtil.clear();
    Model model = ModelUtil.create("Model");
    model.label("prescribed_velocity_heat_" + meshName + ".mph");
    model.param().set("L", "0.30[m]");
    model.param().set("W", "0.082333333[m]");
    model.param().set("Hair", "1[mm]");
    model.param().set("Hs", "10[mm]");
    model.param().set("uin", "4[m/s]");
    model.param().set("Tin", "25[degC]");
    model.param().set("qflux", "2425.5[W/m^2]");
    model.param().set("rhoAir", "1.184[kg/m^3]");
    model.param().set("CpAir", "1007[J/(kg*K)]");
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

    createBoxBoundary(model, "selAirIn", "-tol", "tol", "-tol", "W+tol",
        "-tol", "Hair+tol");
    createBoxBoundary(model, "selAirOut", "L-tol", "L+tol", "-tol", "W+tol",
        "-tol", "Hair+tol");
    createBoxBoundary(model, "selSource", "-tol", "tol", "-tol", "W+tol",
        "-tol", "Hair+Hs+tol");
    createBoxBoundary(model, "selTarget", "L-tol", "L+tol", "-tol", "W+tol",
        "-tol", "Hair+Hs+tol");
    createBoxBoundary(model, "selHeat", "-tol", "L+tol", "-tol", "W+tol",
        "Hair+Hs-tol", "Hair+Hs+tol");
    model.component("comp1").selection().create("selAllDomains", "Explicit");
    model.component("comp1").selection("selAllDomains").geom("geom1", 3);
    model.component("comp1").selection("selAllDomains").all();

    // The inlet cross section is meshed entirely with mapped quadrilaterals.
    // Explicit edge counts make coarse/medium/fine a genuine nested refinement
    // instead of three unrelated free-triangle topologies.
    createBoxEdge(model, "selWidthBottom", "-tol", "tol", "-tol", "W+tol",
        "-tol", "tol");
    createBoxEdge(model, "selWidthInterface", "-tol", "tol", "-tol", "W+tol",
        "Hair-tol", "Hair+tol");
    createBoxEdge(model, "selWidthTop", "-tol", "tol", "-tol", "W+tol",
        "Hair+Hs-tol", "Hair+Hs+tol");
    model.component("comp1").selection().create("selWidthEdges", "Union");
    model.component("comp1").selection("selWidthEdges").set("entitydim", 1);
    model.component("comp1").selection("selWidthEdges").set("input",
        new String[]{"selWidthBottom", "selWidthInterface", "selWidthTop"});

    createBoxEdge(model, "selAirThicknessLeft", "-tol", "tol", "-tol", "tol",
        "-tol", "Hair+tol");
    createBoxEdge(model, "selAirThicknessRight", "-tol", "tol", "W-tol", "W+tol",
        "-tol", "Hair+tol");
    model.component("comp1").selection().create("selAirThicknessEdges", "Union");
    model.component("comp1").selection("selAirThicknessEdges").set("entitydim", 1);
    model.component("comp1").selection("selAirThicknessEdges").set("input",
        new String[]{"selAirThicknessLeft", "selAirThicknessRight"});

    createBoxEdge(model, "selSolidThicknessLeft", "-tol", "tol", "-tol", "tol",
        "Hair-tol", "Hair+Hs+tol");
    createBoxEdge(model, "selSolidThicknessRight", "-tol", "tol", "W-tol", "W+tol",
        "Hair-tol", "Hair+Hs+tol");
    model.component("comp1").selection().create("selSolidThicknessEdges", "Union");
    model.component("comp1").selection("selSolidThicknessEdges").set("entitydim", 1);
    model.component("comp1").selection("selSolidThicknessEdges").set("input",
        new String[]{"selSolidThicknessLeft", "selSolidThicknessRight"});

    if (model.component("comp1").selection("selAirIn").entities().length != 1
        || model.component("comp1").selection("selAirOut").entities().length != 1
        || model.component("comp1").selection("selHeat").entities().length != 1) {
      throw new IllegalStateException("Thermal boundary selection audit failed.");
    }

    model.component("comp1").material().create("matAir", "Common");
    model.component("comp1").material("matAir").selection().named("geom1_air_dom");
    model.component("comp1").material("matAir").propertyGroup("def")
         .set("density", "rhoAir");
    model.component("comp1").material("matAir").propertyGroup("def")
         .set("heatcapacity", "CpAir");
    model.component("comp1").material("matAir").propertyGroup("def")
         .set("thermalconductivity", "0.0257[W/(m*K)]");

    model.component("comp1").material().create("matGraphite", "Common");
    model.component("comp1").material("matGraphite").selection()
         .named("geom1_solid_dom");
    model.component("comp1").material("matGraphite").propertyGroup("def")
         .set("density", "1800[kg/m^3]");
    model.component("comp1").material("matGraphite").propertyGroup("def")
         .set("heatcapacity", "710[J/(kg*K)]");
    model.component("comp1").material("matGraphite").propertyGroup("def")
         .set("thermalconductivity", "150[W/(m*K)]");

    model.component("comp1").physics().create(
        "ht", "HeatTransferInSolidsAndFluids", "geom1");
    // solid1 is the locked default domain model; assigning the air domain to
    // fluid1 automatically leaves the graphite domain on solid1.
    model.component("comp1").physics("ht").feature("fluid1").selection()
         .named("geom1_air_dom");
    model.component("comp1").physics("ht").feature("fluid1")
         .set("u_src", "userdef");
    model.component("comp1").physics("ht").feature("fluid1")
         .set("u", new String[]{"uin", "0", "0"});
    model.component("comp1").physics("ht").feature().create("in1", "Inflow", 2);
    model.component("comp1").physics("ht").feature("in1").selection()
         .named("selAirIn");
    model.component("comp1").physics("ht").feature("in1").set("Tustr", "Tin");
    // The default thermal-insulation condition on the outlet is the zero
    // diffusive-flux condition required by the prescribed-velocity outflow.
    model.component("comp1").physics("ht").feature().create(
        "hf1", "HeatFluxBoundary", 2);
    model.component("comp1").physics("ht").feature("hf1").selection()
         .named("selHeat");
    model.component("comp1").physics("ht").feature("hf1").set("q0", "qflux");
    model.component("comp1").physics("ht").feature("init1").set("Tinit", "Tin");

    model.component("comp1").mesh().create("mesh1");
    model.component("comp1").mesh("mesh1").feature().create("edgWidth", "Edge");
    model.component("comp1").mesh("mesh1").feature("edgWidth").selection()
         .named("selWidthEdges");
    model.component("comp1").mesh("mesh1").feature("edgWidth")
         .create("dis1", "Distribution");
    model.component("comp1").mesh("mesh1").feature("edgWidth").feature("dis1")
         .set("numelem", nWidth);
    model.component("comp1").mesh("mesh1").feature().create("edgAir", "Edge");
    model.component("comp1").mesh("mesh1").feature("edgAir").selection()
         .named("selAirThicknessEdges");
    model.component("comp1").mesh("mesh1").feature("edgAir")
         .create("dis1", "Distribution");
    model.component("comp1").mesh("mesh1").feature("edgAir").feature("dis1")
         .set("numelem", nAir);
    model.component("comp1").mesh("mesh1").feature().create("edgSolid", "Edge");
    model.component("comp1").mesh("mesh1").feature("edgSolid").selection()
         .named("selSolidThicknessEdges");
    model.component("comp1").mesh("mesh1").feature("edgSolid")
         .create("dis1", "Distribution");
    model.component("comp1").mesh("mesh1").feature("edgSolid").feature("dis1")
         .set("numelem", nSolid);
    model.component("comp1").mesh("mesh1").feature().create("map1", "Map");
    model.component("comp1").mesh("mesh1").feature("map1").selection()
         .named("selSource");
    model.component("comp1").mesh("mesh1").feature().create("swe1", "Sweep");
    model.component("comp1").mesh("mesh1").feature("swe1").selection()
         .named("selAllDomains");
    model.component("comp1").mesh("mesh1").feature("swe1").selection("sourceface")
         .named("selSource");
    model.component("comp1").mesh("mesh1").feature("swe1").selection("targetface")
         .named("selTarget");
    model.component("comp1").mesh("mesh1").feature("swe1")
         .create("dis1", "Distribution");
    model.component("comp1").mesh("mesh1").feature("swe1").feature("dis1")
         .set("numelem", nLength);
    model.component("comp1").mesh("mesh1").run();

    model.component("comp1").cpl().create("intSolid", "Integration");
    model.component("comp1").cpl("intSolid").selection().named("geom1_solid_dom");
    model.component("comp1").cpl().create("intHeat", "Integration");
    model.component("comp1").cpl("intHeat").selection().named("selHeat");
    model.component("comp1").cpl().create("aveIn", "Average");
    model.component("comp1").cpl("aveIn").selection().named("selAirIn");
    model.component("comp1").cpl().create("aveOut", "Average");
    model.component("comp1").cpl("aveOut").selection().named("selAirOut");
    model.component("comp1").cpl().create("maxSolid", "Maximum");
    model.component("comp1").cpl("maxSolid").selection().named("geom1_solid_dom");

    model.study().create("std1");
    model.study("std1").feature().create("stat", "Stationary");
    model.study("std1").feature().create("param", "Parametric");
    model.study("std1").feature("param").set("pname", new String[]{"uin"});
    model.study("std1").feature("param").set(
        "plistarr", new String[]{"4 5 6 8 10 12"});
    model.study("std1").feature("param").set("punit", new String[]{"m/s"});
    return model;
  }

  private static double[][] solveAndEvaluate(
      String meshName, int nLength, int nWidth, int nAir, int nSolid,
      boolean exportImages) throws Exception {
    Model model = buildModel(meshName, nLength, nWidth, nAir, nSolid);
    System.out.println("HEAT_SOLVE_BEGIN mesh=" + meshName);
    model.study("std1").run();
    System.out.println("HEAT_SOLVE_CONVERGED mesh=" + meshName);

    String z1 = "intSolid(if(x<L/3,T-273.15[K],0[K]))/intSolid(if(x<L/3,1,0))";
    String z2 = "intSolid(if(x>=L/3&&x<2*L/3,T-273.15[K],0[K]))"
        + "/intSolid(if(x>=L/3&&x<2*L/3,1,0))";
    String z3 = "intSolid(if(x>=2*L/3,T-273.15[K],0[K]))"
        + "/intSolid(if(x>=2*L/3,1,0))";
    String qIn = "intHeat(qflux)";
    String qAir = "rhoAir*CpAir*uin*W*Hair*(aveOut(T)-Tin)";
    String[] expressions = {
      z1, z2, z3,
      "maxSolid(T)-273.15[K]", "maxSolid(T,x)", "maxSolid(T,y)", "maxSolid(T,z)",
      "Tin-273.15[K]", "aveIn(T)-273.15[K]", "aveOut(T)-273.15[K]",
      qIn, qAir, "abs((" + qIn + ")-(" + qAir + "))/abs(" + qIn + ")"
    };
    String[] names = {
      "T_inlet_region_C", "T_middle_region_C", "T_outlet_region_C",
      "Tmax_C", "hotspot_x_m", "hotspot_y_m", "hotspot_z_m",
      "air_inflow_upstream_temperature_C", "air_inlet_boundary_temperature_C",
      "air_outlet_temperature_C",
      "heat_input_W", "air_enthalpy_gain_W", "relative_energy_error"
    };
    model.result().numerical().create("gevHeat", "EvalGlobal");
    model.result().numerical("gevHeat").set("expr", expressions);
    double[][] values = model.result().numerical("gevHeat").getReal();
    for (int j = 0; j < VELOCITIES.length; j++) {
      System.out.printf("HEAT_RESULT mesh=%s velocity_m_s=%.0f", meshName, VELOCITIES[j]);
      for (int i = 0; i < names.length; i++) {
        System.out.printf(" %s=%.12g", names[i], values[i][j]);
      }
      System.out.println();
    }

    model.result().table().create("tblHeat", "Table");
    model.result().numerical("gevHeat").set("table", "tblHeat");
    model.result().numerical("gevHeat").setResult();
    model.result().export().create("tblHeatExport", "Table");
    model.result().export("tblHeatExport").set("table", "tblHeat");
    model.result().export("tblHeatExport").set("filename",
        "comsol/v6_runtime/heat_metrics_" + meshName + ".csv");
    model.result().export("tblHeatExport").run();

    if (exportImages) {
      try {
        model.result().create("pgTemp", "PlotGroup3D");
        model.result("pgTemp").feature().create("surf1", "Surface");
        model.result("pgTemp").feature("surf1").set("expr", "T-273.15[K]");
        model.result("pgTemp").feature("surf1").set("colortable", "ThermalLight");
        model.result().export().create("imgTemp", "pgTemp", "Image");
        model.result().export("imgTemp").set("size", "manualweb");
        model.result().export("imgTemp").set("width", 1200);
        model.result().export("imgTemp").set("height", 700);
        model.result().export("imgTemp").set("zoomextents", "on");
        for (int j = 0; j < VELOCITIES.length; j++) {
          model.result("pgTemp").set("looplevel", new int[]{j + 1});
          model.result("pgTemp").run();
          model.result().export("imgTemp").set("pngfilename", String.format(
              "comsol/v6_runtime/fig7_like_%02.0fms.png", VELOCITIES[j]));
          model.result().export("imgTemp").run();
        }
        System.out.println("HEAT_IMAGE_EXPORT_PASS");
      } catch (Exception imageError) {
        System.out.println("HEAT_IMAGE_EXPORT_FAIL=" + imageError.getMessage());
        imageError.printStackTrace(System.err);
      }
    }
    model.save("comsol/v6_runtime/prescribed_velocity_heat_" + meshName + ".mph");
    return values;
  }

  public static void main(String[] args) throws Exception {
    try {
      double[][] coarse = solveAndEvaluate("coarse", 30, 12, 2, 4, false);
      double[][] medium = solveAndEvaluate("medium", 60, 24, 4, 8, true);
      double[][] fine = solveAndEvaluate("fine", 100, 36, 6, 12, false);
      boolean passed = true;
      for (int j = 0; j < VELOCITIES.length; j++) {
        double tmaxRel = Math.abs(medium[3][j] - fine[3][j]) / Math.abs(fine[3][j]);
        double maxRegionDiff = Math.max(Math.abs(medium[0][j] - fine[0][j]),
            Math.max(Math.abs(medium[1][j] - fine[1][j]),
                Math.abs(medium[2][j] - fine[2][j])));
        double energyError = fine[12][j];
        boolean velocityPass = tmaxRel < 0.01 && maxRegionDiff < 0.2 && energyError < 0.005;
        passed = passed && velocityPass;
        System.out.printf(
            "HEAT_MESH_AUDIT velocity_m_s=%.0f medium_fine_Tmax_relative=%.12g"
            + " medium_fine_max_region_difference_C=%.12g fine_energy_error=%.12g pass=%s%n",
            VELOCITIES[j], tmaxRel, maxRegionDiff, energyError, velocityPass);
      }
      System.out.println("HEAT_ALL_MESHES_CONVERGED=true");
      System.out.println("HEAT_ACCEPTANCE_PASS=" + passed);
      if (!passed) {
        throw new IllegalStateException("Thermal mesh or energy acceptance gate failed.");
      }
    } catch (Exception error) {
      error.printStackTrace(System.err);
      throw error;
    }
  }
}
