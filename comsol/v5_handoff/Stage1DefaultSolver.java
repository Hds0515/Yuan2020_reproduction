import com.comsol.model.*;
import com.comsol.model.util.*;

/** Minimal 0.1 m/s laminar-flow audit using COMSOL's generated default solver. */
public class Stage1DefaultSolver {
  public static void main(String[] args) throws Exception {
    ModelUtil.clear();
    Model model = ModelUtil.create("Model");
    model.label("stage1_default_solver_0p1ms.mph");

    model.param().set("L", "0.30[m]");
    model.param().set("W", "0.082333333[m]");
    model.param().set("Hair", "1[mm]");
    model.param().set("Ain", "3962[mm^2]");
    model.param().set("areaScale", "Ain/(W*Hair)");
    model.param().set("uin", "0.1[m/s]");
    model.param().set("tol", "1e-7[m]");

    model.component().create("comp1", true);
    model.component("comp1").geom().create("geom1", 3);
    model.component("comp1").geom("geom1").lengthUnit("m");
    model.component("comp1").geom("geom1").feature().create("air", "Block");
    model.component("comp1").geom("geom1").feature("air")
         .set("size", new String[]{"L", "W", "Hair"});
    model.component("comp1").geom("geom1").feature("air").set("selresult", "on");
    model.component("comp1").geom("geom1").run();

    String[] tags = {"selInlet", "selOutlet"};
    String[] xmins = {"-tol", "L-tol"};
    String[] xmaxs = {"tol", "L+tol"};
    for (int i = 0; i < tags.length; i++) {
      model.component("comp1").selection().create(tags[i], "Box");
      model.component("comp1").selection(tags[i]).set("entitydim", 2);
      model.component("comp1").selection(tags[i]).set("condition", "inside");
      model.component("comp1").selection(tags[i]).set("xmin", xmins[i]);
      model.component("comp1").selection(tags[i]).set("xmax", xmaxs[i]);
      model.component("comp1").selection(tags[i]).set("ymin", "-tol");
      model.component("comp1").selection(tags[i]).set("ymax", "W+tol");
      model.component("comp1").selection(tags[i]).set("zmin", "-tol");
      model.component("comp1").selection(tags[i]).set("zmax", "Hair+tol");
    }
    int[] inlet = model.component("comp1").selection("selInlet").entities();
    int[] outlet = model.component("comp1").selection("selOutlet").entities();
    if (inlet.length != 1 || outlet.length != 1) {
      throw new IllegalStateException("Selection audit failed before the flow solve.");
    }

    model.component("comp1").selection().create("selEdgeZ0", "Box");
    model.component("comp1").selection("selEdgeZ0").set("entitydim", 1);
    model.component("comp1").selection("selEdgeZ0").set("condition", "inside");
    model.component("comp1").selection("selEdgeZ0").set("xmin", "-tol");
    model.component("comp1").selection("selEdgeZ0").set("xmax", "tol");
    model.component("comp1").selection("selEdgeZ0").set("ymin", "-tol");
    model.component("comp1").selection("selEdgeZ0").set("ymax", "tol");
    model.component("comp1").selection("selEdgeZ0").set("zmin", "-tol");
    model.component("comp1").selection("selEdgeZ0").set("zmax", "Hair+tol");
    model.component("comp1").selection().create("selEdgeZW", "Box");
    model.component("comp1").selection("selEdgeZW").set("entitydim", 1);
    model.component("comp1").selection("selEdgeZW").set("condition", "inside");
    model.component("comp1").selection("selEdgeZW").set("xmin", "-tol");
    model.component("comp1").selection("selEdgeZW").set("xmax", "tol");
    model.component("comp1").selection("selEdgeZW").set("ymin", "W-tol");
    model.component("comp1").selection("selEdgeZW").set("ymax", "W+tol");
    model.component("comp1").selection("selEdgeZW").set("zmin", "-tol");
    model.component("comp1").selection("selEdgeZW").set("zmax", "Hair+tol");
    model.component("comp1").selection().create("selThicknessEdges", "Union");
    model.component("comp1").selection("selThicknessEdges").set("entitydim", 1);
    model.component("comp1").selection("selThicknessEdges")
         .set("input", new String[]{"selEdgeZ0", "selEdgeZW"});

    model.component("comp1").selection().create("selWidthEdge0", "Box");
    model.component("comp1").selection("selWidthEdge0").set("entitydim", 1);
    model.component("comp1").selection("selWidthEdge0").set("condition", "inside");
    model.component("comp1").selection("selWidthEdge0").set("xmin", "-tol");
    model.component("comp1").selection("selWidthEdge0").set("xmax", "tol");
    model.component("comp1").selection("selWidthEdge0").set("ymin", "-tol");
    model.component("comp1").selection("selWidthEdge0").set("ymax", "W+tol");
    model.component("comp1").selection("selWidthEdge0").set("zmin", "-tol");
    model.component("comp1").selection("selWidthEdge0").set("zmax", "tol");
    model.component("comp1").selection().create("selWidthEdgeH", "Box");
    model.component("comp1").selection("selWidthEdgeH").set("entitydim", 1);
    model.component("comp1").selection("selWidthEdgeH").set("condition", "inside");
    model.component("comp1").selection("selWidthEdgeH").set("xmin", "-tol");
    model.component("comp1").selection("selWidthEdgeH").set("xmax", "tol");
    model.component("comp1").selection("selWidthEdgeH").set("ymin", "-tol");
    model.component("comp1").selection("selWidthEdgeH").set("ymax", "W+tol");
    model.component("comp1").selection("selWidthEdgeH").set("zmin", "Hair-tol");
    model.component("comp1").selection("selWidthEdgeH").set("zmax", "Hair+tol");
    model.component("comp1").selection().create("selWidthEdges", "Union");
    model.component("comp1").selection("selWidthEdges").set("entitydim", 1);
    model.component("comp1").selection("selWidthEdges")
         .set("input", new String[]{"selWidthEdge0", "selWidthEdgeH"});

    model.component("comp1").material().create("airmat", "Common");
    model.component("comp1").material("airmat").selection().named("geom1_air_dom");
    model.component("comp1").material("airmat").propertyGroup("def")
         .set("density", "1.184[kg/m^3]");
    model.component("comp1").material("airmat").propertyGroup("def")
         .set("dynamicviscosity", "1.849e-5[Pa*s]");

    model.component("comp1").physics().create("spf", "LaminarFlow", "geom1");
    model.component("comp1").physics("spf").selection().named("geom1_air_dom");
    // COMSOL 6.4 ShapeProperty uses order_fluid; value 2 selects the
    // quadratic-velocity/linear-pressure pair without altering the solver tree.
    model.component("comp1").physics("spf").prop("ShapeProperty")
         .set("order_fluid", 2);
    model.component("comp1").physics("spf").feature().create("inl1", "Inlet", 2);
    model.component("comp1").physics("spf").feature("inl1").selection().named("selInlet");
    model.component("comp1").physics("spf").feature("inl1")
         .set("BoundaryCondition", "Velocity");
    model.component("comp1").physics("spf").feature("inl1").set("U0in", "uin");
    model.component("comp1").physics("spf").feature().create("out1", "Outlet", 2);
    model.component("comp1").physics("spf").feature("out1").selection().named("selOutlet");
    model.component("comp1").physics("spf").feature("out1")
         .set("BoundaryCondition", "Pressure");
    model.component("comp1").physics("spf").feature("out1").set("p0", "0[Pa]");

    model.component("comp1").mesh().create("mesh1");
    model.component("comp1").mesh("mesh1").feature().create("edg1", "Edge");
    model.component("comp1").mesh("mesh1").feature("edg1").selection()
         .named("selThicknessEdges");
    model.component("comp1").mesh("mesh1").feature("edg1")
         .create("dis1", "Distribution");
    model.component("comp1").mesh("mesh1").feature("edg1").feature("dis1")
         .set("numelem", 4);
    model.component("comp1").mesh("mesh1").feature().create("edg2", "Edge");
    model.component("comp1").mesh("mesh1").feature("edg2").selection()
         .named("selWidthEdges");
    model.component("comp1").mesh("mesh1").feature("edg2")
         .create("dis1", "Distribution");
    model.component("comp1").mesh("mesh1").feature("edg2").feature("dis1")
         .set("numelem", 16);
    model.component("comp1").mesh("mesh1").feature().create("map1", "Map");
    model.component("comp1").mesh("mesh1").feature("map1").selection()
         .named("selInlet");
    model.component("comp1").mesh("mesh1").feature().create("swe1", "Sweep");
    model.component("comp1").mesh("mesh1").feature("swe1").selection()
         .named("geom1_air_dom");
    model.component("comp1").mesh("mesh1").feature("swe1").selection("sourceface")
         .named("selInlet");
    model.component("comp1").mesh("mesh1").feature("swe1").selection("targetface")
         .named("selOutlet");
    model.component("comp1").mesh("mesh1").feature("swe1")
         .create("dis1", "Distribution");
    model.component("comp1").mesh("mesh1").feature("swe1").feature("dis1")
         .set("numelem", 60);
    model.component("comp1").mesh("mesh1").run();

    model.component("comp1").cpl().create("intIn", "Integration");
    model.component("comp1").cpl("intIn").selection().named("selInlet");
    model.component("comp1").cpl().create("intOut", "Integration");
    model.component("comp1").cpl("intOut").selection().named("selOutlet");
    model.component("comp1").cpl().create("maxAir", "Maximum");
    model.component("comp1").cpl("maxAir").selection().named("geom1_air_dom");
    model.component("comp1").cpl().create("minAir", "Minimum");
    model.component("comp1").cpl("minAir").selection().named("geom1_air_dom");

    model.study().create("std1");
    model.study("std1").feature().create("stat", "Stationary");
    System.out.println("STAGE1_DEFAULT_SOLVER_BEGIN");
    model.study("std1").run();
    System.out.println("STAGE1_DEFAULT_SOLVER_CONVERGED");

    String flux = "spf.rho*(u*nx+v*ny+w*nz)";
    String[] expressions = {
      "abs(intIn(" + flux + "))",
      "abs(intOut(" + flux + "))",
      "areaScale*abs(intIn(" + flux + "))",
      "areaScale*abs(intOut(" + flux + "))",
      "abs(abs(intIn(" + flux + "))-abs(intOut(" + flux + ")))/abs(intIn(" + flux + "))",
      "minAir(p)",
      "maxAir(p)",
      "maxAir(p)-minAir(p)",
      "maxAir(sqrt(u^2+v^2+w^2))"
    };
    model.result().numerical().create("gevStage1", "EvalGlobal");
    model.result().numerical("gevStage1").set("expr", expressions);
    double[][] values = model.result().numerical("gevStage1").getReal();
    String[] names = {
      "channel_inlet_mass_flow_kg_s", "channel_outlet_mass_flow_kg_s",
      "scaled_inlet_mass_flow_kg_s", "scaled_outlet_mass_flow_kg_s",
      "relative_mass_imbalance", "minimum_pressure_Pa", "maximum_pressure_Pa",
      "pressure_range_Pa", "maximum_speed_m_s"
    };
    for (int i = 0; i < names.length; i++) {
      System.out.printf("%s=%.12g%n", names[i], values[i][0]);
    }
    if (!Double.isFinite(values[4][0]) || values[4][0] >= 0.001) {
      throw new IllegalStateException("Mass imbalance must be below 0.1 percent.");
    }
    model.save("comsol/v6_runtime/stage1_default_solver_0p1ms.mph");
    System.out.println("STAGE1_DEFAULT_SOLVER_PASS");
  }
}
