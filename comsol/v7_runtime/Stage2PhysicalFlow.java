import com.comsol.model.*;
import com.comsol.model.util.*;

/** V7 default-solver flow audit using the same representative-cell scaling as heat. */
public class Stage2PhysicalFlow {
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

  public static void main(String[] args) throws Exception {
    ModelUtil.clear();
    Model model = ModelUtil.create("Model");
    model.label("v7_physical_flow_0p1ms.mph");
    model.param().set("AinStack", "3962[mm^2]");
    model.param().set("nCells", "40");
    model.param().set("AinCell", "AinStack/nCells");
    model.param().set("Hair", "1[mm]");
    model.param().set("W", "AinCell/Hair");
    model.param().set("AheatCell", "0.0247[m^2]");
    model.param().set("L", "AheatCell/W");
    model.param().set("uin", "0.1[m/s]");
    model.param().set("rhoAir", "1.184[kg/m^3]");
    model.param().set("tol", "1e-7[m]");

    model.component().create("comp1", true);
    model.component("comp1").geom().create("geom1", 3);
    model.component("comp1").geom("geom1").lengthUnit("m");
    model.component("comp1").geom("geom1").feature().create("air", "Block");
    model.component("comp1").geom("geom1").feature("air")
        .set("size", new String[]{"L", "W", "Hair"});
    model.component("comp1").geom("geom1").feature("air").set("selresult", "on");
    model.component("comp1").geom("geom1").run();

    box(model, "selInlet", 2, "-tol", "tol", "-tol", "W+tol", "-tol", "Hair+tol");
    box(model, "selOutlet", 2, "L-tol", "L+tol", "-tol", "W+tol", "-tol", "Hair+tol");
    box(model, "selZ0", 1, "-tol", "tol", "-tol", "tol", "-tol", "Hair+tol");
    box(model, "selZW", 1, "-tol", "tol", "W-tol", "W+tol", "-tol", "Hair+tol");
    model.component("comp1").selection().create("selThickness", "Union");
    model.component("comp1").selection("selThickness").set("entitydim", 1);
    model.component("comp1").selection("selThickness").set("input", new String[]{"selZ0", "selZW"});
    box(model, "selWidth0", 1, "-tol", "tol", "-tol", "W+tol", "-tol", "tol");
    box(model, "selWidthH", 1, "-tol", "tol", "-tol", "W+tol", "Hair-tol", "Hair+tol");
    model.component("comp1").selection().create("selWidth", "Union");
    model.component("comp1").selection("selWidth").set("entitydim", 1);
    model.component("comp1").selection("selWidth").set("input", new String[]{"selWidth0", "selWidthH"});
    if (model.component("comp1").selection("selInlet").entities().length != 1
        || model.component("comp1").selection("selOutlet").entities().length != 1) {
      throw new IllegalStateException("V7 flow inlet/outlet selection must each contain one boundary");
    }

    model.component("comp1").material().create("matAir", "Common");
    model.component("comp1").material("matAir").selection().all();
    model.component("comp1").material("matAir").propertyGroup("def").set("density", "rhoAir");
    model.component("comp1").material("matAir").propertyGroup("def")
        .set("dynamicviscosity", "1.849e-5[Pa*s]");
    model.component("comp1").material("matAir").propertyGroup("def")
        .set("thermalconductivity", "0.0251[W/(m*K)]");
    model.component("comp1").material("matAir").propertyGroup("def")
        .set("heatcapacity", "1007[J/(kg*K)]");

    model.component("comp1").physics().create("spf", "LaminarFlow", "geom1");
    model.component("comp1").physics("spf").prop("ShapeProperty").set("order_fluid", 2);
    model.component("comp1").physics("spf").feature("fp1").set("rho_mat", "userdef");
    model.component("comp1").physics("spf").feature("fp1").set("rho", "rhoAir");
    model.component("comp1").physics("spf").feature("fp1").set("mu_mat", "userdef");
    model.component("comp1").physics("spf").feature("fp1").set("mu", "1.849e-5[Pa*s]");
    model.component("comp1").physics("spf").feature().create("inl1", "Inlet", 2);
    model.component("comp1").physics("spf").feature("inl1").selection().named("selInlet");
    model.component("comp1").physics("spf").feature("inl1").set("BoundaryCondition", "MassFlow");
    model.component("comp1").physics("spf").feature("inl1").set("MassFlowType", "MassFlowRate");
    model.component("comp1").physics("spf").feature("inl1").set("mfr", "rhoAir*uin*AinCell");
    model.component("comp1").physics("spf").feature().create("out1", "Outlet", 2);
    model.component("comp1").physics("spf").feature("out1").selection().named("selOutlet");
    model.component("comp1").physics("spf").feature("out1").set("BoundaryCondition", "Pressure");
    model.component("comp1").physics("spf").feature("out1").set("p0", "0[Pa]");

    model.component("comp1").mesh().create("mesh1");
    model.component("comp1").mesh("mesh1").feature().create("edgZ", "Edge");
    model.component("comp1").mesh("mesh1").feature("edgZ").selection().named("selThickness");
    model.component("comp1").mesh("mesh1").feature("edgZ").create("dis1", "Distribution");
    model.component("comp1").mesh("mesh1").feature("edgZ").feature("dis1").set("numelem", 4);
    model.component("comp1").mesh("mesh1").feature().create("edgW", "Edge");
    model.component("comp1").mesh("mesh1").feature("edgW").selection().named("selWidth");
    model.component("comp1").mesh("mesh1").feature("edgW").create("dis1", "Distribution");
    model.component("comp1").mesh("mesh1").feature("edgW").feature("dis1").set("numelem", 16);
    model.component("comp1").mesh("mesh1").feature().create("map1", "Map");
    model.component("comp1").mesh("mesh1").feature("map1").selection().named("selInlet");
    model.component("comp1").mesh("mesh1").feature().create("swe1", "Sweep");
    model.component("comp1").mesh("mesh1").feature("swe1").selection().named("geom1_air_dom");
    model.component("comp1").mesh("mesh1").feature("swe1").selection("sourceface").named("selInlet");
    model.component("comp1").mesh("mesh1").feature("swe1").selection("targetface").named("selOutlet");
    model.component("comp1").mesh("mesh1").feature("swe1").create("dis1", "Distribution");
    model.component("comp1").mesh("mesh1").feature("swe1").feature("dis1").set("numelem", 60);
    model.component("comp1").mesh("mesh1").run();

    model.component("comp1").cpl().create("intIn", "Integration");
    model.component("comp1").cpl("intIn").selection().named("selInlet");
    model.component("comp1").cpl().create("intOut", "Integration");
    model.component("comp1").cpl("intOut").selection().named("selOutlet");
    model.component("comp1").cpl().create("minAir", "Minimum");
    model.component("comp1").cpl("minAir").selection().all();
    model.component("comp1").cpl().create("maxAir", "Maximum");
    model.component("comp1").cpl("maxAir").selection().all();

    model.study().create("std1");
    model.study("std1").feature().create("stat", "Stationary");
    System.out.println("V7_FLOW_DEFAULT_SOLVER_BEGIN");
    model.study("std1").run();
    System.out.println("V7_FLOW_DEFAULT_SOLVER_CONVERGED");
    String flux = "spf.rho*(u*nx+v*ny+w*nz)";
    String[] expressions = {
      "abs(intIn(" + flux + "))", "abs(intOut(" + flux + "))",
      "nCells*abs(intIn(" + flux + "))", "nCells*abs(intOut(" + flux + "))",
      "rhoAir*uin*AinCell", "rhoAir*uin*AinStack",
      "abs(abs(intIn(" + flux + "))-rhoAir*uin*AinCell)/(rhoAir*uin*AinCell)",
      "abs(abs(intIn(" + flux + "))-abs(intOut(" + flux + ")))/abs(intIn(" + flux + "))",
      "minAir(p)", "maxAir(p)", "maxAir(p)-minAir(p)"
    };
    String[] names = {
      "channel_inlet_mdot_kg_s", "channel_outlet_mdot_kg_s",
      "stack_inlet_mdot_kg_s", "stack_outlet_mdot_kg_s",
      "expected_channel_mdot_kg_s", "expected_stack_mdot_kg_s",
      "mdot_identity_relative_error", "mass_imbalance_relative",
      "minimum_pressure_Pa", "maximum_pressure_Pa", "pressure_range_Pa"
    };
    model.result().numerical().create("gev", "EvalGlobal");
    model.result().numerical("gev").set("expr", expressions);
    double[][] values = model.result().numerical("gev").getReal();
    for (int i = 0; i < names.length; i++) {
      System.out.printf("%s=%.12g%n", names[i], values[i][0]);
    }
    boolean pass = values[6][0] < 0.001 && values[7][0] < 0.001;
    System.out.println("V7_FLOW_ACCEPTANCE_PASS=" + pass);
    if (!pass) throw new IllegalStateException("V7 flow mass gate failed");
    model.save("comsol/v7_runtime/physical_flow_0p1ms.mph");
  }
}
