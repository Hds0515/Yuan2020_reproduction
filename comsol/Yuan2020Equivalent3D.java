import com.comsol.model.*;
import com.comsol.model.util.*;

/**
 * Stage 1: isothermal 0.1 m/s flow in one representative parallel channel.
 *
 * Geometry decision B is used. The solved channel has area W*Hair and every
 * extensive total is multiplied by areaScale=Ain/(W*Hair). Therefore the
 * reported total inlet mass flow is exactly rho*uin*Ain while the hydraulic
 * diameter remains that of the declared representative channel.
 */
public class Yuan2020Equivalent3D {

  public static Model run() {
    ModelUtil.clear();
    Model model = ModelUtil.create("Model");
    model.label("flow_0p1ms_converged.mph");

    model.param().set("L", "0.30[m]", "Representative flow length");
    model.param().set("W", "0.082333333[m]", "Representative channel width");
    model.param().set("Hair", "1[mm]", "Declared representative channel height");
    model.param().set("Ain", "3962[mm^2]", "Reported total coolant inlet area");
    model.param().set("areaScale", "Ain/(W*Hair)", "Parallel-channel multiplier");
    model.param().set("uin", "0.1[m/s]", "Stage 1 inlet velocity");
    model.param().set("tol", "1e-7[m]");

    model.component().create("comp1", true);
    model.component("comp1").geom().create("geom1", 3);
    model.component("comp1").geom("geom1").lengthUnit("m");
    model.component("comp1").geom("geom1").feature().create("air", "Block");
    model.component("comp1").geom("geom1").feature("air")
         .set("size", new String[]{"L", "W", "Hair"});
    model.component("comp1").geom("geom1").feature("air").set("selresult", "on");
    model.component("comp1").geom("geom1").run();
    System.out.println("COMSOL_STAGE1 geometry_complete");

    model.component("comp1").selection().create("selInlet", "Box");
    model.component("comp1").selection("selInlet").set("entitydim", 2);
    model.component("comp1").selection("selInlet").set("xmin", "-tol");
    model.component("comp1").selection("selInlet").set("xmax", "tol");
    model.component("comp1").selection("selInlet").set("ymin", "-tol");
    model.component("comp1").selection("selInlet").set("ymax", "W+tol");
    model.component("comp1").selection("selInlet").set("zmin", "-tol");
    model.component("comp1").selection("selInlet").set("zmax", "Hair+tol");

    model.component("comp1").selection().create("selOutlet", "Box");
    model.component("comp1").selection("selOutlet").set("entitydim", 2);
    model.component("comp1").selection("selOutlet").set("xmin", "L-tol");
    model.component("comp1").selection("selOutlet").set("xmax", "L+tol");
    model.component("comp1").selection("selOutlet").set("ymin", "-tol");
    model.component("comp1").selection("selOutlet").set("ymax", "W+tol");
    model.component("comp1").selection("selOutlet").set("zmin", "-tol");
    model.component("comp1").selection("selOutlet").set("zmax", "Hair+tol");

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
    System.out.println("COMSOL_STAGE1 selections_complete");

    model.component("comp1").material().create("matAir", "Common");
    model.component("comp1").material("matAir").selection().named("geom1_air_dom");
    model.component("comp1").material("matAir").propertyGroup("def")
         .set("density", "1.184[kg/m^3]");
    model.component("comp1").material("matAir").propertyGroup("def")
         .set("dynamicviscosity", "1.849e-5[Pa*s]");

    model.component("comp1").physics().create("spf", "LaminarFlow", "geom1");
    model.component("comp1").physics("spf").selection().named("geom1_air_dom");
    model.component("comp1").physics("spf").prop("PhysicalModelProperty")
         .set("StokesFlowProp", "1");
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
    System.out.println("COMSOL_STAGE1 physics_complete");

    model.component("comp1").mesh().create("mesh1");
    model.component("comp1").mesh("mesh1").feature().create("edg1", "Edge");
    model.component("comp1").mesh("mesh1").feature("edg1").selection()
         .named("selThicknessEdges");
    model.component("comp1").mesh("mesh1").feature("edg1").create("dis1", "Distribution");
    model.component("comp1").mesh("mesh1").feature("edg1").feature("dis1")
         .set("numelem", 4);
    model.component("comp1").mesh("mesh1").feature().create("edg2", "Edge");
    model.component("comp1").mesh("mesh1").feature("edg2").selection()
         .named("selWidthEdges");
    model.component("comp1").mesh("mesh1").feature("edg2").create("dis1", "Distribution");
    model.component("comp1").mesh("mesh1").feature("edg2").feature("dis1")
         .set("numelem", 12);
    model.component("comp1").mesh("mesh1").feature().create("map1", "Map");
    model.component("comp1").mesh("mesh1").feature("map1").selection().named("selInlet");
    model.component("comp1").mesh("mesh1").feature("map1").create("size1", "Size");
    model.component("comp1").mesh("mesh1").feature("map1").feature("size1")
         .set("custom", "on");
    model.component("comp1").mesh("mesh1").feature("map1").feature("size1")
         .set("hmax", "W/12");
    model.component("comp1").mesh("mesh1").feature().create("swe1", "Sweep");
    model.component("comp1").mesh("mesh1").feature("swe1").selection()
         .named("geom1_air_dom");
    model.component("comp1").mesh("mesh1").feature("swe1").selection("sourceface")
         .named("selInlet");
    model.component("comp1").mesh("mesh1").feature("swe1").selection("targetface")
         .named("selOutlet");
    model.component("comp1").mesh("mesh1").feature("swe1").create("dis1", "Distribution");
    model.component("comp1").mesh("mesh1").feature("swe1").feature("dis1")
         .set("numelem", 40);
    model.component("comp1").mesh("mesh1").run();
    System.out.println("COMSOL_STAGE1 swept_mesh_complete");

    model.component("comp1").cpl().create("intIn", "Integration");
    model.component("comp1").cpl("intIn").selection().named("selInlet");
    model.component("comp1").cpl().create("intOut", "Integration");
    model.component("comp1").cpl("intOut").selection().named("selOutlet");
    model.component("comp1").cpl().create("maxAir", "Maximum");
    model.component("comp1").cpl("maxAir").selection().named("geom1_air_dom");
    model.component("comp1").cpl().create("minOut", "Minimum");
    model.component("comp1").cpl("minOut").selection().named("selOutlet");

    model.study().create("std1");
    model.study("std1").feature().create("stat", "Stationary");
    model.sol().create("sol1");
    model.sol("sol1").study("std1");
    model.sol("sol1").attach("std1");
    model.sol("sol1").create("st1", "StudyStep");
    model.sol("sol1").feature("st1").set("study", "std1");
    model.sol("sol1").feature("st1").set("studystep", "stat");
    model.sol("sol1").create("v1", "Variables");
    model.sol("sol1").feature("v1").set("control", "stat");
    model.sol("sol1").create("s1", "Stationary");
    model.sol("sol1").feature("s1").set("control", "stat");
    model.sol("sol1").feature("s1").create("fc1", "FullyCoupled");
    model.sol("sol1").feature("s1").feature("fc1").set("maxiter", 50);
    model.sol("sol1").feature("s1").create("d1", "Direct");
    model.sol("sol1").feature("s1").feature("d1").set("linsolver", "pardiso");
    model.sol("sol1").feature("s1").feature("fc1").set("linsolver", "d1");
    model.sol("sol1").feature("s1").feature().remove("fcDef");
    System.out.println("COMSOL_STAGE1 stokes_initialization_fully_coupled_pardiso");
    model.sol("sol1").runAll();
    model.component("comp1").physics("spf").prop("PhysicalModelProperty")
         .set("StokesFlowProp", "0");
    model.sol("sol1").feature("v1").set("initmethod", "sol");
    model.sol("sol1").feature("v1").set("initsol", "sol1");
    System.out.println("COMSOL_STAGE1 inertial_laminar_continuation");
    model.sol("sol1").runAll();
    System.out.println("COMSOL_STAGE1 converged");

    String flux = "spf.rho*(u*nx+v*ny+w*nz)";
    model.result().table().create("tblStage1", "Table");
    model.result().numerical().create("gevStage1", "EvalGlobal");
    model.result().numerical("gevStage1").set("expr", new String[]{
       "areaScale*abs(intIn(" + flux + "))",
       "areaScale*abs(intOut(" + flux + "))",
       "abs(abs(intIn(" + flux + "))-abs(intOut(" + flux + ")))/abs(intIn(" + flux + "))",
       "1.184[kg/m^3]*uin*Ain",
       "maxAir(sqrt(u^2+v^2+w^2))",
       "maxAir(abs(p))",
       "minOut(u*nx+v*ny+w*nz)",
       "areaScale"
    });
    model.result().numerical("gevStage1").set("descr", new String[]{
       "scaled inlet mass flow", "scaled outlet mass flow", "relative mass imbalance",
       "rho u Ain target", "maximum velocity", "maximum absolute pressure",
       "minimum outlet-normal velocity", "area scale"
    });
    model.result().numerical("gevStage1").set("table", "tblStage1");
    model.result().numerical("gevStage1").setResult();
    model.result().export().create("tblExport", "Table");
    model.result().export("tblExport").set("table", "tblStage1");
    model.result().export("tblExport").set("filename", "results_v4/stage1_flow_metrics.csv");
    model.result().export("tblExport").run();
    return model;
  }

  public static void main(String[] args) throws Exception {
    try {
      Model model = run();
      model.save("results_v4/flow_0p1ms_converged.mph");
    } catch (Exception exception) {
      exception.printStackTrace(System.err);
      throw exception;
    }
  }
}
