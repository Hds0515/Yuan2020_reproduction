import com.comsol.model.*;
import com.comsol.model.util.*;

/**
 * Yuan et al. (2020) equivalent 3D coolant-channel model.
 *
 * Target: COMSOL Multiphysics 6.4, Java API.
 * This is an equivalent homogenized stack/channel geometry because the paper
 * does not publish the original CAD/channel dimensions.
 */
public class Yuan2020Equivalent3D {

  public static Model run() {
    ModelUtil.clear();
    Model model = ModelUtil.create("Model");
    model.label("Yuan2020_equivalent_3D.mph");

    model.param().set("L", "0.30[m]", "Flow-direction length");
    model.param().set("W", "0.082333333[m]", "Width; L*W = 247 cm^2");
    model.param().set("Ain", "3962[mm^2]", "Reported total coolant inlet area");
    model.param().set("Hair", "Ain/W", "Equivalent homogenized air gap");
    model.param().set("Hs", "0.010[m]", "Equivalent graphite solid thickness");
    model.param().set("Tin", "25[degC]");
    model.param().set("uin", "4[m/s]");
    model.param().set("qflux", "2425.5[W/m^2]");
    model.param().set("h_ext", "5[W/(m^2*K)]", "Assumed; not reported");
    model.param().set("Tamb", "25[degC]");
    model.param().set("tol", "1e-6[m]");

    model.component().create("comp1", true);
    model.component("comp1").geom().create("geom1", 3);
    model.component("comp1").geom("geom1").lengthUnit("m");

    model.component("comp1").geom("geom1").feature().create("air", "Block");
    model.component("comp1").geom("geom1").feature("air")
         .set("size", new String[]{"L", "W", "Hair"});
    model.component("comp1").geom("geom1").feature("air")
         .set("pos", new String[]{"0", "0", "0"});
    model.component("comp1").geom("geom1").feature("air").set("selresult", "on");

    model.component("comp1").geom("geom1").feature().create("solid", "Block");
    model.component("comp1").geom("geom1").feature("solid")
         .set("size", new String[]{"L", "W", "Hs"});
    model.component("comp1").geom("geom1").feature("solid")
         .set("pos", new String[]{"0", "0", "Hair"});
    model.component("comp1").geom("geom1").feature("solid").set("selresult", "on");

    model.component("comp1").geom("geom1").run();

    // Coordinate-based selections avoid relying on manually entered entity IDs.
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

    model.component("comp1").selection().create("selHeat", "Box");
    model.component("comp1").selection("selHeat").set("entitydim", 2);
    model.component("comp1").selection("selHeat").set("xmin", "-tol");
    model.component("comp1").selection("selHeat").set("xmax", "L+tol");
    model.component("comp1").selection("selHeat").set("ymin", "-tol");
    model.component("comp1").selection("selHeat").set("ymax", "W+tol");
    model.component("comp1").selection("selHeat").set("zmin", "Hair+Hs-tol");
    model.component("comp1").selection("selHeat").set("zmax", "Hair+Hs+tol");

    String[] zones = {"selZone1", "selZone2", "selZone3"};
    String[] xmin = {"0", "L/3", "2*L/3"};
    String[] xmax = {"L/3", "2*L/3", "L"};
    for (int i = 0; i < 3; i++) {
      model.component("comp1").selection().create(zones[i], "Box");
      model.component("comp1").selection(zones[i]).set("entitydim", 3);
      model.component("comp1").selection(zones[i]).set("xmin", xmin[i]);
      model.component("comp1").selection(zones[i]).set("xmax", xmax[i]);
      model.component("comp1").selection(zones[i]).set("ymin", "0");
      model.component("comp1").selection(zones[i]).set("ymax", "W");
      model.component("comp1").selection(zones[i]).set("zmin", "Hair");
      model.component("comp1").selection(zones[i]).set("zmax", "Hair+Hs");
    }

    model.component("comp1").material().create("matAir", "Common");
    model.component("comp1").material("matAir").selection().named("geom1_air_dom");
    model.component("comp1").material("matAir").propertyGroup("def")
         .set("thermalconductivity", new String[]{"0.0251", "0", "0", "0", "0.0251", "0", "0", "0", "0.0251"});
    model.component("comp1").material("matAir").propertyGroup("def")
         .set("heatcapacity", "1007");
    model.component("comp1").material("matAir").propertyGroup("def")
         .set("density", "1.184");
    model.component("comp1").material("matAir").propertyGroup("def")
         .set("dynamicviscosity", "1.849e-5");

    model.component("comp1").material().create("matGraphite", "Common");
    model.component("comp1").material("matGraphite").selection().named("geom1_solid_dom");
    model.component("comp1").material("matGraphite").propertyGroup("def")
         .set("thermalconductivity", new String[]{"24", "0", "0", "0", "24", "0", "0", "0", "24"});
    model.component("comp1").material("matGraphite").propertyGroup("def")
         .set("heatcapacity", "460");
    model.component("comp1").material("matGraphite").propertyGroup("def")
         .set("density", "2250");

    model.component("comp1").physics().create("spf", "LaminarFlow", "geom1");
    model.component("comp1").physics("spf").selection().named("geom1_air_dom");
    model.component("comp1").physics("spf").feature().create("inl1", "Inlet", 2);
    model.component("comp1").physics("spf").feature("inl1").selection().named("selInlet");
    model.component("comp1").physics("spf").feature("inl1").set("BoundaryCondition", "Velocity");
    model.component("comp1").physics("spf").feature("inl1").set("U0in", "uin");
    model.component("comp1").physics("spf").feature().create("out1", "Outlet", 2);
    model.component("comp1").physics("spf").feature("out1").selection().named("selOutlet");
    model.component("comp1").physics("spf").feature("out1").set("p0", "0[Pa]");

    model.component("comp1").physics().create("ht", "HeatTransfer", "geom1");
    model.component("comp1").physics("ht").feature("solid1").selection().named("geom1_solid_dom");
    model.component("comp1").physics("ht").feature().create("fluid1", "Fluid", 3);
    model.component("comp1").physics("ht").feature("fluid1").selection().named("geom1_air_dom");
    model.component("comp1").physics("ht").feature().create("temp1", "TemperatureBoundary", 2);
    model.component("comp1").physics("ht").feature("temp1").selection().named("selInlet");
    model.component("comp1").physics("ht").feature("temp1").set("T0", "Tin");
    model.component("comp1").physics("ht").feature().create("hf1", "HeatFluxBoundary", 2);
    model.component("comp1").physics("ht").feature("hf1").selection().named("selHeat");
    model.component("comp1").physics("ht").feature("hf1").set("q0", "qflux");

    model.component("comp1").multiphysics().create("nitf1", "NonisothermalFlow", 3);
    model.component("comp1").multiphysics("nitf1").selection().named("geom1_air_dom");

    model.component("comp1").mesh().create("mesh1");
    model.component("comp1").mesh("mesh1").autoMeshSize(3);

    for (int i = 0; i < 3; i++) {
      String tag = "ave" + (i + 1);
      model.component("comp1").cpl().create(tag, "Average");
      model.component("comp1").cpl(tag).selection().named(zones[i]);
    }

    model.study().create("std1");
    model.study("std1").feature().create("stat", "Stationary");
    model.study("std1").feature().create("param", "Parametric");
    model.study("std1").feature("param").set("pname", new String[]{"uin"});
    model.study("std1").feature("param").set("plistarr",
         new String[]{"4[m/s] 5[m/s] 6[m/s] 8[m/s] 10[m/s] 12[m/s]"});
    model.study("std1").feature("param").set("punit", new String[]{"m/s"});

    model.result().table().create("tblZones", "Table");
    model.result().numerical().create("gevZones", "EvalGlobal");
    model.result().numerical("gevZones").set("expr", new String[]{
       "ave1(T)-273.15[K]", "ave2(T)-273.15[K]", "ave3(T)-273.15[K]",
       "maxop1(T)-273.15[K]"
    });
    model.result().numerical("gevZones").set("descr", new String[]{
       "Zone 1 mean", "Zone 2 mean", "Zone 3 mean", "Maximum temperature"
    });
    model.result().numerical("gevZones").set("table", "tblZones");

    model.study("std1").run();
    model.result().numerical("gevZones").setResult();

    model.result().export().create("tblExport", "Table");
    model.result().export("tblExport").set("table", "tblZones");
    model.result().export("tblExport").set("filename", "comsol_zone_temperatures_raw.csv");
    model.result().export("tblExport").run();

    return model;
  }

  public static void main(String[] args) {
    Model model = run();
    model.save("Yuan2020_equivalent_3D.mph");
  }
}
