import com.comsol.model.*;
import com.comsol.model.util.*;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/** Geometry/entity audit only; no physics and no solve. Verified with COMSOL 6.4. */
public class Stage1SelectionAudit {
  private static int[] difference(int[] all, int[] inlet, int[] outlet) {
    Set<Integer> excluded = new HashSet<Integer>();
    for (int id : inlet) excluded.add(id);
    for (int id : outlet) excluded.add(id);
    return Arrays.stream(all).filter(id -> !excluded.contains(id)).toArray();
  }

  private static double selectedArea(Model model, int[] entityIds) {
    GeomMeasureFinal measure = model.component("comp1").geom("geom1").measureFinal();
    measure.selection().geom("geom1", 2);
    measure.selection().set(entityIds);
    return measure.getArea();
  }

  public static void main(String[] args) throws Exception {
    ModelUtil.clear();
    Model model = ModelUtil.create("Model");
    model.param().set("L", "0.30[m]");
    model.param().set("W", "0.082333333[m]");
    model.param().set("Hair", "1[mm]");
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

    model.component("comp1").selection().create("selAllBoundaries", "Explicit");
    model.component("comp1").selection("selAllBoundaries").geom("geom1", 2);
    model.component("comp1").selection("selAllBoundaries").all();
    int[] domains = model.component("comp1").selection("geom1_air_dom").entities();
    int[] allBoundaries = model.component("comp1").selection("selAllBoundaries").entities();
    int[] inlet = model.component("comp1").selection("selInlet").entities();
    int[] outlet = model.component("comp1").selection("selOutlet").entities();
    int[] walls = difference(allBoundaries, inlet, outlet);

    System.out.println("fluid_domains=" + Arrays.toString(domains));
    System.out.println("inlet_boundaries=" + Arrays.toString(inlet));
    System.out.println("outlet_boundaries=" + Arrays.toString(outlet));
    System.out.println("wall_boundaries=" + Arrays.toString(walls));
    if (domains.length != 1 || inlet.length != 1 || outlet.length != 1) {
      throw new IllegalStateException(
          "Expected exactly one fluid domain, one inlet boundary, and one outlet boundary.");
    }

    double inletArea = selectedArea(model, inlet);
    double outletArea = selectedArea(model, outlet);
    double expectedArea = 0.082333333 * 0.001;
    double relativeAreaDifference = Math.abs(inletArea - outletArea) / inletArea;
    System.out.printf("inlet_area_m2=%.12g%n", inletArea);
    System.out.printf("outlet_area_m2=%.12g%n", outletArea);
    System.out.printf("expected_channel_area_m2=%.12g%n", expectedArea);
    System.out.printf("relative_inlet_outlet_area_difference=%.12g%n", relativeAreaDifference);

    if (Math.abs(inletArea - expectedArea) / expectedArea > 1e-6
        || relativeAreaDifference > 1e-12) {
      throw new IllegalStateException("Inlet/outlet area audit failed.");
    }

    model.save("comsol/v6_runtime/stage1_selection_audit.mph");
    System.out.println("STAGE1_SELECTION_AUDIT_PASS");
  }
}
