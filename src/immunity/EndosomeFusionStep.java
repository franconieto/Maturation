package immunity;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import repast.simphony.context.Context;
import repast.simphony.query.space.grid.GridCell;
import repast.simphony.query.space.grid.GridCellNgh;
import repast.simphony.space.continuous.ContinuousSpace;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridPoint;
import repast.simphony.util.ContextUtils;

public class EndosomeFusionStep {
	private static ContinuousSpace<Object> space;
	private static Grid<Object> grid;
	
	public static void fusion (Endosome endosome) {
		String maxRab = Collections.max(endosome.rabContent.entrySet(), Map.Entry.comparingByValue()).getKey();
		// Solo fucionan cisternas.  Esto evita fusión entre vesíclas pequeñas
		if ((endosome.area < Cell.minCistern)) return;
		HashMap<String, Double> rabContent = new HashMap<String, Double>(endosome.getRabContent());
		double areaGolgi = 0d;
		for (String rab : rabContent.keySet()){
			String name = CellProperties.getInstance().rabOrganelle.get(rab);
			if (name.contains("Golgi")) {areaGolgi = areaGolgi + rabContent.get(rab);} 
		}
		boolean isGolgi = false;
		if (areaGolgi/endosome.area >= 0.5) {
			isGolgi = true;
		}
		HashMap<String, Double> membraneContent = new HashMap<String, Double>(endosome.getMembraneContent());
		HashMap<String, Double> solubleContent = new HashMap<String, Double>(endosome.getSolubleContent());
		space = endosome.getSpace();
		grid = endosome.getGrid();
		double cellLimit = 3 * Cell.orgScale;
	
		GridPoint pt = grid.getLocation(endosome);
		// I calculated that the 50 x 50 grid is equivalent to a 750 x 750 nm square
		// Hence, size/15 is in grid units
		int gridSize = (int) Math.round(endosome.size*Cell.orgScale / 15d);// para aumentar fusión.  Lo normal es 15d
		GridCellNgh<Endosome> nghCreator = new GridCellNgh<Endosome>(grid, pt,
				Endosome.class, gridSize, gridSize);
		List<GridCell<Endosome>> cellList = nghCreator.getNeighborhood(true);
		List<Endosome> endosomes_to_delete = new ArrayList<Endosome>();
		double vv = endosome.volume;
		double ss = endosome.area;
		boolean isCistern = (ss * ss * ss / (vv * vv) > 36.01 * Math.PI);// is a sphere

		for (GridCell<Endosome> gr : cellList) {

			// include all endosomes
			for (Endosome end : gr.items()) {
				
				vv = end.volume;
				ss = end.area;
				boolean isCistern2 = (ss * ss * ss / (vv * vv) < 36.01 * Math.PI);
				if (isGolgi){
					if (end != endosome  // it is not itself
							&& (end.volume <= 4/3 * Math.PI*Math.pow(Cell.rcyl, 3))// use to be endosome.volume) // the other is smaller
							&& (Math.random() < EndosomeAssessCompatibility.compatibles(endosome, end))) {

						endosomes_to_delete.add(end);
					}
				}
				else { // it is not a Golgi
					if (end != endosome  // it is not itself
						&& (end.volume < endosome.volume)
						&& (Math.random() < EndosomeAssessCompatibility.compatibles(endosome, end))) {

						endosomes_to_delete.add(end);
					}
				}
			}
		}
		for (Endosome endosome2 : endosomes_to_delete) {
			endosome.volume = endosome.volume + endosome2.volume;
			endosome.area = endosome.area + endosome2.area;
			endosome.rabContent = sumRabContent(endosome, endosome2);
			endosome.membraneContent = sumMembraneContent(endosome, endosome2);
			endosome.solubleContent = sumSolubleContent(endosome, endosome2);
			Context<Object> context = ContextUtils.getContext(endosome);
			context.remove(endosome2);
		}
		double rsphere = Math.pow(endosome.volume * 3d / 4d / Math.PI, (1d / 3d));
		double size = rsphere;
		endosome.speed = 1d/ size;
		Endosome.endosomeShape(endosome);
		endosome.getEndosomeTimeSeries().clear();
		endosome.getRabTimeSeries().clear();
//		The time series will be re-calculated by COPASI call in the next tick
//		

	}
	private static HashMap<String, Double> sumRabContent(Endosome endosome1,
			Endosome endosome2) {

		HashMap<String, Double> rabSum = new HashMap<String, Double>();
		for (String key1 : endosome1.rabContent.keySet()) {
			if (endosome2.rabContent.containsKey(key1)) {
				double sum = endosome1.rabContent.get(key1)
						+ endosome2.rabContent.get(key1);
				rabSum.put(key1, sum);
			} else
				rabSum.put(key1, endosome1.rabContent.get(key1));
		}
		for (String key2 : endosome2.rabContent.keySet()) {
			if (!endosome1.rabContent.containsKey(key2)) {
				rabSum.put(key2, endosome2.rabContent.get(key2));
			}
		}

		return rabSum;
	}

	private static HashMap<String, Double> sumMembraneContent(Endosome endosome1,
			Endosome endosome2) {
		HashMap<String, Double> memSum = new HashMap<String, Double>();
		for (String key1 : endosome1.membraneContent.keySet()) {
			if (endosome2.membraneContent.containsKey(key1)) {
				double sum = endosome1.membraneContent.get(key1)
						+ endosome2.membraneContent.get(key1);
				memSum.put(key1, sum);
			} else
				memSum.put(key1, endosome1.membraneContent.get(key1));
		}
		for (String key2 : endosome2.membraneContent.keySet()) {
			if (!endosome1.membraneContent.containsKey(key2)) {
				double sum = endosome2.membraneContent.get(key2);
				memSum.put(key2, sum);
			}
		}
		return memSum;
	}

	private static HashMap<String, Double> sumSolubleContent(Endosome endosome1,
			Endosome endosome2) {
		HashMap<String, Double> solSum = new HashMap<String, Double>();
		for (String key1 : endosome1.solubleContent.keySet()) {
			if (endosome2.solubleContent.containsKey(key1)) {
				double sum = endosome1.solubleContent.get(key1)
						+ endosome2.solubleContent.get(key1);
				solSum.put(key1, sum);
			} else
				solSum.put(key1, endosome1.solubleContent.get(key1));
		}
		for (String key2 : endosome2.solubleContent.keySet()) {
			if (!endosome1.solubleContent.containsKey(key2)) {
				double sum = endosome2.solubleContent.get(key2);
				solSum.put(key2, sum);
			}
		}
		return solSum;
	}


}
