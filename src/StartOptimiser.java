public class StartOptimiser {
  public static int nSc;

  public static void main(String argv[]) throws Exception {
	  
	int numIndividuals = 100;
	int trnSize = 7;
	double mutRate = 0.07;
	  
	//nSc = Integer.parseInt(argv[0]);
	KusiakLayoutEvaluator eval = new KusiakLayoutEvaluator();
	WindScenario sc = new WindScenario("./Scenarios/practice_1.xml");	  
	eval.initialize(sc);
	Solver algorithm = new Solver(eval, numIndividuals, trnSize, mutRate);
	algorithm.run_cw();
  }
}
