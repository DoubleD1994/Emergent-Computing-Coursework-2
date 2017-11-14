import java.util.ArrayList;
import java.util.Random;

public class Solver {

    WindFarmLayoutEvaluator wfle;
    boolean[][] individuals;
    double[] fits;
    Random rand;
    
    int num_individuals;
    int tournamentSize;
    double mutationRate;
    
    ArrayList<double[]> grid;

    public Solver(WindFarmLayoutEvaluator evaluator, int numIndividuals, int tournSize, double mutRate) {
        wfle = evaluator;
        rand = new Random();
        grid = new ArrayList<double[]>();
        
        // set up any parameter here, e.g pop size, cross_rate etc.
        num_individuals = numIndividuals;  // change this to anything you want
        tournamentSize = tournSize;
        mutationRate = mutRate;       
    }
    
    
   
    public void run_cw() {
    

    	
    /************set up grid for scenario chosen  ***************************/	
    // do not change or remove this section of code
    // centers must be > 8*R apart and avoid obstacles
    	
  	double interval = 8.001 * wfle.getTurbineRadius();

  	for (double x=0.0; x<wfle.getFarmWidth(); x+=interval) {
  	    for (double y=0.0; y<wfle.getFarmHeight(); y+=interval) {
                boolean valid = true;
                for (int o=0; o<wfle.getObstacles().length; o++) {
                    double[] obs = wfle.getObstacles()[o];
                    if (x>obs[0] && y>obs[1] && x<obs[2] && y<obs[3]) {
                        valid = false;
                    }
                }

                if (valid) {
                    double[] point = {x, y};
                    grid.add(point);
                }
            }
        }
  	
  	
  	
  	/************initialize a population:*****************************/
  	
        //  the variable grid.size() denotes the 
  		// maximum number of turbines for the given scenario
  	
        individuals = new boolean[num_individuals][grid.size()];
        fits = new double[num_individuals];

        for (int p=0; p<num_individuals; p++) {
            for (int i=0; i<grid.size(); i++) {
                individuals[p][i] = rand.nextBoolean();
                
            }
        }

       /****** evaluate initial population  *************************/
     
        // this populates the fit[] array with the fitness values for each solution
        System.out.println("Initital:");
        double bestFit = evaluate();

        /**** PUT YOUR OPTIMISER CODE HERE ***********/
        
        int parentOne;
    	int parentTwo;
    	
    	boolean[] child;
        
    	double childFitness;
    	
        for (int i=0; i<(9000); i++) {

        	// add some code to evolve a solution
                	
        	//Tournament Selection to get two parents
        	parentOne = tournamentSelect();
        	parentTwo = tournamentSelect();
        	
        	//Uniform Cross-over	
        	child = uniform_crossover(parentOne, parentTwo);
        	
        	// flip mutation
        	flipMutation(child);
        	
        	childFitness = evaluate_individual(child);

        	//Replace worst individual if worse than child
        	worstReplace(child, childFitness);
        	
        	System.out.println("Iteration: " + i + " Child Fitness: " + childFitness );
        	
        	if(childFitness<bestFit)
        	{
        		bestFit = childFitness;
            	System.out.println("Iteration: " + i + " New Best Fitness: " + bestFit);
        	}
        	
        }
        
        System.out.println("Best Fitness Before Hill Climber: " + bestFit);
        
        double bestFitness = fits[0];
        int bestFitId=0;
        
        //get current best solution 
        for(int j=0; j<num_individuals; j++)
    	{
    		double individualFitness = fits[j];
    		
    		if(individualFitness<bestFitness)
    		{
    			bestFitId = j;
    			bestFitness = fits[j];
    			//System.out.println("Best Fitness: " + bestFitness);
    		}
    	}
        
        //run hill climber on current to perfect best solution
        for(int i = 0; i<1000; i++){     	
        	System.out.println("Iteration: " + i);
        	bestFit = swapTurbinePosition(bestFitId);
        }
        
        //Display the best solution with its ID
        System.out.println("Best Fitness: " + bestFit + " With solution: " + bestFitId);
      }
     
    //Tournament selection
    //Selecting best of the individuals to become a parent
    public int tournamentSelect(){
    	
    	int[] potentialParent = new int[tournamentSize];
    	
    	//Selecting random individuals
    	for(int i = 0; i<tournamentSize; i++)
    	{
    		potentialParent[i] = rand.nextInt(num_individuals);
    	}
    	
    	//get fitness of first individual in list
    	int bestParentId = potentialParent[0];
    	double bestFitness = fits[bestParentId];
    	
    	//find best individual
    	for(int i = 1; i<tournamentSize; i++)
    	{
    		int parentId = potentialParent[i];
    		double parentFitness = fits[parentId];   		
    		
    		if(parentFitness < bestFitness)
    		{
    			bestParentId = parentId;
    			bestFitness = parentFitness;
    		}
    	}
    	
    	return bestParentId;
    }
    
    //one-point cross-over
    public boolean[] one_point_crossover(int parentOneId, int parentTwoId)
    {
    	int i;
    	
    	//pick cut point
    	int cutPoint = rand.nextInt(grid.size());
    	
    	//child array
    	boolean[] child = new boolean[grid.size()];
    	

    	
    	// genes from parent one
    	for(i=0; i<cutPoint; i++)
    	{
    		child[i] = individuals[parentOneId][i];
    	}
    	
    	//genes from parent two
    	for(i=cutPoint; i<grid.size(); i++)
    	{
    		child[i] = individuals[parentTwoId][i];
    	}
    	
    	return child;
    }
    
    //uniform cross-over 
    public boolean[] uniform_crossover(int parentOneId, int parentTwoId){
    	
    	//child array
    	boolean[] child = new boolean[grid.size()];
    	
    	double parentChoice;
    	
    	for(int i=0; i<grid.size(); i++)
    	{
    		parentChoice = rand.nextDouble();
    		
    		if(parentChoice<0.5)
    		{
    			child[i] = individuals[parentOneId][i];
    		}
    		else
    		{
    			child[i] = individuals[parentTwoId][i];
    		}
    	}
    	
    	return child;
    }
    
    //flipMutation
    public void flipMutation(boolean[] child){
    	
    	for(int i = 0; i<grid.size(); i++)
    	{
    		if(rand.nextDouble() < mutationRate)
    		{
    			if(child[i])
    			{	
    				
    				child[i] = false;
    			}
    			else
    			{
    				child[i] = true;
    			}
    		}
    	}
    	
    }
    
    //Replacement
    public void worstReplace(boolean[] child, double childFitness){
    	
    	int worstId = 0;
    	double worstFitness = 0;
    	
    	for(int i=0; i<num_individuals; i++)
    	{
    		double individualFitness = fits[i];
    		if(individualFitness>worstFitness)
    		{
    			worstId = i;
    			worstFitness = fits[i];
    		}
    	}
    	
    	if(childFitness<worstFitness)
    	{
    		for(int i=0; i<grid.size();i++)
    		{
    			individuals[worstId][i] = child[i];
    		}
    		fits[worstId] = childFitness;
    	} 	
    }
    
    //Tournament Replace
    //Selecting worst of the individuals to become a parent
    public void tournamentReplace(boolean[] child, double childFitness){
    	
    	int[] potentialWorst = new int[tournamentSize];
    	
    	//Selecting random individuals
    	for(int i = 0; i<tournamentSize; i++)
    	{
    		potentialWorst[i] = rand.nextInt(num_individuals);
    	}
    	
    	//get fitness of first individual in list
    	int worstId = potentialWorst[0];
    	double worstFitness = fits[worstId];
    	
    	//find best individual
    	for(int i = 1; i<tournamentSize; i++)
    	{
    		int parentId = potentialWorst[i];
    		double parentFitness = fits[parentId];   		
    		
    		if(parentFitness > worstFitness)
    		{
    			worstId = parentId;
    			worstFitness = parentFitness;
    		}
    	}
    	
    	if(childFitness<worstFitness)
    	{
    		for(int i=0; i<grid.size();i++)
    		{
    			individuals[worstId][i] = child[i];
    		}
    		fits[worstId] = childFitness;
    	}
    	
    }
    
    //Swap the position of a random turbine to find out if gives a better solution
    public double swapTurbinePosition(int currentBestId){
    	// copy the current solution
    	
    	boolean[] currentBest = new boolean[grid.size()];
    	boolean[] newSoln = new boolean[grid.size()];
    	
        	for (int i=0;i<grid.size();i++){
        		currentBest[i] = individuals[currentBestId][i];
        		newSoln[i]=individuals[currentBestId][i];
        	}

        	// mutate it
    	
        	int place1 = (int)(Math.random()*grid.size());
        	int place2 = (int)(Math.random()*grid.size());
        	
        	while(!currentBest[place1]){
        		place1 = (int)(Math.random()*grid.size());
        	}
        	
        	while(newSoln[place2]){
        		place2 = (int)(Math.random()*grid.size());
        	}
        	
    	    
    	
        	boolean swap = currentBest[place1];
        	newSoln[place1] = currentBest[place2];
        	newSoln[place2] = swap;
        	
        	double newSolnFitness = evaluate_individual(newSoln);
        	double currentBestFitness = fits[currentBestId];
        	
        	if(newSolnFitness<currentBestFitness)
        	{
        		System.out.println("Hill Climber found better solution with fitness: " + newSolnFitness);
        		for(int i=0; i<grid.size();i++)
        		{
        			individuals[currentBestId][i] = newSoln[i];
        		}
        		fits[currentBestId] = newSolnFitness;
        		
        		return newSolnFitness;
        	}    
        	
        	return currentBestFitness;
        }
    
    // evaluate a single chromosome
    private double evaluate_individual(boolean[] child) {
 
       
         int nturbines=0;
         for (int i=0; i<grid.size(); i++) {
                if (child[i]) {
                    nturbines++;
                }
         }
            

          double[][] layout = new double[nturbines][2];
            int l_i = 0;
            for (int i=0; i<grid.size(); i++) {
                if (child[i]) {
                    layout[l_i][0] = grid.get(i)[0];
                    layout[l_i][1] = grid.get(i)[1];
                    l_i++;
                }
            }
	    
	    double coe;
	    if (wfle.checkConstraint(layout)) {
		wfle.evaluate(layout);
		coe = wfle.getEnergyCost();
		//System.out.println("layout valid");
	    } else {
		coe = Double.MAX_VALUE;
	    }

     return coe;
	  
        
    }

    // evaluates the whole population
    private double evaluate() {
        double minfit = Double.MAX_VALUE;
        for (int p=0; p<num_individuals; p++) {
            int nturbines=0;
            for (int i=0; i<grid.size(); i++) {
                if (individuals[p][i]) {
                    nturbines++;
                }
            }

            double[][] layout = new double[nturbines][2];
            int l_i = 0;
            for (int i=0; i<grid.size(); i++) {
                if (individuals[p][i]) {
                    layout[l_i][0] = grid.get(i)[0];
                    layout[l_i][1] = grid.get(i)[1];
                    l_i++;
                }
            }
	    
	    double coe;
	    if (wfle.checkConstraint(layout)) {
		wfle.evaluate(layout);
		coe = wfle.getEnergyCost();
	    } else {
		coe = Double.MAX_VALUE;
	    }

            fits[p] = coe;
            if (fits[p] < minfit) {
                minfit = fits[p];
            }
	   
        }
        System.out.println("min " + minfit);
        
        return minfit;
    }

    
    
}
