package mitmsearch.mitm;

import gurobi.*;

import java.util.List;
import java.util.ArrayList;
import java.util.stream.IntStream;
import java.util.stream.Collectors;
import java.io.IOException;
import java.io.File;
import java.io.FileWriter;
public class Mitm {
  private final int Rounds;
  private final GRBModel model;
  private FileWriter logfile ;
  private final MitmFactory factory;
  private final GRBVar[][][][] Pi_init;
  private final GRBVar[][][] Cond;
  private final GRBVar[][][][][] DA;
  private final GRBVar[][][][] DP;
  private final GRBVar[][][] DC1;
  private final GRBVar[][][][] DP2;
  private final GRBVar[][][] DC12;
  private final GRBVar[][][][][] DB; 
  private final GRBVar[][][][] DC2;   
  private final GRBLinExpr   objective;
  private final GRBVar[] obj;
  private final GRBVar[][] dom;
  private static final int[][] rho = new int[][]{{0,36,3,41,18},{1,44,10,45,2},{62,6,43,15,61},{28,55,25,21,56},{27,20,39,8,14}};

  /**
   * @param env the Gurobi environment
   */
  public Mitm(final GRBEnv env, final int Rounds) throws GRBException {
    model = new GRBModel(env);
    this.Rounds = Rounds;

    factory = new MitmFactory(model);
    Pi_init = new GRBVar[5][5][32][3];
    Cond = new GRBVar[5][32][4];
    DA = new GRBVar[Rounds+1][5][5][32][3];
    DP = new GRBVar[Rounds][5][32][3];
    DP2 = new GRBVar[Rounds][5][32][3];
    DC1 = new GRBVar[Rounds][5][32];
    DC12 = new GRBVar[Rounds][5][32];
    DB = new GRBVar[Rounds][5][5][32][3];
    DC2 = new GRBVar[Rounds][5][5][32];
    // Initialization
    for (int i = 0; i < 5; i++) 
      for (int j = 0; j < 5; j++) 
        for (int k = 0; k < 32; k++) 
          for (int l = 0; l < 3; l++) {     
            Pi_init[i][j][k][l] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "Pi_init_"+i+"_"+j+"_"+k+"_"+l);            
    }

    for (int j = 0; j < 5; j++) 
      for (int k = 0; k < 32; k++) 
        for (int l = 0; l < 4; l++) {     
          Cond[j][k][l] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "Cond_"+j+"_"+k+"_"+l);            
    }

    for (int round = 0; round < Rounds+1; round++)
      for (int i = 0; i < 5; i++) 
        for (int j = 0; j < 5; j++) 
	  for (int k = 0; k < 32; k++) 
            for (int l = 0; l < 3; l++) {     
              DA[round][i][j][k][l] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "DA_"+round+"_"+i+"_"+j+"_"+k+"_"+l);            
    }

    for (int round = 0; round < Rounds; round++)
      for (int i = 0; i < 5; i++) 
        for (int j = 0; j < 5; j++) 
	  for (int k = 0; k < 32; k++)           
            for (int l = 0; l < 3; l++) {     
              DB[round][i][j][k][l] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "DB_"+round+"_"+i+"_"+j+"_"+k+"_"+l);    	          
    }

    for (int round = 0; round < Rounds; round++)
      for (int i = 0; i < 5; i++) 
        for (int j = 0; j < 5; j++) 
	  for (int k = 0; k < 32; k++) {    
            DC2[round][i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "DC2_"+round+"_"+i+"_"+j+"_"+k); 
    	             
    }
    

    for (int round = 0; round < Rounds; round++)
      for (int i = 0; i < 5; i++) 
	for (int k = 0; k < 32; k++) 
	  for (int l = 0; l < 3; l++)  {         
            DP[round][i][k][l] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "DP_"+round+"_"+i+"_"+k+"_"+l);               
    }

    for (int round = 0; round < Rounds; round++)
      for (int i = 0; i < 5; i++) 
	for (int k = 0; k < 32; k++) 
	  for (int l = 0; l < 3; l++)  {         
	    DP2[round][i][k][l] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "DP2_"+round+"_"+i+"_"+k+"_"+l);              
    }
    for (int round = 0; round < Rounds; round++)
      for (int i = 0; i < 5; i++) 
	for (int k = 0; k < 32; k++) {       
          DC1[round][i][k] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "DC1_"+round+"_"+i+"_"+k);                         
    }

    for (int round = 0; round < Rounds; round++)
      for (int i = 0; i < 5; i++) 
	for (int k = 0; k < 32; k++) {      
	  DC12[round][i][k] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "DC12_"+round+"_"+i+"_"+k);               
    }

    dom = new GRBVar[2][32];
    for (int i = 0; i < 2; i++) 
      for (int k = 0; k < 32; k++){
        dom[i][k]   = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "dom_"+i+"_"+k);
    }

    
    //fixed input
    factory.addfixed_in(Pi_init,DA,Cond);
    ArticleTrail();
    factory.addconstr_in(Cond);
    

    // Constraints
    factory.addfivexor_red(DA, DP, DC1);

    factory.addtwoxor_red(DP2, DP, DC12);
    factory.addTheta_red(DA, DP2, DB, DC2);

    factory.addSbox_nc(DB, DA);

    factory.addDoM(DA, dom);
  
    GRBVar[][] beta = new GRBVar[4][32];
    for (int i = 0; i < 4; i++) 
      for (int k = 0; k < 32; k++){
        beta[i][k]   = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "beta_"+i+"_"+k);
    };
    factory.betaConstraints(Pi_init, beta);



    objective = new GRBLinExpr();  
    GRBLinExpr DoF_red = new GRBLinExpr();
    GRBLinExpr DoF_blue = new GRBLinExpr();
    GRBLinExpr DoM = new GRBLinExpr();
    
    GRBVar Obj = model.addVar(0.0, 1600.0, 0.0, GRB.INTEGER, "Obj");
    
    obj = new GRBVar[3];
    obj[0] = model.addVar(0.0, 1600.0, 0.0, GRB.INTEGER, "Obj1"); 
    obj[1] = model.addVar(0.0, 1600.0, 0.0, GRB.INTEGER, "Obj2"); 
    obj[2] = model.addVar(0.0, 1600.0, 0.0, GRB.INTEGER, "Obj3"); 
                
      
    for (int k = 0; k < 32; k++) {    
      DoF_blue.addTerm(1.0, Pi_init[0][0][k][0]);
      DoF_red.addTerm(1.0, Pi_init[0][0][k][2]);
      DoF_blue.addTerm(-1.0, beta[0][k]);
      DoF_red.addTerm(-1.0, beta[0][k]);
      
      DoF_blue.addTerm(1.0, Pi_init[0][1][k][0]);
      DoF_red.addTerm(1.0, Pi_init[0][1][k][2]);
      DoF_blue.addTerm(-1.0, beta[1][k]);
      DoF_red.addTerm(-1.0, beta[1][k]);
       
      DoF_blue.addTerm(1.0, Pi_init[0][2][k][0]);
      DoF_red.addTerm(1.0, Pi_init[0][2][k][2]);
      DoF_blue.addTerm(-1.0, beta[2][k]);
      DoF_red.addTerm(-1.0, beta[2][k]);

      DoF_blue.addTerm(1.0, Pi_init[0][4][k][0]);
      DoF_red.addTerm(1.0, Pi_init[0][4][k][2]);
      DoF_blue.addTerm(-1.0, beta[3][k]);
      DoF_red.addTerm(-1.0, beta[3][k]);
    }

    
    for (int round = 0; round < Rounds; round ++) 
      for (int i = 0; i < 5; i++) 
        for (int k = 0; k < 32; k++) {
          DoF_red.addTerm(-1.0, DC1[round][i][k]);
          DoF_red.addTerm(-1.0, DC12[round][i][k]);
          for (int j = 0; j < 5; j++) {
            DoF_red.addTerm(-1.0, DC2[round][i][j][k]);        
          }
    }
    

   

    for (int i = 0; i < 2; i++) 
      for (int k = 0; k < 32; k++) {
        DoM.addTerm(1.0, dom[i][k]);
    }
    

    objective.addTerm(1.0, Obj);


    model.addConstr(DoF_red, GRB.EQUAL, obj[0], "");
    model.addConstr(DoF_blue, GRB.EQUAL, obj[1], "");
    model.addConstr(DoM, GRB.EQUAL, obj[2], ""); 
    
    model.addConstr(objective, GRB.LESS_EQUAL, DoF_blue, "");
    model.addConstr(objective, GRB.LESS_EQUAL, DoF_red, "");
    model.addConstr(objective, GRB.LESS_EQUAL, DoM, "");
    model.setObjective(objective, GRB.MAXIMIZE);
  }

  public List<MitmSolution> solve(final int nSolutions, final boolean nonOptimalSolutions, final int minObjValue, final int Threads) throws GRBException {
    model.read("tune1.prm");
    model.write("model.lp");
    model.set(GRB.IntParam.Threads, Threads);
    if (minObjValue != -1)
      model.addConstr(objective, GRB.EQUAL, minObjValue, "objectiveFix");
    model.set(GRB.IntParam.DualReductions, 0);
    

    model.optimize();
    model.write("output.sol");
    //model.computeIIS();
    //model.write("model1.ilp");
    return getAllFoundSolutions();
  }

  public void dispose() throws GRBException {
    model.dispose();
  }

  public List<MitmSolution> getAllFoundSolutions() throws GRBException {
    return IntStream.range(0, model.get(GRB.IntAttr.SolCount)).boxed()
      .map(solNb -> getSolution(solNb))
      .collect(Collectors.toList());
  }

  private MitmSolution getSolution(final int solutionNumber) {
    try {
      model.set(GRB.IntParam.SolutionNumber, solutionNumber);
      int[][][][] Pi_initValue     = new int[5][5][32][3];
      int[][][] CondValue     = new int[5][32][4];
      int[][][][][] DAValue     = new int[Rounds+1][5][5][32][3];
      int[][][][] DPValue  = new int[Rounds][5][32][3];
      int[][][][] DP2Value  = new int[Rounds][5][32][3];
      int[][][] DC1Value = new int[Rounds][5][32];
      int[][][] DC12Value = new int[Rounds][5][32];
      int[][][][][] DBValue  = new int[Rounds][5][5][32][3];
      int[][][][] DC2Value  = new int[Rounds][5][5][32];
      int[][] domValue  = new int[2][32];
      int[] objValue  = new int[3];

      for (int round = 0; round < Rounds; round++)
        for (int i = 0; i < 5; i++) 
          for (int j = 0; j < 5; j++) 
	    for (int k = 0; k < 32; k++) {
              DC2Value[round][i][j][k]  = (int) Math.round(DC2[round][i][j][k].get(GRB.DoubleAttr.Xn));
              
	      for (int l = 0; l < 3; l++)  { 
                DAValue[round][i][j][k][l]  = (int) Math.round(DA[round][i][j][k][l].get(GRB.DoubleAttr.Xn));
                DBValue[round][i][j][k][l]  = (int) Math.round(DB[round][i][j][k][l].get(GRB.DoubleAttr.Xn));
              }
      }
      for (int i = 0; i < 5; i++) 
        for (int j = 0; j < 5; j++) 
	  for (int k = 0; k < 32; k++) 
            for (int l = 0; l < 3; l++)  { 
              DAValue[Rounds][i][j][k][l]  = (int) Math.round(DA[Rounds][i][j][k][l].get(GRB.DoubleAttr.Xn));
              Pi_initValue[i][j][k][l]  = (int) Math.round(Pi_init[i][j][k][l].get(GRB.DoubleAttr.Xn));
            }

      for (int j = 0; j < 5; j++) 
	for (int k = 0; k < 32; k++) 
          for (int l = 0; l < 4; l++)  { 
            CondValue[j][k][l]  = (int) Math.round(Cond[j][k][l].get(GRB.DoubleAttr.Xn));
          }
      for (int round = 0; round < Rounds; round++)
        for (int i = 0; i < 5; i++) 
	  for (int k = 0; k < 32; k++) {
            DC1Value[round][i][k]  = (int) Math.round(DC1[round][i][k].get(GRB.DoubleAttr.Xn));
            DC12Value[round][i][k]  = (int) Math.round(DC12[round][i][k].get(GRB.DoubleAttr.Xn));
	    
	    for (int l = 0; l < 3; l++)  { 
              DPValue[round][i][k][l]  = (int) Math.round(DP[round][i][k][l].get(GRB.DoubleAttr.Xn));
	      DP2Value[round][i][k][l]  = (int) Math.round(DP2[round][i][k][l].get(GRB.DoubleAttr.Xn));
	    }
      }
      
      for (int i = 0; i < 2; i++) 
	for (int k = 0; k < 32; k++) 
          domValue[i][k]  = (int) Math.round(dom[i][k].get(GRB.DoubleAttr.Xn));

      for (int i = 0; i < 3; i++) 
        objValue[i]  = (int) Math.round(obj[i].get(GRB.DoubleAttr.Xn));
	   
      return new MitmSolution(Rounds, (int) Math.round(model.get(GRB.DoubleAttr.PoolObjVal)), Pi_initValue, CondValue, DAValue,DBValue, DC2Value, DPValue, DP2Value, DC1Value, DC12Value, domValue, objValue);
    } catch (GRBException e) {
      System.out.println("Error code: " + e.getErrorCode() + ". " + e.getMessage());
      e.printStackTrace();
      System.exit(1);
      return null; // Can't access
    }
  }
  private void ArticleTrail() throws GRBException {
    int[][][][] Pi_article = {{{{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{1,1,0},{0,1,1},{0,1,1},{0,1,1},{1,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{1,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{1,1,1},{0,1,1},{0,1,1},{0,1,1},{1,1,1},{0,1,1},{0,1,1},{0,1,1}},{{0,1,1},{1,1,0},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{1,1,0},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{1,1,1},{0,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1}},{{0,1,1},{1,1,0},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{1,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{1,1,0},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{1,1,1},{0,1,1},{0,1,1},{0,1,1},{1,1,1},{0,1,1},{0,1,1}},{{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1}},{{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{1,1,1},{0,1,1},{0,1,1},{1,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{1,1,0},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{1,1,1},{0,1,1},{0,1,1},{1,1,1},{0,1,1},{0,1,1},{0,1,1}}},{{{0,1,1},{0,1,1},{0,1,1},{0,1,1},{1,1,1},{0,1,1},{0,1,1},{0,1,1},{1,1,1},{0,1,1},{0,1,1},{0,1,1},{1,1,0},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{1,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{1,1,0},{0,1,1},{0,1,1},{0,1,1}},{{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1}},{{0,1,1},{1,1,1},{0,1,1},{0,1,1},{1,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{1,1,1},{0,1,1},{0,1,1},{1,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{1,1,0},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1}},{{1,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{1,1,0},{0,1,1},{0,1,1},{0,1,1},{1,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{1,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{1,1,1},{0,1,1},{0,1,1},{0,1,1}},{{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{1,1,0},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{1,1,1},{0,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{0,1,1},{1,1,0},{0,1,1},{0,1,1},{0,1,1}}},{{{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1}},{{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1}},{{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1}},{{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1}},{{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1}}},{{{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1}},{{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1}},{{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1}},{{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1}},{{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1}}},{{{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1}},{{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1}},{{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1}},{{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1}},{{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1}}}};
    for (int i = 0; i < 5; i++) 
      for (int j = 0; j < 5; j++) 
        for (int k = 0; k < 32; k++) 
          for (int l = 0; l < 3; l++) {     
             model.addConstr(Pi_init[i][j][k][l], GRB.EQUAL, Pi_article[i][j][k][l], "");            
    }
    
  }

}
