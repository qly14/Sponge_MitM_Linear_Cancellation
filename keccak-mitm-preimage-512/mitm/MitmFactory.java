package mitmsearch.mitm;

import gurobi.*;

public class MitmFactory {
  private GRBModel model;
  private static final int[][] rho = new int[][]{{0,36,3,41,18},{1,44,10,45,2},{62,6,43,15,61},{28,55,25,21,56},{27,20,39,8,14}};


  public MitmFactory(final GRBModel model) {
    this.model = model;
  }

  public void addfixed_in(GRBVar[][][][] Pi_init, GRBVar[][][][][] DA, GRBVar[][][] Cond) throws GRBException {
    
    for (int i = 0; i < 5; i++) 
      for (int j = 0; j < 5; j++) {
        for (int k = 0; k < 32; k++) {
	  for (int l = 0; l < 3; l++)  {
            model.addConstr(Pi_init[0][0][k][l], GRB.EQUAL, Pi_init[1][3][(k+rho[0][1])%32][l], "");  
 	    model.addConstr(Pi_init[0][2][(k+rho[1][0])%32][l], GRB.EQUAL, Pi_init[1][0][(k+rho[1][1])%32][l], "");
            model.addConstr(Pi_init[0][4][(k+rho[2][0])%32][l], GRB.EQUAL, Pi_init[1][2][(k+rho[2][1])%32][l], "");
            model.addConstr(Pi_init[0][1][(k+rho[3][0])%32][l], GRB.EQUAL, Pi_init[1][4][(k+rho[3][1])%32][l], "");
            //padding
            if (k >= 32-4) {
              model.addConstr(Pi_init[1][4][(k+rho[3][1])%32][l], GRB.EQUAL, 1, "");
            } 
          }
        }
        if ((i==0&j==0) | (i==0&j==1) | (i==0&j==2) | (i==0&j==4) | (i==1&j==0) | (i==1&j==2) | (i==1&j==3) | (i==1&j==4)) {   
	  for (int k = 0; k < 32; k++) {         
            model.addConstr(Pi_init[i][j][k][1], GRB.EQUAL, 1, ""); 
            GRBLinExpr known = new GRBLinExpr();
            known.addTerm(1, Pi_init[i][j][k][0]);
            known.addTerm(1, Pi_init[i][j][k][2]);
            model.addConstr(known, GRB.GREATER_EQUAL, 1, "");      
	  }
        }
        else {
          for (int k = 0; k < 32; k++) {      
            model.addConstr(Pi_init[i][j][k][0], GRB.EQUAL, 1, "");   
            model.addConstr(Pi_init[i][j][k][1], GRB.EQUAL, 1, "");  
            model.addConstr(Pi_init[i][j][k][2], GRB.EQUAL, 1, ""); 
	  }
        }
      }
    //DA[0]
    for (int j = 0; j < 5; j ++)
      for (int k = 0; k < 32; k ++) {
      //Avoid Red mul Blue
      model.addConstr(linExprOf(Pi_init[0][j][k][0], Pi_init[1][j][k][2]), GRB.GREATER_EQUAL, 1, "");
      model.addConstr(linExprOf(Pi_init[0][j][k][2], Pi_init[1][j][k][0]), GRB.GREATER_EQUAL, 1, "");
      //Out DA0 DA4
      double[] t1 = {0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1};
      model.addConstr(linExprOf(t1, Pi_init[0][j][k][0], Pi_init[0][j][k][2], Pi_init[1][j][k][0], Pi_init[1][j][k][2], DA[0][0][j][k][0], DA[0][0][j][k][2], DA[0][4][j][k][0], DA[0][4][j][k][2], Cond[j][k][0], Cond[j][k][1], Cond[j][k][2]), GRB.EQUAL, 0, ""); 
      double[][] t2 = {{0, 0, 0, 0, 0, 0, 1, 1, -1, -2, -1},
			{1, 1, 0, 0, 0, 0, -1, -1, 0, 1, 0},
			{0, 0, 1, 1, 0, 0, -1, -1, 1, 0, 0},
			{-2, 0, -1, -1, 0, 0, 2, 1, -2, -1, 0},
			{0, -1, 0, -1, 0, 0, 0, 1, 0, 0, 0},
			{0, 1, 0, 0, 1, 0, -1, -1, 1, 1, -1},
			{1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 1, 0, 0, 0, -1, 1, 0, 0},
			{-1, 0, 0, -1, -1, 0, 2, 1, -2, -1, 1},
			{0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 0},
			{-1, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 1, 0, -1, -1, 0},
			{0, -1, 0, -1, -1, 0, 1, 2, -1, -2, 0},
			{1, 0, 0, 0, 0, 0, -1, 0, 0, 1, 0},
			{1, 0, 1, 0, -1, 0, -1, 0, 0, 0, 1}};
      double[] c2 = {0, 0, 0, -2, -1, 0, 0, 0, -1, 0, -1, 0, -1, 0, 0};
      for (int q = 0; q < 15; q ++)   
        model.addConstr(linExprOf(t2[q], Pi_init[0][j][k][0], Pi_init[0][j][k][2], Pi_init[1][j][k][0], Pi_init[1][j][k][2], DA[0][0][j][k][0], DA[0][0][j][k][2], DA[0][4][j][k][0], DA[0][4][j][k][2], Cond[j][k][0], Cond[j][k][1], Cond[j][k][2]), GRB.GREATER_EQUAL, c2[q], "");   

      model.addConstr(DA[0][0][j][k][1], GRB.EQUAL, 1, "");
      model.addConstr(DA[0][4][j][k][1], GRB.EQUAL, 1, "");

      //DA3
      model.addConstr(DA[0][3][j][k][1], GRB.EQUAL, 1, "");
      
      double[] t3 = {1, 1, -1, -1, 1};
      model.addConstr(linExprOf(t3, Pi_init[0][j][k][0], Pi_init[0][j][k][2], DA[0][3][j][k][0], DA[0][3][j][k][2], Cond[j][k][3]), GRB.EQUAL, 0, ""); 
      double[][] t4 = {{0, 0, 1, 1, -1},
			{1, 0, -1, 0, 1},
			{-1, 0, 1, 0, 0}};
      double[] c4 = {1, 0, 0};
      for (int q = 0; q < 3; q ++)   
        model.addConstr(linExprOf(t4[q], Pi_init[0][j][k][0], Pi_init[0][j][k][2], DA[0][3][j][k][0], DA[0][3][j][k][2], Cond[j][k][3]), GRB.GREATER_EQUAL, c4[q], ""); 
      
      
      //DA1
      model.addConstr(DA[0][1][j][k][0], GRB.EQUAL, Pi_init[1][j][k][0], "");
      model.addConstr(DA[0][1][j][k][1], GRB.EQUAL, Pi_init[1][j][k][1], "");
      model.addConstr(DA[0][1][j][k][2], GRB.EQUAL, Pi_init[1][j][k][2], "");
      //DA2
      model.addConstr(DA[0][2][j][k][0], GRB.EQUAL, 1, "");
      model.addConstr(DA[0][2][j][k][1], GRB.EQUAL, 1, "");
      model.addConstr(DA[0][2][j][k][2], GRB.EQUAL, 1, "");
    }

  }  

  public void addconstr_in(GRBVar[][][] Cond) throws GRBException {
    GRBLinExpr Allconstr = new GRBLinExpr();

    for (int i = 0; i < 5; i++)  
      for (int j = 0; j < 4; j++)  
        for (int k = 0; k < 32; k++) {  
          Allconstr.addTerm(1.0, Cond[i][k][j]); 
      }
    
    model.addConstr(Allconstr, GRB.LESS_EQUAL, 400.0, "");
  }  


  public void addfivexor_red(GRBVar[][][][][] DA, GRBVar[][][][] DP, GRBVar[][][] DC1) throws GRBException {
    GRBVar[][][][] DA_Allone = new GRBVar[DP.length][5][32][3];

    
    GRBVar[][] DA_Twored = new GRBVar[5][32];
    double[] t1 = {1.0, 1.0, -1.0};

    for (int round = 0; round < DP.length; round ++) 
      for (int i = 0; i < 5; i ++) 
	for (int k = 0; k < 32; k ++) {
          for (int l = 0; l < 3; l ++) {
            DA_Allone[round][i][k][l] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "DA_Allone_"+round+"_"+i+"_"+k+"_"+l); 
            if (round == 0) {
              if (i==2) 
                model.addConstr(DA_Allone[round][i][k][l], GRB.EQUAL, 1, "");
              else
                Determine_AllOne(DA_Allone[round][i][k][l], DA[round][i][0][k][l], DA[round][i][1][k][l], DA[round][i][2][k][l], DA[round][i][3][k][l], DA[round][i][4][k][l]);           
            }
            else {
              Determine_AllOne(DA_Allone[round][i][k][l],DA[round][i][0][k][l],DA[round][i][1][k][l],DA[round][i][2][k][l],DA[round][i][3][k][l],DA[round][i][4][k][l]);
            }
          }
	  if (round==0) {
	    DA_Twored[i][k] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "DA_Twored_"+i+"_"+k); 
            if (i==2)
              model.addConstr(DA_Twored[i][k], GRB.EQUAL, 0, "");      
            else 
              Determine_twoinfive(DA_Twored[i][k],DA[round][i][0][k][0],DA[round][i][1][k][0],DA[round][i][2][k][0],DA[round][i][3][k][0],DA[round][i][4][k][0]);
          }  
          model.addConstr(DA_Allone[round][i][k][1], GRB.EQUAL, DP[round][i][k][1], "");
          model.addConstr(DA_Allone[round][i][k][1], GRB.GREATER_EQUAL, DP[round][i][k][0], "");
          model.addConstr(DA_Allone[round][i][k][0], GRB.LESS_EQUAL, DP[round][i][k][0], "");
          model.addConstr(linExprOf(t1, DC1[round][i][k], DA_Allone[round][i][k][0], DP[round][i][k][0]), GRB.EQUAL, 0, "");
          model.addConstr(DA_Allone[round][i][k][2], GRB.EQUAL, DP[round][i][k][2], "");

  
          if (round==0)
            model.addConstr(DC1[round][i][k], GRB.LESS_EQUAL, DA_Twored[i][k], "");
          
        }  
  }

 public void addtwoxor_red(GRBVar[][][][] DP2, GRBVar[][][][] DP, GRBVar[][][] DC12) throws GRBException {
    GRBVar[][][][] DP_Allone = new GRBVar[DP.length][5][32][3];
    //two red
    GRBVar[][][] DP_Allzero = new GRBVar[DP.length][5][32];
    for (int round = 0; round < DP.length; round ++)
      for (int i = 0; i < 5; i ++) 
	for (int k = 0; k < 32; k ++) {
          for (int l = 0; l < 3; l ++) {
            DP_Allone[round][i][k][l] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "DP_Allone_"+round+"_"+i+"_"+k+"_"+l); 
            Determine_AllOne(DP_Allone[round][i][k][l],DP[round][(i+4)%5][k][l],DP[round][(i+1)%5][(k+63)%32][l]);
          }

          DP_Allzero[round][i][k] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "DP_Allzero_"+round+"_"+i+"_"+k); 
          Determine_Allzero(DP_Allzero[round][i][k],DP[round][(i+4)%5][k][0],DP[round][(i+1)%5][(k+63)%32][0]);

          double[] t1 = {1.0, 1.0, -1.0};
          model.addConstr(DP_Allone[round][i][k][1], GRB.EQUAL, DP2[round][i][k][1], "");
          model.addConstr(DP_Allone[round][i][k][1], GRB.GREATER_EQUAL, DP2[round][i][k][0], "");
          model.addConstr(DP_Allone[round][i][k][0], GRB.LESS_EQUAL, DP2[round][i][k][0], "");
          model.addConstr(linExprOf(t1, DC12[round][i][k], DP_Allone[round][i][k][0], DP2[round][i][k][0]), GRB.EQUAL, 0, "");
          model.addConstr(DP_Allone[round][i][k][2], GRB.EQUAL, DP2[round][i][k][2], "");

          model.addConstr(DC12[round][i][k], GRB.LESS_EQUAL, DP_Allzero[round][i][k], "");
        }      
  }
  
  
  public void addTheta_red(GRBVar[][][][][] DA, GRBVar[][][][] DP2, GRBVar[][][][][] DB, GRBVar[][][][] DC2) throws GRBException {
    GRBVar[][][][][] DP2_Allone = new GRBVar[DP2.length][5][5][32][3];
    //two red
    GRBVar[][][][] DP2_Allzero = new GRBVar[DP2.length][5][5][32];
    for (int round = 0; round < DP2.length; round ++)
      for (int i = 0; i < 5; i ++) 
        for (int j = 0; j < 5; j ++)
	  for (int k = 0; k < 32; k ++) {
            for (int l = 0; l < 3; l ++) {
              DP2_Allone[round][i][j][k][l] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "DP2_Allone_"+round+"_"+i+"_"+j+"_"+k+"_"+l); 
              Determine_AllOne(DP2_Allone[round][i][j][k][l],DA[round][i][j][k][l], DP2[round][i][k][l]);
            }
          DP2_Allzero[round][i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "DP2_Allzero_"+round+"_"+i+"_"+j+"_"+k); 
          Determine_Allzero(DP2_Allzero[round][i][j][k],DA[round][i][j][k][0], DP2[round][i][k][0]);


          double[] t1 = {1.0, 1.0, -1.0};
          model.addConstr(DP2_Allone[round][i][j][k][1], GRB.EQUAL, DB[round][j][(2*i+3*j)%5][(k+rho[i][j])%32][1], "");
          model.addConstr(DP2_Allone[round][i][j][k][1], GRB.GREATER_EQUAL, DB[round][j][(2*i+3*j)%5][(k+rho[i][j])%32][0], "");
          model.addConstr(DP2_Allone[round][i][j][k][0], GRB.LESS_EQUAL, DB[round][j][(2*i+3*j)%5][(k+rho[i][j])%32][0], "");
          model.addConstr(linExprOf(t1, DC2[round][i][j][k], DP2_Allone[round][i][j][k][0], DB[round][j][(2*i+3*j)%5][(k+rho[i][j])%32][0]), GRB.EQUAL, 0, "");
          model.addConstr(DP2_Allone[round][i][j][k][2], GRB.EQUAL, DB[round][j][(2*i+3*j)%5][(k+rho[i][j])%32][2], "");

          
          model.addConstr(DC2[round][i][j][k], GRB.LESS_EQUAL, DP2_Allzero[round][i][j][k], "");
        } 
  }

  


  public void addSbox_nc(GRBVar[][][][][] DB, GRBVar[][][][][] DA) throws GRBException {

    double[][] t = {{0, 0, 1, 0, 0, 0, 0, 0, -1},
		    {1, 0, 0, 0, 0, 0, -1, 0, 0},
		    {1, 0, 0, 1, 0, 0, -1, -1, 1},
		    {0, 0, 1, 0, 0, 1, 1, -1, -1},
		    {0, -1, 0, -1, 1, -1, 0, 1, 0},
		    {0, 0, 0, 0, 0, 0, -1, 1, 0},
		    {0, 0, 0, 0, 0, 0, 0, 1, -1},
		    {0, 1, -1, -1, 1, 0, 1, -1, 1},
		    {-1, 1, 1, 0, 0, 0, 1, -1, 0},
		    {0, 0, -1, 0, 0, -1, 0, 0, 1},
		    {-1, 0, 0, -1, 0, 0, 1, 0, 0}};
    double[] c = {0, 0, 0, 0, -1, 0, 0, 0, 0, -1, -1};
    GRBVar[][][][][] DB_mul = new GRBVar[DB.length][5][5][32][3];
    for (int round = 0; round < DB.length; round ++)
      for (int i = 0; i < 5; i ++) 
        for (int j = 0; j < 5; j ++) 
          for (int k = 0; k < 32; k ++) {
            for (int l = 0; l < 3; l ++) {
              DB_mul[round][i][j][k][l] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, "DB_mul_"+round+"_"+i+"_"+j+"_"+k+"_"+l); 
            }
            for (int q = 0; q < 11; q++) {
              model.addConstr(linExprOf(t[q], DB[round][(i+1)%5][j][k][0], DB[round][(i+1)%5][j][k][1], DB[round][(i+1)%5][j][k][2], DB[round][(i+2)%5][j][k][0], DB[round][(i+2)%5][j][k][1], DB[round][(i+2)%5][j][k][2], DB_mul[round][i][j][k][0], DB_mul[round][i][j][k][1], DB_mul[round][i][j][k][2]), GRB.GREATER_EQUAL, c[q], "");
     }
    }

    for (int round = 0; round < DB.length; round ++)
      for (int i = 0; i < 5; i ++) 
	for (int j = 0; j < 5; j ++) 
	  for (int k = 0; k < 32; k ++) 
            for (int l = 0; l < 3; l ++) {
              Determine_AllOne(DA[round+1][i][j][k][l],DB[round][i][j][k][l],DB[round][(i+2)%5][j][k][l],DB_mul[round][i][j][k][l]);         
        } 
  }



  public void addDoM(GRBVar[][][][][] DA, GRBVar[][] dom) throws GRBException {
      
    int r=DA.length-1;
    for (int k = 0; k < 32; k ++) {
      Determine_AllOne(dom[0][k],DA[r][3][0][(k-rho[3][0]+32)%32][1],DA[r][3][3][(k-rho[3][0]+32)%32][1],DA[r][0][2][(k-rho[0][2]+32)%32][1],DA[r][0][0][(k-rho[0][2]+32)%32][1]);
      Determine_AllOne(dom[1][k],DA[r][4][1][(k-rho[4][1]+64)%32][1],DA[r][4][4][(k-rho[4][1]+64)%32][1],DA[r][1][3][(k-rho[1][3]+64)%32][1],DA[r][1][1][(k-rho[1][3]+64)%32][1]);
    } 


  }

  public void betaConstraints(GRBVar[][][][] Pi_init, GRBVar[][] beta) throws GRBException {
    double[] t = {1, 1, -2};
    for (int k = 0; k < 32; k ++) {
      model.addConstr(Pi_init[0][0][k][0], GRB.GREATER_EQUAL, beta[0][k], "");
      model.addConstr(Pi_init[0][0][k][2], GRB.GREATER_EQUAL, beta[0][k], "");     
      model.addConstr(linExprOf(t, Pi_init[0][0][k][0], Pi_init[0][0][k][2], beta[0][k]), GRB.LESS_EQUAL, 1, "");
      model.addConstr(Pi_init[0][1][k][0], GRB.GREATER_EQUAL, beta[1][k], "");
      model.addConstr(Pi_init[0][1][k][2], GRB.GREATER_EQUAL, beta[1][k], "");     
      model.addConstr(linExprOf(t, Pi_init[0][1][k][0], Pi_init[0][1][k][2], beta[1][k]), GRB.LESS_EQUAL, 1, "");
      model.addConstr(Pi_init[0][2][k][0], GRB.GREATER_EQUAL, beta[2][k], "");
      model.addConstr(Pi_init[0][2][k][2], GRB.GREATER_EQUAL, beta[2][k], "");     
      model.addConstr(linExprOf(t, Pi_init[0][2][k][0], Pi_init[0][2][k][2], beta[2][k]), GRB.LESS_EQUAL, 1, "");
      model.addConstr(Pi_init[0][4][k][0], GRB.GREATER_EQUAL, beta[3][k], "");
      model.addConstr(Pi_init[0][4][k][2], GRB.GREATER_EQUAL, beta[3][k], "");     
      model.addConstr(linExprOf(t, Pi_init[0][4][k][0], Pi_init[0][4][k][2], beta[3][k]), GRB.LESS_EQUAL, 1, "");
    }
  }

  public void Determine_Allzero(GRBVar mainVar, GRBVar ... vars) throws GRBException {
        GRBLinExpr expr = new GRBLinExpr();
	expr.addTerm(1.0, mainVar);
        for (int i = 0; i < vars.length; i++) {
          expr.addTerm(1.0, vars[i]);
	  model.addConstr(linExprOf(mainVar,vars[i]), GRB.LESS_EQUAL, 1, "");
        }
        model.addConstr(expr, GRB.GREATER_EQUAL, 1, "");
  }

  public void Determine_twointhree(GRBVar mainVar, GRBVar ... vars) throws GRBException {
        GRBLinExpr expr = new GRBLinExpr();
	expr.addTerm(2.0, mainVar);
        for (int i = 0; i < vars.length; i++) {
          expr.addTerm(1.0, vars[i]);
        }
        model.addConstr(expr, GRB.GREATER_EQUAL, 2, "");
        model.addConstr(expr, GRB.LESS_EQUAL, 3, "");
  }

  public void Determine_twoinfive(GRBVar mainVar, GRBVar ... vars) throws GRBException {
        GRBLinExpr expr1 = new GRBLinExpr();
	GRBLinExpr expr2 = new GRBLinExpr();
	expr1.addTerm(4.0, mainVar);
	expr2.addTerm(-2.0, mainVar);
        for (int i = 0; i < vars.length; i++) {
          expr1.addTerm(1.0, vars[i]);
          expr2.addTerm(-1.0, vars[i]);
        }
        model.addConstr(expr1, GRB.GREATER_EQUAL, 4, "");
        model.addConstr(expr2, GRB.GREATER_EQUAL, -5, "");
  }

  public void Determine_twoinfour(GRBVar mainVar, GRBVar ... vars) throws GRBException {
        GRBLinExpr expr1 = new GRBLinExpr();
	GRBLinExpr expr2 = new GRBLinExpr();
	expr1.addTerm(3.0, mainVar);
	expr2.addTerm(-2.0, mainVar);
        for (int i = 0; i < vars.length; i++) {
          expr1.addTerm(1.0, vars[i]);
          expr2.addTerm(-1.0, vars[i]);
        }
        model.addConstr(expr1, GRB.GREATER_EQUAL, 3, "");
        model.addConstr(expr2, GRB.GREATER_EQUAL, -4, "");
  }

  public void Determine_AllOne(GRBVar mainVar, GRBVar ... vars) throws GRBException {
        GRBLinExpr exprm = new GRBLinExpr();
	exprm.addTerm(1.0, mainVar);
	GRBLinExpr exprp = new GRBLinExpr();
	exprp.addTerm(-1.0*vars.length, mainVar);
        for (int i = 0; i < vars.length; i++) {
          exprm.addTerm(-1.0, vars[i]);
	  exprp.addTerm(1.0, vars[i]);
        }
        model.addConstr(exprm, GRB.GREATER_EQUAL, 1-vars.length, "");
	model.addConstr(exprp, GRB.GREATER_EQUAL, 0, "");
  }
  
  public void Determine_ExistOne(GRBVar mainVar, GRBVar ... vars) throws GRBException {
        GRBLinExpr expr1 = new GRBLinExpr();
	expr1.addTerm(1.0, mainVar);
        GRBLinExpr exprm = new GRBLinExpr();
        exprm.addTerm(vars.length, mainVar);
        for (int i = 0; i < vars.length; i++) {
          expr1.addTerm(-1.0, vars[i]);
          exprm.addTerm(-1.0, vars[i]);
        }
        model.addConstr(exprm, GRB.GREATER_EQUAL, 0, "");
	model.addConstr(expr1, GRB.LESS_EQUAL, 0, "");
  }

  public void Determine_Existzero(GRBVar mainVar, GRBVar ... vars) throws GRBException {
        GRBLinExpr expr = new GRBLinExpr();
	expr.addTerm(1.0, mainVar);
        for (int i = 0; i < vars.length; i++) {
          expr.addTerm(1.0, vars[i]);
          model.addConstr(linExprOf(mainVar,vars[i]), GRB.GREATER_EQUAL, 1, "");
        }
	model.addConstr(expr, GRB.LESS_EQUAL, vars.length, "");
  }

  public GRBLinExpr linExprOf(double[] coeffs, GRBVar ... vars) throws GRBException {
    GRBLinExpr ofVars = new GRBLinExpr();
    ofVars.addTerms(coeffs, vars);
    return ofVars;
  }

  public GRBLinExpr linExprOf(double constant, GRBVar ... vars) {
    GRBLinExpr ofVars = linExprOf(vars);
    ofVars.addConstant(constant);
    return ofVars;
  }

  public GRBLinExpr linExprOf(GRBVar ... vars) {
    GRBLinExpr expr = new GRBLinExpr();
    for (GRBVar var : vars)
      expr.addTerm(1.0, var);
    return expr;
  }
 
}
