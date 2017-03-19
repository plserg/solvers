
using MathNet.Numerics.LinearAlgebra;

namespace solvers.solvers
{
   /**
 * Solves
 * 
 * Q.x + [A]T.y = -g, <br>
 * A.x = -h
 * 
 * for diagonal Q.
 * H is expected to be diagonal.
 * Only the subdiagonal elements are relevant.
 * 
 * @see "S.Boyd and L.Vandenberghe, Convex Optimization, p. 542"
 * */
   public class DiagonalKKTSolver : KKTSolver
   {
      public override Vector<double> Solve()
      {
         
         int n = this.QMatrix.ColumnCount;
         int m = this.AMatrix.RowCount;
         var A = this.AMatrix;
         var AT = this.AMatrix.Transpose();
         var h = this.HVector;
         var g = this.GVector;
         //shortcut for solfing the diagonal system Q:

         var Qinv = CreateMatrix.DenseDiagonal<double>(n, n, (i) => 1.00 / this.QMatrix[i,i]);
         var AQinvAT = A.Multiply( Qinv.Multiply( AT ));//Schur compliment  matrix, S = -AQ^{-1}A^T
         var q = A.Multiply(Qinv.Multiply(g)) - h;
         //parts of the solution
         var y = AQinvAT.Solve(q);
         var x = -g - AT.Multiply(y);
         //solve the diagonal system
         for (int i = 0; i < n; ++i)
         {
            x[i] *=  Qinv[i, i];
         }

         var sol = CreateVector.Dense<double>(n+m);
         for(int i = 0; i < n; ++i)
         {
            sol[i] = x[i];
         }
         for(int i = 0; i < m; ++i)
         {
            sol[n + i] = y[i];
         }
         return sol;
      }
   }
}
