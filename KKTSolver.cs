using QP.solvers;
using System;
using MathNet.Numerics.LinearAlgebra;

namespace solvers.solvers
{
   public abstract class KKTSolver : IKKTSolver
   {
      Vector<double> _h;
      Vector<double> _g;
      Matrix<double> _A;
      Matrix<double> _Q;
      double _tol = 1.0E+10 * Double.Epsilon;

     public  Vector<double> GVector
      {
         get
         {
            return this._g;
         }

         set
         {
            this._g = value;
         }
      }

     public  Matrix<double> AMatrix
      {
         get
         {
            return this._A;
         }

         set
         {
            this._A = value;
         }
      }

      public Vector<double> HVector
      {
         get
         {
            return this._h;
         }

         set
         {
            this._h = value;
         }
      }

     public  Matrix<double> QMatrix
      {
         get
         {
            return this._Q;
         }

         set
         {
            this._Q = value;
         }
      }

      public double Tolerance
      {
         get
         {
            return _tol;
         }

         set
         {
            this._tol = value;
         }
      }

      public abstract Vector<double> Solve();
      public Vector<double> residual(Vector<double> sol)
      {
         int n = this.QMatrix.ColumnCount;
         int m = this.AMatrix.RowCount;
         var x = sol.SubVector(0, n);
         var y = sol.SubVector(n, m);
         var res = CreateVector.DenseOfArray<double>(new double[n + m]);

         res.SetSubVector(0, n, QMatrix.Multiply(x) + AMatrix.Transpose().Multiply(y));
         res.SetSubVector(n, m, AMatrix.Multiply(x));

         return res;
      }
     
   }
}
