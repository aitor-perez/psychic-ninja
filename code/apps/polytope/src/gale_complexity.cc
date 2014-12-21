/* Copyright (c) 2014
   Julian Pfeifle & Discrete Geometry class FME/UPC 2014

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version: http://www.gnu.org/licenses/gpl.txt.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
--------------------------------------------------------------------------------
*/

#include "polymake/SparseVector.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/ListMatrix.h"
#include "polymake/Array.h"
//#include "polymake/linalg.h"

namespace polymake { namespace polytope {

class Configuration {
  public:
    int e;
    int n;
    Vector<Rational> v;
    Matrix<Rational> m;
    Matrix<Rational> cocircuits;
    
    Configuration(int _e, int _n, const Vector<Rational>& _v) {
      e = _e;
      n = _n;
      v = _v;
      
      m = Matrix<Rational>(e+1, n);
      for (int j=0; j<n; ++j) {
        m[0][j] = 1;
        for (int i=1; i<e+1; ++i) {
          m[i][j] = v[j*e + i];
        }
      }
      
      // Generar cocircuits
      
    }
};

Matrix<Rational> enumerate_configurations(int e, int n, int m) {
  Matrix<Rational> equations(e+1, e*n+1);
  for (int j=0; j<e; ++j) {
    for (int i=0; i<n; ++i) {
      equations(j, 1 + i*e + j) = 1;
    }
  }
  
  // Adding last equation v_11 = m
  equations(e, 0) = -m;
  equations(e, 1) = 1;
  
  ListMatrix<SparseVector<Rational> > inequalities(0, e*n+1);
  
  // Symmetries
  for (int i=0; i<e-1; ++i) { // leave out last row
    SparseVector<Rational> ineq(e*n+1);
    ineq[1 + i] = 1;
    ineq[1 + i + 1] = -1;
    inequalities /= ineq;
  }
  SparseVector<Rational> last(e*n+1);
  last[e] = 1;
  inequalities /= last;
  
  // (Almost) Lexicographic order
  for (int i=0; i<n-1; ++i) { // leave out last row, since it makes no sense
    SparseVector<Rational> ineq(e*n+1);
    ineq[1 + i*e] = 1;
    ineq[1 + (i+1)*e] = -1;
    inequalities /= ineq;
  }
  
  // Bounding box inequalities
  for (int i=1; i<e*n+1; ++i) {
    SparseVector<Rational> ineq(e*n+1);
    ineq[0] = m;
    ineq[i] = 1;
    inequalities /= ineq;
    ineq[i] = -1;
    inequalities /= ineq;
  }
  
  perl::Object P("Polytope");
  P.take("EQUATIONS") << equations;
  P.take("INEQUALITIES") << inequalities;
  
  const Matrix<Rational> C = P.CallPolymakeMethod("LATTICE_POINTS"); // S'ha de fer així

  return C;
}

bool check_lexicographically_sorted(int e, const Vector<Rational>& v) {
  int n = (v.size()-1)/e;

  int prev = 1;
  for (int i=1; i<n; i++) {
    int next = 1 + e*i;
    
    for (int j=0; j<e; ++j) {
      if (v[prev + j] > v[next + j]) {
        break;
      }
      
      if (v[prev + j] < v[next + j]) {
        return false;
      }
    }
    
    prev = next;
  }
  
  return true;
}

Array<Matrix<Rational> > gale_complexity(int e, int n, int m) {
  Matrix<Rational> configs = enumerate_configurations(e,n,m);
  
  int pos = 0;
  int neg = 0;
  for (int k=0; k<configs.rows(); ++k) {
    if (check_lexicographically_sorted(e, configs[k])) {
      ++pos;
      
      Configuration conf(e, n, configs[k]);
      
      cerr << conf.v << endl;
      cerr << "***" << endl;
      cerr << conf.m << endl;
      cerr << "***" << endl;
      
//      Matrix<Rational> kern = T(null_space(T(conf.m)));
      
  //    cerr << kern << endl;
      break;
    } else {
      ++neg;
    }
  }
  
  cerr << pos << "vs" << neg << endl;
  cerr << "Total: " << pos+neg << endl;
  
  return Array<Matrix<Rational> > (0);
}

UserFunctionTemplate4perl("# @category Computations"
                          "# Computes the lattice points of all polytopes with Gale complexity m."
                          "# @param Int e Dimension of the Gale diagram vector space"
                          "# @param Int n Number of vectors of the Gale diagram (i.e. number of vertices of the polytopes)"
                          "# @param Int m Gale complexity of the polytopes"
                          "# @return Array<Matrix<Rational> > Array of all the lattice points for all the polytopes",
                  			  "gale_complexity(Int, Int, Int)" );

} }

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
