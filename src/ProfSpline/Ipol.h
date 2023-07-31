#ifndef PROF_IPOL_H
#define PROF_IPOL_H

#include "ProfSpline/ParamPoints.h"
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <stdexcept>

namespace Professor {


  /// Throwable error
  struct IpolError : public std::runtime_error {
    IpolError(const std::string& reason) : std::runtime_error(reason) { }
  };


  /// @name Calculator functions for parameterisation elements
  //@{

  /// Calculate the number of coefficients for the given parameter space dimension and polynomial order
  /// @todo Deal in uints
  int calcnumCoeffs(int dim, int order);

  /// Calculate parametrisation coefficients
  /// @note structure is the pre-calculated algebraic structure of the polynomial
  /// @todo Provide an (unscaled) arrays-only version with a vector<vector<double>> in place of ParamPoints
  std::vector<double> calcCoeffs(const ParamPoints& pts, const std::vector<double>& vals, int order,
                                 double threshold, const std::vector<std::vector<int> >& structure);

  /// Calculate an interpolated value
  double calcValue(const std::vector<double>& params,
                   const std::vector<double>& coeffs, int order,
                   const std::vector<std::vector<int> >& structure);

  /// Calculate an interpolated value
  double calcValue(const std::vector<double>& paramslongvector,
                   const std::vector<double>& coeffs);

  /// @brief Make the algebraic coefficient structure
  ///
  /// The structure is a nested vector of ints, with each entry in the outer vector
  /// representing a term in the polynomial and the inner vectors being the powers
  /// to which each param should be raised in that term.
  ///
  /// @note In a given param space, the structure can be reused between any values
  ///   or errors parameterised at the same polynomial order. So pre-computing
  ///   structures at 3rd, 4th and 5th order would cover everything needed...
  std::vector< std::vector<int> > mkStructure(int dim, int order);

  /// Make the vector of polynomial terms to which the coeffs are to be applied, at the given order
  ///
  /// @note The same long vector can be used in any parameterised value calculation
  ///    in the "active" parameter space, between values & errors and multiple bins.
  ///    Just need to "dot" it with the appropriate fixed coefficient vector.
  std::vector<double> mkLongVector(const std::vector<double>& params, int order,
                                   const std::vector< std::vector<int> >& structure);

  // vector<double> mkLongVectorDerivative(const vector<double>& params, int order,
  //                                       const vector<double>& minPV, const vector<double>& maxPV,
  //                                       const vector<vector<int> >& structure);

  // vector<double> mkLongVectorGradient(const vector<double>& params, int coord, int order,
  //                                     const vector<double>& minPV, const vector<double>& maxPV,
  //                                     const vector<vector<int> >& structure);

  //@}



  /// The heart of Professor: the interpolation of a single numerical value through the parameter space
  class Ipol {
  public:

    /// @brief Constructor for calculation of coefficients
    ///
    /// The @a pts list of N-dimensional parameter points must correspond to the @a ptvals
    /// list of values at those points, to be interpolated at the given polynomial @a order.
    /// A name may optionally be given.
    ///
    /// @note Expert settings: The stability of the SVD operation is controlled
    /// by the @a svdthreshold parameter, which should not normally be
    /// touched. The stability is normally ensured by internally scaling
    /// parameter points into unit ranges within the sampled hypercube defined
    /// by @a pts; changing @doscaling to false will disable this scaling, which
    /// simplifies Ipol I/O (no PMin/Max metadata is needed) but risks SVD
    /// instability.
    ///
    Ipol(const ParamPoints& pts, const std::vector<double>& ptvals, int order,
         const std::string& name="", double svdthreshold=1e-20, bool doscaling=true) {
      _dim = pts.dim();
      _order = order;
      _name = name;
      _structure = mkStructure(_dim, _order);
      if (doscaling) {
        _minPV = pts.ptmins();
        _maxPV = pts.ptmaxs();
      }
      _coeffs = calcCoeffs(pts, ptvals, _order, svdthreshold, _structure);
    };

    /// Constructor to read ipol from file (one string for each object)
    /// @todo Also allow optional passing of pmins, pmaxs vectors for the case where the string includes scaling?
    Ipol(const std::string& s) {
      fromString(s);
    };

    ~Ipol() {
      _coeffs.clear();
      _structure.clear();
    };


    /// @name String representations
    //@{

    /// Get a persistent string representation of this Ipol
    std::string toString(const std::string& name="") const;

    /// Read and set coefficients (name), order from string
    /// @todo Also allow optional passing of pmins, pmaxs vectors for the case where the string includes scaling?
    void fromString(const std::string& s);

    /// Get a string representation of the mathematical expression
    std::string exprString() const;

    //@}


    /// @name Calculations
    //@{

    /// Get the value of the parametrisation at point p
    double value(const std::vector<double>& p) const;

    /// Get the value of the derivative of the parametrisation at point p
    /// @todo Expose as a standalone calcDerivative function, cf. calcValue
    double derivative(const std::vector<double>& p) const;

    /// Get the gradient of the parametrisation at point p
    /// @todo Expose as a standalone calcGradient function, cf. calcValue
    std::vector<double> gradient(const std::vector<double>& p) const;

    //@}


    /// @name Coefficient access
    //@{

    /// Get the number of coefficients
    const size_t numCoeffs() const {
      // ::numCoeffs(dim(),order())
      return _coeffs.size();
    }

    /// Get the vector of coefficients by const reference
    const std::vector<double>& coeffs() const { return _coeffs; }

    /// Get a single coefficient
    const double& coeff(size_t i) const { return coeffs()[i]; }

    //@}


    /// @name Polynomial term structure
    //@{

    /// Get the polynomial term exponent structure
    const std::vector< std::vector<int> >& structure() const { return _structure; }

    /// Convert params to scaled params
    std::vector<double> sparams(const std::vector<double>& params) const;

    /// Get a long vector of polynomial terms (sans coefficients) for the given param point
    std::vector<double> longVector(const std::vector<double>& params) const {
      return mkLongVector(sparams(params), order(), structure());
    }

    //@}


    /// @name Basic ipol properties
    //@{

    /// Accessors to the dimension of the param points
    int dim() const { return _dim; }
    int numParams() const { return dim(); }

    /// Get the order of the parametrisation
    int order() const { return _order; }

    /// Get the name of the parametrised object
    std::string name() const { return _name; }

    //@}


    /// @name Limit-setting
    //@{

    void setParamLimits(const std::vector<double>& minpvs, const std::vector<double>& maxpvs) {
      setMinParamVals(minpvs);
      setMaxParamVals(maxpvs);
    }

    const std::vector<double>& minParamVals() { return _minPV; }
    const std::vector<double>& maxParamVals() { return _maxPV; }

    // The if statements here guarantee backwards compatibility with the format 'binned 2' where the
    // ranges were read from the ipol meta. In binned 3 they are read for each ipol from the coefficient
    // output
    void setMinParamVals(const std::vector<double>& minpvs) { if (_minPV.empty()) _minPV = minpvs; }
    void setMaxParamVals(const std::vector<double>& maxpvs) { if (_maxPV.empty()) _maxPV = maxpvs; }

    //@}


  private:

    int _dim, _order;
    std::vector<std::vector<int> > _structure;
    std::string _name;
    std::vector<double> _coeffs, _minPV, _maxPV;

  };


}

#endif
