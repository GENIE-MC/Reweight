#include "ProfSpline/ParamPoints.h"
#include <algorithm>

namespace Professor {

  using namespace std;


  ParamPoints::ParamPoints(const vector< vector<double> >& ppoints) {
    /// @todo Throw a ParamPointsError or similar rather than this assert
    assert(!ppoints.empty());
    // _parampoints.clear();
    // _locked = false;
    // for (size_t i = 0; i < p.size(); ++i) {
    //   _parampoints.push_back(p[i]);
    // }
    _parampoints = ppoints;
    _locked = true;
  }


  vector<double> ParamPoints::ptcenters() const {
    vector<double> temp_max, temp_min;
    for (size_t i = 0; i < dim(); i++) { // iteration over coordinates
      vector<double> temp;
      for (size_t j = 0; j < numPoints(); j++) { // iteration over anchors
        temp.push_back(_parampoints[j][i]);
      }
      temp_max.push_back(*max_element(temp.begin(), temp.end()));
      temp_min.push_back(*min_element(temp.begin(), temp.end()));
    }
    vector<double> center;
    for (size_t i = 0; i < dim(); i++) { // iteration over coordinates
      center.push_back(temp_min[i] + 0.5* (temp_max[i] - temp_min[i]));
    }

    return center;
  }


  vector<double> ParamPoints::ptmins() const {
    vector<double> temp_min;
    for (size_t i = 0; i < dim(); i++) { // iteration over coordinates
      vector<double> temp;
      for (size_t j = 0; j < numPoints(); j++) { // iteration over anchors
        temp.push_back(_parampoints[j][i]);
      }
      temp_min.push_back(*min_element(temp.begin(), temp.end()));
    }
    return temp_min;
  }


  vector<double> ParamPoints::ptmaxs() const {
    vector<double> temp_max;
    for (size_t i = 0; i < dim(); i++) { // iteration over coordinates
      vector<double> temp;
      for (size_t j = 0; j < numPoints(); j++) { // iteration over anchors
        temp.push_back(_parampoints[j][i]);
      }
      temp_max.push_back(*max_element(temp.begin(), temp.end()));
    }
    return temp_max;
  }


  vector< pair<double, double> > ParamPoints::ptedges() const {
    const vector<double> mins = ptmins();
    const vector<double> maxs = ptmaxs();
    vector< pair<double, double> > edge;
    edge.reserve(dim());
    for (size_t i = 0; i < dim(); i++) {
      edge.push_back( pair<double, double>(mins[i], maxs[i]));
    }
    return edge;
  }

  void ParamPoints::setNames(std::vector<std::string > names) {
    if (_names.size() == 0) { // No names set so far
      if (dim() == names.size()) { // Sanity check
        for (size_t i = 0; i < names.size(); i++) {
          _names.push_back(names[i]);
        }
      }
      else {
        stringstream ss;
        ss << "ParamPoints::setNames: dimension mismatch (" << dim() << "dimensions vs. " << names.size() << " names)  ";
        throw ParamPointsError(ss.str());
      }
    }
    else { // Names already set
      stringstream ss;
      ss << "ParamPoints::setNames: Names already set!";
      throw ParamPointsError(ss.str());
    }
  }

  void ParamPoints::printMeta() const {
    cout << "Nr. of points: " << numPoints() << endl;
    cout << "Dimension:     " << dim() << endl;
    cout << "Center:       ";
    for (size_t i = 0; i < dim(); i++) {
      cout << " " << ptcenters()[i];
    }
    cout << endl;
    cout << "Edges:" << endl;
    const vector< pair<double, double> > edges = ptedges();
    for (size_t i = 0; i < dim(); i++) {
      cout << edges[i].first << " < " << edges[i].second << endl;
    }
    cout << endl;
  }


  void ParamPoints::printPoints() const {
    for (size_t i = 0; i < numPoints(); ++i) {
      cout << "Point " << i << ":" << endl;
      for (size_t j = 0; j < dim(); ++j) {
        cout << _parampoints[i][j] << " ";
      }
      cout << endl;
    }
  }


}
