#ifndef PiecewiseCubic_SPLINE_H
#define PiecewiseCubic_SPLINE_H

#include <vector>
#include <iostream>

struct CubicSegment {
    double x0; // interval start
    double x1; // interval end
    std::vector<double> values;
    double xc; // center of interval
};

// Evaluate cubic at x
inline double EvalCubicSegment(const CubicSegment& seg, double x){
    double dx = x - seg.xc;
    return seg.values[0] + seg.values[1]*dx + seg.values[2]*dx*dx + seg.values[3]*dx*dx*dx;
}

class PiecewiseCubicSpline {
public:

    std::vector<CubicSegment> _segs;

    const std::vector<CubicSegment> GetSegments(){
        return _segs;
    }

    void SetSegments(const std::vector<CubicSegment>& s){
        _segs = s;
    }

    void SetParameter(const int n_segment, const int n_parameter, const double value){
         _segs[n_segment].values[n_parameter] = value;
    }

    // default constructor
    PiecewiseCubicSpline() = default;

    PiecewiseCubicSpline(const std::vector<CubicSegment>& s)
    {
        SetSegments(s);
    }

    PiecewiseCubicSpline(const char* file="coeffs_global.txt")
    {

        std::ifstream in(file);
        if(!in.is_open()){
            std::cerr << "Cannot open file: " << file << std::endl;
            return;
        }

        std::string line;
        std::getline(in, line); // skip header

        while(std::getline(in, line)){
            if(line.empty() || line[0]=='#') continue;

            std::istringstream ss(line);

            CubicSegment s;
            ss >> s.xc >> s.x0 >> s.x1
               >> s.values[0] >> s.values[1] >> s.values[2] >> s.values[3];

            if(ss.fail()) continue;

            _segs.push_back(s);
        }

        in.close();

        std::cout << "Loaded segments: " << _segs.size() << std::endl;
    }


    double Eval(double x) const {
        if(_segs.empty()) return 0.0;

        // find interval (linear search, same as your code)
        size_t k = 0;
        for(; k < _segs.size(); ++k){
            if(x >= _segs[k].x0 && x <= _segs[k].x1)
                break;
        }

        if(k == _segs.size()) k = _segs.size() - 1;
        return EvalCubicSegment(_segs[k], x);
    }

    size_t Size() const { return _segs.size(); }
};

#endif
