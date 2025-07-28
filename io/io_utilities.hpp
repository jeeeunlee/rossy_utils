#pragma once

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <algorithm>
#include <list>
#include <Eigen/Dense>
#include <vector>
#include <string_view>
#include <yaml/yaml.h>
#include "rossy_utils/math/typedefs.h"
#include "Configuration.h"

#define ROSSY_DEBUG 1  // or 1
#define ROSSY_PRINT 1  // or 1

enum myColor {
    Red=0,
    BoldRed=1,
    Green=2,
    BoldGreen=3,
    Yellow=4,
    BoldYellow=5,
    Blue=6,
    BoldBlue=7,
    Magneta=8,
    BoldMagneta=9,
    Cyan=10,
    BoldCyan=11
};

namespace rossy_utils
{

    const std::string border = "================================================================================";
    
    // =========================================================================
    // Save Vector
    // =========================================================================
    static std::list< std::string > gs_fileName_string; //global & static
    void cleaningFile(std::string _file_name, std::string& _ret_file, bool b_param);
    
    template <typename Scalar>
    void saveVector(const VectorX<Scalar>& vec_, std::string name_, bool b_param=false) {
        std::string file_name;
        cleaningFile(name_, file_name, b_param);

        std::ofstream savefile(file_name.c_str(), std::ios::app);
        for (int i(0); i < vec_.rows(); ++i) {
            savefile << vec_(i) << "\t";
        }
        savefile << "\n";
        savefile.flush();
    };

    template <typename Scalar>
    void saveVector(const Quaternion<Scalar>& qq, std::string name_, bool b_param=false) {
        std::string file_name;
        cleaningFile(name_, file_name, b_param);

        std::ofstream savefile(file_name.c_str(), std::ios::app);
        
        savefile << qq.w() << "\t" << qq.x() << "\t"
                << qq.y() << "\t" << qq.z() << "\t";  
        
        savefile << "\n";
        savefile.flush();
    };

    template <typename Scalar>
    void saveMatrix(const MatrixX<Scalar>& mtx_, std::string name_, bool b_param=false) {
        for(int j=0; j<mtx_.rows(); ++j)  {
            saveVector<Scalar>(mtx_.row(j), name_ , b_param);
        }
    };

    template <typename Scalar>
    void saveValue(Scalar _value, std::string _name, bool b_param=false) {
        std::string file_name;
        cleaningFile(_name, file_name, b_param);
        std::ofstream savefile(file_name.c_str(), std::ios::app);

        savefile << _value << "\n";
        savefile.flush();
    };

    template <typename Scalar>
    void saveVector(Scalar* _vec, std::string _name, int size, bool b_param=false) {
        std::string file_name;
        cleaningFile(_name, file_name, b_param);
        std::ofstream savefile(file_name.c_str(), std::ios::app);

        for (int i(0); i < size; ++i) {
            savefile << _vec[i] << "\t";
        }
        savefile << "\n";
        savefile.flush();
    };

    template <typename Scalar>
    void saveVector(const std::vector<Scalar>& _vec, std::string _name,
                    bool b_param=false) {
        std::string file_name;
        cleaningFile(_name, file_name, b_param);
        std::ofstream savefile(file_name.c_str(), std::ios::app);
        for (int i(0); i < _vec.size(); ++i) {
            savefile << _vec[i] << "\t";
        }
        savefile << "\n";
        savefile.flush();
    };

    // =========================================================================
    // Read File
    // =========================================================================
    void readFile(std::string file_name_, std::vector<std::string> & _vec);
    void splitString(std::string* str_array, std::string strTarget, std::string strTok );
    
    template <typename YamlType>
    YamlType readParameter(const YAML::Node& node, const std::string& name) {
        try { return node[name.c_str()].as<YamlType>(); }
        catch (...) { throw std::runtime_error(name); }
    };

    template <typename YamlType>
    void readParameter(const YAML::Node& node, const std::string& name, YamlType& parameter) {
        try { parameter = readParameter<YamlType>(node, name); }
        catch (...) { throw std::runtime_error(name); }
    };

    // =========================================================================
    // Pretty Print
    // =========================================================================
    template <typename Scalar>
    void size_print(VectorX<Scalar> const & vv, std::ostream & os,
        std::string const & title, std::string const & prefix="", bool nonl=false){
        char const* nlornot("\n");
        if (nonl) {
            nlornot = "";
        }
        if (!title.empty()) {
            os << title <<"(length: " <<vv.size() <<")"<< nlornot;
        }
    };

    template <typename Scalar>
    void size_print(MatrixX<Scalar> const & vv, std::ostream & os,
        std::string const & title, std::string const & prefix="", bool nonl=false){
        char const* nlornot("\n");
        if (nonl) {
            nlornot = "";
        }
        if (!title.empty()) {
            os << title <<"(rows: " <<vv.rows() <<", cols: "<<vv.cols()<<")"<< nlornot;
        }
    };

    template <typename Scalar>
    std::string pretty_string(Scalar vv) {
        static int const buflen(32);
        static char buf[buflen];
        memset(buf, 0, sizeof(buf));
        snprintf(buf, buflen - 1, "% 6.6f  ", vv);
        std::string str(buf);
        return str;
    };

    template <typename Scalar>
    void pretty_print(const MatrixX<Scalar>& mm, std::ostream& os,
                    const std::string& title, const std::string& prefix="",
                    bool vecmode=false, bool nonl=false) {
        char const* nlornot("\n");
        if (nonl) {
            nlornot = "";
        }
        if (!title.empty()) {
            os << title <<"(rows: " <<mm.rows() <<", cols: "<<mm.cols()<<")"<< nlornot;
        }
        if ((mm.rows() <= 0) || (mm.cols() <= 0)) {
            os << prefix << " (empty)" << nlornot;
        } else {
            // if (mm.cols() == 1) {
            //   vecmode = true;
            // }

            if (vecmode) {
                if (!prefix.empty()) os << prefix;
                for (int ir(0); ir < mm.rows(); ++ir) {
                    os << pretty_string(mm.coeff(ir, 0));
                }
                os << nlornot;

            } else {
                for (int ir(0); ir < mm.rows(); ++ir) {
                    if (!prefix.empty()) os << prefix;
                    for (int ic(0); ic < mm.cols(); ++ic) {
                        os << pretty_string(mm.coeff(ir, ic));
                    }
                    os << nlornot;
                }
            }
        }
    };

    template <typename Scalar>
    void pretty_print(VectorX<Scalar> const& vv, std::ostream& os,
                    std::string const& title, std::string const& prefix="",
                    bool nonl=false) {
        pretty_print((MatrixX<Scalar> const&)vv, os, title, prefix, true, nonl);
    };

    template <typename Scalar>
    void pretty_print(Vector3<Scalar> const& vv, std::ostream& os,
                    std::string const& title, std::string const& prefix="",
                    bool nonl=false) {
        pretty_print((MatrixX<Scalar> const&)vv, os, title, prefix, true, nonl);
    };

    template <typename Scalar>
    void pretty_print(Eigen::Quaternion<Scalar> const& qq, std::ostream& os,
                    std::string const& title, std::string const& prefix="",
                    bool nonl=false) {
        VectorX<Scalar> wxyz = VectorX<Scalar>::Zero(4);
        wxyz << qq.w(), qq.x(), qq.y(), qq.z();
        pretty_print(wxyz, os, title, prefix, true, nonl);
        // pretty_print(qq.coeffs(), os, title, prefix, true, nonl);
    };

    template <typename Scalar>
    void pretty_print(const std::vector<Scalar>& _vec, const char* title) {
        std::printf("%s: ", title);
        for (int i(0); i < _vec.size(); ++i) {
            std::printf("% 6.4f, \t", _vec[i]);
        }
        std::printf("\n");
    };

    template <typename Scalar>
    std::string pretty_string(VectorX<Scalar> const& vv) {
        std::ostringstream os;
        pretty_print(vv, os, "", "", true);
        return os.str();
    };

    template <typename Scalar>
    std::string pretty_string(MatrixX<Scalar> const& mm,
                            std::string const& prefix) {
        std::ostringstream os;
        pretty_print(mm, os, "", prefix);
        return os.str();
    }

    void pretty_constructor(const int& _num_tab, const std::string& _name);
    void color_print(const myColor & _color, const std::string& _name, bool line_change=true);

    // =========================================================================

    inline std::ostream& debug_stream() {
    #if ROSSY_DEBUG
        return std::cout;
    #else
        static std::ofstream null_stream("/dev/null");  // discards output on Linux
        return null_stream;
    #endif
    }

    #if ROSSY_PRINT

    // -------------------------------
    // Eigen Matrix or Vector
    // -------------------------------
    template<typename Derived,
    typename std::enable_if<std::is_base_of<Eigen::EigenBase<Derived>, Derived>::value>::type>
    inline void save_to_bin(const std::string& filename, const Derived& data) {
        std::string full_path = LOG_DIR + filename;
        std::ofstream file(full_path, std::ios::binary);
        file.write(reinterpret_cast<const char*>(data.derived().data()), data.size() * sizeof(float));
        file.close();
    }

    // -------------------------------
    // std::vector<T> (float or uint)
    template<typename T,
    typename std::enable_if<std::is_same<T, float>::value || std::is_same<T, double>::value || std::is_same<T, uint>::value>::type>
    inline void save_to_bin(const std::string& filename, const std::vector<T>& vec) {
        std::string full_path = LOG_DIR + filename;
        std::ofstream file(full_path, std::ios::binary);
        file.write(reinterpret_cast<const char*>(vec.data()), vec.size() * sizeof(T));
        file.close();
    }

    // -------------------------------
    // Scalar float
    template<typename T,
    typename std::enable_if<std::is_same<T, float>::value || std::is_same<T, double>::value>::type>
    inline void save_to_bin(const std::string& filename, T value) {
        std::string full_path = LOG_DIR + filename;
        std::ofstream file(full_path, std::ios::binary);
        file.write(reinterpret_cast<const char*>(&value), sizeof(T));
        file.close();
    }

    #else

    // If ROSSY_PRINT is 0, make them no-ops
    template<typename... Args>
    void save_to_bin(const Args&...) {}

    #endif // ROSSY_PRINT

} /* rossy_utils */



