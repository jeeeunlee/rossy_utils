#include "rossy_utils/io/io_utilities.hpp"

#include <string>
#include <vector>

namespace rossy_utils {

void cleaningFile(std::string _file_name, std::string& _ret_file, bool b_param) {
    if (b_param)
        _ret_file += CURRENT_DIR;
    else
        _ret_file += CURRENT_DIR "experiment_data/";

    _ret_file += _file_name;
    _ret_file += ".txt";

    std::list<std::string>::iterator iter = std::find(
        gs_fileName_string.begin(), gs_fileName_string.end(), _file_name);
    if (gs_fileName_string.end() == iter) {
        gs_fileName_string.push_back(_file_name);
        remove(_ret_file.c_str());
    }
}


void pretty_constructor(const int& _num_tab, const std::string& _name) {
    myColor color;
    color_print(myColor::BoldCyan, "|", false);
    std::string content = " ";
    int space_to_go(0);
    if (_num_tab != 0) {
        for (int i = 0; i < _num_tab; ++i) {
            content += "    ";
        }
        content = content + "||--" + _name;
        switch (_num_tab) {
            case 1:
                color = myColor::BoldGreen;
                break;
            case 2:
                color = myColor::BoldYellow;
                break;
            case 3:
                color = myColor::BoldBlue;
                break;
            case 4:
                color = myColor::BoldMagneta;
                break;
            default:
                std::cout << "no such color in pretty_constructor" << std::endl;
                exit(0);
        }
    } else {
        content += _name;
        color = myColor::BoldRed;
    }
    space_to_go = 78 - content.length();
    // std::cout << space_to_go << std::endl;
    for (int i = 0; i < space_to_go; ++i) {
        content += " ";
    }
    color_print(color, content, false);
    color_print(myColor::BoldCyan, "|");
}

void color_print(const myColor& _color, const std::string& _name,
                 bool line_change) {
    switch (_color) {
        case Red:
            printf("\033[0;31m");
            break;
        case BoldRed:
            printf("\033[1;31m");
            break;
        case Green:
            printf("\033[0;32m");
            break;
        case BoldGreen:
            printf("\033[1;32m");
            break;
        case Yellow:
            printf("\033[0;33m");
            break;
        case BoldYellow:
            printf("\033[1;33m");
            break;
        case Blue:
            printf("\033[0;34m");
            break;
        case BoldBlue:
            printf("\033[1;34m");
            break;
        case Magneta:
            printf("\033[0;35m");
            break;
        case BoldMagneta:
            printf("\033[1;35m");
            break;
        case Cyan:
            printf("\033[0;36m");
            break;
        case BoldCyan:
            printf("\033[1;36m");
            break;
        default:
            std::cout << "No Such Color" << std::endl;
            exit(0);
    }
    if (line_change)
        printf("%s\n", _name.c_str());
    else
        printf("%s", _name.c_str());
    printf("\033[0m");
}




void readFile(std::string _file_name, std::vector<std::string>& _vec) {
    std::ifstream InputFile(_file_name.c_str());
    std::string tempstring;
    if (!InputFile.is_open()) {
        std::cout << "Data file load error... check the data file" << std::endl;
        exit(0);
    } else {
        while (!InputFile.eof()) {
            InputFile.clear();
            std::getline(InputFile, tempstring);
            _vec.push_back(tempstring);
        }
        InputFile.close();
    }
}

void splitString(std::string* str_array, std::string strTarget,
                 std::string strTok) {
    int nCutPos = 0;
    int nIndex = 0;
    while ((nCutPos = strTarget.find_first_of(strTok)) != strTarget.npos) {
        if (nCutPos > 0) {
            str_array[nIndex++] = strTarget.substr(0, nCutPos);
        }
        strTarget = strTarget.substr(nCutPos + 1);
    }
    if (strTarget.length() > 0) {
        str_array[nIndex++] = strTarget.substr(0, nCutPos);
    }
}

// explicitly instatiate
// template void pretty_print<double>(const std::vector<double>&, const char*);
// template void pretty_print<float>(const std::vector<float>&, const char*);

// template void pretty_print<float>(const Eigen::MatrixXf&, std::ostream& ,
//                     const std::string&, const std::string& ,
//                     bool, bool);
// template void pretty_print<double>(const Eigen::MatrixXd&, std::ostream& ,
//                     const std::string&, const std::string& ,
//                     bool, bool);

}  // namespace rossy_utils
