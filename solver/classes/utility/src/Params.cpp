#include "../include/Params.h"

#include <fstream>
#include <iostream>

using namespace std;

rapidjson::Document Params::getAllParams() {
    std::ifstream ifs("config/params.json");
    string jsonContent((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));
    rapidjson::Document d;
    d.Parse(jsonContent.c_str());
    return d;
}
