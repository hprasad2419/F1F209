#include <cstdlib> // std::stod, std::stoi
#include <cstring> // strcmp
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "F1F209Wrapper.hh"

using namespace std;

class InputStructure {
public:
  int _p;
  int _n;
  double _Z;
  double _A;
  double _E;
  double _theta;
  int _nop;

  InputStructure()
      : _p(6), _n(6), _Z(6), _A(12), _E(0.680), _theta(60), _nop(100) {}

  void SetBeamEnergy(double E) { _E = E; }
  void SetLeptonAngle(double theta) { _theta = theta; }

  //----------------//

  std::string stripComments(const std::string &line) {
    size_t pos_hash = line.find('#');
    size_t pos_slash = line.find("//");

    size_t pos = std::min(pos_hash, pos_slash);
    if (pos == std::string::npos)
      pos = std::max(pos_hash, pos_slash);

    return line.substr(0, pos);
  }

  std::string trim(const std::string &s) {
    const char *ws = " \t\n\r";
    size_t start = s.find_first_not_of(ws);
    size_t end = s.find_last_not_of(ws);
    return (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
  }

  void loadFromFile(const std::string &filename) {
    std::ifstream file(filename);
    if (!file)
      throw std::runtime_error("Cannot open input file: " + filename);

    std::string line;
    while (std::getline(file, line)) {
      line = stripComments(line);
      line = trim(line);
      if (line.empty())
        continue;

      auto eq_pos = line.find('=');
      if (eq_pos == std::string::npos)
        throw std::runtime_error("Invalid line (missing '='): " + line);

      std::string key = trim(line.substr(0, eq_pos));
      std::string value = trim(line.substr(eq_pos + 1));

      std::istringstream iss(value);

      if (key == "number_of_points")
        iss >> _nop;
      else if (key == "number_of_protons")
        iss >> _p;
      else if (key == "number_of_neutrons")
        iss >> _n;
      else if (key == "beam_energy")
        iss >> _E;
      else if (key == "outgoing_lepton_angle")
        iss >> _theta;
      else
        throw std::runtime_error("Unknown parameter: " + key);

      _Z = _p;
      _A = _p + _n;
    }
  }

  //---------------- //
};

class Application {
public:
  void parseCommandLine(int argc, char **argv) {
    if (argc < 2) {
      PrintUsage("./Fit");
      std::exit(EXIT_FAILURE);
    }

    for (int i = 1; i < argc; ++i) {
      std::string arg = argv[i];

      if (arg == "--input" && i + 1 < argc) {
        inputFile = argv[++i];
        _input.loadFromFile(inputFile);
      } else if (arg == "--output" && i + 1 < argc) {
        outputFile = argv[++i];
      } else if (arg == "--p" && i + 1 < argc) {
        string line = argv[++i];

        auto eq_pos = line.find('=');
        if (eq_pos == std::string::npos)
          throw std::runtime_error("Invalid line (missing '='): " + line);

        std::string key = trim(line.substr(0, eq_pos));
        std::string value = trim(line.substr(eq_pos + 1));

        std::istringstream iss(value);

        if (key == "beam_energy")
          _input.SetBeamEnergy(std::stod(iss.str()));
        else if (key == "outgoing_lepton_angle")
          _input.SetLeptonAngle(std::stod(iss.str()));

      } else {
        PrintUsage("./Fit");
        throw std::runtime_error("Unknown command-line option: " + arg);
      }
    }

    if (inputFile.empty())
      throw std::runtime_error("Missing --input argument");

    if (outputFile.empty())
      throw std::runtime_error("Missing --output argument");
  }
  inline void Run();

private:
  InputStructure _input;
  string inputFile;
  string outputFile;

  void PrintUsage(const char *progName) {
    cerr << "Usage:\n"
         << progName
         << " --input <input-file.txt> --output <output-file.txt> [--p "
            "beam_energy=<value> --p outgoing_lepton_angle=<value>]\n";
  }
  std::string stripComments(const std::string &line) {
    size_t pos_hash = line.find('#');
    size_t pos_slash = line.find("//");

    size_t pos = std::min(pos_hash, pos_slash);
    if (pos == std::string::npos)
      pos = std::max(pos_hash, pos_slash);

    return line.substr(0, pos);
  }

  std::string trim(const std::string &s) {
    const char *ws = " \t\n\r";
    size_t start = s.find_first_not_of(ws);
    size_t end = s.find_last_not_of(ws);
    return (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
  }
};

void Application::Run() {
  static const double kDEG = 3.14159265358979323846 / 180.0;
  static const double M = 0.9383;              // Mass of proton (in GeV)
  static const double Ml = 0.510998910 * 1e-3; // Mass of electron (in GeV)

  cout << "Output file: " << outputFile << endl;
  cout << "Input file: " << inputFile << endl;
  cout << "Number of points: " << _input._nop << endl;
  cout << "Number of protons: " << _input._p << endl;
  cout << "Number of neutrons: " << _input._n << endl;
  cout << "Beam energy: " << _input._E << endl;
  cout << "Outgoing lepton angle: " << _input._theta << endl;
  cout << endl;

  double E = _input._E;
  double theta = _input._theta;
  int Npoints = _input._nop;
  double Z = _input._Z;
  double A = _input._A;

  F1F209Wrapper *model = new F1F209Wrapper();

  std::vector<std::vector<double>> Data;

  for (int index = 0; index < Npoints; index++) {
    double width = (E - Ml) / Npoints;
    double lowerEdge = width * index;
    double upperEdge = lowerEdge + width;
    double binCenter = 0.5 * (lowerEdge + upperEdge);
    double Eout = binCenter;
    double w = E - Eout;

    double xsec_elastic = model->GetXS_QE(Z, A, E, Eout, theta);
    double xsec_mec = model->GetXS_MEC(Z, A, E, Eout, theta);
    double xsec_inclusive = model->GetXS(Z, A, E, Eout, theta);
    std::vector<double> row{w, xsec_elastic, xsec_mec, xsec_inclusive};
    Data.push_back(row);
  }

  std::ofstream out;
  out.open(outputFile.c_str());
  std::string header = "w(GeV)\t\t simga_{QE} (mb/GeV sr)\t\t sigma_{MEC} "
                       "(mb/GeV sr)\t\t sigma_{inclusive} (mb/GeV sr)\n";
  out << header;
  out.flush();
  for (int index = 0; index < Npoints; index++) {
    for (int col = 0; col < Data[index].size(); col++) {
      if (col == Data[index].size() - 1)
        out << Data[index][col] << "\n";
      else
        out << Data[index][col] << " ";
    }
  }
  out.flush();
  std::cout << "Data written in " << outputFile << " ...\n";
}

extern Application app;
