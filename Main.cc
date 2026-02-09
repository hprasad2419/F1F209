#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>

#include "Application.hh"
#include "F1F209Wrapper.hh"

Application app;

int main(int argc, char *argv[]) {

  cout << R"(
 /$$$$$$$$                                             /$$                     /$$       /$$$$$$$$ /$$   /$$    
| $$_____/                                            |__/                    | $$      | $$_____/|__/  | $$    
| $$       /$$$$$$/$$$$   /$$$$$$   /$$$$$$   /$$$$$$  /$$  /$$$$$$$  /$$$$$$ | $$      | $$       /$$ /$$$$$$  
| $$$$$   | $$_  $$_  $$ /$$__  $$ /$$__  $$ /$$__  $$| $$ /$$_____/ |____  $$| $$      | $$$$$   | $$|_  $$_/  
| $$__/   | $$ \ $$ \ $$| $$  \ $$| $$$$$$$$| $$  \__/| $$| $$        /$$$$$$$| $$      | $$__/   | $$  | $$    
| $$      | $$ | $$ | $$| $$  | $$| $$_____/| $$      | $$| $$       /$$__  $$| $$      | $$      | $$  | $$ /$$
| $$$$$$$$| $$ | $$ | $$| $$$$$$$/|  $$$$$$$| $$      | $$|  $$$$$$$|  $$$$$$$| $$      | $$      | $$  |  $$$$/
|________/|__/ |__/ |__/| $$____/  \_______/|__/      |__/ \_______/ \_______/|__/      |__/      |__/   \___/  
                        | $$                                                                                    
                        | $$                                                                                    
                        |__/                                                                                    
  )" << endl;

  try {
    app.parseCommandLine(argc, argv);
    app.Run();
  } catch (const exception &ex) {
    cout << "Erro : " << ex.what() << endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;

  return 0;
}
