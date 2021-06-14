#include <iostream>
#include <string>

// Function to import the background data from file_name
//  The function will import only the necessary values (from z_in)
bool mg_import(const double a_in, const double fourpiG, mg_cosmology * mg_cosmo, const string file_name) {

    ifstream fin;

    // Check the presence of the file, otherwise throw an error
    try {

      fin.open(file_name, ios::in);
      if (!fin)
        throw std::runtime_error(" The file can not be found or can not be opened.");

    }
    catch (std::exception &ex) {
      std::cout << " Quintessence ERROR :" << ex.what() << endl;
      parallel.abortForce();
    }

    // Counter set to -2 to eliminate header and fist row (bad results)
    int i = -2;

    vector<string> row;
    string line, value;

    do {

        row.clear();

        // Read an entire row and store it in a string variable 'line'
        getline(fin, line);

        // Split the row in words (separated by space)
        stringstream s(line);

        // Populate row vector with different values (eliminating ,)
        // getline(s, value, ',');

        while ( getline(s, value, ',') ) {

            row.push_back(value);
        }

        if(i >= 0){

          if ( stod(row[0]) >= 9. * a_in / 10. ){

            // Populating vectors with file values
            mg_cosmo->a_vec.push_back( stod(row[1]) );
            mg_cosmo->H_vec.push_back( stod(row[2]) * sqrt(2.*fourpiG/3.) );
            mg_cosmo->H_prime_vec.push_back( stod(row[3]) * sqrt(2.*fourpiG/3.) * sqrt(2.*fourpiG/3.) );
            mg_cosmo->w_mg_vec.push_back( stod(row[8]) );
            mg_cosmo->Omega_m_vec.push_back( stod(row[6]) + stod(row[7]) );
            mg_cosmo->Omega_rad_vec.push_back( stod(row[5]) );
            mg_cosmo->Omega_mg_vec.push_back( stod(row[4]) );
            mg_cosmo->particleHorizon_vec.push_back( stod(row[9]) * 2./3./ sqrt(fourpiG) );
            mg_cosmo->mg_field_vec.push_back( stod(row[10]) );
            mg_cosmo->mg_field_p_vec.push_back( stod(row[11]) * sqrt(2.*fourpiG/3.) );

          }

        }

        if (i < 0) i++;

    } while (line != "");

    // Set length of vectors for GSL use
    mg_cosmo->last_int = mg_cosmo->a_vec.size()-1;

    // Defining the pointers to standard arrays
    mg_cosmo->a = &mg_cosmo->a_vec[0];
    mg_cosmo->H = &mg_cosmo->H_vec[0];
    mg_cosmo->H_prime = &mg_cosmo->H_prime_vec[0];
    mg_cosmo->w_mg = &mg_cosmo->w_mg_vec[0];
    mg_cosmo->Omega_m = &mg_cosmo->Omega_m_vec[0];
    mg_cosmo->Omega_rad = &mg_cosmo->Omega_rad_vec[0];
    mg_cosmo->Omega_mg = &mg_cosmo->Omega_mg_vec[0];
    mg_cosmo->particleHorizon = &mg_cosmo->particleHorizon_vec[0];
    mg_cosmo->mg_field = &mg_cosmo->mg_field_vec[0];
    mg_cosmo->mg_field_p = &mg_cosmo->mg_field_p_vec[0];

    // Check if the import is good
    if ( mg_cosmo->a_vec.size() > 1 && mg_cosmo->a_vec.size() == mg_cosmo->H_vec.size() && mg_cosmo->a_vec.size() == mg_cosmo->H_prime_vec.size() && mg_cosmo->a_vec.size() == mg_cosmo->particleHorizon_vec.size() && mg_cosmo->a_vec.size() == mg_cosmo->mg_field_vec.size() && mg_cosmo->a_vec.size() == mg_cosmo->mg_field_p_vec.size()) {
      return true;
    }
    else {
      return false;
    }

}
