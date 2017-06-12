#include <iostream>
#include <cstdlib>
#include <gmpxx.h>
#include <unistd.h>
#include "QuadForm.h"
#include "Genus.h"
#include "Character.h"
#include "Math.h"

typedef Isometry<mpz_class, mpq_class> IsometryZZ;
typedef Genus<mpz_class, mpq_class> GenusZZ;
typedef QuadForm<mpz_class, mpq_class> QuadFormZZ;
typedef Character<mpz_class, mpq_class> CharacterZZ;
typedef Math<mpz_class, mpq_class> MathZZ;

int main(int argc, char** argv)
{
    /*** INITIALIZATION *****************************************************/

    mpz_class a, b, c, f, g, h;
    mpz_class upto = 0;
    std::string eigfilename;
    std::string genfilename;
    std::string outfilename;
    std::vector<mpz_class> conductors;
    std::vector<mpz_class> primes;
    bool allConductors = false;
    bool quadformProvided = false;
    bool helpMessage = false;
    bool helpMessageLong = false;
    int64_t numThreads = 0;
    int64_t maxThreads = std::thread::hardware_concurrency();

    /*** READ QUADRATIC FORM ARGUMENT FROM COMMAND-LINE, IF PROVIDED ********/

    if (argc >= 7)
    {
        bool isQuadForm = true;
        for (int64_t n = argc-6; n < argc; n++)
        {
            std::string s(argv[n]);
            if (s.empty())
            {
                isQuadForm = false;
                continue;
            }
            int64_t pos = s[0] == '-' ? 1 : 0;
            isQuadForm = isQuadForm && (s.substr(pos).find_first_not_of("0123456789") == std::string::npos);
        }

        if (isQuadForm)
        {
            // Assign quadratic form coefficients from end of parameter list.
            a = std::atoi(argv[argc-6]);
            b = std::atoi(argv[argc-5]);
            c = std::atoi(argv[argc-4]);
            f = std::atoi(argv[argc-3]);
            g = std::atoi(argv[argc-2]);
            h = std::atoi(argv[argc-1]);
            argc -= 6;
            quadformProvided = true;
        }
    }

    /*** PARSE COMMAND-LINE ARGUMENTS ***************************************/

    int cc;
    while ((cc = getopt(argc, argv, "he:i:ac:o:p:u:j:x")) != -1)
    {
        switch (cc)
        {
            case 'a':
                allConductors = true;
                break;
            case 'c':
                conductors.push_back(std::move(mpz_class(optarg)));
                break;
            case 'e':
                eigfilename = optarg;
                break;
            case 'h':
                if (helpMessage)
                {
                    helpMessageLong = true;
                }
                helpMessage = true;
                break;
            case 'i':
                genfilename = optarg;
                break;
            case 'j':
                numThreads = std::atoi(optarg);
                if (numThreads < 0)
                {
                    numThreads = 0;
                }
                else if (numThreads > maxThreads)
                {
                    numThreads = maxThreads;
                }
                break;
            case 'o':
                outfilename = optarg;
                break;
            case 'p':
                primes.push_back(mpz_class(optarg));
                break;
            case 'u':
                upto = mpz_class(optarg);
                break;
            case 'x':
                numThreads = maxThreads;
                break;
        }
    }

    /*** DISPLAY HELP MESSAGE IF REQUESTED **********************************/

    if (helpMessage || helpMessageLong)
    {
        std::cout << "Usage: birch [-h[h]] [-x] [-a] [[-c conductor] ...] "
                     "[[-p prime] ...]" << std::endl
                  << "             [-i infile] [-o outfile] [-e eigfile] "
                     "[-j threads] [-u upto]" << std::endl
                  << "             [a b c f g h]" << std::endl;

        std::cout << "Computes the the genus, Hecke operators, and Hecke eigenvalues of integral" << std::endl;
        std::cout << "ternary quadratic forms with shape" << std::endl;
        std::cout << "        ax^2 + by^2 + cz^2 + fyz + gxz + hxy" << std::endl << std::endl;
 
        std::cout << "Options:" << std::endl;
        std::cout << "  -a       attach a character for each squarefree divisor of the discriminant" << std::endl;
        std::cout << "  -c N     attach a character with conductor N" << std::endl;
        std::cout << "  -e FILE  read eigenvector data from FILE" << std::endl;
        std::cout << "  -h       display this message" << std::endl;
        std::cout << "  -hh      displays a longer help message" << std::endl;
        std::cout << "  -i FILE  read genus data from FILE" << std::endl;
        std::cout << "  -j N     launch N threads when performing computations" << std::endl;
        std::cout << "  -o FILE  save genus and Hecke operator data to FILE" << std::endl;
        std::cout << "  -p N     compute Hecke operator or eigenvalues at N; must be a prime" << std::endl;
        std::cout << "  -u N     compute Hecke operators or eigenvalues for all primes <= N" << std::endl;
        std::cout << "  -x       launch the maximum number of threads supported by the system" << std::endl;
        std::cout << std::endl;

        if (helpMessageLong)
        {
            std::cout << "Examples:" << std::endl;
            std::cout << "(1) compute the genus and all Hecke operators at p=2,3" << std::endl;
            std::cout << "       birch -a -o 85085.txt -p 2 -p 3 1 130 172 -65 -1 0" << std::endl << std::endl;
            std::cout << "(2) compute Hecke operator at p=29 after loading data from saved genus file" << std::endl;
            std::cout << "       birch -a -i 85085.txt -o 85085.txt -p 29" << std::endl << std::endl;
            std::cout << "(3) compute the first 1000 Hecke eigenvalues for all eigenvectors loaded from a file" << std::endl;
            std::cout << "       birch -i 85085.txt -e 85085-eigs.txt -u 1000" << std::endl << std::endl;
            std::cout << "(4) compute the first 1000 Hecke eigenvalues using maximum threads" << std::endl;
            std::cout << "       birch -x -i 85085.txt -e 85085-eigs.txt -u 1000" << std::endl;
            std::cout << std::endl;

            std::cout << "Notes:" << std::endl;
            std::cout << "- Eigenvector files should have the following format:" << std::endl << std::endl;
            std::cout << "     M [conductor 1] [conductor 2] ... [conductor M]" << std::endl << std::endl;
            std::cout << "     [vector conductor] [vector coefficients relative to genus reps in genus file]" << std::endl;
            std::cout << "     [num coefficients] [prime 1] [value 1] [prime 2] [value 2] ... [prime N] [value N]" << std::endl;
            std::cout << "     ...repeated for each eigenvector..." << std::endl << std::endl;

            std::cout << "  For example," << std::endl << std::endl;
            std::cout << "     2 1 30" << std::endl << std::endl;
            std::cout << "     1 3 -1" << std::endl;
            std::cout << "     6 2 1 3 1 5 1 7 1 11 -4 13 -2" << std::endl << std::endl; 
            std::cout << "     30 1" << std::endl;
            std::cout << "     6 2 -1 3 -1 5 -1 7 1 11 -4 13 -2" << std::endl;
            std::cout << std::endl;

            std::cout << "  This builds two characters with conductors 1 and 30, then reads two" << std::endl;
            std::cout << "  eigenvectors associated to these conductors, providing their coordinates, and" << std::endl;
            std::cout << "   a few of their Hecke eigenvalues." << std::endl << std::endl;

            std::cout << " - An eigenvector file must be accompanied by a genus files." << std::endl << std::endl;
        }

        return EXIT_SUCCESS;
    }

    /*** QUALITY CONTROL ****************************************************/

    // Handle case where eigenvectors are provided, but no genus provided.
    if (!eigfilename.empty() && genfilename.empty())
    {
        throw std::runtime_error("Must provide input genus file when computing eigenvalues.");
    }

    // Handle case where quadratic form and genus file both provided.
    if (quadformProvided && !genfilename.empty())
    {
        throw std::runtime_error("Ambiguous input: quadratic form and genus file provided.");
    }

    // Handle case where no quadratic form input has been provided.
    if (!quadformProvided && genfilename.empty())
    {
        throw std::runtime_error("No quadratic form input provided.");
    }

    /*** BUILD GENUS ********************************************************/

    // Shared pointers to the genus and initial quadratic form.
    std::shared_ptr<GenusZZ> genus;
    std::shared_ptr<QuadFormZZ> q;

    // Construct the genus.
    if (quadformProvided)
    {
        q = std::make_shared<QuadFormZZ>(a, b, c, f, g, h);
        genus = std::make_shared<GenusZZ>(*q, numThreads);
    }
    else
    {
        genus = std::make_shared<GenusZZ>(genfilename);
        q = genus->quad_form();
    }

    /*** BUILD CHARACTERS ***************************************************/

    // Determine characters to be computed.
    if (eigfilename.empty())
    {
        // If no eigenvectors specified, and the all conductors flag is set,
        // determine squarefree conductors and add their characters.
        if (allConductors)
        {
            auto divs = std::move(MathZZ::squarefree_divisors(q->discriminant()));
            for (auto& d : divs)
            {
                CharacterZZ chi(d);
                genus->add_character(chi);
            }
        }
        else
        {
            // Otherwise, add conductors based upon user input.
            for (mpz_class& cond : conductors)
            {
                std::vector<mpz_class> ps = std::move(MathZZ::prime_divisors_naive(cond));
                CharacterZZ chi(ps);
                genus->add_character(chi);
            }

            // If no user specified conductors, add the trivial character.
            if (conductors.empty())
            {
                CharacterZZ chi(1);
                genus->add_character(chi);
            }
        }
    }
    else
    {
        // If eigenvector file specified, import them. This will implicitly
        // build characters based on the conductors of the eigenvectors being
        // loaded. User input conductors will be ignored to avoid unnecessary
        // computations.
        genus->import_eigenvectors(eigfilename);
    }

    /*** ASSIGN PRIMES IF WE'RE COMPUTING UP TO A LIMIT *********************/

    if (upto > 0)
    {
        primes = std::move(MathZZ::primes_up_to(upto, genus->discriminant()));
    }

    /*** EITHER COMPUTE HECKE OPERATORS OR HECKE EIGENVALUES ****************/

    if (eigfilename.empty())
    {
        // If no eigenvector file specified, we're going to compute Hecke
        // operators instead.
        for (auto& p : primes)
        {
            genus->compute_hecke_operators(p, numThreads);
        }
    }
    else
    {
        // If eigenvector file specified, let's compute Hecke eigenvalues.
        genus->compute_eigenvalues(primes, numThreads);
    }

    /*** SAVE GENUS DATA ****************************************************/

    if (outfilename.empty())
    {
        // If no output filename was specified, and the genus was actually
        // computed, print the genus to stdout.
        if (genus->computed())
        {
            std::cout << *genus;
        }
    }
    else
    {
        // Otherwise, print the genus to the specified file.
        std::fstream f(outfilename, std::ios::out);
        if (!f.is_open())
        {
            std::cout << *genus;
            throw std::runtime_error("Unable to open output file.");
        }

        f << *genus;
        f.close();
    }

    /*** SAVE HECKE EIGENVALUE DATA *****************************************/

    if (!eigfilename.empty())
    {
        genus->export_eigenvectors(eigfilename);
    }

    return EXIT_SUCCESS;
}
