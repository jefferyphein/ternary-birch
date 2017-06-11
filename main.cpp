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
    while ((cc = getopt(argc, argv, "e:g:ac:o:p:u:")) != -1)
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
            case 'g':
                genfilename = optarg;
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
        }
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
        genus = std::make_shared<GenusZZ>(*q);
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
            genus->compute_hecke_operators(p, 8);
        }
    }
    else
    {
        // If eigenvector file specified, let's compute Hecke eigenvalues.
        genus->compute_eigenvalues(primes, 8);
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
