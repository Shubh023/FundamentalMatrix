// Imagine++ project
// Project:  Fundamental
// Author:   Pascal Monasse

#include "./Imagine/Features.h"
#include <Imagine/Graphics.h>
#include <Imagine/LinAlg.h>
#include <vector>
#include <cstdlib>
#include <ctime>
using namespace Imagine;
using namespace std;

static const float BETA = 0.01f; // Probability of failure

struct Match {
    float x1, y1, x2, y2;
};

// Display SIFT points and fill vector of point correspondences
void algoSIFT(Image<Color,2> I1, Image<Color,2> I2,
              vector<Match>& matches) {
    // Find interest points
    SIFTDetector D;
    D.setFirstOctave(-1);
    Array<SIFTDetector::Feature> feats1 = D.run(I1);
    drawFeatures(feats1, Coords<2>(0,0));
    cout << "Im1: " << feats1.size() << flush;
    Array<SIFTDetector::Feature> feats2 = D.run(I2);
    drawFeatures(feats2, Coords<2>(I1.width(),0));
    cout << " Im2: " << feats2.size() << flush;

    const double MAX_DISTANCE = 100.0*100.0;
    for(size_t i=0; i < feats1.size(); i++) {
        SIFTDetector::Feature f1=feats1[i];
        for(size_t j=0; j < feats2.size(); j++) {
            double d = squaredDist(f1.desc, feats2[j].desc);
            if(d < MAX_DISTANCE) {
                Match m;
                m.x1 = f1.pos.x();
                m.y1 = f1.pos.y();
                m.x2 = feats2[j].pos.x();
                m.y2 = feats2[j].pos.y();
                matches.push_back(m);
            }
        }
    }
}

/**
 * Function to check if a value exists in a vector<int>
 */
bool contains(vector<int> vec, const int & elem)
{
    bool result = false;
    for (auto & x : vec)
    {
        if (x == elem)
        {
            result = true;
            break;
        }
    }
    return result;
}

/**
 * Get a random Sample of matches
 * @param matches input vector of matches
 * @param nbSample nb of samples to provide
 * @return return nbSample number of random samples from matches vector
 */
vector<Match> getSample(vector<Match>& matches, size_t nbSample)
{
    vector<int> indexes(nbSample, -1);
    vector<Match> sample;

    for (size_t i = 0; i < nbSample; i++)
    {
        auto index = rand() % matches.size();
        if (!contains(indexes, index)) {
            indexes[i] = index;
            sample.push_back(matches[index]);
        }
    }
    return sample;
}

/**
 * Compute Matrix A
 * @param points : vector containing pairs of matching points
 * @return computed matrix A
 */
FMatrix<float, 9, 9> compute_A(vector<Match> samples) {
    // Create linear system matrix
    FMatrix<float,9,9> A;
    for(int i = 0; i < 8; i++){
        Match m = samples[i];
        // Get points
        DoublePoint3 p1, p2;
        p1[0] = m.x1;
        p1[1] = m.y1;
        p1[2] = 1;
        p2[0] = m.x2;
        p2[1] = m.y2;
        p2[2] = 1;
        A(i,0) = p1[0] * p2[0];
        A(i,1) = p1[0] * p2[1];
        A(i,2) = p1[0];
        A(i,3) = p1[1] * p2[0];
        A(i,4) = p1[1] * p2[1];
        A(i,5) = p1[1];
        A(i,6) = p2[0];
        A(i,7) = p2[1];
        A(i,8) = 1;
    }
    for(int i = 0; i < 9; i++)
        A(8,i) = 0;

    A(8,8) = 0;
    return A;
}


/**
 * Get Epipolar distance between two match Points using the calulated F Matrix
 * @param m : a given match of points
 * @param F : some precomputed Fundamental matrix
 * @return (float) Epipolar distance between point p1 & p2 from the given match m
 */
float compute_EPDistance(Match m, FMatrix<float, 3, 3>& F) {
    DoublePoint3 p1, p2;
    // Assign Matching points values to p1 & p2
    p1[0] = m.x1;
    p1[1] = m.y1;
    p1[2] = 1;
    p2[0] = m.x2;
    p2[1] = m.y2;
    p2[2] = 1;

    FVector<float, 3> tf;
    tf = F * p2;
    float norm = sqrt(pow(tf[0],2.0) + pow(tf[1],2.0));
    tf /= norm;
    auto distance = abs(p1 * tf);
    return distance;
}

/**
 * Normalize the values from all points of a Sample of Matches
 * @param sample
 * @return (vector<Match>) vector with all matche's values normalized
 */
vector<Match> normalizeSample(vector<Match> sample) {
    vector<Match> normalized_sample(sample);
    for (size_t i = 0; i < normalized_sample.size(); i++) {
        // Point 1
        normalized_sample[i].x1 *= 10e-3;
        normalized_sample[i].y1 *= 10e-3;

        // Point 2
        normalized_sample[i].x2 *= 10e-3;
        normalized_sample[i].y2 *= 10e-3;
    }
    return normalized_sample;
}


FMatrix<float,3,3> computeFundamentalMatrix(vector<Match>& samples){


    // Normalization of all matching points
    auto normalized_sample = normalizeSample(samples);

    auto A = compute_A(normalized_sample);

    // Use SVD to solve equations
    FVector<float,9> S;
    FMatrix<float,9,9> U, V;
    svd(A,U,S,V);

    FMatrix<float,3,3> tF;
    for (int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            tF(i, j)= V.getRow(8)[3 * i + j];
     
    // Get the desired rank for F using SVD again to solve with rank as a constraint
    FVector<float,3> S2;
    FMatrix<float,3,3> U2, V2;
    svd(tF, U2, S2, V2);
    S2[2] = 0;
    tF = U2 * Diagonal(S2) * V2;

    // Normalize values of tF
    FMatrix<float,3,3> N(0.f);
    N(0,0) = 10e-3;
    N(1,1) = 10e-3;
    N(2,2) = 1;
    tF = N * tF * N;
    return tF;
}


// RANSAC algorithm to compute F from point matches (8-point algorithm)
// Parameter matches is filtered to keep only inliers as output.
FMatrix<float,3,3> computeF(vector<Match>& matches) {
    const float distMax = 1.5f; // Pixel error for inlier/outlier discrimination
    int Niter=100000; // Adjusted dynamically
    FMatrix<float,3,3> bestF;
    vector<int> bestInliers;
    // --------------- TODO ------------
    // DO NOT FORGET NORMALIZATION OF POINTS

    int nbSamples = 8;
    int ct = 0;
    vector<int> inliers;
    vector<Match> samples;
    FMatrix<float,3,3> F;

    // Iteration on Niter for RANSAC algorithm
    while (ct < Niter)
    {
        // Get some random sample of size nbSamples from all elements of matches
        samples = getSample(matches, nbSamples);

        // Using the extract sample of size 8 compute a fundamental Matrix F
        F = computeFundamentalMatrix(samples);

        // Clear old computed inliers
        inliers.clear();

        // Calulate new inliers using the EpiPolar Distance function
        for(int i = 0; i < matches.size(); i++){
            auto epdist = compute_EPDistance(matches[i], F);

            // If distance is under the max dist bound than keep match m as an inlier
            if(epdist < distMax){
                inliers.push_back(i);
            }
        }

        // Compute the bestF Matrix by keeping the current F Matrix and inliers as best possible result we can get
        if (inliers.size() > bestInliers.size()){
            bestF = F;
            bestInliers = inliers;
            Niter = min((int)(std::log(BETA) / std::log(1 - pow(((float)bestInliers.size() / matches.size()), nbSamples))), 10000);
        }

        // Increment the counter
        ct++;
    }

    cout << "Nb of Iterations needed for RANSAC to get bestF: " << ct << endl;

    // Remaining to add is the Optimization with mean squared method...


    // Updating matches with inliers only
    vector<Match> all=matches;
    matches.clear();
    for(size_t i=0; i<bestInliers.size(); i++)
        matches.push_back(all[bestInliers[i]]);
    return bestF;
}

// Expects clicks in one image and show corresponding line in other image.
// Stop at right-click.
void displayEpipolar(Image<Color> I1, Image<Color> I2,
                     const FMatrix<float,3,3>& F) {
    while(true) {
        int x,y;
        if(getMouse(x,y) == 3)
            break;
        // --------------- TODO ------------
        DoublePoint3 pt1, pt2;
        int w = I1.width();
        // If x belongs to I1
        if(x <= w){
            pt1[0] = x; pt1[1] = y; pt1[2] = 1;
            FVector<float, 3> line;
            line = transpose(F) * pt1;
            // Draw corresponding line and circle in the corresponding images
            fillCircle(x,y,4, BLUE);
            drawLine(w,
                     (-1)*(line[2])/line[1],
                     2*w,
                     (-1)*(line[2]+line[0]*w)/line[1],
                     BLUE);
        }
        // If x belongs to I2
        if(x > w){
            pt2[0] = x - w; pt2[1] = y; pt2[2] = 1;
            FVector<float, 3> line;
            line = F * pt2;            
            // Draw corresponding line and circle in the corresponding images
            fillCircle(x,y,4, RED);
            drawLine(0,
                     (-1)*line[2]/line[1],
                     w,
                     (-1)*(line[2]+line[0]*w)/line[1],
                     RED);
        }
    }
}

int main(int argc, char* argv[])
{
    srand((unsigned int)time(0));

    const char* s1 = argc>1? argv[1]: srcPath("im1.jpg");
    const char* s2 = argc>2? argv[2]: srcPath("im2.jpg");

    // Load and display images
    Image<Color,2> I1, I2;
    if( ! load(I1, s1) ||
        ! load(I2, s2) ) {
        cerr<< "Unable to load images" << endl;
        return 1;
    }
    int w = I1.width();
    openWindow(2*w, I1.height());
    display(I1,0,0);
    display(I2,w,0);

    vector<Match> matches;
    algoSIFT(I1, I2, matches);
    const int n = (int)matches.size();
    cout << " matches: " << n << endl;
    drawString(100,20,std::to_string(n)+ " matches", RED);
    click();

    FMatrix<float,3,3> F = computeF(matches);
    cout << "F="<< endl << F;

    // Redisplay with matches
    display(I1,0,0);
    display(I2,w,0);
    for(size_t i=0; i<matches.size(); i++) {
        Color c(rand()%256,rand()%256,rand()%256);
        fillCircle(matches[i].x1+0, matches[i].y1, 2, c);
        fillCircle(matches[i].x2+w, matches[i].y2, 2, c);
    }
    drawString(100, 20, to_string(matches.size())+"/"+to_string(n)+" inliers", RED);
    click();

    // Redisplay without SIFT points
    display(I1,0,0);
    display(I2,w,0);
    displayEpipolar(I1, I2, F);

    endGraphics();
    return 0;
}
