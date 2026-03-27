/*
 *
 *
 * Macro for performing B feed-down subtraction 
 * and applying correction to BackgroundSubtraction.C
 * and SignalExtraction.C resulting distributions.
 * 
 * 
 * 
 * 
**/


using namespace std;



void FeedDownSubtraction(){

    // Execution time calculation
    time_t start, end;
    time(&start); // initial instant of program execution

    

    time(&end); // end instant of program execution
    
    // Calculating total time taken by the program. 
    double time_taken = double(end - start); 
    cout << "Time taken by program is : " << fixed 
         << time_taken/60 << setprecision(5); 
    cout << " min " << endl; 

    
}

int main(){
    FeedDownSubtraction();
    return 0;
}
