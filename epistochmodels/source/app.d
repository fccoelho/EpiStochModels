import std.stdio;
import std.datetime.stopwatch;
import std.file;
import core.exception : RangeError;

import models: SIR, SIR_Dem;


void main()
{
    double beta = 0.7;
    double gam = 0.3;
    int N = 1000000;
    bool constant = false;
    int I0 = 10;
    double tf = 1000;
    auto sw = StopWatch(AutoStart.no);
    auto model = new SIR(N, beta, gam);
    model.initialize(N-I0, I0, 0);
    sw.start();
    auto sim = model.run(0, tf);
    sw.stop();
    writeln("Time of the SIR run with N=100000: ", sw.peek());
    writefln("Number of steps: %s", sim[0].length);

    File outf = File("sim.csv", "w");
    outf.writeln("time,S,I,dt");
    foreach (i, double t; sim[0])
    {
        try
        {
            outf.writefln("%s,%s,%s,%s", t, sim[1][i][0], sim[1][i][1], sim[1][i][2]);
        }
        catch (RangeError e)
        {
            writefln("t: %s\t S: %s\t I: %s\t dt:%s", t, sim[1][i][0], sim[1][i][1], sim[1][i][2]);
        }
    }
    outf.close();
    // double alpha = 0.1;
    // auto model2 = new SIR_Dem(N, alpha, beta, gam);
    // model2.initialize(N-I0, I0, 0);
    // sw.start();
    // auto sim2 = model2.run(0, tf);
    // sw.stop(); 
    // writeln("Time of the SIR_Dem run with N=100000: ",sw.peek());
    // writefln("Number of steps for SIR_Dem: %s", sim2[0].length);
}
