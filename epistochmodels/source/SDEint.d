/**
Ordinary Stochastic differencial Equations solvers
*/
module sde;
import std.math;
import std.stdio;
// import std.parallelism;
import std.range : enumerate;
import std.typecons : tuple, Tuple;
import core.exception : RangeError;
// import mir.ndslice;
import mir.random;
// import mir.math;
import mir.random.variable : uniformVar, exponentialVar, normalVar;
// import pyd.pyd;

/// Defining deterministic (f) and Stochastic components of a Stochastic diferential equation.
alias Shift = double[] delegate(double, double[], double[]);
alias Drift = double[] delegate(double, double[], double[]);

/**
Euler-Maruyama first-order integrator
*/
class EMaruyama
{
    Shift f;
    Drift g;

    this(Shift f, Drift g)
    {
        this.f = f;
        this.g = g;
    }
    
    double[][] integrate(double[] inits, double[] trange, double[] pars) 
    {
        double[][] sol; /// Solution 
        sol ~= inits;
        foreach (i,t; trange){
            double h = i==0 ? trange[i+1]-t : t-trange[i-1];
            auto nrv = normalVar!double(0,h);
            // writeln(i, "==>", this.g(t, sol[i],pars));
            sol ~= sol[i].dup;
            try{
                sol[$-1][] = f(t,sol[i],pars)[]*h + g(t, sol[i], pars)[] * nrv(rne);
            }
            catch(Exception){
                writeln("error: ", sol[$-1]);
                break;
            }
            writeln("State: ", sol[$-1]);
        }
        return sol;
    }
}

@("Euler-maruyama growth")    
unittest{
    double[] F(double t, double[] y, double[] pars) @safe
    {
        double b = pars[0];
        double d = pars[1];
        auto res = y.dup;
        res[] *= (b-d);
        // writeln("f ==> ",y,res);
        return res;
    }
    double[] G(double t, double[] y, double[] pars){
        double b = pars[0];
        double d = pars[1];
        auto res = y.dup;
        res[] *= (b+d);
        foreach (i,x; res){
            res[i] = sqrt(x);
        }
        // writeln("g ==> ", y,res);
        return res;
    }
    double[1] inits = [20];
    double[100] trange;
    for (uint i=0; i<trange.length; i++){
        trange[i] = i;
    }
    double[2] pars = [0.3,0.27];
    auto EM = new EMaruyama(&F, &G);
    auto res = EM.integrate(inits,trange,pars);
}


    /**
    Milstein second-order integrator
    */
    class Milstein
    {
        Shift f;
        Drift g;

        this(Shift f, Drift g)
        {
            this.f = f;
            this.g = g;
        }

        /// Integrates the system.
        double[][] integrate(double[] inits, double[] trange, double[] pars)
        {
        double [][] sol; /// Solution array
        sol[0] = inits;
        return sol;
        }
    }
