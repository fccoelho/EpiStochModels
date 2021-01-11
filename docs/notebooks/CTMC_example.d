/+dub.json:
{
	"name": "CTMC_example",
	"authors": [
		"Flávio Codeço Coelho"
	],
	"description": "Example usage of SIR model from epistochmodels",
	"copyright": "Copyright © 2020, Flávio Codeço Coelho",
	"license": "GPL-3",
    "dependencies": {
		"epistochmodels": "*"
	}
}
+/
import std.datetime.stopwatch;
import std.stdio;
import std.algorithm;
import models;

void main()
{
    int[][] tmat = [[-1, 1, 0, 0], [0, -1, 1, 0], [0, -1, 0, 1]];
    double[string] pars = ["beta" : 0.3, "gam" : 0.05, "mu" : 0.01];
    alias double delegate(int[], double[string]) dy;
    dy[] props = [
        (v, p) => p["beta"] * v[0] * v[1] / sum(v[0 .. $ - 1]), // infection
        (v, p) => p["gam"] * v[1], // recovery
        (v, p) => p["mu"] * v[1] // death
    ];
    CTMC model = new CTMC(tmat, props);
    model.initialize([995, 5, 0, 0], pars);
    auto sw = StopWatch(AutoStart.no);
    sw.start();
    auto res = model.run(0, 1000);
    sw.stop();
    writeln("Time of the SIR run with N=1000: ", sw.peek());
}
