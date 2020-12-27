/+dub.json:
{
	"name": "sir_example",
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
import models;

void main(){
    double beta = 0.7;
    double gam = 0.3;
    int N = 1_000_000;
    int I0 = 10;
    double tf = 1000;
    auto sw = StopWatch(AutoStart.no);
    auto model = new SIR(N, beta, gam);
    model.initialize(N-I0, I0, 0);
    sw.start();
    auto sim = model.run(0, tf);
    sw.stop();
    writeln("Time of the SIR run with N=1000000: ", sw.peek());
    writefln("Number of steps: %s", sim[0].length);
}