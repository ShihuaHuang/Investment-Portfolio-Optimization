// Written by Prof. Sreenivas for IE523: Financial Computing

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cstdlib>

#include "lp_lib.h"

using namespace std;

const double error = 1e-10;//create a error to use in Newton-Raphson method as a measurement of accuracy
int number_of_cash_flows;
vector <double> price_list;
vector <int> maturity_list;
vector <double> yield_to_maturity;
vector <double> duration;
vector <double> convexity;
double debt_obligation_amount;
double time_when_debt_is_due;
vector <double> percentage_of_cash_flow_to_meet_debt_obligation;

double function(vector <double> cash_flow, double price, int maturity, double rate)
{
	// write a function that computes f(r) in page 2 of lesson 3 of my notes
	double sum, fr;
	sum = 0;
	for (int i = 1; i <= maturity; i++)
		sum = sum + cash_flow[i - 1]*pow((1 + rate), maturity - i);//calculate the term of summation
	return fr = price*pow((1 + rate), maturity) - sum;//use root finding function to find YTM
}

double derivative_function(vector <double> cash_flow, double price, int maturity, double rate)
{
	// write a function that computes f'(r) in the bottom of page 2 of lesson 3 
	// of my notes
	double sum, dfr;
	sum = 0;
	for (int i = 1; i <= maturity - 1; i++)
		sum = sum + cash_flow[i - 1] * (maturity - i)*pow((1 + rate), maturity - i - 1);
	return dfr = maturity*price*pow((1 + rate), maturity - 1) - sum;//calculate the derivative of root finding function
}

double Newton_Raphson(vector <double> cash_flow, double price, int maturity, double rate)
{
	// write a function that finds the (only) +ve root of f(r) of page 2 of 
	// lesson 3 using Newton-Raphson method
	double ratenew;
	for (int i = 0; i < 100; i++)
	{
		ratenew = rate - function(cash_flow, price, maturity, rate) / derivative_function(cash_flow, price, maturity, rate);
		rate = ratenew;//calculate for 100 times to find the closest value of YTM, using recurision
	}
	return rate;
}

double get_duration(vector <double> cash_flow, double price, int maturity, double rate)
{
	// write a function that computes the duration of a cash flow
	double sum, duration;
	sum = 0;
	for (int i = 1; i <= maturity; i++)
		sum = sum + i*cash_flow[i - 1] / pow((1 + rate), i);//the formula for duration on page3 of lecture 4
	return duration = sum / price;//let the variation of discounted value of obligation equal to the variation of bond's present value
}

double get_convexity(vector <double> cash_flow, double price, int maturity, double rate)
{
	// write a function that computes the convexity of a cash flow
	double sum, convexity;
	sum = 0;
	for (int i = 1; i <= maturity; i++)
		sum = sum + i*(i + 1)*cash_flow[i - 1] / pow((1 + rate), i + 2);//the formula for convexity on page8 of lecture 4
	return convexity = sum / price;//bond-convexity is the second-derivative of the bond-price
}

double present_value_of_debt()
{
	// compute PV of future debt obligation
	// using the average-value-of-the-YTMs 
	double sum = 0;
	for (int i = 0; i < number_of_cash_flows; i++)
	{
		sum = sum + yield_to_maturity[i];
	}
	double mean = sum / number_of_cash_flows;//using the mean of bonds' YTM
	double pv = debt_obligation_amount / pow((1 + mean), time_when_debt_is_due);//the formula for present value equation1 on page2 of lecture 4
	return pv;
}

void print_data(char *filename)
{
	cout << "Input File: " << filename << endl;
	cout << "We owe " << debt_obligation_amount << " in " << time_when_debt_is_due << " years" << endl;
	cout << "Number of Cash Flows: " << number_of_cash_flows << endl;
	double sum = 0;
	for (int i = 0; i < number_of_cash_flows; i++)
	{
		cout << "---------------------------" << endl;
		cout << "Cash Flow #" << i + 1 << endl;
		cout << "Price = " << price_list[i] << endl;
		cout << "Maturity = " << maturity_list[i] << endl;
		cout << "Yield to Maturity = " << yield_to_maturity[i] << endl;
		cout << "Duration = " << duration[i] << endl;
		cout << "Convexity = " << convexity[i] << endl;
		cout << "Percentage of Face Value that would meet the obligation = " <<
			percentage_of_cash_flow_to_meet_debt_obligation[i] << endl;
		sum = sum + yield_to_maturity[i];
	}
	cout << "***************************" << endl;
	cout << "Average YTM (which is used to compute PV of Debt) = " << sum / number_of_cash_flows << endl;//print average YTM as we calculated
	cout << "Present value of debt = " << present_value_of_debt() << endl;//print present value of debt as we calculated
	cout << "***************************" << endl;
}

void get_data(char* argv[])
{
	// write the code that reads the data from the file identified 
	// on the command-line. 
	double value;
	vector <double> p;
	vector <double> cash_flow;
	ifstream input_filename(argv[1]);//read file as we input the name in the command line
	if (input_filename.is_open())
	{
		while (input_filename >> value)
		{
			p.push_back(value);//read every value one by one and store them in p vector
		}
	}
	number_of_cash_flows = p[0];//the value of number of cash flow is the first element in p vector
	debt_obligation_amount = p[p.size() - 2];//the value of debt obligation amount is the last but one element in p vector
	time_when_debt_is_due = p[p.size() - 1];//the value of time when debt is due is the last element in p vector
	for (int i = 0; i < number_of_cash_flows; i++)
	{
		double sum = 0;
		for (int j = 0; j < i; j++)
		{
			sum = sum + maturity_list[j];//count total number of cash flows to calculate price and maturity element
		}
		price_list.push_back(p[1 + sum + 2 * i]);//create price vector and store no.1+sum+2*i value in p vector, we add number of cash flows between them
		maturity_list.push_back(p[2 + sum + 2 * i]);//create maturity vector and store no.2+sum+2*i value in p vector, we add number of cash flows between them

		int k = maturity_list[i];
		cash_flow.assign(p.begin() + 3 + sum + 2 * i, p.begin() + 3 + sum + 2 * i + k);//given bond's number, assign cashflow from their right place, k elements in total

		double startingYTM = 0.06;//give a starting number of YTM to calculate
		yield_to_maturity.push_back(Newton_Raphson(cash_flow, price_list[i], maturity_list[i], startingYTM));//use Newton_Raphson method to calculate YTM and 
		                                                                                                     //store them in a vector
		duration.push_back(get_duration(cash_flow, price_list[i], maturity_list[i], yield_to_maturity[i]));//use get_duration function to calculate duration 
		                                                                                                   //and store them in a vector
		convexity.push_back(get_convexity(cash_flow, price_list[i], maturity_list[i], yield_to_maturity[i]));//use get_convexity function to calculate convexity 
		                                                                                                     //and store them in a vector
	}
	for (int i = 0; i < number_of_cash_flows; i++)
	{
		percentage_of_cash_flow_to_meet_debt_obligation.push_back(present_value_of_debt() / price_list[i]);//use present_value_of_debt function to calculate the 
		                                                                                                   //percentage of cash flow and store them in a vector
	}
}

void get_optimal_portfolio()
{
	lprec *lp;
	vector<double>money;
	// write the lp_solve specific material that 
	// computes the optimal_portfolio;
	REAL solution[10];

	lp = make_lp(0, number_of_cash_flows);//set number of variables equal to the number of cash flows
	{
		double row[10];
		for (int i = 1; i <= number_of_cash_flows; i++)
			row[i] = price_list[i - 1];//the cofficient equal to the price of bond
		row[0] = 0;
		add_constraint(lp, row, EQ, present_value_of_debt()); //add constraint, set it equal to present value of debt
	}

	{
		double row[10];
		for (int i = 1; i <= number_of_cash_flows; i++)
			row[i] = duration[i - 1]/ percentage_of_cash_flow_to_meet_debt_obligation[i - 1];//the cofficient equal to the duration divided by percentage of cash flow
		row[0] = 0;
		add_constraint(lp, row, EQ, time_when_debt_is_due);//add constraint, set it equal to time when debt is due
	}

	{
		double row[10];
		for (int i = 1; i <= number_of_cash_flows; i++)
			row[i] = -convexity[i - 1]/ percentage_of_cash_flow_to_meet_debt_obligation[i - 1];//the cofficient equal to the convexity divided by percentage of cash flow
		row[0] = 0;
		set_obj_fn(lp, row);//set objective of minimizing the negative convexity
	}

	double ret = solve(lp);

	// print the optimizing values of the variables
	if (ret == 0)//if there are solutions
	{
		get_variables(lp, solution);//solve lp problem
		cout << "Largest Convexity we can get is:" << -1 * get_objective(lp) << endl;//print the maximium of convexity
		cout << "Optimal portfolio:" << endl;
		for (int i = 0; i < number_of_cash_flows; i++)
		{
			cout << "%Cash Flow" << i + 1 << ": " << solution[i] << endl;//print the solution of percentage of bond
		}
		cout << "That is, buy" << endl;
		for (int i = 0; i < number_of_cash_flows; i++)
		{
			money.push_back(solution[i] * price_list[i]);//print the calculation of value of bond, if solution is zero then print the $0 bond
			cout << "$" << money[i] << " of Cash Flow#" << i + 1 << endl;
		}
	}
	else
	{
		cout << "There is no portfolio that meets the duration constraint of "<< time_when_debt_is_due <<" years" << endl;//if there is no solution, print out accordingly
	}
	delete_lp(lp);
}

int main(int argc, char* argv[])
{
	if (argc == 1) {
		cout << "Input filename missing" << endl;
	}
	else
	{
		get_data(argv);

		print_data(argv[1]);

		get_optimal_portfolio();
	}
	return (0);
}