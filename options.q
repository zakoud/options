// Options calculations library

/ Constants
.opt.pi:                      acos -1;
.opt.inverseSquareRootOf2Pi:  reciprocal sqrt 2 * .opt.pi;
.opt.inverseSquareRootOfTwo:  reciprocal sqrt 2;
.opt.optionModelApplied:      `Black76;
.opt.pow2:                    {x * x};


/ Cumulative Density function (CDF) of a standard normal (Gaussian) random variable. It is closely related to the error function erf
/ New version, more accurate - Version 19/11/2014
/ "Handbook of Mathematical Functions" (1965), by Abramowitz and Stegun
/ Equation 7.1.26 
.opt.cdf:{

    L:abs x;
   
    K : 1.0 % (1.0 + (0.2316419 * L));
    
    a1 : 0.31938153;
    a2 : -0.356563782;
    a3 : 1.781477937;
    a4 : -1.821255978;
    a5 : 1.330274429;
    a12345k:(a1 * K)+ (a2 * K * K) + (a3 * K * K * K) + (a4 * K * K * K * K) + (a5 * K * K * K * K * K);
     
    result: 1.0 - .opt.inverseSquareRootOf2Pi * exp[neg .opt.pow2[L] % 2.0] * a12345k;
    if[x<0;
        result:1 - result];
    result
 }
 
 
/ Black-Scholes Modeling
/ Extension of ql.q from github.com/kimtang/ql.q 

.opt.bls:()!();

.opt.bls[`ParamsTable]:flip `spot`strike`maturity`volatility`yield`rate`pricingStart`averageSpot`direct`type_!"FFFFFFFFSS"$\:();

.opt.bls[`CalculateOptionPrice]:{[params]
    modelApplied:.opt.optionModelApplied;
    if[not `type_ in key params;'`$"Type is missing"];
    if[`American = params`type_; :.opt.bls[`BjerksundStenslandFormula][params]];
    if[`Asian = params`type_; :.opt.bls[`TurnbullWakemanFormula][params]];

    $[(`European=params[`type_]) & (modelApplied=`BlackScholes);
		[params,:.opt.bls[`BlackScholesDistribution] params; :.opt.bls[`BlackScholesFormula] params];
		[params,:.opt.bls[`Black76Distribution] params;:.opt.bls[`Black76Formula] params]];
 }; 


/ Black's 76 Model
.opt.bls[`Black76Formula]:{[params]
	neg params[`sign]*(exp[neg params[`rate]*params[`maturity]]) * ((params[`strike] * .opt.cdf params[`sign]*params[`d2]) - (params[`spot] * .opt.cdf params[`sign]*params[`d1]))
 };

.opt.bls[`Black76Distribution]:{ [params] 
    d1:(log[params[`spot] % params[`strike]] + 0.5 * (params[`maturity]*params[`volatility] xexp 2)) % params[`volatility] * sqrt params[`maturity];
    d2:d1 - params[`volatility] * sqrt params[`maturity];
    sign:neg 1f-2f*`float$(params[`direct]=`call);
    (`d1`d2`sign)!(d1;d2;sign)
 };
 
 
/ Black-Scholes Model 
.opt.bls[`BlackScholesFormula]:{[params]
	neg params[`sign]*(params[`strike]*exp[neg params[`rate] * params[`maturity]] * .opt.cdf params[`sign]*params[`d2] ) - params[`spot] * .opt.cdf params[`sign]*params[`d1]
 };

.opt.bls[`BlackScholesDistribution]:{[params] 
	d1:(log[params[`spot] % params[`strike]] + params[`maturity]*params[`rate]+0.5*params[`volatility] xexp 2) %params[`volatility] * sqrt params[`maturity];
	d2:d1-params[`volatility] *sqrt params[`maturity];
	sign:neg 1f-2f*`float$(params[`direct]=`call );
	(`d1`d2`sign)!(d1;d2;sign)
    };
	

/ Bjerksund - Stensland Approximation formula for pricing American options
/ "The Complete Guide to Option Pricing Formulas " (2006), by Espen Haug
.opt.bls[`BjerksundStenslandFormula]:{[params]
	/ prepare params
	d:params[`yield];
	sigma:params[`volatility];
	r:params[`rate];
	T:params[`maturity];

	/ Drift rate or cost of carry
	b:r-d;

	/ put-call transformation
	$[`put=params[`direct];
		[r:r-b;
		b:neg b;
		K:params[`spot];
		S:params[`strike]];
		[S:params[`spot];
		K:params[`strike]];
	];
	sig2:sigma xexp 2;

	/ Barrier options step 1
	beta:(0.5-(b%sig2))+sqrt[(((b%sig2)-0.5)xexp 2)+((2*r)%sig2)];
        
	/ define boundaries
	b_inf:K*(beta%(beta-1f));
	b0:max[(K;K*(r%(r-b)))];
        
	/ Bjerksund-Stensland closed form boundary
	h:neg(((b*T)+(2*sigma*sqrt[T]))*(K xexp 2)%(b0*(b_inf-b0)));
	X:b0+((1-exp(h))*(b_inf-b0));
	
	alpha:(X-K)*(X xexp neg beta);

	/ Bjerksund-Stensland approximation
	tmp1:alpha * S xexp beta;
	tmp2:alpha * .opt.bls[`BjerksundStenslandPhi] @ (`spot`maturity`gamma`H`X`sigma`rate`b)!(S;T;beta;X;X;sigma;r;b);
	tmp3:.opt.bls[`BjerksundStenslandPhi] @ (`spot`maturity`gamma`H`X`sigma`rate`b)!(S;T;1;X;X;sigma;r;b);
	tmp4:.opt.bls[`BjerksundStenslandPhi] @ (`spot`maturity`gamma`H`X`sigma`rate`b)!(S;T;1;K;X;sigma;r;b);
	tmp5:K*.opt.bls[`BjerksundStenslandPhi] @ (`spot`maturity`gamma`H`X`sigma`rate`b)!(S;T;0;X;X;sigma;r;b);
	tmp6:K*.opt.bls[`BjerksundStenslandPhi] @ (`spot`maturity`gamma`H`X`sigma`rate`b)!(S;T;0;K;X;sigma;r;b);

	/ Equivalent to tmp1-tmp2+tmp3-tmp4-tmp5+tmp6
	:(tmp1-tmp2)+(tmp3-tmp4)+(tmp6-tmp5);
 }

/ Bjerksund - Stensland Phi function
.opt.bls[`BjerksundStenslandPhi]:{[params]
	/ prepare params
	T:params[`maturity];
	S:params[`spot];
	gamma:params[`gamma];
	sigma:params[`sigma];
	H:params[`H];
	X:params[`X];
	r:params[`rate];
	b:params[`b];

	sig2:sigma xexp 2;
	k:((2*b)%sig2)+((2*gamma)-1);
	lambda:(neg r)+(gamma*b)+(0.5*gamma)*(gamma-1)*sig2;
	tmp1:(T*(sig2*(gamma-0.5)+b)+log[S%H])%(sigma*sqrt[T]);
	tmp2:(T*(sig2*(gamma-0.5)+b)+log[(X xexp 2)%(S*H)])%(sigma*sqrt[T]);

	N_tmp1:.opt.cdf[neg 0f^tmp1];
	N_tmp2:.opt.cdf[neg 0f^tmp2];

	:exp[lambda*T]*(S xexp gamma)*(N_tmp1 - (((X%S) xexp k)*N_tmp2));
	
 };
 
 

/ Turnbull - Wakeman Approximation formula for pricing Asian options (specifically futures)
/ "The Complete Guide to Option Pricing Formulas " (2006), by Espen Haug
/ "Asian Options with Cost Of Carry Zero" (2006), by Espen Haug
.opt.bls[`TurnbullWakemanFormula]:{[params]
    / prepare params
    T:params[`maturity];
    spot:params[`spot];
    strike:params[`strike];
    rate:params[`rate];
    yield:params[`yield];
    avgSpot:params[`averageSpot];
    sigma:params[`volatility];
    T2:params[`pricingStart];

    if[0 = T;
        if[params[`direct]=`call;
            :max[(params[`averageSpot] - params[`strike];0f)];
        ];
        :max[(params[`strike] - params[`averageSpot];0f)];
    ];
    

    if[T2>=0;
        if[T<>T2;
            sigma:sqrt[log[(.opt.bls[`EvaluateM2][params])%((.opt.bls[`EvaluateM1][params])xexp 2)]%T];
            yield:rate - log[.opt.bls[`EvaluateM1][params]]%T;
            :.opt.bls[`CalculateOptionPrice][(`type_`spot`strike`maturity`volatility`rate`direct)!(`European;spot;strike;T;sigma;rate;params[`direct])];
        ];
        :.opt.bls[`CalculateOptionPrice][(`type_`spot`strike`maturity`volatility`rate`direct)!(`European;spot;strike;T;sigma;rate;params[`direct])];
    ];

    if[T2<0;
        params[`strike]:((strike*(T-T2)) + T2*avgSpot)%T;
        if[strike>0;
            params[`averageSpot]:0f;
            params[`pricingStart]:0f;
            :(T%(T-T2))*.opt.bls[`TurnbullWakemanFormula][params];
        ];
        if[params[`direct]=`call;
            params[`pricingStart]:0f;
            :(T%(T-T2))*(strike*(.opt.bls[`EvaluateM1][params]-params[`strike]))*exp[neg rate*T];
        ];
        if[params[`direct]=`put;
            :0f;
        ];
    ];
 };
 
.opt.bls[`EvaluateM1]:{[params]
    $[params[`rate]<>params[`yield];
        :(exp[params[`rate]- params[`yield]]*params[`maturity]) - (exp[params[`rate]- params[`yield]]*params[`pricingStart])%((params[`rate]-params[`yield])*(params[`maturity]-params[`pricingStart]));
        :1f;
    ];
 };

.opt.bls[`EvaluateM2]:{[params]
    $[params[`rate]<>params[`yield];
        :(2*exp[(2*(params[`rate]-params[`yield])+(params[`volatility]xexp 2))*params[`maturity]])%((params[`rate]-params[`yield] + params[`volatility]xexp 2)*(2*(params[`rate]-params[`yield]))+((params[`maturity] - params[`pricingStart])xexp 2)) + ((2*exp[(2*(params[`rate]-params[`yield])+(params[`volatility]xexp 2))*params[`pricingStart]])% (((params[`rate]-params[`yield])*(params[`volatility]xexp 2))*(1%(2*(params[`rate]-params[`yield])+(params[`volatility]xexp 2))-exp[(params[`rate]-params[`yield])*(params[`maturity] - params[`pricingStart])]%((params[`rate]-params[`yield])+(params[`volatility]xexp 2)))));
        :(2%((params[`volatility]xexp 4)*((params[`maturity]-params[`pricingStart])xexp 2)))*(exp[(params[`volatility]xexp 2)*params[`maturity]]-exp[(params[`volatility]xexp 2)*params[`pricingStart]]*(1+(params[`volatility]xexp 2)*(params[`maturity]-params[`pricingStart])));
    ];
 };

