# candidate delta and p values to run the model
delta_list = [ 0.00001, 0.00001468, 0.00002154, 0.00003162, 0.00004642, 0.00006813, ...
               0.0001,  0.0001468,  0.0002154,  0.0003162,  0.0004642,  0.0006813, ...
               0.001,   0.001468,   0.002154,   0.003162,   0.004642,   0.006813, ...
               0.01,    0.01468,    0.02154,    0.03162,    0.04642,    0.06813, ...
               0.1,     0.1468,     0.2154,     0.3162,     0.4642,     0.6813];
p_list = 0:0.2:2;

for random_seed = [0]
    for delta = delta_list
        for p = p_list
            test_rcv1(p, delta, random_seed);
	    test_covtype(p, delta, random_seed);
        end
    end
end
