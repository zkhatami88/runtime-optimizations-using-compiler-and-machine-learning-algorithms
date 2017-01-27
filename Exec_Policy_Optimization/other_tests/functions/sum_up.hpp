auto f = [](unsigned i) -> unsigned
{
    unsigned sum = 0;
    for (int j = 0; j < 100000; ++j)
        sum += j * i;
 
    return sum;
};

int num_iters = 10000000;

int hpx_main(int argc, char* argv[])
{
    auto r = boost::irange(0, num_iters);
        
    hpx::util::high_resolution_timer t1;
    
    hpx::parallel::for_each_n(hpx::parallel::seq, r.begin(), num_iters, f);
    
    std::cout << t1.elapsed() << " ";
    
    hpx::util::high_resolution_timer t2;

    hpx::parallel::for_each_n(hpx::parallel::par, r.begin(), num_iters, f);

    std::cout << t2.elapsed() << std::endl;

    return hpx::finalize();
}