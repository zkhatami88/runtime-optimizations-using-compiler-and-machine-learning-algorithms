
//  main.cpp
//  Struct1
//
//  Created by Zahra Khatami.
//  Copyright (c) 2016 Zahra Khatami. All rights reserved.

//plot results and compare them
//if you can plot them on something that shows moving

#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>
#include <hpx/version.hpp>
#include <hpx/include/parallel_algorithm.hpp>
#include <hpx/include/parallel_executors.hpp>
#include <hpx/include/parallel_transform.hpp>
#include <hpx/include/parallel_executor_parameters.hpp>
#include <hpx/include/iostreams.hpp>
#include <hpx/include/threads.hpp>
#include <hpx/util/safe_lexical_cast.hpp>

#include <hpx/parallel/util/numa_allocator.hpp>

#include <boost/format.hpp>
#include <boost/range/functions.hpp>

#include <string>
#include <vector>

///////////////////////////////////////////////////////////////////////////////

double G=6.673*pow(10.0,-11.0);

std::vector<double> a1, b1;
std::vector<double> a2, b2;
std::vector<double> a3, b3;
std::vector<double> a4, b4;
std::vector<double> a5, b5;

///////////////////////////////////////////////////////////////////////////////
struct Body
{
public:
    int ID1,parent,m1;
    std::vector<double> r1, v1, force;

    friend class hpx::serialization::access;
    template<typename Ar> void serialize(Ar &ar, unsigned){
     ar &ID1 &parent &m1 &r1 &v1 &force;}
}body1;

std::vector<Body> b;
template<typename Ar> void serialize(Ar &ar, unsigned){
        ar & b;}

///////////////////////////////////////////////////////////////////////////////
void initialize_body(std::size_t N)
{
    body1.ID1=0; body1.m1=0; body1.parent=0;
    for(int i=0; i<3; ++i){
        body1.r1.push_back(1.0); body1.v1.push_back(1.0); body1.force.push_back(0.0);
    }
    for(int i=0; i<N; ++i)
        b.push_back(body1);
}

///////////////////////////////////////////////////////////////////////////////
void read_Input(std::size_t n)
{
    std::size_t N=n;
    for(std::size_t i=0; i<N; ++i){
        b[i].ID1=i; b[i].m1=1;
        b[i].parent=0;
        for(int j=0; j<3; ++j){
            b[i].r1[j]=float(std::rand() % 100);
            b[i].v1[j]=0;
        }
    }
}
///////////////////////////////////////////////////////////////////////////////

void copy_input(std::size_t n)
{

    std::size_t j, N=n;
    for(std::size_t i=0; i<N; ++i)
    {
        if(i<(N/5))
        {
            for(j=0; j<3; ++j)
            {
                a1.push_back(b[i].force[j]);
                b1.push_back(b[i].r1[j]);
            }
        }

        if(i>=(N/5) && i<(2*N/5))
        {
            for(j=0; j<3; ++j)
            {
                a2.push_back(b[i].force[j]);
                b2.push_back(b[i].r1[j]);
            }
        }

        if(i>=(2*N/5) && i<(3*N/5))
        {
            for(j=0; j<3; ++j)
            {
                a3.push_back(b[i].force[j]);
                b3.push_back(b[i].r1[j]);
            }
        }

        if(i>=(3*N/5) && i<(4*N/5))
        {
            for(j=0; j<3; ++j)
            {
                a4.push_back(b[i].force[j]);
                b4.push_back(b[i].r1[j]);
            }
        }

        if(i>=(4*N/5) && i<N)
        {
            for(j=0; j<3; ++j)
            {
                a5.push_back(b[i].force[j]);
                b5.push_back(b[i].r1[j]);
            }
        }
}}

///////////////////////////////////////////////////////////////////////////////
void N_Body(std::size_t n)
{
    initialize_body(n);
    read_Input(n);
    copy_input(n);
}

///////////////////////////////////////////////////////////////////////////////
hpx::threads::topology& retrieve_topology()
{
    static hpx::threads::topology& topo = hpx::threads::create_topology();
    return topo;
}

///////////////////////////////////////////////////////////////////////////////
double mysecond()
{
    return hpx::util::high_resolution_clock::now() * 1e-9;
}

///////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<double> >
numa_domain_worker(std::size_t domain, hpx::lcos::local::latch& l,
    std::size_t part_size, std::size_t iterations,
    std::size_t prefetch_distance_factor, std::size_t range_size)
{
    l.count_down_and_wait();

    std::vector<std::vector<double> > timing(2, std::vector<double>(iterations));

    //--------------------------------------
    std::vector<std::size_t> range(range_size);
        for(std::size_t i=0; i<range_size; ++i)
            range[i]=i;

    auto ctx = hpx::parallel::util::make_prefetcher_context(range.begin(),
        range.end(), prefetch_distance_factor,
        a1, b1, a2, b2, a3, b3, a4, b4, a5, b5);

    for(std::size_t it=0 ; it!=iterations; ++it)
    {
        //random access iterator
        timing[0][it] = mysecond();
        hpx::parallel::for_each(hpx::parallel::par,
            range.begin(), range.end(),
            [&](std::size_t i)
            {
                a1[6*i] = a1[6*i] + G * (b1[6*i + 3]-b[0].r1[0])/pow((1+pow((b1[6*i + 3]-b[0].r1[0]),2.0)+
                    pow((b1[6*i + 4]-b[0].r1[1]),2) +pow((b1[6*i + 5] -b[0].r1[2]),2.0)),1.5);

                a1[6*i + 1] = a1[6*i + 1] + G * (b1[6*i + 4]-b[0].r1[0])/pow((1+pow((b1[6*i + 4]-b[0].r1[0]),2.0)+
                    pow((b1[6*i + 4]-b[0].r1[1]),2) +pow((b1[6*i + 5] -b[0].r1[2]),2.0)),1.5);

                a1[6*i + 2] = a1[6*i + 2] + G * (b1[6*i + 5]-b[0].r1[0])/pow((1+pow((b1[6*i + 5]-b[0].r1[0]),2.0)+
                    pow((b1[6*i + 4]-b[0].r1[1]),2) +pow((b1[6*i + 5] -b[0].r1[2]),2.0)),1.5);


                a2[6*i] = a2[6*i] + G * (b2[6*i + 3]-b[0].r1[0])/pow((1+pow((b2[6*i + 3]-b[0].r1[0]),2.0)+
                    pow((b2[6*i + 4]-b[0].r1[1]),2) +pow((b2[6*i + 5] -b[0].r1[2]),2.0)),1.5);

                a2[6*i + 1] = a2[6*i + 1] + G * (b2[6*i + 4]-b[0].r1[0])/pow((1+pow((b2[6*i + 4]-b[0].r1[0]),2.0)+
                    pow((b2[6*i + 4]-b[0].r1[1]),2) +pow((b2[6*i + 5] -b[0].r1[2]),2.0)),1.5);

                a2[6*i + 2] = a2[6*i + 2] + G * (b2[6*i + 5]-b[0].r1[0])/pow((1+pow((b2[6*i + 5]-b[0].r1[0]),2.0)+
                    pow((b2[6*i + 4]-b[0].r1[1]),2) +pow((b2[6*i + 5] -b[0].r1[2]),2.0)),1.5);


                a3[6*i] = a3[6*i] + G * (b3[6*i + 3]-b[0].r1[0])/pow((1+pow((b3[6*i + 3]-b[0].r1[0]),2.0)+
                    pow((b3[6*i + 4]-b[0].r1[1]),2) +pow((b3[6*i + 5] -b[0].r1[2]),2.0)),1.5);

                a3[6*i + 1] = a3[6*i + 1] + G * (b3[6*i + 4]-b[0].r1[0])/pow((1+pow((b3[6*i + 4]-b[0].r1[0]),2.0)+
                    pow((b3[6*i + 4]-b[0].r1[1]),2) +pow((b3[6*i + 5] -b[0].r1[2]),2.0)),1.5);

                a3[6*i + 2] = a3[6*i + 2] + G * (b3[6*i + 5]-b[0].r1[0])/pow((1+pow((b3[6*i + 5]-b[0].r1[0]),2.0)+
                    pow((b3[6*i + 4]-b[0].r1[1]),2) +pow((b3[6*i + 5] -b[0].r1[2]),2.0)),1.5);


                a4[6*i] = a4[6*i] + G * (b4[6*i + 3]-b[0].r1[0])/pow((1+pow((b4[6*i + 3]-b[0].r1[0]),2.0)+
                    pow((b4[6*i + 4]-b[0].r1[1]),2) +pow((b4[6*i + 5] -b[0].r1[2]),2.0)),1.5);

                a4[6*i + 1] = a4[6*i + 1] + G * (b4[6*i + 4]-b[0].r1[0])/pow((1+pow((b4[6*i + 4]-b[0].r1[0]),2.0)+
                    pow((b4[6*i + 4]-b[0].r1[1]),2) +pow((b4[6*i + 5] -b[0].r1[2]),2.0)),1.5);

                a4[6*i + 2] = a4[6*i + 2] + G * (b4[6*i + 5]-b[0].r1[0])/pow((1+pow((b4[6*i + 5]-b[0].r1[0]),2.0)+
                    pow((b4[6*i + 4]-b[0].r1[1]),2) +pow((b4[6*i + 5] -b[0].r1[2]),2.0)),1.5);

                a5[6*i] = a5[6*i] + G * (b5[6*i + 3]-b[0].r1[0])/pow((1+pow((b5[6*i + 3]-b[0].r1[0]),2.0)+
                    pow((b5[6*i + 4]-b[0].r1[1]),2) +pow((b5[6*i + 5] -b[0].r1[2]),2.0)),1.5);

                a5[6*i + 1] = a5[6*i + 1] + G * (b5[6*i + 4]-b[0].r1[0])/pow((1+pow((b5[6*i + 4]-b[0].r1[0]),2.0)+
                    pow((b5[6*i + 4]-b[0].r1[1]),2) +pow((b5[6*i + 5] -b[0].r1[2]),2.0)),1.5);

                a5[6*i + 2] = a5[6*i + 2] + G * (b5[6*i + 5]-b[0].r1[0])/pow((1+pow((b5[6*i + 5]-b[0].r1[0]),2.0)+
                    pow((b5[6*i + 4]-b[0].r1[1]),2) +pow((b5[6*i + 5] -b[0].r1[2]),2.0)),1.5);

            }
        );

        timing[0][it] = mysecond() - timing[0][it];


        timing[1][it] = mysecond();

        //prefetching iterator
        hpx::parallel::for_each(hpx::parallel::par,
            ctx.begin(), ctx.end(),
            [&](std::size_t i)
            {
                a1[6*i] = a1[6*i] + G * (b1[6*i + 3]-b[0].r1[0])/pow((1+pow((b1[6*i + 3]-b[0].r1[0]),2.0)+
                    pow((b1[6*i + 4]-b[0].r1[1]),2) +pow((b1[6*i + 5] -b[0].r1[2]),2.0)),1.5);

                a1[6*i + 1] = a1[6*i + 1] + G * (b1[6*i + 4]-b[0].r1[0])/pow((1+pow((b1[6*i + 4]-b[0].r1[0]),2.0)+
                    pow((b1[6*i + 4]-b[0].r1[1]),2) +pow((b1[6*i + 5] -b[0].r1[2]),2.0)),1.5);

                a1[6*i + 2] = a1[6*i + 2] + G * (b1[6*i + 5]-b[0].r1[0])/pow((1+pow((b1[6*i + 5]-b[0].r1[0]),2.0)+
                    pow((b1[6*i + 4]-b[0].r1[1]),2) +pow((b1[6*i + 5] -b[0].r1[2]),2.0)),1.5);


                a2[6*i] = a2[6*i] + G * (b2[6*i + 3]-b[0].r1[0])/pow((1+pow((b2[6*i + 3]-b[0].r1[0]),2.0)+
                    pow((b2[6*i + 4]-b[0].r1[1]),2) +pow((b2[6*i + 5] -b[0].r1[2]),2.0)),1.5);

                a2[6*i + 1] = a2[6*i + 1] + G * (b2[6*i + 4]-b[0].r1[0])/pow((1+pow((b2[6*i + 4]-b[0].r1[0]),2.0)+
                    pow((b2[6*i + 4]-b[0].r1[1]),2) +pow((b2[6*i + 5] -b[0].r1[2]),2.0)),1.5);

                a2[6*i + 2] = a2[6*i + 2] + G * (b2[6*i + 5]-b[0].r1[0])/pow((1+pow((b2[6*i + 5]-b[0].r1[0]),2.0)+
                    pow((b2[6*i + 4]-b[0].r1[1]),2) +pow((b2[6*i + 5] -b[0].r1[2]),2.0)),1.5);


                a3[6*i] = a3[6*i] + G * (b3[6*i + 3]-b[0].r1[0])/pow((1+pow((b3[6*i + 3]-b[0].r1[0]),2.0)+
                    pow((b3[6*i + 4]-b[0].r1[1]),2) +pow((b3[6*i + 5] -b[0].r1[2]),2.0)),1.5);

                a3[6*i + 1] = a3[6*i + 1] + G * (b3[6*i + 4]-b[0].r1[0])/pow((1+pow((b3[6*i + 4]-b[0].r1[0]),2.0)+
                    pow((b3[6*i + 4]-b[0].r1[1]),2) +pow((b3[6*i + 5] -b[0].r1[2]),2.0)),1.5);

                a3[6*i + 2] = a3[6*i + 2] + G * (b3[6*i + 5]-b[0].r1[0])/pow((1+pow((b3[6*i + 5]-b[0].r1[0]),2.0)+
                    pow((b3[6*i + 4]-b[0].r1[1]),2) +pow((b3[6*i + 5] -b[0].r1[2]),2.0)),1.5);


                a4[6*i] = a4[6*i] + G * (b4[6*i + 3]-b[0].r1[0])/pow((1+pow((b4[6*i + 3]-b[0].r1[0]),2.0)+
                    pow((b4[6*i + 4]-b[0].r1[1]),2) +pow((b4[6*i + 5] -b[0].r1[2]),2.0)),1.5);

                a4[6*i + 1] = a4[6*i + 1] + G * (b4[6*i + 4]-b[0].r1[0])/pow((1+pow((b4[6*i + 4]-b[0].r1[0]),2.0)+
                    pow((b4[6*i + 4]-b[0].r1[1]),2) +pow((b4[6*i + 5] -b[0].r1[2]),2.0)),1.5);

                a4[6*i + 2] = a4[6*i + 2] + G * (b4[6*i + 5]-b[0].r1[0])/pow((1+pow((b4[6*i + 5]-b[0].r1[0]),2.0)+
                    pow((b4[6*i + 4]-b[0].r1[1]),2) +pow((b4[6*i + 5] -b[0].r1[2]),2.0)),1.5);

                a5[6*i] = a5[6*i] + G * (b5[6*i + 3]-b[0].r1[0])/pow((1+pow((b5[6*i + 3]-b[0].r1[0]),2.0)+
                    pow((b5[6*i + 4]-b[0].r1[1]),2) +pow((b5[6*i + 5] -b[0].r1[2]),2.0)),1.5);

                a5[6*i + 1] = a5[6*i + 1] + G * (b5[6*i + 4]-b[0].r1[0])/pow((1+pow((b5[6*i + 4]-b[0].r1[0]),2.0)+
                    pow((b5[6*i + 4]-b[0].r1[1]),2) +pow((b5[6*i + 5] -b[0].r1[2]),2.0)),1.5);

                a5[6*i + 2] = a5[6*i + 2] + G * (b5[6*i + 5]-b[0].r1[0])/pow((1+pow((b5[6*i + 5]-b[0].r1[0]),2.0)+
                    pow((b5[6*i + 4]-b[0].r1[1]),2) +pow((b5[6*i + 5] -b[0].r1[2]),2.0)),1.5);

            }
        );


        timing[1][it] = mysecond() - timing[1][it];

    }

    return timing;

}


///////////////////////////////////////////////////////////////////////////////
std::size_t get_num_numa_nodes(hpx::threads::topology const& topo,
    boost::program_options::variables_map& vm)
{
    std::size_t numa_nodes = topo.get_number_of_numa_nodes();
    if (numa_nodes == 0)
        numa_nodes = topo.get_number_of_sockets();

    std::string num_numa_domains_str = vm["numa-domains"].as<std::string>();
    if (num_numa_domains_str != "all")
    {
        numa_nodes = hpx::util::safe_lexical_cast<std::size_t>(num_numa_domains_str);
    }
    return numa_nodes;
}

std::pair<std::size_t, std::size_t> get_num_numa_pus(
    hpx::threads::topology const& topo, std::size_t numa_nodes,
    boost::program_options::variables_map& vm)
{
    std::size_t numa_pus = hpx::threads::hardware_concurrency() / numa_nodes;

    std::string num_threads_str = vm["threads"].as<std::string>();
    std::size_t pus = numa_pus;

    if(num_threads_str != "all")
    {
        pus = hpx::util::safe_lexical_cast<std::size_t>(num_threads_str);
    }

    return std::make_pair(numa_pus, pus);
}

///////////////////////////////////////////////////////////////////////////////

int hpx_main(boost::program_options::variables_map& vm)
{

    // extract hardware topology
    hpx::threads::topology const& topo = retrieve_topology();

    std::size_t numa_nodes = get_num_numa_nodes(topo, vm);
    std::pair<std::size_t, std::size_t> pus =
        get_num_numa_pus(topo, numa_nodes, vm);

    std::size_t iterations = vm["iterations"].as<std::size_t>();
    std::string num_numa_domains_str = vm["numa-domains"].as<std::string>();
    std::size_t prefetch_distance_factor = vm["prefetch_distance_factor"].as<std::size_t>();
    std::size_t range_size = vm["range_size"].as<std::size_t>();
    std::size_t problem_size = vm["problem_size"].as<std::size_t>();

    using namespace hpx::parallel;

    typedef hpx::threads::executors::local_priority_queue_attached_executor
        executor_type;
    typedef std::vector<executor_type> executors_vector;

    executors_vector execs;
    execs.reserve(numa_nodes);

    // creating our executors ....
    for (std::size_t i = 0; i != numa_nodes; ++i)
    {
        // create executor for this NUMA domain
        execs.emplace_back(i * pus.first, pus.second);
    }

    // allocate data
    typedef hpx::parallel::util::numa_allocator<
            double, executors_vector
        > allocator_type;
    allocator_type alloc(execs, retrieve_topology());

    // perform benchmark
    hpx::lcos::local::latch l(numa_nodes);

    double time_total = mysecond();
    std::vector<hpx::future<std::vector<std::vector<double> > > > workers;
    workers.reserve(numa_nodes);

    std::size_t part_size = problem_size / numa_nodes;

    N_Body(problem_size);

    for (std::size_t i = 0; i != numa_nodes; ++i)
    {
        auto policy = par.on(execs[i]);
        workers.push_back(
            hpx::async(execs[i], &numa_domain_worker,
                i, boost::ref(l),
                part_size, iterations, prefetch_distance_factor, range_size)
            );
    }

    std::vector<std::vector<std::vector<double> > >
        timings_all = hpx::util::unwrapped(workers);

    time_total = mysecond() - time_total;


    const char *label[2] = {
        "Force Computation:             ",
        "Force Computation_WPrefetching:",
    };

    const double bytes[2] = {
        15 * 10 * sizeof(double) * static_cast<double>(range_size),
        15 * 10 * sizeof(double) * static_cast<double>(range_size)
    };


    std::vector<std::vector<double> > timing(2, std::vector<double>(iterations, 0.0));

    for(auto const & times : timings_all)
    {
        for(std::size_t iteration = 0; iteration != iterations; ++iteration)
        {
            timing[0][iteration] += times[0][iteration];
            timing[1][iteration] += times[1][iteration];
        }
    }

    for(std::size_t iteration = 0; iteration != iterations; ++iteration)
    {
        timing[0][iteration] /= numa_nodes;
        timing[1][iteration] /= numa_nodes;
    }

    // Note: skip first iteration
    std::vector<double> avgtime(2, 0.0);
    std::vector<double> mintime(2, (std::numeric_limits<double>::max)());
    std::vector<double> maxtime(2, 0.0);
    for(std::size_t iteration = 1; iteration != iterations; ++iteration)
    {
        for (std::size_t j=0; j<2; j++){
            avgtime[j] = avgtime[j] + timing[j][iteration];
            mintime[j] = (std::min)(mintime[j], timing[j][iteration]);
            maxtime[j] = (std::max)(maxtime[j], timing[j][iteration]);
        }
    }


    std::cout<< "-------------------------------------------------------------\n";

    printf("Function                      Best Rate MB/s    Avg time     Min time     Max time\n");
    for (std::size_t j=0; j<2; j++)
    {
        avgtime[j] = avgtime[j]/(double)(iterations-1);

        printf("%s%12.1f  %11.6f  %11.6f  %11.6f\n", label[j],
            1.0E-06 * bytes[j]/mintime[j],
            avgtime[j],
            mintime[j],
            maxtime[j]);
    }


    std::cout
        << "\nTotal time: " << time_total
        << " (per iteration: " << time_total/iterations << ")\n";

    std::cout<< "-------------------------------------------------------------\n";


    return hpx::finalize();
}

int main(int argc, char* argv[])
{

    using namespace boost::program_options;

    options_description cmdline("usage: " HPX_APPLICATION_STRING " [options]");

    cmdline.add_options()
        (   "iterations",
            boost::program_options::value<std::size_t>()->default_value(100),
            "iterations (default: 1000)")
        (   "threads",
            boost::program_options::value<std::string>()->default_value("all"),
            "number of threads per NUMA domain to use. (default: all)")
        (   "numa-domains",
            boost::program_options::value<std::string>()->default_value("all"),
            "number of NUMA domains to use. (default: all)")
        (   "prefetch_distance_factor",
            boost::program_options::value<std::size_t>()->default_value(20),
            "Distance (in chunk_size) between each preteching data. (default: 1)")
        (   "range_size",
            boost::program_options::value<std::size_t>()->default_value(100000),
            "size of range. (default: 100000000)")
        (   "problem_size",
            boost::program_options::value<std::size_t>()->default_value(100000),
            "size of problem. (default: 100000000)")
    ;

    // parse command line here to extract the necessary settings for HPX
    parsed_options opts =
        command_line_parser(argc, argv)
            .allow_unregistered()
            .options(cmdline)
            .style(command_line_style::unix_style)
            .run();

    variables_map vm;
    store(opts, vm);

    hpx::threads::topology const& topo = retrieve_topology();
    std::size_t numa_nodes = get_num_numa_nodes(topo, vm);
    std::pair<std::size_t, std::size_t> pus =
        get_num_numa_pus(topo, numa_nodes, vm);
    std::size_t num_cores = topo.get_number_of_numa_node_cores(0);

    std::vector<std::string> cfg;
    cfg.push_back("hpx.numa_sensitive=2");  // no-cross NUMA stealing

    // block all cores of requested number of NUMA-domains
    cfg.push_back(boost::str(
        boost::format("hpx.cores=%d") % (numa_nodes * num_cores)
    ));
    cfg.push_back(boost::str(
        boost::format("hpx.os_threads=%d") % (numa_nodes * pus.second)
    ));

    std::string node_name("numanode");
    if (topo.get_number_of_numa_nodes() == 0)
        node_name = "socket";

    std::string bind_desc("hpx.bind!=");
    for (std::size_t i = 0; i != numa_nodes; ++i)
    {
        if (i != 0)
            bind_desc += ";";

        std::size_t base_thread = i * pus.second;
        bind_desc += boost::str(
            boost::format("thread:%d-%d=%s:%d.core:0-%d.pu:0")
              % base_thread % (base_thread+pus.second-1)  // thread:%d-%d
              % node_name % i                             // %s:%d
              % (pus.second-1)                            // core:0-%d
        );
    }
    cfg.push_back(bind_desc);

    return hpx::init(cmdline, argc, argv, cfg);

}
