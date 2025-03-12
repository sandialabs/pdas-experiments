
#ifndef PDAS_EXPERIMENTS_OBSERVER_HPP_
#define PDAS_EXPERIMENTS_OBSERVER_HPP_

class StateObserver
{
public:
    StateObserver(const std::string & f0, int freq)
        : myfile_(f0,  std::ios::out | std::ios::binary),
        sampleFreq_(freq){}

    explicit StateObserver(int freq)
        : myfile_("state_snapshots.bin",  std::ios::out | std::ios::binary),
        sampleFreq_(freq){}

    StateObserver() = default;
    StateObserver & operator=(StateObserver &) = default;
    StateObserver & operator=(StateObserver &&) = default;

    ~StateObserver(){ myfile_.close(); }

    template<typename TimeType, typename ObservableType>
    std::enable_if_t< pressio::is_vector_eigen<ObservableType>::value >
    operator()(pressio::ode::StepCount step,
            const TimeType /*timeIn*/,
            const ObservableType & state)
    {
        if (step.get() % sampleFreq_ == 0) {
            const std::size_t ext = state.size()*sizeof(typename ObservableType::Scalar);
            myfile_.write(reinterpret_cast<const char*>(&state(0)), ext);
        }
    }

private:
    std::ofstream myfile_;
    int sampleFreq_ = {};
};

class RuntimeObserver
{
public:
    RuntimeObserver(const std::string & f0)
        : timeFile_(f0, std::ios::out | std::ios::binary)
    {}

    ~RuntimeObserver() { timeFile_.close(); }

    void operator() (double runtimeVal, int itersVal)
    {
        // number of subiterations, total runtime for single outer loop iteration
        std::size_t iters_st = static_cast<std::size_t>(itersVal);
        timeFile_.write(reinterpret_cast<const char*>(&iters_st), sizeof(std::size_t));
        timeFile_.write(reinterpret_cast<const char*>(&runtimeVal), sizeof(double));
    }

    void operator() (double runtimeVal)
    {
        // Overload for handling non-Schwarz input, have to manually time it and pass total runtime
        // Just write it as if there's one iteration, one subiteration
        std::size_t one = static_cast<std::size_t>(1);
        timeFile_.write(reinterpret_cast<const char*>(&one), sizeof(std::size_t));
        timeFile_.write(reinterpret_cast<const char*>(&runtimeVal), sizeof(double));
    }

private:
    std::ofstream timeFile_;
};

#endif
