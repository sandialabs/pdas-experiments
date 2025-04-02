#ifndef PDAS_EXPERIMENTS_FOM_HPP_
#define PDAS_EXPERIMENTS_FOM_HPP_

#include "pressio/ode_steppers_implicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "observer.hpp"
#include <chrono>

template<class AppType, class ParserType>
void run_mono_fom(AppType & system, ParserType & parser)
{
    PRESSIOLOG_INITIALIZE(parser.loglevel(), parser.logtarget(), parser.logfile());

    using app_t = AppType;
    using scalar_t = typename app_t::scalar_type;
    using state_t = typename app_t::state_type;
    using jacob_t = typename app_t::jacobian_type;

    // initial condition
    state_t state = system.initialCondition();
    std::string icFile = parser.icFile();
    if (!(icFile.empty())) {
        // load from file
        auto instate = pschwarz::read_vector_from_binary<scalar_t>(icFile);
        int nrows = instate.rows();
        if (nrows == state.rows()) {
            state = instate;
        }
        else {
            throw std::runtime_error("Invalid icFile dimensions: " + std::to_string(nrows));
        }
    }

    const auto odeScheme = parser.odeScheme();
    auto stepperObj = pressio::ode::create_implicit_stepper(odeScheme, system);

    using lin_solver_t = pressio::linearsolvers::Solver<
        pressio::linearsolvers::iterative::Bicgstab, jacob_t>;
    lin_solver_t linSolverObj;
    auto NonLinSolver = pressio::create_newton_solver(stepperObj, linSolverObj);
    NonLinSolver.setStopCriterion(pressio::nonlinearsolvers::Stop::WhenAbsolutel2NormOfCorrectionBelowTolerance);
    NonLinSolver.setStopTolerance(1e-5);

    StateObserver Obs(parser.stateSamplingFreq());
    RuntimeObserver Obs_run("runtime.bin");

    const auto startTime = static_cast<scalar_t>(0.0);
    auto runtimeStart = std::chrono::high_resolution_clock::now();
    pressio::ode::advance_n_steps(
        stepperObj, state, startTime,
        parser.timeStepSize(),
        pressio::ode::StepCount(parser.numSteps()),
        Obs, NonLinSolver);
    auto runtimeEnd = std::chrono::high_resolution_clock::now();
    auto nsElapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(runtimeEnd - runtimeStart).count();
    double secElapsed = static_cast<double>(nsElapsed) * 1e-9;
    Obs_run(secElapsed);

    PRESSIOLOG_FINALIZE();

}

#endif