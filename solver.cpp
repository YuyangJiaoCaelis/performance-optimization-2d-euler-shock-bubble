#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

#if defined(USE_VECTOR_STATE_STORAGE)
constexpr bool kUseVectorStateStorage = true;
constexpr const char* kStateStorageLabel = "std::vector";
#else
constexpr bool kUseVectorStateStorage = false;
constexpr const char* kStateStorageLabel = "std::array";
#endif

#if defined(PASS_STATE_BY_VALUE)
constexpr bool kPassStateByValue = true;
constexpr const char* kStatePassingLabel = "value";
#else
constexpr bool kPassStateByValue = false;
constexpr const char* kStatePassingLabel = "reference";
#endif

#if defined(PASS_GRID_BY_VALUE)
constexpr bool kPassGridByValue = true;
constexpr const char* kGridPassingLabel = "value";
#else
constexpr bool kPassGridByValue = false;
constexpr const char* kGridPassingLabel = "reference";
#endif

#if defined(NON_CACHE_FRIENDLY_ACCESS)
constexpr bool kUseCacheUnfriendlyAccess = true;
constexpr const char* kAccessPatternLabel = "cache-unfriendly";
#else
constexpr bool kUseCacheUnfriendlyAccess = false;
constexpr const char* kAccessPatternLabel = "cache-friendly";
#endif

#if defined(DISABLE_RECONSTRUCTION_CACHE)
constexpr bool kUseReconstructionCache = false;
constexpr const char* kReconstructionCacheLabel = "disabled";
#else
constexpr bool kUseReconstructionCache = true;
constexpr const char* kReconstructionCacheLabel = "enabled";
#endif

template <std::size_t Count>
using FixedStorage = std::conditional_t<kUseVectorStateStorage, std::vector<double>, std::array<double, Count>>;

template <std::size_t Count>
FixedStorage<Count> makeStorage() {
    if constexpr (kUseVectorStateStorage) {
        return FixedStorage<Count>(Count, 0.0);
    } else {
        return {};
    }
}

class Primitive {
public:
    Primitive()
        : values_(makeStorage<4>()) {
    }

    Primitive(double density, double velocity_x, double velocity_y, double pressure)
        : Primitive() {
        this->density() = density;
        this->velocity_x() = velocity_x;
        this->velocity_y() = velocity_y;
        this->pressure() = pressure;
    }

    double& density() {
        return values_[0];
    }

    const double& density() const {
        return values_[0];
    }

    double& velocity_x() {
        return values_[1];
    }

    const double& velocity_x() const {
        return values_[1];
    }

    double& velocity_y() {
        return values_[2];
    }

    const double& velocity_y() const {
        return values_[2];
    }

    double& pressure() {
        return values_[3];
    }

    const double& pressure() const {
        return values_[3];
    }

private:
    FixedStorage<4> values_;
};

class Conserved {
public:
    Conserved()
        : values_(makeStorage<4>()) {
    }

    Conserved(double density, double momentum_x, double momentum_y, double energy)
        : Conserved() {
        this->density() = density;
        this->momentum_x() = momentum_x;
        this->momentum_y() = momentum_y;
        this->energy() = energy;
    }

    double& density() {
        return values_[0];
    }

    const double& density() const {
        return values_[0];
    }

    double& momentum_x() {
        return values_[1];
    }

    const double& momentum_x() const {
        return values_[1];
    }

    double& momentum_y() {
        return values_[2];
    }

    const double& momentum_y() const {
        return values_[2];
    }

    double& energy() {
        return values_[3];
    }

    const double& energy() const {
        return values_[3];
    }

private:
    FixedStorage<4> values_;
};

Primitive operator+(const Primitive& left, const Primitive& right) {
    return {
        left.density() + right.density(),
        left.velocity_x() + right.velocity_x(),
        left.velocity_y() + right.velocity_y(),
        left.pressure() + right.pressure()
    };
}

Primitive operator-(const Primitive& left, const Primitive& right) {
    return {
        left.density() - right.density(),
        left.velocity_x() - right.velocity_x(),
        left.velocity_y() - right.velocity_y(),
        left.pressure() - right.pressure()
    };
}

Primitive operator*(double scalar, const Primitive& state) {
    return {
        scalar * state.density(),
        scalar * state.velocity_x(),
        scalar * state.velocity_y(),
        scalar * state.pressure()
    };
}

Conserved operator+(const Conserved& left, const Conserved& right) {
    return {
        left.density() + right.density(),
        left.momentum_x() + right.momentum_x(),
        left.momentum_y() + right.momentum_y(),
        left.energy() + right.energy()
    };
}

Conserved operator-(const Conserved& left, const Conserved& right) {
    return {
        left.density() - right.density(),
        left.momentum_x() - right.momentum_x(),
        left.momentum_y() - right.momentum_y(),
        left.energy() - right.energy()
    };
}

Conserved operator*(double scalar, const Conserved& state) {
    return {
        scalar * state.density(),
        scalar * state.momentum_x(),
        scalar * state.momentum_y(),
        scalar * state.energy()
    };
}

using PrimitiveArgument = std::conditional_t<kPassStateByValue, Primitive, const Primitive&>;
using ConservedArgument = std::conditional_t<kPassStateByValue, Conserved, const Conserved&>;

enum class Direction {
    X,
    Y
};

enum class InitialCase {
    ShockBubble,
    PlanarShock
};

struct Settings {
    int nx = 500;
    int ny = 197;
    double final_time = 3.0e-4;
    double cfl = 0.45;
    double snapshot_interval = 0.0;
    std::string output_prefix = "output/shock_bubble";
    InitialCase initial_case = InitialCase::ShockBubble;
};

struct LineWorkspace {
    std::vector<Primitive> primitive_values;
    std::vector<Conserved> left_half_step;
    std::vector<Conserved> right_half_step;
    std::vector<Conserved> interface_fluxes;

    void resize(std::size_t length) {
        primitive_values.resize(length);
        left_half_step.resize(length);
        right_half_step.resize(length);
        interface_fluxes.resize(length);
    }
};

double minmod(double left, double right) {
    if (left * right <= 0.0) {
        return 0.0;
    }
    return std::copysign(std::min(std::abs(left), std::abs(right)), left);
}

int chooseScrambledStride(int count) {
    if (count <= 2) {
        return 1;
    }

    int stride = count / 2 + 1;
    while (std::gcd(stride, count) != 1) {
        ++stride;
    }
    return stride;
}

std::vector<int> buildIndexOrder(int first_index, int count) {
    std::vector<int> order;
    order.reserve(static_cast<std::size_t>(count));

    if (!kUseCacheUnfriendlyAccess) {
        for (int offset = 0; offset < count; ++offset) {
            order.push_back(first_index + offset);
        }
        return order;
    }

    const int stride = chooseScrambledStride(count);
    int position = 0;
    for (int step = 0; step < count; ++step) {
        order.push_back(first_index + position);
        position = (position + stride) % count;
    }

    return order;
}

std::string nextArgument(int& index, int argc, char** argv) {
    if (index + 1 >= argc) {
        throw std::runtime_error("Missing value after " + std::string(argv[index]) + ".");
    }
    ++index;
    return argv[index];
}

std::string snapshotFileName(const std::string& prefix, int snapshot_index) {
    std::ostringstream builder;
    builder << prefix << "_" << std::setw(4) << std::setfill('0') << snapshot_index << ".csv";
    return builder.str();
}

std::string caseName(InitialCase initial_case) {
    switch (initial_case) {
        case InitialCase::ShockBubble:
            return "shock-bubble";
        case InitialCase::PlanarShock:
            return "planar-shock";
    }

    throw std::runtime_error("Unsupported initial case.");
}

InitialCase parseCaseName(const std::string& value) {
    if (value == "shock-bubble") {
        return InitialCase::ShockBubble;
    }
    if (value == "planar-shock") {
        return InitialCase::PlanarShock;
    }

    throw std::runtime_error(
        "Unknown case '" + value + "'. Expected 'shock-bubble' or 'planar-shock'."
    );
}

std::string variantSummary() {
    std::ostringstream summary;
    summary << "state-storage=" << kStateStorageLabel
            << ", state-passing=" << kStatePassingLabel
            << ", grid-passing=" << kGridPassingLabel
            << ", access-pattern=" << kAccessPatternLabel
            << ", reconstruction-cache=" << kReconstructionCacheLabel;
    return summary.str();
}

void printUsage(const char* program_name) {
    std::cout
        << "Usage: " << program_name << " [options]\n"
        << "Options:\n"
        << "  --nx <int>                 Number of active cells in x (default: 500)\n"
        << "  --ny <int>                 Number of active cells in y (default: 197)\n"
        << "  --final-time <double>      Final simulation time in seconds (default: 3.0e-4)\n"
        << "  --cfl <double>             CFL number (default: 0.45)\n"
        << "  --snapshot-interval <double>  Write intermediate CSV snapshots every given time interval\n"
        << "  --output-prefix <string>   Prefix for CSV output files (default: output/shock_bubble)\n"
        << "  --case <name>              Initial condition: shock-bubble or planar-shock\n"
        << "  --help                     Show this message\n";
}

Settings parseCommandLine(int argc, char** argv) {
    Settings settings;

    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];

        if (argument == "--help") {
            printUsage(argv[0]);
            std::exit(0);
        }
        if (argument == "--nx") {
            settings.nx = std::stoi(nextArgument(index, argc, argv));
            continue;
        }
        if (argument == "--ny") {
            settings.ny = std::stoi(nextArgument(index, argc, argv));
            continue;
        }
        if (argument == "--final-time") {
            settings.final_time = std::stod(nextArgument(index, argc, argv));
            continue;
        }
        if (argument == "--cfl") {
            settings.cfl = std::stod(nextArgument(index, argc, argv));
            continue;
        }
        if (argument == "--snapshot-interval") {
            settings.snapshot_interval = std::stod(nextArgument(index, argc, argv));
            continue;
        }
        if (argument == "--output-prefix") {
            settings.output_prefix = nextArgument(index, argc, argv);
            continue;
        }
        if (argument == "--case") {
            settings.initial_case = parseCaseName(nextArgument(index, argc, argv));
            continue;
        }

        throw std::runtime_error("Unknown argument: " + argument);
    }

    if (settings.nx < 4 || settings.ny < 4) {
        throw std::runtime_error("Both nx and ny must be at least 4.");
    }
    if (settings.final_time <= 0.0) {
        throw std::runtime_error("Final time must be positive.");
    }
    if (settings.cfl <= 0.0 || settings.cfl >= 1.0) {
        throw std::runtime_error("CFL should lie between 0 and 1.");
    }
    if (settings.snapshot_interval < 0.0) {
        throw std::runtime_error("Snapshot interval cannot be negative.");
    }

    return settings;
}

class EulerSolver2D {
public:
    explicit EulerSolver2D(Settings settings)
        : settings_(std::move(settings)),
          nx_total_(settings_.nx + 2 * ghost_cells_),
          ny_total_(settings_.ny + 2 * ghost_cells_),
          y_min_(-0.5 * domain_width_),
          y_max_(0.5 * domain_width_),
          dx_((x_max_ - x_min_) / static_cast<double>(settings_.nx)),
          dy_((y_max_ - y_min_) / static_cast<double>(settings_.ny)),
          density_floor_(1.0e-10),
          pressure_floor_(ambient_pressure_ * 1.0e-8),
          x_full_access_order_(buildIndexOrder(0, nx_total_)),
          y_full_access_order_(buildIndexOrder(0, ny_total_)),
          x_active_access_order_(buildIndexOrder(ghost_cells_, settings_.nx)),
          y_active_access_order_(buildIndexOrder(ghost_cells_, settings_.ny)),
          cells_(static_cast<std::size_t>(nx_total_ * ny_total_)) {
        initialiseCase();
        applyTransmissiveBoundaries();
    }

    void run() {
        createOutputDirectory();

        const Primitive post_shock = postShockState();
        std::cout << std::setprecision(8)
                  << "Running Euler case with " << settings_.nx << " x " << settings_.ny
                  << " active cells.\n"
                  << "Initial condition: " << caseName(settings_.initial_case) << '\n'
                  << "Variant: " << variantSummary() << '\n'
                  << "Final time: " << settings_.final_time << " s, CFL: " << settings_.cfl << '\n'
                  << "Post-shock state: rho = " << post_shock.density()
                  << ", u = " << post_shock.velocity_x()
                  << ", p = " << post_shock.pressure() << '\n';

        double time = 0.0;
        int step = 0;
        int snapshot_index = 0;
        double last_snapshot_time = -1.0;
        const double epsilon = 1.0e-14;
        double next_snapshot_time = settings_.snapshot_interval > 0.0
            ? settings_.snapshot_interval
            : std::numeric_limits<double>::infinity();

        writeSnapshot(snapshotFileName(settings_.output_prefix, snapshot_index++), time, step);
        last_snapshot_time = time;

        while (time + epsilon < settings_.final_time) {
            applyTransmissiveBoundaries();
            double dt = computeStableTimeStep();
            if (time + dt > settings_.final_time) {
                dt = settings_.final_time - time;
            }

            const bool x_first = (step % 2 == 0);
            takeStrangStep(dt, x_first);

            time += dt;
            ++step;

            if (step == 1 || step % 25 == 0 || time + epsilon >= settings_.final_time) {
                std::cout << "Step " << step
                          << ", time = " << time
                          << ", dt = " << dt << '\n';
            }

            if (time + epsilon >= next_snapshot_time && time + epsilon < settings_.final_time) {
                writeSnapshot(snapshotFileName(settings_.output_prefix, snapshot_index++), time, step);
                last_snapshot_time = time;
                next_snapshot_time += settings_.snapshot_interval;
            }
        }

        if (std::abs(last_snapshot_time - time) > epsilon) {
            writeSnapshot(snapshotFileName(settings_.output_prefix, snapshot_index), time, step);
        }

        std::cout << "Finished. Snapshots written with prefix '" << settings_.output_prefix << "'.\n";
    }

private:
    using Grid = std::vector<Conserved>;

    static constexpr int ghost_cells_ = 2;
    static constexpr double gamma_gas_ = 1.4;
    static constexpr double x_min_ = 0.0;
    static constexpr double x_max_ = 0.225;
    static constexpr double domain_width_ = 0.089;
    static constexpr double bubble_center_x_ = 0.035;
    static constexpr double bubble_center_y_ = 0.0;
    static constexpr double bubble_radius_ = 0.025;
    static constexpr double shock_location_ = 0.005;
    static constexpr double ambient_pressure_ = 1.01325e5;
    static constexpr double air_density_ = 1.29;
    static constexpr double helium_density_ = 0.214;
    static constexpr double shock_mach_number_ = 1.22;

    Settings settings_;
    int nx_total_ = 0;
    int ny_total_ = 0;
    double y_min_ = 0.0;
    double y_max_ = 0.0;
    double dx_ = 0.0;
    double dy_ = 0.0;
    double density_floor_ = 0.0;
    double pressure_floor_ = 0.0;
    std::vector<int> x_full_access_order_;
    std::vector<int> y_full_access_order_;
    std::vector<int> x_active_access_order_;
    std::vector<int> y_active_access_order_;
    Grid cells_;

    int index(int i, int j) const {
        return j * nx_total_ + i;
    }

    double cellCenterX(int i) const {
        return x_min_ + (static_cast<double>(i - ghost_cells_) + 0.5) * dx_;
    }

    double cellCenterY(int j) const {
        return y_min_ + (static_cast<double>(j - ghost_cells_) + 0.5) * dy_;
    }

    void createOutputDirectory() const {
        const std::filesystem::path prefix_path(settings_.output_prefix);
        const std::filesystem::path directory = prefix_path.parent_path();
        if (!directory.empty()) {
            std::filesystem::create_directories(directory);
        }
    }

    void clampPrimitive(Primitive& state) const {
        state.density() = std::max(state.density(), density_floor_);
        state.pressure() = std::max(state.pressure(), pressure_floor_);
    }

    void clampConserved(Conserved& state) const {
        state.density() = std::max(state.density(), density_floor_);
        const double inverse_density = 1.0 / state.density();
        const double velocity_x = state.momentum_x() * inverse_density;
        const double velocity_y = state.momentum_y() * inverse_density;
        const double kinetic_energy = 0.5 * state.density()
            * (velocity_x * velocity_x + velocity_y * velocity_y);
        const double minimum_energy = pressure_floor_ / (gamma_gas_ - 1.0) + kinetic_energy;
        if (state.energy() < minimum_energy) {
            state.energy() = minimum_energy;
        }
    }

    Primitive conservedToPrimitive(Conserved state) const {
        clampConserved(state);
        const double inverse_density = 1.0 / state.density();
        const double velocity_x = state.momentum_x() * inverse_density;
        const double velocity_y = state.momentum_y() * inverse_density;
        const double kinetic_energy = 0.5 * state.density()
            * (velocity_x * velocity_x + velocity_y * velocity_y);
        const double pressure = (gamma_gas_ - 1.0) * (state.energy() - kinetic_energy);
        return {state.density(), velocity_x, velocity_y, std::max(pressure, pressure_floor_)};
    }

    Conserved primitiveToConserved(Primitive state) const {
        clampPrimitive(state);
        const double kinetic_energy = 0.5 * state.density()
            * (state.velocity_x() * state.velocity_x() + state.velocity_y() * state.velocity_y());
        return {
            state.density(),
            state.density() * state.velocity_x(),
            state.density() * state.velocity_y(),
            state.pressure() / (gamma_gas_ - 1.0) + kinetic_energy
        };
    }

    Conserved physicalFlux(ConservedArgument state, Direction direction) const {
        const Primitive primitive = conservedToPrimitive(state);

        if (direction == Direction::X) {
            return {
                state.momentum_x(),
                state.momentum_x() * primitive.velocity_x() + primitive.pressure(),
                state.momentum_x() * primitive.velocity_y(),
                primitive.velocity_x() * (state.energy() + primitive.pressure())
            };
        }

        return {
            state.momentum_y(),
            state.momentum_y() * primitive.velocity_x(),
            state.momentum_y() * primitive.velocity_y() + primitive.pressure(),
            primitive.velocity_y() * (state.energy() + primitive.pressure())
        };
    }

    Conserved forceFlux(
        ConservedArgument left_state,
        ConservedArgument right_state,
        double cell_width,
        double dt,
        Direction direction
    ) const {
        const Conserved left_flux = physicalFlux(left_state, direction);
        const Conserved right_flux = physicalFlux(right_state, direction);

        const Conserved lax_friedrichs_flux =
            0.5 * (left_flux + right_flux) - 0.5 * (cell_width / dt) * (right_state - left_state);

        Conserved richtmyer_state =
            0.5 * (left_state + right_state) - 0.5 * (dt / cell_width) * (right_flux - left_flux);
        clampConserved(richtmyer_state);
        const Conserved richtmyer_flux = physicalFlux(richtmyer_state, direction);

        return 0.5 * (lax_friedrichs_flux + richtmyer_flux);
    }

    Primitive limitedSlope(
        PrimitiveArgument left,
        PrimitiveArgument centre,
        PrimitiveArgument right
    ) const {
        return {
            minmod(centre.density() - left.density(), right.density() - centre.density()),
            minmod(
                centre.velocity_x() - left.velocity_x(),
                right.velocity_x() - centre.velocity_x()
            ),
            minmod(
                centre.velocity_y() - left.velocity_y(),
                right.velocity_y() - centre.velocity_y()
            ),
            minmod(centre.pressure() - left.pressure(), right.pressure() - centre.pressure())
        };
    }

    void reconstructHalfStepStates(
        PrimitiveArgument left_cell,
        PrimitiveArgument centre_cell,
        PrimitiveArgument right_cell,
        double cell_width,
        double dt,
        Direction direction,
        Conserved& left_half_step,
        Conserved& right_half_step
    ) const {
        const Primitive slope = limitedSlope(left_cell, centre_cell, right_cell);

        Primitive left_state = centre_cell - 0.5 * slope;
        Primitive right_state = centre_cell + 0.5 * slope;
        clampPrimitive(left_state);
        clampPrimitive(right_state);

        const Conserved left_conserved = primitiveToConserved(left_state);
        const Conserved right_conserved = primitiveToConserved(right_state);
        const Conserved left_flux = physicalFlux(left_conserved, direction);
        const Conserved right_flux = physicalFlux(right_conserved, direction);

        left_half_step = left_conserved - 0.5 * (dt / cell_width) * (right_flux - left_flux);
        right_half_step = right_conserved - 0.5 * (dt / cell_width) * (right_flux - left_flux);

        clampConserved(left_half_step);
        clampConserved(right_half_step);
    }

    Primitive postShockState() const {
        const double sound_speed = std::sqrt(gamma_gas_ * ambient_pressure_ / air_density_);
        const double shock_speed = shock_mach_number_ * sound_speed;
        const double mach_squared = shock_mach_number_ * shock_mach_number_;

        const double density_ratio =
            ((gamma_gas_ + 1.0) * mach_squared) / ((gamma_gas_ - 1.0) * mach_squared + 2.0);
        const double pressure_ratio =
            1.0 + (2.0 * gamma_gas_ / (gamma_gas_ + 1.0)) * (mach_squared - 1.0);

        const double post_shock_density = air_density_ * density_ratio;
        const double post_shock_pressure = ambient_pressure_ * pressure_ratio;
        const double post_shock_velocity = shock_speed * (1.0 - 1.0 / density_ratio);

        return {post_shock_density, post_shock_velocity, 0.0, post_shock_pressure};
    }

    void initialiseShockBubble() {
        const Primitive shocked_air = postShockState();

        for (int j : y_active_access_order_) {
            for (int i : x_active_access_order_) {
                const double x = cellCenterX(i);
                const double y = cellCenterY(j);

                Primitive state;
                if (x < shock_location_) {
                    state = shocked_air;
                } else {
                    const double distance_to_bubble =
                        std::hypot(x - bubble_center_x_, y - bubble_center_y_);
                    const bool inside_bubble = distance_to_bubble <= bubble_radius_;
                    state = {
                        inside_bubble ? helium_density_ : air_density_,
                        0.0,
                        0.0,
                        ambient_pressure_
                    };
                }

                cells_[index(i, j)] = primitiveToConserved(state);
            }
        }
    }

    void initialisePlanarShock() {
        const Primitive shocked_air = postShockState();

        for (int j : y_active_access_order_) {
            for (int i : x_active_access_order_) {
                const double x = cellCenterX(i);
                const Primitive state = (x < shock_location_)
                    ? shocked_air
                    : Primitive{air_density_, 0.0, 0.0, ambient_pressure_};
                cells_[index(i, j)] = primitiveToConserved(state);
            }
        }
    }

    void initialiseCase() {
        if (settings_.initial_case == InitialCase::ShockBubble) {
            initialiseShockBubble();
            return;
        }

        initialisePlanarShock();
    }

    void applyTransmissiveBoundaries() {
        for (int j = ghost_cells_; j < ny_total_ - ghost_cells_; ++j) {
            for (int ghost = 0; ghost < ghost_cells_; ++ghost) {
                cells_[index(ghost, j)] = cells_[index(ghost_cells_, j)];
                cells_[index(nx_total_ - 1 - ghost, j)] =
                    cells_[index(nx_total_ - 1 - ghost_cells_, j)];
            }
        }

        for (int i = 0; i < nx_total_; ++i) {
            for (int ghost = 0; ghost < ghost_cells_; ++ghost) {
                cells_[index(i, ghost)] = cells_[index(i, ghost_cells_)];
                cells_[index(i, ny_total_ - 1 - ghost)] =
                    cells_[index(i, ny_total_ - 1 - ghost_cells_)];
            }
        }
    }

    double computeStableTimeStep() const {
        double max_x_speed = 0.0;
        double max_y_speed = 0.0;

        for (int j : y_active_access_order_) {
            for (int i : x_active_access_order_) {
                const Primitive state = conservedToPrimitive(cells_[index(i, j)]);
                const double sound_speed = std::sqrt(gamma_gas_ * state.pressure() / state.density());
                max_x_speed = std::max(max_x_speed, std::abs(state.velocity_x()) + sound_speed);
                max_y_speed = std::max(max_y_speed, std::abs(state.velocity_y()) + sound_speed);
            }
        }

        const double denominator = max_x_speed / dx_ + max_y_speed / dy_;
        if (denominator <= 0.0) {
            throw std::runtime_error("Computed a non-positive stable time step.");
        }

        return settings_.cfl / denominator;
    }

    void advanceLine(
        std::vector<Conserved>& line,
        double cell_width,
        double dt,
        Direction direction,
        LineWorkspace& workspace
    ) const {
        const int length = static_cast<int>(line.size());

        for (int cell = 0; cell < length; ++cell) {
            workspace.primitive_values[static_cast<std::size_t>(cell)] =
                conservedToPrimitive(line[static_cast<std::size_t>(cell)]);
        }

        if constexpr (kUseReconstructionCache) {
            for (int cell = 1; cell < length - 1; ++cell) {
                reconstructHalfStepStates(
                    workspace.primitive_values[static_cast<std::size_t>(cell - 1)],
                    workspace.primitive_values[static_cast<std::size_t>(cell)],
                    workspace.primitive_values[static_cast<std::size_t>(cell + 1)],
                    cell_width,
                    dt,
                    direction,
                    workspace.left_half_step[static_cast<std::size_t>(cell)],
                    workspace.right_half_step[static_cast<std::size_t>(cell)]
                );
            }

            for (int interface_index = ghost_cells_; interface_index <= length - ghost_cells_; ++interface_index) {
                const Conserved& left_state =
                    workspace.right_half_step[static_cast<std::size_t>(interface_index - 1)];
                const Conserved& right_state =
                    workspace.left_half_step[static_cast<std::size_t>(interface_index)];
                workspace.interface_fluxes[static_cast<std::size_t>(interface_index)] =
                    forceFlux(left_state, right_state, cell_width, dt, direction);
            }
        } else {
            for (int interface_index = ghost_cells_; interface_index <= length - ghost_cells_; ++interface_index) {
                Conserved left_left_half_step;
                Conserved left_right_half_step;
                Conserved right_left_half_step;
                Conserved right_right_half_step;

                reconstructHalfStepStates(
                    workspace.primitive_values[static_cast<std::size_t>(interface_index - 2)],
                    workspace.primitive_values[static_cast<std::size_t>(interface_index - 1)],
                    workspace.primitive_values[static_cast<std::size_t>(interface_index)],
                    cell_width,
                    dt,
                    direction,
                    left_left_half_step,
                    left_right_half_step
                );

                reconstructHalfStepStates(
                    workspace.primitive_values[static_cast<std::size_t>(interface_index - 1)],
                    workspace.primitive_values[static_cast<std::size_t>(interface_index)],
                    workspace.primitive_values[static_cast<std::size_t>(interface_index + 1)],
                    cell_width,
                    dt,
                    direction,
                    right_left_half_step,
                    right_right_half_step
                );

                workspace.interface_fluxes[static_cast<std::size_t>(interface_index)] =
                    forceFlux(left_right_half_step, right_left_half_step, cell_width, dt, direction);
            }
        }

        for (int cell = ghost_cells_; cell < length - ghost_cells_; ++cell) {
            line[static_cast<std::size_t>(cell)] = line[static_cast<std::size_t>(cell)] - (dt / cell_width)
                * (
                    workspace.interface_fluxes[static_cast<std::size_t>(cell + 1)]
                    - workspace.interface_fluxes[static_cast<std::size_t>(cell)]
                );
            clampConserved(line[static_cast<std::size_t>(cell)]);
        }
    }

    void sweepXInPlace(Grid& grid, double dt) const {
        std::vector<Conserved> line(static_cast<std::size_t>(nx_total_));
        LineWorkspace workspace;
        workspace.resize(line.size());

        for (int j : y_active_access_order_) {
            for (int i : x_full_access_order_) {
                line[static_cast<std::size_t>(i)] = grid[static_cast<std::size_t>(index(i, j))];
            }

            advanceLine(line, dx_, dt, Direction::X, workspace);

            for (int i : x_active_access_order_) {
                grid[static_cast<std::size_t>(index(i, j))] = line[static_cast<std::size_t>(i)];
            }
        }
    }

    void sweepYInPlace(Grid& grid, double dt) const {
        std::vector<Conserved> line(static_cast<std::size_t>(ny_total_));
        LineWorkspace workspace;
        workspace.resize(line.size());

        for (int i : x_active_access_order_) {
            for (int j : y_full_access_order_) {
                line[static_cast<std::size_t>(j)] = grid[static_cast<std::size_t>(index(i, j))];
            }

            advanceLine(line, dy_, dt, Direction::Y, workspace);

            for (int j : y_active_access_order_) {
                grid[static_cast<std::size_t>(index(i, j))] = line[static_cast<std::size_t>(j)];
            }
        }
    }

    Grid sweepXByValue(Grid grid, double dt) const {
        sweepXInPlace(grid, dt);
        return grid;
    }

    Grid sweepYByValue(Grid grid, double dt) const {
        sweepYInPlace(grid, dt);
        return grid;
    }

    void takeStrangStep(double dt, bool x_first) {
        if constexpr (kPassGridByValue) {
            if (x_first) {
                cells_ = sweepXByValue(cells_, 0.5 * dt);
                applyTransmissiveBoundaries();
                cells_ = sweepYByValue(cells_, dt);
                applyTransmissiveBoundaries();
                cells_ = sweepXByValue(cells_, 0.5 * dt);
            } else {
                cells_ = sweepYByValue(cells_, 0.5 * dt);
                applyTransmissiveBoundaries();
                cells_ = sweepXByValue(cells_, dt);
                applyTransmissiveBoundaries();
                cells_ = sweepYByValue(cells_, 0.5 * dt);
            }
        } else {
            if (x_first) {
                sweepXInPlace(cells_, 0.5 * dt);
                applyTransmissiveBoundaries();
                sweepYInPlace(cells_, dt);
                applyTransmissiveBoundaries();
                sweepXInPlace(cells_, 0.5 * dt);
            } else {
                sweepYInPlace(cells_, 0.5 * dt);
                applyTransmissiveBoundaries();
                sweepXInPlace(cells_, dt);
                applyTransmissiveBoundaries();
                sweepYInPlace(cells_, 0.5 * dt);
            }
        }

        applyTransmissiveBoundaries();
    }

    void writeSnapshot(const std::string& file_name, double time, int step) const {
        std::ofstream output(file_name);
        if (!output) {
            throw std::runtime_error("Failed to open output file: " + file_name);
        }

        output << std::setprecision(12);
        output << "# case=" << caseName(settings_.initial_case) << '\n';
        output << "# variant=" << variantSummary() << '\n';
        output << "# state_storage=" << kStateStorageLabel << '\n';
        output << "# state_passing=" << kStatePassingLabel << '\n';
        output << "# grid_passing=" << kGridPassingLabel << '\n';
        output << "# access_pattern=" << kAccessPatternLabel << '\n';
        output << "# reconstruction_cache=" << kReconstructionCacheLabel << '\n';
        output << "# nx=" << settings_.nx << '\n';
        output << "# ny=" << settings_.ny << '\n';
        output << "# dx=" << dx_ << '\n';
        output << "# dy=" << dy_ << '\n';
        output << "# time=" << time << '\n';
        output << "# step=" << step << '\n';
        output << "x,y,density,pressure,velocity_x,velocity_y\n";

        for (int j = ghost_cells_; j < ny_total_ - ghost_cells_; ++j) {
            for (int i = ghost_cells_; i < nx_total_ - ghost_cells_; ++i) {
                const Primitive state = conservedToPrimitive(cells_[static_cast<std::size_t>(index(i, j))]);
                output << cellCenterX(i) << ','
                       << cellCenterY(j) << ','
                       << state.density() << ','
                       << state.pressure() << ','
                       << state.velocity_x() << ','
                       << state.velocity_y() << '\n';
            }
        }
    }
};

}  // namespace

int main(int argc, char** argv) {
    try {
        const Settings settings = parseCommandLine(argc, argv);
        EulerSolver2D solver(settings);
        solver.run();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
