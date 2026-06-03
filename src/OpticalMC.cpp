#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <random>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct RunConfig {
    std::string ratio_tag;
    double thickness_um = 0.0;
    std::string readout_surface = "back";
    int samples_per_step = 16;
    std::string xy_boundary = "infinite";
    uint64_t random_seed = 12345;
    int num_threads = 1;
    int max_steps = 10000;
    double roulette_threshold = 0.0;
    double roulette_survival_probability = 1.0;
    std::string front_reflection_model = "aluminum_fresnel";
    double front_reflectance = 0.0;
    std::string front_reflection_mode = "specular";
    double front_aluminum_n = 0.65;
    double front_aluminum_k = 5.3;
    std::string back_reflection_model = "air_fresnel";
    double back_air_n = 1.000293;
    bool output_detected_photons = false;
    double psf_bin_size_um = 5.0;
    double psf_range_um = 500.0;
    std::string optical_properties_csv = "optical_properties.csv";
    std::string source_steps_csv;
    std::string event_sources_csv;
    std::string output_dir = ".";
    double wavelength_nm = std::numeric_limits<double>::quiet_NaN();
    double incident_event_count = 0.0;
};

struct OpticalProperties {
    std::string ratio_tag;
    double wavelength_nm = 0.0;
    double mu_a_per_um = 0.0;
    double mu_s_per_um = 0.0;
    double g = 0.0;
    double n_eff = 1.0;
    double mu_s_prime_per_um = std::numeric_limits<double>::quiet_NaN();
    std::string phase_function_csv;
    std::vector<double> phase_mu_min;
    std::vector<double> phase_mu_max;
    std::vector<double> phase_cdf;
};

struct SourceStep {
    int source_id = 0;
    std::string source_event_uid;
    std::string eventID;
    std::string ratio_tag;
    double wavelength_nm = 0.0;
    Vec3 p0;
    Vec3 p1;
    double n_photon_step = 0.0;
    bool has_psf_anchor = false;
    double psf_anchor_x = 0.0;
    double psf_anchor_y = 0.0;
};

struct EventSource {
    std::string source_event_uid;
    std::string eventID;
    double depth_um = 0.0;
    double total_n_photon = 0.0;
    bool has_psf_anchor = false;
    double psf_anchor_x = 0.0;
    double psf_anchor_y = 0.0;
};

struct PhotonRecord {
    std::string source_event_uid;
    int source_id = 0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double ux = 0.0;
    double uy = 0.0;
    double uz = 0.0;
    double weight = 0.0;
    double path_length_um = 0.0;
    int scatter_count = 0;
    std::string surface;
};

struct StepStats {
    double detected_weight = 0.0;
    double absorbed_weight = 0.0;
    double front_escape_weight = 0.0;
    double back_escape_weight = 0.0;
    double lost_weight = 0.0;
    double weighted_path_length = 0.0;
    double weighted_scatter_count = 0.0;
    double launched_weight = 0.0;
    int samples = 0;
};

struct EventStats {
    double detected_weight = 0.0;
    double sum_x = 0.0;
    double sum_y = 0.0;
    double sum_x2 = 0.0;
    double sum_y2 = 0.0;
};

struct Accumulators {
    std::vector<StepStats> step_stats;
    std::unordered_map<std::string, EventStats> event_stats;
    std::vector<double> psf_2d;
    std::vector<double> lsf_x;
    std::vector<double> lsf_y;
    std::vector<PhotonRecord> detected;
    double detected_weight = 0.0;
    double front_escape_weight = 0.0;
    double back_escape_weight = 0.0;
    double absorbed_weight = 0.0;
    double lost_weight = 0.0;
    double total_source_weight = 0.0;
    double detected_sum_x = 0.0;
    double detected_sum_y = 0.0;
    double detected_sum_x2 = 0.0;
    double detected_sum_y2 = 0.0;
};

std::string trim(const std::string& s) {
    const auto first = s.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return "";
    }
    const auto last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, last - first + 1);
}

std::string unquote(std::string s) {
    s = trim(s);
    if (s.size() >= 2 && s.front() == '"' && s.back() == '"') {
        s = s.substr(1, s.size() - 2);
    }
    std::string out;
    out.reserve(s.size());
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '\\' && i + 1 < s.size()) {
            ++i;
            switch (s[i]) {
                case 'n': out.push_back('\n'); break;
                case 'r': out.push_back('\r'); break;
                case 't': out.push_back('\t'); break;
                default: out.push_back(s[i]); break;
            }
        } else {
            out.push_back(s[i]);
        }
    }
    return out;
}

std::vector<std::string> splitCsvLine(const std::string& line) {
    std::vector<std::string> fields;
    std::string current;
    bool in_quotes = false;
    for (size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (c == '"') {
            if (in_quotes && i + 1 < line.size() && line[i + 1] == '"') {
                current.push_back('"');
                ++i;
            } else {
                in_quotes = !in_quotes;
            }
        } else if (c == ',' && !in_quotes) {
            fields.push_back(current);
            current.clear();
        } else {
            current.push_back(c);
        }
    }
    fields.push_back(current);
    return fields;
}

std::string csvEscape(const std::string& s) {
    if (s.find_first_of(",\"\r\n") == std::string::npos) {
        return s;
    }
    std::string out = "\"";
    for (char c : s) {
        if (c == '"') {
            out += "\"\"";
        } else {
            out.push_back(c);
        }
    }
    out.push_back('"');
    return out;
}

std::string num(double value) {
    if (!std::isfinite(value)) {
        return "";
    }
    std::ostringstream os;
    os << std::setprecision(12) << value;
    return os.str();
}

double toDouble(const std::unordered_map<std::string, std::string>& row, const std::string& key,
                double default_value = 0.0) {
    const auto it = row.find(key);
    if (it == row.end() || trim(it->second).empty()) {
        return default_value;
    }
    return std::stod(it->second);
}

bool optionalDouble(const std::unordered_map<std::string, std::string>& row,
                    const std::string& key, double& value) {
    const auto it = row.find(key);
    if (it == row.end() || trim(it->second).empty()) {
        return false;
    }
    try {
        value = std::stod(it->second);
    } catch (...) {
        return false;
    }
    return std::isfinite(value);
}

int toInt(const std::unordered_map<std::string, std::string>& row, const std::string& key,
          int default_value = 0) {
    const auto it = row.find(key);
    if (it == row.end() || trim(it->second).empty()) {
        return default_value;
    }
    return std::stoi(it->second);
}

std::string toString(const std::unordered_map<std::string, std::string>& row,
                     const std::string& key) {
    const auto it = row.find(key);
    return it == row.end() ? "" : it->second;
}

std::vector<std::unordered_map<std::string, std::string>> readCsvRows(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("cannot open CSV: " + path);
    }
    std::string header_line;
    if (!std::getline(in, header_line)) {
        throw std::runtime_error("empty CSV: " + path);
    }
    if (!header_line.empty() && static_cast<unsigned char>(header_line[0]) == 0xEF) {
        header_line.erase(0, 3);
    }
    auto headers = splitCsvLine(header_line);
    std::vector<std::unordered_map<std::string, std::string>> rows;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        auto values = splitCsvLine(line);
        std::unordered_map<std::string, std::string> row;
        for (size_t i = 0; i < headers.size(); ++i) {
            row[headers[i]] = i < values.size() ? values[i] : "";
        }
        rows.push_back(std::move(row));
    }
    return rows;
}

bool pathExists(const std::string& path) {
    std::ifstream in(path);
    return static_cast<bool>(in);
}

bool isAbsolutePath(const std::string& path) {
    if (path.empty()) {
        return false;
    }
    if (path.front() == '/' || path.front() == '\\') {
        return true;
    }
    return path.size() >= 3 && std::isalpha(static_cast<unsigned char>(path[0])) &&
           path[1] == ':' && (path[2] == '\\' || path[2] == '/');
}

std::string directoryName(const std::string& path) {
    const auto pos = path.find_last_of("/\\");
    if (pos == std::string::npos) {
        return "";
    }
    return path.substr(0, pos);
}

std::string joinPath(const std::string& dir, const std::string& file) {
    if (dir.empty() || file.empty() || isAbsolutePath(file)) {
        return file;
    }
    const char last = dir.back();
    if (last == '/' || last == '\\') {
        return dir + file;
    }
    return dir + "/" + file;
}

std::string resolveDataPath(const std::string& candidate, const std::string& base_file) {
    if (candidate.empty() || isAbsolutePath(candidate) || pathExists(candidate)) {
        return candidate;
    }
    const auto sibling = joinPath(directoryName(base_file), candidate);
    return pathExists(sibling) ? sibling : candidate;
}

void loadPhaseFunction(OpticalProperties& op, const std::string& optical_properties_path) {
    const auto resolved = resolveDataPath(trim(op.phase_function_csv), optical_properties_path);
    double cumulative = 0.0;
    op.phase_mu_min.clear();
    op.phase_mu_max.clear();
    op.phase_cdf.clear();

    for (const auto& row : readCsvRows(resolved)) {
        double probability = toDouble(row, "probability", -1.0);
        if (probability < 0.0) {
            probability = toDouble(row, "probability_mean", 0.0);
        }
        if (probability <= 0.0) {
            continue;
        }
        const double mu_min = toDouble(row, "cos_theta_min", 0.0);
        const double mu_max = toDouble(row, "cos_theta_max", 0.0);
        if (!std::isfinite(mu_min) || !std::isfinite(mu_max) || mu_max < mu_min) {
            throw std::runtime_error("invalid cos(theta) bin in phase function: " + resolved);
        }
        cumulative += probability;
        op.phase_mu_min.push_back(std::max(-1.0, std::min(1.0, mu_min)));
        op.phase_mu_max.push_back(std::max(-1.0, std::min(1.0, mu_max)));
        op.phase_cdf.push_back(cumulative);
    }
    if (op.phase_cdf.empty() || cumulative <= 0.0) {
        throw std::runtime_error("phase function has no positive-probability bins: " + resolved);
    }
    for (double& value : op.phase_cdf) {
        value /= cumulative;
    }
    op.phase_cdf.back() = 1.0;
    op.phase_function_csv = resolved;
}

std::string readText(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("cannot open file: " + path);
    }
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

bool jsonString(const std::string& text, const std::string& key, std::string& value) {
    const std::regex re("\"" + key + "\"\\s*:\\s*\"((?:\\\\.|[^\"])*)\"");
    std::smatch match;
    if (std::regex_search(text, match, re)) {
        value = unquote("\"" + match[1].str() + "\"");
        return true;
    }
    return false;
}

bool jsonNumber(const std::string& text, const std::string& key, double& value) {
    const std::regex re("\"" + key + "\"\\s*:\\s*(-?(?:\\d+\\.?\\d*|\\.\\d+)(?:[eE][+-]?\\d+)?)");
    std::smatch match;
    if (std::regex_search(text, match, re)) {
        value = std::stod(match[1].str());
        return true;
    }
    return false;
}

bool jsonBool(const std::string& text, const std::string& key, bool& value) {
    const std::regex re("\"" + key + "\"\\s*:\\s*(true|false|1|0)");
    std::smatch match;
    if (std::regex_search(text, match, re)) {
        const auto v = match[1].str();
        value = (v == "true" || v == "1");
        return true;
    }
    return false;
}

RunConfig readConfig(const std::string& path) {
    RunConfig cfg;
    const auto text = readText(path);
    double d = 0.0;
    std::string s;
    bool b = false;
    if (jsonString(text, "ratio_tag", s)) cfg.ratio_tag = s;
    if (jsonNumber(text, "thickness_um", d)) cfg.thickness_um = d;
    if (jsonString(text, "readout_surface", s)) cfg.readout_surface = s;
    if (jsonNumber(text, "samples_per_step", d)) cfg.samples_per_step = static_cast<int>(d);
    if (jsonString(text, "xy_boundary", s)) cfg.xy_boundary = s;
    if (jsonNumber(text, "random_seed", d)) cfg.random_seed = static_cast<uint64_t>(d);
    if (jsonNumber(text, "num_threads", d)) cfg.num_threads = static_cast<int>(d);
    if (jsonNumber(text, "max_steps", d)) cfg.max_steps = static_cast<int>(d);
    if (jsonNumber(text, "roulette_threshold", d)) cfg.roulette_threshold = d;
    if (jsonNumber(text, "roulette_survival_probability", d)) cfg.roulette_survival_probability = d;
    if (jsonString(text, "front_reflection_model", s)) cfg.front_reflection_model = s;
    if (jsonNumber(text, "front_reflectance", d)) cfg.front_reflectance = d;
    if (jsonString(text, "front_reflection_mode", s)) cfg.front_reflection_mode = s;
    if (jsonNumber(text, "front_aluminum_n", d)) cfg.front_aluminum_n = d;
    if (jsonNumber(text, "front_aluminum_k", d)) cfg.front_aluminum_k = d;
    if (jsonString(text, "back_reflection_model", s)) cfg.back_reflection_model = s;
    if (jsonNumber(text, "back_air_n", d)) cfg.back_air_n = d;
    if (jsonBool(text, "output_detected_photons", b)) cfg.output_detected_photons = b;
    if (jsonNumber(text, "psf_bin_size_um", d)) cfg.psf_bin_size_um = d;
    if (jsonNumber(text, "psf_range_um", d)) cfg.psf_range_um = d;
    if (jsonString(text, "optical_properties_csv", s)) cfg.optical_properties_csv = s;
    if (jsonString(text, "source_steps_csv", s)) cfg.source_steps_csv = s;
    if (jsonString(text, "event_sources_csv", s)) cfg.event_sources_csv = s;
    if (jsonString(text, "output_dir", s)) cfg.output_dir = s;
    if (jsonNumber(text, "wavelength_nm", d)) cfg.wavelength_nm = d;
    if (jsonNumber(text, "incident_event_count", d)) cfg.incident_event_count = d;

    if (cfg.thickness_um <= 0.0) {
        throw std::runtime_error("run_config.json must provide positive thickness_um");
    }
    if (cfg.samples_per_step <= 0) {
        throw std::runtime_error("samples_per_step must be positive");
    }
    if (cfg.max_steps <= 0) {
        throw std::runtime_error("max_steps must be positive");
    }
    if (cfg.num_threads <= 0) {
        cfg.num_threads = 1;
    }
    if (cfg.psf_bin_size_um <= 0.0 || cfg.psf_range_um <= 0.0) {
        throw std::runtime_error("psf_bin_size_um and psf_range_um must be positive");
    }
    if (cfg.front_reflectance < 0.0 || cfg.front_reflectance > 1.0) {
        throw std::runtime_error("front_reflectance must be between 0 and 1");
    }
    if (cfg.front_reflection_mode != "none" && cfg.front_reflection_mode != "specular" &&
        cfg.front_reflection_mode != "diffuse" && cfg.front_reflection_mode != "lambertian") {
        throw std::runtime_error(
            "front_reflection_mode must be none, specular, diffuse, or lambertian");
    }
    if (cfg.front_reflection_model != "effective" &&
        cfg.front_reflection_model != "aluminum_fresnel") {
        throw std::runtime_error("front_reflection_model must be effective or aluminum_fresnel");
    }
    if (cfg.front_aluminum_n < 0.0 || cfg.front_aluminum_k < 0.0) {
        throw std::runtime_error("front_aluminum_n and front_aluminum_k must be non-negative");
    }
    if (cfg.front_reflection_model == "aluminum_fresnel" &&
        cfg.front_aluminum_n == 0.0 && cfg.front_aluminum_k == 0.0) {
        throw std::runtime_error(
            "front_aluminum_n and front_aluminum_k cannot both be zero for aluminum_fresnel");
    }
    if (cfg.back_reflection_model != "none" && cfg.back_reflection_model != "air_fresnel") {
        throw std::runtime_error("back_reflection_model must be none or air_fresnel");
    }
    if (cfg.back_air_n <= 0.0) {
        throw std::runtime_error("back_air_n must be positive");
    }
    return cfg;
}

std::vector<SourceStep> readSources(const std::string& path) {
    std::vector<SourceStep> sources;
    for (const auto& row : readCsvRows(path)) {
        SourceStep s;
        s.source_id = toInt(row, "source_id", static_cast<int>(sources.size()));
        s.source_event_uid = toString(row, "source_event_uid");
        s.eventID = toString(row, "eventID");
        s.ratio_tag = toString(row, "ratio_tag");
        s.wavelength_nm = toDouble(row, "wavelength_nm", 0.0);
        s.p0 = {toDouble(row, "src_x0_um"), toDouble(row, "src_y0_um"),
                toDouble(row, "src_z0_um")};
        s.p1 = {toDouble(row, "src_x1_um"), toDouble(row, "src_y1_um"),
                toDouble(row, "src_z1_um")};
        s.n_photon_step = toDouble(row, "n_photon_step");
        double anchor_x = 0.0;
        double anchor_y = 0.0;
        if (optionalDouble(row, "macro_anchor_x_um", anchor_x) &&
            optionalDouble(row, "macro_anchor_y_um", anchor_y)) {
            s.has_psf_anchor = true;
            s.psf_anchor_x = anchor_x;
            s.psf_anchor_y = anchor_y;
        }
        if (s.n_photon_step > 0.0) {
            sources.push_back(std::move(s));
        }
    }
    return sources;
}

std::vector<EventSource> readEvents(const std::string& path) {
    std::vector<EventSource> events;
    for (const auto& row : readCsvRows(path)) {
        EventSource e;
        e.source_event_uid = toString(row, "source_event_uid");
        e.eventID = toString(row, "eventID");
        e.depth_um = toDouble(row, "depth_um", 0.0);
        e.total_n_photon = toDouble(row, "total_n_photon", 0.0);
        double anchor_x = 0.0;
        double anchor_y = 0.0;
        if (optionalDouble(row, "macro_anchor_x_um", anchor_x) &&
            optionalDouble(row, "macro_anchor_y_um", anchor_y)) {
            e.has_psf_anchor = true;
            e.psf_anchor_x = anchor_x;
            e.psf_anchor_y = anchor_y;
        }
        events.push_back(std::move(e));
    }
    return events;
}

void applyEventAnchors(std::vector<SourceStep>& sources,
                       const std::vector<EventSource>& events) {
    std::unordered_map<std::string, const EventSource*> by_uid;
    by_uid.reserve(events.size());
    for (const auto& event : events) {
        if (event.has_psf_anchor) {
            by_uid[event.source_event_uid] = &event;
        }
    }
    for (auto& source : sources) {
        if (source.has_psf_anchor) {
            continue;
        }
        const auto it = by_uid.find(source.source_event_uid);
        if (it == by_uid.end()) {
            continue;
        }
        source.has_psf_anchor = true;
        source.psf_anchor_x = it->second->psf_anchor_x;
        source.psf_anchor_y = it->second->psf_anchor_y;
    }
}

OpticalProperties selectOpticalProperties(const std::string& path, const RunConfig& cfg,
                                          const std::vector<SourceStep>& sources) {
    const double target_wavelength = std::isfinite(cfg.wavelength_nm)
                                         ? cfg.wavelength_nm
                                         : (sources.empty() ? std::numeric_limits<double>::quiet_NaN()
                                                            : sources.front().wavelength_nm);
    std::vector<OpticalProperties> candidates;
    for (const auto& row : readCsvRows(path)) {
        OpticalProperties op;
        op.ratio_tag = toString(row, "ratio_tag");
        op.wavelength_nm = toDouble(row, "wavelength_nm", 0.0);
        op.mu_a_per_um = toDouble(row, "mu_a_per_um", 0.0);
        op.mu_s_per_um = toDouble(row, "mu_s_per_um", 0.0);
        op.g = toDouble(row, "g", 0.0);
        op.n_eff = toDouble(row, "n_eff", 1.0);
        op.mu_s_prime_per_um = toDouble(row, "mu_s_prime_per_um",
                                        std::numeric_limits<double>::quiet_NaN());
        op.phase_function_csv = toString(row, "phase_function_csv");
        if (!cfg.ratio_tag.empty() && op.ratio_tag != cfg.ratio_tag) {
            continue;
        }
        candidates.push_back(op);
    }
    if (candidates.empty()) {
        throw std::runtime_error("no optical_properties.csv row matches ratio_tag=" + cfg.ratio_tag);
    }
    auto finalize = [&](OpticalProperties op) {
        if (!trim(op.phase_function_csv).empty()) {
            loadPhaseFunction(op, path);
        }
        return op;
    };
    if (!std::isfinite(target_wavelength)) {
        return finalize(candidates.front());
    }
    auto best = std::min_element(candidates.begin(), candidates.end(), [&](const auto& a, const auto& b) {
        return std::abs(a.wavelength_nm - target_wavelength) <
               std::abs(b.wavelength_nm - target_wavelength);
    });
    return finalize(*best);
}

double uniformOpen(std::mt19937_64& rng) {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    double x = 0.0;
    do {
        x = dist(rng);
    } while (x <= 0.0);
    return x;
}

Vec3 isotropicDirection(std::mt19937_64& rng) {
    const double mu = 2.0 * uniformOpen(rng) - 1.0;
    const double phi = 2.0 * kPi * uniformOpen(rng);
    const double sin_theta = std::sqrt(std::max(0.0, 1.0 - mu * mu));
    return {sin_theta * std::cos(phi), sin_theta * std::sin(phi), mu};
}

Vec3 scatterByCosTheta(const Vec3& u, double cos_theta, std::mt19937_64& rng) {
    cos_theta = std::min(1.0, std::max(-1.0, cos_theta));
    const double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
    const double phi = 2.0 * kPi * uniformOpen(rng);
    const double cos_phi = std::cos(phi);
    const double sin_phi = std::sin(phi);

    Vec3 out;
    if (std::abs(u.z) > 0.99999) {
        out.x = sin_theta * cos_phi;
        out.y = sin_theta * sin_phi;
        out.z = cos_theta * (u.z >= 0.0 ? 1.0 : -1.0);
    } else {
        const double denom = std::sqrt(1.0 - u.z * u.z);
        out.x = sin_theta * (u.x * u.z * cos_phi - u.y * sin_phi) / denom + u.x * cos_theta;
        out.y = sin_theta * (u.y * u.z * cos_phi + u.x * sin_phi) / denom + u.y * cos_theta;
        out.z = -sin_theta * cos_phi * denom + u.z * cos_theta;
    }
    const double norm = std::sqrt(out.x * out.x + out.y * out.y + out.z * out.z);
    return {out.x / norm, out.y / norm, out.z / norm};
}

double sampleHGCostheta(double g, std::mt19937_64& rng) {
    const double xi = uniformOpen(rng);
    if (std::abs(g) < 1.0e-12) {
        return 2.0 * xi - 1.0;
    }
    const double term = (1.0 - g * g) / (1.0 - g + 2.0 * g * xi);
    return std::min(1.0, std::max(-1.0, (1.0 + g * g - term * term) / (2.0 * g)));
}

double sampleTabulatedCostheta(const OpticalProperties& op, std::mt19937_64& rng) {
    const double xi = uniformOpen(rng);
    const auto it = std::lower_bound(op.phase_cdf.begin(), op.phase_cdf.end(), xi);
    const size_t index = static_cast<size_t>(
        std::min<std::ptrdiff_t>(std::distance(op.phase_cdf.begin(), it),
                                 static_cast<std::ptrdiff_t>(op.phase_cdf.size() - 1)));
    const double a = op.phase_mu_min[index];
    const double b = op.phase_mu_max[index];
    return a + (b - a) * uniformOpen(rng);
}

Vec3 scatterPhoton(const Vec3& u, const OpticalProperties& op, std::mt19937_64& rng) {
    const double cos_theta = op.phase_cdf.empty() ? sampleHGCostheta(op.g, rng)
                                                  : sampleTabulatedCostheta(op, rng);
    return scatterByCosTheta(u, cos_theta, rng);
}

bool isReadoutSurface(const std::string& readout, const std::string& surface) {
    return readout == "both" || readout == surface;
}

bool hasFrontReflection(const RunConfig& cfg) {
    return cfg.front_reflection_mode != "none";
}

bool hasBackReflection(const RunConfig& cfg) {
    return cfg.back_reflection_model != "none";
}

Vec3 diffuseFrontReflection(std::mt19937_64& rng) {
    const double cos_theta = std::sqrt(uniformOpen(rng));
    const double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
    const double phi = 2.0 * kPi * uniformOpen(rng);
    return {sin_theta * std::cos(phi), sin_theta * std::sin(phi), cos_theta};
}

void reflectAtFront(Vec3& u, const RunConfig& cfg, std::mt19937_64& rng) {
    if (cfg.front_reflection_mode == "specular") {
        u.z = std::abs(u.z);
    } else {
        u = diffuseFrontReflection(rng);
    }
}

double aluminumFresnelReflectance(const Vec3& u, const RunConfig& cfg,
                                  const OpticalProperties& op) {
    const double n1 = std::isfinite(op.n_eff) && op.n_eff > 0.0 ? op.n_eff : 1.0;
    const double cos_i = std::min(1.0, std::max(0.0, -u.z));
    const double sin2_i = std::max(0.0, 1.0 - cos_i * cos_i);
    const std::complex<double> n2(cfg.front_aluminum_n, cfg.front_aluminum_k);
    const std::complex<double> sin_t = (n1 / n2) * std::sqrt(sin2_i);
    const std::complex<double> cos_t = std::sqrt(1.0 - sin_t * sin_t);
    const std::complex<double> rs =
        (n1 * cos_i - n2 * cos_t) / (n1 * cos_i + n2 * cos_t);
    const std::complex<double> rp =
        (n2 * cos_i - n1 * cos_t) / (n2 * cos_i + n1 * cos_t);
    const double reflectance = 0.5 * (std::norm(rs) + std::norm(rp));
    return std::min(1.0, std::max(0.0, reflectance));
}

double frontReflectionProbability(const Vec3& u, const RunConfig& cfg,
                                  const OpticalProperties& op) {
    if (cfg.front_reflection_model == "aluminum_fresnel") {
        return aluminumFresnelReflectance(u, cfg, op);
    }
    return cfg.front_reflectance;
}

double dielectricFresnelReflectance(double n1, double n2, double cos_i) {
    cos_i = std::min(1.0, std::max(0.0, cos_i));
    const double sin2_i = std::max(0.0, 1.0 - cos_i * cos_i);
    const double eta = n1 / n2;
    const double sin2_t = eta * eta * sin2_i;
    if (sin2_t >= 1.0) {
        return 1.0;
    }
    const double cos_t = std::sqrt(std::max(0.0, 1.0 - sin2_t));
    const double rs_den = n1 * cos_i + n2 * cos_t;
    const double rp_den = n2 * cos_i + n1 * cos_t;
    const double rs = rs_den != 0.0 ? (n1 * cos_i - n2 * cos_t) / rs_den : 1.0;
    const double rp = rp_den != 0.0 ? (n2 * cos_i - n1 * cos_t) / rp_den : 1.0;
    const double reflectance = 0.5 * (rs * rs + rp * rp);
    return std::min(1.0, std::max(0.0, reflectance));
}

double backReflectionProbability(const Vec3& u, const RunConfig& cfg,
                                 const OpticalProperties& op) {
    if (cfg.back_reflection_model == "none") {
        return 0.0;
    }
    const double n1 = std::isfinite(op.n_eff) && op.n_eff > 0.0 ? op.n_eff : 1.0;
    return dielectricFresnelReflectance(n1, cfg.back_air_n, u.z);
}

int binIndex(double value, double range, double bin_size) {
    if (value < -range || value >= range) {
        return -1;
    }
    return static_cast<int>(std::floor((value + range) / bin_size));
}

double binCenter(int index, double range, double bin_size) {
    return -range + (static_cast<double>(index) + 0.5) * bin_size;
}

double psfAnchorX(const SourceStep& source) {
    return source.has_psf_anchor ? source.psf_anchor_x : 0.5 * (source.p0.x + source.p1.x);
}

double psfAnchorY(const SourceStep& source) {
    return source.has_psf_anchor ? source.psf_anchor_y : 0.5 * (source.p0.y + source.p1.y);
}

void recordDetected(Accumulators& acc, const SourceStep& source, int source_index,
                    const PhotonRecord& rec, int n_bins, const RunConfig& cfg) {
    const double dx = rec.x - psfAnchorX(source);
    const double dy = rec.y - psfAnchorY(source);

    acc.detected_weight += rec.weight;
    acc.detected_sum_x += rec.weight * dx;
    acc.detected_sum_y += rec.weight * dy;
    acc.detected_sum_x2 += rec.weight * dx * dx;
    acc.detected_sum_y2 += rec.weight * dy * dy;

    auto& ev = acc.event_stats[source.source_event_uid];
    ev.detected_weight += rec.weight;
    ev.sum_x += rec.weight * dx;
    ev.sum_y += rec.weight * dy;
    ev.sum_x2 += rec.weight * dx * dx;
    ev.sum_y2 += rec.weight * dy * dy;

    const int bx = binIndex(dx, cfg.psf_range_um, cfg.psf_bin_size_um);
    const int by = binIndex(dy, cfg.psf_range_um, cfg.psf_bin_size_um);
    if (bx >= 0) {
        acc.lsf_x[bx] += rec.weight;
    }
    if (by >= 0) {
        acc.lsf_y[by] += rec.weight;
    }
    if (bx >= 0 && by >= 0) {
        acc.psf_2d[static_cast<size_t>(by) * n_bins + bx] += rec.weight;
    }
    if (cfg.output_detected_photons) {
        PhotonRecord copy = rec;
        copy.source_id = source.source_id >= 0 ? source.source_id : source_index;
        copy.source_event_uid = source.source_event_uid;
        acc.detected.push_back(copy);
    }
}

bool handleBoundaryHit(Accumulators& acc, StepStats& ss, const SourceStep& source,
                       int source_index, const RunConfig& cfg, int n_bins,
                       const OpticalProperties& op,
                       double photon_weight, Vec3& p, Vec3& u, double path_length,
                       int scatters, const std::string& boundary_surface,
                       std::mt19937_64& rng) {
    if (boundary_surface == "front" && hasFrontReflection(cfg)) {
        p.z = 0.0;
        if (uniformOpen(rng) < frontReflectionProbability(u, cfg, op)) {
            reflectAtFront(u, cfg, rng);
            return true;
        }
        ss.absorbed_weight += photon_weight;
        acc.absorbed_weight += photon_weight;
        ss.weighted_path_length += photon_weight * path_length;
        ss.weighted_scatter_count += photon_weight * scatters;
        return false;
    }

    if (boundary_surface == "back" && hasBackReflection(cfg)) {
        p.z = cfg.thickness_um;
        if (uniformOpen(rng) < backReflectionProbability(u, cfg, op)) {
            u.z = -std::abs(u.z);
            return true;
        }
    }

    if (boundary_surface == "front") {
        ss.front_escape_weight += photon_weight;
        acc.front_escape_weight += photon_weight;
    } else {
        ss.back_escape_weight += photon_weight;
        acc.back_escape_weight += photon_weight;
    }
    if (isReadoutSurface(cfg.readout_surface, boundary_surface)) {
        ss.detected_weight += photon_weight;
        PhotonRecord rec;
        rec.x = p.x;
        rec.y = p.y;
        rec.z = p.z;
        rec.ux = u.x;
        rec.uy = u.y;
        rec.uz = u.z;
        rec.weight = photon_weight;
        rec.path_length_um = path_length;
        rec.scatter_count = scatters;
        rec.surface = boundary_surface;
        recordDetected(acc, source, source_index, rec, n_bins, cfg);
    }
    ss.weighted_path_length += photon_weight * path_length;
    ss.weighted_scatter_count += photon_weight * scatters;
    return false;
}

Accumulators runRange(const std::vector<SourceStep>& sources, size_t begin, size_t end,
                      const RunConfig& cfg, const OpticalProperties& op, int n_bins,
                      uint64_t seed_offset) {
    Accumulators acc;
    acc.step_stats.resize(sources.size());
    acc.psf_2d.assign(static_cast<size_t>(n_bins) * n_bins, 0.0);
    acc.lsf_x.assign(n_bins, 0.0);
    acc.lsf_y.assign(n_bins, 0.0);
    std::mt19937_64 rng(cfg.random_seed + 0x9E3779B97F4A7C15ull * (seed_offset + 1));
    const double mu_t = op.mu_a_per_um + op.mu_s_per_um;
    if (mu_t < 0.0 || op.mu_a_per_um < 0.0 || op.mu_s_per_um < 0.0) {
        throw std::runtime_error("optical coefficients must be non-negative");
    }

    for (size_t si = begin; si < end; ++si) {
        const auto& source = sources[si];
        auto& ss = acc.step_stats[si];
        if (source.n_photon_step <= 0.0) {
            continue;
        }
        const double photon_weight = source.n_photon_step / static_cast<double>(cfg.samples_per_step);
        ss.launched_weight = source.n_photon_step;
        ss.samples = cfg.samples_per_step;
        acc.total_source_weight += source.n_photon_step;

        for (int sample = 0; sample < cfg.samples_per_step; ++sample) {
            const double t = uniformOpen(rng);
            Vec3 p{
                source.p0.x + t * (source.p1.x - source.p0.x),
                source.p0.y + t * (source.p1.y - source.p0.y),
                source.p0.z + t * (source.p1.z - source.p0.z),
            };
            Vec3 u = isotropicDirection(rng);
            double path_length = 0.0;
            int scatters = 0;
            bool alive = true;

            for (int step = 0; alive && step < cfg.max_steps; ++step) {
                double boundary_distance = std::numeric_limits<double>::infinity();
                std::string boundary_surface;
                if (u.z > 0.0) {
                    boundary_distance = (cfg.thickness_um - p.z) / u.z;
                    boundary_surface = "back";
                } else if (u.z < 0.0) {
                    boundary_distance = -p.z / u.z;
                    boundary_surface = "front";
                }
                if (boundary_distance < 0.0) {
                    boundary_distance = 0.0;
                }

                if (mu_t <= 0.0) {
                    p.x += boundary_distance * u.x;
                    p.y += boundary_distance * u.y;
                    p.z += boundary_distance * u.z;
                    path_length += boundary_distance;
                    alive = handleBoundaryHit(acc, ss, source, static_cast<int>(si), cfg, n_bins,
                                              op,
                                              photon_weight, p, u, path_length, scatters,
                                              boundary_surface, rng);
                    continue;
                }

                const double free_path = -std::log(uniformOpen(rng)) / mu_t;
                if (boundary_distance <= free_path) {
                    p.x += boundary_distance * u.x;
                    p.y += boundary_distance * u.y;
                    p.z += boundary_distance * u.z;
                    path_length += boundary_distance;
                    alive = handleBoundaryHit(acc, ss, source, static_cast<int>(si), cfg, n_bins,
                                              op,
                                              photon_weight, p, u, path_length, scatters,
                                              boundary_surface, rng);
                } else {
                    p.x += free_path * u.x;
                    p.y += free_path * u.y;
                    p.z += free_path * u.z;
                    path_length += free_path;
                    if (uniformOpen(rng) < op.mu_a_per_um / mu_t) {
                        ss.absorbed_weight += photon_weight;
                        acc.absorbed_weight += photon_weight;
                        ss.weighted_path_length += photon_weight * path_length;
                        ss.weighted_scatter_count += photon_weight * scatters;
                        alive = false;
                    } else {
                        u = scatterPhoton(u, op, rng);
                        ++scatters;
                    }
                }
            }

            if (alive) {
                ss.lost_weight += photon_weight;
                acc.lost_weight += photon_weight;
                ss.weighted_path_length += photon_weight * path_length;
                ss.weighted_scatter_count += photon_weight * scatters;
            }
        }
    }
    return acc;
}

void mergeAcc(Accumulators& dst, const Accumulators& src) {
    if (dst.step_stats.empty()) {
        dst.step_stats = src.step_stats;
        dst.psf_2d = src.psf_2d;
        dst.lsf_x = src.lsf_x;
        dst.lsf_y = src.lsf_y;
    } else {
        for (size_t i = 0; i < dst.step_stats.size(); ++i) {
            dst.step_stats[i].detected_weight += src.step_stats[i].detected_weight;
            dst.step_stats[i].absorbed_weight += src.step_stats[i].absorbed_weight;
            dst.step_stats[i].front_escape_weight += src.step_stats[i].front_escape_weight;
            dst.step_stats[i].back_escape_weight += src.step_stats[i].back_escape_weight;
            dst.step_stats[i].lost_weight += src.step_stats[i].lost_weight;
            dst.step_stats[i].weighted_path_length += src.step_stats[i].weighted_path_length;
            dst.step_stats[i].weighted_scatter_count += src.step_stats[i].weighted_scatter_count;
            dst.step_stats[i].launched_weight += src.step_stats[i].launched_weight;
            dst.step_stats[i].samples += src.step_stats[i].samples;
        }
        for (size_t i = 0; i < dst.psf_2d.size(); ++i) {
            dst.psf_2d[i] += src.psf_2d[i];
        }
        for (size_t i = 0; i < dst.lsf_x.size(); ++i) {
            dst.lsf_x[i] += src.lsf_x[i];
            dst.lsf_y[i] += src.lsf_y[i];
        }
    }
    for (const auto& kv : src.event_stats) {
        auto& e = dst.event_stats[kv.first];
        e.detected_weight += kv.second.detected_weight;
        e.sum_x += kv.second.sum_x;
        e.sum_y += kv.second.sum_y;
        e.sum_x2 += kv.second.sum_x2;
        e.sum_y2 += kv.second.sum_y2;
    }
    dst.detected.insert(dst.detected.end(), src.detected.begin(), src.detected.end());
    dst.detected_weight += src.detected_weight;
    dst.front_escape_weight += src.front_escape_weight;
    dst.back_escape_weight += src.back_escape_weight;
    dst.absorbed_weight += src.absorbed_weight;
    dst.lost_weight += src.lost_weight;
    dst.total_source_weight += src.total_source_weight;
    dst.detected_sum_x += src.detected_sum_x;
    dst.detected_sum_y += src.detected_sum_y;
    dst.detected_sum_x2 += src.detected_sum_x2;
    dst.detected_sum_y2 += src.detected_sum_y2;
}

double rmsFromSums(double w, double sx, double sx2) {
    if (w <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double mean = sx / w;
    return std::sqrt(std::max(0.0, sx2 / w - mean * mean));
}

double inferIncidentEventCount(const std::vector<EventSource>& events) {
    int64_t max_event_id = -1;
    for (const auto& event : events) {
        try {
            size_t consumed = 0;
            const auto parsed = std::stoll(event.eventID, &consumed);
            if (consumed == event.eventID.size() && parsed > max_event_id) {
                max_event_id = parsed;
            }
        } catch (...) {
        }
    }
    return max_event_id >= 0 ? static_cast<double>(max_event_id + 1) : 0.0;
}

double fwhmFromLsf(const std::vector<double>& lsf, double range, double bin_size) {
    const auto max_it = std::max_element(lsf.begin(), lsf.end());
    if (max_it == lsf.end() || *max_it <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double half = 0.5 * (*max_it);
    const int peak = static_cast<int>(std::distance(lsf.begin(), max_it));

    int left = peak;
    while (left > 0 && lsf[left] >= half) {
        --left;
    }
    int right = peak;
    while (right + 1 < static_cast<int>(lsf.size()) && lsf[right] >= half) {
        ++right;
    }
    const bool left_clipped = (left == 0 && lsf[left] >= half);
    const bool right_clipped =
        (right + 1 >= static_cast<int>(lsf.size()) && lsf[right] >= half);
    if (left_clipped || right_clipped) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    auto interp = [&](int i0, int i1) {
        const double x0 = binCenter(i0, range, bin_size);
        const double x1 = binCenter(i1, range, bin_size);
        const double y0 = lsf[i0];
        const double y1 = lsf[i1];
        if (std::abs(y1 - y0) < 1.0e-15) {
            return 0.5 * (x0 + x1);
        }
        return x0 + (half - y0) * (x1 - x0) / (y1 - y0);
    };

    const double x_left = (left == peak) ? binCenter(peak, range, bin_size)
                                         : interp(left, left + 1);
    double x_right = binCenter(right, range, bin_size);
    if (right > peak && right < static_cast<int>(lsf.size()) && lsf[right] < half) {
        x_right = interp(right - 1, right);
    }
    return std::max(0.0, x_right - x_left);
}

std::string pathJoin(const std::string& dir, const std::string& file) {
    if (dir.empty() || dir == ".") {
        return file;
    }
    const char last = dir.back();
    if (last == '/' || last == '\\') {
        return dir + file;
    }
    return dir + "/" + file;
}

void writeSummary(const std::string& path, const RunConfig& cfg, const OpticalProperties& op,
                  const std::vector<EventSource>& events, const Accumulators& acc) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("cannot write " + path);
    }
    const double total_capture_light = std::accumulate(
        events.begin(), events.end(), 0.0,
        [](double s, const EventSource& e) { return s + e.total_n_photon; });
    const double n_events = static_cast<double>(events.size());
    const double incident_events = cfg.incident_event_count > 0.0
                                       ? cfg.incident_event_count
                                       : inferIncidentEventCount(events);
    const double rms_x = rmsFromSums(acc.detected_weight, acc.detected_sum_x, acc.detected_sum_x2);
    const double rms_y = rmsFromSums(acc.detected_weight, acc.detected_sum_y, acc.detected_sum_y2);
    const double rms_r = std::isfinite(rms_x) && std::isfinite(rms_y)
                             ? std::sqrt(rms_x * rms_x + rms_y * rms_y)
                             : std::numeric_limits<double>::quiet_NaN();
    const double fwhm_x = fwhmFromLsf(acc.lsf_x, cfg.psf_range_um, cfg.psf_bin_size_um);
    const double fwhm_y = fwhmFromLsf(acc.lsf_y, cfg.psf_range_um, cfg.psf_bin_size_um);
    const double mu_s_prime_from_g = op.mu_s_per_um * (1.0 - op.g);
    const double mu_s_prime = std::isfinite(op.mu_s_prime_per_um) ? op.mu_s_prime_per_um
                                                                  : mu_s_prime_from_g;
    const double mu_tr = op.mu_a_per_um + mu_s_prime;
    const double absorption_length_um =
        op.mu_a_per_um > 0.0 ? 1.0 / op.mu_a_per_um : std::numeric_limits<double>::infinity();
    const double scattering_length_um =
        op.mu_s_per_um > 0.0 ? 1.0 / op.mu_s_per_um : std::numeric_limits<double>::infinity();
    const double transport_mfp_um =
        mu_tr > 0.0 ? 1.0 / mu_tr : std::numeric_limits<double>::infinity();
    const double diffusion_length_um =
        (op.mu_a_per_um > 0.0 && mu_tr > 0.0)
            ? std::sqrt(1.0 / (3.0 * op.mu_a_per_um * mu_tr))
            : std::numeric_limits<double>::infinity();

    out << "ratio_tag,wavelength_nm,thickness_um,mu_a_per_um,mu_s_per_um,g,n_eff,"
           "phase_function_mode,phase_function_csv,"
           "front_reflection_model,front_reflectance,front_reflection_mode,"
           "front_aluminum_n,front_aluminum_k,back_reflection_model,back_air_n,"
           "mu_s_prime_per_um,mu_s_prime_from_g_per_um,mu_tr_per_um,absorption_length_um,scattering_length_um,transport_mfp_um,diffusion_length_um,"
           "n_events,incident_event_count,capture_fraction,n_source_steps,samples_per_step,total_source_weight,total_detected_weight,"
           "front_escape_weight,back_escape_weight,absorbed_weight,lost_weight,"
           "mean_light_per_capture,mean_detected_light_per_capture,mean_light_per_incident,mean_detected_light_per_incident,detection_efficiency,"
           "spot_rms_x,spot_rms_y,spot_rms_r,fwhm_x,fwhm_y\n";
    out << csvEscape(cfg.ratio_tag.empty() ? op.ratio_tag : cfg.ratio_tag) << ','
        << num(op.wavelength_nm) << ',' << num(cfg.thickness_um) << ','
        << num(op.mu_a_per_um) << ',' << num(op.mu_s_per_um) << ','
        << num(op.g) << ',' << num(op.n_eff) << ','
        << csvEscape(op.phase_cdf.empty() ? "HG" : "tabulated_mu") << ','
        << csvEscape(op.phase_function_csv) << ','
        << csvEscape(cfg.front_reflection_model) << ',' << num(cfg.front_reflectance) << ','
        << csvEscape(cfg.front_reflection_mode) << ','
        << num(cfg.front_aluminum_n) << ',' << num(cfg.front_aluminum_k) << ','
        << csvEscape(cfg.back_reflection_model) << ',' << num(cfg.back_air_n) << ','
        << num(mu_s_prime) << ',' << num(mu_s_prime_from_g) << ',' << num(mu_tr) << ','
        << num(absorption_length_um) << ',' << num(scattering_length_um) << ','
        << num(transport_mfp_um) << ',' << num(diffusion_length_um) << ','
        << events.size() << ',' << num(incident_events) << ','
        << num(incident_events > 0.0 ? n_events / incident_events : 0.0) << ','
        << acc.step_stats.size() << ',' << cfg.samples_per_step << ','
        << num(acc.total_source_weight) << ',' << num(acc.detected_weight) << ','
        << num(acc.front_escape_weight) << ',' << num(acc.back_escape_weight) << ','
        << num(acc.absorbed_weight) << ',' << num(acc.lost_weight) << ','
        << num(n_events > 0.0 ? total_capture_light / n_events : 0.0) << ','
        << num(n_events > 0.0 ? acc.detected_weight / n_events : 0.0) << ','
        << num(incident_events > 0.0 ? total_capture_light / incident_events : 0.0) << ','
        << num(incident_events > 0.0 ? acc.detected_weight / incident_events : 0.0) << ','
        << num(acc.total_source_weight > 0.0 ? acc.detected_weight / acc.total_source_weight : 0.0)
        << ',' << num(rms_x) << ',' << num(rms_y) << ',' << num(rms_r) << ','
        << num(fwhm_x) << ',' << num(fwhm_y) << '\n';
}

void writeEventSummary(const std::string& path, const std::vector<EventSource>& events,
                       const Accumulators& acc) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("cannot write " + path);
    }
    out << "source_event_uid,eventID,depth_um,total_n_photon,detected_weight,"
           "detection_efficiency,centroid_x,centroid_y,spot_rms_x,spot_rms_y,spot_rms_r\n";
    for (const auto& event : events) {
        EventStats stats;
        const auto it = acc.event_stats.find(event.source_event_uid);
        if (it != acc.event_stats.end()) {
            stats = it->second;
        }
        const double cx = stats.detected_weight > 0.0 ? stats.sum_x / stats.detected_weight
                                                      : std::numeric_limits<double>::quiet_NaN();
        const double cy = stats.detected_weight > 0.0 ? stats.sum_y / stats.detected_weight
                                                      : std::numeric_limits<double>::quiet_NaN();
        const double rx = rmsFromSums(stats.detected_weight, stats.sum_x, stats.sum_x2);
        const double ry = rmsFromSums(stats.detected_weight, stats.sum_y, stats.sum_y2);
        const double rr = std::isfinite(rx) && std::isfinite(ry)
                              ? std::sqrt(rx * rx + ry * ry)
                              : std::numeric_limits<double>::quiet_NaN();
        out << csvEscape(event.source_event_uid) << ',' << csvEscape(event.eventID) << ','
            << num(event.depth_um) << ',' << num(event.total_n_photon) << ','
            << num(stats.detected_weight) << ','
            << num(event.total_n_photon > 0.0 ? stats.detected_weight / event.total_n_photon : 0.0)
            << ',' << num(cx) << ',' << num(cy) << ',' << num(rx) << ',' << num(ry) << ','
            << num(rr) << '\n';
    }
}

void writeStepSummary(const std::string& path, const std::vector<SourceStep>& sources,
                      const Accumulators& acc) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("cannot write " + path);
    }
    out << "source_id,source_event_uid,eventID,detected_weight,absorbed_weight,escape_weight,"
           "front_escape_weight,back_escape_weight,lost_weight,mean_path_length,mean_scatter_count,"
           "samples\n";
    for (size_t i = 0; i < sources.size(); ++i) {
        const auto& s = sources[i];
        const auto& st = acc.step_stats[i];
        const double denom = st.launched_weight > 0.0 ? st.launched_weight : 0.0;
        out << s.source_id << ',' << csvEscape(s.source_event_uid) << ',' << csvEscape(s.eventID)
            << ',' << num(st.detected_weight) << ',' << num(st.absorbed_weight) << ','
            << num(st.front_escape_weight + st.back_escape_weight) << ','
            << num(st.front_escape_weight) << ',' << num(st.back_escape_weight) << ','
            << num(st.lost_weight) << ','
            << num(denom > 0.0 ? st.weighted_path_length / denom : 0.0) << ','
            << num(denom > 0.0 ? st.weighted_scatter_count / denom : 0.0) << ','
            << st.samples << '\n';
    }
}

void writeDetectedPhotons(const std::string& path, const std::vector<PhotonRecord>& photons) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("cannot write " + path);
    }
    out << "source_event_uid,source_id,readout_surface,readout_x_um,readout_y_um,readout_z_um,"
           "dir_x,dir_y,dir_z,weight,path_length_um,scatter_count\n";
    for (const auto& p : photons) {
        out << csvEscape(p.source_event_uid) << ',' << p.source_id << ','
            << csvEscape(p.surface) << ',' << num(p.x) << ',' << num(p.y) << ','
            << num(p.z) << ',' << num(p.ux) << ',' << num(p.uy) << ',' << num(p.uz)
            << ',' << num(p.weight) << ',' << num(p.path_length_um) << ','
            << p.scatter_count << '\n';
    }
}

void writePsfLsf(const RunConfig& cfg, const Accumulators& acc, int n_bins) {
    {
        std::ofstream out(pathJoin(cfg.output_dir, "psf_2d.csv"));
        if (!out) throw std::runtime_error("cannot write psf_2d.csv");
        out << "x_bin_center_um,y_bin_center_um,weight\n";
        for (int y = 0; y < n_bins; ++y) {
            for (int x = 0; x < n_bins; ++x) {
                out << num(binCenter(x, cfg.psf_range_um, cfg.psf_bin_size_um)) << ','
                    << num(binCenter(y, cfg.psf_range_um, cfg.psf_bin_size_um)) << ','
                    << num(acc.psf_2d[static_cast<size_t>(y) * n_bins + x]) << '\n';
            }
        }
    }
    {
        std::ofstream out(pathJoin(cfg.output_dir, "lsf_x.csv"));
        if (!out) throw std::runtime_error("cannot write lsf_x.csv");
        out << "x_bin_center_um,weight\n";
        for (int x = 0; x < n_bins; ++x) {
            out << num(binCenter(x, cfg.psf_range_um, cfg.psf_bin_size_um)) << ','
                << num(acc.lsf_x[x]) << '\n';
        }
    }
    {
        std::ofstream out(pathJoin(cfg.output_dir, "lsf_y.csv"));
        if (!out) throw std::runtime_error("cannot write lsf_y.csv");
        out << "y_bin_center_um,weight\n";
        for (int y = 0; y < n_bins; ++y) {
            out << num(binCenter(y, cfg.psf_range_um, cfg.psf_bin_size_um)) << ','
                << num(acc.lsf_y[y]) << '\n';
        }
    }
}

void ensureOutputDirectory(const std::string& dir) {
#ifdef _WIN32
    const std::string command = "if not exist \"" + dir + "\" mkdir \"" + dir + "\"";
#else
    const std::string command = "mkdir -p \"" + dir + "\"";
#endif
    if (!dir.empty() && dir != ".") {
        std::system(command.c_str());
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 2) {
            std::cerr << "usage: OpticalMC run_config.json [source_steps.csv event_sources.csv "
                         "optical_properties.csv output_dir]\n";
            return 2;
        }
        RunConfig cfg = readConfig(argv[1]);
        if (argc > 2) cfg.source_steps_csv = argv[2];
        if (argc > 3) cfg.event_sources_csv = argv[3];
        if (argc > 4) cfg.optical_properties_csv = argv[4];
        if (argc > 5) cfg.output_dir = argv[5];
        if (cfg.source_steps_csv.empty() || cfg.event_sources_csv.empty()) {
            throw std::runtime_error(
                "source_steps_csv and event_sources_csv must be in run_config.json or CLI args");
        }
        if (cfg.readout_surface != "front" && cfg.readout_surface != "back" &&
            cfg.readout_surface != "both") {
            throw std::runtime_error("readout_surface must be front, back, or both");
        }
        if (cfg.xy_boundary != "infinite") {
            std::cerr << "warning: only xy_boundary=\"infinite\" is currently implemented; "
                         "continuing with infinite x/y extent\n";
        }
        ensureOutputDirectory(cfg.output_dir);

        auto sources = readSources(cfg.source_steps_csv);
        auto events = readEvents(cfg.event_sources_csv);
        applyEventAnchors(sources, events);
        auto op = selectOpticalProperties(cfg.optical_properties_csv, cfg, sources);
        if (cfg.ratio_tag.empty()) {
            cfg.ratio_tag = op.ratio_tag;
        }
        const int n_bins = static_cast<int>(std::ceil(2.0 * cfg.psf_range_um / cfg.psf_bin_size_um));
        if (n_bins <= 0) {
            throw std::runtime_error("invalid PSF bin configuration");
        }

        std::cout << "sources," << sources.size() << "\n";
        std::cout << "events," << events.size() << "\n";
        std::cout << "optical_mu_a_per_um," << op.mu_a_per_um << "\n";
        std::cout << "optical_mu_s_per_um," << op.mu_s_per_um << "\n";
        std::cout << "optical_phase_function_mode,"
                  << (op.phase_cdf.empty() ? "HG" : "tabulated_mu") << "\n";

        const int n_threads = std::max(1, std::min<int>(cfg.num_threads, sources.empty() ? 1 : static_cast<int>(sources.size())));
        std::vector<std::thread> threads;
        std::vector<Accumulators> partials(static_cast<size_t>(n_threads));
        std::exception_ptr thread_error = nullptr;
        std::mutex error_mutex;

        for (int tid = 0; tid < n_threads; ++tid) {
            const size_t begin = sources.size() * static_cast<size_t>(tid) / static_cast<size_t>(n_threads);
            const size_t end = sources.size() * static_cast<size_t>(tid + 1) / static_cast<size_t>(n_threads);
            threads.emplace_back([&, tid, begin, end]() {
                try {
                    partials[tid] = runRange(sources, begin, end, cfg, op, n_bins, static_cast<uint64_t>(tid));
                } catch (...) {
                    std::lock_guard<std::mutex> lock(error_mutex);
                    if (!thread_error) {
                        thread_error = std::current_exception();
                    }
                }
            });
        }
        for (auto& t : threads) {
            t.join();
        }
        if (thread_error) {
            std::rethrow_exception(thread_error);
        }

        Accumulators acc;
        for (const auto& part : partials) {
            mergeAcc(acc, part);
        }

        writeSummary(pathJoin(cfg.output_dir, "optical_mc_summary.csv"), cfg, op, events, acc);
        writeEventSummary(pathJoin(cfg.output_dir, "optical_mc_event_summary.csv"), events, acc);
        writeStepSummary(pathJoin(cfg.output_dir, "optical_mc_source_step_summary.csv"), sources, acc);
        writePsfLsf(cfg, acc, n_bins);
        if (cfg.output_detected_photons) {
            writeDetectedPhotons(pathJoin(cfg.output_dir, "detected_photons.csv"), acc.detected);
        }

        std::cout << "detected_weight," << acc.detected_weight << "\n";
        std::cout << "wrote," << pathJoin(cfg.output_dir, "optical_mc_summary.csv") << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
