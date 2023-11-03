#include "solver/nlopt_solve.h"

const double FUNCTION_TOLERANCE = 1e-7;
const double SMALL_CONSTRAINT_TOLERANCE = 1e-6;
const double LARGE_CONSTRAINT_TOLERANCE = 1e-2;
const double MAX_FUNC_EVALS = 100000;

QPair<double, std::vector<double>> Nlopt_Solve::solve(optimization_data_s optimization_data)
{
    std::vector<double> solution = optimization_data.initial_guess;
    size_t n = solution.size();

    nlopt::opt opt(nlopt::LD_MMA, n);

    std::vector<double> lb;
    std::vector<double> ub;

    lb.resize(n);
    ub.resize(n);

    double ub_val = 20.0;
    double lb_val = -20.0;

    for (size_t i = 0; i < n; ++i) {
        ub[i] = ub_val;
        lb[i] = lb_val;
    }

    opt.set_upper_bounds(ub);
    opt.set_lower_bounds(lb);

    minimization_data_s* min_data = new minimization_data_s(optimization_data.intersection_datas);
    opt.set_min_objective(minimization_function, min_data);
    opt.set_ftol_abs(FUNCTION_TOLERANCE);
    opt.set_maxeval(MAX_FUNC_EVALS);

    std::vector<constraint_orthogonality_data_s*> constraint_orthogonality_datas;
    for (size_t i = 0; i < optimization_data.intersection_datas.size(); ++i) {
        constraint_orthogonality_data_s* constraint_orthogonality_data = new constraint_orthogonality_data_s(optimization_data.intersection_datas[i]);
        constraint_orthogonality_datas.push_back(constraint_orthogonality_data);

        opt.add_inequality_constraint(constraint_function_tangent_orthogonality_pos, constraint_orthogonality_data, SMALL_CONSTRAINT_TOLERANCE);
        opt.add_inequality_constraint(constraint_function_tangent_orthogonality_neg, constraint_orthogonality_data, SMALL_CONSTRAINT_TOLERANCE);
        opt.add_inequality_constraint(constraint_function_plane_orthogonality_pos, constraint_orthogonality_data, SMALL_CONSTRAINT_TOLERANCE);
        opt.add_inequality_constraint(constraint_function_plane_orthogonality_neg, constraint_orthogonality_data, SMALL_CONSTRAINT_TOLERANCE);
    }

    std::vector<constraint_offset_data_s*> constraint_offset_datas;
    for (size_t i = 1; i < optimization_data.intersection_datas.size(); ++i) {
        constraint_offset_data_s* constraint_offset_data = new constraint_offset_data_s(optimization_data.intersection_datas[0], optimization_data.intersection_datas[i]);
        constraint_offset_datas.push_back(constraint_offset_data);
        opt.add_inequality_constraint(constraint_function_offset_pos, constraint_offset_data, SMALL_CONSTRAINT_TOLERANCE);
        opt.add_inequality_constraint(constraint_function_offset_neg, constraint_offset_data, SMALL_CONSTRAINT_TOLERANCE);
    }

    // TODO this is technically not ortho, adjust structs
    opt.add_inequality_constraint(constraint_function_isec0_ca_cb_pos, constraint_orthogonality_datas[0], SMALL_CONSTRAINT_TOLERANCE);
    opt.add_inequality_constraint(constraint_function_isec0_ca_cb_neg, constraint_orthogonality_datas[0], SMALL_CONSTRAINT_TOLERANCE);

    double minf;
    nlopt::result result;

    try {
        result = opt.optimize(solution, minf);
    }
    catch(std::exception &e) {
        std::cout << ">> nlopt failed: " << e.what() << std::endl;
    }

    delete min_data;
    constraint_orthogonality_datas.clear();
    constraint_offset_datas.clear();
    return qMakePair(minf, solution);
}

double Nlopt_Solve::minimization_function(const std::vector<double> &x, std::vector<double> &grad, void *min_data)
{
    minimization_data_s *data = reinterpret_cast<minimization_data_s*>(min_data);
    std::vector<intersection_data_s> intersection_datas = data->intersection_datas;

    double out = 0.0; // minfunc

    for (size_t i = 0; i < intersection_datas.size(); ++i) {
        intersection_data_s intersection_data = intersection_datas[i];

        vector_t atan = intersection_data.tangent_a;
        vector_t btan = intersection_data.tangent_b;

        unsigned int an_x = intersection_data.polybezier_a_index * 3 + 0;
        unsigned int an_y = intersection_data.polybezier_a_index * 3 + 1;
        unsigned int an_c = intersection_data.polybezier_a_index * 3 + 2;
        unsigned int bn_x = intersection_data.polybezier_b_index * 3 + 0;
        unsigned int bn_y = intersection_data.polybezier_b_index * 3 + 1;
        unsigned int bn_c = intersection_data.polybezier_b_index * 3 + 2;

        // -------------------
        // MINIMIZATION ENERGY
        // ---------------------------------------------------------------------------
        // (6) MIN(ni) = SUM[(||tji x ni||^2 + ||tij x nj||^2) + (tij.z + tji.z)^2
        // ---------------------------------------------------------------------------

        // --------------------------------------
        // Minimal Foreshortening
        // --------------------------------------
        // (4) (tij.z)^2 + (tij.z)^2
        // (7) t1z = -(t1.x * n1.x + t1.y * n1.y)
        // (7) t2z = -(t2.x * n2.x + t2.y * n2.y)
        // --------------------------------------
        // --------------------------------------------------------------------
        // Cross-sections as Local Geodesics
        // --------------------------------------------------------------------
        // (3) ||tij x nj||^2
        // (7) t1z = -(t1.x * n1.x + t1.y * n1.y)
        // norm(cross([t2.x, t2.y, -(t2.x*n2.x + t2.y*n2.y)], [n1.x, n1.y, 1]))
        // --------------------------------------------------------------------
        // --------------------------------------------------------------------
        // (3) ||tji x ni||^2
        // (7) t2z = -(t2.x * n2.x + t2.y * n2.y)
        // norm(cross([t1.x, t1.y, -(t1.x*n1.x + t1.y*n1.y)], [n2.x, n2.y, 1]))
        // --------------------------------------------------------------------
        out += (btan.y - x[an_y] * (-1.0 * btan.x * x[bn_x] - btan.y * x[bn_y])) * (btan.y - x[an_y] * (-1.0 * btan.x * x[bn_x] - btan.y * x[bn_y])) + (x[an_x] * (-btan.x * x[bn_x] - btan.y * x[bn_y]) - btan.x) * (x[an_x] * (-btan.x * x[bn_x] - btan.y * x[bn_y]) - btan.x) + (btan.x * x[an_y] - btan.y * x[an_x]) * (btan.x * x[an_y] - btan.y * x[an_x]) + (atan.y - x[bn_y] * (-1.0 * atan.x * x[an_x] - atan.y * x[an_y])) * (atan.y - x[bn_y] * (-1.0 * atan.x * x[an_x] - atan.y * x[an_y])) + (x[bn_x] * (-1.0 * atan.x * x[an_x] - atan.y * x[an_y]) - atan.x) * (x[bn_x] * (-1.0 * atan.x * x[an_x] - atan.y * x[an_y]) - atan.x) + (atan.x * x[bn_y] - atan.y * x[bn_x]) * (atan.x * x[bn_y] - atan.y * x[bn_x]) + (-1.0 * atan.x * x[an_x] - atan.y * x[an_y] - btan.x * x[bn_x] - btan.y * x[bn_y]) * (-1.0 * atan.x * x[an_x] - atan.y * x[an_y] - btan.x * x[bn_x] - btan.y * x[bn_y]);

        if (!grad.empty()) {
            grad[an_x] += 2.0 * x[bn_y] * atan.x * (atan.y - x[bn_y] * (-1.0 * x[an_y] * atan.y - x[an_x] * atan.x)) - 2.0 * x[bn_x] * atan.x * (x[bn_x] * (-1.0 * x[an_y] * atan.y - x[an_x] * atan.x) - atan.x) - 2.0 * atan.x * (-1.0 * btan.y * x[bn_y] - x[an_y] * atan.y - btan.x * x[bn_x] - x[an_x] * atan.x) + 2.0 * (-1.0 * btan.y * x[bn_y] - btan.x * x[bn_x]) * (x[an_x] * (-1.0 * btan.y * x[bn_y] - btan.x * x[bn_x]) - btan.x) - 2.0 * btan.y * (x[an_y] * btan.x - btan.y * x[an_x]);
            grad[an_y] += 2.0 * x[bn_y] * atan.y * (atan.y - x[bn_y] * (-1.0 * x[an_y] * atan.y - x[an_x] * atan.x)) - 2.0 * x[bn_x] * atan.y * (x[bn_x] * (-1.0 * x[an_y] * atan.y - x[an_x] * atan.x) - atan.x) - 2.0 * atan.y * (-1.0 * btan.y * x[bn_y] - x[an_y] * atan.y - btan.x * x[bn_x] - x[an_x] * atan.x) + 2.0 * btan.x * (x[an_y] * btan.x - btan.y * x[an_x]) + 2.0 * (btan.y * x[bn_y] + btan.x * x[bn_x]) * (btan.y - x[an_y] * (-1.0 * btan.y * x[bn_y] - btan.x * x[bn_x]));
            grad[an_c] += 0.0;
            grad[bn_x] += 2.0 * x[an_y] * btan.x * (btan.y - x[an_y] * (-1.0 * btan.y * x[bn_y] - btan.x * x[bn_x])) - 2.0 * btan.x * x[an_x] * (x[an_x] * (-1.0 * btan.y * x[bn_y] - btan.x * x[bn_x]) - btan.x) - 2.0 * btan.x * (-1.0 * btan.y * x[bn_y] - x[an_y] * atan.y - btan.x * x[bn_x] - x[an_x] * atan.x) + 2.0 * (-1.0 * x[an_y] * atan.y - x[an_x] * atan.x) * (x[bn_x] * (-1.0 * x[an_y] * atan.y - x[an_x] * atan.x) - atan.x) - 2.0 * atan.y * (x[bn_y] * atan.x - x[bn_x] * atan.y);
            grad[bn_y] += 2.0 * btan.y * x[an_y] * (btan.y - x[an_y] * (-1.0 * btan.y * x[bn_y] - btan.x * x[bn_x])) - 2.0 * btan.y * x[an_x] * (x[an_x] * (-1.0 * btan.y * x[bn_y] - btan.x * x[bn_x]) - btan.x) - 2.0 * btan.y * (-1.0 * btan.y * x[bn_y] - x[an_y] * atan.y - btan.x * x[bn_x] - x[an_x] * atan.x) + 2.0 * atan.x * (x[bn_y] * atan.x - x[bn_x] * atan.y) + 2.0 * (x[an_y] * atan.y + x[an_x] * atan.x) * (atan.y - x[bn_y] * (-1.0 * x[an_y] * atan.y - x[an_x] * atan.x));
            grad[bn_c] += 0.0;
        }
    }

    return out; //minfunc
}

double Nlopt_Solve::constraint_function_tangent_orthogonality_pos(const std::vector<double> &x, std::vector<double> &grad, void *con_data)
{
    constraint_orthogonality_data_s *data = reinterpret_cast<constraint_orthogonality_data_s*>(con_data);
    intersection_data_s intersection_data = data->intersection_data;

    double out = 0.0; // confunc

    vector_t atan = intersection_data.tangent_a;
    vector_t btan = intersection_data.tangent_b;

    unsigned int an_x = intersection_data.polybezier_a_index * 3 + 0;
    unsigned int an_y = intersection_data.polybezier_a_index * 3 + 1;
    unsigned int an_c = intersection_data.polybezier_a_index * 3 + 2;
    unsigned int bn_x = intersection_data.polybezier_b_index * 3 + 0;
    unsigned int bn_y = intersection_data.polybezier_b_index * 3 + 1;
    unsigned int bn_c = intersection_data.polybezier_b_index * 3 + 2;

    // --------------------------------
    // TANGENT ORTHOGONALITY CONSTRAINT
    // -----------------------------------------------------------------------------------------
    // (2) tij . tji < 0.1
    // (7) t1z = -(t1.x * n1.x + t1.y * n1.y)
    // (7) t2z = -(t2.x * n2.x + t2.y * n2.y)
    // tij.x * tji.x + tij.y * tji.y + tij.z * tji.z
    // tij.x * tji.x + tij.y * tji.y + -(t1.x * n1.x + t1.y * n1.y) * -(t2.x * n2.x + t2.y * n2.y)
    // -----------------------------------------------------------------------------------------
    out += ((atan.x) * (btan.x) + (atan.y) * (btan.y) + (atan.x) * (x[an_x]) * btan.x * x[bn_x] + atan.x * x[an_x] * btan.y * x[bn_y] + atan.y * x[an_y] * btan.x * x[bn_x] + atan.y * x[an_y] * btan.y * x[bn_y]) - 0.01;

    if (!grad.empty()) {
        grad[an_x] += atan.x * btan.x * x[bn_x] + atan.x * btan.y * x[bn_y];
        grad[an_y] += atan.y * btan.x * x[bn_x] + atan.y * btan.y * x[bn_y];
        grad[an_c] += 0.0;
        grad[bn_x] += atan.x * btan.x * x[an_x] + atan.y * btan.x * x[an_y];
        grad[bn_y] += atan.x * btan.y * x[an_x] + atan.y * btan.y * x[an_y];
        grad[bn_c] += 0.0;
    }

    return out; //confunc
}

double Nlopt_Solve::constraint_function_tangent_orthogonality_neg(const std::vector<double> &x, std::vector<double> &grad, void *con_data)
{
    constraint_orthogonality_data_s *data = reinterpret_cast<constraint_orthogonality_data_s*>(con_data);
    intersection_data_s intersection_data = data->intersection_data;

    double out = 0.0; // confunc

    vector_t atan = intersection_data.tangent_a;
    vector_t btan = intersection_data.tangent_b;

    unsigned int an_x = intersection_data.polybezier_a_index * 3 + 0;
    unsigned int an_y = intersection_data.polybezier_a_index * 3 + 1;
    unsigned int an_c = intersection_data.polybezier_a_index * 3 + 2;
    unsigned int bn_x = intersection_data.polybezier_b_index * 3 + 0;
    unsigned int bn_y = intersection_data.polybezier_b_index * 3 + 1;
    unsigned int bn_c = intersection_data.polybezier_b_index * 3 + 2;

    // --------------------------------
    // TANGENT ORTHOGONALITY CONSTRAINT
    // -----------------------------------------------------------------------------------------
    // (2) tij . tji > -0.1
    // (7) t1z = -(t1.x * n1.x + t1.y * n1.y)
    // (7) t2z = -(t2.x * n2.x + t2.y * n2.y)
    // tij.x * tji.x + tij.y * tji.y + tij.z * tji.z
    // tij.x * tji.x + tij.y * tji.y + -(t1.x * n1.x + t1.y * n1.y) * -(t2.x * n2.x + t2.y * n2.y)
    // -----------------------------------------------------------------------------------------
    out += -1.0 * ((atan.x) * (btan.x) + (atan.y) * (btan.y) + (atan.x) * (x[an_x]) * btan.x * x[bn_x] + atan.x * x[an_x] * btan.y * x[bn_y] + atan.y * x[an_y] * btan.x * x[bn_x] + atan.y * x[an_y] * btan.y * x[bn_y]) - 0.01;

    if (!grad.empty()) {
        grad[an_x] += -1.0 * atan.x * btan.x * x[bn_x] - atan.x * btan.y * x[bn_y];
        grad[an_y] += -1.0 * atan.y * btan.x * x[bn_x] - atan.y * btan.y * x[bn_y];
        grad[an_c] += 0.0;
        grad[bn_x] += -1.0 * atan.x * btan.x * x[an_x] - atan.y * btan.x * x[an_y];
        grad[bn_y] += -1.0 * atan.x * btan.y * x[an_x] - atan.y * btan.y * x[an_y];
        grad[bn_c] += 0.0;
    }

    return out; //confunc
}

double Nlopt_Solve::constraint_function_plane_orthogonality_pos(const std::vector<double> &x, std::vector<double> &grad, void *con_data)
{
    constraint_orthogonality_data_s *data = reinterpret_cast<constraint_orthogonality_data_s*>(con_data);
    intersection_data_s intersection_data = data->intersection_data;

    double out = 0.0; // confunc

    vector_t atan = intersection_data.tangent_a;
    vector_t btan = intersection_data.tangent_b;

    unsigned int an_x = intersection_data.polybezier_a_index * 3 + 0;
    unsigned int an_y = intersection_data.polybezier_a_index * 3 + 1;
    unsigned int an_c = intersection_data.polybezier_a_index * 3 + 2;
    unsigned int bn_x = intersection_data.polybezier_b_index * 3 + 0;
    unsigned int bn_y = intersection_data.polybezier_b_index * 3 + 1;
    unsigned int bn_c = intersection_data.polybezier_b_index * 3 + 2;

    // ------------------------------
    // PLANE ORTHOGONALITY CONSTRAINT
    // -------------------------------------
    // (2) ni . nj < 0.1
    // ni = [ni.x, ni.y, 1.0]
    // nj = [nj.x, nj.y, 1.0]
    // ni.x * nj.x + ni.y * nj.y + 1.0 * 1.0
    // -------------------------------------
    out += ((x[an_x]) * (x[bn_x]) + (x[an_y]) * (x[bn_y]) + 1) - 0.1;

    if (!grad.empty()) {
        grad[an_x] += x[bn_x];
        grad[an_y] += x[bn_y];
        grad[an_c] += 0.0;
        grad[bn_x] += x[an_x];
        grad[bn_y] += x[an_y];
        grad[bn_c] += 0.0;
    }

    return out; // confunc
}

double Nlopt_Solve::constraint_function_plane_orthogonality_neg(const std::vector<double> &x, std::vector<double> &grad, void *con_data)
{
    constraint_orthogonality_data_s *data = reinterpret_cast<constraint_orthogonality_data_s*>(con_data);
    intersection_data_s intersection_data = data->intersection_data;

    double out = 0.0; // confunc

    vector_t atan = intersection_data.tangent_a;
    vector_t btan = intersection_data.tangent_b;

    unsigned int an_x = intersection_data.polybezier_a_index * 3 + 0;
    unsigned int an_y = intersection_data.polybezier_a_index * 3 + 1;
    unsigned int an_c = intersection_data.polybezier_a_index * 3 + 2;
    unsigned int bn_x = intersection_data.polybezier_b_index * 3 + 0;
    unsigned int bn_y = intersection_data.polybezier_b_index * 3 + 1;
    unsigned int bn_c = intersection_data.polybezier_b_index * 3 + 2;

    // ------------------------------
    // PLANE ORTHOGONALITY CONSTRAINT
    // -------------------------------------
    // (2) ni . nj > -0.1
    // ni = [ni.x, ni.y, 1.0]
    // nj = [nj.x, nj.y, 1.0]
    // ni.x * nj.x + ni.y * nj.y + 1.0 * 1.0
    // -------------------------------------
    out += -1.0 * ((x[an_x]) * (x[bn_x]) + (x[an_y]) * (x[bn_y]) + 1) - 0.1;

    if (!grad.empty()) {
        grad[an_x] += -1.0 * x[bn_x];
        grad[an_y] += -1.0 * x[bn_y];
        grad[an_c] += 0.0;
        grad[bn_x] += -1.0 * x[an_x];
        grad[bn_y] += -1.0 * x[an_y];
        grad[bn_c] += 0.0;
    }

    return out; //confunc
}

double Nlopt_Solve::constraint_function_offset_pos(const std::vector<double> &x, std::vector<double> &grad, void *con_data)
{
    constraint_offset_data_s *data = reinterpret_cast<constraint_offset_data_s*>(con_data);
    intersection_data_s intersection_data_fixed = data->intersection_data_fixed;
    intersection_data_s intersection_data = data->intersection_data;

    double out = 0.0; // confunc

    unsigned int an_x = intersection_data.polybezier_a_index * 3 + 0;
    unsigned int an_y = intersection_data.polybezier_a_index * 3 + 1;
    unsigned int an_c = intersection_data.polybezier_a_index * 3 + 2;
    unsigned int bn_x = intersection_data.polybezier_b_index * 3 + 0;
    unsigned int bn_y = intersection_data.polybezier_b_index * 3 + 1;
    unsigned int bn_c = intersection_data.polybezier_b_index * 3 + 2;

    point_t p = intersection_data_fixed.position - intersection_data.position;

    // ------------------------------
    // CROSSHAIR 3D OFFSET CONSTRAINT
    // ------------------------------------------------------------------------
    // (8) xi * ni + ci = 0 && xj * nj + cj = 0
    // (9) xi *(ni.x − nj.x) + xj * (ni.y − nj.y) + ( ci − cj ) = 0
    // ni = [ni.x, ni.y, 1.0]
    // nj = [nj.x, nj.y, 1.0]
    // (*) where xi and xj are any points on the ni and nj planes, respectively
    // ------------------------------------------------------------------------
    out += p.x * ((x[an_x]) - (x[bn_x])) + p.y * ((x[an_y]) - (x[bn_y])) + (x[an_c] - x[bn_c]);

    if (!grad.empty()) {
        grad[an_x] += p.x;
        grad[an_y] += p.y;
        grad[an_c] += 1.0;
        grad[bn_x] += -1.0 * p.x;
        grad[bn_y] += -1.0 * p.y;
        grad[bn_c] += -1.0;
    }

    return out; // confunc
}

double Nlopt_Solve::constraint_function_offset_neg(const std::vector<double> &x, std::vector<double> &grad, void *con_data)
{
    constraint_offset_data_s *data = reinterpret_cast<constraint_offset_data_s*>(con_data);
    intersection_data_s intersection_data_fixed = data->intersection_data_fixed;
    intersection_data_s intersection_data = data->intersection_data;

    double out = 0.0; // confunc

    unsigned int an_x = intersection_data.polybezier_a_index * 3 + 0;
    unsigned int an_y = intersection_data.polybezier_a_index * 3 + 1;
    unsigned int an_c = intersection_data.polybezier_a_index * 3 + 2;
    unsigned int bn_x = intersection_data.polybezier_b_index * 3 + 0;
    unsigned int bn_y = intersection_data.polybezier_b_index * 3 + 1;
    unsigned int bn_c = intersection_data.polybezier_b_index * 3 + 2;

    point_t p = intersection_data_fixed.position - intersection_data.position;

    // ------------------------------
    // CROSSHAIR 3D OFFSET CONSTRAINT
    // ------------------------------------------------------------------------
    // (8) xi * ni + ci = 0 && xj * nj + cj = 0
    // (9) xi *(ni.x − nj.x) + xj * (ni.y − nj.y) + ( ci − cj ) = 0
    // ni = [ni.x, ni.y, 1.0]
    // nj = [nj.x, nj.y, 1.0]
    // (*) where xi and xj are any points on the ni and nj planes, respectively
    // ------------------------------------------------------------------------
    out += -1.0 * (p.x * ((x[an_x]) - (x[bn_x])) + p.y * ((x[an_y]) - (x[bn_y])) + (x[an_c] - x[bn_c]));

    if (!grad.empty()) {
        grad[an_x] += -1.0 * p.x;
        grad[an_y] += -1.0 * p.y;
        grad[an_c] += -1.0;
        grad[bn_x] += p.x;
        grad[bn_y] += p.y;
        grad[bn_c] += 1.0;
    }

    return out; // confunc
}

double Nlopt_Solve::constraint_function_isec0_ca_cb_pos(const std::vector<double> &x, std::vector<double> &grad, void *con_data)
{
    constraint_orthogonality_data_s *data = reinterpret_cast<constraint_orthogonality_data_s*>(con_data);
    intersection_data_s intersection_data = data->intersection_data;

    double out = 0.0; // confunc

    unsigned int an_x = intersection_data.polybezier_a_index * 3 + 0;
    unsigned int an_y = intersection_data.polybezier_a_index * 3 + 1;
    unsigned int an_c = intersection_data.polybezier_a_index * 3 + 2;
    unsigned int bn_x = intersection_data.polybezier_b_index * 3 + 0;
    unsigned int bn_y = intersection_data.polybezier_b_index * 3 + 1;
    unsigned int bn_c = intersection_data.polybezier_b_index * 3 + 2;

    // -----------------------------
    // CROSSHAIR-0 OFFSET CONSTRAINT
    // -----------------------------------------------
    // (9.1) Spline i and j at crosshair ij has depth of 0.0
    // (*) necessary for convergence
    // -----------------------------------------------
    out += x[an_c] - x[bn_c];

    if (!grad.empty()) {
        grad[an_x] += 0.0;
        grad[an_y] += 0.0;
        grad[an_c] += 1.0;
        grad[bn_x] += 0.0;
        grad[bn_y] += 0.0;
        grad[bn_c] += -1.0;
    }

    return out; // confunc
}

double Nlopt_Solve::constraint_function_isec0_ca_cb_neg(const std::vector<double> &x, std::vector<double> &grad, void *con_data)
{
    constraint_orthogonality_data_s *data = reinterpret_cast<constraint_orthogonality_data_s*>(con_data);
    intersection_data_s intersection_data = data->intersection_data;

    double out = 0.0; // confunc

    unsigned int an_x = intersection_data.polybezier_a_index * 3 + 0;
    unsigned int an_y = intersection_data.polybezier_a_index * 3 + 1;
    unsigned int an_c = intersection_data.polybezier_a_index * 3 + 2;
    unsigned int bn_x = intersection_data.polybezier_b_index * 3 + 0;
    unsigned int bn_y = intersection_data.polybezier_b_index * 3 + 1;
    unsigned int bn_c = intersection_data.polybezier_b_index * 3 + 2;

    // -----------------------------
    // CROSSHAIR-0 OFFSET CONSTRAINT
    // -----------------------------------------------
    // (9.1) Spline i and j at crosshair ij has depth of 0.0
    // (*) necessary for convergence
    // -----------------------------------------------
    out += x[bn_c] - x[an_c];

    if (!grad.empty()) {
        grad[an_x] += 0.0;
        grad[an_y] += 0.0;
        grad[an_c] += -1.0;
        grad[bn_x] += 0.0;
        grad[bn_y] += 0.0;
        grad[bn_c] += 1.0;
    }

    return out; // confunc
}
