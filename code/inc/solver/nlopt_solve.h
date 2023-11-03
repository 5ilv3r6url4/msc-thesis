#ifndef NLOPT_SOLVE_H
#define NLOPT_SOLVE_H

#include "model/polybezier.h"
#include "model/intersection.h"

#include "utils/math/util_math_defs.h"

#include "nlopt.hpp"

#include <vector>
#include <iomanip>

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QTextStream>

class Nlopt_Solve : public QObject
{
    Q_OBJECT

public slots:

signals:

public:

    typedef struct intersection_data_s {
        int polybezier_a_index;
        int polybezier_b_index;
        vector_t tangent_a;
        vector_t tangent_b;
        point_t position;

        intersection_data_s() {
            this->polybezier_a_index = -1;
            this->polybezier_b_index = -1;
            this->tangent_a = { 0.0, 0.0, 0.0 };
            this->tangent_b = { 0.0, 0.0, 0.0 };
            this->position = { 0.0, 0.0, 0.0 };
        }
        intersection_data_s(int polybezier_a_index,
                            int polybezier_b_index,
                            vector_t tangent_a,
                            vector_t tangent_b,
                            point_t position) {
            this->polybezier_a_index = polybezier_a_index;
            this->polybezier_b_index = polybezier_b_index;
            this->tangent_a = tangent_a;
            this->tangent_b = tangent_b;
            this->position = position;
        }
    } intersection_data_s;

    typedef struct optimization_data_s {
        bool silhouette_set;
        std::vector<double> initial_guess;
        std::vector<intersection_data_s> intersection_datas;

        optimization_data_s() {
            this->silhouette_set = false;
            this->initial_guess.resize(0);
            this->intersection_datas.resize(0);
        }
        optimization_data_s(std::vector<double> initial_guess,
                            std::vector<intersection_data_s> intersection_datas) {
            this->silhouette_set = false;
            this->initial_guess = initial_guess;
            this->intersection_datas = intersection_datas;
        }
    } optimization_data_s;

    typedef struct minimization_data_s {
        std::vector<intersection_data_s> intersection_datas;

        minimization_data_s() {
            this->intersection_datas.resize(0);
        }
        minimization_data_s(std::vector<intersection_data_s> intersection_datas) {
            this->intersection_datas = intersection_datas;
        }
    } minimization_data_s;

    typedef struct constraint_orthogonality_data_s {
        intersection_data_s intersection_data;

        constraint_orthogonality_data_s() {
            this->intersection_data = intersection_data_s();
        }
        constraint_orthogonality_data_s(intersection_data_s intersection_data) {
            this->intersection_data = intersection_data;
        }
    } constraint_orthogonality_data_s;

    typedef struct constraint_offset_data_s {
        intersection_data_s intersection_data_fixed;
        intersection_data_s intersection_data;

        constraint_offset_data_s() {
            this->intersection_data_fixed = intersection_data_s();
            this->intersection_data = intersection_data_s();
        }
        constraint_offset_data_s(intersection_data_s intersection_data_fixed, intersection_data_s intersection_data) {
            this->intersection_data_fixed = intersection_data_fixed;
            this->intersection_data = intersection_data;
        }
    } constraint_offset_data_s;

    static QPair<double, std::vector<double>> solve(optimization_data_s optimization_data);

    static double minimization_function(const std::vector<double> &x, std::vector<double> &grad, void *min_data);
    static double constraint_function_tangent_orthogonality_pos(const std::vector<double> &x, std::vector<double> &grad, void *con_data);
    static double constraint_function_tangent_orthogonality_neg(const std::vector<double> &x, std::vector<double> &grad, void *con_data);
    static double constraint_function_plane_orthogonality_pos(const std::vector<double> &x, std::vector<double> &grad, void *con_data);
    static double constraint_function_plane_orthogonality_neg(const std::vector<double> &x, std::vector<double> &grad, void *con_data);
    static double constraint_function_offset_pos(const std::vector<double> &x, std::vector<double> &grad, void *con_data);
    static double constraint_function_offset_neg(const std::vector<double> &x, std::vector<double> &grad, void *con_data);
    static double constraint_function_isec0_ca_cb_pos(const std::vector<double> &x, std::vector<double> &grad, void *con_data);
    static double constraint_function_isec0_ca_cb_neg(const std::vector<double> &x, std::vector<double> &grad, void *con_data);

    enum NLOPT_RETURN_CODE { NLOPT_SUCCESS          =  1,
                             NLOPT_STOPVAL_REACHED  =  2,
                             NLOPT_FTOL_REACHED     =  3,
                             NLOPT_XTOL_REACHED     =  4,
                             NLOPT_MAXEVAL_REACHED  =  5,
                             NLOPT_MAXTIME_REACHED  =  6,
                             NLOPT_FAILURE          = -1,
                             NLOPT_INVALID_ARGS     = -2,
                             NLOPT_OUT_OF_MEMORY    = -3,
                             NLOPT_ROUNDOFF_LIMITED = -4,
                             NLOPT_FORCED_STOP      = -5 };
};

#endif // NLOPT_SOLVE_H
