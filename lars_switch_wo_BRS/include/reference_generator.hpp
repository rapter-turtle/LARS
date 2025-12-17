#pragma once
#include <vector>
#include <cmath>
#include <algorithm>

struct RefPoint {
    double x;
    double y;
    double yaw;
    double u;
    double r;
};

class ReferenceGenerator
{
public:

    // Figure-eight trajectory
    static std::vector<RefPoint> generateFigureEight(
        double tfinal,
        double dt,
        double tx,
        double ty,
        double theta,
        double A = 80.0,
        double B = 80.0,
        double C = 60.0)
    {
        int N = int(tfinal / dt);
        std::vector<RefPoint> ref;
        ref.reserve(N);

        std::vector<double> x(N), y(N);

        // Position
        for (int i = 0; i < N; i++) {
            double t = i * dt;
            x[i] = A * std::sin(t / C);
            y[i] = B * std::sin(t / C) * std::cos(t / C);
        }

        // Velocity
        std::vector<double> dx(N), dy(N);
        for (int i = 1; i < N - 1; i++) {
            dx[i] = (x[i + 1] - x[i - 1]) / (2 * dt);
            dy[i] = (y[i + 1] - y[i - 1]) / (2 * dt);
        }
        dx[0] = dx[1];
        dy[0] = dy[1];
        dx[N - 1] = dx[N - 2];
        dy[N - 1] = dy[N - 2];

        // Heading
        std::vector<double> yaw(N);
        yaw[0] = std::atan2(dy[0], dx[0]);
        for (int i = 1; i < N; i++) {
            double raw = std::atan2(dy[i], dx[i]);
            yaw[i] = unwrapAngle(yaw[i - 1], raw);
        }

        // Rotational speed
        std::vector<double> r(N);
        for (int i = 1; i < N - 1; i++) {
            r[i] = (yaw[i + 1] - yaw[i - 1]) / (2 * dt);
        }
        r[0] = r[1];
        r[N - 1] = r[N - 2];

        // Output trajectory
        for (int i = 0; i < N; i++) {
            RefPoint pt;

            // Rotation
            double xr = x[i] * std::cos(theta) - y[i] * std::sin(theta) + tx;
            double yr = x[i] * std::sin(theta) + y[i] * std::cos(theta) + ty;

            pt.x = xr;
            pt.y = yr;
            pt.yaw = yaw[i] + theta;
            pt.u = std::sqrt(dx[i] * dx[i] + dy[i] * dy[i]);
            pt.r = r[i];

            ref.push_back(pt);
        }

        return ref;
    }


    // Linear trajectory (LARS)
    static std::vector<RefPoint> generateLinear(
        double tfinal,
        double dt,
        double tx,
        double ty,
        double theta,
        double A = 20.0)
    {
        int N = int(tfinal / dt);
        std::vector<RefPoint> ref(N);

        for (int i = 0; i < N; i++) {
            double t = i * dt;

            double x = A * t;
            double y = 0.0;

            double xr = x * std::cos(theta) - y * std::sin(theta) + tx;
            double yr = x * std::sin(theta) + y * std::cos(theta) + ty;

            RefPoint r;
            r.x = xr;
            r.y = yr;
            r.yaw = theta;
            r.u = A;
            r.r = 0.0;

            ref[i] = r;
        }

        return ref;
    }


    // Continuous-arc-length figure-eight (optional)
    static std::vector<RefPoint> generateFigureEightCon(
        double tfinal,
        double dt,
        double tx,
        double ty,
        double theta,
        double A = 80.0,
        double B = 80.0,
        double C = 30.0)
    {
        int N = int(tfinal / dt);

        std::vector<double> x(N), y(N);

        for (int i = 0; i < N; i++) {
            double t = i * dt;
            x[i] = A * std::sin(t / C);
            y[i] = B * std::sin(t / C) * std::cos(t / C);
        }

        // arc length derivative
        std::vector<double> dx(N), dy(N), ds(N);
        for (int i = 1; i < N - 1; i++) {
            dx[i] = (x[i + 1] - x[i - 1]) / 2.0;
            dy[i] = (y[i + 1] - y[i - 1]) / 2.0;
            ds[i] = std::sqrt(dx[i] * dx[i] + dy[i] * dy[i]);
        }

        // cumulative arc length
        std::vector<double> s(N);
        s[0] = 0;
        for (int i = 1; i < N; i++)
            s[i] = s[i - 1] + ds[i];

        // uniform s
        std::vector<double> s_uniform(N);
        double s_final = s.back();
        for (int i = 0; i < N; i++)
            s_uniform[i] = (double)i / (N - 1) * s_final;

        // map s → t using linear interpolation
        std::vector<double> t_uniform(N);
        for (int i = 0; i < N; i++)
            t_uniform[i] = interp1D(s, t_uniform[i], s_uniform);

        // recompute uniform trajectory
        std::vector<RefPoint> ref(N);
        for (int i = 0; i < N; i++)
        {
            double t = t_uniform[i];

            double xx = A * std::sin(t / C);
            double yy = B * std::sin(t / C) * std::cos(t / C);

            double dxu = (A / C) * std::cos(t / C);
            double dyu = (B / C) * (std::cos(t / C)*std::cos(t / C) -
                                    std::sin(t / C)*std::sin(t / C));

            double yaw = std::atan2(dyu, dxu);

            double vel = std::sqrt(dxu*dxu + dyu*dyu);

            // rotation
            double xr = xx * std::cos(theta) - yy * std::sin(theta) + tx;
            double yr = xx * std::sin(theta) + yy * std::cos(theta) + ty;

            RefPoint rp;
            rp.x = xr;
            rp.y = yr;
            rp.yaw = yaw + theta;
            rp.u = vel;
            rp.r = 0.0;

            ref[i] = rp;
        }

        return ref;
    }


private:

    // unwrap logic (same as numpy.unwrap)
    static double unwrapAngle(double prev, double now)
    {
        double a = now - prev;

        if (a > M_PI) now -= 2 * M_PI;
        else if (a < -M_PI) now += 2 * M_PI;

        return now;
    }


    // simple 1D interpolation function (linear)
    static double interp1D(const std::vector<double> &x, double xq,
                           const std::vector<double> &y)
    {
        auto it = std::lower_bound(x.begin(), x.end(), xq);
        int idx = std::max(1, int(it - x.begin()));

        double x1 = x[idx - 1], x2 = x[idx];
        double y1 = y[idx - 1], y2 = y[idx];

        double ratio = (xq - x1) / (x2 - x1);
        return y1 + ratio * (y2 - y1);
    }
};
