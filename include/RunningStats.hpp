#ifndef __RUNNINGSTATS_HPP__
#define __RUNNINGSTATS_HPP__

/**
 * Shamelessly stolen from http://www.johndcook.com/blog/standard_deviation/
 */
class RunningStats {
private:
    int m_n;
    double m_oldM, m_newM, m_oldS, m_newS;

public:

    RunningStats() : m_n(0) {
    }

    void clear() {
        m_n = 0;
    }

    void push(double x) {
        m_n++;

        // See Knuth TAOCP vol 2, 3rd edition, page 232
        if (m_n == 1) {
            m_oldM = m_newM = x;
            m_oldS = 0.0;
        } else {
            m_newM = m_oldM + (x - m_oldM) / m_n;
            m_newS = m_oldS + (x - m_oldM)*(x - m_newM);

            // set up for next iteration
            m_oldM = m_newM;
            m_oldS = m_newS;
        }
    }

    int count() const {
        return m_n;
    }

    double mean() const {
        return (m_n > 0) ? m_newM : 0.0;
    }

    double variance() const {
        return ( (m_n > 1) ? m_newS / (m_n - 1) : 0.0);
    }

    double stdDev() const {
        return sqrt(variance());
    }

};

#endif