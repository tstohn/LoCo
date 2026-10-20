#include <unordered_set>
#include <random>
#include <algorithm> 
#include <limits>
#include <set>
//#include <execution> // C++17 can improve runtime of std::sort with std::sort(std::execution::par, begin, end, ...), but we disabled it for R weaknesses to handle tbb
#include <atomic>
#include <chrono>
#include <iomanip>
#include <map>
#include <functional>
#include <tuple>

#include "Neighborhood.hpp"

namespace neighborhoodCalculations
{

    // random engine for the neighbourhood construction: seed 0 = a random seed (not reproducible)
    inline std::mt19937 make_rng(const unsigned int seed)
    {
        return std::mt19937(seed == 0 ? std::random_device{}() : seed);
    }

    std::vector<int> get_random_elements(unsigned int numbers, unsigned int maxNum, const unsigned int seed = 0) 
    {
        if (numbers > maxNum)
            throw std::invalid_argument("Creating more neighbourhoods that cells is not allowed!");

        std::vector<int> pool(maxNum);
        std::iota(pool.begin(), pool.end(), 0);

        std::mt19937 gen = make_rng(seed);

        std::shuffle(pool.begin(), pool.end(), gen);

        return std::vector<int>(pool.begin(), pool.begin() + numbers);
    }

    // per-cell status of the neighbourhood construction (kmeans)
    enum KmediansCellStatus : char { KM_PLACED = 0, KM_NO_ROOM = 1, KM_DISSOLVED = 2, KM_DROPPED = 3 };
    struct KmediansResult
    {
        std::vector<std::pair<int, std::vector<int>>> groups; // anchor (first member) + all members
        std::vector<char> status;                              // KmediansCellStatus per cell
    };

    // Disjoint neighbourhoods by a size-constrained k-medians clustering (Manhattan distance) in the cell state space
    // (= the nodes of scData: the cellStateGeneFile features, or all features). Neighbourhood sizes may vary by
    // `slack` around s (sMin = ceil(s(1-slack)) ... sMax = floor(s(1+slack)), e.g. 95-105 for s = 100, slack = 0.05).
    // Not every cell has to be in a neighbourhood.
    //   1.) k-means++ seeding: first centre a random cell, every further centre a random cell drawn with probability
    //       proportional to its distance to the nearest centre so far
    //   2.) k-medians iterations (at most maxIterations, until at most 0.1 % of the cells change):
    //         a.) ALL cells are assigned from scratch, cells with the largest regret (distance to the 2nd nearest minus
    //             distance to the nearest centre) first: each cell goes to the nearest of its R nearest centres that has
    //             fewer than s cells; cells that did not fit then go to the nearest of their R nearest centres with fewer
    //             than sMax cells (the slack is a tolerance, not the target: otherwise all neighbourhoods fill up to sMax
    //             and ~slack of them starve); if all R are full, the cell is not placed in this round (it is considered
    //             again in the next round). Neighbourhoods with fewer than sMin cells then take their nearest not-placed cells
    //             (even outside these cells' R nearest) up to sMin, so that no centre starves. From the 3rd iteration on,
    //             neighbourhoods that still have fewer than sMin cells are starving (per iteration only the smallest of
    //             neighbouring ones): their centre is re-seeded at the worst served cell (furthest from its centre), so the
    //             other centres of a small population take over its cells and the capacity for all cells is kept; in the
    //             last 5 iterations starving neighbourhoods are dissolved instead
    //         b.) centre = coordinate-wise median of its placed cells
    //   3.) final neighbourhoods:
    //         a.) neighbourhoods with fewer than sMin cells are dissolved (they do not exist), their cells are not placed
    //         b.) not-placed cells join the nearest of their R nearest neighbourhoods with room (closest cells first)
    //         c.) median update
    //         d.) swap phase for every placed cell whose own neighbourhood is NOT among its R nearest (worst first):
    //             move it into the nearest of its R nearest neighbourhoods with room (if its own keeps >= sMin), else swap
    //             it with the member of one of its R nearest neighbourhoods that fits best into its own neighbourhood
    //             (if the total distance decreases), else leave it unplaced (if its own keeps >= sMin), else it stays
    //             ("unresolved")
    //         e.) not-placed cells join the nearest of their R nearest neighbourhoods with room
    //         f.) final median update, anchor = member closest to the median (first in the member list)
    // Every final neighbourhood has sMin..sMax cells. Statistics are printed.
    KmediansResult kmedians_with_repair(const std::shared_ptr<const GraphData>& scData, const unsigned int kIn,
                                        const unsigned int s, const unsigned int Rin, const double slack,
                                        const unsigned int seed, const int threads, const int maxIterations = 30)
    {
        const auto t0 = std::chrono::steady_clock::now();
        const size_t n = scData->number_of_nodes();
        const size_t k = kIn;
        if (k < 1 || k > n)
            throw std::invalid_argument("kmeans sampling: number of neighbourhoods must be between 1 and the number of cells!");
        const size_t dim = scData->get_node_at(0)->dimensions();
        const size_t sMax = static_cast<size_t>(std::floor(s * (1.0 + slack) + 1e-9));
        const size_t sMin = std::max<size_t>(1, static_cast<size_t>(std::ceil(s * (1.0 - slack) - 1e-9)));
        const size_t R = std::max<size_t>(1, std::min<size_t>(Rin, k));

        // flat copy of the cell state coordinates
        std::vector<double> X(n * dim);
        for (size_t i = 0; i < n; ++i)
        {
            const nodePtr node = scData->get_node_at(i);
            for (size_t d = 0; d < dim; ++d) X[i * dim + d] = node->value_at(static_cast<int>(d));
        }
        auto manhattan = [&](const double* a, const double* b)
        {
            double dist = 0.0;
            for (size_t d = 0; d < dim; ++d) dist += std::abs(a[d] - b[d]);
            return dist;
        };

        // body(begin, end) over [0, count), split over the threads of one pool (reused for all loops)
        ThreadPool pool(std::max(1, threads));
        auto parallel_for = [&](size_t count, const std::function<void(size_t, size_t)>& body)
        {
            const size_t chunks = std::min<size_t>(std::max(1, threads), std::max<size_t>(1, count / 512));
            if (chunks <= 1) { body(0, count); return; }
            const size_t step = (count + chunks - 1) / chunks;
            for (size_t b = 0; b < count; b += step)
            {
                const size_t e = std::min(count, b + step);
                pool.enqueue([&body, b, e]() { body(b, e); });
            }
            pool.wait_for_tasks();
        };

        std::mt19937 gen = make_rng(seed);

        // ---- 1.) k-means++ seeding (weights = Manhattan distance to the nearest centre) ----
        std::vector<double> centres(k * dim);
        {
            std::vector<double> minDist(n, std::numeric_limits<double>::max());
            size_t next = std::uniform_int_distribution<size_t>(0, n - 1)(gen);
            for (size_t c = 0; c < k; ++c)
            {
                std::copy(&X[next * dim], &X[next * dim] + dim, &centres[c * dim]);
                double sum = 0.0;
                for (size_t i = 0; i < n; ++i)
                {
                    minDist[i] = std::min(minDist[i], manhattan(&X[i * dim], &centres[c * dim]));
                    sum += minDist[i];
                }
                if (c + 1 == k) break;
                if (sum <= 0.0) { next = std::uniform_int_distribution<size_t>(0, n - 1)(gen); continue; }
                double target = std::uniform_real_distribution<double>(0.0, sum)(gen);
                next = n - 1;
                for (size_t i = 0; i < n; ++i)
                {
                    target -= minDist[i];
                    if (target <= 0.0) { next = i; break; }
                }
            }
        }

        std::vector<char> alive(k, 1);
        std::vector<int> aliveIdx;
        std::vector<int> assign(n, -1);
        std::vector<int> prefs(n * R, -1);     // R nearest alive centres of every cell (nearest first)
        std::vector<double> dNearest(n), regret(n), dOwn(n);
        std::vector<int> rankOwn(n, 0);        // rank of the own centre among all alive centres (1 = nearest), 0 = not placed
        size_t Ra = R;

        // distances of every cell to all alive centres: R nearest, regret, rank of the own centre
        auto compute_pass = [&]()
        {
            aliveIdx.clear();
            for (size_t c = 0; c < k; ++c) if (alive[c]) aliveIdx.push_back(static_cast<int>(c));
            Ra = std::min(R, aliveIdx.size());
            const size_t m = std::min(aliveIdx.size(), std::max<size_t>(Ra, 2));
            parallel_for(n, [&](size_t b, size_t e)
            {
                std::vector<double> d(aliveIdx.size());
                std::vector<int> ord(aliveIdx.size());
                for (size_t i = b; i < e; ++i)
                {
                    const int own = assign[i];
                    double dOwnI = 0.0;
                    for (size_t q = 0; q < aliveIdx.size(); ++q)
                    {
                        d[q] = manhattan(&X[i * dim], &centres[static_cast<size_t>(aliveIdx[q]) * dim]);
                        if (aliveIdx[q] == own) dOwnI = d[q];
                    }
                    std::iota(ord.begin(), ord.end(), 0);
                    std::partial_sort(ord.begin(), ord.begin() + m, ord.end(),
                                      [&](int x, int y) { return d[x] < d[y] || (d[x] == d[y] && x < y); });
                    for (size_t p = 0; p < Ra; ++p) prefs[i * R + p] = aliveIdx[ord[p]];
                    dNearest[i] = d[ord[0]];
                    regret[i] = (aliveIdx.size() > 1) ? d[ord[1]] - d[ord[0]] : 0.0;
                    if (own >= 0)
                    {
                        int r = 1;
                        for (size_t q = 0; q < aliveIdx.size(); ++q) if (d[q] < dOwnI) ++r;
                        rankOwn[i] = r;
                        dOwn[i] = dOwnI;
                    }
                    else { rankOwn[i] = 0; dOwn[i] = 0.0; }
                }
            });
        };

        // centre = coordinate-wise median of its placed cells (the minimiser of the sum of Manhattan distances)
        auto update_medians = [&]()
        {
            std::vector<std::vector<int>> mem(k);
            for (size_t i = 0; i < n; ++i) if (assign[i] >= 0) mem[assign[i]].push_back(static_cast<int>(i));
            std::vector<double> values;
            for (size_t c = 0; c < k; ++c)
            {
                if (!alive[c] || mem[c].empty()) continue;
                values.resize(mem[c].size());
                for (size_t d = 0; d < dim; ++d)
                {
                    for (size_t m = 0; m < mem[c].size(); ++m) values[m] = X[static_cast<size_t>(mem[c][m]) * dim + d];
                    const size_t mid = values.size() / 2;
                    std::nth_element(values.begin(), values.begin() + mid, values.end());
                    double median = values[mid];
                    if (values.size() % 2 == 0) median = 0.5 * (median + *std::max_element(values.begin(), values.begin() + mid));
                    centres[c * dim + d] = median;
                }
            }
        };

        auto dist_to_centre = [&](int i, int c) { return manhattan(&X[static_cast<size_t>(i) * dim], &centres[static_cast<size_t>(c) * dim]); };

        std::vector<size_t> sizes(k, 0);
        auto count_sizes = [&]()
        {
            std::fill(sizes.begin(), sizes.end(), 0);
            for (size_t i = 0; i < n; ++i) if (assign[i] >= 0) ++sizes[assign[i]];
        };

        // ---- 2.) k-medians iterations ----
        std::vector<int> order(n), newAssign(n);
        std::vector<size_t> notPlacedPerIteration;
        int iterations = 0;
        bool converged = false;
        size_t lastChanged = n;
        size_t dissolvedDuringIterations = 0;
        size_t reseeded = 0;
        std::vector<char> reseedUsed(n, 0);
        for (int it = 0; it < maxIterations; ++it)
        {
            ++iterations;
            compute_pass();
            std::iota(order.begin(), order.end(), 0);
            std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return regret[a] > regret[b]; });
            std::fill(sizes.begin(), sizes.end(), 0);
            std::fill(newAssign.begin(), newAssign.end(), -1);
            // pass 1: fill up to the normal size s; pass 2: cells that did not fit may use the slack (up to sMax)
            for (const size_t cap : {static_cast<size_t>(s), sMax})
            {
                for (const int i : order)
                {
                    if (newAssign[i] >= 0) continue;
                    for (size_t p = 0; p < Ra; ++p)
                    {
                        const int c = prefs[static_cast<size_t>(i) * R + p];
                        if (sizes[c] < cap) { newAssign[i] = c; ++sizes[c]; break; }
                    }
                }
            }
            size_t notPlaced = 0;
            for (size_t i = 0; i < n; ++i) if (newAssign[i] < 0) ++notPlaced;
            // under-full neighbourhoods (< sMin) take their nearest not-placed cells (even outside the cells' R nearest)
            // so that no centre starves: its median then moves towards the cells nobody has room for. Cells placed far
            // away this way are repaired in the final swap phase.
            {
                std::vector<int> underFull;
                for (const int c : aliveIdx) if (sizes[c] < sMin) underFull.push_back(c);
                if (!underFull.empty() && notPlaced > 0)
                {
                    std::vector<std::tuple<double, int, int>> candidates;   // (distance, cell, centre)
                    for (size_t i = 0; i < n; ++i)
                    {
                        if (newAssign[i] >= 0) continue;
                        for (const int c : underFull) candidates.emplace_back(dist_to_centre(static_cast<int>(i), c), static_cast<int>(i), c);
                    }
                    std::sort(candidates.begin(), candidates.end());
                    for (const auto& [d, i, c] : candidates)
                    {
                        if (newAssign[i] >= 0 || sizes[c] >= sMin) continue;
                        newAssign[i] = c;
                        ++sizes[c];
                        --notPlaced;
                    }
                }
            }
            // starving neighbourhoods (still < sMin, no free cells to fill them) are dissolved during the iterations
            // (from the 3rd on, the first ones only follow the seeding): e.g. 4 centres sharing a population of 315 cells
            // (~79 each) would otherwise all be dissolved at the end and the population lost. Per iteration only the
            // smallest of neighbouring starving neighbourhoods is dissolved (its R nearest centres are skipped), so the
            // remaining centres of such a population take over its cells in the next iteration.
            size_t dissolvedNow = 0;
            if (it >= 2)
            {
                std::vector<int> starving;
                for (const int c : aliveIdx) if (sizes[c] < sMin) starving.push_back(c);
                std::stable_sort(starving.begin(), starving.end(), [&](int x, int y) { return sizes[x] < sizes[y]; });
                std::vector<char> skip(k, 0);
                for (const int c : starving)
                {
                    if (skip[c]) continue;
                    // its R nearest alive centres are not dissolved in this iteration
                    std::vector<std::pair<double, int>> near;
                    for (const int o : aliveIdx)
                    {
                        if (o == c || !alive[o]) continue;
                        near.emplace_back(manhattan(&centres[static_cast<size_t>(c) * dim], &centres[static_cast<size_t>(o) * dim]), o);
                    }
                    const size_t m = std::min(R, near.size());
                    std::partial_sort(near.begin(), near.begin() + m, near.end());
                    for (size_t q = 0; q < m; ++q) skip[near[q].second] = 1;
                    for (size_t i = 0; i < n; ++i) if (newAssign[i] == c) newAssign[i] = -1;
                    if (it < maxIterations - 5)
                    {
                        // RE-SEED instead of losing the centre (keeps the capacity for all cells): move it to the worst
                        // served cell, i.e. the cell furthest from its own centre (a not-placed cell: furthest from its
                        // nearest centre) that is not already used as a re-seed point in this iteration
                        int best = -1;
                        double bestD = -1.0;
                        for (size_t i = 0; i < n; ++i)
                        {
                            if (reseedUsed[i]) continue;
                            const double d = (newAssign[i] >= 0) ? dist_to_centre(static_cast<int>(i), newAssign[i]) : dNearest[i];
                            if (d > bestD) { bestD = d; best = static_cast<int>(i); }
                        }
                        if (best >= 0)
                        {
                            reseedUsed[best] = 1;
                            std::copy(&X[static_cast<size_t>(best) * dim], &X[static_cast<size_t>(best) * dim] + dim, &centres[static_cast<size_t>(c) * dim]);
                            ++reseeded;
                            ++dissolvedNow;   // the assignment changed
                            continue;
                        }
                    }
                    alive[c] = 0;
                    ++dissolvedNow;
                    ++dissolvedDuringIterations;
                }
            }
            notPlacedPerIteration.push_back(notPlaced);
            // stable: at most 0.1 % of the cells changed (a few cells can flip between two equally close centres forever)
            lastChanged = 0;
            for (size_t i = 0; i < n; ++i) if (newAssign[i] != assign[i]) ++lastChanged;
            converged = (dissolvedNow == 0 && lastChanged <= n / 1000);
            assign.swap(newAssign);
            if (converged) break;
            {
                // re-seeded centres have no cells yet: keep their new position (update_medians skips empty centres)
                update_medians();
            }
        }

        // ---- 3.) final neighbourhoods ----
        std::vector<char> status(n, KM_PLACED);
        for (size_t i = 0; i < n; ++i) if (assign[i] < 0) status[i] = KM_NO_ROOM;

        // a.) dissolve neighbourhoods with fewer than sMin cells
        count_sizes();
        size_t dissolved = 0;
        for (size_t c = 0; c < k; ++c)
        {
            if (alive[c] && sizes[c] < sMin) { alive[c] = 0; ++dissolved; }
        }
        for (size_t i = 0; i < n; ++i)
        {
            if (assign[i] >= 0 && !alive[assign[i]]) { assign[i] = -1; status[i] = KM_DISSOLVED; }
        }
        if (dissolved + dissolvedDuringIterations == k)
            throw std::invalid_argument("kmeans sampling: no neighbourhood reaches the minimum size (neighbourhoodSize too large?)");
        count_sizes();

        // place not-placed cells (closest to a centre first) into the nearest of their R nearest neighbourhoods with room
        auto readmit = [&]()
        {
            std::vector<int> unplaced;
            for (size_t i = 0; i < n; ++i) if (assign[i] < 0) unplaced.push_back(static_cast<int>(i));
            std::stable_sort(unplaced.begin(), unplaced.end(), [&](int a, int b) { return dNearest[a] < dNearest[b]; });
            size_t placed = 0;
            for (const int i : unplaced)
            {
                for (size_t p = 0; p < Ra; ++p)
                {
                    const int c = prefs[static_cast<size_t>(i) * R + p];
                    if (sizes[c] < sMax) { assign[i] = c; ++sizes[c]; status[i] = KM_PLACED; ++placed; break; }
                }
            }
            return placed;
        };

        // b.) + c.)
        compute_pass();
        const size_t readmittedAfterDissolve = readmit();
        update_medians();

        // d.) swap phase for the placed cells whose own neighbourhood is not among their R nearest
        compute_pass();
        std::vector<std::vector<int>> members(k);
        std::vector<size_t> posInMembers(n, 0);
        for (size_t i = 0; i < n; ++i)
        {
            if (assign[i] < 0) continue;
            posInMembers[i] = members[assign[i]].size();
            members[assign[i]].push_back(static_cast<int>(i));
        }
        auto relocate = [&](int i, int to)   // to = -1: not placed
        {
            const int from = assign[i];
            if (from >= 0)
            {
                std::vector<int>& mf = members[from];
                const int last = mf.back();
                mf[posInMembers[i]] = last;
                posInMembers[last] = posInMembers[i];
                mf.pop_back();
                --sizes[from];
            }
            if (to >= 0)
            {
                posInMembers[i] = members[to].size();
                members[to].push_back(i);
                ++sizes[to];
            }
            assign[i] = to;
        };
        auto dist_to = [&](int i, int c) { return manhattan(&X[static_cast<size_t>(i) * dim], &centres[static_cast<size_t>(c) * dim]); };

        std::vector<int> badly;
        for (size_t i = 0; i < n; ++i) if (assign[i] >= 0 && rankOwn[i] > static_cast<int>(R)) badly.push_back(static_cast<int>(i));
        std::stable_sort(badly.begin(), badly.end(), [&](int a, int b)
        {
            if (rankOwn[a] != rankOwn[b]) return rankOwn[a] > rankOwn[b];
            return (dOwn[a] - dNearest[a]) > (dOwn[b] - dNearest[b]);
        });
        const size_t badlyBefore = badly.size();
        size_t moved = 0, swapped = 0, dropped = 0, unresolved = 0;
        for (const int a : badly)
        {
            const int A = assign[a];
            if (A < 0) continue;
            const double dA = dist_to(a, A);
            // still outside its R nearest? (it may have been moved as a swap partner meanwhile)
            int rank = 1;
            for (const int c : aliveIdx) if (dist_to(a, c) < dA) ++rank;
            if (rank <= static_cast<int>(R)) continue;

            bool done = false;
            // 1.) move into the nearest of its R nearest neighbourhoods with room
            if (sizes[A] > sMin)
            {
                for (size_t p = 0; p < Ra && !done; ++p)
                {
                    const int B = prefs[static_cast<size_t>(a) * R + p];
                    if (sizes[B] < sMax) { relocate(a, B); ++moved; done = true; }
                }
            }
            // 2.) swap with the member of one of its R nearest neighbourhoods that fits best into A
            for (size_t p = 0; p < Ra && !done; ++p)
            {
                const int B = prefs[static_cast<size_t>(a) * R + p];
                const double gainA = dA - dist_to(a, B);
                int bestB = -1;
                double bestGainB = -std::numeric_limits<double>::infinity();
                for (const int b : members[B])
                {
                    const double g = dist_to(b, B) - dist_to(b, A);
                    if (g > bestGainB) { bestGainB = g; bestB = b; }
                }
                if (bestB >= 0 && gainA + bestGainB > 1e-12 * std::max(1.0, dA))
                {
                    relocate(a, B);
                    relocate(bestB, A);
                    ++swapped;
                    done = true;
                }
            }
            // 3.) leave it unplaced, 4.) otherwise it stays
            if (!done)
            {
                if (sizes[A] > sMin) { relocate(a, -1); status[a] = KM_DROPPED; ++dropped; }
                else ++unresolved;
            }
        }

        // e.) not-placed cells into their R nearest neighbourhoods with room (centres unchanged since the last pass)
        const size_t readmittedAfterSwaps = readmit();

        // f.) final medians, statistics and anchors
        update_medians();
        compute_pass();
        count_sizes();

        KmediansResult result;
        result.status = status;
        std::vector<size_t> finalSizes;
        std::vector<double> radii;
        std::vector<std::vector<int>> finalMembers(k);
        for (size_t i = 0; i < n; ++i) if (assign[i] >= 0) finalMembers[assign[i]].push_back(static_cast<int>(i));
        for (size_t c = 0; c < k; ++c)
        {
            if (!alive[c]) continue;
            std::vector<int>& mem = finalMembers[c];
            if (mem.size() < sMin || mem.size() > sMax)
                throw std::runtime_error("kmeans sampling: internal error, a neighbourhood has " + std::to_string(mem.size()) + " cells");
            size_t anchorPos = 0;
            double radius = 0.0;
            for (size_t m = 0; m < mem.size(); ++m)
            {
                radius += dOwn[mem[m]];
                if (dOwn[mem[m]] < dOwn[mem[anchorPos]]) anchorPos = m;
            }
            radii.push_back(radius / static_cast<double>(mem.size()));
            finalSizes.push_back(mem.size());
            std::swap(mem[0], mem[anchorPos]);
            result.groups.emplace_back(mem[0], mem);
        }

        // ---- short summary ----
        size_t placedCells = 0, rank1 = 0, rank2toR = 0;
        for (size_t i = 0; i < n; ++i)
        {
            if (assign[i] < 0) continue;
            ++placedCells;
            if (rankOwn[i] == 1) ++rank1;
            else if (rankOwn[i] <= static_cast<int>(R)) ++rank2toR;
        }
        // counters of the individual steps (not printed)
        (void)iterations; (void)converged; (void)lastChanged; (void)reseeded; (void)dissolvedDuringIterations;
        (void)readmittedAfterDissolve; (void)readmittedAfterSwaps; (void)badlyBefore; (void)moved; (void)swapped;
        (void)dropped; (void)unresolved; (void)t0;
        auto pct = [](double x, double total) { return total > 0 ? 100.0 * x / total : 0.0; };
        LOCO_OUT << std::fixed << std::setprecision(1);
        LOCO_OUT << "\tBuilding neighbourhoods with kmeans: " << result.groups.size() << " neighbourhoods of " << sMin << "-" << sMax << " cells\n";
        LOCO_OUT << "\t" << pct(placedCells, n) << " % of cells assigned to a neighbourhood\n";
        LOCO_OUT << "\tof these: " << pct(rank1, placedCells) << " % in their closest neighbourhood, " << pct(rank2toR, placedCells)
                 << " % in their 2nd-" << R << "th closest\n";
        LOCO_OUT.unsetf(std::ios_base::floatfield);
        LOCO_OUT << std::setprecision(6);
        return result;
    }

    //return LOWEST values first: we want the values with LWOEST laplacian first
    bool sort_corr(const std::pair<int, int>& a, const std::pair<int, int>& b,
                   std::unordered_map< std::pair<int, int>, const double, pair_hash> correlationLaplacian)
    {
        return(correlationLaplacian.at(a) < correlationLaplacian.at(b));
    }

    std::vector<std::vector<int>> generate_subsets(
        const std::vector<int>& set,
        int minSize)
    {
        std::vector<std::vector<int>> subsets;

        int n = set.size();
        int maxMask = 1 << n;

        for(int mask = 0; mask < maxMask; ++mask)
        {
            int count = __builtin_popcount(mask);
            if(count < minSize) continue;

            std::vector<int> subset;
            subset.reserve(count);

            for(int i = 0; i < n; ++i)
            {
                if(mask & (1 << i))
                    subset.push_back(set[i]);
            }

            std::sort(subset.begin(), subset.end());
            subsets.push_back(subset);
        }

        return subsets;
    }


    bool contains_subset(const std::vector<int>& superset,
                        const std::vector<int>& subset)
    {
        return std::includes(
            superset.begin(), superset.end(),
            subset.begin(), subset.end());
    }
}

//it filters the x top pairs for CORRELATIONS and SLOPES
std::vector<std::pair<int, int>> Neighborhood::filter_best_pairs(size_t numberGenes)
{

    std::vector<std::pair<int, int>> correlations;
    for(auto& correlation : correlationLaplacian)
    {
        correlations.push_back(correlation.first);
    }
    //if the correlation filter is bigger than actual found correlation pairs reduce correlation filter to prevent segfault
    if(numberGenes > correlations.size() || numberGenes==0){numberGenes = correlations.size();}

    std::sort(correlations.begin(), correlations.end(),
    [&](const std::pair<int,int>& a, const std::pair<int,int>& b)
    {
        return neighborhoodCalculations::sort_corr(a, b, correlationLaplacian);
    });
    std::vector<std::pair<int, int>> filteredCorrelations(correlations.begin(), correlations.begin() + numberGenes);

    std::sort(correlations.begin(), correlations.end(),
    [&](const std::pair<int,int>& a, const std::pair<int,int>& b)
    {
        return neighborhoodCalculations::sort_corr(a, b, slopeLaplacian);
    });
    std::vector<std::pair<int, int>> filteredCorrelations_2(correlations.begin(), correlations.begin() + numberGenes);

    //insert all new correlatiopns from filteredCorrelations_2 into first vector
    for (const auto& correlation : filteredCorrelations_2) {
        // Check if element is not yet in vector b
        if (std::find(filteredCorrelations.begin(), filteredCorrelations.end(), correlation) == filteredCorrelations.end()) 
        {
            // Insert element into vector b
            filteredCorrelations.push_back(correlation);
        }
    }

    return(filteredCorrelations);
}


// old function to write all the permuted laplacian scores: not logner done
void Neighborhood::write_shuffled_laplacians(const std::string& outFile, const std::string& prefix)
{
    std::ofstream outputFile;
    std::string folder_prefix = outFile + "/" + prefix;

    // CORRELATION LAPLACIANS
    std::string outputCorr = folder_prefix + "_shuffledCorr.tsv" ;
    std::remove(outputCorr.c_str());
    outputFile.open(outputCorr, std::ofstream::app);

    // collect all localPairs
    std::vector<std::pair<int,int>> localPairs;
    for(const auto& it : shuffledCorrLaplacians)
    {
        localPairs.push_back(it.first);
    }

    if(localPairs.empty())
    {
        outputFile.close();
        return;
    }

    int permutations = shuffledCorrLaplacians.at(localPairs[0]).size();

    // HEADER
    outputFile << "Permutation";
    for(const std::pair<int,int>& pair : localPairs)
    {
        int name_a_idx = pair.first;
        int name_b_idx = pair.second;
        if(!corrStateGenes.empty())
        {
            name_a_idx = corrStateGenes.at(name_a_idx);
            name_b_idx = corrStateGenes.at(name_b_idx);
        }
        outputFile << "\t"
                   << inputDataOrigional.geneNames.at(name_a_idx)
                   << "_"
                   << inputDataOrigional.geneNames.at(name_b_idx);
    }
    outputFile << "\n";

    // rows = permutations
    for(int p = 0; p < permutations; ++p)
    {
        outputFile << p;

        for(const std::pair<int,int>& pair : localPairs)
        {
            outputFile << "\t" << shuffledCorrLaplacians.at(pair).at(p);
        }

        outputFile << "\n";
    }

    outputFile.close();

}

void Neighborhood::write_results_to_file(const std::string& outFile, const std::string& prefix, bool calcSets)
{
    std::ofstream outputFile;   

    /////////////////////////////
    // 1.) LAPLCACIAN SCORES FILE
    /////////////////////////////
    // Create an array of indices [0, 1, 2, ..., N-1]
    std::vector<size_t> indices(laplacianScores.pairNames.size());
    std::iota(indices.begin(), indices.end(), 0);

    // Sort the indices based on the Laplacian values (L) in ascending order
    std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
        
        bool isNanA = std::isnan(laplacianScores.L[a]);
        bool isNanB = std::isnan(laplacianScores.L[b]);
        
        // Safety: Push NaNs to the very bottom of the file
        if (isNanA && isNanB) return false; 
        if (isNanA) return false;
        if (isNanB) return true;
        
        // Ascending order: smallest L goes first
        return laplacianScores.L[a] < laplacianScores.L[b];
    });
    //1FILE with annotations: corrpair - LaplacianCorrscore - LaplacianSlopeScore - origionalCliqueItIsFrom
    std::string outputLaplace = outFile + "/" + prefix + "_laplacian.tsv";
    std::remove(outputLaplace.c_str());
    outputFile.open (outputLaplace, std::ofstream::app);
    //write HEADER
    outputFile << "FeaturePair" << "\t" << "CorrelationScore" << "\t" << "p_CorrelationScore"; 
    if(calcSets){ outputFile << "\t" << "Feauture Set";} 
    outputFile << "\n";
    //write LINES
    for(size_t i = 0; i < indices.size(); ++i)
    {
        size_t idx = indices[i]; // Map loop iteration to the correct sorted row

        // write PAIR << Laplacian << p_value
        outputFile << laplacianScores.pairNames.at(idx) << "\t" << 
        std::to_string(laplacianScores.L[idx]) << "\t" << 
        std::to_string(laplacianScores.p_values[idx]);

        if(calcSets)
        {
            outputFile << "\t" << neighbourhoodCorrs.featureSetString[idx];
        }

        outputFile << "\n";
    }
    outputFile.close();

    /////////////////////////////
    // CORRELATIONS:  for NeighborhoodID * CorrpAir: value == corr
    /////////////////////////////
    //CORRELATIONS FILE
    std::string outputCorr = outFile + "/" + prefix + "_correlations.tsv";
    std::remove(outputCorr.c_str());
    outputFile.open (outputCorr, std::ofstream::app);
    //write HEADER: first col has all featurePair names and then every col-header is a neighbourhood (inverse to origional version)
    outputFile << "FeaturePair";
    for(const auto& neighbourhoddPtr : neighbourhoodCorrs.nodePtrList)
    {
        outputFile << "\t" << neighbourhoddPtr->get_name();
    }
    outputFile << "\n";

    for(int featurePairIdx =0; featurePairIdx < neighbourhoodCorrs.pairNames.size(); ++featurePairIdx) 
    {

        std::string featurePair = neighbourhoodCorrs.pairNames.at(featurePairIdx);
        //get all correlation values for all neighbourhood for this featurePair
        const double* pairCorrs = neighbourhoodCorrs.pairToListOfAllNCorrs.get_row(featurePairIdx);

        outputFile << featurePair;
        for(int Nidx = 0; Nidx < neighbourhoodNum; ++Nidx)
        {
            outputFile << "\t" << pairCorrs[Nidx];
        }
        outputFile << "\n";
    }
    outputFile << "\n";
    outputFile.close();
    
    /////////////////////////////
    // 3.)  File with annotations: NeighborhoodID - centralNodeDimensions[start ... end]
    /////////////////////////////
    std::string outputNeighborhoodCoordinates = outFile + "/" + prefix + "_coord.tsv";
    std::remove(outputNeighborhoodCoordinates.c_str());
    outputFile.open (outputNeighborhoodCoordinates, std::ofstream::app);
    //write HEADER
    outputFile << "Neighborhood";
    for(const std::string& feature : get_feaure_names()) 
    {
        outputFile << "\t" << feature;
    }
    outputFile << "\n";
    //write lines
    for(const nodePtr& node : neighborhoodGraph->get_all_nodes())
    {
        outputFile << node->get_name();
        for(const double& dim : node->all_values())
        {
            outputFile << "\t" << dim;
        }
        outputFile << "\n";
    }
    outputFile.close();

    /////////////////////////////
    // 4.) File with cellIDs in each neighborhood: NeighborhoodID in header - all cellIDs in rows for each neighborhood columns
    /////////////////////////////
    std::string outputNeighborhoodCells = outFile + "/" + prefix + "_cells.tsv";
    std::remove(outputNeighborhoodCells.c_str());
    outputFile.open (outputNeighborhoodCells, std::ofstream::app);
    //write HEADER (all neigborhood names)
    const std::vector<nodePtr> neighborHoodPtrVector = neighborhoodGraph->get_all_nodes();
    for(size_t neiborhoodID = 0; neiborhoodID < neighborHoodPtrVector.size(); ++neiborhoodID)
    {
        outputFile << neighborHoodPtrVector.at(neiborhoodID)->get_name();
        if(neiborhoodID < (neighborHoodPtrVector.size()-1)){outputFile << "\t";}
    }  
    outputFile << "\n";
    //write lines
    //iterate over cellIDs (0 - <number cells in neighborhood>), for every id write out the cell with this id for every neighborhood
    size_t maxNeighbourhoodCells = neighborhoodSize;
    for(const auto& nb : neighborhoods) { maxNeighbourhoodCells = std::max(maxNeighbourhoodCells, nb.second.size()); }
    for(size_t cellID = 0; cellID < maxNeighbourhoodCells; ++cellID)
    {
        //for every neighborhood (in columns)
        for(size_t neiborhoodID = 0; neiborhoodID < neighborHoodPtrVector.size(); ++neiborhoodID)
        {
            const std::vector<int>& nbCells = neighborhoods.at(neighborHoodPtrVector.at(neiborhoodID));
            if(cellID < nbCells.size()) { outputFile << nbCells.at(cellID); } else { outputFile << "NA"; }
            if(neiborhoodID < (neighborHoodPtrVector.size()-1)){outputFile << "\t";}
        }
        outputFile << "\n";
    }
    outputFile.close();
}

void Neighborhood::calculate_correlations_variance(unsigned int numberNodes, 
                                          const std::pair<int, int>& correlationpair,
                                          std::unordered_map<const std::pair<int, int>, double, pair_hash>& corrVariance,
                                          int totalCount, double& currentCount)
{
    unsigned int nonNanNodes = 0;
    //calcualte mean first
    double mean = 0;
    for(const nodePtr& node: neighborhoodGraph->get_all_nodes())
    {
        double correlationTmp = corrResult.at(node).correlationResult.at(correlationpair);
        if(!std::isnan(correlationTmp))
        {
            mean += correlationTmp;
            nonNanNodes++;
        }
    }

    double variance = 0;

    //if non-nan nodes are in less than 10% of the neighborhoods (>90% are nan-nodes), 
    // then assign the whole variance adn thereby laplacian to nan
    if(nonNanNodes < numberNodes/10)
    {
       variance = 0; //leave variance explicitely at zero. Laplacian score will therefore also be zero due to division by zero
    }
    else
    {
        mean /= nonNanNodes;
        //calcualte deviations from mean
        for(const nodePtr& node: neighborhoodGraph->get_all_nodes())
        {
            double correlationTmp = corrResult.at(node).correlationResult.at(correlationpair);
            if(!std::isnan(correlationTmp))
            {
                variance += pow((correlationTmp-mean), 2);
            }
        }

        variance /= nonNanNodes;
    }

    threadLock.lock();
    corrVariance.insert(std::pair<const std::pair<int, int>, double>(correlationpair, variance));
    ++currentCount; 
    double perc = currentCount / totalCount;
    printProgress(perc);
    threadLock.unlock();
}

void Neighborhood::calculate_slopes(unsigned int numberNodes,
                                    const std::pair<int, int>& correlationpair,
                                    std::unordered_map<const std::pair<int, int>, double, pair_hash>& slopeVariance,
                                    int totalCount, double& currentCount)
{
        //calcualte mean first
    double mean = 0;
    unsigned int nonNanNodes = 0;
    for(const nodePtr& node: neighborhoodGraph->get_all_nodes())
    {
        double slopeTmp = corrResult.at(node).slopeResult.at(correlationpair);
        if(!std::isnan(slopeTmp))
        {
            mean += slopeTmp;
            nonNanNodes++;
        }
    }

    double variance = 0;
    if(nonNanNodes < numberNodes/10)
    {
       variance = 0; //enforce nan values further down by 0-division
    }
    else
    {
        mean /= nonNanNodes;
        //calcualte deviations from mean
        for(const nodePtr& node: neighborhoodGraph->get_all_nodes())
        {
            double slopeTmp = corrResult.at(node).slopeResult.at(correlationpair);
            if(!std::isnan(slopeTmp))
            {
                variance += pow((slopeTmp-mean), 2);
            }
        }
        variance /= nonNanNodes;
    }

    threadLock.lock();
    slopeVariance.insert(std::pair<const std::pair<int, int>, double>(correlationpair, variance));
    ++currentCount; 
    double perc = currentCount / totalCount;
    printProgress(perc);
    threadLock.unlock();
}

void Neighborhood::laplacian_score(const std::pair<int, int>& pair,
                                   const std::unordered_map<const std::pair<int, int>, double, pair_hash>& corrVariance,
                                   const std::unordered_map<const std::pair<int, int>, double, pair_hash>& slopeVariance,
                                   int totalCount, double& currentCount)
{

    //iterate over adjacency matrix (its a double** weightMatrix in GraphHandler)
    unsigned int nodeNumber = neighborhoodGraph->get_all_nodes().size();

    double corrWeightSum = 0;
    double slopeWeightSum = 0;

    for(unsigned int i = 0; i < (nodeNumber-1); ++i)
    {
        for(unsigned int j = i+1; j < nodeNumber; ++j)
        {
            nodePtr nodeA = neighborhoodGraph->get_node_at(i);
            nodePtr nodeB = neighborhoodGraph->get_node_at(j);

            //edge weight for closeness of values
            double weight = neighborhoodGraph->get_edge_weight_between_nodes(nodeA, nodeB);
            if(weight == 0){continue;}

            //actual feature values for nodes
            double featureNodeACorr = corrResult.at(nodeA).correlationResult.at(pair);
            double featureNodeBCorr = corrResult.at(nodeB).correlationResult.at(pair);
            double featureNodeASlope = corrResult.at(nodeA).slopeResult.at(pair);
            double featureNodeBSlope = corrResult.at(nodeB).slopeResult.at(pair);

            if( !std::isnan(featureNodeACorr) && !std::isnan(featureNodeBCorr))
            {
                corrWeightSum += weight * pow((featureNodeACorr - featureNodeBCorr), 2); //square the result
            }
            if( !std::isnan(featureNodeASlope) && !std::isnan(featureNodeBSlope))
            {
                slopeWeightSum += weight * pow((featureNodeASlope - featureNodeBSlope), 2); //square the result
            }
        }
    }

    //substract variance of the correlations
    corrWeightSum /= corrVariance.at(pair);
    slopeWeightSum /= slopeVariance.at(pair);

    //store the result, for every pair of nodes (here used pair in indices) - indices are the feature indices as in origional data/ and all node feature of Neighborhood Graph
    threadLock.lock();
    correlationLaplacian.insert(std::pair<const std::pair<int, int>, const double>(pair, corrWeightSum));
    slopeLaplacian.insert(std::pair<const std::pair<int, int>, const double>(pair, slopeWeightSum));
    ++currentCount; 
    double perc = currentCount / totalCount;
    printProgress(perc);
    threadLock.unlock();
}

//similar to laplacian score but with randomly shuffled correlations
//we shuffle the whole group of corelations for all nieghborhoods, that means calculated correlations per neighborhood stay together
//but get assigned to a new neighborhood
void Neighborhood::laplacian_significance(
    const std::pair<int,int>& pair,
    const std::unordered_map<const std::pair<int,int>,double,pair_hash>& corrVariance,
    const std::unordered_map<const std::pair<int,int>,double,pair_hash>& slopeVariance,
    const std::vector<CorrelationPropagationResult>& vectorizedResults,
    int totalCount,
    double& currentCount)
{
    thread_local std::mt19937 rng(std::random_device{}());
    const unsigned int nodeNumber = neighborhoodGraph->get_all_nodes().size();

    // ----- cache variances (1 map lookup only) -----
    const double corrVar  = corrVariance.at(pair);
    const double slopeVar = slopeVariance.at(pair);

    // ----- cache observed laplacians -----
    const double corrObserved  = correlationLaplacian.at(pair);
    const double slopeObserved = slopeLaplacian.at(pair);

    // ----- cache feature values (avoid map lookup in loops) -----
    std::vector<double> corrValues(nodeNumber);
    std::vector<double> slopeValues(nodeNumber);

    for(unsigned int i=0;i<nodeNumber;i++)
    {
        corrValues[i]  = vectorizedResults[i].correlationResult.at(pair);
        slopeValues[i] = vectorizedResults[i].slopeResult.at(pair);
    }

    // ----- build edge list once -----
    struct Edge { int i,j; double w; };
    std::vector<Edge> edges;
    edges.reserve(nodeNumber * 5); // guess due to KNN graph, will be more in reality since its not a directed graph

    for(unsigned int i=0;i<nodeNumber-1;i++)
    {
        nodePtr nodeA = neighborhoodGraph->get_node_at(i);

        for(unsigned int j=i+1;j<nodeNumber;j++)
        {
            nodePtr nodeB = neighborhoodGraph->get_node_at(j);

            double w = neighborhoodGraph->get_edge_weight_between_nodes(nodeA,nodeB);

            if(w!=0.0)
                edges.push_back({(int)i,(int)j,w});
        }
    }

    // ----- permutation vector -----
    std::vector<int> perm(nodeNumber);
    std::iota(perm.begin(),perm.end(),0);

    // ----- result storage -----
    std::vector<double> shuffledCorrLaplacianVector(permutations);
    std::vector<double> shuffledSlopeLaplacianVector(permutations);

    // ----- permutation loop -----
    for(int p=0;p<permutations;p++)
    {
        std::shuffle(perm.begin(),perm.end(),rng);

        double corrSum  = 0.0;
        double slopeSum = 0.0;

        for(const Edge& e : edges)
        {
            int pi = perm[e.i];
            int pj = perm[e.j];

            double corrA = corrValues[pi];
            double corrB = corrValues[pj];

            double slopeA = slopeValues[pi];
            double slopeB = slopeValues[pj];

            if(!std::isnan(corrA) && !std::isnan(corrB))
            {
                double d = corrA - corrB;
                corrSum += e.w * d * d;
            }

            if(!std::isnan(slopeA) && !std::isnan(slopeB))
            {
                double d = slopeA - slopeB;
                slopeSum += e.w * d * d;
            }
        }

        shuffledCorrLaplacianVector[p]  = corrSum  / corrVar;
        shuffledSlopeLaplacianVector[p] = slopeSum / slopeVar;
    }

    // ----- p-value calculation -----
    double p_corr_tmp  = 0.0;
    double p_slope_tmp = 0.0;

    for(int i=0;i<permutations;i++)
    {
        if(shuffledCorrLaplacianVector[i] <= corrObserved)
            p_corr_tmp++;

        if(shuffledSlopeLaplacianVector[i] <= slopeObserved)
            p_slope_tmp++;
    }

    p_corr_tmp  /= permutations;
    p_slope_tmp /= permutations;

    // ----- store results -----
    threadLock.lock();

    p_corr.insert({pair,p_corr_tmp});
    p_slope.insert({pair,p_slope_tmp});

    shuffledCorrLaplacians.insert({pair,shuffledCorrLaplacianVector});
    shuffledSlopeLaplacians.insert({pair,shuffledSlopeLaplacianVector});

    ++currentCount;
    double perc = currentCount / totalCount;
    printProgress(perc);

    threadLock.unlock();
}

//for now simply calculate sum of all edge weights
// -> smaller values mean tiny changes between nodes
//in origional paper is a matrix-multiplication formular
void Neighborhood::calculate_laplacian_score(int threads)
{   

    //calculate once the variance that we observe in all correlations/ slopes
    unsigned int numberNodes = neighborhoodGraph->get_all_nodes().size();

    //CORRELATION VARIANCE
    LOCO_OUT << "\tSTEP[6a]:\tCalculate variance of smooth Correlations\n";
    std::unordered_map<const std::pair<int, int>, double, pair_hash> corrVariance;
    //boost::asio::thread_pool pool_corr(threads);
    ThreadPool pool_corr(threads);

    double count = 0;
    for(const std::pair<int, int>& correlationpair : pairs)
    {
        //boost::asio::post(pool_corr, std::bind(&Neighborhood::calculate_correlations_variance, this, 
        //                  numberNodes, print, std::cref(correlationpair), std::ref(corrVariance), pairs.size(), std::ref(count)));
        pool_corr.enqueue([
            this,                       
            numberNodes,                
            &correlationpair,           
            &corrVariance,              
            pairCount = pairs.size(),   
            &count                      
        ]() {
            calculate_correlations_variance(
                numberNodes,
                correlationpair,
                corrVariance,
                pairCount,
                count
            );
        });
    }
    pool_corr.wait_for_tasks();
    printProgress(1);
    LOCO_OUT << "\n";

    //SLOPE VARIANCE
    LOCO_OUT << "\tSTEP[6b]:\tCalculate variance of smooth Slopes\n";
    std::unordered_map<const std::pair<int, int>, double, pair_hash> slopeVariance;
    //boost::asio::thread_pool pool_slope(threads);
    ThreadPool pool_slope(threads);
    count = 0;
    for(const std::pair<int, int>& correlationpair : pairs)
    {
        //boost::asio::post(pool_slope, std::bind(&Neighborhood::calculate_slopes, this, 
        //                  numberNodes, std::cref(correlationpair), std::ref(slopeVariance), pairs.size(), std::ref(count)));
        pool_slope.enqueue([
            this,                       
            numberNodes,                
            &correlationpair,           
            &slopeVariance,              
            pairCount = pairs.size(),   
            &count                      
        ]() {
            calculate_slopes(
                numberNodes,
                correlationpair,
                slopeVariance,
                pairCount,
                count
            );
        });
    }
    pool_slope.wait_for_tasks();
    printProgress(1);
    LOCO_OUT << "\n";

    //FULL LAPLACIAN SCORE
    LOCO_OUT << "\tSTEP[6c]:\tCalculate Laplacian Scores\n";
    //boost::asio::thread_pool pool_lapl(threads);
    ThreadPool pool_lapl(threads);

    count = 0;
    for(const std::pair<int, int>& pair : pairs)
    {
        //boost::asio::post(pool_lapl, std::bind(&Neighborhood::laplacian_score, this, 
        //                  std::cref(pair), std::cref(corrVariance), std::cref(slopeVariance),
        //                  pairs.size(), std::ref(count)));
        pool_lapl.enqueue([
            this,                       
            &pair,                
            &corrVariance,           
            &slopeVariance,              
            pairCount = pairs.size(),   
            &count                      
        ]() {
            laplacian_score(
                pair,                
                corrVariance,           
                slopeVariance,              
                pairCount,   
                count  
            );
        });
    }
    pool_lapl.wait_for_tasks();
    printProgress(1);
    LOCO_OUT << "\n";

    //calcualte significance: calcualte laplacian score for randomly shuffled features across N
    LOCO_OUT << "\tSTEP[6d]:\tCalculate Significance: Laplacian Scores after shuffling correlations\n";
    LOCO_OUT << "\t\t We keep the KNN graph and randomly reassign all correlations/slops of a neighborhood to a different neighborhood\n";
    count = 0;
    //firstly write the correlation/slope results to a vector that is easily shufflable
    std::vector<CorrelationPropagationResult> vectorizedResults;
    vectorizedResults.reserve(corrResult.size());
    for (const auto& [node, result] : corrResult)
    {
        vectorizedResults.push_back(result);
    }
    ThreadPool pool_shuffle(threads);
    for(const std::pair<int, int>& pair : pairs)
    {
        //boost::asio::post(pool_shuffle, std::bind(&Neighborhood::laplacian_significance, this, 
        //                  std::cref(pair), std::cref(corrVariance), std::cref(slopeVariance),
        //                 std::cref(vectorizedResults),
        //                  pairs.size(), std::ref(count)));
        pool_shuffle.enqueue([
            this,                       
            &pair,                
            &corrVariance,           
            &slopeVariance,   
            &vectorizedResults,           
            pairCount = pairs.size(),   
            &count                      
        ]() {
            laplacian_significance(
                pair,                
                corrVariance,           
                slopeVariance, 
                vectorizedResults,             
                pairCount,   
                count  
            );
        });
    }
    pool_shuffle.wait_for_tasks();
    printProgress(1);
    LOCO_OUT << "\n";

}

//fills the reuslt per neighborhood:
//we take correlations from already calcualted neighborhood data
//slopes have to be calculated again

//only do it for interesting correlations within cliquesVector
void Neighborhood::extract_pairs_from_correlation_sets(std::unordered_map<nodePtr, std::shared_ptr<GraphData>>& neighborhoodCorrelations)
{

    //iterte through all cliques & insert all pairs in pair vector
    for(const std::vector<int>& cliques : cliquesVector)
    {
        //for all pairs of nodes from the clique
        for(size_t i = 0; i < (cliques.size()-1); ++i)
        {
            for(size_t j = (i+1); j < cliques.size(); ++j)
            {
                std::pair<int, int> tmpPair(cliques.at(i), cliques.at(j));
                //if the pair has a reported clique but needs to be updated
                if(pairToClique.find(tmpPair) != pairToClique.end())
                {
                    std::vector<std::string> features;
                    for(int cIdx : cliques)
                    {
                        std::string featureName = inputDataOrigional.geneNames.at(cIdx);
                        if(!corrStateGenes.empty())
                        {
                            int correctedIndex = corrStateGenes.at(cIdx);
                            featureName =  inputDataOrigional.geneNames.at(correctedIndex);
                        }
                        features.push_back(featureName);
                    }
                    pairToClique.at(tmpPair).push_back(features);
                }
                else 
                {
                    //if the pair has not been reported so far
                    std::vector<std::string> features;
                    for(int cIdx : cliques)
                    {
                        std::string featureName = inputDataOrigional.geneNames.at(cIdx);
                        if(!corrStateGenes.empty())
                        {
                            int correctedIndex = corrStateGenes.at(cIdx);
                            featureName =  inputDataOrigional.geneNames.at(correctedIndex);
                        }
                        features.push_back(featureName);
                    }
                    std::vector<std::vector<std::string>> stringCliquesVector;
                    stringCliquesVector.push_back(features);
                    pairToClique.insert(std::pair<const std::pair<int, int>, std::vector<std::vector<std::string>>>(tmpPair, stringCliquesVector));
                }
                if(std::find(pairs.begin(), pairs.end(), tmpPair) == pairs.end()) //if the pair was not stored yet
                {
                    pairs.push_back(tmpPair);
                }
            }
        }
    }

    //pre-compute ranks for spearman correlation
    if(correlationType == "spearman")
    {
        for(auto const& [center, tmpData] : neighborhoodCorrelations) 
        {
            for(int i = 0; i < tmpData->number_of_nodes(); ++i) 
            {
                tmpData->get_node_at(i)->compute_ranks(); 
            }
        }
    }

    //PROCESS ALL THOSE PAIRS
    for(const nodePtr& neighborhoodCenter : centralNeighborhoodPtrs)
    {
        //temporary result for this neighborhood
        CorrelationPropagationResult tmpResult;

        //all feature-pairs (correlations between two features)
        for(const std::pair<int, int>& pair : pairs)
        {
            //get current neighborhood Correlation Data
            std::shared_ptr<GraphData> tmpData = neighborhoodCorrelations.at(neighborhoodCenter);

            //store all CORRELATIONS for those pairs (get them from already calculated neighborhoodCorrelations)
            nodePtr featureNodeA = tmpData->get_node_at(pair.first);
            nodePtr featureNodeB = tmpData->get_node_at(pair.second);
            //re-calcualte correlations (before we had absolute values)
            //const double corr = tmpData->get_distance_between_nodes(featureNodeA, featureNodeB);

            double corr;
            if(correlationType == "spearman")
            {
                corr = calculate_correlation_coefficient(featureNodeA->all_ranked_values(), featureNodeB->all_ranked_values());
            }
            else if(correlationType == "pearson")
            {
                corr = calculate_correlation_coefficient(featureNodeA->all_values(), featureNodeB->all_values());
            }

            //if (std::isnan(corr)) 
            //{
            //    std::cout << "WARNING: Correlation value is NaN" << " between features: " << featureNodeA->get_name() << " and " << featureNodeB->get_name() << "\n";    
            //}
                
            tmpResult.correlationResult.insert(std::pair< std::pair<int, int>, double >(pair, corr));

            //store all SLOPES for those pairs (have to calculate new)
            const double slope = calculate_slope(featureNodeA->all_values(), featureNodeB->all_values());
            tmpResult.slopeResult.insert(std::pair< std::pair<int, int>, double >(pair, slope));
        }

        corrResult.insert(std::pair<nodePtr, CorrelationPropagationResult>(neighborhoodCenter, tmpResult));
    }
    
}

void Neighborhood::detect_cliques_in_neighborhood(nodePtr neighborhoodCenter, const double& correlationStrengthCutoff, int minCliqueSize,
                                                  std::unordered_map<nodePtr, std::vector<std::vector<int>>>& cliquesPerNeighborhood, 
                                                  std::unordered_map<nodePtr, std::shared_ptr<GraphData>>& neighborhoodCorrelations,
                                                  int totalCount, double& currentCount)
{
        //calculate all correlations/ slopes of them
        SingleCellData inputDataTmp = filter_singleCelldata(inputDataOrigional, neighborhoods.at(neighborhoodCenter));

        //for the protein graph calculate brute force all distances instead of reading form a KD-tree
        // this way we can filter correlations based on value
        unsigned int proteinGraphKnn = 0;
        std::shared_ptr<GraphData> correlationData;
        if(correlationType == "spearman")
        {
            correlationData = std::make_shared<GraphData>(inputDataTmp, corrStateGenes, proteinGraphKnn, &GraphIni::protein_correlation_graph_spearman);
        }
        else if(correlationType == "pearson")
        {
            correlationData = std::make_shared<GraphData>(inputDataTmp, corrStateGenes, proteinGraphKnn, &GraphIni::protein_correlation_graph_pearson);
        }
        
        GraphHandler corrGraphBuilder = GraphHandler(correlationData, 0, correlationStrengthCutoff, -1);
        corrGraphBuilder.create_graph();

        //results are stored here afterwards
        std::vector<std::vector<int>> cliqueVectorRaw;
        if(correlatedSetMode == 1)
        {
            corrGraphBuilder.calculate_connected_components(cliqueVectorRaw, minCliqueSize);
        }
        else if(correlatedSetMode == 0)
        {
            corrGraphBuilder.calculate_fully_connected_components(cliqueVectorRaw, minCliqueSize);
        }
        else if(correlatedSetMode >= 2)
        {
            corrGraphBuilder.calculate_min_edge_connected_components(cliqueVectorRaw, minCliqueSize, correlatedSetMode);
        }
        else
        {
            LOCO_OUT << "Invalid mode for detection of correlated sets: fallback to connected component -> mode 0\n";
            corrGraphBuilder.calculate_connected_components(cliqueVectorRaw, minCliqueSize);
        }

        //insert cliques in map<neighborhoodID -> cliques>
        threadLock.lock();
        cliquesPerNeighborhood.insert(std::make_pair(neighborhoodCenter, cliqueVectorRaw));
        neighborhoodCorrelations.insert(std::pair<nodePtr, std::shared_ptr<GraphData>>(neighborhoodCenter, correlationData));
        ++currentCount; 
        double perc = currentCount / totalCount;
        printProgress(perc);
        threadLock.unlock();
}

//filters cliques that do not 'randomly' occur in one neighborhood, but the ones that are shared also among the neighbors of this neighborhood
//(NOT IUSED AT THIS POINT: Since we do not look at subsets of cliques it almost never happens that the EXACT CLIQUE of several proteins is the same...)
void Neighborhood::filter_cliques_present(nodePtr neighborhoodCenter, std::unordered_map<nodePtr, std::vector<std::vector<int>>>& cliquesPerNeighborhood)
{
        //get closest neighborhoods (now only knn neighbors in graph)
        //TODO: connect via graph of whole data and BFS through it for neighbors (igraph_bfs_callback.c with own callback function)
        //for now just get 'numberNeighbors' closest neighbors for our main node (maybe later use graph and BFS)
        int numberNeighbors = 5;
        std::vector<nodePtr> neighborIDs = neighborhoodGraph->get_data()->get_adjacent_nodes_knn(neighborhoodCenter, numberNeighbors);

        //find common cliques in all neighborhoods
        //for every clique in actual nHood
        for(std::vector<int> clique : cliquesPerNeighborhood.at(neighborhoodCenter))
        {

            bool cliqueInNeighbors = true;
            //check if it is contained in all neighboring neighborhoods
            for(nodePtr neighboringNHood : neighborIDs)
            {
                std::vector<std::vector<int>> nCliques = cliquesPerNeighborhood.at(neighboringNHood);

                //the cliques are all sorted in ascending order
                if(std::find(nCliques.begin(), nCliques.end(), clique) == nCliques.end()) //if the clique is not in the neighboring cliques
                {
                    cliqueInNeighbors = false;
                    break;
                }
            }

            //if clique is not existing in neighbors reject it (filtering of cliques that are only in ONE neighborhood and probably outliers)
            if(!cliqueInNeighbors){continue;}

            threadLock.lock();
            //if we managed to run through for loop, this clique is in all neighboring-neighborHoods
            if(std::find(cliquesVector.begin(), cliquesVector.end(), clique) == cliquesVector.end()) //if the clique does not exist yet
            {
                cliquesVector.push_back(clique);
            }
            threadLock.unlock();
        }
}

// Helper to create a fast bitmask for pruning
uint64_t create_signature(const std::vector<int>& vec) {
    uint64_t mask = 0;
    for (int x : vec) mask |= (1ULL << (static_cast<uint32_t>(x) % 64));
    return mask;
}

void Neighborhood::filter_consistent_correlation_sets_sota(
        const std::unordered_map<nodePtr, std::vector<std::vector<int>>>& cliquesPerNeighborhood,
        int minCliqueSize)
{
    // 1. GENERATE CANDIDATES: all possible subsets of found sets of allowed size
    std::set<std::vector<int>> candidateSets;
    for(const auto& clique : cliquesVector) {
        std::vector<int> sortedClique = clique;
        std::sort(sortedClique.begin(), sortedClique.end());
        auto subsets = neighborhoodCalculations::generate_subsets(sortedClique, minCliqueSize);
        for(auto& s : subsets){ candidateSets.insert(std::move(s));}
    }

    // 2. PRE-PROCESS NEIGHBORHOODS (The Speed Boost)
    // We sort everything once and pre-calculate bitmasks
    struct IndexedNeighborhood {
        std::vector<std::vector<int>> sortedCliques;
        std::vector<uint64_t> masks;
    };
    
    std::vector<IndexedNeighborhood> processedData;
    processedData.reserve(cliquesPerNeighborhood.size());

    for (auto const& [node, cliques] : cliquesPerNeighborhood) {
        IndexedNeighborhood idx;
        for (auto c : cliques) {
            std::sort(c.begin(), c.end());
            idx.masks.push_back(create_signature(c));
            idx.sortedCliques.push_back(std::move(c));
        }
        processedData.push_back(std::move(idx));
    }

    const int neighborhoodCount = processedData.size();
    const int required = std::ceil(minimumCorrSetAbundance * neighborhoodCount);
    std::vector<std::vector<int>> result;

    // check all candidates
    for(const auto& candidate : candidateSets) {
        uint64_t candMask = create_signature(candidate);
        int count = 0;

        for(const auto& neighborhood : processedData) {
            bool found = false;
            for(size_t i = 0; i < neighborhood.sortedCliques.size(); ++i) {
                // PRUNING STEP: Bitmask check is nearly free
                if ((candMask & neighborhood.masks[i]) != candMask) continue;
                
                // Size check: Candidate can't be a subset if it's bigger
                if (candidate.size() > neighborhood.sortedCliques[i].size()) continue;

                // FINAL CHECK: std::includes is the fastest way to check subset on sorted vectors
                if (std::includes(neighborhood.sortedCliques[i].begin(), neighborhood.sortedCliques[i].end(),
                                  candidate.begin(), candidate.end())) {
                    found = true;
                    break;
                }
            }
            if(found) count++;
            if(count >= required) break; // Early exit for this candidate
        }

        if(count >= required) result.push_back(candidate);
    }

    cliquesVector = std::move(result);
}

//calcualte all apirwise corrs within a N and push it into shared tempCalculationMatrix
void Neighborhood::calculate_correlations_for_N(nodePtr neighborhoodCenter, size_t neighborhood_idx, 
                                                const double& corrThreshold, FlatMatrix& tempCalculationMatrix,
                                                const std::vector<std::pair<int, int>>& tmpAllPairs , 
                                                std::vector<std::atomic<int>>& tmpCorrelationCountAboveThreshold, 
                                                std::atomic<int>& currentCount)
{
    //1.) subset the data used for this neighbourhood

    // add filtering step to choose only corrGEnes
    std::vector<std::vector<double>> pointCloud = return_filtered_point_cloud(inputDataOrigional, neighborhoods.at(neighborhoodCenter), corrStateGenes);
    
    // make rank computation independent of nodes: do it on the transpose of raw data to get all features quickly
    //2.) compute ranks of values for spearman correlation
    std::vector<std::vector<double>> featureMajorData = transpose_to_feature_major(pointCloud);

    //1b.) features with too little variation in this neighbourhood: fewer than minVaryingFraction of the cells differ from
    // the most frequent value (e.g. all zeroes and a single 1) -> all their correlations are NaN in this neighbourhood
    // (done BEFORE ranking, the ties are the same but this avoids depending on rank precision)
    std::vector<bool> featureVaries(featureMajorData.size(), true);
    if(minVaryingFraction > 0.0)
    {
        for (size_t f = 0; f < featureMajorData.size(); ++f)
        {
            const double n = static_cast<double>(featureMajorData[f].size());
            // small tolerance: e.g. 0.07 * 100 = 7.0000000000000009 in floating point, which would require 8 cells
            featureVaries[f] = static_cast<double>(count_non_modal_values(featureMajorData[f])) >= minVaryingFraction * n - 1e-9;
        }
    }

    if(correlationType == "spearman")
    {
        for (size_t f = 0; f < featureMajorData.size(); ++f) 
        {
            rankify(featureMajorData[f]); 
        }
    }

    // 3.) Calculate correlation coefficients and update FlatMatrix safely!
    for(size_t p = 0; p < tmpAllPairs.size(); ++p)
    {
        int featA = tmpAllPairs[p].first;
        int featB = tmpAllPairs[p].second;

        // Calculate correlation (assuming your function takes two vectors)
        double corr = std::numeric_limits<double>::quiet_NaN();
        if(featureVaries[featA] && featureVaries[featB])
        {
            corr = calculate_correlation_coefficient(featureMajorData[featA], featureMajorData[featB]);
        }

        // Write directly to the matrix lock-free! (Using the neighborhood_idx as the row)
        tempCalculationMatrix(neighborhood_idx, p) = corr;

        // If above threshold, safely increment the global counter without a lock
        if (std::abs(corr) >= corrThreshold) 
        {
            // std::atomic safely handles multiple threads incrementing this at once
            tmpCorrelationCountAboveThreshold[p].fetch_add(1, std::memory_order_relaxed);
        }
    }

    //update all pairwise corrs for this N: add it to FlatMatrix
    currentCount++;
}

void Neighborhood::step_2_calculate_correlation(const double& corrThreshold, const int threads)
{

    //1.) INITIALIZE DATA STRUCTURES

    // a) Determine your feature size
    size_t geneSize = inputDataOrigional.geneNames.size();
    if (!corrStateGenes.empty()) 
    {
        geneSize = corrStateGenes.size();
    }
    size_t totalPairs = (geneSize * (geneSize - 1)) / 2;

    //b) temporary TRANSPOSE matrix to quickly iterae over N and calc corrs for every pair
    //(later transpose: neighbourhoodCorrs has outer vector for corr-pairs to quickly calc laplacian for corr-pairs across N)
    std::vector<std::pair<int, int>> tmpAllPairs; //later keep only the pairs with x% above corr-threshold
    tmpAllPairs.reserve(totalPairs);
    for (int i = 0; i < static_cast<int>(geneSize); ++i) 
    {
    for (int j = i + 1; j < static_cast<int>(geneSize); ++j) 
        {
            // emplace_back constructs the pair directly in place, 
            // which is slightly faster than push_back({i, j})
            tmpAllPairs.emplace_back(i, j); 
        }
    }
    FlatMatrix tempCalculationMatrix(neighbourhoodNum, totalPairs);

    //c) create structures
    // temporary structures to count correlation-pair with >threshold correlations
    //count in which feature-pair correlation is above threshold
    std::vector<std::atomic<int>> tmpCorrelationCountAboveThreshold(totalPairs); 

    //3.) calcualte all correlations in threads
    std::atomic<int> completedNCount{0};
    size_t neighborhood_idx = 0;
    {
        ThreadPool pool_correlations(threads);
        for( nodePtr neighborhoodCenter : centralNeighborhoodPtrs)
        {
            pool_correlations.enqueue([
                this,                              
                neighborhoodCenter,
                neighborhood_idx,                        
                &corrThreshold, 
                &tempCalculationMatrix,  
                &tmpAllPairs,                                                
                &tmpCorrelationCountAboveThreshold,    
                &completedNCount                        
            ]() {
                calculate_correlations_for_N(
                    neighborhoodCenter,      
                    neighborhood_idx,                  
                    corrThreshold,                              
                    tempCalculationMatrix,
                    tmpAllPairs,     
                    tmpCorrelationCountAboveThreshold,               
                    completedNCount  
                );
            });

            neighborhood_idx++; // Increment for the next loop iteration
        }

        // 2. The Main R Thread prints the progress while workers work!
        int last_val = -1;

        while (completedNCount < neighbourhoodNum) 
        {
            double percentage = (double)completedNCount.load() / neighbourhoodNum;
            int val = (int)(percentage * 100);

            // Only call the print function if the integer percentage actually changed
            if (val != last_val) 
            {
                printProgressMainThread(percentage);
                last_val = val;
            }

            // Sleep briefly so the main thread doesn't hog the CPU while waiting
            std::this_thread::sleep_for(std::chrono::milliseconds(50));
        }

        // 3. Guarantee it prints 100% and a newline at the very end
        pool_correlations.wait_for_tasks();
    }
    printProgressMainThread(1.0);
    LOCO_OUT << "\n";

    //4.) Create result for correlations (per N)
    // get only corrs that are x% over threshood and safe them in neighbourhoodCorrs
    // do this inside the init function of the neighbourhoodCorrelations and instead of the matrix transpose do a 
    // in-place writing into the new amtrix (transposed) only for coors > threshold
    const int minRequired = std::ceil(minimumCorrSetAbundance * neighbourhoodNum);    
    neighbourhoodCorrs.init(tempCalculationMatrix, 
                      tmpAllPairs, 
                      tmpCorrelationCountAboveThreshold, 
                      minRequired, 
                      centralNeighborhoodPtrs, 
                      inputDataOrigional.geneNames,
                      corrStateGenes);

}

// ---------------------------------------------------------------------------
// LAPLACIAN SCORE (He, Cai & Niyogi 2005) of the correlations of one feature pair across the neighbourhood graph
//
// Let c_i be the correlation of the pair in neighbourhood i and w_ij the edge weight between neighbourhoods i,j.
// Neighbourhoods where the correlation is NOT defined (NaN) are removed from the graph (induced subgraph on the
// valid neighbourhoods V, same edges and weights, nothing is reconnected or imputed). In this subgraph:
//
//   d_i  = sum_{j in V} w_ij                      (degree)
//   mu_w = sum_i d_i c_i / sum_i d_i              (degree-weighted mean)
//   L    = sum_{i<j in V} w_ij (c_i - c_j)^2  /  sum_{i in V} d_i (c_i - mu_w)^2
//
// L is ~0 for a smooth signal, ~1 for a signal without relation to the graph and up to 2 for alternating signals.
// ---------------------------------------------------------------------------
namespace
{
    // the induced subgraph of the neighbourhood graph on all neighbourhoods with a defined correlation
    struct ValidSubgraph
    {
        std::vector<int> validNodes;      // neighbourhoods with a defined correlation
        std::vector<size_t> edgeIdx;      // indices into the global edge list: edges with both endpoints valid
        std::vector<double> degree;       // degree within the valid subgraph (size: number of neighbourhoods)
        double totalDegree = 0.0;         // sum of all degrees (= 2 * sum of valid edge weights)
    };

    // a correlation is defined if it is not NaN (-2.0 is the 'not set' default of the FlatMatrix)
    inline bool is_defined_correlation(const double c)
    {
        return !std::isnan(c) && c != -2.0;
    }

    ValidSubgraph build_valid_subgraph(const double* pairCorrs, const size_t numberNodes, const std::vector<Edge>& edges)
    {
        ValidSubgraph graph;
        graph.degree.assign(numberNodes, 0.0);
        for(size_t n = 0; n < numberNodes; ++n)
        {
            if(is_defined_correlation(pairCorrs[n])) { graph.validNodes.push_back(static_cast<int>(n)); }
        }
        for(size_t e = 0; e < edges.size(); ++e)
        {
            const Edge& edge = edges[e];
            if(is_defined_correlation(pairCorrs[edge.i]) && is_defined_correlation(pairCorrs[edge.j]))
            {
                graph.edgeIdx.push_back(e);
                graph.degree[edge.i] += edge.w;
                graph.degree[edge.j] += edge.w;
                graph.totalDegree += 2.0 * edge.w;
            }
        }
        return graph;
    }

    // Laplacian score for the values c (indexed by neighbourhood, only valid neighbourhoods are read).
    // Returns NaN if the score is undefined (no edges or no variation among the connected neighbourhoods).
    // denominator: degree-weighted variance (out)
    double degree_weighted_laplacian(const double* c, const ValidSubgraph& graph, const std::vector<Edge>& edges, double& denominator)
    {
        denominator = 0.0;
        if(graph.totalDegree <= 0.0) { return std::numeric_limits<double>::quiet_NaN(); }

        double numerator = 0.0;
        for(const size_t e : graph.edgeIdx)
        {
            const double diff = c[edges[e].i] - c[edges[e].j];
            numerator += edges[e].w * diff * diff;
        }

        double weightedSum = 0.0;
        for(const int n : graph.validNodes) { weightedSum += graph.degree[n] * c[n]; }
        const double weightedMean = weightedSum / graph.totalDegree;

        for(const int n : graph.validNodes)
        {
            const double dev = c[n] - weightedMean;
            denominator += graph.degree[n] * dev * dev;
        }

        // no variation among the connected neighbourhoods
        if(denominator < 1e-18) { return std::numeric_limits<double>::quiet_NaN(); }
        return numerator / denominator;
    }
}

void Neighborhood::calculate_laplacian_score_for_pair(const int featurePairIdx, const std::vector<Edge>& edges)
{
    const double* pairCorrs = neighbourhoodCorrs.pairToListOfAllNCorrs.get_row(featurePairIdx);
    const ValidSubgraph graph = build_valid_subgraph(pairCorrs, neighbourhoodNum, edges);

    // no minimum number of defined neighbourhoods here: pairs with too few were already removed by minNonNAFraction
    // (step_2b_filter_pairs_by_non_na). The score is NaN only if it is undefined (no edges or no variation among the
    // connected neighbourhoods); these pairs are removed in step_3_calculate_laplacian_score
    double denominator = 0.0;
    double score = std::numeric_limits<double>::quiet_NaN();
    if(!graph.validNodes.empty())
    {
        score = degree_weighted_laplacian(pairCorrs, graph, edges, denominator);
    }

    laplacianScores.degreeWeightedVariance[featurePairIdx] = denominator;
    laplacianScores.L[featurePairIdx] = score;
}

// significance of the Laplacian score by permutations:
// the mask of neighbourhoods with a defined correlation stays FIXED. Only the defined correlations are shuffled
// among the valid neighbourhoods, so every permutation is scored on exactly the same valid subgraph
// (same edges, same degrees) as the observed score. NaNs never move.
void Neighborhood::laplacian_significance_for_pair(size_t pair_idx, const std::vector<Edge>& edges, std::atomic<int>& currentCount)
{
    const double observedL = laplacianScores.L[pair_idx];

    // Safety Check: If L is NaN (undefined score), p-value is also NaN
    if (std::isnan(observedL))
    {
        laplacianScores.p_values[pair_idx] = std::numeric_limits<double>::quiet_NaN();
        currentCount++;
        return;
    }

    const double* pairCorrs = neighbourhoodCorrs.pairToListOfAllNCorrs.get_row(pair_idx);
    const ValidSubgraph graph = build_valid_subgraph(pairCorrs, neighbourhoodNum, edges);

    // defined correlations (shuffled) and the vector with the permuted values at the valid neighbourhoods
    std::vector<double> validValues;
    validValues.reserve(graph.validNodes.size());
    for(const int n : graph.validNodes) { validValues.push_back(pairCorrs[n]); }
    std::vector<double> permutedCorrs(pairCorrs, pairCorrs + neighbourhoodNum);

    thread_local std::mt19937 rng(std::random_device{}());
    double p_count = 0.0;

    for(int p = 0; p < permutations; ++p)
    {
        std::shuffle(validValues.begin(), validValues.end(), rng);
        for(size_t k = 0; k < graph.validNodes.size(); ++k)
        {
            permutedCorrs[graph.validNodes[k]] = validValues[k];
        }

        double denominatorPermuted = 0.0;
        const double shuffledL = degree_weighted_laplacian(permutedCorrs.data(), graph, edges, denominatorPermuted);

        // an undefined shuffled score means no variation among the connected neighbourhoods = perfectly smooth:
        // count it as at least as smooth as the observed score
        if(std::isnan(shuffledL) || shuffledL <= observedL)
        {
            p_count++;
        }
    }

    // Lock-free write to pre-allocated results
    // add a pseudocount (observed statistic counts as one permutation): p can never be exactly 0
    laplacianScores.p_values[pair_idx] = (permutations > 0) ?
        (p_count + 1.0) / (permutations + 1.0) : std::numeric_limits<double>::quiet_NaN();
    currentCount++;
}

void Neighborhood::step_3_calculate_laplacian_score(const int threads)
{
    laplacianScores.pairNames = neighbourhoodCorrs.pairNames;
    size_t totalPairs = laplacianScores.pairNames.size();

    // 1.) edge list of the neighbourhood graph, built exactly ONCE (used for scores AND permutations)
    std::vector<Edge> globalEdges;
    globalEdges.reserve(neighbourhoodNum * 5); // guess due to KNN graph, will be more in reality since its not a directed graph
    for(unsigned int i = 0; i + 1 < neighbourhoodNum; ++i)
    {
        for(unsigned int j = i + 1; j < neighbourhoodNum; ++j)
        {
            double w = neighborhoodGraph->get_edge_weight_between_nodes(i, j);
            if(w != 0.0)
            {
                globalEdges.push_back({(int)i, (int)j, w});
            }
        }
    }

    // 2.) calcualte laplacian score (degree-weighted, on the subgraph of neighbourhoods with a defined correlation)
    LOCO_OUT << "\t Calculate Laplacian Score for found correlations\n";
    // pre-allocate result vectors for thread-safe writing
    laplacianScores.degreeWeightedVariance.resize(totalPairs, 0.0);
    laplacianScores.L.resize(totalPairs, 0.0);         // Holds final Laplacian scores
    {
        ThreadPool pool_L(threads);
        for (size_t p = 0; p < totalPairs; ++p)
        {
            pool_L.enqueue([this, p, &globalEdges]() {
                this->calculate_laplacian_score_for_pair(p, globalEdges);
            });
        }
        pool_L.wait_for_tasks();
    }

    // pairs with an undefined score (no edges or no variation among the connected neighbourhoods) are removed:
    // they get no permutations/ p-value and are not part of the result
    {
        std::vector<size_t> definedIdx;
        for(size_t p = 0; p < totalPairs; ++p)
        {
            if(!std::isnan(laplacianScores.L[p])) { definedIdx.push_back(p); }
        }
        LOCO_OUT << "\t " << definedIdx.size() << " of " << totalPairs << " pairs have a defined Laplacian score\n";
        if(definedIdx.size() < totalPairs) { keep_only_pairs(definedIdx); }
        totalPairs = laplacianScores.L.size();
    }

    // 3.) calcualte significance values (only print status here)
    LOCO_OUT << "\t Calculate significance for Laplacian Score across " << std::to_string(permutations) << " permutations\n";
    // Pre-allocate p-values
    laplacianScores.p_values.resize(totalPairs, 0.0);
    // Initialize the atomic counter
    std::atomic<int> completedPairsCount{0};
    {
        ThreadPool pool_signifiance(threads);
        // Enqueue all tasks
        for (size_t p = 0; p < totalPairs; ++p)
        {
            // Capture the atomic counter by reference
            pool_signifiance.enqueue([this, p, &globalEdges, &completedPairsCount]() 
            {
                this->laplacian_significance_for_pair(p, globalEdges, completedPairsCount);
            });
        }
        // MAIN THREAD PROGRESS LOOP
        // The main thread waits here and updates the UI instead of blocking entirely
        while (completedPairsCount < totalPairs)
        {
            // Calculate percentage
            double perc = static_cast<double>(completedPairsCount) / totalPairs;
            printProgressMainThread(perc);
            // Put the main thread to sleep for 100ms so it doesn't waste CPU cycles 
            // while checking the counter.
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        // 6. Ensure everything is cleanly wrapped up before destroying the pool
        pool_signifiance.wait_for_tasks(); 
    }
    // Print final 100% state and a newline
    printProgressMainThread(1.0);
    LOCO_OUT << "\n";

}


/*
OLD CODE AS EXMAPLE

      //calculate all correlations/ slopes of them
        SingleCellData inputDataTmp = filter_singleCelldata(inputDataOrigional, neighborhoods.at(neighborhoodCenter));

        //for the protein graph calculate brute force all distances instead of reading form a KD-tree
        // this way we can filter correlations based on value
        unsigned int proteinGraphKnn = 0;
        std::shared_ptr<GraphData> correlationData;
        if(correlationType == "spearman")
        {
            correlationData = std::make_shared<GraphData>(inputDataTmp, corrStateGenes, proteinGraphKnn, &GraphIni::protein_correlation_graph_spearman);
        }
        else if(correlationType == "pearson")
        {
            correlationData = std::make_shared<GraphData>(inputDataTmp, corrStateGenes, proteinGraphKnn, &GraphIni::protein_correlation_graph_pearson);
        }
        
        GraphHandler corrGraphBuilder = GraphHandler(correlationData, 0, correlationStrengthCutoff, -1);
        corrGraphBuilder.create_graph();

        //results are stored here afterwards
        std::vector<std::vector<int>> cliqueVectorRaw;
        if(correlatedSetMode == 1)
        {
            corrGraphBuilder.calculate_connected_components(cliqueVectorRaw, minCliqueSize);
        }
        else if(correlatedSetMode == 0)
        {
            corrGraphBuilder.calculate_fully_connected_components(cliqueVectorRaw, minCliqueSize);
        }
        else if(correlatedSetMode >= 2)
        {
            corrGraphBuilder.calculate_min_edge_connected_components(cliqueVectorRaw, minCliqueSize, correlatedSetMode);
        }
        else

*/

// Benjamini-Hochberg adjusted p-values over ALL pairs with a defined p-value: exactly R's
// p.adjust(p, method = "BH") -- NaN p-values stay NaN and do not count as tests,
// q = min(1, cummin((n / i) * p)) over the p-values sorted in decreasing order
void Neighborhood::step_3b_adjust_p_values()
{
    const std::vector<double>& p = laplacianScores.p_values;
    laplacianScores.p_adj.assign(p.size(), std::numeric_limits<double>::quiet_NaN());

    std::vector<size_t> defined;
    defined.reserve(p.size());
    for(size_t i = 0; i < p.size(); ++i)
    {
        if(!std::isnan(p[i])) { defined.push_back(i); }
    }
    const size_t n = defined.size();
    if(n == 0) { return; }

    // decreasing p (ties keep their index order, like R's order(); tied p-values get the same q anyway)
    std::stable_sort(defined.begin(), defined.end(), [&p](size_t a, size_t b) { return p[a] > p[b]; });

    double runningMin = std::numeric_limits<double>::infinity();
    for(size_t k = 0; k < n; ++k)
    {
        const size_t rank = n - k; // rank of this p-value in increasing order (n ... 1)
        const double q = (static_cast<double>(n) / static_cast<double>(rank)) * p[defined[k]];
        runningMin = std::min(runningMin, q);
        laplacianScores.p_adj[defined[k]] = std::min(1.0, runningMin);
    }
}

// keep only the pairs at keepIdx (indices into the current pair vectors, increasing): all per-pair data
// (pair list/ names, correlation matrix rows, Laplacian results) is subset, the original pair order is kept.
// Laplacian results that are not calculated yet (empty vectors) stay empty, so this can be used at every step
void Neighborhood::keep_only_pairs(const std::vector<size_t>& keepIdx)
{
    const size_t numNeighbourhoods = neighbourhoodCorrs.pairToListOfAllNCorrs.cols;

    std::vector<std::pair<int, int>> pairList;
    std::vector<std::string> pairNames;
    FlatMatrix corrs(keepIdx.size(), numNeighbourhoods);

    for(size_t k = 0; k < keepIdx.size(); ++k)
    {
        const size_t p = keepIdx[k];
        pairList.push_back(neighbourhoodCorrs.pairList[p]);
        pairNames.push_back(neighbourhoodCorrs.pairNames[p]);
        const double* row = neighbourhoodCorrs.pairToListOfAllNCorrs.get_row(p);
        std::copy(row, row + numNeighbourhoods, &corrs.data[k * numNeighbourhoods]);
    }

    auto subset = [&keepIdx](auto& values) {
        if(values.empty()) { return; }
        std::remove_reference_t<decltype(values)> kept;
        kept.reserve(keepIdx.size());
        for(const size_t p : keepIdx) { kept.push_back(values[p]); }
        values = std::move(kept);
    };
    subset(laplacianScores.L);
    subset(laplacianScores.p_values);
    subset(laplacianScores.p_adj);
    subset(laplacianScores.degreeWeightedVariance);
    subset(laplacianScores.sortRank);

    neighbourhoodCorrs.pairList = std::move(pairList);
    neighbourhoodCorrs.pairNames = pairNames;
    neighbourhoodCorrs.pairToListOfAllNCorrs = std::move(corrs);
    laplacianScores.pairNames = std::move(pairNames);
}

// minNonNAFraction > 0: keep only pairs with a defined (non-NaN) correlation in MORE than this fraction of all
// neighbourhoods. Applied right after the correlations, BEFORE the Laplacian score: removed pairs get no score,
// no p-value and are not part of the result (and do not count as tests for p_adj)
void Neighborhood::step_2b_filter_pairs_by_non_na(double minNonNAFraction)
{
    const size_t numPairs = neighbourhoodCorrs.pairNames.size();
    if(minNonNAFraction <= 0.0 || numPairs == 0) { return; }

    const size_t numNeighbourhoods = neighbourhoodCorrs.pairToListOfAllNCorrs.cols;
    std::vector<size_t> keepIdx;
    for(size_t p = 0; p < numPairs; ++p)
    {
        const double* row = neighbourhoodCorrs.pairToListOfAllNCorrs.get_row(p);
        size_t defined = 0;
        for(size_t n = 0; n < numNeighbourhoods; ++n)
        {
            if(!std::isnan(row[n])) { ++defined; }
        }
        // fraction computed like R's mean() of a logical vector (long double sum / n)
        const double fraction = static_cast<double>(static_cast<long double>(defined) / static_cast<long double>(numNeighbourhoods));
        if(fraction > minNonNAFraction) { keepIdx.push_back(p); }
    }
    LOCO_OUT << "minNonNAFraction: " << keepIdx.size() << " of " << numPairs << " pairs have a correlation in > "
             << minNonNAFraction * 100 << "% of neighbourhoods\n";
    // only copy the pairs if any were removed (keepIdx is then all pairs in the original order)
    if(keepIdx.size() < numPairs) { keep_only_pairs(keepIdx); }
}

// pair filter, applied AFTER p-values and p_adj are calculated for all tested pairs (p_adj corrects for all tested
// pairs) and BEFORE the feature sets are calculated (feature sets only use the kept pairs):
//   numberCorrelations > 0: keep the significant pairs (p_adj <= significanceCutoff) and of those the
//   <numberCorrelations> with the lowest Laplacian score (ties: lower p_adj, then output order)
// (minNonNAFraction is applied before the Laplacian score, see step_2b_filter_pairs_by_non_na)
void Neighborhood::step_3b_filter_pairs(unsigned int numberCorrelations, double significanceCutoff)
{
    const size_t numPairs = laplacianScores.L.size();

    // output order: rank of every pair when ALL pairs are sorted by L (same sort as for the output before)
    {
        std::vector<size_t> indices(numPairs);
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(),
            [this](size_t a, size_t b) { return laplacianScores.L[a] < laplacianScores.L[b]; });
        laplacianScores.sortRank.assign(numPairs, 0);
        for(size_t r = 0; r < numPairs; ++r) { laplacianScores.sortRank[indices[r]] = r; }
    }

    const size_t remainingPairs = laplacianScores.L.size();
    if(numberCorrelations > 0 && remainingPairs > 0)
    {
        std::vector<size_t> sigIdx;
        for(size_t p = 0; p < remainingPairs; ++p)
        {
            const double q = laplacianScores.p_adj[p];
            if(!std::isnan(q) && q <= significanceCutoff) { sigIdx.push_back(p); }
        }
        std::sort(sigIdx.begin(), sigIdx.end(), [this](size_t a, size_t b) {
            if(laplacianScores.L[a] != laplacianScores.L[b]) { return laplacianScores.L[a] < laplacianScores.L[b]; }
            if(laplacianScores.p_adj[a] != laplacianScores.p_adj[b]) { return laplacianScores.p_adj[a] < laplacianScores.p_adj[b]; }
            return laplacianScores.sortRank[a] < laplacianScores.sortRank[b];
        });
        const size_t numberSig = sigIdx.size();
        if(sigIdx.size() > numberCorrelations) { sigIdx.resize(numberCorrelations); }
        std::sort(sigIdx.begin(), sigIdx.end()); // back to the original pair order
        LOCO_OUT << "numberCorrelations: " << numberSig << " of " << remainingPairs << " pairs are significant (p_adj <= "
                 << significanceCutoff << "), keeping " << sigIdx.size() << " with the lowest Laplacian score\n";
        if(sigIdx.size() < remainingPairs) { keep_only_pairs(sigIdx); }
    }
}

void Neighborhood::step_4_calculate_feature_sets(int minFeatureSetSize)
{
    // create graph with features as nodes
    // use the second GraphHandler initializer: it creates directly a graph from nodes/ edges
    // input are the feature names and the pairList
    GraphHandler consensusCorrelationGraph = GraphHandler(neighbourhoodCorrs.featureNames, neighbourhoodCorrs.pairList);
    std::vector<std::vector<int>> featureSetsRaw;

    //call subgraph algorithm on this data
    if(correlatedSetMode == 1)
    {
        consensusCorrelationGraph.calculate_connected_components(featureSetsRaw, minFeatureSetSize);
    }
    else if(correlatedSetMode == 0)
    {
        consensusCorrelationGraph.calculate_fully_connected_components(featureSetsRaw, minFeatureSetSize);
    }
    else if(correlatedSetMode >= 2)
    {
        consensusCorrelationGraph.calculate_min_edge_connected_components(featureSetsRaw, minFeatureSetSize, correlatedSetMode);
    }
    else
    {
        LOCO_OUT << "Invalid mode for detection of correlated sets: fallback to connected component -> mode 0\n";
        consensusCorrelationGraph.calculate_connected_components(featureSetsRaw, minFeatureSetSize);
    } 

    //iterate through the list of feature-pairs: for each pair check in which featureSetsRaw they are and create the string for it

    size_t numPairs = neighbourhoodCorrs.pairList.size();

    // 1. Setup fast lookup map (Pair -> Index in pairList)
    std::unordered_map<std::pair<int, int>, int, pair_hash> pairToIndexMap;
    pairToIndexMap.reserve(numPairs);    
    for (size_t i = 0; i < numPairs; ++i) {
        pairToIndexMap[neighbourhoodCorrs.pairList[i]] = i;
    }

    // 2. Storage for mapping pair index -> lists of feature sets that contain this pair
    std::vector<std::vector<std::vector<int>>> pairToSets(numPairs);

    // 3. The Smart Iteration: Loop sets, build pairs, map to index
    for (const auto& featureSet : featureSetsRaw) {
        
        if (featureSet.size() < 2) continue; // Skip sets too small to form a pair

        for (size_t i = 0; i < featureSet.size(); ++i) {
            for (size_t j = i + 1; j < featureSet.size(); ++j) {
                
                int u = featureSet[i];
                int v = featureSet[j];

                // Ensure the pair is always ordered [smaller, larger] to match pairList
                if (u > v) {
                    std::swap(u, v);
                }

                // Check if this generated pair exists in our main pairList
                auto it = pairToIndexMap.find({u, v});
                if (it != pairToIndexMap.end()) {
                    int pairIndex = it->second;
                    // Add this entire featureSet to the pair's storage
                    pairToSets[pairIndex].push_back(featureSet);
                }
            }
        }
    }

    // 4. Create the final vector of strings aligned exactly with pairList
    neighbourhoodCorrs.featureSetString.resize(numPairs);
    
    for (size_t i = 0; i < numPairs; ++i) {
        // Call your custom function, passing the collected sets for this specific pair
        neighbourhoodCorrs.featureSetString[i] = consensusCorrelationGraph.return_named_feature_set(pairToSets[i]);
    }

    // Now finalFeatureSetStrings[i] contains the correctly formatted string 
    // for the pair at neighbourhoodCorrs.pairList[i]
    
}

//centralneighbourhoodPtrs: all the neighbourhoodsPTrs, sotring all cells in enighbourhoods etc.
void Neighborhood::calculate_correlation_propagation(double correlationStrengthCutoff, int minFeatureSetSize, const bool calcSets, int threads,
                                                     double minNonNAFraction, unsigned int numberCorrelations, double significanceCutoff,
                                                     double minVaryingFraction)
{
    this->minVaryingFraction = minVaryingFraction;
    if(minVaryingFraction > 0.0)
    {
        LOCO_OUT << "minVaryingFraction: correlations of a feature are NA in neighbourhoods where less than "
                 << minVaryingFraction * 100 << "% of cells differ from its most frequent value\n";
    }

    LOCO_OUT << "STEP 2";
    LOCO_OUT << "\tCalculate correlations in all neighbourhoods\n";
    //output: neighbourhoodCorrs now stores all correlations across neighbourhoods (outer vector is feature-pairs for quick L calculation)
    step_2_calculate_correlation(correlationStrengthCutoff, threads);
    // pairs with too few defined correlations are removed BEFORE scoring (no score, not in the result)
    step_2b_filter_pairs_by_non_na(minNonNAFraction);

    LOCO_OUT << "STEP 3";
    LOCO_OUT << "\tCalculate Laplacian score\n";
    // also removes the pairs with an undefined score
    step_3_calculate_laplacian_score(threads);

    // p_adj over all tested pairs, then the numberCorrelations filter: feature sets are only calculated for the kept pairs
    step_3b_adjust_p_values();
    step_3b_filter_pairs(numberCorrelations, significanceCutoff);

    if(calcSets)
    {
        LOCO_OUT << "STEP 4";
        LOCO_OUT << "\tCalcualte correlated feature sets (biological programs)";
        std::vector<std::vector<int>> featureSetsRaw;
        step_4_calculate_feature_sets(minFeatureSetSize);
    }
}

void Neighborhood::create_neighborhood_graph(int knn)
{
    //create graph of neighborhoods: input is only the central nodes of neighborhoods
    std::cout << "Create neighbourhood graph\n";
    std::shared_ptr<GraphData> scGraphData = std::make_shared<GraphData>(centralNeighborhoodPtrs, cellStateGenes, knn, &GraphIni::cell_similarity_graph_manhattan_nodes);
    //inout radius does not matter, bandwidth=0 triggers auto-estimation via calc_bandwidth
    neighborhoodGraph = std::make_shared<GraphHandler>(scGraphData, knn, 0, 0);
    neighborhoodGraph->create_graph();
}

Neighborhood::Neighborhood(const std::shared_ptr<const GraphData> scData, unsigned int neighborhoodNumber, 
                           unsigned int neighborhoodSize, int neighborhoodKNN,
                           const SingleCellData& inputData,
                           const std::vector<int>& cellStateGenes, const std::vector<int>& corrStateGenes, int permutations,
                           const double& corrSetAbundance, const unsigned int correlatedSetMode, const std::string& correlationType,
                           const std::string& nSampling,
                           const int threads, const unsigned int kmeansR, const double kmeansSlack,
                           const unsigned int seed) :
                           neighborhoodSize(neighborhoodSize), neighbourhoodNum(neighborhoodNumber), inputDataOrigional(inputData),
                           cellStateGenes(cellStateGenes), corrStateGenes(corrStateGenes), permutations(permutations),
                           minimumCorrSetAbundance(corrSetAbundance), correlatedSetMode(correlatedSetMode), correlationType(correlationType)
{
    if(nSampling == "kmeans")
    {
        // disjoint neighbourhoods: size-constrained k-medians clusters (Manhattan) in the cell state space.
        // anchor = cell closest to the cluster median, neighbourhood = all cells of the cluster
        const neighborhoodCalculations::KmediansResult clusters = neighborhoodCalculations::kmedians_with_repair(
            scData, neighborhoodNumber, neighborhoodSize, kmeansR, kmeansSlack, seed, threads);
        // why a cell is in no neighbourhood
        const char* labels[] = {"placed", "no_room", "dissolved", "dropped"};
        cellNeighbourhoodStatus.resize(clusters.status.size());
        for(size_t i = 0; i < clusters.status.size(); ++i) cellNeighbourhoodStatus[i] = labels[static_cast<int>(clusters.status[i])];
        for(const auto& cluster : clusters.groups)
        {
            const nodePtr centerNode = scData->get_node_at(cluster.first);
            centralNeighborhoodPtrs.push_back(centerNode);
            neighborhoods.insert(std::make_pair(centerNode, cluster.second));
        }
    }
    else if(nSampling == "random")
    {
        int cellIDRange = scData->number_of_nodes();
        //save all neighborhood IDs & node IDs making up neighborhoods
        std::vector<int> centralNodeIDs = neighborhoodCalculations::get_random_elements(neighborhoodNumber, cellIDRange, seed);
        for(int centerNodeID : centralNodeIDs)
        {
            const nodePtr centerNode = scData->get_node_at(centerNodeID);
            centralNeighborhoodPtrs.push_back(centerNode);
            //std::vector<int> value = scData->get_adjacent_node_ids_knn(centerNode, neighborhoodSize);
            std::vector<int> value = scData->get_adjacent_node_ids_knn_kdsearch(centerNode);

            neighborhoods.insert(std::make_pair(centerNode, value));
        }
    }
    else
    {
        throw std::invalid_argument("Nsampling must be 'random' or 'kmeans'");
    }
    // the number of neighbourhoods actually built (kmeans can dissolve neighbourhoods)
    neighbourhoodNum = centralNeighborhoodPtrs.size();

    // per cell: "placed" (in at least one neighbourhood) or why not (kmeans: no_room/ dissolved/ dropped, set above)
    if(cellNeighbourhoodStatus.empty())
    {
        cellNeighbourhoodStatus.assign(scData->number_of_nodes(), "not_in_neighbourhood");
        for(const auto& nb : neighborhoods)
        {
            for(const int cell : nb.second) cellNeighbourhoodStatus.at(cell) = "placed";
        }
    }

    //create the neighborhood graph (how to neighborhoods connect)
    // creates single-cell Graph and neighborhoods from that
    create_neighborhood_graph(neighborhoodKNN);
}

//the final R structure needs:
//return R-object of LoCo results
// create a List of several dataframes
// 1.a) raw data table
// 1.b) UMAP data table

// 2.a) neighbourhood - cells
// 2.b) laplacian scores
// 2c) shuffled results

// 3a) neighbourhood - coords
// 3b) neighbourhood correlations data table
void Neighborhood::fill_result_data(

    std::vector<std::string>& nIDs, // all neighborhoods IDs
    std::vector<std::string>& nID_anchorCellID, //achnor cell IDs for neighborhoods
    std::vector<std::vector<std::string>>& nID_allCellIDs, //vector off all cellIDs for all neighborhoods (same order as nIDs)

    std::vector<std::string>& correlation_pairs_sorted, //all names of the correlation pairs

    std::vector<std::string>& correlation_pairs_origional, //all names of the correlation pairs (order of get_pair_correlations)

    std::vector<double>& corrL, 
    std::vector<double>& pCorrL, 
    std::vector<double>& pAdjL,
    std::vector<std::string>& featureSetString,
    const bool calcFeatureSets
)
{

    // =========================
    // anchor cell IDs for neighbourhood IDs
    // =========================
    const std::vector<nodePtr> nodes = neighborhoodGraph->get_all_nodes();
    for(const auto& node : nodes)
    {
        std::string cellID = node->get_name();
        nIDs.push_back("N_" + cellID);
        nID_anchorCellID.push_back(cellID);
    }

    // =========================
    // LAPLACIAN RESULTS
    // =========================
    size_t numPairs = neighbourhoodCorrs.pairNames.size();
    std::vector<size_t> indices(numPairs);
    std::iota(indices.begin(), indices.end(), 0);
    // Sort by the rank of the pair when ALL pairs (before filtering) were sorted by laplacianScores.L in
    // increasing order (see step_3b_filter_pairs): unfiltered output = sorted by L, filtering keeps the order
    std::sort(indices.begin(), indices.end(), 
        [this](size_t a, size_t b) {
            return laplacianScores.sortRank[a] < laplacianScores.sortRank[b];
        }
    );
    // Use the sorted indices to populate your lists
    for(size_t i = 0; i < numPairs; ++i)
    {
        size_t idx = indices[i]; // Grab the index that belongs in this sorted position

        //store correlations in order for Laplacian
        std::string featurePairTmp = neighbourhoodCorrs.pairNames[idx];
        correlation_pairs_sorted.push_back(featurePairTmp);

        // values
        corrL.push_back(laplacianScores.L[idx]);
        pCorrL.push_back(laplacianScores.p_values[idx]);
        pAdjL.push_back(laplacianScores.p_adj[idx]);

        // get feature sets
        if(calcFeatureSets)
        {
            featureSetString.push_back(neighbourhoodCorrs.featureSetString[idx]);
        }
    }

    // =========================
    // CORRELATION MATRIX
    // =========================
    size_t numNeighborhoods = nIDs.size();

    //store correlation in order of the correlation matrix
    correlation_pairs_origional = neighbourhoodCorrs.pairNames;

    //check order is as expected
    for (size_t i = 0; i < numNeighborhoods; ++i)
    {
        std::string graphName = "N_" + neighborhoodGraph->get_node_at(i)->get_name();
        
        if (nIDs[i] != graphName)
        {
            throw std::runtime_error(
                "Order assumption failed! Matrix column " + std::to_string(i) + 
                " is '" + graphName + "', but nIDs expects '" + nIDs[i] + "'."
            );
        }
    }

    // the correlations themselves are not copied here: get_pair_correlations() returns them (row = pair in the
    // order of correlation_pairs_origional, col = neighbourhood in the order of nIDs, checked above)
    if(neighbourhoodCorrs.pairToListOfAllNCorrs.rows != correlation_pairs_origional.size() ||
       neighbourhoodCorrs.pairToListOfAllNCorrs.cols != numNeighborhoods)
    {
        throw std::runtime_error("Correlation matrix does not match the number of pairs/ neighbourhoods.");
    }

    // =========================
    // CELL IDs PER NEIGHBORHOOD
    // =========================
    std::vector<std::string> nIDs_cell_order;
    for(const auto& neighbourhoodCellList : neighborhoods)
    {
        std::vector<std::string> cellIDList;
        //every cellID is an int from the vector of ints for cell-positions in raw data
        for(const auto& cellID : neighbourhoodCellList.second)
        {
            //get the cellID at position cell and push into cellIDList 
            cellIDList.push_back(inputDataOrigional.cellIDs.at(cellID));
        }
        nID_allCellIDs.push_back(cellIDList);
        nIDs_cell_order.push_back("N_" + neighbourhoodCellList.first->get_name());
    }

    //reorder according to nIDs:
    // Build lookup: current order -> index
    std::unordered_map<std::string, size_t> indexMap_allCellNids;
    for (size_t i = 0; i < nIDs_cell_order.size(); ++i)
    {
        indexMap_allCellNids[nIDs_cell_order[i]] = i;
    }
    // Prepare reordered containers
    std::vector<std::vector<std::string>> reorderedCellIDs;
    reorderedCellIDs.reserve(nIDs.size());
    // Reorder according to nIDs
    for (const auto& nID : nIDs)
    {
        auto it = indexMap_allCellNids.find(nID);
        if (it == indexMap_allCellNids.end())
        {
            throw std::runtime_error("nID not found in cellID data: " + nID);
        }
        size_t idx = it->second;
        reorderedCellIDs.push_back(nID_allCellIDs[idx]);
    }
    // Replace originals
    nID_allCellIDs = std::move(reorderedCellIDs);
}
