#include "identification.h"

#include <algorithm>
#include <chrono>
#include <deque>
#include <memory>
#include <stdexcept>

namespace fasttreeid {

std::vector<std::pair<size_t, size_t>> missingEdges (size_t n,
		const std::vector<std::pair<size_t, size_t>> &original) {

	std::vector<std::vector<bool>> adj(n, std::vector<bool>(n, false));
	for (auto [u,v] : original) {
		adj[u][v] = true;
	}
	std::vector<std::pair<size_t, size_t>> result;
	for (size_t i = 0; i < n; i++) {
		for (size_t j = i + 1; j < n; j++) {
			if (!adj[i][j]) {
				result.emplace_back(i, j);
			}
		}
	}
	return result;
}

void generateTreks (const std::vector<std::vector<std::shared_ptr<Lambda>>> &lambda,
		const std::vector<std::vector<std::shared_ptr<Omega>>> &omega,
		const std::vector<std::vector<std::shared_ptr<Sigma>>> &sigma,
		const std::vector<size_t> &directed,
		int64_t prime) {

	size_t n = directed.size() + 1;

	/* First compute σ_{0,i} and σ_{i,0} for all i */
	sigma[0][0]->setInTermsOfLambdaAndOmega(omega[0][0]->getValue());

	for (size_t i = 1; i < n; i++) {
		size_t p = directed[i-1];
		int64_t result = multiply(sigma[0][p]->getInTermsOfLambdaAndOmega(), lambda[p][i]->getValue(), prime);
		if (omega[0][i]) {
			result += omega[0][i]->getValue();
			if (result >= prime) result -= prime;
		}
		sigma[0][i]->setInTermsOfLambdaAndOmega(result);
		sigma[i][0]->setInTermsOfLambdaAndOmega(result);
	}

	/* Now compute σ_{i,j} for 0 < i,j */
	for (size_t i = 1; i < n; i++) {
		size_t p = directed[i-1];
		for (size_t j = i; j < n; j++) {
			size_t q = directed[j-1];
			int64_t result = multiply(sigma[p][j]->getInTermsOfLambdaAndOmega(), lambda[p][i]->getValue(), prime);
			result += multiply(sigma[i][q]->getInTermsOfLambdaAndOmega(), lambda[q][j]->getValue(), prime);
			if (result >= prime) result -= prime;
			result += prime - multiply(sigma[p][q]->getInTermsOfLambdaAndOmega(), multiply(lambda[p][i]->getValue(), lambda[q][j]->getValue(), prime), prime);
			if (result >= prime) result -= prime;
			if (omega[i][j]) {
				result += omega[i][j]->getValue();
				if (result >= prime) result -= prime;
			}
			sigma[i][j]->setInTermsOfLambdaAndOmega(result);
			sigma[j][i]->setInTermsOfLambdaAndOmega(result);
		}
	}
}

int64_t computeRank (size_t i, size_t j, size_t p, size_t q,
		const std::shared_ptr<RandomTool> &randomTool,
		const std::vector<std::vector<std::shared_ptr<Sigma>>> &sigma) {

	/* Compute the rank of the edge {i,j} where p is the parent of i and q is the parent of j */

	int64_t prime = randomTool->getPrime();
	int64_t sigma_pq = sigma[p][q]->getInTermsOfLambdaAndOmega();
	int64_t sigma_iq = sigma[i][q]->getInTermsOfLambdaAndOmega();
	int64_t sigma_pj = sigma[p][j]->getInTermsOfLambdaAndOmega();
	int64_t sigma_ij = sigma[i][j]->getInTermsOfLambdaAndOmega();

	if (sigma_pq == 0 && sigma_iq == 0 && sigma_pj == 0 && sigma_ij == 0) {
		return 0;
	}
	if ((multiply(sigma_pq, sigma_ij, prime) - multiply(sigma_pj, sigma_iq, prime) + prime) % prime == 0) {
		return 1;
	}
	return 2;
}

void floodComponent (size_t u, std::vector<bool> &visited, std::vector<size_t> &component,
		const std::vector<std::vector<size_t>> &graph) {

	if (!visited[u]) {
		visited[u] = true;
		component.push_back(u);
		for (size_t v : graph[u]) {
			floodComponent(v, visited, component, graph);
		}
	}
}

void bfsPropagate (size_t u, size_t n, const std::vector<std::vector<size_t>> &graph,
		std::vector<std::optional<std::variant<Fraction, Cycle, Path>>> &identification) {

	if (u == 0 || identification[u]) {
		return;
	}
	std::vector<std::optional<size_t>> parent(n);
	std::deque<size_t> q;
	q.push_back(u);
	while (!q.empty()) {
		size_t cur = q.front();
		q.pop_front();
		if (identification[cur] && !std::holds_alternative<Path>(*identification[cur])) {
			std::vector<size_t> path;
			while (cur != u) {
				path.push_back(cur);
				cur = *parent[cur];
			}
			std::reverse(path.begin(), path.end());
			identification[u] = Path(path);
			return;
		}
		for (size_t v : graph[cur]) {
			if (!parent[v]) {
				parent[v] = cur;
				q.push_back(v);
			}
		}
	}
	throw std::logic_error("bfsPropagate: Failure to find a path. Please report this error.");
}

std::vector<std::vector<std::vector<std::vector<int64_t>>>> multiplyComponentwise
		(const std::vector<std::vector<std::vector<std::vector<int64_t>>>> &initial,
		const std::shared_ptr<RandomTool> &randomTool) {

	auto result = initial;
	for (auto &t0 : result) {
		for (auto &t1 : t0) {
			int64_t val = randomTool->getRandomSmallerPrime();
			for (auto &t2 : t1) {
				for (auto &t3 : t2){
					t3 = multiply(t3, val, randomTool->getPrime());
				}
			}
		}
	}
	return result;
}

std::vector<std::vector<int64_t>> matrixMultiply (const std::vector<std::vector<int64_t>> &left,
		const std::vector<std::vector<int64_t>> &right,
		int64_t prime) {

	std::vector<std::vector<int64_t>> result(2, std::vector<int64_t>(2));
	for (size_t i = 0; i < 2; i++) {
		for (size_t j = 0; j < 2; j++) {
			for (size_t k = 0; k < 2; k++) {
				result[i][j] += multiply(left[i][k], right[k][j], prime);
			}
			if (result[i][j] >= prime){
				result[i][j] -= prime;
			}
		}
	}
	return result;
}

std::vector<std::vector<std::vector<std::vector<int64_t>>>> matrixMultiply
		(const std::vector<std::vector<std::vector<std::vector<int64_t>>>> &left,
		const std::vector<std::vector<std::vector<std::vector<int64_t>>>> &right,
		int64_t prime) {

	size_t n = left.size(), m = right.size(), l = right[0].size();
	std::vector<std::vector<std::vector<std::vector<int64_t>>>> result(n,
			std::vector<std::vector<std::vector<int64_t>>>(l,
			std::vector<std::vector<int64_t>>(2, std::vector<int64_t>(2))));

	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < l; j++) {
			for (size_t k = 0; k < m; k++) {
				result[i][j][0][0] += multiply(left[i][k][0][0], right[k][j][0][0], prime) + multiply(left[i][k][0][1], right[k][j][1][0], prime);
				result[i][j][0][1] += multiply(left[i][k][0][0], right[k][j][0][1], prime) + multiply(left[i][k][0][1], right[k][j][1][1], prime);
				result[i][j][1][0] += multiply(left[i][k][1][0], right[k][j][0][0], prime) + multiply(left[i][k][1][1], right[k][j][1][0], prime);
				result[i][j][1][1] += multiply(left[i][k][1][0], right[k][j][0][1], prime) + multiply(left[i][k][1][1], right[k][j][1][1], prime);

				if (result[i][j][0][0] >= prime) {
					result[i][j][0][0] -= prime;
				}
				if (result[i][j][0][0] >= prime) {
					result[i][j][0][0] -= prime;
				}
				if (result[i][j][0][1] >= prime) {
					result[i][j][0][1] -= prime;
				}
				if (result[i][j][0][1] >= prime) {
					result[i][j][0][1] -= prime;
				}
				if (result[i][j][1][0] >= prime) {
					result[i][j][1][0] -= prime;
				}
				if (result[i][j][1][0] >= prime) {
					result[i][j][1][0] -= prime;
				}
				if (result[i][j][1][1] >= prime) {
					result[i][j][1][1] -= prime;
				}
				if (result[i][j][1][1] >= prime) {
					result[i][j][1][1] -= prime;
				}
			}
		}
	}
	return result;
}

std::vector<size_t> selfReducibility (size_t start, size_t n, size_t length, int64_t prime,
		const std::vector<std::vector<std::vector<std::vector<int64_t>>>> &initial) {

	std::vector<std::vector<std::vector<std::vector<std::vector<int64_t>>>>> suffix(length);

	/* suffix[length - 1] = e_start (start-th unit vector) */
	suffix[length - 1] = std::vector<std::vector<std::vector<std::vector<int64_t>>>>(n,
			std::vector<std::vector<std::vector<int64_t>>>(1,
				std::vector<std::vector<int64_t>>(2, std::vector<int64_t>(2))));
	suffix[length - 1][start][0] = {{1, 0}, {0, 1}};

	for (int i = length - 2; i >= 0; i--) { /* suffix products in O(n^3) */
		suffix[i] = matrixMultiply(initial, suffix[i + 1], prime);
	}

	std::vector<size_t> result;
	result.push_back(start);
	std::vector<std::vector<int64_t>> left = {{1, 0}, {0, 1}};
	size_t cur = start;

	/* Construct the path step by step, currently, the path goes start..cur, its weight is left */
	for (size_t t = 0; t < length - 1; t++) {
		int found = -1;
		for (size_t v = 0; v < n; v++) {
			/* It doesn't matter if the edge doesn't exist, then, the result will be 0 */
			auto tmp = matrixMultiply(left, matrixMultiply(initial[cur][v], suffix[t][v][0], prime), prime);
			if (tmp[0][1] != 0 || tmp[0][0] != tmp[1][1]) {
				found = v;
				break;
			} else if (tmp[1][0]) {
				throw std::logic_error("selfReducibility: Unexpected result. Please report this error.");
			}
		}
		if (found < 0) {
			throw std::logic_error("selfReducibility: Identifying cycle not found. Please report this error.");
		}
		left = matrixMultiply(left, initial[cur][found], prime);
		result.push_back(found);
		cur = found;
	}
	return result;
}

void checkSolutions (int64_t a, int64_t b, int64_t c, int64_t discriminant,
		int64_t prime, size_t u, size_t n, const std::vector<size_t> &component,
		const std::vector<std::vector<size_t>> &graph,
		const std::vector<size_t> &directed,
		const std::vector<std::vector<std::shared_ptr<Sigma>>> &sigma,
		std::vector<std::optional<std::variant<Fraction, Cycle, Path>>> &identification) {
	
	std::vector<std::optional<FASTP>> fastp(n);
	fastp[u] = FASTP((b - c + prime) % prime, 1, 2 * a % prime, discriminant, 0);
	std::deque<size_t> queue;
	queue.push_back(u);
	while (!queue.empty()) {
		size_t j = queue.front();
		queue.pop_front();
		for (size_t i : graph[j]) {
			if (!fastp[i]) {
				size_t p = directed[i-1], q = directed[j-1];
				int64_t sigma_pq = sigma[p][q]->getInTermsOfLambdaAndOmega();
				int64_t sigma_iq = sigma[i][q]->getInTermsOfLambdaAndOmega();
				int64_t sigma_pj = sigma[p][j]->getInTermsOfLambdaAndOmega();
				int64_t sigma_ij = sigma[i][j]->getInTermsOfLambdaAndOmega();
				fastp[i] = FASTP((multiply(fastp[j]->getP(), sigma_iq, prime) - multiply(fastp[j]->getR(), sigma_ij, prime) + prime) % prime,
							(multiply(fastp[j]->getQ(), sigma_iq, prime) - multiply(fastp[j]->getT(), sigma_ij, prime) + prime) % prime,
							(multiply(fastp[j]->getP(), sigma_pq, prime) - multiply(fastp[j]->getR(), sigma_pj, prime) + prime) % prime,
							fastp[j]->getS(),
							(multiply(fastp[j]->getQ(), sigma_pq, prime) - multiply(fastp[j]->getT(), sigma_pj, prime) + prime) % prime);

				if (!(fastp[i]->getR() != 0 || fastp[i]->getT() != 0)) {/* Otherwise, divisor = 0 */
					throw std::logic_error("checkSolutions: Divisor of FASTP is 0. Please report this error.");
				}

				/* First check whether the divisor is 0 for one of the options for s */
				if (multiply(fastp[i]->getR(), fastp[i]->getR(), prime) == multiply(multiply(fastp[i]->getT(), fastp[i]->getT(), prime), fastp[i]->getS(), prime)) {
					std::get<Cycle>(*identification[u]).setOneIdentifiableAZero(p, i);
					return;
				}
				queue.push_back(i);
			}
		}
	}
	for (size_t i : component) {
		for (size_t j : graph[i]) {
			size_t p = directed[i-1], q = directed[j-1];
			int64_t sigma_pq = sigma[p][q]->getInTermsOfLambdaAndOmega();
			int64_t sigma_iq = sigma[i][q]->getInTermsOfLambdaAndOmega();
			int64_t sigma_pj = sigma[p][j]->getInTermsOfLambdaAndOmega();
			int64_t sigma_ij = sigma[i][j]->getInTermsOfLambdaAndOmega();

			int64_t A = getA(*fastp[i], *fastp[j], prime, sigma_pq, sigma_iq, sigma_pj, sigma_ij);
			int64_t B = getB(*fastp[i], *fastp[j], prime, sigma_pq, sigma_iq, sigma_pj, sigma_ij);

			/* It must hold that A + B * sqrt(s) = 0 */
			if (A != 0 || B != 0) {
				if (!(A != 0 && B != 0)) { /* Otherwise no solution */
					throw std::logic_error("checkSolutions: No solution exists. Please report this error.");
				}
				if (!(multiply(A, A, prime) == multiply(multiply(B, B, prime), fastp[i]->getS(), prime))) {
					throw std::logic_error("checkSolutions: FASTPs are not equal. Please report this error.");
				}
				std::get<Cycle>(*identification[u]).setOneIdentifiableOneOption(i, j);
				return;
			}
		}
	}
}

void findIdentifyingCycle (size_t n, int64_t prime,
		const std::vector<size_t> &component,
		const std::vector<std::vector<std::vector<std::vector<int64_t>>>> &initial,
		const std::vector<std::vector<size_t>> &graph,
		const std::vector<size_t> &directed,
		const std::vector<std::vector<std::shared_ptr<Sigma>>> &sigma,
		std::vector<std::optional<std::variant<Fraction, Cycle, Path>>> &identification) {

	auto result = initial;
	for (size_t t = 2; t <= component.size(); t++) {
		result = matrixMultiply(result, initial, prime);
		for (size_t i = 0; i < component.size(); i++) {
			if (result[i][i][0][1] != 0 || result[i][i][0][0] != result[i][i][1][1]) {
				/* No need to transpose here, doesn't change whether it's a multiple of Id */
				auto cycle = selfReducibility(i, component.size(), t, prime, initial);
				for (size_t &v : cycle) {
					v = component[v];
				}

				identification[component[i]] = Cycle(cycle);
				for (size_t v : component) {
					bfsPropagate(v, n, graph, identification);
				}

				int64_t a = result[i][i][0][1], b = result[i][i][0][0], c = result[i][i][1][1],
					d = result[i][i][1][0];
				/* Lemma 8, a and d swapped because of transposition */

				if (a == 0) { /* Case 2 in lemma 8, a = 0 */
					std::get<Cycle>(*identification[component[i]]).setOneIdentifiableAZero(directed[component[i]-1], component[i]);
				} else {
					int64_t discriminant = (multiply((c - b + prime) % prime, (c - b + prime) % prime, prime) + 4 * multiply(a, d, prime) % prime) % prime;
					if (discriminant == 0) {
						std::get<Cycle>(*identification[component[i]]).setOneIdentifiableDiscriminantZero();
					} else {
						checkSolutions(a, b, c, discriminant, prime, component[i], n, component, graph, directed, sigma, identification);
					}
				}
				return;
			} else if (result[i][i][1][0]) {
				throw std::logic_error("findIdentifyingCycle: No solution. Please report this error.");
				/* Otherwise, there would be no solution (case 3 in lemma 8) */
			}
		}
	}
}

std::vector<std::vector<std::vector<std::vector<int64_t>>>> matrixExponentiate (const std::vector<std::vector<std::vector<std::vector<int64_t>>>> &matrix, size_t exponent, int64_t prime) {
	if (exponent == 1) {
		return matrix;
	}
	auto tmp = matrixExponentiate (matrix, exponent / 2, prime);
	tmp = matrixMultiply(tmp, tmp, prime);
	if (exponent % 2 == 1) {
		tmp = matrixMultiply(tmp, matrix, prime);
	}
	return tmp;
}

bool hasNonIdentityEntry (const std::vector<std::vector<std::vector<std::vector<int64_t>>>> &matrix, size_t &which) {
	for (size_t i = 0; i < matrix.size(); i++) {
		/* No need to transpose here, doesn't change whether it's a multiple of Id */
		if (matrix[i][i][0][1] != 0 || matrix[i][i][0][0] != matrix[i][i][1][1]) {
			which = i;
			return true;
		} else if (matrix[i][i][1][0]) {
			throw std::logic_error("hasNonIdentityEntry: No solution. Please report this error.");
			/* Otherwise, there would be no solution (case 3 in lemma 8) */
		}
	}
	return false;
}

void findIdentifyingCycleFaster (size_t n, int64_t prime,
		const std::vector<size_t> &component,
		const std::vector<std::vector<std::vector<std::vector<int64_t>>>> &initial,
		const std::vector<std::vector<size_t>> &graph,
		const std::vector<size_t> &directed,
		const std::vector<std::vector<std::shared_ptr<Sigma>>> &sigma,
		std::vector<std::optional<std::variant<Fraction, Cycle, Path>>> &identification) {

	size_t which = component.size();
	if (!hasNonIdentityEntry(matrixExponentiate(initial, component.size(), prime), which)) {
		return;
	}
	size_t low = 2, high = component.size(); /* TODO: can we start this at 3? */
	while (low < high) {
		size_t mid = (low + high) / 2;
		if (hasNonIdentityEntry(matrixExponentiate(initial, mid, prime), which)) {
			high = mid;
		} else {
			low = mid + 1;
		}
	}
	auto result = matrixExponentiate(initial, low, prime);
	if (!hasNonIdentityEntry(result, which)) {
		throw std::logic_error("findIdentifyingCycleFaster: No solution. Please report this error.");
	}
	auto cycle = selfReducibility(which, component.size(), low, prime, initial);
	for (size_t &v : cycle) {
		v = component[v];
	}

	identification[component[which]] = Cycle(cycle);
	for (size_t v : component) {
		bfsPropagate(v, n, graph, identification);
	}

	int64_t a = result[which][which][0][1], b = result[which][which][0][0], c = result[which][which][1][1],
	d = result[which][which][1][0];
	/* Lemma 8, a and d swapped because of transposition */

	if (a == 0) { /* Case 2 in lemma 8, a = 0 */
		std::get<Cycle>(*identification[component[which]]).setOneIdentifiableAZero(directed[component[which]-1], component[which]);
	} else {
		int64_t discriminant = (multiply((c - b + prime) % prime, (c - b + prime) % prime, prime) + 4 * multiply(a, d, prime) % prime) % prime;
		if (discriminant == 0) {
			std::get<Cycle>(*identification[component[which]]).setOneIdentifiableDiscriminantZero();
		} else {
			checkSolutions(a, b, c, discriminant, prime, component[which], n, component, graph, directed, sigma, identification);
		}
	}
}

void identifyComponent (size_t u, size_t n, std::vector<bool> &visited,
		const std::shared_ptr<RandomTool> &randomTool,
		const std::vector<std::vector<size_t>> &graph,
		const std::vector<size_t> &directed,
		const std::vector<std::vector<std::shared_ptr<Sigma>>> &sigma,
		std::vector<std::optional<std::variant<Fraction, Cycle, Path>>> &identification) {

	int64_t prime = randomTool->getPrime();
	std::vector<size_t> component;
	floodComponent(u, visited, component, graph);
	bool hasMarked = false;
	for (size_t v : component) {
		if (identification[v]) {
			hasMarked = true;
		}
	}
	if (hasMarked) {
		for (size_t v : component) {
			bfsPropagate(v, n, graph, identification);
		}
		return;
	}

	std::vector<size_t> rev_component(n);
	for (size_t i = 0; i < component.size(); i++) {
		rev_component[component[i]] = i;
	}
	std::vector<std::vector<std::vector<std::vector<int64_t>>>> initial(component.size(),
			std::vector<std::vector<std::vector<int64_t>>>(component.size()));

	for (size_t i = 0; i < component.size(); i++) {
		for (size_t j = 0; j < component.size(); j++) {
			initial[i][j] = {{0, 0}, {0, 0}};
		}
		for (size_t j : graph[component[i]]) {
			size_t p = directed[component[i]-1], q = directed[j-1];
			int64_t sigma_pq = sigma[p][q]->getInTermsOfLambdaAndOmega();
			int64_t sigma_iq = sigma[component[i]][q]->getInTermsOfLambdaAndOmega();
			int64_t sigma_pj = sigma[p][j]->getInTermsOfLambdaAndOmega();
			int64_t sigma_ij = sigma[component[i]][j]->getInTermsOfLambdaAndOmega();

			initial[i][rev_component[j]] = {{sigma_pj, sigma_pq},
				{(prime - sigma_ij) % prime, (prime - sigma_iq) % prime}};
			/* The matrix is transposed to multiply from the left instead of from the right */

			/* Multiply with a random value (indeterminate for this edge) */
			int64_t val = randomTool->getRandomSmallerPrime();
			initial[i][rev_component[j]][0][0] = multiply(initial[i][rev_component[j]][0][0], val, prime);
			initial[i][rev_component[j]][0][1] = multiply(initial[i][rev_component[j]][0][1], val, prime);
			initial[i][rev_component[j]][1][0] = multiply(initial[i][rev_component[j]][1][0], val, prime);
			initial[i][rev_component[j]][1][1] = multiply(initial[i][rev_component[j]][1][1], val, prime);
		}
		/* Now add the self loop */
		int64_t val = randomTool->getRandomSmallerPrime();
		initial[i][i][0][0] = val;
		initial[i][i][1][1] = val;
	}

	if (component.size() < 50) {
		findIdentifyingCycle(n, prime, component, initial, graph, directed, sigma, identification);
	} else {
		findIdentifyingCycleFaster(n, prime, component, initial, graph, directed, sigma, identification);
	}
}

uint64_t generateSeed () {
	std::random_device rd;
	uint64_t rd_seed = (static_cast<uint64_t>(rd()) << 32) | rd(); /* may be deterministic on some platforms */
	uint64_t time_seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	return rd_seed ^ time_seed; /* so xor with the current time to get something "non-deterministic" on all platforms */
}

IdentificationResult identify (const std::vector<std::pair<size_t, size_t>> &bidirected,
		const std::vector<size_t> &directed, std::optional<uint64_t> seed_opt, std::optional<int64_t> prime_opt) {

	uint64_t seed = seed_opt.has_value() ? seed_opt.value() : generateSeed();
	auto randomTool = std::make_shared<RandomTool>(seed, prime_opt);
	int64_t prime = randomTool->getPrime();

	if (!((1ll << 58) <= prime && prime < (1ll << 59))) {
		throw std::invalid_argument("The prime should be between 2^58 and 2^59");
	}

	size_t n = directed.size() + 1;
	std::vector<std::vector<std::shared_ptr<Lambda>>> lambda(n, std::vector<std::shared_ptr<Lambda>>(n));
	std::vector<std::vector<std::shared_ptr<Omega>>> omega(n, std::vector<std::shared_ptr<Omega>>(n));
	std::vector<std::vector<std::shared_ptr<Sigma>>> sigma(n, std::vector<std::shared_ptr<Sigma>>(n));

	for (size_t i = 1; i < n; i++) {
		lambda[directed[i - 1]][i] = std::make_shared<Lambda>(randomTool, directed[i - 1], i);
	}
	for (auto [u, v] : bidirected) {
		omega[u][v] = std::make_shared<Omega>(randomTool, u, v);
		omega[v][u] = omega[u][v];
	}
	for (size_t i = 0; i < n; i++) {
		omega[i][i] = std::make_shared<Omega>(randomTool, i, i);
	}
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			sigma[i][j] = std::make_shared<Sigma>(i, j);
		}
	}
	generateTreks(lambda, omega, sigma, directed, prime);

	auto missing = missingEdges(n, bidirected);
	std::vector<std::pair<size_t, size_t>> involvingRoot, rank2;
	for (auto [i,j] : missing) {
		if (i == 0) {
			involvingRoot.emplace_back(i, j);
		} else if (computeRank(i, j, directed[i - 1], directed[j - 1], randomTool, sigma) == 2) {
			rank2.emplace_back(i, j);
		}
	}

	IdentificationResult identificationResult(n, seed, prime);
	for (auto [i,j] : involvingRoot) {
		identificationResult.identification[j] = Fraction(sigma[0][j], sigma[0][directed[j - 1]]);
	}

	std::vector<std::vector<size_t>> graph(n);
	for (auto [i,j] : rank2) {
		graph[i].push_back(j);
		graph[j].push_back(i);
	}
	std::vector<bool> visited(n);
	for (size_t i = 1; i < n; i++) {
		if (!visited[i]) {
			identifyComponent(i, n, visited, randomTool, graph, directed, sigma, identificationResult.identification);
		}
	}

	return identificationResult;
}

}
