import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import matplotlib.ticker as ticker

# --- CẤU HÌNH BÀI TOÁN ---
DIM = 5                    
POP_SIZE = 100              
MAX_GEN = 1000              
BOUNDS = [-5.12, 5.12]      

# --- HÀM MỤC TIÊU (RASTRIGIN) ---
def rastrigin(X):
    A = 10
    return A * DIM + np.sum(X**2 - A * np.cos(2 * np.pi * X), axis=1)

# --- CÁC HÀM THUẬT TOÁN (GIỮ NGUYÊN LOGIC CŨ) ---
# Mình gom gọn lại để tập trung vào phần vẽ biểu đồ
def run_de_rand_1_bin(max_gen, pop_size, dim, bounds):
    F = 0.5; CR = 0.9
    pop = np.random.uniform(bounds[0], bounds[1], (pop_size, dim))
    fitness = rastrigin(pop)
    best_fitness_history = [np.min(fitness)]
    for g in range(max_gen):
        new_pop = np.copy(pop)
        for i in range(pop_size):
            idxs = [idx for idx in range(pop_size) if idx != i]
            r1, r2, r3 = np.random.choice(idxs, 3, replace=False)
            v = pop[r1] + F * (pop[r2] - pop[r3])
            cross_points = np.random.rand(dim) < CR
            if not np.any(cross_points): cross_points[np.random.randint(0, dim)] = True
            u = np.where(cross_points, v, pop[i]); u = np.clip(u, bounds[0], bounds[1])
            f_u = rastrigin(u.reshape(1, -1))[0]
            if f_u < fitness[i]: new_pop[i] = u; fitness[i] = f_u
        pop = new_pop; best_fitness_history.append(np.min(fitness))
    return best_fitness_history

def run_de_current_to_best_1(max_gen, pop_size, dim, bounds):
    F = 0.5; CR = 0.9
    pop = np.random.uniform(bounds[0], bounds[1], (pop_size, dim))
    fitness = rastrigin(pop)
    best_fitness_history = [np.min(fitness)]
    for g in range(max_gen):
        new_pop = np.copy(pop)
        best_idx = np.argmin(fitness); x_best = pop[best_idx]
        for i in range(pop_size):
            idxs = [idx for idx in range(pop_size) if idx != i]
            r1, r2 = np.random.choice(idxs, 2, replace=False)
            v = pop[i] + F * (x_best - pop[i]) + F * (pop[r1] - pop[r2])
            cross_points = np.random.rand(dim) < CR
            if not np.any(cross_points): cross_points[np.random.randint(0, dim)] = True
            u = np.where(cross_points, v, pop[i]); u = np.clip(u, bounds[0], bounds[1])
            f_u = rastrigin(u.reshape(1, -1))[0]
            if f_u < fitness[i]: new_pop[i] = u; fitness[i] = f_u
        pop = new_pop; best_fitness_history.append(np.min(fitness))
    return best_fitness_history

def run_shade(max_gen, pop_size, dim, bounds):
    H_SIZE = 100; M_CR = np.ones(H_SIZE) * 0.5; M_F = np.ones(H_SIZE) * 0.5; k = 0
    pop = np.random.uniform(bounds[0], bounds[1], (pop_size, dim))
    fitness = rastrigin(pop)
    best_fitness_history = [np.min(fitness)]
    for g in range(max_gen):
        new_pop = np.copy(pop); new_fitness = np.copy(fitness)
        S_CR = []; S_F = []; S_diff = []; sorted_indices = np.argsort(fitness)
        for i in range(pop_size):
            r_idx = np.random.randint(0, H_SIZE)
            cr_i = np.clip(np.random.normal(M_CR[r_idx], 0.1), 0, 1)
            while True:
                f_i = stats.cauchy.rvs(loc=M_F[r_idx], scale=0.1); 
                if f_i > 0: break
            f_i = min(f_i, 1.0)
            p_limit = max(2, int(pop_size * 0.1)); pbest_idx = np.random.choice(sorted_indices[:p_limit]); x_pbest = pop[pbest_idx]
            idxs = [idx for idx in range(pop_size) if idx != i]; r1, r2 = np.random.choice(idxs, 2, replace=False)
            v = pop[i] + f_i * (x_pbest - pop[i]) + f_i * (pop[r1] - pop[r2])
            cross_points = np.random.rand(dim) < cr_i; cross_points[np.random.randint(0, dim)] = True
            u = np.where(cross_points, v, pop[i]); u = np.clip(u, bounds[0], bounds[1])
            f_u = rastrigin(u.reshape(1, -1))[0]
            if f_u < fitness[i]:
                diff = fitness[i] - f_u; new_pop[i] = u; new_fitness[i] = f_u
                S_CR.append(cr_i); S_F.append(f_i); S_diff.append(diff)
        pop = new_pop; fitness = new_fitness
        if len(S_CR) > 0:
            weights = np.array(S_diff) / np.sum(S_diff)
            M_CR[k] = np.sum(weights * np.array(S_CR))
            M_F[k] = np.sum(weights * (np.array(S_F)**2)) / np.sum(weights * np.array(S_F))
            k = (k + 1) % H_SIZE
        best_fitness_history.append(np.min(fitness))
    return best_fitness_history

# --- MAIN PROGRAM ---
print(f"Running comparison (D={DIM})...")
loss_rand = run_de_rand_1_bin(MAX_GEN, POP_SIZE, DIM, BOUNDS)
loss_ctb = run_de_current_to_best_1(MAX_GEN, POP_SIZE, DIM, BOUNDS)
loss_shade = run_shade(MAX_GEN, POP_SIZE, DIM, BOUNDS)

# --- XỬ LÝ DỮ LIỆU ĐỂ GIẢ LẬP SỐ 0 ---
ZERO_THRESHOLD = 1e-15 # Ngưỡng cực nhỏ để coi là 0 trên biểu đồ

def clip_for_log_plot(data, threshold):
    # Thay thế các giá trị quá nhỏ bằng threshold để vẽ không bị lỗi
    return np.maximum(np.array(data), threshold)

loss_rand_plot = clip_for_log_plot(loss_rand, ZERO_THRESHOLD)
loss_ctb_plot = clip_for_log_plot(loss_ctb, ZERO_THRESHOLD)
loss_shade_plot = clip_for_log_plot(loss_shade, ZERO_THRESHOLD)

# --- VẼ BIỂU ĐỒ ---
plt.figure(figsize=(10, 6))

# Vẽ các đường
plt.semilogy(loss_rand_plot, label='DE/rand/1/bin', linewidth=1.5, linestyle='--', color='tab:blue')
plt.semilogy(loss_ctb_plot, label='DE/current-to-best/1/bin', linewidth=1.5, linestyle='-.', color='tab:orange')
plt.semilogy(loss_shade_plot, label='SHADE (Adaptive)', linewidth=2.0, color='red')

# Cấu hình trục
plt.title(f'So sánh tốc độ hội tụ trên hàm Rastrigin (D={DIM})', fontsize=14)
plt.xlabel('Thế hệ (Generation)', fontsize=12)
plt.ylabel('Giá trị Fitness tối ưu', fontsize=12)

# Grid
plt.grid(True, which="major", ls="-", alpha=0.4)
plt.grid(True, which="minor", ls=":", alpha=0.2)

# --- TÙY CHỈNH TRỤC Y ĐỂ HIỆN SỐ 0 GIẢ ---
ax = plt.gca()

# Thiết lập các mốc hiển thị (Ticks) thủ công
# Bạn muốn hiện: 10^2, 10^1, 10^0 (1), và số 0 (thực tế là ZERO_THRESHOLD)
major_ticks = [100, 10, 1, ZERO_THRESHOLD] 
ax.set_yticks(major_ticks)

# Đặt lại nhãn (Labels) cho các mốc đó. Chỗ ZERO_THRESHOLD ta ghi là "0"
ax.set_yticklabels(['$10^2$', '$10^1$', '$10^0$', '0'])

# Giới hạn trục Y để đẹp
plt.ylim(bottom=ZERO_THRESHOLD * 0.5) 

plt.legend(fontsize=12)
plt.tight_layout()
plt.show()

print(f"Final Rand: {loss_rand[-1]:.5e}")
print(f"Final C-to-B: {loss_ctb[-1]:.5e}")
print(f"Final SHADE: {loss_shade[-1]:.5e}")