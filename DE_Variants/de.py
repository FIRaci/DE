import numpy as np

# --- ĐỊNH NGHĨA BÀI TOÁN ---
# Đây là hàm mục tiêu (fitness function)
# Chính là hàm f(x) = x1^2 + x2^2 từ file PDF
def fitness_function(solution):
    """Tính toán giá trị fitness của một cá thể (giải pháp)."""
    return np.sum(solution**2)

# --- THUẬT TOÁN DI TRUYỀN VI SAI (DE) ---
def differential_evolution(fitness_func, lower_bound, upper_bound, pop_size, dim, F, CR, max_gen):
    """
    Thực thi thuật toán DE/rand/1/bin.
    
    Tham số:
    - fitness_func: Hàm mục tiêu cần tối ưu (minimize)
    - lower_bound: Biên dưới của không gian tìm kiếm
    - upper_bound: Biên trên của không gian tìm kiếm
    - pop_size (N): Kích thước quần thể
    - dim (D): Số chiều của bài toán
    - F: Hệ số co giãn (Scaling Factor)
    - CR (P_cr): Xác suất lai ghép (Crossover Probability)
    - max_gen: Số thế hệ tối đa (điều kiện dừng)
    """
    
    # --- Bước 1: Khởi tạo Quần thể (Initialization) ---
    # Tạo ngẫu nhiên quần thể P gồm N cá thể D chiều
    # Giá trị nằm trong khoảng [lower_bound, upper_bound]
    pop = lower_bound + (upper_bound - lower_bound) * np.random.rand(pop_size, dim)
    
    # Tính toán fitness cho từng cá thể trong quần thể ban đầu
    fitness = np.array([fitness_func(ind) for ind in pop])
    
    # Tìm cá thể tốt nhất ban đầu
    best_idx = np.argmin(fitness)
    best_solution = pop[best_idx].copy()
    best_fitness = fitness[best_idx]
    
    print(f"--- BẮT ĐẦU THUẬT TOÁN DE ---")
    print(f"Quần thể ban đầu, Best Fitness: {best_fitness:.6f}")
    
    # --- Bắt đầu vòng lặp các Thế hệ (Generations) ---
    for g in range(max_gen):
        
        # Lặp qua từng cá thể (Parent Vector) trong quần thể
        for i in range(pop_size):
            parent_vector = pop[i]
            
            # --- Bước 2: Đột biến (Mutation) ---
            # Chọn ngẫu nhiên 3 chỉ số r1, r2, r3
            # Chúng phải khác nhau và khác với chỉ số i hiện tại
            indices = [idx for idx in range(pop_size) if idx != i]
            r1, r2, r3 = np.random.choice(indices, 3, replace=False)
            
            # Lấy 3 vector tương ứng
            x_r1 = pop[r1]
            x_r2 = pop[r2]
            x_r3 = pop[r3]
            
            # Áp dụng công thức đột biến (DE/rand/1)
            # v_i = x_r1 + F * (x_r2 - x_r3)
            # (File PDF gọi đây là Trial Vector u_i)
            trial_vector = x_r1 + F * (x_r2 - x_r3)
            
            # Xử lý biên (Bounds Handling) - Rất quan trọng!
            # Nếu giá trị nào vượt biên, "kéo" nó về biên gần nhất
            trial_vector = np.clip(trial_vector, lower_bound, upper_bound)
            
            # --- Bước 3: Lai ghép (Crossover - Binomial) ---
            # Tạo Cá thể con (Offspring)
            offspring_vector = parent_vector.copy()
            
            # Chọn 1 chiều j_rand BẮT BUỘC thay đổi
            j_rand = np.random.randint(0, dim)
            
            for j in range(dim):
                # Kiểm tra điều kiện lai ghép
                if np.random.rand() < CR or j == j_rand:
                    # Lấy gen từ Trial Vector
                    offspring_vector[j] = trial_vector[j]
            
            # --- Bước 4: Lựa chọn (Selection) ---
            # Tính fitness của "con" (Offspring)
            offspring_fitness = fitness_func(offspring_vector)
            
            # So sánh "con" và "cha" (Parent)
            if offspring_fitness < fitness[i]:
                # "Con" tốt hơn -> thay thế "cha"
                pop[i] = offspring_vector
                fitness[i] = offspring_fitness
                
                # Cập nhật giải pháp tốt nhất (best_solution) nếu cần
                if offspring_fitness < best_fitness:
                    best_fitness = offspring_fitness
                    best_solution = offspring_vector.copy()
        
        # Thông báo kết quả cuối mỗi thế hệ
        print(f"Thế hệ {g+1}/{max_gen}, Best Fitness: {best_fitness:.6f}")

    print(f"--- KẾT THÚC ---")
    return best_solution, best_fitness

# --- CHẠY THUẬT TOÁN ---
if __name__ == "__main__":
    # Đặt các tham số giống như trong file DE.pdf
    POP_SIZE = 5       # N (Kích thước quần thể)
    DIMENSION = 2      # D (Số chiều)
    LOWER_BOUND = -5   # Biên dưới
    UPPER_BOUND = 5    # Biên trên
    F_FACTOR = 0.5     # F (Hệ số co giãn)
    CR_RATE = 0.7      # P_cr (Xác suất lai ghép)
    MAX_GEN = 20        # Số thế hệ tối đa
    
    # Chạy thuật toán
    final_solution, final_fitness = differential_evolution(
        fitness_function, 
        LOWER_BOUND, 
        UPPER_BOUND, 
        POP_SIZE, 
        DIMENSION, 
        F_FACTOR, 
        CR_RATE, 
        MAX_GEN
    )
    
    print(f"\nGiải pháp tốt nhất tìm được: {final_solution}")
    print(f"Giá trị Fitness tốt nhất: {final_fitness:.6f}")
    print(f"(Đáp án lý tưởng là [0, 0] với fitness = 0.0)")