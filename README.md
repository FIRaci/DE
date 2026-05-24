# Giải thuật Tiến hóa Vi phân (Differential Evolution - DE) & Ứng dụng trong Cell-free ISAC

Dự án này là tập hợp hai nghiên cứu chính liên quan đến giải thuật Tiến hóa Vi phân (Differential Evolution - DE).

## Cấu trúc Dự án

### 1. `DE_Variants/` (So sánh các biến thể của DE)
Phần này giải quyết một bài toán thử nghiệm đơn giản nhằm mục đích đánh giá và so sánh hiệu năng hội tụ giữa các biến thể khác nhau của giải thuật Tiến hóa Vi phân (DE).
- **`de.py`**: Chứa các lớp/hàm triển khai giải thuật DE.
- **`compare_de.py`**: Kịch bản chạy mô phỏng và vẽ đồ thị so sánh giữa các biến thể DE.

### 2. `Cell_Free_ISAC/` (Ứng dụng DE vào Cell-free ISAC)
Phần này ứng dụng giải thuật DE để giải quyết bài toán thiết kế Beamforming trong hệ thống Cell-free ISAC (Integrated Sensing and Communication). Kết quả thu được từ thuật toán DE sẽ được đem so sánh với thuật toán tối ưu lồi đã được sử dụng trong bài toán gốc.
- **`Cell-free-ISAC-beamforming-main/`**: Mã nguồn MATLAB/Python mô phỏng bài toán gốc.
- **`cvx-w64/`**: Bộ giải bài toán tối ưu lồi (dùng để đối chiếu).
- **Hình ảnh kết quả (`Figure_*.png`)**: Đồ thị biểu diễn kết quả mô phỏng và so sánh hiệu năng.

## Hướng dẫn chạy
- Đối với `DE_Variants`, bạn có thể chạy file `compare_de.py` bằng Python:
  ```bash
  cd DE_Variants
  python compare_de.py
  ```
- Đối với phần `Cell_Free_ISAC`, tuỳ thuộc vào ngôn ngữ của bài toán gốc (Python hoặc MATLAB), bạn hãy chạy kịch bản mô phỏng tương ứng bên trong thư mục `Cell_Free_ISAC`.

## Tài liệu báo cáo
Các tài liệu phân tích chi tiết, đánh giá hiệu năng thuật toán và báo cáo bài tập lớn của nhóm (môn IT4593 - ĐHBK Hà Nội) có thể được tìm thấy trong thư mục `Docs/`:
- **`Nhóm 38 - Báo cáo bài tập lớn.pdf`**: Báo cáo chi tiết về cơ sở lý thuyết, mô hình hệ thống Cell-free ISAC và kết quả mô phỏng thuật toán DE so với phương pháp gốc.
- **`Nhóm 38 - Slide báo cáo BTL.pdf`**: Slide thuyết trình tóm tắt các nội dung chính của bài tập lớn.
