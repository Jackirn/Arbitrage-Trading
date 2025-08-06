# Energy Trading – Mean-Reversion Strategy on Heating Oil & Gasoil Futures

This repository presents a comprehensive analysis of **mean-reverting statistical arbitrage strategies** applied to a futures pair:  
**Heating Oil (HO)** and **Low-Sulfur Gasoil (LGO)**, using high-frequency intraday data.

The study includes:
- Data loading, filtering, and preprocessing
- Ornstein-Uhlenbeck (OU) process calibration using bootstrap methods
- Optimal trading band selection under transaction costs and stop-loss constraints
- In-sample (IS) and out-of-sample (OS) performance assessment

Two parallel implementations are provided:
- `Energy-Trading-Python/`: Complete Python pipeline
- `Energy-Trading-MATLAB/`: MATLAB modules and scripts

---

## 📁 Repository Structure

### Root Files
- `README.md`: This file
- `HO-LGO.xlsm`: Raw futures data in Excel with macros
- `requirements.txt`: Python dependencies
- `Final_Project_Group_8B_Report.pdf`: Final report (PDF-style)
- `.gitignore`: Git tracking exclusions
- `venv/`: Python virtual environment (not tracked)

---

### 🔧 Python Directory — `Energy-Trading-Python/`

Python implementation of the project with utilities and crypto extensions:

| File/Folder                          | Description |
|-------------------------------------|-------------|
| `Final_Project_Group_*.ipynb`       | 📓 Jupyter notebooks with full pipeline |
| `Final_Project_Crypto_*.ipynb`      | Crypto adaptation of the strategy |
| `crypto_dataset/`                   | Additional datasets for crypto testing |
| `utilities/`                        | Modular utility functions (data, plots, stats, etc.) |
| `path/`                             | File and directory path utilities |
| `PDF/`                              | Generated PDF outputs |

---

### 🧮 MATLAB Directory — `Energy-Trading-MATLAB/`

Modular MATLAB codebase structured into functional blocks:

| Folder                        | Description |
|-------------------------------|-------------|
| `Cost & Analysis Utilities/`  | Transaction cost handling and performance metrics |
| `Data/`                       | Data files (converted/cleaned) |
| `Data Processing/`            | Trimming and filtering scripts |
| `OU model Estimation .../`    | OU model calibration and bootstrap |
| `Trading Strategy & Bands/`   | Trading rules and band computation |
| `Prints and Plots/`           | Visual results |
| `Final_Project_Group_*.mlx`   | Main live script notebook (MATLAB) |

---

## ⚙️ How to Run (Python)

1. **Install dependencies** (preferably in a virtual environment):

```bash
pip install -r requirements.txt
```

2. **Run the main Jupyter Notebooks**:

```bash
jupyter notebook Energy-Trading-Python/Final_Project_Group_*.ipynb
```

---

## 📚 Academic Context

This project was developed as part of a **university coursework** on quantitative trading and statistical arbitrage.  
The methodology is based on the paper by **Baviera (2019)**, available in the `PDF/` folder.

---

## ✅ Evaluation

This project received the **highest grade** in the course, validating both the **correctness of the implementation** and the **depth of the analysis**.  
All results and methods were thoroughly reviewed by the academic supervisors.

---

## 📄 License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.

---

## 🔒 Note on Repository Privacy

Please note that the **original repository**, which contains the **complete commit history and development process**, is private.  
This is due to the presence of **academic material and internal coursework content** that cannot be publicly shared.

This public version has been curated for **portfolio and educational purposes**, and faithfully reflects the structure, results, and methodology of the original work.