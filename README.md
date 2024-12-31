# Adaptive and Wiener Filtering for ECG Signal Processing

## Project Overview
This project investigates the application of Wiener and adaptive filters for ECG signal processing. 
The focus is on noise reduction and the comparison of filter designs, including the Least Mean Squares (LMS) and Recursive Least Squares (RLS) methods, 
for handling stationary and non-stationary noise in bio-signals.

---

## Contents

### **1. Wiener Filtering**
- **Time-Domain Implementation**:
  - Determined the optimum filter order by minimizing the Mean Square Error (MSE).
  - Observed successful attenuation of 50 Hz noise and high-frequency components beyond 150 Hz.
  - ![image](https://github.com/user-attachments/assets/f4421a27-1f05-49e0-87c5-7a8fb43c97dd)
  - ![image](https://github.com/user-attachments/assets/62b21250-e0c8-4266-92f7-997eb1c2530d)

- **Frequency Domain Implementation**:
  - Implemented Wiener filtering in the frequency domain.
    ![image](https://github.com/user-attachments/assets/7b76cbaf-26ea-45c3-89e3-af5df9de297c)

Used an ideal ECG template and a manually modeled ECG signal as the desired signal and analysed the difference; ideal ECG slightly outperforming the modeled ECG due to its closer representation of real-world characteristics.

- **Effect of Non-Stationary Noise**:
  - Wiener filtering is less effective for non-stationary noise as it assumes stationary noise and signal characteristics.
    ![image](https://github.com/user-attachments/assets/b2f3fbec-dd0f-4645-8884-520497d209ae)

---

### **2. Adaptive Filtering**
- **Methodology**:
  - Primary signal: ECG signal with added noise.
  - Reference signal: A signal closely related to the noise, uncorrelated with the desired signal.
  - Adaptive filters dynamically update coefficients to minimize error between the filtered and desired signals.
    ![image](https://github.com/user-attachments/assets/1dc71b5a-2b5c-415e-baa3-c9d36c717b2f)

- **LMS Method**:
  - Utilized gradient descent to iteratively update filter weights.
  - Performance evaluated for a sawtooth ECG signal with added non-stationary noise.
    ![image](https://github.com/user-attachments/assets/3f8411b0-771d-4782-abd3-76ba5e519a9c)
    ** Note:- Here error signal is the filtered signal.


- **RLS Method**:
  - Introduced a forgetting factor to account for previous errors in weight updates.
    ![image](https://github.com/user-attachments/assets/bc3722e0-673d-4618-9166-7b3581d3ad9b)
     ** Note:- Here error signal is the filtered signal.
    

- **Comparison**:
  - Both LMS and RLS filters effectively attenuated the noise, outperforming Wiener filtering for non-stationary noise scenarios. RLS provides slightly better adaptation due to its incorporation of historical error data.

---

## Acknowledgments
This project was completed as part of the BM4152 Bio-Signal Processing course.

