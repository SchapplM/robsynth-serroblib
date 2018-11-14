% Calculate inertial parameters regressor of joint inertia matrix for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:46
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function MM_reg = S4RPPP1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_inertiaJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t22 = sin(pkin(4));
t34 = 0.2e1 * t22;
t19 = t22 ^ 2;
t23 = cos(pkin(6));
t33 = t19 * t23;
t21 = sin(pkin(6));
t15 = t21 * t22;
t24 = cos(pkin(4));
t32 = t21 * t24;
t18 = t22 * t23;
t31 = t23 * t24;
t30 = qJ(2) * t22;
t8 = pkin(1) * t32 + t23 * t30;
t29 = t21 * t33;
t28 = t22 * t32;
t27 = t22 * t31;
t26 = -pkin(1) * t23 - pkin(2);
t25 = -qJ(3) * t21 - pkin(1);
t4 = -t24 * qJ(3) - t8;
t20 = t24 ^ 2;
t17 = t19 * t23 ^ 2;
t16 = t19 * t21 ^ 2;
t12 = t21 * t30;
t11 = -0.2e1 * t27;
t10 = 0.2e1 * t28;
t9 = 0.2e1 * t29;
t7 = pkin(1) * t31 - t12;
t6 = (-pkin(2) * t23 + t25) * t22;
t5 = t26 * t24 + t12;
t3 = ((-pkin(2) - qJ(4)) * t23 + t25) * t22;
t2 = pkin(3) * t18 - t4;
t1 = pkin(3) * t15 + t12 + (-qJ(4) + t26) * t24;
t13 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t16, t9, t10, t17, 0.2e1 * t27, t20, 0.2e1 * pkin(1) * t33 + 0.2e1 * t7 * t24, -0.2e1 * t19 * pkin(1) * t21 - 0.2e1 * t8 * t24 (-t21 * t7 + t23 * t8) * t34, t19 * pkin(1) ^ 2 + t7 ^ 2 + t8 ^ 2, t20, -0.2e1 * t28, t11, t16, t9, t17 (t21 * t5 - t23 * t4) * t34, 0.2e1 * t6 * t18 + 0.2e1 * t5 * t24, -0.2e1 * t6 * t15 - 0.2e1 * t4 * t24, t4 ^ 2 + t5 ^ 2 + t6 ^ 2, t20, t11, t10, t17, -0.2e1 * t29, t16 (t1 * t21 + t2 * t23) * t34, -0.2e1 * t3 * t15 + 0.2e1 * t2 * t24, -0.2e1 * t1 * t24 - 0.2e1 * t3 * t18, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t15, 0, -t22 * pkin(1), 0, 0, 0, 0, 0, 0, 0, t18, -t15, t6, 0, 0, 0, 0, 0, 0, 0, -t15, -t18, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, t24, 0, t5, 0, 0, 0, 0, 0, 0, t15, 0, -t24, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t24, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t13;
