% Calculate inertial parameters regressor of joint inertia matrix for
% S4RPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RPRR9_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_inertiaJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t15 = cos(qJ(4));
t37 = 0.2e1 * t15;
t36 = 2 * qJ(2);
t16 = cos(qJ(3));
t35 = t16 * pkin(3);
t10 = t15 ^ 2;
t13 = sin(qJ(4));
t8 = t13 ^ 2;
t34 = t10 + t8;
t11 = t16 ^ 2;
t14 = sin(qJ(3));
t9 = t14 ^ 2;
t5 = t9 + t11;
t17 = -pkin(1) - pkin(5);
t33 = t11 * t17;
t32 = t13 * t14;
t31 = t13 * t15;
t30 = t13 * t16;
t29 = t14 * t17;
t28 = t15 * t14;
t6 = t15 * t16;
t27 = t15 * t17;
t26 = t16 * t14;
t25 = t16 * t17;
t24 = -0.2e1 * t26;
t23 = t13 * t6;
t22 = t34 * t14;
t21 = -pkin(6) * t14 - t35;
t4 = t14 * pkin(3) - t16 * pkin(6) + qJ(2);
t1 = -t13 * t29 + t15 * t4;
t2 = t13 * t4 + t14 * t27;
t20 = -t1 * t13 + t2 * t15;
t18 = qJ(2) ^ 2;
t12 = t17 ^ 2;
t7 = t11 * t12;
t3 = t5 * t17;
t19 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -2 * pkin(1), t36, pkin(1) ^ 2 + t18, t11, t24, 0, t9, 0, 0, t14 * t36, t16 * t36, -0.2e1 * t3, t9 * t12 + t18 + t7, t10 * t11, -0.2e1 * t11 * t31, t26 * t37, t8 * t11, t13 * t24, t9, 0.2e1 * t1 * t14 - 0.2e1 * t13 * t33, -0.2e1 * t11 * t27 - 0.2e1 * t2 * t14, 0.2e1 * (-t1 * t15 - t13 * t2) * t16, t1 ^ 2 + t2 ^ 2 + t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, -t5, t3, 0, 0, 0, 0, 0, 0, -t5 * t13, -t5 * t15, 0, t20 * t14 + t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34 * t9 + t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, -t14, 0, t25, -t29, 0, 0, t23, (t10 - t8) * t16, t32, -t23, t28, 0, t21 * t13 + t15 * t25, -t13 * t25 + t21 * t15, t20, pkin(3) * t25 + t20 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t14, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t30, t22, pkin(6) * t22 + t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t8, 0.2e1 * t31, 0, t10, 0, 0, pkin(3) * t37, -0.2e1 * pkin(3) * t13, 0.2e1 * t34 * pkin(6), t34 * pkin(6) ^ 2 + pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, -t30, t14, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t28, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, t15, 0, -t13 * pkin(6), -t15 * pkin(6), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t19;
