% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRPRR3_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t26 = sin(pkin(9));
t27 = cos(pkin(9));
t30 = sin(qJ(2));
t33 = cos(qJ(2));
t10 = t26 * t30 - t27 * t33;
t8 = t10 ^ 2;
t43 = t27 * pkin(2);
t22 = -pkin(3) - t43;
t32 = cos(qJ(4));
t18 = -t32 * pkin(4) + t22;
t46 = 0.2e1 * t18;
t29 = sin(qJ(4));
t45 = 0.2e1 * t29;
t44 = t26 * pkin(2);
t28 = sin(qJ(5));
t42 = t28 * pkin(4);
t31 = cos(qJ(5));
t41 = t31 * pkin(4);
t21 = pkin(6) + t44;
t40 = pkin(7) + t21;
t12 = t26 * t33 + t27 * t30;
t39 = t29 * t12;
t38 = t32 * t12;
t24 = t29 ^ 2;
t25 = t32 ^ 2;
t37 = t24 + t25;
t36 = t37 * t21;
t17 = t28 * t32 + t31 * t29;
t15 = t28 * t29 - t31 * t32;
t14 = t17 ^ 2;
t13 = t15 ^ 2;
t9 = t12 ^ 2;
t7 = t40 * t32;
t6 = t40 * t29;
t4 = -t28 * t6 + t31 * t7;
t3 = -t28 * t7 - t31 * t6;
t2 = -t28 * t39 + t31 * t38;
t1 = t17 * t12;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30 ^ 2 + t33 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9 + t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t9 + t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 ^ 2 + t2 ^ 2 + t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t30, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t12, 0, (-t10 * t27 + t12 * t26) * pkin(2), 0, 0, 0, 0, 0, 0, -t10 * t32, t10 * t29, t37 * t12, t10 * t22 + t12 * t36, 0, 0, 0, 0, 0, 0, t10 * t15, t10 * t17, t1 * t17 - t2 * t15, -t1 * t3 + t10 * t18 + t2 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t43, -0.2e1 * t44, 0, (t26 ^ 2 + t27 ^ 2) * pkin(2) ^ 2, t24, t32 * t45, 0, t25, 0, 0, -0.2e1 * t22 * t32, t22 * t45, 0.2e1 * t36, t37 * t21 ^ 2 + t22 ^ 2, t14, -0.2e1 * t17 * t15, 0, t13, 0, 0, t15 * t46, t17 * t46, -0.2e1 * t4 * t15 - 0.2e1 * t3 * t17, t18 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * t15 + t2 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3 * t15 + t4 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 + t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, -t38, 0, 0, 0, 0, 0, 0, 0, 0, -t1, -t2, 0, (-t1 * t31 + t2 * t28) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, t32, 0, -t29 * t21, -t32 * t21, 0, 0, 0, 0, t17, 0, -t15, 0, t3, -t4, (-t15 * t28 - t17 * t31) * pkin(4), (t28 * t4 + t3 * t31) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t29, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t17, 0, (-t15 * t31 + t17 * t28) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t41, -0.2e1 * t42, 0, (t28 ^ 2 + t31 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, -t15, 0, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t17, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t41, -t42, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t5;
