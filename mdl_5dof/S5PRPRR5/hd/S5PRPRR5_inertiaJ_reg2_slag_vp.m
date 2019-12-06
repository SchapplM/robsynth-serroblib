% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRPRR5
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
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRPRR5_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t28 = sin(pkin(9));
t29 = cos(pkin(9));
t31 = sin(qJ(4));
t34 = cos(qJ(4));
t17 = t31 * t28 - t34 * t29;
t23 = -t29 * pkin(3) - pkin(2);
t12 = t17 * pkin(4) + t23;
t45 = 0.2e1 * t12;
t44 = 0.2e1 * t23;
t43 = 0.2e1 * t29;
t30 = sin(qJ(5));
t42 = t30 * pkin(4);
t33 = cos(qJ(5));
t41 = t33 * pkin(4);
t40 = pkin(6) + qJ(3);
t24 = t28 ^ 2;
t25 = t29 ^ 2;
t39 = t24 + t25;
t20 = t40 * t28;
t21 = t40 * t29;
t10 = -t34 * t20 - t31 * t21;
t38 = t39 * qJ(3);
t11 = -t31 * t20 + t34 * t21;
t19 = t34 * t28 + t31 * t29;
t35 = cos(qJ(2));
t32 = sin(qJ(2));
t27 = t35 ^ 2;
t26 = t32 ^ 2;
t14 = t17 * t32;
t13 = t19 * t32;
t9 = -t30 * t17 + t33 * t19;
t7 = t33 * t17 + t30 * t19;
t6 = -t17 * pkin(7) + t11;
t5 = -t19 * pkin(7) + t10;
t4 = -t30 * t13 - t33 * t14;
t3 = -t33 * t13 + t30 * t14;
t2 = t30 * t5 + t33 * t6;
t1 = -t30 * t6 + t33 * t5;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26 + t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39 * t26 + t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 ^ 2 + t14 ^ 2 + t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 ^ 2 + t4 ^ 2 + t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t32, 0, 0, 0, 0, 0, 0, 0, 0, t35 * t29, -t35 * t28, t39 * t32, t35 * pkin(2) + t32 * t38, 0, 0, 0, 0, 0, 0, -t35 * t17, -t35 * t19, t13 * t19 + t14 * t17, -t13 * t10 - t14 * t11 - t35 * t23, 0, 0, 0, 0, 0, 0, -t35 * t7, -t35 * t9, -t3 * t9 - t4 * t7, t3 * t1 - t35 * t12 + t4 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t24, t28 * t43, 0, t25, 0, 0, pkin(2) * t43, -0.2e1 * pkin(2) * t28, 0.2e1 * t38, t39 * qJ(3) ^ 2 + pkin(2) ^ 2, t19 ^ 2, -0.2e1 * t19 * t17, 0, t17 ^ 2, 0, 0, t17 * t44, t19 * t44, -0.2e1 * t10 * t19 - 0.2e1 * t11 * t17, t10 ^ 2 + t11 ^ 2 + t23 ^ 2, t9 ^ 2, -0.2e1 * t9 * t7, 0, t7 ^ 2, 0, 0, t7 * t45, t9 * t45, -0.2e1 * t1 * t9 - 0.2e1 * t2 * t7, t1 ^ 2 + t12 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, t28, 0, -pkin(2), 0, 0, 0, 0, 0, 0, t17, t19, 0, t23, 0, 0, 0, 0, 0, 0, t7, t9, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, t14, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4, 0, (t3 * t33 + t30 * t4) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, -t17, 0, t10, -t11, 0, 0, 0, 0, t9, 0, -t7, 0, t1, -t2, (-t30 * t7 - t33 * t9) * pkin(4), (t1 * t33 + t2 * t30) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t41, -0.2e1 * t42, 0, (t30 ^ 2 + t33 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, -t7, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t41, -t42, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t8;
