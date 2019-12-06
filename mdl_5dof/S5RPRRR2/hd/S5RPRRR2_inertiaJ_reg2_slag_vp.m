% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRR2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t38 = sin(pkin(9));
t39 = cos(pkin(9));
t42 = sin(qJ(3));
t51 = cos(qJ(3));
t25 = t42 * t38 - t51 * t39;
t27 = t51 * t38 + t42 * t39;
t41 = sin(qJ(4));
t44 = cos(qJ(4));
t15 = t44 * t25 + t41 * t27;
t32 = -t39 * pkin(2) - pkin(1);
t20 = t25 * pkin(3) + t32;
t10 = t15 * pkin(4) + t20;
t57 = 0.2e1 * t10;
t56 = 0.2e1 * t20;
t55 = 0.2e1 * t32;
t54 = 0.2e1 * t39;
t40 = sin(qJ(5));
t53 = t40 * pkin(4);
t52 = t41 * pkin(3);
t50 = pkin(6) + qJ(2);
t36 = t38 ^ 2;
t37 = t39 ^ 2;
t49 = t36 + t37;
t43 = cos(qJ(5));
t48 = t43 * t52;
t28 = t50 * t38;
t29 = t50 * t39;
t18 = -t51 * t28 - t42 * t29;
t12 = -t27 * pkin(7) + t18;
t19 = -t42 * t28 + t51 * t29;
t13 = -t25 * pkin(7) + t19;
t5 = t44 * t12 - t41 * t13;
t35 = t44 * pkin(3);
t33 = t35 + pkin(4);
t21 = t43 * t33 - t40 * t52;
t6 = t41 * t12 + t44 * t13;
t34 = t43 * pkin(4);
t22 = t40 * t33 + t48;
t17 = -t41 * t25 + t44 * t27;
t9 = -t40 * t15 + t43 * t17;
t7 = t43 * t15 + t40 * t17;
t4 = -t15 * pkin(8) + t6;
t3 = -t17 * pkin(8) + t5;
t2 = t40 * t3 + t43 * t4;
t1 = t43 * t3 - t40 * t4;
t8 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t36, t38 * t54, 0, t37, 0, 0, pkin(1) * t54, -0.2e1 * pkin(1) * t38, 0.2e1 * t49 * qJ(2), t49 * qJ(2) ^ 2 + pkin(1) ^ 2, t27 ^ 2, -0.2e1 * t27 * t25, 0, t25 ^ 2, 0, 0, t25 * t55, t27 * t55, -0.2e1 * t18 * t27 - 0.2e1 * t19 * t25, t18 ^ 2 + t19 ^ 2 + t32 ^ 2, t17 ^ 2, -0.2e1 * t17 * t15, 0, t15 ^ 2, 0, 0, t15 * t56, t17 * t56, -0.2e1 * t6 * t15 - 0.2e1 * t5 * t17, t20 ^ 2 + t5 ^ 2 + t6 ^ 2, t9 ^ 2, -0.2e1 * t9 * t7, 0, t7 ^ 2, 0, 0, t7 * t57, t9 * t57, -0.2e1 * t1 * t9 - 0.2e1 * t2 * t7, t1 ^ 2 + t10 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t38, 0, -pkin(1), 0, 0, 0, 0, 0, 0, t25, t27, 0, t32, 0, 0, 0, 0, 0, 0, t15, t17, 0, t20, 0, 0, 0, 0, 0, 0, t7, t9, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, -t25, 0, t18, -t19, 0, 0, 0, 0, t17, 0, -t15, 0, t5, -t6, (-t15 * t41 - t17 * t44) * pkin(3), (t41 * t6 + t44 * t5) * pkin(3), 0, 0, t9, 0, -t7, 0, t1, -t2, -t21 * t9 - t22 * t7, t1 * t21 + t2 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t35, -0.2e1 * t52, 0, (t41 ^ 2 + t44 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t21, -0.2e1 * t22, 0, t21 ^ 2 + t22 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, -t15, 0, t5, -t6, 0, 0, 0, 0, t9, 0, -t7, 0, t1, -t2, (-t40 * t7 - t43 * t9) * pkin(4), (t1 * t43 + t2 * t40) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t35, -t52, 0, 0, 0, 0, 0, 0, 0, 1, t21 + t34, -t48 + (-pkin(4) - t33) * t40, 0, (t21 * t43 + t22 * t40) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t34, -0.2e1 * t53, 0, (t40 ^ 2 + t43 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, -t7, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t21, -t22, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t34, -t53, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t8;
