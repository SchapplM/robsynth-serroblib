% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRPR1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t34 = sin(pkin(9));
t32 = t34 ^ 2;
t35 = cos(pkin(9));
t33 = t35 ^ 2;
t41 = t32 + t33;
t42 = t41 * qJ(4);
t27 = -t35 * pkin(4) - pkin(3);
t38 = cos(qJ(3));
t47 = t38 * pkin(2);
t18 = t27 - t47;
t53 = 0.2e1 * t18;
t52 = 0.2e1 * t27;
t51 = 0.2e1 * t35;
t36 = sin(qJ(5));
t45 = cos(qJ(5));
t15 = t36 * t34 - t45 * t35;
t17 = t45 * t34 + t36 * t35;
t37 = sin(qJ(3));
t48 = t37 * pkin(2);
t26 = qJ(4) + t48;
t11 = (-pkin(7) - t26) * t34;
t31 = t35 * pkin(7);
t12 = t35 * t26 + t31;
t6 = t45 * t11 - t36 * t12;
t7 = t36 * t11 + t45 * t12;
t50 = -t7 * t15 - t6 * t17;
t19 = (-pkin(7) - qJ(4)) * t34;
t20 = t35 * qJ(4) + t31;
t10 = t36 * t19 + t45 * t20;
t9 = t45 * t19 - t36 * t20;
t49 = -t10 * t15 - t9 * t17;
t28 = -pkin(3) - t47;
t46 = pkin(3) - t28;
t44 = t18 + t27;
t43 = t41 * t26;
t23 = t34 * t51;
t14 = t17 ^ 2;
t13 = t15 ^ 2;
t8 = -0.2e1 * t17 * t15;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 + t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15 * t6 + t17 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t47, -0.2e1 * t48, 0, (t37 ^ 2 + t38 ^ 2) * pkin(2) ^ 2, t32, t23, 0, t33, 0, 0, -0.2e1 * t28 * t35, 0.2e1 * t28 * t34, 0.2e1 * t43, t41 * t26 ^ 2 + t28 ^ 2, t14, t8, 0, t13, 0, 0, t15 * t53, t17 * t53, 0.2e1 * t50, t18 ^ 2 + t6 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 * t10 - t15 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t47, -t48, 0, 0, t32, t23, 0, t33, 0, 0, t46 * t35, -t46 * t34, t42 + t43, -t28 * pkin(3) + t26 * t42, t14, t8, 0, t13, 0, 0, t44 * t15, t44 * t17, t49 + t50, t7 * t10 + t18 * t27 + t6 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t32, t23, 0, t33, 0, 0, pkin(3) * t51, -0.2e1 * pkin(3) * t34, 0.2e1 * t42, t41 * qJ(4) ^ 2 + pkin(3) ^ 2, t14, t8, 0, t13, 0, 0, t15 * t52, t17 * t52, 0.2e1 * t49, t10 ^ 2 + t27 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, t34, 0, t28, 0, 0, 0, 0, 0, 0, t15, t17, 0, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, t34, 0, -pkin(3), 0, 0, 0, 0, 0, 0, t15, t17, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t17, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, -t15, 0, t6, -t7, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, -t15, 0, t9, -t10, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t1;
