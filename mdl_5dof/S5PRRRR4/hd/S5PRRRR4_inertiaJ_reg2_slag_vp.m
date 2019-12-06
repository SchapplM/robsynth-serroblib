% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRR4_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t37 = sin(qJ(3));
t53 = t37 * pkin(2);
t27 = pkin(7) + t53;
t36 = sin(qJ(4));
t33 = t36 ^ 2;
t39 = cos(qJ(4));
t34 = t39 ^ 2;
t44 = t33 + t34;
t46 = t44 * t27;
t29 = -pkin(4) * t39 - pkin(3);
t40 = cos(qJ(3));
t50 = t40 * pkin(2);
t19 = t29 - t50;
t59 = 0.2e1 * t19;
t58 = 0.2e1 * t29;
t57 = 0.2e1 * t39;
t35 = sin(qJ(5));
t38 = cos(qJ(5));
t16 = t35 * t36 - t38 * t39;
t18 = t35 * t39 + t36 * t38;
t12 = (-pkin(8) - t27) * t36;
t32 = t39 * pkin(8);
t48 = t39 * t27;
t13 = t32 + t48;
t6 = t12 * t38 - t13 * t35;
t7 = t12 * t35 + t13 * t38;
t56 = -t7 * t16 - t6 * t18;
t20 = (-pkin(8) - pkin(7)) * t36;
t51 = t39 * pkin(7);
t21 = t32 + t51;
t10 = t20 * t38 - t21 * t35;
t11 = t20 * t35 + t21 * t38;
t55 = -t10 * t18 - t11 * t16;
t54 = t35 * pkin(4);
t52 = t38 * pkin(4);
t28 = -pkin(3) - t50;
t49 = pkin(3) - t28;
t47 = t19 + t29;
t45 = t44 * pkin(7);
t24 = t36 * t57;
t15 = t18 ^ 2;
t14 = t16 ^ 2;
t9 = -0.2e1 * t18 * t16;
t8 = (-t16 * t35 - t18 * t38) * pkin(4);
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 + t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16 * t6 + t18 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t50, -0.2e1 * t53, 0, (t37 ^ 2 + t40 ^ 2) * pkin(2) ^ 2, t33, t24, 0, t34, 0, 0, -0.2e1 * t28 * t39, 0.2e1 * t28 * t36, 0.2e1 * t46, t44 * t27 ^ 2 + t28 ^ 2, t15, t9, 0, t14, 0, 0, t16 * t59, t18 * t59, 0.2e1 * t56, t19 ^ 2 + t6 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10 * t16 + t11 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t50, -t53, 0, 0, t33, t24, 0, t34, 0, 0, t49 * t39, -t49 * t36, t45 + t46, -t28 * pkin(3) + pkin(7) * t46, t15, t9, 0, t14, 0, 0, t47 * t16, t47 * t18, t55 + t56, t10 * t6 + t11 * t7 + t19 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t33, t24, 0, t34, 0, 0, pkin(3) * t57, -0.2e1 * pkin(3) * t36, 0.2e1 * t45, t44 * pkin(7) ^ 2 + pkin(3) ^ 2, t15, t9, 0, t14, 0, 0, t16 * t58, t18 * t58, 0.2e1 * t55, t10 ^ 2 + t11 ^ 2 + t29 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t36, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t18, 0, (-t16 * t38 + t18 * t35) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, t39, 0, -t36 * t27, -t48, 0, 0, 0, 0, t18, 0, -t16, 0, t6, -t7, t8, (t35 * t7 + t38 * t6) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, t39, 0, -t36 * pkin(7), -t51, 0, 0, 0, 0, t18, 0, -t16, 0, t10, -t11, t8, (t10 * t38 + t11 * t35) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t52, -0.2e1 * t54, 0, (t35 ^ 2 + t38 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t18, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, -t16, 0, t6, -t7, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, -t16, 0, t10, -t11, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t52, -t54, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t1;
