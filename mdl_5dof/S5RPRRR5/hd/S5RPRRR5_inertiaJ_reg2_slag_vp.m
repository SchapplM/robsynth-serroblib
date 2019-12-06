% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRRR5
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
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRR5_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t40 = cos(pkin(9));
t60 = t40 * pkin(1);
t31 = pkin(2) + t60;
t43 = sin(qJ(3));
t46 = cos(qJ(3));
t39 = sin(pkin(9));
t61 = t39 * pkin(1);
t21 = t43 * t31 + t46 * t61;
t19 = pkin(7) + t21;
t42 = sin(qJ(4));
t37 = t42 ^ 2;
t45 = cos(qJ(4));
t38 = t45 ^ 2;
t50 = t37 + t38;
t52 = t50 * t19;
t20 = t46 * t31 - t43 * t61;
t18 = -pkin(3) - t20;
t57 = t45 * pkin(4);
t14 = t18 - t57;
t66 = 0.2e1 * t14;
t33 = -pkin(3) - t57;
t65 = 0.2e1 * t33;
t64 = 0.2e1 * t45;
t41 = sin(qJ(5));
t44 = cos(qJ(5));
t24 = t41 * t42 - t44 * t45;
t26 = t41 * t45 + t44 * t42;
t10 = (-pkin(8) - t19) * t42;
t36 = t45 * pkin(8);
t54 = t45 * t19;
t11 = t36 + t54;
t3 = t44 * t10 - t41 * t11;
t4 = t41 * t10 + t44 * t11;
t63 = -t4 * t24 - t3 * t26;
t27 = (-pkin(8) - pkin(7)) * t42;
t56 = t45 * pkin(7);
t28 = t36 + t56;
t12 = t44 * t27 - t41 * t28;
t13 = t41 * t27 + t44 * t28;
t62 = -t12 * t26 - t13 * t24;
t59 = t41 * pkin(4);
t58 = t44 * pkin(4);
t55 = pkin(3) - t18;
t53 = t14 + t33;
t51 = t50 * pkin(7);
t30 = t42 * t64;
t23 = t26 ^ 2;
t22 = t24 ^ 2;
t9 = -0.2e1 * t26 * t24;
t8 = (-t24 * t41 - t26 * t44) * pkin(4);
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t60, -0.2e1 * t61, 0, (t39 ^ 2 + t40 ^ 2) * pkin(1) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t20, -0.2e1 * t21, 0, t20 ^ 2 + t21 ^ 2, t37, t30, 0, t38, 0, 0, -0.2e1 * t18 * t45, 0.2e1 * t18 * t42, 0.2e1 * t52, t50 * t19 ^ 2 + t18 ^ 2, t23, t9, 0, t22, 0, 0, t24 * t66, t26 * t66, 0.2e1 * t63, t14 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3 * t24 + t4 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 + t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t20, -t21, 0, 0, t37, t30, 0, t38, 0, 0, t55 * t45, -t55 * t42, t51 + t52, -t18 * pkin(3) + pkin(7) * t52, t23, t9, 0, t22, 0, 0, t53 * t24, t53 * t26, t62 + t63, t3 * t12 + t4 * t13 + t14 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24 * t12 + t26 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t37, t30, 0, t38, 0, 0, pkin(3) * t64, -0.2e1 * pkin(3) * t42, 0.2e1 * t51, t50 * pkin(7) ^ 2 + pkin(3) ^ 2, t23, t9, 0, t22, 0, 0, t24 * t65, t26 * t65, 0.2e1 * t62, t12 ^ 2 + t13 ^ 2 + t33 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, t45, 0, -t42 * t19, -t54, 0, 0, 0, 0, t26, 0, -t24, 0, t3, -t4, t8, (t3 * t44 + t4 * t41) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t42, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t26, 0, (-t24 * t44 + t26 * t41) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, t45, 0, -t42 * pkin(7), -t56, 0, 0, 0, 0, t26, 0, -t24, 0, t12, -t13, t8, (t12 * t44 + t13 * t41) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t58, -0.2e1 * t59, 0, (t41 ^ 2 + t44 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, -t24, 0, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t26, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, -t24, 0, t12, -t13, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t58, -t59, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t1;
