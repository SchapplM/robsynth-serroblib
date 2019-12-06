% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRRRR6
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
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRR6_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t50 = sin(qJ(3));
t72 = t50 * pkin(2);
t40 = pkin(7) + t72;
t49 = sin(qJ(4));
t46 = t49 ^ 2;
t53 = cos(qJ(4));
t47 = t53 ^ 2;
t60 = t46 + t47;
t79 = t60 * t40;
t51 = sin(qJ(2));
t54 = cos(qJ(3));
t55 = cos(qJ(2));
t27 = t50 * t51 - t54 * t55;
t22 = t27 ^ 2;
t42 = -pkin(4) * t53 - pkin(3);
t69 = t54 * pkin(2);
t31 = t42 - t69;
t78 = 0.2e1 * t31;
t77 = 0.2e1 * t42;
t76 = 0.2e1 * t53;
t48 = sin(qJ(5));
t52 = cos(qJ(5));
t25 = t48 * t49 - t52 * t53;
t29 = t48 * t53 + t49 * t52;
t19 = (-pkin(8) - t40) * t49;
t45 = t53 * pkin(8);
t64 = t53 * t40;
t20 = t45 + t64;
t8 = t19 * t52 - t20 * t48;
t9 = t19 * t48 + t20 * t52;
t75 = -t9 * t25 - t8 * t29;
t32 = (-pkin(8) - pkin(7)) * t49;
t70 = t53 * pkin(7);
t33 = t45 + t70;
t15 = t32 * t52 - t33 * t48;
t16 = t32 * t48 + t33 * t52;
t74 = -t15 * t29 - t16 * t25;
t73 = t48 * pkin(4);
t71 = t52 * pkin(4);
t41 = -pkin(3) - t69;
t68 = pkin(3) - t41;
t67 = t27 * t53;
t30 = t50 * t55 + t51 * t54;
t66 = t49 * t30;
t65 = t53 * t30;
t63 = t31 + t42;
t61 = pkin(7) * t60;
t11 = t60 * t30;
t36 = t49 * t76;
t24 = t30 ^ 2;
t23 = t29 ^ 2;
t21 = t25 ^ 2;
t18 = t27 * t49;
t14 = t27 * t29;
t13 = t27 * t25;
t12 = -0.2e1 * t29 * t25;
t10 = (-t25 * t48 - t29 * t52) * pkin(4);
t7 = -t48 * t66 + t52 * t65;
t6 = t29 * t30;
t1 = -t25 * t7 + t29 * t6;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51 ^ 2 + t55 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24 + t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24 * t60 + t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 ^ 2 + t7 ^ 2 + t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t51, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t30, 0, (-t27 * t54 + t30 * t50) * pkin(2), 0, 0, 0, 0, 0, 0, -t67, t18, t11, t27 * t41 + t30 * t79, 0, 0, 0, 0, 0, 0, t13, t14, t1, t27 * t31 - t6 * t8 + t7 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t69, -0.2e1 * t72, 0, (t50 ^ 2 + t54 ^ 2) * pkin(2) ^ 2, t46, t36, 0, t47, 0, 0, -0.2e1 * t41 * t53, 0.2e1 * t41 * t49, 0.2e1 * t79, t40 ^ 2 * t60 + t41 ^ 2, t23, t12, 0, t21, 0, 0, t25 * t78, t29 * t78, 0.2e1 * t75, t31 ^ 2 + t8 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t30, 0, 0, 0, 0, 0, 0, 0, 0, -t67, t18, t11, -t27 * pkin(3) + pkin(7) * t11, 0, 0, 0, 0, 0, 0, t13, t14, t1, -t15 * t6 + t16 * t7 + t27 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t69, -t72, 0, 0, t46, t36, 0, t47, 0, 0, t68 * t53, -t68 * t49, t61 + t79, -t41 * pkin(3) + pkin(7) * t79, t23, t12, 0, t21, 0, 0, t63 * t25, t63 * t29, t74 + t75, t15 * t8 + t16 * t9 + t31 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t46, t36, 0, t47, 0, 0, pkin(3) * t76, -0.2e1 * pkin(3) * t49, 0.2e1 * t61, pkin(7) ^ 2 * t60 + pkin(3) ^ 2, t23, t12, 0, t21, 0, 0, t25 * t77, t29 * t77, 0.2e1 * t74, t15 ^ 2 + t16 ^ 2 + t42 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, -t65, 0, 0, 0, 0, 0, 0, 0, 0, -t6, -t7, 0, (t48 * t7 - t52 * t6) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, t53, 0, -t49 * t40, -t64, 0, 0, 0, 0, t29, 0, -t25, 0, t8, -t9, t10, (t48 * t9 + t52 * t8) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, t53, 0, -t49 * pkin(7), -t70, 0, 0, 0, 0, t29, 0, -t25, 0, t15, -t16, t10, (t15 * t52 + t16 * t48) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t71, -0.2e1 * t73, 0, (t48 ^ 2 + t52 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, -t7, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, -t25, 0, t8, -t9, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, -t25, 0, t15, -t16, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t71, -t73, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t2;
