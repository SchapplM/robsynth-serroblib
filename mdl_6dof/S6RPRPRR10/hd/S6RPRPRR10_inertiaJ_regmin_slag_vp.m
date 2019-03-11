% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x31]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPRR10_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t51 = sin(pkin(10));
t52 = cos(pkin(10));
t54 = sin(qJ(5));
t57 = cos(qJ(5));
t37 = t54 * t51 - t57 * t52;
t46 = -t52 * pkin(4) - pkin(3);
t27 = t37 * pkin(5) + t46;
t84 = 0.2e1 * t27;
t83 = 0.2e1 * t46;
t55 = sin(qJ(3));
t82 = -0.2e1 * t55;
t81 = 0.2e1 * t55;
t80 = 2 * qJ(2);
t53 = sin(qJ(6));
t79 = t53 * pkin(5);
t78 = t55 * pkin(5);
t56 = cos(qJ(6));
t77 = t56 * pkin(5);
t38 = t57 * t51 + t54 * t52;
t58 = cos(qJ(3));
t29 = t58 * t38;
t39 = t55 * pkin(3) - t58 * qJ(4) + qJ(2);
t34 = t52 * t39;
t44 = t52 * t58;
t59 = -pkin(1) - pkin(7);
t70 = t51 * t59;
t17 = -pkin(8) * t44 + t34 + (pkin(4) - t70) * t55;
t69 = t55 * t59;
t25 = t51 * t39 + t52 * t69;
t71 = t51 * t58;
t21 = -pkin(8) * t71 + t25;
t9 = t54 * t17 + t57 * t21;
t5 = -t29 * pkin(9) + t9;
t76 = t56 * t5;
t75 = t58 * pkin(3);
t74 = t37 * t55;
t73 = t38 * t55;
t50 = t58 ^ 2;
t72 = t50 * t59;
t68 = t58 * t37;
t67 = t58 * t59;
t66 = pkin(8) + qJ(4);
t65 = t51 ^ 2 + t52 ^ 2;
t49 = t55 ^ 2;
t64 = -t49 - t50;
t8 = t57 * t17 - t54 * t21;
t4 = pkin(9) * t68 + t78 + t8;
t1 = t56 * t4 - t53 * t5;
t40 = t66 * t51;
t41 = t66 * t52;
t22 = -t57 * t40 - t54 * t41;
t35 = pkin(4) * t71 - t67;
t63 = t65 * qJ(4);
t62 = -qJ(4) * t55 - t75;
t24 = -t51 * t69 + t34;
t61 = -t24 * t51 + t25 * t52;
t23 = -t54 * t40 + t57 * t41;
t20 = t29 * pkin(5) + t35;
t19 = -t53 * t37 + t56 * t38;
t18 = t56 * t37 + t53 * t38;
t15 = -t37 * pkin(9) + t23;
t14 = -t38 * pkin(9) + t22;
t13 = -t53 * t29 - t56 * t68;
t12 = -t53 * t73 - t56 * t74;
t11 = t56 * t29 - t53 * t68;
t10 = t53 * t74 - t56 * t73;
t7 = t53 * t14 + t56 * t15;
t6 = t56 * t14 - t53 * t15;
t2 = t53 * t4 + t76;
t3 = [1, 0, 0, -2 * pkin(1), t80, pkin(1) ^ 2 + qJ(2) ^ 2, t50, t58 * t82, 0, 0, 0, t55 * t80, t58 * t80, 0.2e1 * t24 * t55 - 0.2e1 * t50 * t70, -0.2e1 * t25 * t55 - 0.2e1 * t52 * t72, 0.2e1 * (-t24 * t52 - t25 * t51) * t58, t50 * t59 ^ 2 + t24 ^ 2 + t25 ^ 2, t68 ^ 2, 0.2e1 * t68 * t29, -t68 * t81, t29 * t82, t49, 0.2e1 * t35 * t29 + 0.2e1 * t8 * t55, -0.2e1 * t35 * t68 - 0.2e1 * t9 * t55, t13 ^ 2, -0.2e1 * t13 * t11, t13 * t81, t11 * t82, t49, 0.2e1 * t1 * t55 + 0.2e1 * t20 * t11, 0.2e1 * t20 * t13 - 0.2e1 * t2 * t55; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, t64 * t51, t64 * t52, 0, t61 * t55 + t72, 0, 0, 0, 0, 0, -t58 * t29 - t55 * t73, t55 * t74 + t58 * t68, 0, 0, 0, 0, 0, t10 * t55 - t58 * t11, -t12 * t55 - t58 * t13; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65 * t49 + t50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t58, -t55, 0, t67, -t69, t62 * t51 + t52 * t67, -t51 * t67 + t62 * t52, t61, pkin(3) * t67 + t61 * qJ(4), -t68 * t38, -t38 * t29 + t37 * t68, t73, -t74, 0, t22 * t55 + t46 * t29 + t35 * t37, -t23 * t55 + t35 * t38 - t46 * t68, t13 * t19, -t19 * t11 - t13 * t18, t19 * t55, -t18 * t55, 0, t27 * t11 + t20 * t18 + t6 * t55, t27 * t13 + t20 * t19 - t7 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, -t55, t44, -t71, t65 * t55, t55 * t63 + t75, 0, 0, 0, 0, 0, -t68, -t29, 0, 0, 0, 0, 0, -t58 * t18, -t58 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3) * t52, -0.2e1 * pkin(3) * t51, 0.2e1 * t63, t65 * qJ(4) ^ 2 + pkin(3) ^ 2, t38 ^ 2, -0.2e1 * t38 * t37, 0, 0, 0, t37 * t83, t38 * t83, t19 ^ 2, -0.2e1 * t19 * t18, 0, 0, 0, t18 * t84, t19 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, t44, 0, -t67, 0, 0, 0, 0, 0, t29, -t68, 0, 0, 0, 0, 0, t11, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, t51, 0, -pkin(3), 0, 0, 0, 0, 0, t37, t38, 0, 0, 0, 0, 0, t18, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, -t29, t55, t8, -t9, 0, 0, t13, -t11, t55, t55 * t77 + t1, -t76 + (-t4 - t78) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, t74, 0, 0, 0, 0, 0, t10, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t37, 0, t22, -t23, 0, 0, t19, -t18, 0, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t77, -0.2e1 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t11, t55, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, 0, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t77, -t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
