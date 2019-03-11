% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRPR8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t56 = sin(pkin(10));
t57 = cos(pkin(10));
t59 = sin(qJ(4));
t62 = cos(qJ(4));
t39 = t56 * t62 + t57 * t59;
t63 = cos(qJ(3));
t31 = t39 * t63;
t50 = t62 * t63;
t78 = t59 * t63;
t33 = t57 * t50 - t56 * t78;
t58 = sin(qJ(6));
t61 = cos(qJ(6));
t11 = t61 * t31 + t58 * t33;
t85 = -0.2e1 * t11;
t38 = -t56 * t59 + t57 * t62;
t51 = -t62 * pkin(4) - pkin(3);
t29 = -t38 * pkin(5) + t51;
t84 = 0.2e1 * t29;
t83 = 0.2e1 * t62;
t82 = 2 * qJ(2);
t81 = pkin(4) * t56;
t60 = sin(qJ(3));
t80 = t59 * t60;
t79 = t59 * t62;
t64 = -pkin(1) - pkin(7);
t77 = t59 * t64;
t76 = t60 * t64;
t75 = t62 * t60;
t74 = t62 * t64;
t73 = t63 * t60;
t72 = t63 * t64;
t71 = -qJ(5) - pkin(8);
t44 = t60 * pkin(3) - t63 * pkin(8) + qJ(2);
t37 = t62 * t44;
t69 = qJ(5) * t63;
t18 = -t62 * t69 + t37 + (pkin(4) - t77) * t60;
t67 = t60 * t74;
t23 = t67 + (t44 - t69) * t59;
t9 = t56 * t18 + t57 * t23;
t45 = t71 * t59;
t46 = t71 * t62;
t25 = t56 * t45 - t57 * t46;
t53 = t60 ^ 2;
t55 = t63 ^ 2;
t70 = -t53 - t55;
t68 = -0.2e1 * t73;
t8 = t57 * t18 - t56 * t23;
t4 = t60 * pkin(5) - t33 * pkin(9) + t8;
t5 = -t31 * pkin(9) + t9;
t1 = t61 * t4 - t58 * t5;
t24 = t57 * t45 + t56 * t46;
t40 = pkin(4) * t78 - t72;
t66 = -pkin(3) * t63 - pkin(8) * t60;
t2 = t58 * t4 + t61 * t5;
t54 = t62 ^ 2;
t52 = t59 ^ 2;
t49 = t57 * pkin(4) + pkin(5);
t35 = t58 * t49 + t61 * t81;
t34 = t61 * t49 - t58 * t81;
t32 = t38 * t60;
t30 = t39 * t60;
t27 = t59 * t44 + t67;
t26 = -t59 * t76 + t37;
t22 = t31 * pkin(5) + t40;
t20 = t58 * t38 + t61 * t39;
t19 = -t61 * t38 + t58 * t39;
t15 = t38 * pkin(9) + t25;
t14 = -t39 * pkin(9) + t24;
t13 = -t58 * t31 + t61 * t33;
t12 = -t58 * t30 + t61 * t32;
t10 = -t61 * t30 - t58 * t32;
t7 = t58 * t14 + t61 * t15;
t6 = t61 * t14 - t58 * t15;
t3 = [1, 0, 0, -2 * pkin(1), t82, pkin(1) ^ 2 + qJ(2) ^ 2, t55, t68, 0, 0, 0, t60 * t82, t63 * t82, t54 * t55, -0.2e1 * t55 * t79, t73 * t83, t59 * t68, t53, 0.2e1 * t26 * t60 - 0.2e1 * t55 * t77, -0.2e1 * t27 * t60 - 0.2e1 * t55 * t74, -0.2e1 * t9 * t31 - 0.2e1 * t8 * t33, t40 ^ 2 + t8 ^ 2 + t9 ^ 2, t13 ^ 2, t13 * t85, 0.2e1 * t13 * t60, t60 * t85, t53, 0.2e1 * t1 * t60 + 0.2e1 * t22 * t11, 0.2e1 * t22 * t13 - 0.2e1 * t2 * t60; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70 * t59, t70 * t62, t30 * t33 - t32 * t31, -t8 * t30 + t9 * t32 - t40 * t63, 0, 0, 0, 0, 0, t10 * t60 - t63 * t11, -t12 * t60 - t63 * t13; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30 ^ 2 + t32 ^ 2 + t55, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t63, -t60, 0, t72, -t76, t59 * t50 (-t52 + t54) * t63, t80, t75, 0, t66 * t59 + t62 * t72, -t59 * t72 + t66 * t62, -t24 * t33 - t25 * t31 + t9 * t38 - t8 * t39, t8 * t24 + t9 * t25 + t40 * t51, t13 * t20, -t20 * t11 - t13 * t19, t20 * t60, -t19 * t60, 0, t29 * t11 + t22 * t19 + t6 * t60, t29 * t13 + t22 * t20 - t7 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -t60, 0, 0, 0, 0, 0, t50, -t78, t30 * t39 + t32 * t38, -t30 * t24 + t32 * t25 - t63 * t51, 0, 0, 0, 0, 0, -t63 * t19, -t63 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t52, 0.2e1 * t79, 0, 0, 0, pkin(3) * t83, -0.2e1 * pkin(3) * t59, -0.2e1 * t24 * t39 + 0.2e1 * t25 * t38, t24 ^ 2 + t25 ^ 2 + t51 ^ 2, t20 ^ 2, -0.2e1 * t20 * t19, 0, 0, 0, t19 * t84, t20 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -t78, t60, t26, -t27 (-t31 * t56 - t33 * t57) * pkin(4) (t56 * t9 + t57 * t8) * pkin(4), 0, 0, t13, -t11, t60, t34 * t60 + t1, -t35 * t60 - t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, -t75, 0 (-t30 * t57 + t32 * t56) * pkin(4), 0, 0, 0, 0, 0, t10, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t62, 0, -t59 * pkin(8), -t62 * pkin(8) (t38 * t56 - t39 * t57) * pkin(4) (t24 * t57 + t25 * t56) * pkin(4), 0, 0, t20, -t19, 0, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t56 ^ 2 + t57 ^ 2) * pkin(4) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t34, -0.2e1 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, 0, 0, 0, 0, t11, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, 0, 0, 0, 0, t19, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t11, t60, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, 0, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t34, -t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
