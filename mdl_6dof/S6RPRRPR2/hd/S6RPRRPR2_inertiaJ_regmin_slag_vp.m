% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRPR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t54 = sin(pkin(11));
t56 = cos(pkin(11));
t59 = sin(qJ(4));
t61 = cos(qJ(4));
t36 = -t54 * t59 + t56 * t61;
t49 = -t61 * pkin(4) - pkin(3);
t27 = -t36 * pkin(5) + t49;
t83 = 0.2e1 * t27;
t62 = cos(qJ(3));
t82 = -0.2e1 * t62;
t81 = 0.2e1 * t62;
t80 = pkin(3) * t61;
t79 = pkin(4) * t54;
t78 = cos(qJ(6));
t37 = t54 * t61 + t56 * t59;
t58 = sin(qJ(6));
t21 = t58 * t36 + t78 * t37;
t77 = t21 * t62;
t55 = sin(pkin(10));
t46 = t55 * pkin(1) + pkin(7);
t76 = t46 * t59;
t60 = sin(qJ(3));
t75 = t59 * t60;
t74 = t59 * t61;
t73 = t59 * t62;
t72 = t61 * t60;
t71 = t61 * t62;
t20 = -t78 * t36 + t58 * t37;
t70 = t62 * t20;
t69 = t62 * t46;
t68 = -qJ(5) - pkin(8);
t57 = cos(pkin(10));
t48 = -t57 * pkin(1) - pkin(2);
t35 = -t62 * pkin(3) - t60 * pkin(8) + t48;
t30 = t61 * t35;
t67 = qJ(5) * t60;
t14 = -t61 * t67 + t30 + (-pkin(4) - t76) * t62;
t65 = t61 * t69;
t18 = t65 + (t35 - t67) * t59;
t9 = t54 * t14 + t56 * t18;
t41 = t68 * t59;
t42 = t68 * t61;
t25 = t54 * t41 - t56 * t42;
t43 = t60 * t46;
t33 = pkin(4) * t75 + t43;
t66 = t60 * t81;
t29 = -t54 * t75 + t56 * t72;
t8 = t56 * t14 - t54 * t18;
t4 = -t62 * pkin(5) - t29 * pkin(9) + t8;
t28 = t37 * t60;
t5 = -t28 * pkin(9) + t9;
t1 = t78 * t4 - t58 * t5;
t24 = t56 * t41 + t54 * t42;
t2 = t58 * t4 + t78 * t5;
t53 = t62 ^ 2;
t52 = t61 ^ 2;
t51 = t60 ^ 2;
t50 = t59 ^ 2;
t47 = t56 * pkin(4) + pkin(5);
t32 = t58 * t47 + t78 * t79;
t31 = t78 * t47 - t58 * t79;
t23 = t59 * t35 + t65;
t22 = -t59 * t69 + t30;
t19 = t28 * pkin(5) + t33;
t16 = t36 * pkin(9) + t25;
t15 = -t37 * pkin(9) + t24;
t13 = -t58 * t28 + t78 * t29;
t12 = t78 * t28 + t58 * t29;
t7 = t58 * t15 + t78 * t16;
t6 = t78 * t15 - t58 * t16;
t3 = [1, 0, 0 (t55 ^ 2 + t57 ^ 2) * pkin(1) ^ 2, t51, t66, 0, 0, 0, t48 * t82, 0.2e1 * t48 * t60, t52 * t51, -0.2e1 * t51 * t74, -0.2e1 * t60 * t71, t59 * t66, t53, -0.2e1 * t22 * t62 + 0.2e1 * t51 * t76, 0.2e1 * t51 * t46 * t61 + 0.2e1 * t23 * t62, -0.2e1 * t9 * t28 - 0.2e1 * t8 * t29, t33 ^ 2 + t8 ^ 2 + t9 ^ 2, t13 ^ 2, -0.2e1 * t13 * t12, t13 * t82, t12 * t81, t53, -0.2e1 * t1 * t62 + 0.2e1 * t19 * t12, 0.2e1 * t19 * t13 + 0.2e1 * t2 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8 * t28 + t9 * t29 - t33 * t62, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 ^ 2 + t29 ^ 2 + t53, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t60, t62, 0, -t43, -t69, t59 * t72 (-t50 + t52) * t60, -t73, -t71, 0, -t46 * t72 + (-pkin(3) * t60 + pkin(8) * t62) * t59, pkin(8) * t71 + (t76 - t80) * t60, -t24 * t29 - t25 * t28 + t9 * t36 - t8 * t37, t8 * t24 + t9 * t25 + t33 * t49, t13 * t21, -t21 * t12 - t13 * t20, -t77, t70, 0, t27 * t12 + t19 * t20 - t6 * t62, t27 * t13 + t19 * t21 + t7 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, -t60, 0, 0, 0, 0, 0, t71, -t73, t28 * t37 + t29 * t36, -t28 * t24 + t29 * t25 - t62 * t49, 0, 0, 0, 0, 0, -t70, -t77; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t50, 0.2e1 * t74, 0, 0, 0, 0.2e1 * t80, -0.2e1 * pkin(3) * t59, -0.2e1 * t24 * t37 + 0.2e1 * t25 * t36, t24 ^ 2 + t25 ^ 2 + t49 ^ 2, t21 ^ 2, -0.2e1 * t21 * t20, 0, 0, 0, t20 * t83, t21 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, -t75, -t62, t22, -t23 (-t28 * t54 - t29 * t56) * pkin(4) (t54 * t9 + t56 * t8) * pkin(4), 0, 0, t13, -t12, -t62, -t31 * t62 + t1, t32 * t62 - t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, -t72, 0 (-t28 * t56 + t29 * t54) * pkin(4), 0, 0, 0, 0, 0, -t12, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t61, 0, -t59 * pkin(8), -t61 * pkin(8) (t36 * t54 - t37 * t56) * pkin(4) (t24 * t56 + t25 * t54) * pkin(4), 0, 0, t21, -t20, 0, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t54 ^ 2 + t56 ^ 2) * pkin(4) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t31, -0.2e1 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, 0, 0, 0, 0, t12, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, 0, 0, 0, 0, t20, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t12, -t62, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t20, 0, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t31, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
