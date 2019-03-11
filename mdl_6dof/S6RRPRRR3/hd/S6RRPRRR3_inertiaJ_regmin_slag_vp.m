% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x33]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRRR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t68 = sin(qJ(5));
t69 = sin(qJ(4));
t72 = cos(qJ(5));
t73 = cos(qJ(4));
t49 = t68 * t69 - t72 * t73;
t66 = cos(pkin(11));
t57 = -t66 * pkin(2) - pkin(3);
t51 = -t73 * pkin(4) + t57;
t35 = t49 * pkin(5) + t51;
t97 = 0.2e1 * t35;
t65 = sin(pkin(11));
t70 = sin(qJ(2));
t87 = cos(qJ(2));
t46 = t65 * t70 - t66 * t87;
t96 = -0.2e1 * t46;
t95 = 0.2e1 * t46;
t94 = 0.2e1 * t51;
t93 = t46 * pkin(4);
t92 = t46 * pkin(5);
t67 = sin(qJ(6));
t91 = t67 * pkin(5);
t90 = t68 * pkin(4);
t71 = cos(qJ(6));
t61 = t71 * pkin(5);
t47 = t65 * t87 + t66 * t70;
t50 = t68 * t73 + t72 * t69;
t23 = t50 * t47;
t60 = -t87 * pkin(2) - pkin(1);
t28 = t46 * pkin(3) - t47 * pkin(8) + t60;
t52 = (-qJ(3) - pkin(7)) * t70;
t76 = t87 * pkin(7);
t53 = t87 * qJ(3) + t76;
t33 = t65 * t52 + t66 * t53;
t16 = t73 * t28 - t69 * t33;
t79 = t73 * t47;
t11 = -pkin(9) * t79 + t16 + t93;
t80 = t73 * t33;
t12 = t80 + (-pkin(9) * t47 + t28) * t69;
t81 = t72 * t12;
t7 = t68 * t11 + t81;
t5 = -t23 * pkin(10) + t7;
t89 = t71 * t5;
t62 = t72 * pkin(4);
t56 = t65 * pkin(2) + pkin(8);
t88 = pkin(9) + t56;
t30 = -t67 * t49 + t71 * t50;
t86 = t30 * t46;
t85 = t50 * t46;
t84 = t69 * t46;
t83 = t69 * t47;
t82 = t69 * t73;
t78 = 0.2e1 * t87;
t77 = t71 * t90;
t24 = t49 * t47;
t6 = t72 * t11 - t68 * t12;
t4 = t24 * pkin(10) + t6 + t92;
t1 = t71 * t4 - t67 * t5;
t42 = t88 * t69;
t43 = t88 * t73;
t26 = -t72 * t42 - t68 * t43;
t31 = -t66 * t52 + t65 * t53;
t59 = t62 + pkin(5);
t38 = t71 * t59 - t67 * t90;
t21 = pkin(4) * t83 + t31;
t2 = t67 * t4 + t89;
t27 = -t68 * t42 + t72 * t43;
t75 = -t46 * t56 + t47 * t57;
t64 = t73 ^ 2;
t63 = t69 ^ 2;
t45 = t47 ^ 2;
t44 = t46 ^ 2;
t40 = t73 * t46;
t39 = t67 * t59 + t77;
t34 = t49 * t46;
t29 = t71 * t49 + t67 * t50;
t20 = t29 * t46;
t19 = -t49 * pkin(10) + t27;
t18 = -t50 * pkin(10) + t26;
t17 = t69 * t28 + t80;
t15 = t23 * pkin(5) + t21;
t14 = -t67 * t23 - t71 * t24;
t13 = t71 * t23 - t67 * t24;
t9 = t67 * t18 + t71 * t19;
t8 = t71 * t18 - t67 * t19;
t3 = [1, 0, 0, t70 ^ 2, t70 * t78, 0, 0, 0, pkin(1) * t78, -0.2e1 * pkin(1) * t70, 0.2e1 * t31 * t47 - 0.2e1 * t33 * t46, t31 ^ 2 + t33 ^ 2 + t60 ^ 2, t64 * t45, -0.2e1 * t45 * t82, t79 * t95, t83 * t96, t44, 0.2e1 * t16 * t46 + 0.2e1 * t31 * t83, -0.2e1 * t17 * t46 + 0.2e1 * t31 * t79, t24 ^ 2, 0.2e1 * t24 * t23, -t24 * t95, t23 * t96, t44, 0.2e1 * t21 * t23 + 0.2e1 * t6 * t46, -0.2e1 * t21 * t24 - 0.2e1 * t7 * t46, t14 ^ 2, -0.2e1 * t14 * t13, t14 * t95, t13 * t96, t44, 0.2e1 * t1 * t46 + 0.2e1 * t15 * t13, 0.2e1 * t15 * t14 - 0.2e1 * t2 * t46; 0, 0, 0, 0, 0, t70, t87, 0, -t70 * pkin(7), -t76 (-t46 * t65 - t47 * t66) * pkin(2) (-t31 * t66 + t33 * t65) * pkin(2), t69 * t79 (-t63 + t64) * t47, t84, t40, 0, -t31 * t73 + t69 * t75, t31 * t69 + t73 * t75, -t24 * t50, -t50 * t23 + t24 * t49, t85, -t34, 0, t21 * t49 + t51 * t23 + t26 * t46, t21 * t50 - t51 * t24 - t27 * t46, t14 * t30, -t30 * t13 - t14 * t29, t86, -t20, 0, t35 * t13 + t15 * t29 + t8 * t46, t35 * t14 + t15 * t30 - t9 * t46; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t65 ^ 2 + t66 ^ 2) * pkin(2) ^ 2, t63, 0.2e1 * t82, 0, 0, 0, -0.2e1 * t57 * t73, 0.2e1 * t57 * t69, t50 ^ 2, -0.2e1 * t50 * t49, 0, 0, 0, t49 * t94, t50 * t94, t30 ^ 2, -0.2e1 * t30 * t29, 0, 0, 0, t29 * t97, t30 * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, 0, 0, 0, 0, t40, -t84, 0, 0, 0, 0, 0, -t34, -t85, 0, 0, 0, 0, 0, -t20, -t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, -t83, t46, t16, -t17, 0, 0, -t24, -t23, t46, t46 * t62 + t6, -t81 + (-t11 - t93) * t68, 0, 0, t14, -t13, t46, t38 * t46 + t1, -t39 * t46 - t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t73, 0, -t69 * t56, -t73 * t56, 0, 0, t50, -t49, 0, t26, -t27, 0, 0, t30, -t29, 0, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, -t69, 0, 0, 0, 0, 0, -t49, -t50, 0, 0, 0, 0, 0, -t29, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t62, -0.2e1 * t90, 0, 0, 0, 0, 1, 0.2e1 * t38, -0.2e1 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t23, t46, t6, -t7, 0, 0, t14, -t13, t46, t46 * t61 + t1, -t89 + (-t4 - t92) * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -t49, 0, t26, -t27, 0, 0, t30, -t29, 0, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t50, 0, 0, 0, 0, 0, -t29, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t62, -t90, 0, 0, 0, 0, 1, t38 + t61, -t77 + (-pkin(5) - t59) * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t61, -0.2e1 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, t46, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t29, 0, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t38, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t61, -t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
