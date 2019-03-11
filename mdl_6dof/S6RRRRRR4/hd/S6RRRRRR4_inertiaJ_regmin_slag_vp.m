% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% MM_reg [((6+1)*6/2)x38]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRRR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR4_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t80 = sin(qJ(4));
t81 = sin(qJ(3));
t85 = cos(qJ(4));
t86 = cos(qJ(3));
t56 = t80 * t81 - t85 * t86;
t57 = t80 * t86 + t85 * t81;
t79 = sin(qJ(5));
t84 = cos(qJ(5));
t35 = t84 * t56 + t79 * t57;
t69 = -t86 * pkin(3) - pkin(2);
t46 = t56 * pkin(4) + t69;
t24 = t35 * pkin(5) + t46;
t116 = 0.2e1 * t24;
t115 = 0.2e1 * t46;
t114 = 0.2e1 * t69;
t82 = sin(qJ(2));
t113 = -0.2e1 * t82;
t87 = cos(qJ(2));
t112 = -0.2e1 * t87;
t111 = 0.2e1 * t87;
t110 = pkin(8) + pkin(9);
t109 = pkin(2) * t86;
t108 = pkin(7) * t81;
t78 = sin(qJ(6));
t107 = t78 * pkin(5);
t106 = t79 * pkin(4);
t105 = t80 * pkin(3);
t47 = t57 * t82;
t48 = t56 * t82;
t26 = t84 * t47 - t79 * t48;
t102 = t87 * pkin(4);
t60 = -t87 * pkin(2) - t82 * pkin(8) - pkin(1);
t55 = t86 * t60;
t94 = t86 * t82;
t37 = -pkin(9) * t94 + t55 + (-pkin(3) - t108) * t87;
t93 = t86 * t87;
t89 = pkin(7) * t93;
t39 = t89 + (-pkin(9) * t82 + t60) * t81;
t22 = t85 * t37 - t80 * t39;
t15 = t48 * pkin(10) - t102 + t22;
t95 = t85 * t39;
t23 = t80 * t37 + t95;
t18 = -t47 * pkin(10) + t23;
t96 = t84 * t18;
t9 = t79 * t15 + t96;
t7 = -t26 * pkin(11) + t9;
t83 = cos(qJ(6));
t104 = t83 * t7;
t103 = t87 * pkin(3);
t101 = t87 * pkin(5);
t100 = t81 * t82;
t99 = t81 * t86;
t98 = t81 * t87;
t73 = t85 * pkin(3);
t68 = t73 + pkin(4);
t90 = t84 * t105;
t53 = t79 * t68 + t90;
t97 = t83 * t53;
t70 = t82 * pkin(7);
t59 = pkin(3) * t100 + t70;
t92 = t82 * t111;
t91 = t83 * t106;
t27 = -t79 * t47 - t84 * t48;
t8 = t84 * t15 - t79 * t18;
t6 = -t27 * pkin(11) - t101 + t8;
t1 = t83 * t6 - t78 * t7;
t88 = -t53 - t106;
t61 = t110 * t81;
t62 = t110 * t86;
t40 = -t85 * t61 - t80 * t62;
t31 = -t57 * pkin(10) + t40;
t41 = -t80 * t61 + t85 * t62;
t32 = -t56 * pkin(10) + t41;
t16 = t84 * t31 - t79 * t32;
t51 = -t79 * t105 + t84 * t68;
t49 = pkin(5) + t51;
t45 = t83 * t49;
t29 = -t78 * t53 + t45;
t72 = t84 * pkin(4);
t67 = t72 + pkin(5);
t63 = t83 * t67;
t50 = -t78 * t106 + t63;
t38 = t47 * pkin(4) + t59;
t2 = t78 * t6 + t104;
t17 = t79 * t31 + t84 * t32;
t77 = t87 ^ 2;
t76 = t86 ^ 2;
t75 = t82 ^ 2;
t74 = t81 ^ 2;
t71 = t83 * pkin(5);
t52 = t78 * t67 + t91;
t44 = t81 * t60 + t89;
t43 = -pkin(7) * t98 + t55;
t36 = -t79 * t56 + t84 * t57;
t30 = t78 * t49 + t97;
t21 = -t78 * t35 + t83 * t36;
t20 = t83 * t35 + t78 * t36;
t19 = t26 * pkin(5) + t38;
t14 = -t78 * t26 + t83 * t27;
t13 = t83 * t26 + t78 * t27;
t11 = -t35 * pkin(11) + t17;
t10 = -t36 * pkin(11) + t16;
t4 = t78 * t10 + t83 * t11;
t3 = t83 * t10 - t78 * t11;
t5 = [1, 0, 0, t75, t92, 0, 0, 0, pkin(1) * t111, pkin(1) * t113, t76 * t75, -0.2e1 * t75 * t99, t93 * t113, t81 * t92, t77, 0.2e1 * t75 * t108 - 0.2e1 * t43 * t87, 0.2e1 * t75 * pkin(7) * t86 + 0.2e1 * t44 * t87, t48 ^ 2, 0.2e1 * t48 * t47, -t48 * t112, t47 * t111, t77, -0.2e1 * t22 * t87 + 0.2e1 * t59 * t47, 0.2e1 * t23 * t87 - 0.2e1 * t59 * t48, t27 ^ 2, -0.2e1 * t27 * t26, t27 * t112, t26 * t111, t77, 0.2e1 * t38 * t26 - 0.2e1 * t8 * t87, 0.2e1 * t38 * t27 + 0.2e1 * t9 * t87, t14 ^ 2, -0.2e1 * t14 * t13, t14 * t112, t13 * t111, t77, -0.2e1 * t1 * t87 + 0.2e1 * t19 * t13, 0.2e1 * t19 * t14 + 0.2e1 * t2 * t87; 0, 0, 0, 0, 0, t82, t87, 0, -t70, -t87 * pkin(7), t81 * t94 (-t74 + t76) * t82, -t98, -t93, 0, -pkin(7) * t94 + (-pkin(2) * t82 + pkin(8) * t87) * t81, pkin(8) * t93 + (t108 - t109) * t82, -t48 * t57, -t57 * t47 + t48 * t56, -t57 * t87, t56 * t87, 0, -t40 * t87 + t69 * t47 + t59 * t56, t41 * t87 - t69 * t48 + t59 * t57, t27 * t36, -t36 * t26 - t27 * t35, -t36 * t87, t35 * t87, 0, -t16 * t87 + t46 * t26 + t38 * t35, t17 * t87 + t46 * t27 + t38 * t36, t14 * t21, -t21 * t13 - t14 * t20, -t21 * t87, t20 * t87, 0, t24 * t13 + t19 * t20 - t3 * t87, t24 * t14 + t19 * t21 + t4 * t87; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t74, 0.2e1 * t99, 0, 0, 0, 0.2e1 * t109, -0.2e1 * pkin(2) * t81, t57 ^ 2, -0.2e1 * t57 * t56, 0, 0, 0, t56 * t114, t57 * t114, t36 ^ 2, -0.2e1 * t36 * t35, 0, 0, 0, t35 * t115, t36 * t115, t21 ^ 2, -0.2e1 * t21 * t20, 0, 0, 0, t20 * t116, t21 * t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, -t100, -t87, t43, -t44, 0, 0, -t48, -t47, -t87, -t85 * t103 + t22, -t95 + (-t37 + t103) * t80, 0, 0, t27, -t26, -t87, -t51 * t87 + t8, t53 * t87 - t9, 0, 0, t14, -t13, -t87, -t29 * t87 + t1, t30 * t87 - t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t86, 0, -t81 * pkin(8), -t86 * pkin(8), 0, 0, t57, -t56, 0, t40, -t41, 0, 0, t36, -t35, 0, t16, -t17, 0, 0, t21, -t20, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t73, -0.2e1 * t105, 0, 0, 0, 0, 1, 0.2e1 * t51, -0.2e1 * t53, 0, 0, 0, 0, 1, 0.2e1 * t29, -0.2e1 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t47, -t87, t22, -t23, 0, 0, t27, -t26, -t87, -t84 * t102 + t8, -t96 + (-t15 + t102) * t79, 0, 0, t14, -t13, -t87, -t50 * t87 + t1, t52 * t87 - t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t56, 0, t40, -t41, 0, 0, t36, -t35, 0, t16, -t17, 0, 0, t21, -t20, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t73, -t105, 0, 0, 0, 0, 1, t51 + t72, -t90 + (-pkin(4) - t68) * t79, 0, 0, 0, 0, 1, t88 * t78 + t45 + t63, t88 * t83 + (-t49 - t67) * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t72, -0.2e1 * t106, 0, 0, 0, 0, 1, 0.2e1 * t50, -0.2e1 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t26, -t87, t8, -t9, 0, 0, t14, -t13, -t87, -t83 * t101 + t1, -t104 + (-t6 + t101) * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t35, 0, t16, -t17, 0, 0, t21, -t20, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t51, -t53, 0, 0, 0, 0, 1, t29 + t71, -t97 + (-pkin(5) - t49) * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t72, -t106, 0, 0, 0, 0, 1, t50 + t71, -t91 + (-pkin(5) - t67) * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t71, -0.2e1 * t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, -t87, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t20, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t29, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t50, -t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t71, -t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t5;
