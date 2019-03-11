% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x33]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPRR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t113 = cos(qJ(5));
t79 = sin(pkin(6));
t85 = sin(qJ(2));
t110 = t79 * t85;
t81 = cos(pkin(6));
t84 = sin(qJ(3));
t87 = cos(qJ(3));
t60 = t110 * t84 - t81 * t87;
t61 = t110 * t87 + t81 * t84;
t78 = sin(pkin(12));
t80 = cos(pkin(12));
t40 = -t60 * t78 + t61 * t80;
t83 = sin(qJ(5));
t92 = t60 * t80 + t61 * t78;
t26 = t113 * t92 + t40 * t83;
t121 = -0.2e1 * t26;
t74 = -t87 * pkin(3) - pkin(2);
t91 = t78 * t84 - t80 * t87;
t51 = pkin(4) * t91 + t74;
t120 = 0.2e1 * t51;
t119 = 0.2e1 * t87;
t118 = pkin(1) * t85;
t88 = cos(qJ(2));
t109 = t79 * t88;
t98 = pkin(8) * t109;
t53 = t98 + (pkin(9) + t118) * t81;
t54 = (-pkin(2) * t88 - pkin(9) * t85 - pkin(1)) * t79;
t35 = -t53 * t84 + t54 * t87;
t32 = -pkin(3) * t109 - qJ(4) * t61 + t35;
t36 = t53 * t87 + t54 * t84;
t34 = -qJ(4) * t60 + t36;
t15 = t32 * t80 - t34 * t78;
t11 = -pkin(4) * t109 - pkin(10) * t40 + t15;
t16 = t32 * t78 + t34 * t80;
t14 = -pkin(10) * t92 + t16;
t95 = -t11 * t113 + t14 * t83;
t4 = pkin(5) * t109 + t95;
t86 = cos(qJ(6));
t117 = t4 * t86;
t116 = t78 * pkin(3);
t115 = t88 * pkin(1);
t73 = pkin(3) * t80 + pkin(4);
t57 = t113 * t73 - t116 * t83;
t55 = -pkin(5) - t57;
t114 = pkin(5) - t55;
t27 = t113 * t40 - t83 * t92;
t82 = sin(qJ(6));
t19 = -t109 * t82 + t27 * t86;
t17 = t19 * t82;
t102 = -qJ(4) - pkin(9);
t67 = t102 * t84;
t68 = t102 * t87;
t49 = t67 * t78 - t68 * t80;
t39 = -pkin(10) * t91 + t49;
t48 = t67 * t80 + t68 * t78;
t64 = t78 * t87 + t80 * t84;
t90 = -pkin(10) * t64 + t48;
t22 = -t113 * t90 + t39 * t83;
t112 = t22 * t86;
t75 = t79 ^ 2;
t111 = t75 * t88;
t108 = t81 * t85;
t24 = t82 * t26;
t46 = t113 * t91 + t64 * t83;
t42 = t82 * t46;
t47 = t113 * t64 - t83 * t91;
t107 = t82 * t47;
t58 = t113 * t116 + t73 * t83;
t56 = pkin(11) + t58;
t106 = t82 * t56;
t105 = t82 * t86;
t25 = t86 * t26;
t104 = t86 * t47;
t103 = t86 * t56;
t101 = -0.2e1 * t47 * t46;
t100 = -0.2e1 * t109;
t99 = 0.2e1 * t109;
t97 = t84 * t109;
t96 = t87 * t109;
t94 = -pkin(5) * t47 - pkin(11) * t46;
t93 = -t46 * t56 + t47 * t55;
t69 = pkin(8) * t110;
t52 = t69 + (-pkin(2) - t115) * t81;
t7 = t11 * t83 + t113 * t14;
t44 = t60 * pkin(3) + t52;
t28 = pkin(4) * t92 + t44;
t77 = t86 ^ 2;
t76 = t82 ^ 2;
t71 = t75 * t88 ^ 2;
t70 = 0.2e1 * t105;
t63 = pkin(1) * t108 + t98;
t62 = t115 * t81 - t69;
t45 = t47 ^ 2;
t43 = t86 * t46;
t41 = t82 * t104;
t29 = (-t76 + t77) * t47;
t23 = t113 * t39 + t83 * t90;
t21 = t46 * pkin(5) - t47 * pkin(11) + t51;
t20 = t22 * t82;
t18 = t109 * t86 + t27 * t82;
t13 = t21 * t82 + t23 * t86;
t12 = t21 * t86 - t23 * t82;
t9 = -t18 * t82 + t19 * t86;
t8 = t26 * pkin(5) - t27 * pkin(11) + t28;
t5 = -pkin(11) * t109 + t7;
t3 = t4 * t82;
t2 = t5 * t86 + t8 * t82;
t1 = -t5 * t82 + t8 * t86;
t6 = [1, 0, 0, t75 * t85 ^ 2, 0.2e1 * t85 * t111, 0.2e1 * t79 * t108, t81 * t99, t81 ^ 2, 0.2e1 * pkin(1) * t111 + 0.2e1 * t62 * t81, -0.2e1 * t118 * t75 - 0.2e1 * t63 * t81, t61 ^ 2, -0.2e1 * t61 * t60, t61 * t100, t60 * t99, t71, -0.2e1 * t109 * t35 + 0.2e1 * t52 * t60, 0.2e1 * t109 * t36 + 0.2e1 * t52 * t61, -0.2e1 * t15 * t40 - 0.2e1 * t16 * t92, t15 ^ 2 + t16 ^ 2 + t44 ^ 2, t27 ^ 2, t27 * t121, t27 * t100, t26 * t99, t71, 0.2e1 * t109 * t95 + 0.2e1 * t26 * t28, 0.2e1 * t109 * t7 + 0.2e1 * t27 * t28, t19 ^ 2, -0.2e1 * t19 * t18, 0.2e1 * t19 * t26, t18 * t121, t26 ^ 2, 0.2e1 * t1 * t26 + 0.2e1 * t18 * t4, 0.2e1 * t19 * t4 - 0.2e1 * t2 * t26; 0, 0, 0, 0, 0, t110, t109, t81, t62, -t63, t61 * t84, -t60 * t84 + t61 * t87, -t97, -t96, 0, -pkin(2) * t60 + pkin(9) * t97 - t52 * t87, -pkin(2) * t61 + pkin(9) * t96 + t52 * t84, -t15 * t64 - t16 * t91 - t48 * t40 - t49 * t92, t15 * t48 + t16 * t49 + t44 * t74, t27 * t47, -t26 * t47 - t27 * t46, -t47 * t109, t46 * t109, 0, t109 * t22 + t26 * t51 + t28 * t46, t109 * t23 + t27 * t51 + t28 * t47, t19 * t104 (-t18 * t86 - t17) * t47, t104 * t26 + t19 * t46, -t107 * t26 - t18 * t46, t26 * t46, t1 * t46 + t107 * t4 + t12 * t26 + t18 * t22, t104 * t4 - t13 * t26 + t19 * t22 - t2 * t46; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t84 ^ 2, t84 * t119, 0, 0, 0, pkin(2) * t119, -0.2e1 * pkin(2) * t84, -0.2e1 * t48 * t64 - 0.2e1 * t49 * t91, t48 ^ 2 + t49 ^ 2 + t74 ^ 2, t45, t101, 0, 0, 0, t46 * t120, t47 * t120, t77 * t45, -0.2e1 * t45 * t105, 0.2e1 * t46 * t104, t82 * t101, t46 ^ 2, 0.2e1 * t107 * t22 + 0.2e1 * t12 * t46, 0.2e1 * t104 * t22 - 0.2e1 * t13 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t60, -t109, t35, -t36 (-t80 * t40 - t78 * t92) * pkin(3) (t15 * t80 + t16 * t78) * pkin(3), 0, 0, t27, -t26, -t109, -t109 * t57 - t95, t109 * t58 - t7, t17, t9, t24, t25, 0, -t106 * t26 + t18 * t55 - t117, -t103 * t26 + t19 * t55 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, t87, 0, -t84 * pkin(9), -t87 * pkin(9) (-t80 * t64 - t78 * t91) * pkin(3) (t48 * t80 + t49 * t78) * pkin(3), 0, 0, t47, -t46, 0, -t22, -t23, t41, t29, t42, t43, 0, t82 * t93 - t112, t86 * t93 + t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t78 ^ 2 + t80 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t57, -0.2e1 * t58, t76, t70, 0, 0, 0, -0.2e1 * t55 * t86, 0.2e1 * t55 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, 0, 0, 0, 0, t26, t27, 0, 0, 0, 0, 0, t25, -t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0, 0, 0, 0, 0, t46, t47, 0, 0, 0, 0, 0, t43, -t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t26, -t109, -t95, -t7, t17, t9, t24, t25, 0, -pkin(5) * t18 - pkin(11) * t24 - t117, -pkin(5) * t19 - pkin(11) * t25 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t46, 0, -t22, -t23, t41, t29, t42, t43, 0, t82 * t94 - t112, t86 * t94 + t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t57, -t58, t76, t70, 0, 0, 0, t114 * t86, -t114 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t76, t70, 0, 0, 0, 0.2e1 * pkin(5) * t86, -0.2e1 * pkin(5) * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, t26, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, -t107, t46, t12, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, t86, 0, -t106, -t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, -t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, t86, 0, -t82 * pkin(11), -t86 * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t6;
