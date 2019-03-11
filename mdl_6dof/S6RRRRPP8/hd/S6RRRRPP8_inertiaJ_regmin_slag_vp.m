% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRPP8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t66 = sin(qJ(4));
t72 = pkin(4) + pkin(5);
t69 = cos(qJ(4));
t95 = t69 * qJ(5);
t131 = t72 * t66 - t95;
t96 = t66 * qJ(5);
t130 = t69 * t72 + t96;
t64 = sin(pkin(6));
t71 = cos(qJ(2));
t108 = t64 * t71;
t68 = sin(qJ(2));
t109 = t64 * t68;
t65 = cos(pkin(6));
t67 = sin(qJ(3));
t70 = cos(qJ(3));
t38 = t70 * t109 + t65 * t67;
t21 = t69 * t108 + t38 * t66;
t129 = -0.2e1 * t21;
t128 = -0.2e1 * t38;
t41 = pkin(3) + t130;
t127 = 0.2e1 * t41;
t82 = -t69 * pkin(4) - t96;
t44 = -pkin(3) + t82;
t126 = -0.2e1 * t44;
t125 = -0.2e1 * t67;
t124 = 0.2e1 * t67;
t123 = 0.2e1 * t70;
t122 = pkin(1) * t68;
t121 = pkin(1) * t71;
t120 = pkin(3) * t66;
t119 = pkin(3) * t69;
t118 = pkin(9) * t66;
t117 = pkin(9) * t69;
t37 = t67 * t109 - t65 * t70;
t36 = t37 * pkin(4);
t116 = t66 * pkin(10);
t115 = pkin(5) + t72;
t50 = pkin(8) * t109;
t30 = t50 + (-pkin(2) - t121) * t65;
t13 = t37 * pkin(3) - t38 * pkin(10) + t30;
t88 = pkin(8) * t108;
t31 = t88 + (pkin(9) + t122) * t65;
t32 = (-pkin(2) * t71 - pkin(9) * t68 - pkin(1)) * t64;
t17 = t70 * t31 + t67 * t32;
t15 = -pkin(10) * t108 + t17;
t7 = t66 * t13 + t69 * t15;
t16 = -t67 * t31 + t70 * t32;
t14 = pkin(3) * t108 - t16;
t114 = t14 * t66;
t113 = t14 * t69;
t112 = t21 * t69;
t22 = -t66 * t108 + t38 * t69;
t111 = t22 * t66;
t60 = t64 ^ 2;
t110 = t60 * t71;
t107 = t65 * t68;
t106 = t66 * t37;
t105 = t66 * t67;
t104 = t66 * t69;
t103 = t66 * t70;
t102 = t67 * t37;
t101 = t69 * t37;
t56 = t69 * t67;
t100 = t69 * t70;
t6 = t69 * t13 - t66 * t15;
t45 = -t70 * pkin(3) - t67 * pkin(10) - pkin(2);
t99 = pkin(9) * t103 - t69 * t45;
t29 = pkin(9) * t100 + t66 * t45;
t61 = t66 ^ 2;
t63 = t69 ^ 2;
t98 = t61 + t63;
t97 = qJ(5) * t21;
t94 = t69 * qJ(6);
t93 = t70 * qJ(5);
t92 = 0.2e1 * t108;
t91 = t67 * t123;
t90 = pkin(10) * t106;
t89 = pkin(10) * t101;
t35 = t37 * qJ(5);
t4 = t35 + t7;
t87 = 0.2e1 * t35 + t7;
t86 = t67 * t108;
t85 = t70 * t108;
t5 = -t36 - t6;
t59 = t70 * pkin(4);
t25 = t59 + t99;
t84 = -0.2e1 * t93 + t29;
t24 = -t93 + t29;
t83 = t4 * t69 + t5 * t66;
t81 = -pkin(4) * t66 + t95;
t80 = t24 * t69 + t25 * t66;
t79 = t22 * qJ(5) - t14;
t78 = t22 * qJ(6) - t5;
t77 = t67 * t94 - t25;
t75 = qJ(5) ^ 2;
t74 = 0.2e1 * qJ(5);
t62 = t67 ^ 2;
t58 = t69 * pkin(10);
t53 = pkin(10) * t103;
t49 = qJ(6) * t105;
t47 = t58 - t94;
t46 = (pkin(10) - qJ(6)) * t66;
t40 = pkin(1) * t107 + t88;
t39 = t65 * t121 - t50;
t33 = (pkin(9) - t81) * t67;
t23 = (-pkin(9) - t131) * t67;
t20 = t21 * qJ(6);
t19 = t24 + t49;
t18 = t70 * pkin(5) - t77;
t8 = t21 * pkin(4) - t79;
t3 = -t72 * t21 + t79;
t2 = t20 + t4;
t1 = -t37 * pkin(5) - t78;
t9 = [1, 0, 0, t60 * t68 ^ 2, 0.2e1 * t68 * t110, 0.2e1 * t64 * t107, t65 * t92, t65 ^ 2, 0.2e1 * pkin(1) * t110 + 0.2e1 * t39 * t65, -0.2e1 * t60 * t122 - 0.2e1 * t40 * t65, t38 ^ 2, t37 * t128, t108 * t128, t37 * t92, t60 * t71 ^ 2, -0.2e1 * t16 * t108 + 0.2e1 * t30 * t37, 0.2e1 * t17 * t108 + 0.2e1 * t30 * t38, t22 ^ 2, t22 * t129, 0.2e1 * t22 * t37, t37 * t129, t37 ^ 2, 0.2e1 * t14 * t21 + 0.2e1 * t6 * t37, 0.2e1 * t14 * t22 - 0.2e1 * t7 * t37, 0.2e1 * t8 * t21 - 0.2e1 * t5 * t37, -0.2e1 * t4 * t21 + 0.2e1 * t5 * t22, -0.2e1 * t8 * t22 + 0.2e1 * t4 * t37, t4 ^ 2 + t5 ^ 2 + t8 ^ 2, -0.2e1 * t1 * t37 - 0.2e1 * t3 * t21, 0.2e1 * t2 * t37 + 0.2e1 * t3 * t22, -0.2e1 * t1 * t22 + 0.2e1 * t2 * t21, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, t109, t108, t65, t39, -t40, t38 * t67, t38 * t70 - t102, -t86, -t85, 0, -pkin(2) * t37 + pkin(9) * t86 - t30 * t70, -pkin(2) * t38 + pkin(9) * t85 + t30 * t67, t22 * t56 (-t111 - t112) * t67, -t22 * t70 + t37 * t56, -t102 * t66 + t21 * t70, -t37 * t70, -t99 * t37 - t6 * t70 + (pkin(9) * t21 + t114) * t67, -t29 * t37 + t7 * t70 + (pkin(9) * t22 + t113) * t67, t105 * t8 + t33 * t21 - t25 * t37 + t5 * t70, -t24 * t21 + t25 * t22 + (-t4 * t66 + t5 * t69) * t67, -t33 * t22 + t24 * t37 - t4 * t70 - t56 * t8, t4 * t24 + t5 * t25 + t8 * t33, t1 * t70 - t105 * t3 - t18 * t37 - t23 * t21, t19 * t37 - t2 * t70 + t23 * t22 + t3 * t56, -t18 * t22 + t19 * t21 + (-t1 * t69 + t2 * t66) * t67, t1 * t18 + t2 * t19 + t3 * t23; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t62, t91, 0, 0, 0, pkin(2) * t123, pkin(2) * t125, t63 * t62, -0.2e1 * t62 * t104, t100 * t125, t66 * t91, t70 ^ 2, 0.2e1 * t118 * t62 + 0.2e1 * t70 * t99, 0.2e1 * t117 * t62 + 0.2e1 * t29 * t70, 0.2e1 * t105 * t33 + 0.2e1 * t25 * t70 (-t24 * t66 + t25 * t69) * t124, -0.2e1 * t24 * t70 - 0.2e1 * t33 * t56, t24 ^ 2 + t25 ^ 2 + t33 ^ 2, -0.2e1 * t105 * t23 + 0.2e1 * t18 * t70, -0.2e1 * t19 * t70 + 0.2e1 * t23 * t56 (-t18 * t69 + t19 * t66) * t124, t18 ^ 2 + t19 ^ 2 + t23 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t37, -t108, t16, -t17, t111, -t66 * t21 + t22 * t69, t106, t101, 0, -pkin(3) * t21 - t113 - t90, -pkin(3) * t22 + t114 - t89, t44 * t21 - t8 * t69 - t90 (t111 - t112) * pkin(10) + t83, -t44 * t22 - t8 * t66 + t89, pkin(10) * t83 + t8 * t44, -t41 * t21 + t3 * t69 - t46 * t37, t41 * t22 + t3 * t66 + t47 * t37, -t1 * t66 - t2 * t69 + t47 * t21 - t46 * t22, t1 * t46 + t2 * t47 + t3 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t70, 0, -t67 * pkin(9), -t70 * pkin(9), t66 * t56 (-t61 + t63) * t67, -t103, -t100, 0, t53 + (-t117 - t120) * t67, pkin(10) * t100 + (t118 - t119) * t67, t105 * t44 - t33 * t69 + t53, t80, -t33 * t66 + (-pkin(10) * t70 - t44 * t67) * t69, pkin(10) * t80 + t33 * t44, -t105 * t41 + t23 * t69 + t46 * t70, t23 * t66 + t41 * t56 - t47 * t70 (-t46 * t67 - t19) * t69 + (t47 * t67 - t18) * t66, t18 * t46 + t19 * t47 + t23 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t61, 0.2e1 * t104, 0, 0, 0, 0.2e1 * t119, -0.2e1 * t120, t69 * t126, 0.2e1 * t98 * pkin(10), t66 * t126, pkin(10) ^ 2 * t98 + t44 ^ 2, t69 * t127, t66 * t127, -0.2e1 * t46 * t66 - 0.2e1 * t47 * t69, t41 ^ 2 + t46 ^ 2 + t47 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t21, t37, t6, -t7, -t5 + t36, -t22 * pkin(4) - t97, t87, -t5 * pkin(4) + t4 * qJ(5), t115 * t37 + t78, t20 + t87, t72 * t22 + t97, t2 * qJ(5) - t1 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, -t105, -t70, -t99, -t29, -0.2e1 * t59 - t99, t82 * t67, t84, -t25 * pkin(4) + t24 * qJ(5), -t115 * t70 + t77, t49 + t84, t130 * t67, t19 * qJ(5) - t18 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t69, 0, -t116, -t58, -t116, t81, t58, t81 * pkin(10), -t46, t47, t131, t47 * qJ(5) - t46 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, t74, pkin(4) ^ 2 + t75, 0.2e1 * t72, t74, 0, t72 ^ 2 + t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t22, 0, t5, -t37, 0, -t22, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, t56, 0, t25, t70, 0, -t56, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, t116, 0, 0, -t66, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4), -1, 0, 0, -t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, t22, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, t56, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t66, 0, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t9;
