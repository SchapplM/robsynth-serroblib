% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRRPP9
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
% Datum: 2019-03-09 21:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRPP9_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t69 = cos(qJ(4));
t53 = qJ(5) * t69;
t66 = sin(qJ(4));
t130 = -pkin(4) * t66 + t53;
t64 = pkin(4) + qJ(6);
t94 = t66 * qJ(5);
t129 = -t64 * t69 - t94;
t62 = sin(pkin(6));
t71 = cos(qJ(2));
t109 = t62 * t71;
t68 = sin(qJ(2));
t110 = t62 * t68;
t63 = cos(pkin(6));
t67 = sin(qJ(3));
t70 = cos(qJ(3));
t34 = t70 * t110 + t63 * t67;
t20 = t69 * t109 + t34 * t66;
t128 = -0.2e1 * t20;
t127 = -0.2e1 * t34;
t126 = -0.2e1 * t66;
t125 = -0.2e1 * t67;
t124 = 0.2e1 * t67;
t123 = 0.2e1 * t70;
t122 = pkin(1) * t68;
t121 = pkin(1) * t71;
t120 = pkin(3) * t69;
t118 = pkin(9) * t66;
t117 = pkin(10) * t70;
t33 = t67 * t110 - t63 * t70;
t116 = t33 * pkin(4);
t45 = pkin(8) * t110;
t28 = t45 + (-pkin(2) - t121) * t63;
t13 = t33 * pkin(3) - t34 * pkin(10) + t28;
t87 = pkin(8) * t109;
t29 = t87 + (pkin(9) + t122) * t63;
t30 = (-pkin(2) * t71 - pkin(9) * t68 - pkin(1)) * t62;
t17 = t70 * t29 + t67 * t30;
t15 = -pkin(10) * t109 + t17;
t7 = t66 * t13 + t69 * t15;
t16 = -t67 * t29 + t70 * t30;
t14 = pkin(3) * t109 - t16;
t115 = t14 * t66;
t114 = t14 * t69;
t113 = t20 * t69;
t21 = -t66 * t109 + t34 * t69;
t112 = t21 * t66;
t58 = t62 ^ 2;
t111 = t58 * t71;
t108 = t63 * t68;
t106 = t66 * t33;
t105 = t66 * t67;
t104 = t66 * t69;
t103 = t66 * t70;
t102 = t67 * t33;
t101 = t69 * t33;
t51 = t69 * t67;
t100 = t69 * t70;
t99 = -t69 * t13 + t66 * t15;
t55 = t67 * pkin(9);
t98 = pkin(4) * t105 + t55;
t41 = -t70 * pkin(3) - t67 * pkin(10) - pkin(2);
t97 = pkin(9) * t103 - t69 * t41;
t27 = pkin(9) * t100 + t66 * t41;
t59 = t66 ^ 2;
t61 = t69 ^ 2;
t96 = t59 + t61;
t95 = qJ(5) * t20;
t32 = t33 * qJ(5);
t93 = t70 * qJ(5);
t92 = -qJ(6) - t64;
t91 = 0.2e1 * t109;
t90 = t67 * t123;
t89 = pkin(10) * t106;
t88 = pkin(10) * t101;
t86 = t32 + t7;
t85 = t67 * t109;
t84 = t70 * t109;
t57 = t70 * pkin(4);
t24 = t57 + t97;
t83 = -t20 * pkin(5) + t7;
t82 = -t21 * pkin(5) - t99;
t5 = t99 - t116;
t81 = t5 * t66 + t69 * t86;
t79 = -t69 * pkin(4) - t94;
t40 = -pkin(3) + t79;
t80 = -t40 * t67 - t117;
t23 = t93 - t27;
t78 = -t23 * t69 + t24 * t66;
t77 = -pkin(5) * t105 + t27;
t76 = -pkin(5) * t51 - t24;
t75 = -t21 * qJ(5) + t14;
t73 = qJ(5) ^ 2;
t72 = 0.2e1 * qJ(5);
t60 = t67 ^ 2;
t56 = t69 * pkin(10);
t54 = t66 * pkin(10);
t52 = -0.2e1 * t93;
t43 = t69 * pkin(5) + t56;
t42 = t66 * pkin(5) + t54;
t37 = -pkin(3) + t129;
t36 = pkin(1) * t108 + t87;
t35 = t63 * t121 - t45;
t31 = -t67 * t53 + t98;
t22 = (qJ(6) * t66 - t53) * t67 + t98;
t19 = t77 - t93;
t18 = t70 * qJ(6) - t76;
t8 = t20 * pkin(4) + t75;
t3 = t64 * t20 + t75;
t2 = t32 + t83;
t1 = -t64 * t33 - t82;
t4 = [1, 0, 0, t58 * t68 ^ 2, 0.2e1 * t68 * t111, 0.2e1 * t62 * t108, t63 * t91, t63 ^ 2, 0.2e1 * pkin(1) * t111 + 0.2e1 * t35 * t63, -0.2e1 * t58 * t122 - 0.2e1 * t36 * t63, t34 ^ 2, t33 * t127, t109 * t127, t33 * t91, t58 * t71 ^ 2, -0.2e1 * t16 * t109 + 0.2e1 * t28 * t33, 0.2e1 * t17 * t109 + 0.2e1 * t28 * t34, t21 ^ 2, t21 * t128, 0.2e1 * t21 * t33, t33 * t128, t33 ^ 2, 0.2e1 * t14 * t20 - 0.2e1 * t33 * t99, 0.2e1 * t14 * t21 - 0.2e1 * t7 * t33, -0.2e1 * t20 * t86 + 0.2e1 * t5 * t21, -0.2e1 * t8 * t20 + 0.2e1 * t5 * t33, -0.2e1 * t8 * t21 + 0.2e1 * t33 * t86, t5 ^ 2 + t8 ^ 2 + t86 ^ 2, 0.2e1 * t1 * t21 - 0.2e1 * t2 * t20, 0.2e1 * t2 * t33 - 0.2e1 * t3 * t21, -0.2e1 * t1 * t33 + 0.2e1 * t3 * t20, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, t110, t109, t63, t35, -t36, t34 * t67, t34 * t70 - t102, -t85, -t84, 0, -pkin(2) * t33 + pkin(9) * t85 - t28 * t70, -pkin(2) * t34 + pkin(9) * t84 + t28 * t67, t21 * t51 (-t112 - t113) * t67, -t21 * t70 + t33 * t51, -t66 * t102 + t20 * t70, -t33 * t70, -t97 * t33 + t99 * t70 + (pkin(9) * t20 + t115) * t67, -t27 * t33 + t7 * t70 + (pkin(9) * t21 + t114) * t67, t23 * t20 + t24 * t21 + (t5 * t69 - t66 * t86) * t67, -t8 * t105 - t31 * t20 + t24 * t33 - t5 * t70, -t31 * t21 - t23 * t33 - t8 * t51 - t70 * t86, -t23 * t86 + t5 * t24 + t8 * t31, t18 * t21 - t19 * t20 + (t1 * t69 - t2 * t66) * t67, t19 * t33 - t2 * t70 - t22 * t21 - t3 * t51, t1 * t70 + t3 * t105 - t18 * t33 + t22 * t20, t1 * t18 + t2 * t19 + t3 * t22; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t60, t90, 0, 0, 0, pkin(2) * t123, pkin(2) * t125, t61 * t60, -0.2e1 * t60 * t104, t100 * t125, t66 * t90, t70 ^ 2, 0.2e1 * t60 * t118 + 0.2e1 * t70 * t97, 0.2e1 * t60 * pkin(9) * t69 + 0.2e1 * t27 * t70 (t23 * t66 + t24 * t69) * t124, -0.2e1 * t31 * t105 - 0.2e1 * t24 * t70, 0.2e1 * t23 * t70 - 0.2e1 * t31 * t51, t23 ^ 2 + t24 ^ 2 + t31 ^ 2 (t18 * t69 - t19 * t66) * t124, -0.2e1 * t19 * t70 - 0.2e1 * t22 * t51, 0.2e1 * t22 * t105 + 0.2e1 * t18 * t70, t18 ^ 2 + t19 ^ 2 + t22 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t33, -t109, t16, -t17, t112, -t66 * t20 + t21 * t69, t106, t101, 0, -pkin(3) * t20 - t114 - t89, -pkin(3) * t21 + t115 - t88 (t112 - t113) * pkin(10) + t81, -t40 * t20 + t8 * t69 + t89, -t40 * t21 - t8 * t66 + t88, pkin(10) * t81 + t8 * t40, t1 * t66 + t2 * t69 - t43 * t20 + t42 * t21, -t37 * t21 - t3 * t66 + t43 * t33, t37 * t20 - t3 * t69 - t42 * t33, t1 * t42 + t2 * t43 + t3 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t70, 0, -t55, -t70 * pkin(9), t66 * t51 (-t59 + t61) * t67, -t103, -t100, 0, -pkin(9) * t51 + (-pkin(3) * t67 + t117) * t66, pkin(10) * t100 + (t118 - t120) * t67, t78, t31 * t69 + t66 * t80, -t31 * t66 + t69 * t80, pkin(10) * t78 + t31 * t40 (t42 * t67 + t19) * t69 + (-t43 * t67 + t18) * t66, -t22 * t66 - t37 * t51 - t43 * t70, t37 * t105 - t22 * t69 + t42 * t70, t18 * t42 + t19 * t43 + t22 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t59, 0.2e1 * t104, 0, 0, 0, 0.2e1 * t120, pkin(3) * t126, 0.2e1 * t96 * pkin(10), 0.2e1 * t40 * t69, t40 * t126, t96 * pkin(10) ^ 2 + t40 ^ 2, 0.2e1 * t42 * t66 + 0.2e1 * t43 * t69, t37 * t126, -0.2e1 * t37 * t69, t37 ^ 2 + t42 ^ 2 + t43 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t20, t33, -t99, -t7, -t21 * pkin(4) - t95, t99 - 0.2e1 * t116, t86 + t32, -t5 * pkin(4) + qJ(5) * t86, -t64 * t21 - t95, 0.2e1 * t32 + t83 (pkin(4) - t92) * t33 + t82, t2 * qJ(5) - t1 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t105, -t70, -t97, -t27, t79 * t67, 0.2e1 * t57 + t97, t52 + t27, -t24 * pkin(4) - t23 * qJ(5), t129 * t67, t52 + t77, t70 * t92 + t76, t19 * qJ(5) - t18 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t69, 0, -t54, -t56, t130, t54, t56, t130 * pkin(10), -t64 * t66 + t53, t43, -t42, t43 * qJ(5) - t42 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(4), t72, pkin(4) ^ 2 + t73, 0, t72, 0.2e1 * t64, t64 ^ 2 + t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t33, 0, t5, t21, 0, -t33, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t70, 0, t24, t51, 0, t70, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, 0, t54, t66, 0, 0, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(4), 0, 0, -1, -t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, t33, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, -t70, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, 0, 0, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t4;
