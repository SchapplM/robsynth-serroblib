% Calculate inertial parameters regressor of joint inertia matrix for
% S6PRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRRPR6_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_inertiaJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t75 = sin(qJ(4));
t68 = t75 ^ 2;
t79 = cos(qJ(4));
t70 = t79 ^ 2;
t127 = t68 + t70;
t72 = sin(pkin(6));
t81 = cos(qJ(2));
t108 = t72 * t81;
t77 = sin(qJ(2));
t109 = t72 * t77;
t73 = cos(pkin(6));
t76 = sin(qJ(3));
t80 = cos(qJ(3));
t33 = t80 * t109 + t73 * t76;
t13 = t79 * t108 + t33 * t75;
t15 = -t75 * t108 + t33 * t79;
t93 = t13 * t75 + t15 * t79;
t31 = t76 * t109 - t73 * t80;
t126 = t31 ^ 2;
t78 = cos(qJ(6));
t103 = t78 * t75;
t60 = t79 * t76;
t74 = sin(qJ(6));
t26 = -t76 * t103 + t74 * t60;
t125 = -0.2e1 * t26;
t119 = pkin(4) + pkin(5);
t99 = t75 * qJ(5);
t34 = t119 * t79 + pkin(3) + t99;
t124 = 0.2e1 * t34;
t89 = -t79 * pkin(4) - t99;
t46 = -pkin(3) + t89;
t123 = -0.2e1 * t46;
t122 = -0.2e1 * t76;
t121 = 0.2e1 * t76;
t120 = 0.2e1 * t80;
t118 = pkin(3) * t75;
t117 = pkin(3) * t79;
t69 = t76 ^ 2;
t116 = t69 * pkin(8);
t115 = t75 * pkin(9);
t114 = t76 * pkin(8);
t112 = t31 * t75;
t111 = t31 * t76;
t110 = t31 * t79;
t107 = t75 * t76;
t106 = t75 * t79;
t105 = t75 * t80;
t104 = t76 * t80;
t102 = t79 * t80;
t47 = -t80 * pkin(3) - t76 * pkin(9) - pkin(2);
t101 = pkin(8) * t105 - t79 * t47;
t24 = pkin(8) * t102 + t75 * t47;
t100 = t127 * pkin(9) ^ 2;
t98 = t79 * qJ(5);
t97 = t80 * qJ(5);
t96 = t75 * t104;
t95 = t69 * t106;
t66 = t80 * pkin(4);
t21 = t66 + t101;
t94 = (pkin(9) - pkin(10)) * t75;
t92 = t31 * t107 + t13 * t80;
t91 = t13 ^ 2 + t15 ^ 2 + t126;
t90 = t93 * pkin(9);
t20 = -t97 + t24;
t7 = t80 * pkin(5) - pkin(10) * t60 + t21;
t9 = pkin(10) * t107 + t20;
t1 = t78 * t7 - t74 * t9;
t2 = t74 * t7 + t78 * t9;
t88 = -pkin(4) * t75 + t98;
t87 = t20 * t79 + t21 * t75;
t86 = t101 * t75 + t24 * t79;
t85 = t33 * t80 + t111;
t37 = t74 * t75 + t78 * t79;
t84 = (t13 * t79 - t15 * t75) * t76;
t83 = pkin(8) ^ 2;
t71 = t80 ^ 2;
t67 = t72 ^ 2;
t65 = t79 * pkin(9);
t63 = t69 * t83;
t59 = t70 * t69;
t58 = t68 * t69;
t55 = t67 * t81 ^ 2;
t53 = pkin(9) * t105;
t51 = t75 * t60;
t49 = -t79 * pkin(10) + t65;
t48 = t102 * t122;
t45 = t78 * qJ(5) - t74 * t119;
t43 = t74 * qJ(5) + t78 * t119;
t42 = 0.2e1 * t127 * pkin(9);
t40 = -t74 * t79 + t103;
t39 = (t68 - t70) * t76;
t28 = t37 * t76;
t25 = (pkin(8) - t88) * t76;
t19 = t78 * t49 + t74 * t94;
t17 = t74 * t49 - t78 * t94;
t16 = (-t119 * t75 - pkin(8) + t98) * t76;
t6 = t15 * t80 + t31 * t60;
t5 = t13 * t74 + t15 * t78;
t3 = -t13 * t78 + t15 * t74;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67 * t77 ^ 2 + t73 ^ 2 + t55, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33 ^ 2 + t126 + t55, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 ^ 2 + t5 ^ 2 + t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, -t109, 0, 0, 0, 0, 0, 0, 0, 0, t80 * t108, -t76 * t108, t85, pkin(2) * t108 + t85 * pkin(8), 0, 0, 0, 0, 0, 0, t92, t6, t84, pkin(8) * t111 + t101 * t13 + t15 * t24, 0, 0, 0, 0, 0, 0, t92, t84, -t6, t13 * t21 + t15 * t20 + t31 * t25, 0, 0, 0, 0, 0, 0, -t31 * t26 - t3 * t80, -t31 * t28 - t5 * t80, -t5 * t26 + t3 * t28, -t3 * t1 - t31 * t16 + t5 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t69, 0.2e1 * t104, 0, t71, 0, 0, pkin(2) * t120, pkin(2) * t122, 0.2e1 * (t69 + t71) * pkin(8), pkin(2) ^ 2 + t71 * t83 + t63, t59, -0.2e1 * t95, t48, t58, 0.2e1 * t96, t71, 0.2e1 * t101 * t80 + 0.2e1 * t75 * t116, 0.2e1 * t79 * t116 + 0.2e1 * t24 * t80 (t101 * t79 - t24 * t75) * t121, t101 ^ 2 + t24 ^ 2 + t63, t59, t48, 0.2e1 * t95, t71, -0.2e1 * t96, t58, 0.2e1 * t25 * t107 + 0.2e1 * t21 * t80 (-t20 * t75 + t21 * t79) * t121, -0.2e1 * t20 * t80 - 0.2e1 * t25 * t60, t20 ^ 2 + t21 ^ 2 + t25 ^ 2, t28 ^ 2, t28 * t125, t28 * t120, t26 ^ 2, t80 * t125, t71, 0.2e1 * t1 * t80 + 0.2e1 * t16 * t26, 0.2e1 * t16 * t28 - 0.2e1 * t2 * t80, -0.2e1 * t1 * t28 - 0.2e1 * t2 * t26, t1 ^ 2 + t16 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t33, 0, 0, 0, 0, 0, 0, 0, 0, -t110, t112, t93, -t31 * pkin(3) + t90, 0, 0, 0, 0, 0, 0, -t110, t93, -t112, t31 * t46 + t90, 0, 0, 0, 0, 0, 0, -t31 * t37, -t31 * t40, t3 * t40 - t5 * t37, t3 * t17 + t5 * t19 - t31 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, 0, t80, 0, -t114, -t80 * pkin(8), 0, 0, t51, -t39, -t105, -t51, -t102, 0, t53 + (-pkin(8) * t79 - t118) * t76, pkin(9) * t102 + (pkin(8) * t75 - t117) * t76, t86, -pkin(3) * t114 + t86 * pkin(9), t51, -t105, t39, 0, t102, -t51, t46 * t107 - t25 * t79 + t53, t87, -t25 * t75 + (-pkin(9) * t80 - t46 * t76) * t79, pkin(9) * t87 + t25 * t46, t28 * t40, -t40 * t26 - t28 * t37, t40 * t80, t26 * t37, -t37 * t80, 0, t16 * t37 - t17 * t80 + t34 * t26, t16 * t40 - t19 * t80 + t34 * t28, -t1 * t40 + t17 * t28 - t19 * t26 - t2 * t37, -t1 * t17 + t16 * t34 + t2 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t68, 0.2e1 * t106, 0, t70, 0, 0, 0.2e1 * t117, -0.2e1 * t118, t42, pkin(3) ^ 2 + t100, t68, 0, -0.2e1 * t106, 0, 0, t70, t79 * t123, t42, t75 * t123, t46 ^ 2 + t100, t40 ^ 2, -0.2e1 * t40 * t37, 0, t37 ^ 2, 0, 0, t37 * t124, t40 * t124, 0.2e1 * t17 * t40 - 0.2e1 * t19 * t37, t17 ^ 2 + t19 ^ 2 + t34 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t15, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, t15, -t13 * pkin(4) + t15 * qJ(5), 0, 0, 0, 0, 0, 0, t3, t5, 0, t3 * t43 + t5 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, -t107, -t80, -t101, -t24, 0, 0, 0, t60, 0, -t80, t107, 0, -0.2e1 * t66 - t101, t89 * t76, -0.2e1 * t97 + t24, -t21 * pkin(4) + t20 * qJ(5), 0, 0, -t28, 0, t26, -t80, -t43 * t80 - t1, -t45 * t80 + t2, -t45 * t26 + t43 * t28, -t1 * t43 + t2 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, 0, t79, 0, -t115, -t65, 0, 0, 0, t75, 0, 0, -t79, 0, -t115, t88, t65, t88 * pkin(9), 0, 0, -t40, 0, t37, 0, t17, t19, -t45 * t37 + t43 * t40, t17 * t43 + t19 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t43, 0.2e1 * t45, 0, t43 ^ 2 + t45 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3 * t78 + t5 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t60, 0, t21, 0, 0, 0, 0, 0, 0, t78 * t80, -t74 * t80, -t74 * t26 - t78 * t28, t1 * t78 + t2 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, 0, t115, 0, 0, 0, 0, 0, 0, 0, 0, -t74 * t37 - t78 * t40, -t17 * t78 + t19 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4), 0, 0, 0, 0, 0, 0, -t78, t74, 0, -t43 * t78 + t45 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74 ^ 2 + t78 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t5, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, -t26, t80, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, -t37, 0, -t17, -t19, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t43, -t45, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, -t74, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t4;
