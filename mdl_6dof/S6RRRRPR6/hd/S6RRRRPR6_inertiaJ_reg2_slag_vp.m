% Calculate inertial parameters regressor of joint inertia matrix for
% S6RRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRPR6_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_inertiaJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t90 = sin(qJ(4));
t91 = sin(qJ(3));
t93 = cos(qJ(4));
t94 = cos(qJ(3));
t63 = t90 * t91 - t93 * t94;
t65 = t90 * t94 + t93 * t91;
t87 = sin(pkin(11));
t88 = cos(pkin(11));
t40 = t88 * t63 + t87 * t65;
t78 = -t94 * pkin(3) - pkin(2);
t52 = t63 * pkin(4) + t78;
t28 = t40 * pkin(5) + t52;
t128 = 0.2e1 * t28;
t127 = 0.2e1 * t52;
t126 = 0.2e1 * t78;
t92 = sin(qJ(2));
t125 = -0.2e1 * t92;
t95 = cos(qJ(2));
t124 = -0.2e1 * t95;
t123 = 0.2e1 * t95;
t122 = -pkin(9) - pkin(8);
t121 = pkin(2) * t94;
t120 = pkin(7) * t91;
t84 = t92 ^ 2;
t119 = t84 * pkin(7);
t118 = t87 * pkin(4);
t117 = t90 * pkin(3);
t81 = t92 * pkin(7);
t116 = t95 * pkin(3);
t115 = t95 * pkin(4);
t114 = cos(qJ(6));
t113 = t91 * t92;
t112 = t91 * t94;
t111 = t91 * t95;
t108 = t94 * t95;
t104 = pkin(7) * t108;
t68 = -t95 * pkin(2) - t92 * pkin(8) - pkin(1);
t45 = t104 + (-pkin(9) * t92 + t68) * t91;
t110 = t93 * t45;
t109 = t94 * t92;
t62 = t94 * t68;
t43 = -pkin(9) * t109 + t62 + (-pkin(3) - t120) * t95;
t26 = t93 * t43 - t90 * t45;
t56 = t93 * t109 - t90 * t113;
t17 = -t56 * qJ(5) - t115 + t26;
t27 = t90 * t43 + t110;
t54 = t65 * t92;
t21 = -t54 * qJ(5) + t27;
t9 = t87 * t17 + t88 * t21;
t67 = pkin(3) * t113 + t81;
t83 = t91 ^ 2;
t85 = t94 ^ 2;
t107 = t83 + t85;
t106 = t92 * t123;
t105 = t88 * t117;
t103 = t91 * t109;
t32 = -t87 * t54 + t88 * t56;
t8 = t88 * t17 - t87 * t21;
t6 = -t95 * pkin(5) - t32 * pkin(10) + t8;
t30 = t88 * t54 + t87 * t56;
t7 = -t30 * pkin(10) + t9;
t89 = sin(qJ(6));
t1 = t114 * t6 - t89 * t7;
t82 = t93 * pkin(3);
t77 = t82 + pkin(4);
t59 = t87 * t77 + t105;
t102 = t114 * t59;
t71 = t122 * t91;
t72 = t122 * t94;
t46 = t93 * t71 + t90 * t72;
t36 = -t65 * qJ(5) + t46;
t47 = t90 * t71 - t93 * t72;
t37 = -t63 * qJ(5) + t47;
t19 = t88 * t36 - t87 * t37;
t101 = t114 * t118;
t57 = -t87 * t117 + t88 * t77;
t44 = t54 * pkin(4) + t67;
t20 = t87 * t36 + t88 * t37;
t50 = -pkin(7) * t111 + t62;
t51 = t91 * t68 + t104;
t100 = -t50 * t91 + t51 * t94;
t2 = t114 * t7 + t89 * t6;
t97 = pkin(7) ^ 2;
t86 = t95 ^ 2;
t80 = t84 * t97;
t79 = t88 * pkin(4);
t75 = t79 + pkin(5);
t70 = t114 * t75;
t60 = t89 * t75 + t101;
t58 = -t89 * t118 + t70;
t53 = pkin(5) + t57;
t49 = t114 * t53;
t42 = -t87 * t63 + t88 * t65;
t35 = t89 * t53 + t102;
t34 = -t89 * t59 + t49;
t25 = t114 * t42 - t89 * t40;
t23 = t114 * t40 + t89 * t42;
t22 = t30 * pkin(5) + t44;
t14 = t114 * t32 - t89 * t30;
t12 = t114 * t30 + t89 * t32;
t11 = -t40 * pkin(10) + t20;
t10 = -t42 * pkin(10) + t19;
t4 = t89 * t10 + t114 * t11;
t3 = t114 * t10 - t89 * t11;
t5 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t84, t106, 0, t86, 0, 0, pkin(1) * t123, pkin(1) * t125, 0.2e1 * (t84 + t86) * pkin(7), pkin(1) ^ 2 + t86 * t97 + t80, t85 * t84, -0.2e1 * t84 * t112, t108 * t125, t83 * t84, t91 * t106, t86, 0.2e1 * t91 * t119 - 0.2e1 * t50 * t95, 0.2e1 * t94 * t119 + 0.2e1 * t51 * t95, 0.2e1 * (-t50 * t94 - t51 * t91) * t92, t50 ^ 2 + t51 ^ 2 + t80, t56 ^ 2, -0.2e1 * t56 * t54, t56 * t124, t54 ^ 2, -t54 * t124, t86, -0.2e1 * t26 * t95 + 0.2e1 * t67 * t54, 0.2e1 * t27 * t95 + 0.2e1 * t67 * t56, -0.2e1 * t26 * t56 - 0.2e1 * t27 * t54, t26 ^ 2 + t27 ^ 2 + t67 ^ 2, t32 ^ 2, -0.2e1 * t32 * t30, t32 * t124, t30 ^ 2, t30 * t123, t86, 0.2e1 * t44 * t30 - 0.2e1 * t8 * t95, 0.2e1 * t44 * t32 + 0.2e1 * t9 * t95, -0.2e1 * t9 * t30 - 0.2e1 * t8 * t32, t44 ^ 2 + t8 ^ 2 + t9 ^ 2, t14 ^ 2, -0.2e1 * t14 * t12, t14 * t124, t12 ^ 2, t12 * t123, t86, -0.2e1 * t1 * t95 + 0.2e1 * t22 * t12, 0.2e1 * t22 * t14 + 0.2e1 * t2 * t95, -0.2e1 * t1 * t14 - 0.2e1 * t2 * t12, t1 ^ 2 + t2 ^ 2 + t22 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, 0, t95, 0, -t81, -t95 * pkin(7), 0, 0, t103 (-t83 + t85) * t92, -t111, -t103, -t108, 0, -pkin(7) * t109 + (-pkin(2) * t92 + pkin(8) * t95) * t91, pkin(8) * t108 + (t120 - t121) * t92, t100, -pkin(2) * t81 + t100 * pkin(8), t56 * t65, -t65 * t54 - t56 * t63, -t65 * t95, t54 * t63, t63 * t95, 0, -t46 * t95 + t78 * t54 + t67 * t63, t47 * t95 + t78 * t56 + t67 * t65, -t26 * t65 - t27 * t63 - t46 * t56 - t47 * t54, t26 * t46 + t27 * t47 + t67 * t78, t32 * t42, -t42 * t30 - t32 * t40, -t42 * t95, t30 * t40, t40 * t95, 0, -t19 * t95 + t52 * t30 + t44 * t40, t20 * t95 + t52 * t32 + t44 * t42, -t19 * t32 - t20 * t30 - t9 * t40 - t8 * t42, t8 * t19 + t9 * t20 + t44 * t52, t14 * t25, -t25 * t12 - t14 * t23, -t25 * t95, t12 * t23, t23 * t95, 0, t28 * t12 + t22 * t23 - t3 * t95, t28 * t14 + t22 * t25 + t4 * t95, -t1 * t25 - t4 * t12 - t3 * t14 - t2 * t23, t1 * t3 + t2 * t4 + t22 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t83, 0.2e1 * t112, 0, t85, 0, 0, 0.2e1 * t121, -0.2e1 * pkin(2) * t91, 0.2e1 * t107 * pkin(8), t107 * pkin(8) ^ 2 + pkin(2) ^ 2, t65 ^ 2, -0.2e1 * t65 * t63, 0, t63 ^ 2, 0, 0, t63 * t126, t65 * t126, -0.2e1 * t46 * t65 - 0.2e1 * t47 * t63, t46 ^ 2 + t47 ^ 2 + t78 ^ 2, t42 ^ 2, -0.2e1 * t42 * t40, 0, t40 ^ 2, 0, 0, t40 * t127, t42 * t127, -0.2e1 * t19 * t42 - 0.2e1 * t20 * t40, t19 ^ 2 + t20 ^ 2 + t52 ^ 2, t25 ^ 2, -0.2e1 * t25 * t23, 0, t23 ^ 2, 0, 0, t23 * t128, t25 * t128, -0.2e1 * t4 * t23 - 0.2e1 * t3 * t25, t28 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, 0, -t113, -t95, t50, -t51, 0, 0, 0, 0, t56, 0, -t54, -t95, -t93 * t116 + t26, -t110 + (-t43 + t116) * t90 (-t54 * t90 - t56 * t93) * pkin(3) (t26 * t93 + t27 * t90) * pkin(3), 0, 0, t32, 0, -t30, -t95, -t57 * t95 + t8, t59 * t95 - t9, -t59 * t30 - t57 * t32, t8 * t57 + t9 * t59, 0, 0, t14, 0, -t12, -t95, -t34 * t95 + t1, t35 * t95 - t2, -t35 * t12 - t34 * t14, t1 * t34 + t2 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, 0, t94, 0, -t91 * pkin(8), -t94 * pkin(8), 0, 0, 0, 0, t65, 0, -t63, 0, t46, -t47 (-t63 * t90 - t65 * t93) * pkin(3) (t46 * t93 + t47 * t90) * pkin(3), 0, 0, t42, 0, -t40, 0, t19, -t20, -t59 * t40 - t57 * t42, t19 * t57 + t20 * t59, 0, 0, t25, 0, -t23, 0, t3, -t4, -t35 * t23 - t34 * t25, t3 * t34 + t4 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t82, -0.2e1 * t117, 0 (t90 ^ 2 + t93 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t57, -0.2e1 * t59, 0, t57 ^ 2 + t59 ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t34, -0.2e1 * t35, 0, t34 ^ 2 + t35 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, 0, -t54, -t95, t26, -t27, 0, 0, 0, 0, t32, 0, -t30, -t95, -t115 * t88 + t8, t115 * t87 - t9 (-t30 * t87 - t32 * t88) * pkin(4) (t8 * t88 + t87 * t9) * pkin(4), 0, 0, t14, 0, -t12, -t95, -t58 * t95 + t1, t60 * t95 - t2, -t60 * t12 - t58 * t14, t1 * t58 + t2 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, -t63, 0, t46, -t47, 0, 0, 0, 0, t42, 0, -t40, 0, t19, -t20 (-t40 * t87 - t42 * t88) * pkin(4) (t19 * t88 + t20 * t87) * pkin(4), 0, 0, t25, 0, -t23, 0, t3, -t4, -t60 * t23 - t58 * t25, t3 * t58 + t4 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t82, -t117, 0, 0, 0, 0, 0, 0, 0, 1, t57 + t79, -t105 + (-pkin(4) - t77) * t87, 0 (t57 * t88 + t59 * t87) * pkin(4), 0, 0, 0, 0, 0, 1, t49 + t70 + (-t59 - t118) * t89, -t101 - t102 + (-t53 - t75) * t89, 0, t34 * t58 + t35 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t79, -0.2e1 * t118, 0 (t87 ^ 2 + t88 ^ 2) * pkin(4) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t58, -0.2e1 * t60, 0, t58 ^ 2 + t60 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t32, 0, t44, 0, 0, 0, 0, 0, 0, t12, t14, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t42, 0, t52, 0, 0, 0, 0, 0, 0, t23, t25, 0, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, -t12, -t95, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, -t23, 0, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t34, -t35, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t58, -t60, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t5;