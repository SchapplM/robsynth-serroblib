% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPRPP2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:52:29
% EndTime: 2019-03-09 09:52:34
% DurationCPUTime: 1.66s
% Computational Cost: add. (2349->210), mult. (5159->363), div. (0->0), fcn. (4654->6), ass. (0->116)
t80 = cos(qJ(4));
t121 = t80 * qJD(5);
t78 = sin(qJ(4));
t126 = t78 * qJ(5);
t138 = pkin(4) + pkin(5);
t143 = t138 * t80;
t98 = t126 + t143;
t147 = -t98 * qJD(4) + t121;
t76 = sin(pkin(9));
t77 = cos(pkin(9));
t79 = sin(qJ(2));
t81 = cos(qJ(2));
t61 = t76 * t81 + t77 * t79;
t51 = t61 * qJD(2);
t60 = t76 * t79 - t77 * t81;
t144 = t51 * qJ(5) + t60 * qJD(5);
t128 = qJ(5) * t80;
t94 = -t138 * t78 + t128;
t145 = 0.2e1 * t144;
t120 = t80 * qJD(6);
t119 = t81 * qJD(2);
t122 = t79 * qJD(2);
t52 = t77 * t119 - t76 * t122;
t132 = t80 * t52;
t71 = qJD(4) * t78;
t28 = t61 * t71 - t132;
t48 = t51 * pkin(4);
t131 = -qJ(3) - pkin(7);
t109 = qJD(2) * t131;
t49 = t81 * qJD(3) + t79 * t109;
t92 = -t79 * qJD(3) + t81 * t109;
t20 = t77 * t49 + t76 * t92;
t70 = pkin(2) * t122;
t21 = pkin(3) * t51 - pkin(8) * t52 + t70;
t112 = -pkin(2) * t81 - pkin(1);
t31 = pkin(3) * t60 - pkin(8) * t61 + t112;
t63 = t131 * t79;
t64 = t131 * t81;
t36 = t63 * t76 - t64 * t77;
t72 = qJD(4) * t80;
t7 = -t78 * t20 + t80 * t21 - t31 * t71 - t36 * t72;
t4 = -t48 - t7;
t142 = -t28 * qJ(6) + t61 * t120 - t4;
t130 = t78 * t31 + t80 * t36;
t11 = t60 * qJ(5) + t130;
t33 = t78 * t36;
t110 = t80 * t31 - t33;
t12 = -t60 * pkin(4) - t110;
t6 = -t80 * t20 - t78 * t21 - t31 * t72 + t36 * t71;
t3 = -t6 + t144;
t141 = t3 * t78 - t4 * t80 + (t11 * t80 + t12 * t78) * qJD(4);
t137 = pkin(4) * t80;
t104 = t126 + t137;
t140 = t104 * qJD(4) - t121;
t139 = 0.2e1 * qJD(4);
t83 = 0.2e1 * qJD(5);
t67 = pkin(2) * t76 + pkin(8);
t136 = t51 * t67;
t135 = t60 * t67;
t134 = t61 * t78;
t133 = t61 * t80;
t74 = t78 ^ 2;
t75 = t80 ^ 2;
t129 = t74 - t75;
t127 = qJ(6) * t78;
t125 = qJ(6) - t67;
t124 = t78 * qJD(5);
t123 = t78 * qJD(6);
t118 = -0.2e1 * pkin(1) * qJD(2);
t68 = -t77 * pkin(2) - pkin(3);
t117 = t68 * t139;
t116 = t61 * t72;
t115 = t67 * t71;
t114 = t67 * t72;
t113 = t78 * t72;
t111 = -0.4e1 * t78 * t133;
t19 = t76 * t49 - t77 * t92;
t35 = -t77 * t63 - t76 * t64;
t58 = t125 * t80;
t108 = t129 * qJD(4);
t107 = -pkin(4) * t71 + t124;
t97 = -t68 + t126;
t44 = t97 + t143;
t5 = t147 * t61 + t94 * t52 - t19;
t106 = -qJD(4) * t44 * t61 + t5;
t103 = pkin(4) * t78 - t128;
t101 = -t11 * t78 + t12 * t80;
t100 = t52 * t68 - t136;
t99 = -t61 * t68 + t135;
t29 = t51 * t78 + t60 * t72;
t96 = t51 * t80 - t60 * t71;
t95 = t52 * t78 + t116;
t13 = t94 * t61 - t35;
t37 = (-pkin(5) * t78 + t128) * qJD(4) + t107;
t93 = qJD(4) * t13 + t37 * t61 + t44 * t52;
t56 = -t97 - t137;
t8 = t103 * t52 + t140 * t61 + t19;
t91 = -t8 + (t56 * t61 - t135) * qJD(4);
t89 = qJ(6) * t116 + t61 * t123 + t52 * t127 - t6;
t14 = t103 * t61 + t35;
t50 = qJ(5) * t72 + t107;
t87 = qJD(4) * t14 - t50 * t61 + t52 * t56 - t136;
t86 = t94 * qJD(4) + t124;
t1 = -pkin(5) * t51 - t142;
t10 = t61 * t127 + t11;
t2 = t89 + t144;
t9 = t33 + (-qJ(6) * t61 - t31) * t80 - t138 * t60;
t85 = -t1 * t80 + t2 * t78 + (t10 * t80 + t78 * t9) * qJD(4);
t84 = t101 * qJD(4) + t3 * t80 + t4 * t78;
t73 = qJ(5) * t83;
t59 = t61 ^ 2;
t57 = t125 * t78;
t39 = -qJD(4) * t58 - t123;
t38 = t125 * t71 - t120;
t24 = (t74 + t75) * t52;
t15 = [0, 0, 0, 0.2e1 * t79 * t119, 0.2e1 * (-t79 ^ 2 + t81 ^ 2) * qJD(2), 0, 0, 0, t79 * t118, t81 * t118, 0.2e1 * t19 * t61 - 0.2e1 * t20 * t60 + 0.2e1 * t35 * t52 - 0.2e1 * t36 * t51, 0.2e1 * t112 * t70 + 0.2e1 * t35 * t19 + 0.2e1 * t36 * t20, 0.2e1 * t52 * t61 * t75 - 0.2e1 * t59 * t113, t129 * t59 * t139 + t52 * t111, 0.2e1 * t51 * t133 - 0.2e1 * t28 * t60, -0.2e1 * t51 * t134 - 0.2e1 * t95 * t60, 0.2e1 * t60 * t51, 0.2e1 * t110 * t51 + 0.2e1 * t19 * t134 + 0.2e1 * t95 * t35 + 0.2e1 * t7 * t60, -0.2e1 * t130 * t51 + 0.2e1 * t19 * t133 - 0.2e1 * t28 * t35 + 0.2e1 * t6 * t60, -0.2e1 * t12 * t51 + 0.2e1 * t8 * t134 + 0.2e1 * t95 * t14 - 0.2e1 * t4 * t60, 0.2e1 * t101 * t52 - 0.2e1 * t141 * t61, 0.2e1 * t11 * t51 - 0.2e1 * t8 * t133 + 0.2e1 * t28 * t14 + 0.2e1 * t3 * t60, 0.2e1 * t11 * t3 + 0.2e1 * t12 * t4 + 0.2e1 * t14 * t8, -0.2e1 * t1 * t60 - 0.2e1 * t13 * t95 - 0.2e1 * t5 * t134 - 0.2e1 * t9 * t51, 0.2e1 * t10 * t51 - 0.2e1 * t13 * t28 + 0.2e1 * t5 * t133 + 0.2e1 * t2 * t60, 0.2e1 * (t10 * t78 - t80 * t9) * t52 + 0.2e1 * t85 * t61, 0.2e1 * t1 * t9 + 0.2e1 * t10 * t2 + 0.2e1 * t13 * t5; 0, 0, 0, 0, 0, t119, -t122, 0, -pkin(7) * t119, pkin(7) * t122 (-t51 * t76 - t52 * t77) * pkin(2) (-t19 * t77 + t20 * t76) * pkin(2), -t61 * t108 + t78 * t132, qJD(4) * t111 - t129 * t52, t29, t96, 0, -t19 * t80 + t100 * t78 + (t35 * t78 - t99 * t80) * qJD(4), t19 * t78 + t100 * t80 + (t35 * t80 + t99 * t78) * qJD(4), t87 * t78 + t91 * t80, t84, t91 * t78 - t87 * t80, -t14 * t50 + t8 * t56 + t84 * t67, t106 * t80 - t39 * t60 + t51 * t57 - t78 * t93, t106 * t78 + t38 * t60 - t51 * t58 + t80 * t93 (-t39 * t61 + t52 * t57 - t2 + (-t58 * t61 - t9) * qJD(4)) * t80 + (t38 * t61 - t52 * t58 - t1 + (-t57 * t61 + t10) * qJD(4)) * t78, -t1 * t57 + t10 * t38 + t13 * t37 - t2 * t58 + t39 * t9 + t44 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t113, -0.2e1 * t108, 0, 0, 0, t78 * t117, t80 * t117, 0.2e1 * t50 * t80 + 0.2e1 * t56 * t71, 0, 0.2e1 * t50 * t78 - 0.2e1 * t56 * t72, -0.2e1 * t56 * t50, 0.2e1 * t37 * t80 - 0.2e1 * t44 * t71, 0.2e1 * t37 * t78 + 0.2e1 * t44 * t72, -0.2e1 * t38 * t80 - 0.2e1 * t39 * t78 + 0.2e1 * (t57 * t80 - t58 * t78) * qJD(4), 0.2e1 * t37 * t44 - 0.2e1 * t38 * t58 - 0.2e1 * t39 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, 0, 0, 0, 0, 0, t96, -t29, t96, -t24, t29, t141, t96, t29, t24, t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38 * t78 - t39 * t80 + (-t57 * t78 - t58 * t80) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t95, t51, t7, t6, -t4 + t48, -t104 * t52 + (t103 * qJD(4) - t124) * t61, -t6 + t145, -pkin(4) * t4 + qJ(5) * t3 + qJD(5) * t11 (pkin(5) + t138) * t51 + t142, t89 + t145, t52 * t98 + t61 * t86, qJ(5) * t2 + qJD(5) * t10 - t1 * t138; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, -t71, 0, -t114, t115, -t114, -t140, -t115, -t140 * t67, -t39, t38, -t147, qJ(5) * t38 - qJD(5) * t58 - t138 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, -t72, -t71, 0, t72, t50, -t71, t72, 0, t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, t73, 0, t83, 0, t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t28, 0, t4, -t51, 0, t28, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, 0, t114, 0, 0, -t72, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, 0, 0, 0, t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, -t28, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, t72, 0, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t15;
