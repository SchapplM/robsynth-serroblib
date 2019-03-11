% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPRPP5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP5_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:06:17
% EndTime: 2019-03-09 10:06:22
% DurationCPUTime: 1.60s
% Computational Cost: add. (1295->216), mult. (2693->345), div. (0->0), fcn. (1898->4), ass. (0->121)
t75 = sin(qJ(4));
t121 = t75 * qJD(5);
t77 = cos(qJ(4));
t128 = t77 * qJ(5);
t141 = pkin(4) + pkin(5);
t143 = t141 * t75 - t128;
t27 = -t143 * qJD(4) + t121;
t149 = pkin(3) + pkin(7);
t116 = qJ(5) * qJD(2);
t76 = sin(qJ(2));
t78 = cos(qJ(2));
t144 = t76 * qJD(5) + t78 * t116;
t129 = t76 * qJ(3);
t80 = -pkin(2) - pkin(8);
t136 = t78 * t80;
t148 = t129 - t136;
t147 = t141 * qJD(2);
t145 = 0.2e1 * t144;
t74 = t78 ^ 2;
t100 = qJD(2) * (t76 ^ 2 - t74);
t71 = t75 ^ 2;
t134 = -t77 ^ 2 + t71;
t99 = t134 * qJD(4);
t131 = qJ(5) * t75;
t87 = -t141 * t77 - t131;
t142 = 2 * qJD(3);
t81 = 0.2e1 * qJD(5);
t120 = t76 * qJD(2);
t91 = -t75 * pkin(4) + t128;
t31 = t91 * qJD(4) + t121;
t59 = pkin(7) * t120;
t92 = pkin(4) * t77 + t131;
t11 = -t59 + t31 * t78 + (-pkin(3) - t92) * t120;
t140 = t11 * t75;
t95 = -t77 * qJD(5) + qJD(3);
t28 = t92 * qJD(4) + t95;
t139 = t28 * t78;
t44 = -pkin(3) * t120 - t59;
t138 = t44 * t75;
t137 = t76 * t80;
t40 = -pkin(1) - t148;
t51 = t149 * t76;
t135 = t77 * t40 + t75 * t51;
t52 = t149 * t78;
t132 = qJ(3) * t78;
t130 = qJ(6) * t78;
t127 = qJ(6) + t80;
t126 = qJD(2) * t75;
t125 = qJD(2) * t77;
t46 = t127 * t75;
t124 = qJD(4) * t46;
t65 = qJD(4) * t75;
t66 = qJD(4) * t77;
t123 = qJD(4) * t78;
t122 = qJD(4) * t80;
t119 = t77 * qJD(6);
t64 = t78 * qJD(2);
t118 = t78 * qJD(3);
t117 = qJ(3) * qJD(4);
t115 = qJD(2) * qJ(3);
t114 = -0.2e1 * pkin(1) * qJD(2);
t101 = -t75 * t40 + t77 * t51;
t14 = t75 * t130 - t141 * t76 - t101;
t16 = t76 * qJ(5) + t135;
t15 = t77 * t130 + t16;
t106 = t77 * t120;
t111 = t75 * t123;
t34 = t106 + t111;
t97 = pkin(2) * t120 - t76 * qJD(3);
t23 = (pkin(8) * t76 - t132) * qJD(2) + t97;
t60 = pkin(7) * t64;
t45 = pkin(3) * t64 + t60;
t8 = -t77 * t23 + t40 * t65 - t75 * t45 - t51 * t66;
t82 = -t34 * qJ(6) + t78 * t119 - t8;
t4 = t82 + t144;
t113 = t14 * t65 + t15 * t66 + t4 * t75;
t112 = pkin(4) * t64;
t110 = t75 * t122;
t109 = t77 * t123;
t108 = t75 * t120;
t107 = t76 * t64;
t105 = t77 * t64;
t104 = t75 * t66;
t103 = t76 * t116;
t102 = qJ(5) * t123;
t98 = t75 * t23 + t40 * t66 - t77 * t45 + t51 * t65;
t96 = t75 * t106;
t10 = t59 - t27 * t78 + (pkin(3) - t87) * t120;
t39 = -qJ(3) - t143;
t94 = t39 * t123 - t10;
t93 = -t78 * pkin(2) - t129;
t17 = -t76 * pkin(4) - t101;
t90 = t16 * t77 + t17 * t75;
t48 = qJ(3) - t91;
t89 = t48 * t78 + t137;
t88 = -t132 - t137;
t86 = qJ(6) * t66 + t75 * qJD(6);
t85 = -qJ(6) * t108 + t98;
t18 = t87 * t78 - t52;
t21 = t87 * qJD(4) - t95;
t84 = -qJD(4) * t18 + t39 * t120 - t21 * t78;
t83 = t93 * qJD(2) + t118;
t5 = -t8 + t144;
t6 = t98 - t112;
t1 = t90 * qJD(4) + t5 * t75 - t6 * t77;
t70 = qJ(5) * t81;
t56 = t77 * t122;
t53 = 0.2e1 * t107;
t50 = t80 * t105;
t49 = -pkin(1) + t93;
t47 = t127 * t77;
t36 = -t108 + t109;
t35 = t75 * t64 + t76 * t66;
t33 = -t76 * t65 + t105;
t32 = -t78 * t115 + t97;
t26 = t56 + t86;
t25 = -t119 + t124;
t22 = t92 * t78 + t52;
t7 = -t25 * t77 + t26 * t75 + (t46 * t77 - t47 * t75) * qJD(4);
t3 = (t86 - t147) * t78 + t85;
t2 = [0, 0, 0, t53, -0.2e1 * t100, 0, 0, 0, t76 * t114, t78 * t114, 0, -0.2e1 * t49 * t120 + 0.2e1 * t32 * t78, -0.2e1 * t32 * t76 - 0.2e1 * t49 * t64, 0.2e1 * t49 * t32, 0.2e1 * t74 * t104 - 0.2e1 * t71 * t107, -0.2e1 * t74 * t99 - 0.4e1 * t78 * t96, 0.2e1 * t75 * t100 - 0.2e1 * t76 * t109, 0.2e1 * t77 * t100 + 0.2e1 * t76 * t111, t53, 0.2e1 * (-t52 * t125 - t98) * t76 + 0.2e1 * (t101 * qJD(2) + t44 * t77 - t52 * t65) * t78, 0.2e1 * (t52 * t126 + t8) * t76 + 0.2e1 * (-t135 * qJD(2) - t52 * t66 - t138) * t78, 0.2e1 * (-t22 * t125 - t6) * t76 + 0.2e1 * (-qJD(2) * t17 + t11 * t77 - t22 * t65) * t78, 0.2e1 * t90 * t120 + 0.2e1 * (-t5 * t77 - t6 * t75 + (t16 * t75 - t17 * t77) * qJD(4)) * t78, 0.2e1 * (-t22 * t126 + t5) * t76 + 0.2e1 * (qJD(2) * t16 + t22 * t66 + t140) * t78, 0.2e1 * t22 * t11 + 0.2e1 * t16 * t5 + 0.2e1 * t17 * t6, 0.2e1 * (t18 * t125 - t3) * t76 + 0.2e1 * (-qJD(2) * t14 - t10 * t77 + t18 * t65) * t78, 0.2e1 * (t18 * t126 + t4) * t76 + 0.2e1 * (qJD(2) * t15 - t10 * t75 - t18 * t66) * t78, 0.2e1 * (-t14 * t75 - t15 * t77) * t120 + 0.2e1 * (t3 * t75 + t4 * t77 + (t14 * t77 - t15 * t75) * qJD(4)) * t78, 0.2e1 * t18 * t10 + 0.2e1 * t14 * t3 + 0.2e1 * t15 * t4; 0, 0, 0, 0, 0, t64, -t120, 0, -t60, t59, t83, t60, -t59, t83 * pkin(7), t78 * t99 + t96, 0.4e1 * t78 * t104 - t134 * t120, t33, -t35, 0, t138 + t50 + (-t76 * t115 + t118) * t77 + (t52 * t77 + t75 * t88) * qJD(4) (qJD(4) * t88 + t44) * t77 + (qJD(2) * t148 - qJD(4) * t52 - t118) * t75, t140 + t50 + (-t48 * t120 + t139) * t77 + (t22 * t77 - t75 * t89) * qJD(4), -t1 (qJD(4) * t89 - t11) * t77 + (qJD(4) * t22 + t139 + (-t48 * t76 + t136) * qJD(2)) * t75, t1 * t80 + t11 * t48 + t22 * t28, -t25 * t76 + t47 * t64 + t75 * t94 + t77 * t84, t26 * t76 + t46 * t64 + t75 * t84 - t77 * t94 (t47 * t120 + (t25 - t124) * t78) * t75 + (-t46 * t120 - t3 + (-qJD(4) * t47 + t26) * t78) * t77 + t113, t10 * t39 + t14 * t25 + t15 * t26 + t18 * t21 - t3 * t47 + t4 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t142, qJ(3) * t142, -0.2e1 * t104, 0.2e1 * t99, 0, 0, 0, 0.2e1 * qJD(3) * t75 + 0.2e1 * t77 * t117, 0.2e1 * qJD(3) * t77 - 0.2e1 * t75 * t117, 0.2e1 * t28 * t75 + 0.2e1 * t48 * t66, 0, -0.2e1 * t28 * t77 + 0.2e1 * t48 * t65, 0.2e1 * t48 * t28, -0.2e1 * t21 * t75 - 0.2e1 * t39 * t66, 0.2e1 * t21 * t77 - 0.2e1 * t39 * t65, 0.2e1 * t7, 0.2e1 * t39 * t21 - 0.2e1 * t47 * t25 + 0.2e1 * t46 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, 0, 0, t60, 0, 0, 0, 0, 0, t33, -t35, t33, 0, t35, t1, t33, t35, 0, -t3 * t77 + t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, t34, t64, -t98, t8, -t98 + 0.2e1 * t112 (-pkin(4) * t120 + t102) * t75 + (t103 + (pkin(4) * qJD(4) - qJD(5)) * t78) * t77, -t8 + t145, -t6 * pkin(4) + t5 * qJ(5) + t16 * qJD(5) (-t86 + 0.2e1 * t147) * t78 - t85, t82 + t145 (t120 * t141 - t102) * t75 + (-t103 + (-qJD(4) * t141 + qJD(5)) * t78) * t77, t4 * qJ(5) + t15 * qJD(5) - t141 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, -t66, 0, -t110, -t56, -t110, -t31, t56, t31 * t80, -t25, t26, t27, t26 * qJ(5) + t46 * qJD(5) - t141 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, -t66, -t65, 0, t66, t31, -t65, t66, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t70, 0, t81, 0, t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t36, 0, t6, -t64, 0, t36, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, 0, t110, 0, 0, t65, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, 0, 0, t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t36, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, -t65, 0, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t2;
