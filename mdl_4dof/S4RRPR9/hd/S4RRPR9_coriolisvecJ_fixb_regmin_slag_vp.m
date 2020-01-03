% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% tauc_reg [4x21]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRPR9_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR9_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:10:04
% EndTime: 2019-12-31 17:10:07
% DurationCPUTime: 0.95s
% Computational Cost: add. (871->192), mult. (2370->303), div. (0->0), fcn. (1607->6), ass. (0->121)
t89 = sin(qJ(2));
t131 = qJD(1) * t89;
t86 = sin(pkin(7));
t118 = t86 * t131;
t87 = cos(pkin(7));
t125 = t87 * qJD(2);
t57 = t118 - t125;
t117 = t87 * t131;
t126 = t86 * qJD(2);
t59 = t117 + t126;
t88 = sin(qJ(4));
t90 = cos(qJ(4));
t102 = t88 * t57 - t90 * t59;
t91 = cos(qJ(2));
t124 = t91 * qJD(1);
t78 = -qJD(4) + t124;
t149 = t102 * t78;
t122 = qJD(1) * qJD(2);
t148 = -0.2e1 * t122;
t147 = -qJD(4) - t78;
t146 = t102 * qJD(4);
t145 = t86 * t91;
t144 = t87 * t89;
t143 = t87 * t91;
t142 = t88 * t86;
t141 = t90 * t87;
t93 = qJD(1) ^ 2;
t140 = t91 * t93;
t92 = qJD(2) ^ 2;
t139 = t92 * t89;
t138 = t92 * t91;
t137 = pkin(6) + qJ(3);
t61 = -t141 + t142;
t99 = t61 * t91;
t136 = qJD(1) * t99 - t61 * qJD(4);
t62 = t90 * t86 + t88 * t87;
t100 = t62 * t91;
t135 = -qJD(1) * t100 + t62 * qJD(4);
t105 = pkin(2) * t89 - qJ(3) * t91;
t45 = t105 * qJD(2) - t89 * qJD(3);
t36 = t45 * qJD(1);
t81 = pkin(5) * t131;
t67 = (qJD(3) - t81) * qJD(2);
t14 = t86 * t36 + t87 * t67;
t69 = -t91 * pkin(2) - t89 * qJ(3) - pkin(1);
t51 = t69 * qJD(1);
t123 = qJD(2) * qJ(3);
t82 = pkin(5) * t124;
t73 = t82 + t123;
t23 = t86 * t51 + t87 * t73;
t114 = t91 * t122;
t127 = qJD(4) * t90;
t134 = t114 * t141 - t57 * t127;
t64 = t105 * qJD(1);
t28 = pkin(5) * t118 + t87 * t64;
t130 = qJD(2) * t89;
t120 = pkin(5) * t130;
t24 = t86 * t120 + t87 * t45;
t76 = pkin(5) * t143;
t33 = t86 * t69 + t76;
t133 = t89 ^ 2 - t91 ^ 2;
t132 = qJD(2) * pkin(2);
t129 = qJD(2) * t91;
t128 = qJD(4) * t89;
t121 = pkin(5) * t145;
t119 = pkin(3) * t124;
t116 = pkin(3) * t86 + pkin(5);
t13 = t87 * t36 - t86 * t67;
t101 = pkin(3) * t89 - pkin(6) * t143;
t95 = t101 * qJD(2);
t6 = qJD(1) * t95 + t13;
t108 = t86 * t114;
t8 = -pkin(6) * t108 + t14;
t115 = t90 * t6 - t88 * t8;
t113 = t89 * t122;
t22 = t87 * t51 - t86 * t73;
t112 = t57 + t125;
t111 = -t59 + t126;
t110 = pkin(1) * t148;
t109 = qJD(3) - t132;
t107 = t88 * t6 + t90 * t8;
t7 = -t59 * pkin(6) - t119 + t22;
t9 = -t57 * pkin(6) + t23;
t1 = t90 * t7 - t88 * t9;
t2 = t88 * t7 + t90 * t9;
t68 = t109 + t81;
t106 = t109 - t68;
t56 = t87 * t69;
t21 = -pkin(6) * t144 + t56 + (-pkin(5) * t86 - pkin(3)) * t91;
t26 = -t86 * t89 * pkin(6) + t33;
t104 = t90 * t21 - t88 * t26;
t103 = t88 * t21 + t90 * t26;
t98 = -pkin(5) * t144 - pkin(6) * t145;
t72 = t137 * t87;
t97 = t101 * qJD(1) + qJD(3) * t86 + qJD(4) * t72 + t28;
t49 = t86 * t64;
t71 = t137 * t86;
t96 = -t98 * qJD(1) + qJD(3) * t87 - qJD(4) * t71 - t49;
t94 = qJD(2) * t100;
t4 = (-qJD(4) * t59 - t108) * t88 + t134;
t5 = qJD(1) * t94 - t146;
t80 = -t87 * pkin(3) - pkin(2);
t77 = pkin(5) * t114;
t65 = t116 * t89;
t53 = t116 * t129;
t52 = t86 * t119 + t82;
t46 = t90 * t57;
t44 = pkin(3) * t108 + t77;
t42 = t61 * t89;
t41 = t62 * t89;
t37 = t86 * t45;
t32 = t56 - t121;
t29 = -pkin(5) * t117 + t49;
t27 = t57 * pkin(3) + t68;
t25 = -t87 * t120 + t37;
t17 = t88 * t59 + t46;
t16 = t98 * qJD(2) + t37;
t12 = t127 * t144 - t128 * t142 + t94;
t11 = -qJD(2) * t99 - t62 * t128;
t10 = t95 + t24;
t3 = [0, 0, 0, 0.2e1 * t91 * t113, t133 * t148, t138, -t139, 0, -pkin(5) * t138 + t89 * t110, pkin(5) * t139 + t91 * t110, (-qJD(1) * t24 - t13) * t91 + ((pkin(5) * t57 + t68 * t86) * t91 + (t22 + (t32 + 0.2e1 * t121) * qJD(1)) * t89) * qJD(2), (qJD(1) * t25 + t14) * t91 + ((pkin(5) * t59 + t68 * t87) * t91 + (-t23 + (-t33 + 0.2e1 * t76) * qJD(1)) * t89) * qJD(2), -t24 * t59 - t25 * t57 + (-t13 * t87 - t14 * t86) * t89 + (-t22 * t87 - t23 * t86 + (-t32 * t87 - t33 * t86) * qJD(1)) * t129, t13 * t32 + t14 * t33 + t22 * t24 + t23 * t25 + (t68 + t81) * pkin(5) * t129, -t102 * t11 - t4 * t42, t102 * t12 - t11 * t17 - t4 * t41 + t42 * t5, -t11 * t78 - t4 * t91 + (-qJD(1) * t42 - t102) * t130, t12 * t78 + t5 * t91 + (-qJD(1) * t41 - t17) * t130, (-t78 - t124) * t130, -(t90 * t10 - t88 * t16) * t78 - t115 * t91 + t53 * t17 + t65 * t5 + t44 * t41 + t27 * t12 + (t103 * t78 + t2 * t91) * qJD(4) + (t104 * qJD(1) + t1) * t130, (t88 * t10 + t90 * t16) * t78 + t107 * t91 - t53 * t102 + t65 * t4 - t44 * t42 + t27 * t11 + (t1 * t91 + t104 * t78) * qJD(4) + (-t103 * qJD(1) - t2) * t130; 0, 0, 0, -t89 * t140, t133 * t93, 0, 0, 0, t93 * pkin(1) * t89, pkin(1) * t140, ((-t86 * t123 - t22) * t89 + (-t112 * pkin(5) + t106 * t86 + t28) * t91) * qJD(1), ((-t87 * t123 + t23) * t89 + (t111 * pkin(5) + t106 * t87 - t29) * t91) * qJD(1), t28 * t59 + t29 * t57 + (-qJD(3) * t57 + t22 * t124 + t14) * t87 + (qJD(3) * t59 + t23 * t124 - t13) * t86, -t22 * t28 - t23 * t29 + (-t22 * t86 + t23 * t87) * qJD(3) + (-t13 * t86 + t14 * t87) * qJ(3) + (-t68 - t132) * t82, -t102 * t136 + t4 * t62, t102 * t135 - t136 * t17 - t4 * t61 - t62 * t5, -t136 * t78 + (qJD(2) * t62 + t102) * t131, t135 * t78 + (-qJD(2) * t61 + t17) * t131, t78 * t131, -t52 * t17 + t44 * t61 + t80 * t5 + (t96 * t88 + t97 * t90) * t78 + t135 * t27 + ((-t90 * t71 - t88 * t72) * qJD(2) - t1) * t131, t52 * t102 + t80 * t4 + t44 * t62 + (-t97 * t88 + t96 * t90) * t78 + t136 * t27 + (-(-t88 * t71 + t90 * t72) * qJD(2) + t2) * t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111 * t124, t112 * t124, -t57 ^ 2 - t59 ^ 2, t22 * t59 + t23 * t57 + t77, 0, 0, 0, 0, 0, t5 + t149, t46 * t78 + (-t108 + (-qJD(4) + t78) * t59) * t88 + t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102 * t17, t102 ^ 2 - t17 ^ 2, -t17 * t78 + t4, -t62 * t114 + t146 + t149, t113, t102 * t27 + t147 * t2 + t115, t147 * t1 + t27 * t17 - t107;];
tauc_reg = t3;
