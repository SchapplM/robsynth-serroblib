% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% tauc_reg [5x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP13_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP13_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:59:37
% EndTime: 2019-12-31 18:59:41
% DurationCPUTime: 1.27s
% Computational Cost: add. (1436->234), mult. (2976->332), div. (0->0), fcn. (1627->4), ass. (0->121)
t59 = sin(qJ(3));
t106 = t59 * qJD(1);
t53 = qJD(4) + t106;
t60 = cos(qJ(4));
t105 = t60 * qJD(3);
t58 = sin(qJ(4));
t109 = qJD(4) * t58;
t61 = cos(qJ(3));
t93 = t61 * t109;
t67 = t59 * t105 + t93;
t19 = t67 * qJD(1) - qJD(4) * t105;
t123 = t61 * t19;
t124 = t60 * t61;
t130 = t58 * t53;
t107 = t58 * qJD(3);
t113 = qJD(1) * t61;
t90 = t60 * t113;
t43 = t90 + t107;
t83 = qJD(4) * t59 + qJD(1);
t149 = (t59 * (-t43 + t90) + t53 * t124) * qJD(3) - t83 * t130 - t123;
t97 = 0.2e1 * qJD(1);
t102 = qJD(1) * qJD(3);
t87 = t61 * t102;
t82 = pkin(4) * t87;
t108 = qJD(4) * t60;
t80 = pkin(3) * t61 + pkin(7) * t59;
t40 = t80 * qJD(3) + qJD(2);
t27 = t40 * qJD(1);
t47 = t59 * pkin(3) - t61 * pkin(7) + qJ(2);
t33 = t47 * qJD(1);
t62 = -pkin(1) - pkin(6);
t146 = qJD(1) * t62;
t52 = qJD(2) + t146;
t46 = t59 * t52;
t35 = qJD(3) * pkin(7) + t46;
t95 = t61 * t107;
t86 = t35 * t108 + t33 * t109 - t60 * t27 + t52 * t95;
t2 = -t82 + t86;
t12 = t58 * t33 + t60 * t35;
t8 = t53 * qJ(5) + t12;
t148 = -t8 * t53 + t2;
t128 = t59 * t60;
t147 = t62 * t128 + t58 * t47;
t110 = qJD(4) * t43;
t129 = t58 * t59;
t50 = t102 * t129;
t20 = -t50 + t110;
t143 = t43 ^ 2;
t112 = qJD(3) * t59;
t3 = t20 * pkin(4) + t19 * qJ(5) - t43 * qJD(5) + t52 * t112;
t142 = t3 * t58;
t141 = t3 * t60;
t122 = t61 * t52;
t36 = -qJD(3) * pkin(3) - t122;
t41 = t58 * t113 - t105;
t10 = t41 * pkin(4) - t43 * qJ(5) + t36;
t139 = t10 * t43;
t138 = t19 * t58;
t137 = t19 * t59;
t136 = t20 * t59;
t135 = t36 * t58;
t134 = t41 * t53;
t133 = t43 * t41;
t132 = t43 * t53;
t131 = t53 * t60;
t127 = t60 * t40;
t45 = t80 * qJD(1);
t126 = t60 * t45;
t125 = t60 * t47;
t63 = qJD(3) ^ 2;
t121 = t63 * t59;
t120 = t63 * t61;
t75 = pkin(4) * t58 - qJ(5) * t60;
t119 = t58 * qJD(5) - t53 * t75 + t46;
t118 = t60 * t122 + t58 * t45;
t57 = t61 ^ 2;
t116 = t59 ^ 2 - t57;
t64 = qJD(1) ^ 2;
t115 = -t63 - t64;
t114 = t64 * qJ(2);
t111 = qJD(3) * t61;
t11 = t60 * t33 - t58 * t35;
t104 = qJD(5) - t11;
t103 = qJ(2) * qJD(3);
t101 = pkin(7) * t130;
t100 = pkin(7) * t131;
t94 = t61 * t105;
t98 = t47 * t108 + t58 * t40 + t62 * t94;
t96 = pkin(7) * t111;
t92 = t62 * t109;
t91 = t53 * t108;
t89 = qJD(2) * t97;
t88 = t58 * t62 - pkin(4);
t85 = t41 + t105;
t84 = -t43 + t107;
t81 = qJ(5) * t87;
t7 = -t53 * pkin(4) + t104;
t78 = t58 * t8 - t60 * t7;
t77 = t58 * t7 + t60 * t8;
t76 = t60 * pkin(4) + t58 * qJ(5);
t74 = qJD(1) * t57 - t53 * t59;
t72 = -t62 + t75;
t71 = -t10 * t59 + t96;
t70 = t36 * t59 - t96;
t69 = t12 * t53 - t86;
t68 = -t33 * t108 + t35 * t109 - t58 * t27 - t52 * t94;
t1 = t53 * qJD(5) - t68 + t81;
t66 = -t78 * qJD(4) + t1 * t60 + t2 * t58;
t65 = t41 * t112 + (-t20 - t50) * t61 + (-t83 * t60 - t95) * t53;
t48 = -pkin(3) - t76;
t23 = t72 * t61;
t18 = t88 * t59 - t125;
t17 = t59 * qJ(5) + t147;
t15 = t43 * pkin(4) + t41 * qJ(5);
t14 = -t126 + (-pkin(4) * qJD(1) + t52 * t58) * t61;
t13 = qJ(5) * t113 + t118;
t9 = t134 - t19;
t6 = (t76 * qJD(4) - qJD(5) * t60) * t61 - t72 * t112;
t5 = qJD(4) * t147 + t88 * t111 - t127;
t4 = qJ(5) * t111 + (qJD(5) - t92) * t59 + t98;
t16 = [0, 0, 0, 0, t89, qJ(2) * t89, -0.2e1 * t59 * t87, 0.2e1 * t116 * t102, -t121, -t120, 0, -t62 * t121 + (qJD(2) * t59 + t61 * t103) * t97, -t62 * t120 + (qJD(2) * t61 - t59 * t103) * t97, -t60 * t123 - t67 * t43, (t41 * t60 + t43 * t58) * t112 + (t138 - t20 * t60 + (t41 * t58 - t43 * t60) * qJD(4)) * t61, -t53 * t93 - t137 + (t43 * t61 + t74 * t60) * qJD(3), -t61 * t91 - t136 + (-t41 * t61 - t74 * t58) * qJD(3), (t53 + t106) * t111, t53 * t127 - t86 * t59 - t61 * t62 * t20 + (t36 * t124 - t147 * t53) * qJD(4) + ((t62 * t41 - t135) * t59 + (qJD(1) * t125 + t11 + (-t62 * t53 + (t52 - t146) * t59) * t58) * t61) * qJD(3), -(-t59 * t92 + t98) * t53 + t68 * t59 + (-t36 * t109 + t62 * t19) * t61 + ((-qJD(1) * t147 - t12) * t61 + (t62 * t43 + (-t36 + t122) * t60) * t59) * qJD(3), t23 * t20 + t6 * t41 - t5 * t53 + (-t10 * t107 - t2) * t59 + (t10 * t108 + t142 + (-qJD(1) * t18 - t7) * qJD(3)) * t61, -t17 * t20 - t18 * t19 - t4 * t41 + t5 * t43 + t78 * t112 + (-qJD(4) * t77 - t1 * t58 + t2 * t60) * t61, t23 * t19 + t4 * t53 - t6 * t43 + (t10 * t105 + t1) * t59 + (t10 * t109 - t141 + (qJD(1) * t17 + t8) * qJD(3)) * t61, t1 * t17 + t10 * t6 + t2 * t18 + t3 * t23 + t8 * t4 + t7 * t5; 0, 0, 0, 0, -t64, -t114, 0, 0, 0, 0, 0, t115 * t59, t115 * t61, 0, 0, 0, 0, 0, t65, -t149, t65, (-t41 * t111 + t43 * t83 - t136) * t60 + (t43 * t111 + t41 * t83 - t137) * t58, t149, -t78 * qJD(1) + (qJD(3) * t77 - t3) * t61 + (qJD(3) * t10 + t66) * t59; 0, 0, 0, 0, 0, 0, t61 * t64 * t59, -t116 * t64, 0, 0, 0, -t61 * t114, t59 * t114, t43 * t131 - t138, (-t19 - t134) * t60 + (-t20 - t132) * t58, t91 + (t53 * t128 + t84 * t61) * qJD(1), -t53 * t109 + (-t53 * t129 + t61 * t85) * qJD(1), -t53 * t113, -t53 * t126 - pkin(3) * t20 + (t61 * t130 - t59 * t85) * t52 + (-t100 + t135) * qJD(4) + (-t11 * t61 + t58 * t70) * qJD(1), pkin(3) * t19 + t118 * t53 + t84 * t46 + (t36 * t60 + t101) * qJD(4) + (t12 * t61 + t60 * t70) * qJD(1), t14 * t53 + t48 * t20 - t141 - t119 * t41 + (t10 * t58 - t100) * qJD(4) + (-t58 * t71 + t61 * t7) * qJD(1), t13 * t41 - t14 * t43 + (t1 + t53 * t7 + (-t20 + t110) * pkin(7)) * t60 + ((qJD(4) * t41 - t19) * pkin(7) + t148) * t58, -t13 * t53 + t48 * t19 - t142 + t119 * t43 + (-t10 * t60 - t101) * qJD(4) + (t60 * t71 - t61 * t8) * qJD(1), t66 * pkin(7) - t119 * t10 - t8 * t13 - t7 * t14 + t3 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133, -t41 ^ 2 + t143, t9, t132 - t20, t87, -t36 * t43 + t69, t11 * t53 + t36 * t41 + t68, -t15 * t41 - t139 + t69 + 0.2e1 * t82, pkin(4) * t19 - t20 * qJ(5) + (-t12 + t8) * t43 + (t7 - t104) * t41, 0.2e1 * t81 - t10 * t41 + t15 * t43 + (0.2e1 * qJD(5) - t11) * t53 - t68, -t2 * pkin(4) + t1 * qJ(5) - t10 * t15 + t104 * t8 - t7 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87 + t133, t9, -t53 ^ 2 - t143, t139 + t148;];
tauc_reg = t16;
