% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% tauc_reg [5x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPPR7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:36:22
% EndTime: 2019-12-31 19:36:26
% DurationCPUTime: 1.13s
% Computational Cost: add. (1377->212), mult. (3622->288), div. (0->0), fcn. (2505->6), ass. (0->128)
t91 = sin(pkin(8));
t92 = cos(pkin(8));
t94 = sin(qJ(2));
t96 = cos(qJ(2));
t72 = t91 * t96 + t92 * t94;
t64 = t72 * qJD(1);
t154 = qJD(5) + t64;
t93 = sin(qJ(5));
t127 = t93 * qJD(2);
t137 = t92 * t96;
t121 = qJD(1) * t137;
t126 = t94 * qJD(1);
t61 = t91 * t126 - t121;
t95 = cos(qJ(5));
t39 = -t95 * t61 + t127;
t155 = t154 * t39;
t115 = t154 ^ 2;
t122 = qJD(1) * qJD(2);
t119 = t94 * t122;
t54 = qJD(2) * t121 - t91 * t119;
t47 = t95 * t54;
t107 = -t93 * t115 + t47;
t153 = -0.2e1 * t122;
t60 = t64 ^ 2;
t152 = -t61 ^ 2 - t60;
t132 = -qJ(3) - pkin(6);
t78 = t132 * t96;
t75 = qJD(1) * t78;
t67 = t91 * t75;
t77 = t132 * t94;
t74 = qJD(1) * t77;
t36 = t92 * t74 + t67;
t125 = -qJD(4) + t36;
t151 = -qJD(5) + t154;
t146 = t61 * pkin(4);
t138 = t92 * t75;
t70 = qJD(2) * pkin(2) + t74;
t34 = t91 * t70 - t138;
t31 = -qJD(2) * qJ(4) - t34;
t15 = -t31 - t146;
t35 = t91 * t74 - t138;
t86 = -t92 * pkin(2) - pkin(3);
t82 = -pkin(7) + t86;
t150 = t82 * t54 + (t15 - t35 + t146) * t154;
t148 = pkin(3) + pkin(7);
t63 = t72 * qJD(2);
t53 = qJD(1) * t63;
t147 = t53 * pkin(3);
t145 = t64 * pkin(4);
t120 = -t96 * pkin(2) - pkin(1);
t109 = -t72 * qJ(4) + t120;
t71 = t91 * t94 - t137;
t17 = t148 * t71 + t109;
t144 = t17 * t54;
t116 = qJD(2) * t132;
t58 = t96 * qJD(3) + t94 * t116;
t48 = t58 * qJD(1);
t59 = -t94 * qJD(3) + t96 * t116;
t49 = t59 * qJD(1);
t18 = t91 * t48 - t92 * t49;
t37 = -t92 * t77 - t91 * t78;
t143 = t18 * t37;
t128 = qJD(5) * t95;
t22 = -qJD(5) * t127 + t61 * t128 + t93 * t53;
t142 = t22 * t95;
t41 = t95 * qJD(2) + t93 * t61;
t141 = t41 * t61;
t140 = t61 * t39;
t139 = t71 * t93;
t136 = t93 * t54;
t98 = qJD(1) ^ 2;
t135 = t96 * t98;
t97 = qJD(2) ^ 2;
t134 = t97 * t94;
t133 = t97 * t96;
t19 = t92 * t48 + t91 * t49;
t131 = t94 ^ 2 - t96 ^ 2;
t130 = qJD(2) * t94;
t129 = qJD(5) * t71;
t124 = t145 - t125;
t88 = pkin(2) * t130;
t26 = t91 * t58 - t92 * t59;
t33 = t92 * t70 + t67;
t83 = pkin(2) * t119;
t118 = -t54 * qJ(4) + t83;
t117 = pkin(2) * t126 + t61 * qJ(4);
t114 = t154 * t41;
t113 = pkin(1) * t153;
t112 = qJD(4) - t33;
t10 = -t148 * qJD(2) + t112 + t145;
t111 = t120 * qJD(1);
t76 = qJD(3) + t111;
t100 = -t64 * qJ(4) + t76;
t9 = t148 * t61 + t100;
t1 = t95 * t10 - t93 * t9;
t2 = t93 * t10 + t95 * t9;
t27 = t92 * t58 + t91 * t59;
t38 = t91 * t77 - t92 * t78;
t16 = -qJD(2) * qJD(4) - t19;
t24 = t61 * pkin(3) + t100;
t108 = t24 * t64 + t18;
t106 = t71 * t128 + t93 * t63;
t105 = -t64 * qJD(4) + t118;
t66 = qJD(2) * t137 - t91 * t130;
t104 = -t66 * qJ(4) - t72 * qJD(4) + t88;
t5 = -t53 * pkin(4) - t16;
t103 = t5 + (-qJD(5) * t82 + t148 * t64 + t117) * t154;
t28 = t72 * pkin(4) + t37;
t102 = t15 * t63 - t28 * t54 + t5 * t71;
t101 = -t95 * t115 - t136;
t99 = t18 * t72 + t26 * t64 - t27 * t61 + t37 * t54 - t38 * t53;
t84 = t91 * pkin(2) + qJ(4);
t56 = qJD(2) * t61;
t46 = t95 * t53;
t32 = t71 * pkin(3) + t109;
t30 = -qJD(2) * pkin(3) + t112;
t29 = -t71 * pkin(4) + t38;
t25 = t64 * pkin(3) + t117;
t23 = t41 * qJD(5) - t46;
t14 = t63 * pkin(3) + t104;
t13 = -t63 * pkin(4) + t27;
t12 = t66 * pkin(4) + t26;
t8 = t105 + t147;
t7 = t54 * pkin(4) + t18;
t6 = t95 * t7;
t4 = t148 * t63 + t104;
t3 = t148 * t53 + t105;
t11 = [0, 0, 0, 0.2e1 * t96 * t119, t131 * t153, t133, -t134, 0, -pkin(6) * t133 + t94 * t113, pkin(6) * t134 + t96 * t113, -t19 * t71 - t33 * t66 - t34 * t63 + t99, t143 + t19 * t38 - t33 * t26 + t34 * t27 + (t76 + t111) * t88, t16 * t71 + t30 * t66 + t31 * t63 + t99, t26 * qJD(2) - t14 * t61 - t24 * t63 - t32 * t53 - t8 * t71, t27 * qJD(2) - t14 * t64 - t24 * t66 - t32 * t54 - t8 * t72, t24 * t14 - t16 * t38 + t30 * t26 - t31 * t27 + t8 * t32 + t143, t106 * t41 + t22 * t139, (-t39 * t93 + t41 * t95) * t63 + (t142 - t23 * t93 + (-t39 * t95 - t41 * t93) * qJD(5)) * t71, t106 * t154 + t71 * t136 + t22 * t72 + t41 * t66, t71 * t47 - t23 * t72 - t39 * t66 + (-t93 * t129 + t95 * t63) * t154, t154 * t66 + t54 * t72, t1 * t66 + t13 * t39 + t29 * t23 + t6 * t72 + (-t154 * t4 - t3 * t72 - t144) * t93 + (t12 * t154 - t102) * t95 + ((-t95 * t17 - t93 * t28) * t154 - t2 * t72 + t15 * t139) * qJD(5), t13 * t41 - t2 * t66 + t29 * t22 + (-(qJD(5) * t28 + t4) * t154 - t144 - (qJD(5) * t10 + t3) * t72 + t15 * t129) * t95 + (-(-qJD(5) * t17 + t12) * t154 - (-qJD(5) * t9 + t7) * t72 + t102) * t93; 0, 0, 0, -t94 * t135, t131 * t98, 0, 0, 0, t98 * pkin(1) * t94, pkin(1) * t135, (t34 - t35) * t64 + (-t33 + t36) * t61 + (-t53 * t91 - t54 * t92) * pkin(2), t33 * t35 - t34 * t36 + (-t76 * t126 - t18 * t92 + t19 * t91) * pkin(2), -t84 * t53 + t86 * t54 + (-t31 - t35) * t64 + (t30 + t125) * t61, -t35 * qJD(2) + t25 * t61 + t108, -t24 * t61 + t25 * t64 + (0.2e1 * qJD(4) - t36) * qJD(2) + t19, t125 * t31 - t16 * t84 + t18 * t86 - t24 * t25 - t30 * t35, -t93 * t114 + t142, (-t23 - t114) * t95 + (-t22 + t155) * t93, t107 + t141, t101 - t140, t154 * t61, t1 * t61 + t103 * t93 + t124 * t39 + t150 * t95 + t84 * t23, t103 * t95 + t124 * t41 - t150 * t93 - t2 * t61 + t84 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152, t33 * t64 + t34 * t61 + t83, t152, -0.2e1 * t64 * qJD(2), -t54 + t56, t147 - t31 * t61 + (-qJD(4) - t30) * t64 + t118, 0, 0, 0, 0, 0, t101 + t140, -t107 + t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54 + t56, -t64 * t61, -t60 - t97, t31 * qJD(2) + t108, 0, 0, 0, 0, 0, -qJD(2) * t39 + t107, -qJD(2) * t41 + t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41 * t39, -t39 ^ 2 + t41 ^ 2, t22 + t155, t151 * t41 + t46, t54, -t15 * t41 + t151 * t2 - t93 * t3 + t6, t151 * t1 + t15 * t39 - t95 * t3 - t93 * t7;];
tauc_reg = t11;
