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
% tauc_reg [5x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:00
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% StartTime: 2021-01-15 19:59:01
% EndTime: 2021-01-15 19:59:06
% DurationCPUTime: 1.27s
% Computational Cost: add. (1445->227), mult. (3816->307), div. (0->0), fcn. (2625->6), ass. (0->134)
t92 = sin(pkin(8));
t93 = cos(pkin(8));
t95 = sin(qJ(2));
t97 = cos(qJ(2));
t72 = t92 * t97 + t93 * t95;
t104 = qJD(1) * t72;
t159 = qJD(5) + t104;
t94 = sin(qJ(5));
t127 = t94 * qJD(2);
t139 = t93 * t97;
t122 = qJD(1) * t139;
t131 = qJD(1) * t95;
t61 = t131 * t92 - t122;
t96 = cos(qJ(5));
t39 = -t61 * t96 + t127;
t160 = t159 * t39;
t116 = t159 ^ 2;
t123 = qJD(1) * qJD(2);
t120 = t97 * t123;
t121 = t95 * t123;
t79 = t92 * t121;
t54 = t120 * t93 - t79;
t47 = t96 * t54;
t110 = -t116 * t94 + t47;
t158 = -0.2e1 * t123;
t60 = t104 ^ 2;
t157 = -t61 ^ 2 - t60;
t156 = 0.2e1 * t104;
t134 = -qJ(3) - pkin(6);
t78 = t134 * t97;
t75 = qJD(1) * t78;
t67 = t92 * t75;
t77 = t134 * t95;
t74 = qJD(1) * t77;
t36 = t74 * t93 + t67;
t126 = -qJD(4) + t36;
t155 = -qJD(5) + t159;
t149 = t61 * pkin(4);
t140 = t93 * t75;
t70 = qJD(2) * pkin(2) + t74;
t34 = t70 * t92 - t140;
t31 = -qJD(2) * qJ(4) - t34;
t15 = -t31 - t149;
t35 = t74 * t92 - t140;
t86 = -pkin(2) * t93 - pkin(3);
t82 = -pkin(7) + t86;
t154 = t82 * t54 + (t15 - t35 + t149) * t159;
t152 = pkin(3) + pkin(7);
t151 = pkin(2) * t95;
t63 = t72 * qJD(2);
t53 = qJD(1) * t63;
t150 = t53 * pkin(3);
t148 = t104 * pkin(4);
t87 = -t97 * pkin(2) - pkin(1);
t111 = -t72 * qJ(4) + t87;
t71 = t92 * t95 - t139;
t17 = t152 * t71 + t111;
t147 = t17 * t54;
t117 = qJD(2) * t134;
t58 = t97 * qJD(3) + t117 * t95;
t48 = t58 * qJD(1);
t59 = -t95 * qJD(3) + t117 * t97;
t49 = t59 * qJD(1);
t18 = t48 * t92 - t49 * t93;
t37 = -t77 * t93 - t78 * t92;
t146 = t18 * t37;
t128 = qJD(5) * t96;
t22 = -qJD(5) * t127 + t128 * t61 + t53 * t94;
t145 = t22 * t96;
t132 = qJD(1) * t87;
t76 = qJD(3) + t132;
t101 = -qJ(4) * t104 + t76;
t24 = t61 * pkin(3) + t101;
t144 = t24 * t104;
t41 = qJD(2) * t96 + t61 * t94;
t143 = t41 * t61;
t142 = t61 * t39;
t141 = t71 * t94;
t138 = t94 * t54;
t99 = qJD(1) ^ 2;
t137 = t97 * t99;
t98 = qJD(2) ^ 2;
t136 = t98 * t95;
t135 = t98 * t97;
t19 = t48 * t93 + t49 * t92;
t133 = t95 ^ 2 - t97 ^ 2;
t130 = qJD(2) * t95;
t129 = qJD(5) * t71;
t125 = t148 - t126;
t89 = pkin(2) * t130;
t88 = pkin(2) * t131;
t26 = t58 * t92 - t59 * t93;
t33 = t93 * t70 + t67;
t83 = pkin(2) * t121;
t119 = -t54 * qJ(4) + t83;
t118 = t61 * qJ(4) + t88;
t115 = t159 * t41;
t114 = pkin(1) * t158;
t113 = qJD(4) - t33;
t10 = -qJD(2) * t152 + t113 + t148;
t9 = t152 * t61 + t101;
t1 = t10 * t96 - t9 * t94;
t2 = t10 * t94 + t9 * t96;
t27 = t58 * t93 + t59 * t92;
t38 = t77 * t92 - t78 * t93;
t16 = -qJD(2) * qJD(4) - t19;
t109 = t128 * t71 + t63 * t94;
t108 = -qJD(2) * t35 + t18;
t107 = -qJD(4) * t104 + t119;
t66 = qJD(2) * t139 - t130 * t92;
t106 = -t66 * qJ(4) - t72 * qJD(4) + t89;
t5 = -pkin(4) * t53 - t16;
t105 = t5 + (-qJD(5) * t82 + t104 * t152 + t118) * t159;
t28 = pkin(4) * t72 + t37;
t103 = t15 * t63 - t28 * t54 + t5 * t71;
t102 = -t116 * t96 - t138;
t100 = t104 * t26 + t18 * t72 - t27 * t61 + t37 * t54 - t38 * t53;
t84 = pkin(2) * t92 + qJ(4);
t56 = qJD(2) * t61;
t46 = t96 * t53;
t32 = pkin(3) * t71 + t111;
t30 = -qJD(2) * pkin(3) + t113;
t29 = -pkin(4) * t71 + t38;
t25 = pkin(3) * t104 + t118;
t23 = qJD(5) * t41 - t46;
t14 = pkin(3) * t63 + t106;
t13 = -pkin(4) * t63 + t27;
t12 = pkin(4) * t66 + t26;
t8 = t107 + t150;
t7 = pkin(4) * t54 + t18;
t6 = t96 * t7;
t4 = t152 * t63 + t106;
t3 = t152 * t53 + t107;
t11 = [0, 0, 0, 0.2e1 * t95 * t120, t133 * t158, t135, -t136, 0, -pkin(6) * t135 + t114 * t95, pkin(6) * t136 + t114 * t97, t87 * t53 + t76 * t63 + (-t26 + (qJD(1) * t71 + t61) * t151) * qJD(2), t87 * t54 + t76 * t66 + (t151 * t156 - t27) * qJD(2), -t19 * t71 - t33 * t66 - t34 * t63 + t100, t146 + t19 * t38 - t33 * t26 + t34 * t27 + (t76 + t132) * t89, t16 * t71 + t30 * t66 + t31 * t63 + t100, t26 * qJD(2) - t14 * t61 - t24 * t63 - t32 * t53 - t71 * t8, t27 * qJD(2) - t104 * t14 - t24 * t66 - t32 * t54 - t72 * t8, t14 * t24 - t16 * t38 + t26 * t30 - t27 * t31 + t32 * t8 + t146, t109 * t41 + t141 * t22, (-t39 * t94 + t41 * t96) * t63 + (t145 - t23 * t94 + (-t39 * t96 - t41 * t94) * qJD(5)) * t71, t109 * t159 + t138 * t71 + t22 * t72 + t41 * t66, t71 * t47 - t23 * t72 - t39 * t66 + (-t129 * t94 + t63 * t96) * t159, t159 * t66 + t54 * t72, t1 * t66 + t13 * t39 + t29 * t23 + t6 * t72 + (-t159 * t4 - t3 * t72 - t147) * t94 + (t12 * t159 - t103) * t96 + ((-t17 * t96 - t28 * t94) * t159 - t2 * t72 + t15 * t141) * qJD(5), t13 * t41 - t2 * t66 + t29 * t22 + (-(qJD(5) * t28 + t4) * t159 - t147 - (qJD(5) * t10 + t3) * t72 + t15 * t129) * t96 + (-(-qJD(5) * t17 + t12) * t159 - (-qJD(5) * t9 + t7) * t72 + t103) * t94; 0, 0, 0, -t95 * t137, t133 * t99, 0, 0, 0, t99 * pkin(1) * t95, pkin(1) * t137, -t104 * t76 - t61 * t88 - t108, qJD(2) * t36 - t104 * t88 + t61 * t76 - t19, (t34 - t35) * t104 + (-t33 + t36) * t61 + (-t53 * t92 - t54 * t93) * pkin(2), t33 * t35 - t34 * t36 + (-t131 * t76 - t18 * t93 + t19 * t92) * pkin(2), -t84 * t53 + t86 * t54 + (-t31 - t35) * t104 + (t30 + t126) * t61, t25 * t61 + t108 + t144, -t24 * t61 + t25 * t104 + (0.2e1 * qJD(4) - t36) * qJD(2) + t19, t126 * t31 - t16 * t84 + t18 * t86 - t24 * t25 - t30 * t35, -t115 * t94 + t145, (-t23 - t115) * t96 + (-t22 + t160) * t94, t110 + t143, t102 - t142, t159 * t61, t1 * t61 + t105 * t94 + t125 * t39 + t154 * t96 + t84 * t23, t105 * t96 + t125 * t41 - t154 * t94 - t2 * t61 + t84 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t156 * qJD(2), -t79 + (-t61 + t122) * qJD(2), t157, t104 * t33 + t34 * t61 + t83, t157, -0.2e1 * t104 * qJD(2), -t54 + t56, t150 - t31 * t61 + (-qJD(4) - t30) * t104 + t119, 0, 0, 0, 0, 0, t102 + t142, -t110 + t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54 + t56, -t104 * t61, -t60 - t98, qJD(2) * t31 + t144 + t18, 0, 0, 0, 0, 0, -qJD(2) * t39 + t110, -qJD(2) * t41 + t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41 * t39, -t39 ^ 2 + t41 ^ 2, t22 + t160, t155 * t41 + t46, t54, -t15 * t41 + t155 * t2 - t94 * t3 + t6, t1 * t155 + t15 * t39 - t96 * t3 - t94 * t7;];
tauc_reg = t11;
