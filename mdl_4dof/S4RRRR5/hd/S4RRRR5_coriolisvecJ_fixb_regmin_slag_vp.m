% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% tauc_reg [4x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRRR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:28:13
% EndTime: 2019-12-31 17:28:18
% DurationCPUTime: 1.46s
% Computational Cost: add. (1217->232), mult. (3192->362), div. (0->0), fcn. (2168->6), ass. (0->130)
t93 = sin(qJ(2));
t141 = qJD(1) * t93;
t92 = sin(qJ(3));
t125 = t92 * t141;
t95 = cos(qJ(3));
t133 = t95 * qJD(2);
t58 = t125 - t133;
t134 = t92 * qJD(2);
t60 = t95 * t141 + t134;
t91 = sin(qJ(4));
t94 = cos(qJ(4));
t108 = t91 * t58 - t94 * t60;
t19 = t94 * t58 + t91 * t60;
t175 = t108 * t19;
t174 = t108 ^ 2 - t19 ^ 2;
t135 = qJD(4) * t94;
t136 = qJD(4) * t91;
t131 = qJD(1) * qJD(2);
t96 = cos(qJ(2));
t121 = t96 * t131;
t130 = qJD(2) * qJD(3);
t33 = -qJD(3) * t125 + (t121 + t130) * t95;
t126 = t96 * t134;
t137 = qJD(3) * t95;
t127 = t93 * t137;
t103 = t126 + t127;
t34 = t103 * qJD(1) + t92 * t130;
t5 = -t58 * t135 - t60 * t136 + t94 * t33 - t91 * t34;
t132 = t96 * qJD(1);
t80 = -qJD(3) + t132;
t77 = -qJD(4) + t80;
t173 = -t19 * t77 + t5;
t68 = -t96 * pkin(2) - t93 * pkin(6) - pkin(1);
t52 = t68 * qJD(1);
t158 = t92 * t52;
t87 = pkin(5) * t132;
t73 = qJD(2) * pkin(6) + t87;
t26 = t95 * t73 + t158;
t15 = -t58 * pkin(7) + t26;
t13 = t15 * t136;
t72 = -qJD(2) * pkin(2) + pkin(5) * t141;
t35 = t58 * pkin(3) + t72;
t172 = t35 * t19 + t13;
t100 = t108 * qJD(4) - t91 * t33 - t94 * t34;
t171 = t108 * t77 + t100;
t82 = t93 * t131;
t114 = pkin(5) * t82;
t112 = pkin(2) * t93 - pkin(6) * t96;
t66 = t112 * qJD(2);
t53 = qJD(1) * t66;
t147 = -t92 * t114 - t95 * t53;
t102 = -t26 * qJD(3) - t147;
t4 = pkin(3) * t82 - t33 * pkin(7) + t102;
t138 = qJD(3) * t92;
t105 = t52 * t137 - t73 * t138 + t92 * t53;
t99 = -t95 * t114 + t105;
t7 = -t34 * pkin(7) + t99;
t123 = t94 * t4 - t91 * t7;
t25 = t95 * t52 - t92 * t73;
t14 = -t60 * pkin(7) + t25;
t12 = -t80 * pkin(3) + t14;
t154 = t94 * t15;
t2 = t91 * t12 + t154;
t170 = -t2 * qJD(4) + t35 * t108 + t123;
t169 = -0.2e1 * t131;
t62 = t91 * t95 + t94 * t92;
t43 = t62 * t93;
t129 = t93 * t138;
t168 = t96 * t133 - t129;
t167 = qJD(3) + qJD(4);
t166 = pkin(6) + pkin(7);
t165 = pkin(3) * t92;
t164 = t33 * t92;
t163 = t58 * t80;
t162 = t60 * t80;
t161 = t72 * t92;
t160 = t72 * t95;
t159 = t80 * t95;
t157 = t92 * t93;
t156 = t92 * t96;
t155 = t93 * t95;
t153 = t95 * t96;
t98 = qJD(1) ^ 2;
t152 = t96 * t98;
t97 = qJD(2) ^ 2;
t151 = t97 * t93;
t150 = t97 * t96;
t61 = t91 * t92 - t94 * t95;
t104 = t61 * t96;
t149 = qJD(1) * t104 - t167 * t61;
t148 = (-t132 + t167) * t62;
t146 = t68 * t137 + t92 * t66;
t145 = t93 * pkin(5) * t134 + t95 * t66;
t63 = t112 * qJD(1);
t144 = pkin(5) * t125 + t95 * t63;
t81 = pkin(5) * t153;
t143 = t92 * t68 + t81;
t89 = t93 ^ 2;
t142 = -t96 ^ 2 + t89;
t140 = qJD(2) * t93;
t139 = qJD(2) * t96;
t128 = t96 * t138;
t122 = qJD(3) * t166;
t118 = qJD(4) * t12 + t7;
t117 = t58 + t133;
t116 = -t60 + t134;
t115 = pkin(1) * t169;
t113 = pkin(3) * t138 - t132 * t165 - t87;
t106 = pkin(3) * t93 - pkin(7) * t153;
t75 = t166 * t95;
t111 = t106 * qJD(1) + qJD(4) * t75 + t95 * t122 + t144;
t48 = t92 * t63;
t74 = t166 * t92;
t110 = -qJD(4) * t74 - t48 - (-pkin(5) * t155 - pkin(7) * t156) * qJD(1) - t92 * t122;
t57 = t95 * t68;
t24 = -pkin(7) * t155 + t57 + (-pkin(5) * t92 - pkin(3)) * t96;
t29 = -pkin(7) * t157 + t143;
t109 = t91 * t24 + t94 * t29;
t107 = qJD(1) * t89 - t80 * t96;
t85 = -t95 * pkin(3) - pkin(2);
t67 = (pkin(5) + t165) * t93;
t44 = t61 * t93;
t36 = t103 * pkin(3) + pkin(5) * t139;
t17 = t34 * pkin(3) + pkin(5) * t121;
t11 = -t136 * t157 + (t167 * t155 + t126) * t94 + t168 * t91;
t10 = -qJD(2) * t104 - t167 * t43;
t9 = -t103 * pkin(7) + (-t93 * t133 - t128) * pkin(5) + t146;
t8 = t106 * qJD(2) + (-t81 + (pkin(7) * t93 - t68) * t92) * qJD(3) + t145;
t1 = t94 * t12 - t91 * t15;
t3 = [0, 0, 0, 0.2e1 * t96 * t82, t142 * t169, t150, -t151, 0, -pkin(5) * t150 + t93 * t115, pkin(5) * t151 + t96 * t115, t33 * t155 + t168 * t60, (-t58 * t95 - t60 * t92) * t139 + (-t164 - t34 * t95 + (t58 * t92 - t60 * t95) * qJD(3)) * t93, t80 * t129 - t33 * t96 + (t107 * t95 + t60 * t93) * qJD(2), t80 * t127 + t34 * t96 + (-t107 * t92 - t58 * t93) * qJD(2), (-t80 - t132) * t140, -(-t68 * t138 + t145) * t80 + (t72 * t137 + pkin(5) * t34 + (qJD(1) * t57 + t25) * qJD(2)) * t93 + ((pkin(5) * t58 + t161) * qJD(2) + (t158 + (pkin(5) * t80 + t73) * t95) * qJD(3) + t147) * t96, (-pkin(5) * t128 + t146) * t80 + t105 * t96 + (pkin(5) * t33 - t72 * t138) * t93 + ((pkin(5) * t60 + t160) * t96 + (-pkin(5) * t159 - t143 * qJD(1) - t26) * t93) * qJD(2), -t10 * t108 - t5 * t44, -t10 * t19 - t100 * t44 + t108 * t11 - t5 * t43, -t10 * t77 - t5 * t96 + (-qJD(1) * t44 - t108) * t140, t11 * t77 - t100 * t96 + (-qJD(1) * t43 - t19) * t140, (-t77 - t132) * t140, -(t94 * t8 - t91 * t9) * t77 - t123 * t96 + t36 * t19 - t67 * t100 + t17 * t43 + t35 * t11 + (t109 * t77 + t2 * t96) * qJD(4) + ((t94 * t24 - t91 * t29) * qJD(1) + t1) * t140, t35 * t10 - t13 * t96 - t17 * t44 - t36 * t108 + t67 * t5 + ((-qJD(4) * t29 + t8) * t77 + t4 * t96) * t91 + ((qJD(4) * t24 + t9) * t77 + t118 * t96) * t94 + (-qJD(1) * t109 - t2) * t140; 0, 0, 0, -t93 * t152, t142 * t98, 0, 0, 0, t98 * pkin(1) * t93, pkin(1) * t152, -t60 * t159 + t164, (t33 + t163) * t95 + (-t34 + t162) * t92, -t80 * t137 + (t116 * t93 + t80 * t153) * qJD(1), t80 * t138 + (t117 * t93 - t80 * t156) * qJD(1), t80 * t141, -pkin(2) * t34 + t144 * t80 + (pkin(6) * t159 + t161) * qJD(3) + ((-pkin(6) * t134 - t25) * t93 + (-pkin(5) * t117 - t161) * t96) * qJD(1), -pkin(2) * t33 - t48 * t80 + (-t92 * pkin(6) * t80 + t160) * qJD(3) + (-t72 * t153 + (-pkin(6) * t133 + t26) * t93 + (t116 * t96 + t80 * t155) * pkin(5)) * qJD(1), -t108 * t149 + t5 * t62, t100 * t62 + t108 * t148 - t149 * t19 - t5 * t61, -t149 * t77 + (qJD(2) * t62 + t108) * t141, t148 * t77 + (-qJD(2) * t61 + t19) * t141, t77 * t141, t17 * t61 - t85 * t100 + (t110 * t91 + t111 * t94) * t77 + t148 * t35 + t113 * t19 + ((-t94 * t74 - t91 * t75) * qJD(2) - t1) * t141, t17 * t62 + t85 * t5 + (t110 * t94 - t111 * t91) * t77 + t149 * t35 - t113 * t108 + (-(-t91 * t74 + t94 * t75) * qJD(2) + t2) * t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60 * t58, -t58 ^ 2 + t60 ^ 2, t33 - t163, -t162 - t34, t82, -t26 * t80 - t72 * t60 + t102, -t25 * t80 + t72 * t58 - t99, -t175, t174, t173, t171, t82, (-t91 * t14 - t154) * t77 + (t77 * t136 - t60 * t19 + t94 * t82) * pkin(3) + t170, (t15 * t77 - t4) * t91 + (-t14 * t77 - t118) * t94 + (t108 * t60 + t77 * t135 - t91 * t82) * pkin(3) + t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t175, t174, t173, t171, t82, -t2 * t77 + t170, -t1 * t77 - t118 * t94 - t91 * t4 + t172;];
tauc_reg = t3;
