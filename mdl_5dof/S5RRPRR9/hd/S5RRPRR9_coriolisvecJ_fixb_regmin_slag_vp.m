% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauc_reg [5x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRR9_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR9_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:21:58
% EndTime: 2019-12-31 20:22:06
% DurationCPUTime: 2.41s
% Computational Cost: add. (2894->292), mult. (7541->421), div. (0->0), fcn. (5685->8), ass. (0->164)
t134 = sin(pkin(9));
t138 = sin(qJ(2));
t181 = qJD(1) * t138;
t135 = cos(pkin(9));
t141 = cos(qJ(2));
t190 = t135 * t141;
t101 = qJD(1) * t190 - t134 * t181;
t214 = qJD(4) + qJD(5);
t223 = t101 - t214;
t136 = sin(qJ(5));
t137 = sin(qJ(4));
t139 = cos(qJ(5));
t140 = cos(qJ(4));
t116 = t136 * t140 + t137 * t139;
t204 = t223 * t116;
t113 = t134 * t141 + t135 * t138;
t103 = t113 * qJD(1);
t176 = t140 * qJD(2);
t82 = t103 * t137 - t176;
t84 = qJD(2) * t137 + t103 * t140;
t153 = t136 * t82 - t139 * t84;
t26 = t136 * t84 + t139 * t82;
t222 = t153 * t26;
t221 = t153 ^ 2 - t26 ^ 2;
t177 = qJD(5) * t139;
t178 = qJD(5) * t136;
t180 = qJD(4) * t137;
t175 = qJD(1) * qJD(2);
t167 = t141 * t175;
t168 = t138 * t175;
t94 = -t134 * t168 + t135 * t167;
t40 = qJD(4) * t176 - t103 * t180 + t140 * t94;
t41 = qJD(4) * t84 + t137 * t94;
t8 = -t136 * t41 + t139 * t40 - t82 * t177 - t84 * t178;
t96 = qJD(4) - t101;
t92 = qJD(5) + t96;
t220 = t26 * t92 + t8;
t171 = -pkin(2) * t141 - pkin(1);
t155 = t171 * qJD(1);
t119 = qJD(3) + t155;
t45 = -pkin(3) * t101 - pkin(7) * t103 + t119;
t208 = -qJ(3) - pkin(6);
t121 = t208 * t138;
t117 = qJD(1) * t121;
t203 = qJD(2) * pkin(2);
t109 = t117 + t203;
t122 = t208 * t141;
t118 = qJD(1) * t122;
t191 = t135 * t118;
t69 = t134 * t109 - t191;
t64 = qJD(2) * pkin(7) + t69;
t21 = t137 * t45 + t140 * t64;
t16 = -pkin(8) * t82 + t21;
t13 = t16 * t178;
t106 = t134 * t118;
t68 = t109 * t135 + t106;
t63 = -qJD(2) * pkin(3) - t68;
t24 = pkin(4) * t82 + t63;
t219 = t24 * t26 + t13;
t144 = qJD(5) * t153 - t136 * t40 - t139 * t41;
t218 = -t153 * t92 + t144;
t125 = pkin(2) * t168;
t102 = t113 * qJD(2);
t93 = qJD(1) * t102;
t44 = pkin(3) * t93 - pkin(7) * t94 + t125;
t36 = t140 * t44;
t164 = qJD(2) * t208;
t100 = -t138 * qJD(3) + t141 * t164;
t147 = t100 * qJD(1);
t99 = qJD(3) * t141 + t138 * t164;
t88 = t99 * qJD(1);
t39 = t134 * t147 + t135 * t88;
t145 = -qJD(4) * t21 - t137 * t39 + t36;
t2 = pkin(4) * t93 - pkin(8) * t40 + t145;
t179 = qJD(4) * t140;
t151 = t137 * t44 + t140 * t39 + t45 * t179 - t64 * t180;
t5 = -pkin(8) * t41 + t151;
t170 = -t136 * t5 + t139 * t2;
t20 = -t137 * t64 + t140 * t45;
t15 = -pkin(8) * t84 + t20;
t10 = pkin(4) * t96 + t15;
t194 = t139 * t16;
t4 = t136 * t10 + t194;
t217 = -qJD(5) * t4 + t24 * t153 + t170;
t216 = -0.2e1 * t175;
t215 = t140 * t93 - t96 * t180;
t189 = t136 * t137;
t115 = -t139 * t140 + t189;
t205 = t223 * t115;
t213 = -t116 * t93 - t205 * t92;
t212 = t214 * t113;
t211 = t82 * t96;
t210 = t84 * t96;
t127 = pkin(2) * t134 + pkin(7);
t209 = pkin(8) + t127;
t56 = pkin(2) * t181 + pkin(3) * t103 - pkin(7) * t101;
t72 = t117 * t135 + t106;
t207 = t137 * t56 + t140 * t72;
t112 = t134 * t138 - t190;
t67 = pkin(3) * t112 - pkin(7) * t113 + t171;
t81 = t121 * t134 - t122 * t135;
t74 = t140 * t81;
t206 = t137 * t67 + t74;
t202 = t103 * t26;
t201 = t103 * t153;
t200 = t103 * t82;
t199 = t103 * t84;
t197 = t137 * t40;
t196 = t137 * t93;
t193 = t113 * t137;
t192 = t113 * t140;
t188 = t137 * t101;
t105 = t112 * qJD(2);
t187 = t137 * t105;
t186 = t140 * t105;
t143 = qJD(1) ^ 2;
t185 = t141 * t143;
t142 = qJD(2) ^ 2;
t184 = t142 * t138;
t183 = t142 * t141;
t182 = t138 ^ 2 - t141 ^ 2;
t174 = t138 * t203;
t128 = -pkin(2) * t135 - pkin(3);
t169 = t113 * t179;
t166 = qJD(5) * t10 + t5;
t38 = t134 * t88 - t135 * t147;
t54 = -t135 * t100 + t134 * t99;
t163 = qJD(4) * t209;
t162 = pkin(1) * t216;
t71 = t117 * t134 - t191;
t80 = -t135 * t121 - t122 * t134;
t161 = t140 * t96;
t160 = -t115 * t93 + t204 * t92;
t159 = -t71 + (t180 - t188) * pkin(4);
t158 = qJD(4) * t127 * t96 + t38;
t110 = t209 * t137;
t157 = -pkin(8) * t188 + qJD(5) * t110 + t137 * t163 + t207;
t111 = t209 * t140;
t49 = t140 * t56;
t156 = pkin(4) * t103 + qJD(5) * t111 - t137 * t72 + t49 + (-pkin(8) * t101 + t163) * t140;
t154 = t38 * t113 - t81 * t93;
t152 = t96 * t188 + t215;
t55 = t100 * t134 + t135 * t99;
t57 = pkin(3) * t102 + pkin(7) * t105 + t174;
t150 = t137 * t57 + t140 * t55 + t67 * t179 - t81 * t180;
t149 = t169 - t187;
t148 = -t127 * t93 + t96 * t63;
t120 = -pkin(4) * t140 + t128;
t70 = t93 * t112;
t62 = t140 * t67;
t60 = t115 * t113;
t59 = t116 * t113;
t53 = pkin(4) * t193 + t80;
t50 = t140 * t57;
t23 = pkin(4) * t149 + t54;
t22 = -pkin(8) * t193 + t206;
t18 = pkin(4) * t41 + t38;
t17 = pkin(4) * t112 - pkin(8) * t192 - t137 * t81 + t62;
t12 = -t136 * t186 + (t214 * t192 - t187) * t139 - t189 * t212;
t11 = t115 * t105 - t116 * t212;
t7 = -pkin(8) * t149 + t150;
t6 = pkin(8) * t186 + pkin(4) * t102 - t137 * t55 + t50 + (-t74 + (pkin(8) * t113 - t67) * t137) * qJD(4);
t3 = t139 * t10 - t136 * t16;
t1 = [0, 0, 0, 0.2e1 * t138 * t167, t182 * t216, t183, -t184, 0, -pkin(6) * t183 + t138 * t162, pkin(6) * t184 + t141 * t162, t101 * t55 - t102 * t69 + t103 * t54 + t105 * t68 - t39 * t112 + t80 * t94 + t154, t38 * t80 + t39 * t81 - t68 * t54 + t69 * t55 + (t119 + t155) * t174, -t84 * t186 + (t140 * t40 - t84 * t180) * t113, -(-t137 * t84 - t140 * t82) * t105 + (-t197 - t140 * t41 + (t137 * t82 - t140 * t84) * qJD(4)) * t113, t102 * t84 + t112 * t40 + t215 * t113 - t96 * t186, t96 * t187 - t102 * t82 - t112 * t41 + (-t96 * t179 - t196) * t113, t102 * t96 + t70, (-t81 * t179 + t50) * t96 + t62 * t93 + (-t64 * t179 + t36) * t112 + t20 * t102 + t54 * t82 + t80 * t41 + t63 * t169 + ((-qJD(4) * t67 - t55) * t96 + (-qJD(4) * t45 - t39) * t112 - t63 * t105 + t154) * t137, -t150 * t96 - t206 * t93 - t151 * t112 - t21 * t102 + t54 * t84 + t80 * t40 - t63 * t186 + (t38 * t140 - t63 * t180) * t113, -t11 * t153 - t60 * t8, -t11 * t26 + t12 * t153 - t144 * t60 - t59 * t8, -t102 * t153 + t11 * t92 + t112 * t8 - t60 * t93, -t102 * t26 + t112 * t144 - t12 * t92 - t59 * t93, t102 * t92 + t70, (-t136 * t7 + t139 * t6) * t92 + (-t136 * t22 + t139 * t17) * t93 + t170 * t112 + t3 * t102 + t23 * t26 - t53 * t144 + t18 * t59 + t24 * t12 + ((-t136 * t17 - t139 * t22) * t92 - t4 * t112) * qJD(5), -t4 * t102 + t24 * t11 + t13 * t112 - t18 * t60 - t23 * t153 + t53 * t8 + (-(-qJD(5) * t22 + t6) * t92 - t17 * t93 - t2 * t112) * t136 + (-(qJD(5) * t17 + t7) * t92 - t22 * t93 - t166 * t112) * t139; 0, 0, 0, -t138 * t185, t182 * t143, 0, 0, 0, t143 * pkin(1) * t138, pkin(1) * t185, (t69 - t71) * t103 + (t68 - t72) * t101 + (-t134 * t93 - t135 * t94) * pkin(2), t68 * t71 - t69 * t72 + (-t119 * t181 + t134 * t39 - t135 * t38) * pkin(2), t161 * t84 + t197, (t40 - t211) * t140 + (-t41 - t210) * t137, t161 * t96 + t196 - t199, t152 + t200, -t96 * t103, -t20 * t103 + t128 * t41 - t49 * t96 - t71 * t82 - t158 * t140 + (t72 * t96 + t148) * t137, t21 * t103 + t128 * t40 + t158 * t137 + t148 * t140 + t207 * t96 - t71 * t84, t116 * t8 - t153 * t205, -t115 * t8 + t116 * t144 - t153 * t204 - t205 * t26, t201 - t213, t160 + t202, -t92 * t103, (-t110 * t139 - t111 * t136) * t93 - t120 * t144 + t18 * t115 - t3 * t103 + (t136 * t157 - t139 * t156) * t92 + t159 * t26 - t204 * t24, -(-t110 * t136 + t111 * t139) * t93 + t120 * t8 + t18 * t116 + t4 * t103 + (t136 * t156 + t139 * t157) * t92 - t159 * t153 + t205 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101 ^ 2 - t103 ^ 2, -t101 * t69 + t103 * t68 + t125, 0, 0, 0, 0, 0, t152 - t200, -t96 ^ 2 * t140 - t196 - t199, 0, 0, 0, 0, 0, t160 - t202, t201 + t213; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84 * t82, -t82 ^ 2 + t84 ^ 2, t40 + t211, t210 - t41, t93, t21 * t96 - t63 * t84 + t145, t20 * t96 + t63 * t82 - t151, -t222, t221, t220, t218, t93, -(-t136 * t15 - t194) * t92 + (t139 * t93 - t92 * t178 - t84 * t26) * pkin(4) + t217, (-t16 * t92 - t2) * t136 + (t15 * t92 - t166) * t139 + (-t136 * t93 + t153 * t84 - t92 * t177) * pkin(4) + t219; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t222, t221, t220, t218, t93, t4 * t92 + t217, -t136 * t2 - t166 * t139 + t3 * t92 + t219;];
tauc_reg = t1;
