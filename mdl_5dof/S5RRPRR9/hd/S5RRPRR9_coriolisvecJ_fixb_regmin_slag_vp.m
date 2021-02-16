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
% tauc_reg [5x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% StartTime: 2021-01-15 21:47:33
% EndTime: 2021-01-15 21:47:46
% DurationCPUTime: 2.62s
% Computational Cost: add. (2962->307), mult. (7735->439), div. (0->0), fcn. (5805->8), ass. (0->168)
t142 = cos(qJ(2));
t197 = cos(pkin(9));
t165 = t197 * t142;
t125 = qJD(1) * t165;
t136 = sin(pkin(9));
t139 = sin(qJ(2));
t184 = qJD(1) * t139;
t101 = t136 * t184 - t125;
t219 = qJD(4) + qJD(5);
t228 = t101 + t219;
t137 = sin(qJ(5));
t138 = sin(qJ(4));
t140 = cos(qJ(5));
t141 = cos(qJ(4));
t117 = t137 * t141 + t138 * t140;
t208 = t228 * t117;
t179 = t141 * qJD(2);
t114 = t136 * t142 + t197 * t139;
t186 = qJD(1) * t114;
t82 = t138 * t186 - t179;
t84 = qJD(2) * t138 + t141 * t186;
t155 = t137 * t82 - t140 * t84;
t26 = t137 * t84 + t140 * t82;
t227 = t155 * t26;
t226 = t155 ^ 2 - t26 ^ 2;
t180 = qJD(5) * t140;
t181 = qJD(5) * t137;
t183 = qJD(4) * t138;
t178 = qJD(1) * qJD(2);
t171 = t139 * t178;
t124 = t136 * t171;
t95 = qJD(2) * t125 - t124;
t40 = qJD(4) * t179 + t141 * t95 - t183 * t186;
t41 = t84 * qJD(4) + t138 * t95;
t8 = -t137 * t41 + t140 * t40 - t82 * t180 - t84 * t181;
t97 = qJD(4) + t101;
t93 = qJD(5) + t97;
t225 = t26 * t93 + t8;
t132 = -pkin(2) * t142 - pkin(1);
t185 = qJD(1) * t132;
t120 = qJD(3) + t185;
t45 = pkin(3) * t101 - pkin(7) * t186 + t120;
t212 = -qJ(3) - pkin(6);
t122 = t212 * t139;
t118 = qJD(1) * t122;
t207 = qJD(2) * pkin(2);
t110 = t118 + t207;
t123 = t212 * t142;
t119 = qJD(1) * t123;
t166 = t197 * t119;
t69 = t136 * t110 - t166;
t64 = qJD(2) * pkin(7) + t69;
t21 = t138 * t45 + t141 * t64;
t16 = -pkin(8) * t82 + t21;
t13 = t16 * t181;
t107 = t136 * t119;
t68 = t197 * t110 + t107;
t63 = -qJD(2) * pkin(3) - t68;
t24 = t82 * pkin(4) + t63;
t224 = t24 * t26 + t13;
t145 = t155 * qJD(5) - t137 * t40 - t140 * t41;
t223 = -t155 * t93 + t145;
t126 = pkin(2) * t171;
t103 = t114 * qJD(2);
t94 = qJD(1) * t103;
t44 = pkin(3) * t94 - pkin(7) * t95 + t126;
t36 = t141 * t44;
t168 = qJD(2) * t212;
t100 = qJD(3) * t142 + t139 * t168;
t88 = t100 * qJD(1);
t148 = -qJD(3) * t139 + t142 * t168;
t89 = t148 * qJD(1);
t39 = t136 * t89 + t197 * t88;
t146 = -t21 * qJD(4) - t138 * t39 + t36;
t2 = pkin(4) * t94 - pkin(8) * t40 + t146;
t182 = qJD(4) * t141;
t153 = t138 * t44 + t141 * t39 + t45 * t182 - t64 * t183;
t5 = -pkin(8) * t41 + t153;
t173 = -t137 * t5 + t140 * t2;
t20 = -t138 * t64 + t141 * t45;
t15 = -pkin(8) * t84 + t20;
t10 = pkin(4) * t97 + t15;
t198 = t140 * t16;
t4 = t10 * t137 + t198;
t222 = -t4 * qJD(5) + t24 * t155 + t173;
t221 = -0.2e1 * t178;
t220 = t141 * t94 - t97 * t183;
t194 = t137 * t138;
t116 = -t140 * t141 + t194;
t209 = t228 * t116;
t218 = -t117 * t94 + t209 * t93;
t217 = t114 * t219;
t216 = pkin(2) * t139;
t215 = t82 * t97;
t214 = t84 * t97;
t128 = pkin(2) * t136 + pkin(7);
t213 = pkin(8) + t128;
t177 = pkin(2) * t184;
t56 = pkin(3) * t186 + pkin(7) * t101 + t177;
t72 = t197 * t118 + t107;
t211 = t138 * t56 + t141 * t72;
t149 = -t136 * t139 + t165;
t67 = -pkin(3) * t149 - pkin(7) * t114 + t132;
t81 = t136 * t122 - t197 * t123;
t74 = t141 * t81;
t210 = t138 * t67 + t74;
t206 = t186 * t26;
t205 = t186 * t155;
t204 = t186 * t82;
t203 = t186 * t84;
t201 = t138 * t40;
t200 = t138 * t94;
t196 = t114 * t138;
t195 = t114 * t141;
t193 = t138 * t101;
t106 = t149 * qJD(2);
t192 = t138 * t106;
t191 = t141 * t106;
t144 = qJD(1) ^ 2;
t190 = t142 * t144;
t143 = qJD(2) ^ 2;
t189 = t143 * t139;
t188 = t143 * t142;
t187 = t139 ^ 2 - t142 ^ 2;
t176 = t139 * t207;
t172 = t114 * t182;
t170 = qJD(5) * t10 + t5;
t38 = t136 * t88 - t197 * t89;
t167 = qJD(4) * t213;
t54 = t100 * t136 - t197 * t148;
t164 = pkin(1) * t221;
t71 = t118 * t136 - t166;
t80 = -t197 * t122 - t123 * t136;
t163 = t141 * t97;
t162 = 0.2e1 * t186;
t161 = -t116 * t94 - t208 * t93;
t160 = -t71 + (t183 + t193) * pkin(4);
t129 = -t197 * pkin(2) - pkin(3);
t159 = qJD(4) * t128 * t97 + t38;
t111 = t213 * t138;
t158 = pkin(8) * t193 + qJD(5) * t111 + t138 * t167 + t211;
t112 = t213 * t141;
t49 = t141 * t56;
t157 = pkin(4) * t186 + qJD(5) * t112 - t138 * t72 + t49 + (pkin(8) * t101 + t167) * t141;
t156 = t38 * t114 - t81 * t94;
t154 = -t97 * t193 + t220;
t55 = t197 * t100 + t136 * t148;
t57 = pkin(3) * t103 - pkin(7) * t106 + t176;
t152 = t138 * t57 + t141 * t55 + t67 * t182 - t81 * t183;
t151 = t172 + t192;
t150 = -t128 * t94 + t97 * t63;
t121 = -t141 * pkin(4) + t129;
t70 = t94 * t149;
t62 = t141 * t67;
t60 = t116 * t114;
t59 = t117 * t114;
t53 = pkin(4) * t196 + t80;
t50 = t141 * t57;
t23 = t151 * pkin(4) + t54;
t22 = -pkin(8) * t196 + t210;
t18 = pkin(4) * t41 + t38;
t17 = -pkin(4) * t149 - pkin(8) * t195 - t138 * t81 + t62;
t12 = t137 * t191 + (t219 * t195 + t192) * t140 - t194 * t217;
t11 = -t116 * t106 - t117 * t217;
t7 = -t151 * pkin(8) + t152;
t6 = -pkin(8) * t191 + pkin(4) * t103 - t138 * t55 + t50 + (-t74 + (pkin(8) * t114 - t67) * t138) * qJD(4);
t3 = t140 * t10 - t137 * t16;
t1 = [0, 0, 0, 0.2e1 * t142 * t171, t187 * t221, t188, -t189, 0, -pkin(6) * t188 + t139 * t164, pkin(6) * t189 + t142 * t164, t103 * t120 + t132 * t94 + (-t54 + (-qJD(1) * t149 + t101) * t216) * qJD(2), t106 * t120 + t132 * t95 + (t162 * t216 - t55) * qJD(2), -t101 * t55 - t103 * t69 - t106 * t68 + t149 * t39 + t186 * t54 + t80 * t95 + t156, t38 * t80 + t39 * t81 - t54 * t68 + t55 * t69 + (t120 + t185) * t176, t84 * t191 + (t141 * t40 - t183 * t84) * t114, (-t138 * t84 - t141 * t82) * t106 + (-t201 - t141 * t41 + (t138 * t82 - t141 * t84) * qJD(4)) * t114, t103 * t84 + t220 * t114 - t149 * t40 + t97 * t191, -t97 * t192 - t103 * t82 + t149 * t41 + (-t182 * t97 - t200) * t114, t103 * t97 - t70, (-t182 * t81 + t50) * t97 + t62 * t94 - (-t182 * t64 + t36) * t149 + t20 * t103 + t54 * t82 + t80 * t41 + t63 * t172 + ((-qJD(4) * t67 - t55) * t97 - (-qJD(4) * t45 - t39) * t149 + t63 * t106 + t156) * t138, -t152 * t97 - t210 * t94 + t153 * t149 - t21 * t103 + t54 * t84 + t80 * t40 + t63 * t191 + (t38 * t141 - t183 * t63) * t114, -t11 * t155 - t60 * t8, -t11 * t26 + t12 * t155 - t145 * t60 - t59 * t8, -t103 * t155 + t11 * t93 - t149 * t8 - t60 * t94, -t103 * t26 - t12 * t93 - t145 * t149 - t59 * t94, t103 * t93 - t70, (-t137 * t7 + t140 * t6) * t93 + (-t137 * t22 + t140 * t17) * t94 - t173 * t149 + t3 * t103 + t23 * t26 - t53 * t145 + t18 * t59 + t24 * t12 + ((-t137 * t17 - t140 * t22) * t93 + t4 * t149) * qJD(5), -t4 * t103 + t24 * t11 - t13 * t149 - t18 * t60 - t23 * t155 + t53 * t8 + (-(-qJD(5) * t22 + t6) * t93 - t17 * t94 + t2 * t149) * t137 + (-(qJD(5) * t17 + t7) * t93 - t22 * t94 + t170 * t149) * t140; 0, 0, 0, -t139 * t190, t187 * t144, 0, 0, 0, t144 * pkin(1) * t139, pkin(1) * t190, qJD(2) * t71 - t101 * t177 - t120 * t186 - t38, t72 * qJD(2) + t120 * t101 - t177 * t186 - t39, (t69 - t71) * t186 + (-t68 + t72) * t101 + (-t136 * t94 - t197 * t95) * pkin(2), t68 * t71 - t69 * t72 + (-t120 * t184 + t136 * t39 - t197 * t38) * pkin(2), t163 * t84 + t201, (t40 - t215) * t141 + (-t41 - t214) * t138, t163 * t97 + t200 - t203, t154 + t204, -t97 * t186, -t20 * t186 + t129 * t41 - t49 * t97 - t71 * t82 - t159 * t141 + (t72 * t97 + t150) * t138, t129 * t40 + t159 * t138 + t150 * t141 + t186 * t21 + t211 * t97 - t71 * t84, t117 * t8 + t155 * t209, -t116 * t8 + t117 * t145 + t155 * t208 + t209 * t26, t205 - t218, t161 + t206, -t93 * t186, (-t111 * t140 - t112 * t137) * t94 - t121 * t145 + t18 * t116 - t3 * t186 + (t137 * t158 - t140 * t157) * t93 + t160 * t26 + t208 * t24, -(-t111 * t137 + t112 * t140) * t94 + t121 * t8 + t18 * t117 + t4 * t186 + (t137 * t157 + t140 * t158) * t93 - t160 * t155 - t209 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t162 * qJD(2), -t124 + (t125 - t101) * qJD(2), -t101 ^ 2 - t186 ^ 2, t101 * t69 + t186 * t68 + t126, 0, 0, 0, 0, 0, t154 - t204, -t141 * t97 ^ 2 - t200 - t203, 0, 0, 0, 0, 0, t161 - t206, t205 + t218; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84 * t82, -t82 ^ 2 + t84 ^ 2, t40 + t215, t214 - t41, t94, t21 * t97 - t63 * t84 + t146, t20 * t97 + t63 * t82 - t153, -t227, t226, t225, t223, t94, -(-t137 * t15 - t198) * t93 + (t140 * t94 - t181 * t93 - t84 * t26) * pkin(4) + t222, (-t16 * t93 - t2) * t137 + (t15 * t93 - t170) * t140 + (-t137 * t94 + t155 * t84 - t180 * t93) * pkin(4) + t224; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t227, t226, t225, t223, t94, t4 * t93 + t222, -t137 * t2 - t140 * t170 + t3 * t93 + t224;];
tauc_reg = t1;
