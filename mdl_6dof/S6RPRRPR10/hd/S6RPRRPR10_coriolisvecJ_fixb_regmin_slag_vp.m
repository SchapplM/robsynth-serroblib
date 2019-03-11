% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
% 
% Output:
% tauc_reg [6x31]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRPR10_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR10_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:37:10
% EndTime: 2019-03-09 05:37:20
% DurationCPUTime: 3.57s
% Computational Cost: add. (2894->384), mult. (6116->535), div. (0->0), fcn. (3769->6), ass. (0->205)
t134 = sin(qJ(3));
t221 = qJD(1) * t134;
t122 = qJD(4) + t221;
t139 = -pkin(1) - pkin(7);
t133 = sin(qJ(4));
t136 = cos(qJ(4));
t239 = qJ(5) * t136;
t259 = pkin(4) + pkin(5);
t156 = t133 * t259 - t239;
t276 = t139 - t156;
t214 = qJD(4) * t134;
t137 = cos(qJ(3));
t177 = pkin(3) * t137 + pkin(8) * t134;
t92 = qJD(3) * t177 + qJD(2);
t275 = -t139 * t214 + t92;
t265 = qJD(4) - qJD(6);
t132 = sin(qJ(6));
t135 = cos(qJ(6));
t207 = t136 * qJD(3);
t220 = qJD(1) * t137;
t93 = t133 * t220 - t207;
t196 = t136 * t220;
t218 = qJD(3) * t133;
t95 = t196 + t218;
t169 = t132 * t95 - t135 * t93;
t46 = t132 * t93 + t135 * t95;
t274 = -t169 ^ 2 + t46 ^ 2;
t217 = qJD(3) * t134;
t190 = t133 * t217;
t111 = qJD(1) * t190;
t237 = qJD(4) * t95;
t58 = -t111 + t237;
t110 = t122 * qJD(5);
t205 = qJD(1) * qJD(3);
t123 = t137 * t205;
t119 = qJ(5) * t123;
t213 = qJD(4) * t136;
t215 = qJD(4) * t133;
t72 = t92 * qJD(1);
t103 = pkin(3) * t134 - pkin(8) * t137 + qJ(2);
t82 = t103 * qJD(1);
t118 = qJD(1) * t139 + qJD(2);
t102 = t134 * t118;
t84 = qJD(3) * pkin(8) + t102;
t229 = t136 * t137;
t90 = t118 * t229;
t154 = -qJD(3) * t90 - t133 * t72 - t213 * t82 + t215 * t84;
t9 = t110 + t119 - t154;
t6 = pkin(9) * t58 + t9;
t216 = qJD(3) * t137;
t195 = t133 * t216;
t187 = t118 * t195 - t136 * t72 + t213 * t84 + t215 * t82;
t212 = qJD(4) * t137;
t193 = t133 * t212;
t57 = qJD(1) * (t134 * t207 + t193) - qJD(4) * t207;
t7 = pkin(9) * t57 - t123 * t259 + t187;
t197 = -t132 * t6 + t135 * t7;
t236 = t118 * t137;
t85 = -qJD(3) * pkin(3) - t236;
t148 = qJ(5) * t95 - t85;
t23 = -t259 * t93 + t148;
t273 = t23 * t46 - t197;
t179 = qJD(1) + t214;
t235 = t122 * t133;
t242 = t137 * t57;
t272 = ((-t95 + t196) * t134 + t122 * t229) * qJD(3) - t179 * t235 - t242;
t203 = 0.2e1 * qJD(1);
t270 = t23 * t169;
t269 = t46 * t169;
t178 = pkin(4) * t123;
t12 = -t178 + t187;
t113 = t122 * qJ(5);
t37 = t133 * t82 + t136 * t84;
t30 = t113 + t37;
t268 = -t122 * t30 + t12;
t36 = -t133 * t84 + t136 * t82;
t224 = qJD(5) - t36;
t267 = qJD(5) * t133 + t102;
t206 = qJD(6) - t122;
t266 = -qJD(6) + t206;
t208 = qJD(6) * t135;
t209 = qJD(6) * t132;
t10 = t132 * t58 - t135 * t57 + t208 * t93 - t209 * t95;
t264 = t169 * t206 + t10;
t240 = qJ(5) * t133;
t262 = -t136 * t259 - t240;
t260 = t95 ^ 2;
t258 = pkin(8) - pkin(9);
t33 = pkin(4) * t93 - t148;
t257 = t33 * t95;
t256 = t95 * t93;
t97 = t132 * t133 + t135 * t136;
t151 = qJD(1) * t97;
t47 = t265 * t97;
t255 = t134 * t151 + t47;
t164 = t132 * t136 - t133 * t135;
t153 = t164 * t134;
t254 = qJD(1) * t153 + t132 * t213 + t133 * t208 - t135 * t215 - t136 * t209;
t172 = pkin(4) * t133 - t239;
t253 = -t122 * t172 + t267;
t252 = -t122 * t156 + t267;
t99 = t177 * qJD(1);
t251 = t133 * t99 + t90;
t250 = qJ(5) * t93;
t248 = t122 * t93;
t247 = t122 * t95;
t145 = -qJ(5) * t57 + qJD(5) * t95 - t118 * t217;
t13 = pkin(4) * t58 - t145;
t246 = t13 * t133;
t245 = t13 * t136;
t26 = pkin(9) * t93 + t37;
t20 = t113 + t26;
t244 = t132 * t20;
t243 = t133 * t57;
t230 = t134 * t139;
t241 = t103 * t133 + t136 * t230;
t74 = t97 * t137;
t238 = qJD(1) * t74;
t234 = t122 * t134;
t233 = t122 * t136;
t232 = t133 * t137;
t231 = t134 * t136;
t140 = qJD(3) ^ 2;
t228 = t140 * t134;
t227 = t140 * t137;
t141 = qJD(1) ^ 2;
t226 = t141 * qJ(2);
t225 = -pkin(9) * t95 + t224;
t131 = t137 ^ 2;
t223 = t134 ^ 2 - t131;
t222 = -t140 - t141;
t219 = qJD(3) * t206;
t210 = qJD(5) * t136;
t17 = -t122 * t259 + t225;
t204 = t132 * t7 + t135 * t6 + t17 * t208;
t202 = pkin(8) * t235;
t201 = pkin(8) * t233;
t39 = qJ(5) * t220 + t251;
t200 = pkin(8) * t216;
t199 = qJD(2) * t203;
t55 = qJ(5) * t134 + t241;
t109 = t258 * t136;
t194 = t139 * t216;
t191 = t136 * t212;
t189 = -t132 * t57 - t135 * t58;
t188 = -t118 * t232 + t136 * t99;
t186 = qJD(4) * t93 - t57;
t185 = -t58 + t237;
t184 = -t85 + t236;
t183 = -t95 + t218;
t182 = t93 + t207;
t115 = t133 * t230;
t181 = t103 * t136 - t115;
t180 = t206 ^ 2;
t108 = t258 * t133;
t175 = pkin(9) * t133 * t221 + qJD(6) * t108 - t215 * t258 - t39;
t150 = pkin(9) * t231 - t137 * t259;
t174 = qJD(1) * t150 - t109 * t265 - t188;
t173 = pkin(4) * t136 + t240;
t2 = t132 * t17 + t135 * t20;
t35 = t115 + (-pkin(9) * t137 - t103) * t136 - t259 * t134;
t41 = pkin(9) * t232 + t55;
t171 = -t132 * t41 + t135 * t35;
t170 = t132 * t35 + t135 * t41;
t29 = -pkin(4) * t122 + t224;
t168 = t133 * t30 - t136 * t29;
t167 = t133 * t29 + t136 * t30;
t166 = qJ(5) * t135 - t132 * t259;
t165 = -qJ(5) * t132 - t135 * t259;
t163 = qJD(1) * t131 - t234;
t162 = -t103 * t215 - t133 * t194 + t136 * t275;
t161 = -t134 * t33 + t200;
t160 = t134 * t85 - t200;
t159 = -t139 + t172;
t158 = -t20 * t209 + t204;
t157 = t122 * t37 - t187;
t152 = qJD(1) * t164;
t147 = t103 * t213 + t133 * t275 + t136 * t194;
t146 = t122 * t36 + t154;
t19 = qJ(5) * t216 + qJD(5) * t134 + t147;
t11 = qJD(6) * t46 + t189;
t144 = -qJD(4) * t168 + t12 * t133 + t136 * t9;
t143 = t265 * t164;
t142 = t93 * t217 + (-t58 - t111) * t137 + (-t136 * t179 - t195) * t122;
t112 = t134 * t123;
t104 = -pkin(3) - t173;
t88 = pkin(3) - t262;
t73 = t132 * t229 - t135 * t232;
t64 = t159 * t137;
t56 = -pkin(4) * t134 - t181;
t54 = t276 * t137;
t49 = pkin(4) * t95 + t250;
t40 = -pkin(4) * t220 - t188;
t34 = -t259 * t95 - t250;
t31 = t248 - t57;
t28 = (qJD(4) * t173 - t210) * t137 - t159 * t217;
t24 = -pkin(4) * t216 - t162;
t22 = t137 * t143 - t217 * t97;
t21 = qJD(3) * t153 + t137 * t47;
t18 = (qJD(4) * t262 + t210) * t137 - t276 * t217;
t15 = (-t190 + t191) * pkin(9) + t19;
t14 = pkin(9) * t193 + qJD(3) * t150 - t162;
t8 = -t259 * t58 + t145;
t1 = t135 * t17 - t244;
t3 = [0, 0, 0, 0, t199, qJ(2) * t199, -0.2e1 * t112, 0.2e1 * t223 * t205, -t228, -t227, 0, -t139 * t228 + (qJ(2) * t216 + qJD(2) * t134) * t203, -t139 * t227 + (-qJ(2) * t217 + qJD(2) * t137) * t203, -t95 * t193 + (-t217 * t95 - t242) * t136 (t133 * t95 + t136 * t93) * t217 + (t243 - t136 * t58 + (t133 * t93 - t136 * t95) * qJD(4)) * t137, -t122 * t193 - t134 * t57 + (t136 * t163 + t137 * t95) * qJD(3), -t122 * t191 - t134 * t58 + (-t133 * t163 - t137 * t93) * qJD(3), t122 * t216 + t112, t162 * t122 - t187 * t134 + (-t139 * t58 + t85 * t213) * t137 + ((qJD(1) * t181 + t36) * t137 + (t133 * t184 + t139 * t93) * t134) * qJD(3), -t147 * t122 + t154 * t134 + (t139 * t57 - t85 * t215) * t137 + ((-qJD(1) * t241 - t37) * t137 + (t136 * t184 + t139 * t95) * t134) * qJD(3), -t122 * t24 + t28 * t93 + t58 * t64 + (-t218 * t33 - t12) * t134 + (t33 * t213 + t246 + (-qJD(1) * t56 - t29) * qJD(3)) * t137, -t19 * t93 + t24 * t95 - t55 * t58 - t56 * t57 + t168 * t217 + (-qJD(4) * t167 + t12 * t136 - t133 * t9) * t137, t122 * t19 - t28 * t95 + t57 * t64 + (t207 * t33 + t9) * t134 + (t33 * t215 - t245 + (qJD(1) * t55 + t30) * qJD(3)) * t137, t12 * t56 + t13 * t64 + t19 * t30 + t24 * t29 + t28 * t33 + t55 * t9, t10 * t74 + t22 * t46, -t10 * t73 - t11 * t74 - t169 * t22 + t21 * t46, -t10 * t134 + t206 * t22 + (-t46 - t238) * t216, t11 * t134 + t206 * t21 + (qJD(1) * t73 + t169) * t216, -t206 * t216 + t112 (-t132 * t15 + t135 * t14) * t206 - t197 * t134 + t18 * t169 + t54 * t11 + t8 * t73 - t23 * t21 + (t134 * t2 - t170 * t206) * qJD(6) + (-qJD(1) * t171 - t1) * t216 -(qJD(6) * t171 + t132 * t14 + t135 * t15) * t206 + t158 * t134 + t18 * t46 + t54 * t10 + t8 * t74 + t23 * t22 + (qJD(1) * t170 + t2) * t216; 0, 0, 0, 0, -t141, -t226, 0, 0, 0, 0, 0, t222 * t134, t222 * t137, 0, 0, 0, 0, 0, t142, -t272, t142 (qJD(1) * t95 + t134 * t185 - t216 * t93) * t136 + (qJD(1) * t93 + t134 * t186 + t216 * t95) * t133, t272, -t168 * qJD(1) + (qJD(3) * t167 - t13) * t137 + (qJD(3) * t33 + t144) * t134, 0, 0, 0, 0, 0, t206 * t151 + (-t164 * t219 + t11) * t137 + (t206 * t47 + (t137 * t152 - t169) * qJD(3)) * t134, -t206 * t152 + (-t219 * t97 + t10) * t137 + (-t143 * t206 + (-t46 + t238) * qJD(3)) * t134; 0, 0, 0, 0, 0, 0, t137 * t141 * t134, -t223 * t141, 0, 0, 0, -t137 * t226, t134 * t226, t233 * t95 - t243 (-t57 - t248) * t136 + (-t58 - t247) * t133, t122 * t213 + (t122 * t231 + t137 * t183) * qJD(1), -t122 * t215 + (-t133 * t234 + t137 * t182) * qJD(1), -t122 * t220, -pkin(3) * t58 - t188 * t122 - t182 * t102 + (t133 * t85 - t201) * qJD(4) + (t133 * t160 - t36 * t137) * qJD(1), pkin(3) * t57 + t251 * t122 + t183 * t102 + (t136 * t85 + t202) * qJD(4) + (t136 * t160 + t37 * t137) * qJD(1), t104 * t58 + t122 * t40 - t245 - t253 * t93 + (t133 * t33 - t201) * qJD(4) + (-t133 * t161 + t137 * t29) * qJD(1), t39 * t93 - t40 * t95 + (pkin(8) * t185 + t122 * t29 + t9) * t136 + (pkin(8) * t186 + t268) * t133, t104 * t57 - t122 * t39 - t246 + t253 * t95 + (-t136 * t33 - t202) * qJD(4) + (t136 * t161 - t137 * t30) * qJD(1), pkin(8) * t144 + t104 * t13 - t253 * t33 - t29 * t40 - t30 * t39, -t10 * t164 + t255 * t46, -t10 * t97 + t11 * t164 - t169 * t255 - t254 * t46, t255 * t206 + (qJD(3) * t164 + t46) * t220, -t254 * t206 + (qJD(3) * t97 - t169) * t220, t206 * t220, t88 * t11 + t8 * t97 + t252 * t169 + t254 * t23 - (t132 * t175 + t135 * t174) * t206 + (-(t108 * t135 - t109 * t132) * qJD(3) + t1) * t220, t88 * t10 - t8 * t164 + t252 * t46 + t255 * t23 - (-t132 * t174 + t135 * t175) * t206 + ((t108 * t132 + t109 * t135) * qJD(3) - t2) * t220; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t256, -t93 ^ 2 + t260, t31, t247 - t58, t123, -t85 * t95 + t157, t85 * t93 + t146, -t49 * t93 + t157 + 0.2e1 * t178 - t257, pkin(4) * t57 - qJ(5) * t58 + (t30 - t37) * t95 + (t29 - t224) * t93, -t33 * t93 + t49 * t95 + 0.2e1 * t110 + 0.2e1 * t119 - t146, -pkin(4) * t12 + qJ(5) * t9 + t224 * t30 - t29 * t37 - t33 * t49, -t269, -t274, -t264, -t206 * t46 + t11, t123, -t165 * t123 - t34 * t169 - (t132 * t225 + t135 * t26) * t206 + (-t166 * t206 + t2) * qJD(6) + t273, t166 * t123 - t34 * t46 - t270 - (-t132 * t26 + t135 * t225) * t206 + (-t165 * t206 - t244) * qJD(6) + t204; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t123 + t256, t31, -t122 ^ 2 - t260, t257 + t268, 0, 0, 0, 0, 0, -t123 * t135 - t132 * t180 - t169 * t95, t123 * t132 - t135 * t180 - t46 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t269, t274, t264, t266 * t46 - t189, -t123, t2 * t266 - t273, t1 * t206 - t158 + t270;];
tauc_reg  = t3;
