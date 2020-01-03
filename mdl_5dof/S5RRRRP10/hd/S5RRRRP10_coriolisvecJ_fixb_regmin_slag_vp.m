% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% 
% Output:
% tauc_reg [5x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRP10_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP10_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:12:21
% EndTime: 2019-12-31 22:12:34
% DurationCPUTime: 4.40s
% Computational Cost: add. (4269->380), mult. (11490->550), div. (0->0), fcn. (8665->8), ass. (0->180)
t151 = sin(pkin(5));
t158 = cos(qJ(2));
t218 = qJD(1) * t158;
t143 = t151 * t218;
t174 = t143 - qJD(3);
t155 = sin(qJ(2));
t152 = cos(pkin(5));
t210 = t152 * qJD(1);
t206 = pkin(1) * t210;
t141 = t155 * t206;
t115 = pkin(7) * t143 + t141;
t154 = sin(qJ(3));
t157 = cos(qJ(3));
t259 = -t115 - t174 * (pkin(3) * t154 - pkin(9) * t157);
t184 = qJD(2) + t210;
t219 = qJD(1) * t155;
t200 = t151 * t219;
t161 = -t154 * t200 + t157 * t184;
t208 = qJD(1) * qJD(2);
t195 = t151 * t208;
t177 = t158 * t195;
t171 = t157 * t177;
t160 = qJD(3) * t161 + t171;
t258 = -qJD(4) * t174 + t160;
t245 = -qJ(5) - pkin(9);
t257 = qJ(5) * t161 + qJD(4) * t245;
t153 = sin(qJ(4));
t156 = cos(qJ(4));
t215 = qJD(3) * t154;
t205 = pkin(8) * t215;
t256 = t153 * t205 + t259 * t156;
t255 = t155 * t158;
t181 = t157 * t143;
t213 = qJD(3) * t157;
t254 = t181 - t213;
t136 = -t157 * pkin(3) - t154 * pkin(9) - pkin(2);
t211 = qJD(4) * t156;
t112 = -pkin(7) * t200 + t158 * t206;
t168 = (pkin(2) * t155 - pkin(8) * t158) * t151;
t113 = qJD(1) * t168;
t225 = t157 * t112 + t154 * t113;
t55 = pkin(9) * t200 + t225;
t253 = -t136 * t211 - t259 * t153 + t156 * t55;
t99 = t154 * t184 + t157 * t200;
t66 = -t153 * t174 + t156 * t99;
t252 = t66 ^ 2;
t159 = qJD(1) ^ 2;
t81 = -pkin(2) * t184 - t112;
t39 = -pkin(3) * t161 - t99 * pkin(9) + t81;
t82 = pkin(8) * t184 + t115;
t111 = (-pkin(2) * t158 - pkin(8) * t155 - pkin(1)) * t151;
t92 = qJD(1) * t111;
t48 = t154 * t92 + t157 * t82;
t42 = -t174 * pkin(9) + t48;
t13 = -t153 * t42 + t156 * t39;
t9 = -t66 * qJ(5) + t13;
t93 = qJD(4) - t161;
t7 = t93 * pkin(4) + t9;
t251 = t7 - t9;
t250 = pkin(1) * t155;
t249 = pkin(1) * t158;
t248 = pkin(4) * t153;
t64 = t153 * t99 + t156 * t174;
t247 = t64 * t93;
t246 = t66 * t93;
t227 = t156 * t157;
t146 = pkin(8) * t227;
t182 = t154 * t143;
t209 = t156 * qJD(5);
t232 = qJ(5) * t154;
t220 = qJD(1) * t151;
t226 = t156 * t158;
t87 = (t153 * t155 + t157 * t226) * t220;
t244 = -pkin(4) * t182 + t87 * qJ(5) + t153 * t55 - t154 * t209 + (pkin(4) * t154 - qJ(5) * t227) * qJD(3) + (-t146 + (-t136 + t232) * t153) * qJD(4) + t256;
t228 = t154 * t156;
t86 = t153 * t181 - t156 * t200;
t243 = t86 * qJ(5) + (-pkin(8) * qJD(3) - qJ(5) * qJD(4)) * t228 + (-qJD(5) * t154 + (-pkin(8) * qJD(4) - qJ(5) * qJD(3)) * t157) * t153 - t253;
t47 = -t154 * t82 + t157 * t92;
t58 = t99 * pkin(3) - pkin(9) * t161;
t242 = t153 * t58 + t156 * t47;
t230 = t151 * t155;
t144 = pkin(7) * t230;
t109 = t144 + (-pkin(2) - t249) * t152;
t120 = -t152 * t157 + t154 * t230;
t121 = t152 * t154 + t157 * t230;
t51 = t120 * pkin(3) - t121 * pkin(9) + t109;
t229 = t151 * t158;
t207 = pkin(7) * t229;
t110 = t207 + (pkin(8) + t250) * t152;
t239 = t157 * t110 + t154 * t111;
t53 = -pkin(9) * t229 + t239;
t241 = t153 * t51 + t156 * t53;
t172 = t154 * t177;
t68 = t99 * qJD(3) + t172;
t237 = t153 * t68;
t236 = t156 * t68;
t178 = t155 * t195;
t212 = qJD(4) * t153;
t29 = -t153 * t178 - t258 * t156 + t99 * t212;
t235 = t29 * t153;
t234 = t257 * t153 + t209 - t242;
t57 = t156 * t58;
t233 = -t99 * pkin(4) - t57 + t257 * t156 + (-qJD(5) + t47) * t153;
t148 = t151 ^ 2;
t231 = t148 * t159;
t222 = t153 * t136 + t146;
t221 = t155 ^ 2 - t158 ^ 2;
t217 = qJD(2) * t155;
t216 = qJD(3) * t153;
t214 = qJD(3) * t156;
t204 = t93 * t214;
t203 = t93 * t212;
t202 = t153 * t229;
t201 = pkin(8) + t248;
t199 = t151 * t217;
t198 = qJD(2) * t229;
t197 = t151 * t152 * t159;
t196 = t148 * t208;
t194 = -t153 * t53 + t156 * t51;
t114 = qJD(2) * t168;
t106 = qJD(1) * t114;
t116 = (t152 * t249 - t144) * qJD(2);
t107 = qJD(1) * t116;
t192 = t157 * t106 - t154 * t107 - t82 * t213 - t92 * t215;
t191 = -t154 * t110 + t157 * t111;
t190 = t156 * t93;
t189 = t158 * t174;
t188 = t174 * t151;
t187 = qJD(3) * t174;
t185 = 0.2e1 * t196;
t183 = qJD(2) + 0.2e1 * t210;
t21 = -pkin(3) * t178 - t192;
t179 = pkin(9) * qJD(4) * t93 + t21;
t52 = pkin(3) * t229 - t191;
t176 = -0.2e1 * pkin(1) * t196;
t14 = t153 * t39 + t156 * t42;
t173 = pkin(7) * t177;
t170 = -t110 * t213 - t111 * t215 + t157 * t114 - t154 * t116;
t169 = -t93 * t211 - t237;
t41 = t174 * pkin(3) - t47;
t167 = -pkin(9) * t68 + t93 * t41;
t75 = t121 * t153 + t151 * t226;
t166 = -t154 * t106 - t157 * t107 - t92 * t213 + t82 * t215;
t20 = pkin(9) * t178 - t166;
t28 = t68 * pkin(3) - pkin(9) * t160 + qJD(2) * t141 + t173;
t5 = t153 * t28 + t156 * t20 + t39 * t211 - t42 * t212;
t164 = -t110 * t215 + t111 * t213 + t154 * t114 + t157 * t116;
t24 = pkin(9) * t199 + t164;
t117 = (t152 * t250 + t207) * qJD(2);
t73 = qJD(3) * t121 + t154 * t198;
t74 = -qJD(3) * t120 + t157 * t198;
t33 = t73 * pkin(3) - t74 * pkin(9) + t117;
t165 = t153 * t33 + t156 * t24 + t51 * t211 - t53 * t212;
t163 = pkin(1) * (-t152 * t208 + t231);
t100 = t154 * t112;
t54 = -pkin(3) * t200 - t157 * t113 + t100;
t25 = -pkin(3) * t199 - t170;
t6 = -qJD(4) * t14 - t153 * t20 + t156 * t28;
t162 = -t241 * qJD(4) - t153 * t24 + t156 * t33;
t30 = t258 * t153 - t156 * t178 + t99 * t211;
t8 = t30 * pkin(4) + t21;
t138 = t245 * t156;
t137 = t245 * t153;
t127 = t156 * t136;
t108 = qJD(1) * t117;
t79 = -t153 * t232 + t222;
t76 = t121 * t156 - t202;
t69 = -qJ(5) * t228 + t127 + (-pkin(8) * t153 - pkin(4)) * t157;
t63 = t64 ^ 2;
t38 = -qJD(4) * t75 + t153 * t199 + t74 * t156;
t37 = -qJD(4) * t202 + t121 * t211 + t74 * t153 - t156 * t199;
t22 = t64 * pkin(4) + qJD(5) + t41;
t15 = -t75 * qJ(5) + t241;
t11 = t120 * pkin(4) - t76 * qJ(5) + t194;
t10 = -t64 * qJ(5) + t14;
t4 = -t37 * qJ(5) - t75 * qJD(5) + t165;
t3 = t73 * pkin(4) - t38 * qJ(5) - t76 * qJD(5) + t162;
t2 = -t30 * qJ(5) - t64 * qJD(5) + t5;
t1 = t68 * pkin(4) + t29 * qJ(5) - t66 * qJD(5) + t6;
t12 = [0, 0, 0, t185 * t255, -t221 * t185, t183 * t198, -t183 * t199, 0, -t108 * t152 - t117 * t184 + t155 * t176, -t107 * t152 - t116 * t184 + t158 * t176, t121 * t160 + t99 * t74, -t120 * t160 - t121 * t68 + t161 * t74 - t99 * t73, t121 * t178 - t160 * t229 - t74 * t174 + t99 * t199, t73 * t174 + (t68 * t158 + (-qJD(1) * t120 + t161) * t217) * t151, (-t148 * t218 - t188) * t217, -t170 * t174 - t117 * t161 + t109 * t68 + t108 * t120 + t81 * t73 + (-t192 * t158 + (qJD(1) * t191 + t47) * t217) * t151, t108 * t121 + t109 * t160 + t117 * t99 + t164 * t174 - t166 * t229 - t239 * t178 - t48 * t199 + t81 * t74, -t29 * t76 + t66 * t38, t29 * t75 - t76 * t30 - t66 * t37 - t38 * t64, -t29 * t120 + t38 * t93 + t66 * t73 + t76 * t68, -t30 * t120 - t37 * t93 - t64 * t73 - t75 * t68, t68 * t120 + t93 * t73, t6 * t120 + t13 * t73 + t162 * t93 + t194 * t68 + t21 * t75 + t25 * t64 + t52 * t30 + t41 * t37, -t5 * t120 - t14 * t73 - t165 * t93 + t21 * t76 - t241 * t68 + t25 * t66 - t52 * t29 + t41 * t38, -t1 * t76 - t10 * t37 + t11 * t29 - t15 * t30 - t2 * t75 - t3 * t66 - t38 * t7 - t4 * t64, t2 * t15 + t10 * t4 + t1 * t11 + t7 * t3 + t8 * (pkin(4) * t75 + t52) + t22 * (pkin(4) * t37 + t25); 0, 0, 0, -t231 * t255, t221 * t231, -t158 * t197, t155 * t197, 0, t115 * t184 + t155 * t163 - t173, pkin(7) * t178 + t112 * t184 + t158 * t163, -qJD(3) * t154 ^ 2 * t200 + ((qJD(3) * t184 + t177) * t154 - t174 * t99) * t157, -t154 * t68 + t157 * t160 + (t182 - t215) * t99 - t254 * t161, -t157 * t187 + (t157 * t189 + (t154 * qJD(2) - t99) * t155) * t220, t154 * t187 + (-t154 * t189 + (qJD(2) * t157 - t161) * t155) * t220, t188 * t219, -pkin(2) * t68 + t81 * t215 - t100 * t174 + t115 * t161 + (pkin(8) * t187 + t113 * t174 - t108) * t157 + (-t47 * t155 + (-pkin(8) * t217 - t158 * t81) * t154) * t220, -t157 * pkin(8) * t178 - pkin(2) * t160 + t108 * t154 - t115 * t99 + t48 * t200 - t254 * t81 + (-t205 - t225) * t174, -t29 * t228 + (-t154 * t212 + t156 * t213 - t87) * t66, t87 * t64 + t66 * t86 + (-t153 * t66 - t156 * t64) * t213 + (t235 - t156 * t30 + (t153 * t64 - t156 * t66) * qJD(4)) * t154, -t87 * t93 + (t29 + t204) * t157 + (-t174 * t66 - t203 + t236) * t154, t86 * t93 + (-t216 * t93 + t30) * t157 + (t174 * t64 + t169) * t154, -t154 * t174 * t93 - t68 * t157, t127 * t68 - t41 * t86 - t54 * t64 + ((-qJD(4) * t136 + t55) * t153 + t256) * t93 + (t41 * t216 - t6 + (qJD(3) * t64 + t169) * pkin(8)) * t157 + (pkin(8) * t30 - t13 * t174 + t21 * t153 + t211 * t41) * t154, -t222 * t68 - t54 * t66 - t41 * t87 + t253 * t93 + (t41 * t214 + t5 + (qJD(3) * t66 + t203) * pkin(8)) * t157 + (-t41 * t212 + t21 * t156 + t174 * t14 + (-t29 + t204) * pkin(8)) * t154, t10 * t86 + t69 * t29 - t79 * t30 + t7 * t87 - t244 * t66 - t243 * t64 + (-t10 * t153 - t156 * t7) * t213 + (-t1 * t156 - t153 * t2 + (-t10 * t156 + t153 * t7) * qJD(4)) * t154, t201 * t8 * t154 + t1 * t69 + t243 * t10 + t2 * t79 + t244 * t7 + (t201 * t213 - t54 + (t211 * t154 - t86) * pkin(4)) * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99 * t161, -t161 ^ 2 + t99 ^ 2, t161 * t143 + t171, -t99 * t143 - t172, t178, -t48 * t174 - t81 * t99 + t192, -t161 * t81 - t47 * t174 + t166, t190 * t66 - t235, (-t29 - t247) * t156 + (-t30 - t246) * t153, t190 * t93 - t66 * t99 + t237, -t93 ^ 2 * t153 + t64 * t99 + t236, -t93 * t99, -pkin(3) * t30 - t13 * t99 - t48 * t64 - t57 * t93 - t179 * t156 + (t47 * t93 + t167) * t153, pkin(3) * t29 + t14 * t99 + t153 * t179 + t156 * t167 + t242 * t93 - t48 * t66, t137 * t29 + t138 * t30 - t233 * t66 - t234 * t64 + (-t7 * t93 + t2) * t156 + (-t10 * t93 - t1) * t153, -t2 * t138 + t1 * t137 + t8 * (-pkin(4) * t156 - pkin(3)) + t233 * t7 + (t93 * t248 - t48) * t22 + t234 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66 * t64, -t63 + t252, -t29 + t247, -t30 + t246, t68, t14 * t93 - t41 * t66 + t6, t13 * t93 + t41 * t64 - t5, pkin(4) * t29 - t251 * t64, t251 * t10 + (-t22 * t66 + t1) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63 - t252, t10 * t64 + t7 * t66 + t8;];
tauc_reg = t12;
