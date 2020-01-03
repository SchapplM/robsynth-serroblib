% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauc_reg [5x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPR9_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR9_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:24:38
% EndTime: 2019-12-31 21:24:50
% DurationCPUTime: 3.64s
% Computational Cost: add. (3585->335), mult. (9171->499), div. (0->0), fcn. (6532->8), ass. (0->190)
t180 = sin(qJ(5));
t183 = cos(qJ(5));
t232 = qJD(5) * t180;
t181 = sin(qJ(3));
t182 = sin(qJ(2));
t238 = qJD(1) * t182;
t217 = t181 * t238;
t184 = cos(qJ(3));
t229 = t184 * qJD(2);
t145 = t217 - t229;
t230 = t181 * qJD(2);
t147 = t184 * t238 + t230;
t178 = sin(pkin(9));
t179 = cos(pkin(9));
t199 = -t179 * t145 - t178 * t147;
t264 = t183 * t199;
t185 = cos(qJ(2));
t226 = qJD(1) * qJD(2);
t214 = t185 * t226;
t225 = qJD(2) * qJD(3);
t109 = -qJD(3) * t217 + (t214 + t225) * t184;
t218 = t185 * t230;
t233 = qJD(3) * t184;
t219 = t182 * t233;
t278 = t218 + t219;
t110 = qJD(1) * t278 + t181 * t225;
t60 = -t178 * t109 - t179 * t110;
t61 = t179 * t109 - t178 * t110;
t94 = t178 * t145 - t179 * t147;
t10 = qJD(5) * t264 + t180 * t60 + t183 * t61 + t232 * t94;
t227 = t185 * qJD(1);
t165 = -qJD(3) + t227;
t158 = -qJD(5) + t165;
t46 = t180 * t94 + t264;
t262 = t46 * t158;
t290 = t10 + t262;
t277 = t180 * t199 - t183 * t94;
t289 = t277 * t46;
t288 = t277 ^ 2 - t46 ^ 2;
t172 = pkin(6) * t227;
t156 = qJD(2) * pkin(7) + t172;
t151 = -t185 * pkin(2) - t182 * pkin(7) - pkin(1);
t135 = t151 * qJD(1);
t254 = t181 * t135;
t103 = t184 * t156 + t254;
t71 = -t145 * qJ(4) + t103;
t266 = t179 * t71;
t102 = t184 * t135 - t181 * t156;
t70 = -t147 * qJ(4) + t102;
t62 = -t165 * pkin(3) + t70;
t25 = t178 * t62 + t266;
t275 = pkin(8) * t199;
t16 = t25 + t275;
t15 = t16 * t232;
t169 = t182 * t226;
t203 = pkin(2) * t182 - pkin(7) * t185;
t149 = t203 * qJD(2);
t136 = qJD(1) * t149;
t206 = pkin(6) * t169;
t246 = -t184 * t136 - t181 * t206;
t190 = -qJD(3) * t103 - t246;
t23 = pkin(3) * t169 - t109 * qJ(4) - t147 * qJD(4) + t190;
t235 = qJD(3) * t181;
t196 = t135 * t233 + t181 * t136 - t156 * t235;
t188 = -t184 * t206 + t196;
t28 = -t110 * qJ(4) - t145 * qJD(4) + t188;
t6 = -t178 * t28 + t179 * t23;
t4 = pkin(4) * t169 - t61 * pkin(8) + t6;
t155 = -qJD(2) * pkin(2) + pkin(6) * t238;
t107 = t145 * pkin(3) + qJD(4) + t155;
t52 = -pkin(4) * t199 + t107;
t287 = -t180 * t4 - t52 * t46 + t15;
t11 = qJD(5) * t277 + t180 * t61 - t183 * t60;
t263 = t277 * t158;
t285 = -t11 - t263;
t148 = t203 * qJD(1);
t131 = t181 * t148;
t271 = -qJ(4) - pkin(7);
t210 = qJD(3) * t271;
t228 = t184 * qJD(4);
t251 = t182 * t184;
t252 = t181 * t185;
t284 = t181 * t210 + t228 - t131 - (-pkin(6) * t251 - qJ(4) * t252) * qJD(1);
t250 = t184 * t185;
t195 = pkin(3) * t182 - qJ(4) * t250;
t243 = pkin(6) * t217 + t184 * t148;
t283 = -qJD(1) * t195 - t181 * qJD(4) + t184 * t210 - t243;
t7 = t178 * t23 + t179 * t28;
t5 = t60 * pkin(8) + t7;
t221 = -t180 * t5 + t183 * t4;
t282 = -t52 * t277 + t221;
t281 = t94 * pkin(8);
t139 = t178 * t184 + t179 * t181;
t192 = t139 * t185;
t280 = qJD(1) * t192 - t139 * qJD(3);
t198 = t178 * t181 - t179 * t184;
t279 = t165 * t198;
t276 = -0.2e1 * t226;
t268 = -t284 * t178 + t283 * t179;
t267 = t283 * t178 + t284 * t179;
t273 = -t172 + (-t181 * t227 + t235) * pkin(3);
t272 = pkin(3) * t178;
t175 = t182 * pkin(6);
t200 = -t180 * t139 - t183 * t198;
t270 = qJD(5) * t200 + t280 * t180 + t279 * t183;
t92 = t183 * t139 - t180 * t198;
t269 = qJD(5) * t92 + t279 * t180 - t280 * t183;
t167 = pkin(6) * t250;
t244 = t184 * t149 + t230 * t175;
t39 = -t182 * t228 + t195 * qJD(2) + (-t167 + (qJ(4) * t182 - t151) * t181) * qJD(3) + t244;
t245 = t181 * t149 + t151 * t233;
t48 = (-pkin(6) * qJD(2) - qJ(4) * qJD(3)) * t251 + (-qJD(4) * t182 + (-pkin(6) * qJD(3) - qJ(4) * qJD(2)) * t185) * t181 + t245;
t13 = t178 * t39 + t179 * t48;
t66 = t178 * t71;
t31 = t179 * t70 - t66;
t24 = t179 * t62 - t66;
t14 = -t165 * pkin(4) + t24 + t281;
t265 = t183 * t14;
t242 = t181 * t151 + t167;
t253 = t181 * t182;
t104 = -qJ(4) * t253 + t242;
t141 = t184 * t151;
t99 = -qJ(4) * t251 + t141 + (-pkin(6) * t181 - pkin(3)) * t185;
t50 = t179 * t104 + t178 * t99;
t261 = -t280 * pkin(4) + t273;
t260 = t109 * t181;
t259 = t145 * t165;
t258 = t147 * t165;
t257 = t155 * t181;
t256 = t155 * t184;
t255 = t165 * t184;
t187 = qJD(1) ^ 2;
t249 = t185 * t187;
t186 = qJD(2) ^ 2;
t248 = t186 * t182;
t247 = t186 * t185;
t153 = t271 * t181;
t154 = t271 * t184;
t106 = t178 * t153 - t179 * t154;
t240 = pkin(3) * t253 + t175;
t176 = t182 ^ 2;
t239 = -t185 ^ 2 + t176;
t237 = qJD(2) * t182;
t236 = qJD(2) * t185;
t234 = qJD(3) * t182;
t231 = t155 * qJD(3);
t223 = t278 * pkin(3) + pkin(6) * t236;
t222 = -t184 * pkin(3) - pkin(2);
t220 = t181 * t234;
t216 = t185 * t229;
t93 = t110 * pkin(3) + pkin(6) * t214;
t212 = qJD(5) * t14 + t5;
t12 = -t178 * t48 + t179 * t39;
t30 = -t178 * t70 - t266;
t49 = -t178 * t104 + t179 * t99;
t209 = pkin(1) * t276;
t105 = t179 * t153 + t178 * t154;
t208 = t145 + t229;
t207 = -t147 + t230;
t82 = -pkin(8) * t198 + t106;
t205 = pkin(4) * t238 + t279 * pkin(8) + qJD(5) * t82 - t268;
t81 = -t139 * pkin(8) + t105;
t204 = t280 * pkin(8) + qJD(5) * t81 + t267;
t2 = t180 * t14 + t183 * t16;
t122 = t198 * t182;
t32 = -t185 * pkin(4) + t122 * pkin(8) + t49;
t121 = t139 * t182;
t34 = -t121 * pkin(8) + t50;
t202 = t180 * t32 + t183 * t34;
t201 = -t183 * t121 + t180 * t122;
t74 = -t180 * t121 - t183 * t122;
t197 = qJD(1) * t176 - t165 * t185;
t168 = t179 * pkin(3) + pkin(4);
t194 = t180 * t168 + t183 * t272;
t193 = t183 * t168 - t180 * t272;
t115 = pkin(4) * t198 + t222;
t100 = t121 * pkin(4) + t240;
t76 = t139 * t234 + t178 * t218 - t179 * t216;
t75 = -qJD(2) * t192 + t198 * t234;
t63 = t147 * pkin(3) - t94 * pkin(4);
t51 = -t75 * pkin(4) + t223;
t33 = -t60 * pkin(4) + t93;
t20 = t31 + t281;
t19 = t30 - t275;
t18 = qJD(5) * t74 - t180 * t76 - t183 * t75;
t17 = qJD(5) * t201 + t180 * t75 - t183 * t76;
t9 = t75 * pkin(8) + t13;
t8 = pkin(4) * t237 + t76 * pkin(8) + t12;
t1 = -t180 * t16 + t265;
t3 = [0, 0, 0, 0.2e1 * t185 * t169, t239 * t276, t247, -t248, 0, -pkin(6) * t247 + t182 * t209, pkin(6) * t248 + t185 * t209, t109 * t251 + (t216 - t220) * t147, (-t145 * t184 - t147 * t181) * t236 + (-t260 - t110 * t184 + (t145 * t181 - t147 * t184) * qJD(3)) * t182, t165 * t220 - t109 * t185 + (t147 * t182 + t184 * t197) * qJD(2), t165 * t219 + t110 * t185 + (-t145 * t182 - t181 * t197) * qJD(2), (-t165 - t227) * t237, -(-t151 * t235 + t244) * t165 + (t184 * t231 + pkin(6) * t110 + (qJD(1) * t141 + t102) * qJD(2)) * t182 + ((pkin(6) * t145 + t257) * qJD(2) + (t254 + (pkin(6) * t165 + t156) * t184) * qJD(3) + t246) * t185, (-pkin(6) * t185 * t235 + t245) * t165 + t196 * t185 + (pkin(6) * t109 - t181 * t231) * t182 + ((pkin(6) * t147 + t256) * t185 + (-pkin(6) * t255 - qJD(1) * t242 - t103) * t182) * qJD(2), t12 * t94 - t7 * t121 + t6 * t122 + t13 * t199 + t24 * t76 + t25 * t75 - t49 * t61 + t50 * t60, t107 * t223 + t24 * t12 + t25 * t13 + t240 * t93 + t6 * t49 + t7 * t50, t10 * t74 + t17 * t277, t10 * t201 - t74 * t11 + t17 * t46 - t18 * t277, -t10 * t185 - t17 * t158 + (qJD(1) * t74 + t277) * t237, t11 * t185 + t18 * t158 + (qJD(1) * t201 + t46) * t237, (-t158 - t227) * t237, -(-t180 * t9 + t183 * t8) * t158 - t221 * t185 - t51 * t46 + t100 * t11 - t33 * t201 + t52 * t18 + (t158 * t202 + t185 * t2) * qJD(5) + ((-t180 * t34 + t183 * t32) * qJD(1) + t1) * t237, t100 * t10 - t15 * t185 + t52 * t17 + t33 * t74 + t51 * t277 + ((-qJD(5) * t34 + t8) * t158 + t4 * t185) * t180 + ((qJD(5) * t32 + t9) * t158 + t212 * t185) * t183 + (-qJD(1) * t202 - t2) * t237; 0, 0, 0, -t182 * t249, t239 * t187, 0, 0, 0, t187 * pkin(1) * t182, pkin(1) * t249, -t147 * t255 + t260, (t109 + t259) * t184 + (-t110 + t258) * t181, -t165 * t233 + (t165 * t250 + t182 * t207) * qJD(1), t165 * t235 + (-t165 * t252 + t182 * t208) * qJD(1), t165 * t238, -pkin(2) * t110 + t243 * t165 + (pkin(7) * t255 + t257) * qJD(3) + ((-pkin(7) * t230 - t102) * t182 + (-pkin(6) * t208 - t257) * t185) * qJD(1), -pkin(2) * t109 - t131 * t165 + (-t181 * pkin(7) * t165 + t256) * qJD(3) + (-t155 * t250 + (-pkin(7) * t229 + t103) * t182 + (t165 * t251 + t185 * t207) * pkin(6)) * qJD(1), -t105 * t61 + t106 * t60 - t6 * t139 - t7 * t198 + t267 * t199 - t279 * t24 + t280 * t25 + t268 * t94, t6 * t105 + t7 * t106 + t273 * t107 + t93 * t222 + t268 * t24 + t267 * t25, t10 * t92 + t270 * t277, t10 * t200 - t92 * t11 - t269 * t277 + t270 * t46, -t270 * t158 + (qJD(2) * t92 - t277) * t238, t269 * t158 + (qJD(2) * t200 - t46) * t238, t158 * t238, t115 * t11 - t33 * t200 + t269 * t52 - t261 * t46 + (t180 * t204 + t183 * t205) * t158 + ((-t180 * t82 + t183 * t81) * qJD(2) - t1) * t238, t115 * t10 + t33 * t92 + t270 * t52 + t261 * t277 + (-t180 * t205 + t183 * t204) * t158 + (-(t180 * t81 + t183 * t82) * qJD(2) + t2) * t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147 * t145, -t145 ^ 2 + t147 ^ 2, t109 - t259, -t110 - t258, t169, -t103 * t165 - t155 * t147 + t190, -t102 * t165 + t155 * t145 - t188, (t178 * t60 - t179 * t61) * pkin(3) + (t24 - t31) * t199 + (-t25 - t30) * t94, -t24 * t30 - t25 * t31 + (-t107 * t147 + t178 * t7 + t179 * t6) * pkin(3), -t289, t288, t290, t285, t169, t193 * t169 + (-t180 * t20 + t183 * t19) * t158 + t63 * t46 + (t158 * t194 - t2) * qJD(5) + t282, -t194 * t169 - t183 * t5 - (t180 * t19 + t183 * t20) * t158 - t63 * t277 + (t158 * t193 - t265) * qJD(5) + t287; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t199 ^ 2 - t94 ^ 2, -t199 * t25 - t24 * t94 + t93, 0, 0, 0, 0, 0, t11 - t263, t10 - t262; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t289, t288, t290, t285, t169, (-qJD(5) - t158) * t2 + t282, -t1 * t158 - t183 * t212 + t287;];
tauc_reg = t3;
