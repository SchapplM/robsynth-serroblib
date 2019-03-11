% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRRP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:05:16
% EndTime: 2019-03-09 06:05:25
% DurationCPUTime: 3.23s
% Computational Cost: add. (5364->395), mult. (12308->506), div. (0->0), fcn. (8132->8), ass. (0->197)
t167 = sin(qJ(4));
t170 = cos(qJ(3));
t237 = qJD(1) * t170;
t220 = t167 * t237;
t279 = -pkin(9) - pkin(8);
t221 = qJD(4) * t279;
t168 = sin(qJ(3));
t196 = pkin(3) * t168 - pkin(8) * t170;
t124 = t196 * qJD(1);
t169 = cos(qJ(4));
t152 = sin(pkin(10)) * pkin(1) + pkin(7);
t133 = t152 * qJD(1);
t287 = qJD(2) * t170 - t168 * t133;
t254 = t167 * t124 + t169 * t287;
t300 = pkin(9) * t220 + t167 * t221 - t254;
t245 = t169 * t170;
t191 = pkin(4) * t168 - pkin(9) * t245;
t205 = t169 * t124 - t167 * t287;
t299 = -qJD(1) * t191 + t169 * t221 - t205;
t238 = qJD(1) * t168;
t214 = t167 * t238;
t233 = qJD(3) * t169;
t119 = -t214 + t233;
t235 = qJD(3) * t167;
t120 = t169 * t238 + t235;
t166 = sin(qJ(5));
t278 = cos(qJ(5));
t182 = t166 * t119 + t278 * t120;
t265 = qJD(3) * pkin(3);
t88 = -t287 - t265;
t65 = -pkin(4) * t119 + t88;
t68 = -t278 * t119 + t120 * t166;
t23 = pkin(5) * t68 - qJ(6) * t182 + t65;
t298 = t23 * t68;
t297 = t65 * t68;
t248 = t166 * t167;
t180 = t278 * t169 - t248;
t101 = t180 * t168;
t150 = -qJD(4) + t237;
t140 = -qJD(5) + t150;
t226 = qJD(1) * qJD(3);
t155 = t168 * t226;
t232 = qJD(3) * t170;
t199 = t278 * t232;
t219 = t167 * t232;
t122 = t166 * t169 + t278 * t167;
t284 = qJD(4) + qJD(5);
t74 = t284 * t122;
t46 = t166 * t219 + t168 * t74 - t169 * t199;
t296 = t101 * t155 + t46 * t140;
t212 = t278 * qJD(5);
t285 = t278 * qJD(4) + t212;
t267 = -t169 * t285 + t180 * t237 + t248 * t284;
t266 = -t122 * t237 + t74;
t274 = t182 * t68;
t234 = qJD(3) * t168;
t230 = qJD(4) * t168;
t218 = t167 * t230;
t176 = t169 * t232 - t218;
t225 = qJD(3) * qJD(4);
t175 = qJD(1) * t176 + t169 * t225;
t211 = qJD(1) * t230;
t222 = qJD(1) * t219 + t167 * t225 + t169 * t211;
t228 = qJD(5) * t166;
t25 = -t119 * t212 + t120 * t228 + t166 * t222 - t278 * t175;
t295 = t170 * t25 + t182 * t234;
t229 = qJD(4) * t169;
t294 = t168 * t229 + t219;
t280 = t182 ^ 2;
t293 = -t68 ^ 2 + t280;
t13 = -t140 * t68 - t25;
t39 = pkin(5) * t182 + qJ(6) * t68;
t275 = t23 * t182;
t291 = t65 * t182;
t136 = t279 * t167;
t137 = t279 * t169;
t181 = t278 * t136 + t166 * t137;
t290 = t181 * qJD(5) + t299 * t166 + t278 * t300;
t81 = t166 * t136 - t278 * t137;
t289 = t81 * qJD(5) + t166 * t300 - t299 * t278;
t154 = -cos(pkin(10)) * pkin(1) - pkin(2);
t114 = -pkin(3) * t170 - pkin(8) * t168 + t154;
t92 = t114 * qJD(1);
t262 = t167 * t92;
t161 = t168 * qJD(2);
t97 = t170 * t133 + t161;
t89 = qJD(3) * pkin(8) + t97;
t54 = t169 * t89 + t262;
t45 = pkin(9) * t119 + t54;
t26 = qJD(5) * t182 + t166 * t175 + t278 * t222;
t208 = -t170 * t26 + t68 * t234;
t231 = qJD(4) * t167;
t198 = -t97 + (-t220 + t231) * pkin(4);
t288 = t170 * t222;
t100 = t122 * t168;
t247 = t167 * t168;
t47 = t167 * t199 - t166 * t218 - t228 * t247 + (t166 * t232 + t168 * t285) * t169;
t269 = -t100 * t155 + t47 * t140;
t286 = t208 + t269;
t283 = t295 - t296;
t282 = -t140 * t182 - t26;
t103 = t169 * t114;
t246 = t168 * t169;
t249 = t152 * t167;
t60 = -pkin(9) * t246 + t103 + (-pkin(4) - t249) * t170;
t123 = t152 * t245;
t240 = t167 * t114 + t123;
t64 = -pkin(9) * t247 + t240;
t185 = t166 * t60 + t278 * t64;
t177 = t191 * qJD(3);
t127 = t196 * qJD(3);
t241 = t169 * t127 + t234 * t249;
t32 = t177 + (-t123 + (pkin(9) * t168 - t114) * t167) * qJD(4) + t241;
t217 = t170 * t231;
t256 = t114 * t229 + t167 * t127;
t34 = (-t168 * t233 - t217) * t152 - t294 * pkin(9) + t256;
t281 = -qJD(5) * t185 - t166 * t34 + t278 * t32;
t277 = pkin(8) * t150;
t273 = qJ(6) * t238 - t290;
t272 = pkin(5) * t238 + t289;
t271 = -t101 * t26 + t46 * t68;
t270 = t266 * pkin(5) + t267 * qJ(6) - qJD(6) * t122 + t198;
t264 = t166 * t45;
t263 = t167 * t88;
t91 = qJD(3) * t161 + t133 * t232;
t259 = t91 * t167;
t258 = t91 * t169;
t53 = -t167 * t89 + t169 * t92;
t44 = -pkin(9) * t120 + t53;
t37 = -pkin(4) * t150 + t44;
t9 = t278 * t37 - t264;
t257 = qJD(6) - t9;
t12 = t278 * t44 - t264;
t255 = -pkin(4) * t212 - qJD(6) + t12;
t253 = qJD(3) * t181;
t252 = qJD(3) * t81;
t251 = t119 * t168;
t250 = t120 * t150;
t244 = t170 * t150;
t171 = qJD(3) ^ 2;
t243 = t171 * t168;
t242 = t171 * t170;
t106 = pkin(4) * t247 + t168 * t152;
t162 = t168 ^ 2;
t239 = -t170 ^ 2 + t162;
t134 = qJD(1) * t154;
t227 = t162 * qJD(1);
t223 = t278 * t45;
t75 = t294 * pkin(4) + t152 * t232;
t160 = -pkin(4) * t169 - pkin(3);
t215 = t150 * t229;
t113 = qJD(1) * t127;
t90 = t287 * qJD(3);
t210 = -t169 * t113 + t167 * t90;
t16 = qJD(1) * t177 - t45 * qJD(4) - t210;
t184 = t167 * t113 + t169 * t90 + t92 * t229 - t89 * t231;
t19 = -t222 * pkin(9) + t184;
t207 = -t166 * t16 - t278 * t19 - t37 * t212 + t45 * t228;
t206 = -t278 * t16 + t166 * t19 + t45 * t212 + t37 * t228;
t204 = t150 + t237;
t203 = -t119 + t233;
t202 = qJD(4) + t237;
t201 = pkin(5) * t155;
t200 = t168 * t215;
t11 = t166 * t44 + t223;
t197 = pkin(4) * t228 - t11;
t195 = -t100 * t25 + t182 * t47;
t194 = -t150 + t202;
t192 = 0.2e1 * qJD(3) * t134;
t129 = t140 * qJD(6);
t146 = qJ(6) * t155;
t1 = t146 - t129 - t207;
t190 = -t140 * t9 + t207;
t10 = t166 * t37 + t223;
t189 = -t10 * t140 - t206;
t187 = -t166 * t64 + t278 * t60;
t183 = t166 * t32 + t60 * t212 - t64 * t228 + t278 * t34;
t179 = t202 * t235;
t178 = (-t227 + t244) * t167;
t56 = t222 * pkin(4) + t91;
t2 = -t201 + t206;
t172 = qJD(1) ^ 2;
t159 = -t278 * pkin(4) - pkin(5);
t153 = pkin(4) * t166 + qJ(6);
t104 = t120 * t234;
t66 = -pkin(5) * t180 - qJ(6) * t122 + t160;
t50 = pkin(5) * t100 - qJ(6) * t101 + t106;
t31 = pkin(4) * t120 + t39;
t28 = t170 * pkin(5) - t187;
t27 = -qJ(6) * t170 + t185;
t8 = -t140 * qJ(6) + t10;
t7 = t140 * pkin(5) + t257;
t6 = pkin(5) * t47 + qJ(6) * t46 - qJD(6) * t101 + t75;
t5 = -pkin(5) * t234 - t281;
t4 = t26 * pkin(5) + t25 * qJ(6) - qJD(6) * t182 + t56;
t3 = qJ(6) * t234 - qJD(6) * t170 + t183;
t14 = [0, 0, 0, 0, 0.2e1 * t170 * t155, -0.2e1 * t239 * t226, t242, -t243, 0, -t152 * t242 + t168 * t192, t152 * t243 + t170 * t192, t120 * t176 + t175 * t246 (t119 * t169 - t120 * t167) * t232 + ((-t119 + t214) * t231 + (-t120 * qJD(4) - t179 - t222) * t169) * t168, t104 + t204 * t218 + (t227 + (-t150 - t202) * t170) * t233, t200 + t288 + (t178 + t251) * qJD(3), -t204 * t234 -(-t114 * t231 + t241) * t150 + ((-t119 * t152 + t263) * qJD(3) + (t262 + (t150 * t152 + t89) * t169) * qJD(4) + t210) * t170 + (t152 * t222 + t259 + t88 * t229 + ((-t170 * t249 + t103) * qJD(1) + t53) * qJD(3)) * t168 (-t152 * t217 + t256) * t150 + t184 * t170 + (t258 + (-t152 * t238 - t88) * t231) * t168 + ((t120 * t152 + t88 * t169) * t170 + (t152 * t169 * t194 - qJD(1) * t240 - t54) * t168) * qJD(3), -t101 * t25 - t182 * t46, -t195 + t271, t295 + t296, t269 - t208 (-t140 - t237) * t234, t56 * t100 + t106 * t26 - t140 * t281 + t187 * t155 + t206 * t170 + t9 * t234 + t65 * t47 + t75 * t68, t183 * t140 - t207 * t170 + t75 * t182 - t106 * t25 + t56 * t101 - t65 * t46 + (-t185 * qJD(1) - t10) * t234, t100 * t4 + t140 * t5 + t170 * t2 + t23 * t47 + t26 * t50 + t6 * t68 + (-qJD(1) * t28 - t7) * t234, -t1 * t100 + t101 * t2 + t182 * t5 - t25 * t28 - t26 * t27 - t3 * t68 - t46 * t7 - t47 * t8, -t1 * t170 - t101 * t4 - t140 * t3 + t23 * t46 + t25 * t50 - t6 * t182 + (qJD(1) * t27 + t8) * t234, t1 * t27 + t2 * t28 + t23 * t6 + t3 * t8 + t4 * t50 + t5 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t243, -t242, 0, 0, 0, 0, 0, t200 - t288 + (t178 - t251) * qJD(3), t104 + (-t150 + t237) * t218 + (-t170 * t194 - t227) * t233, 0, 0, 0, 0, 0, t286, t283, t286, t195 + t271, -t283, t1 * t101 + t100 * t2 - t170 * t4 + t23 * t234 - t46 * t8 + t47 * t7; 0, 0, 0, 0, -t170 * t172 * t168, t239 * t172, 0, 0, 0, qJD(3) * t97 - t134 * t238 - t91, -t134 * t237, -t167 ^ 2 * t211 + (t179 - t250) * t169 (-t222 + t250) * t167 + ((t119 + t233) * qJD(4) + (t170 * t203 - t218) * qJD(1)) * t169, -t215 + (t169 * t244 + (-t120 + t235) * t168) * qJD(1), t150 * t231 + (-t167 * t244 + t168 * t203) * qJD(1), t150 * t238, -pkin(3) * t222 - t258 + t205 * t150 + t97 * t119 + (t169 * t277 + t263) * qJD(4) + (-t53 * t168 + (-pkin(8) * t234 - t88 * t170) * t167) * qJD(1), t259 - t254 * t150 - t97 * t120 + (-t167 * t277 + (t88 - t265) * t169) * qJD(4) + ((-t88 - t265) * t245 + (pkin(3) * t231 - pkin(8) * t233 + t54) * t168) * qJD(1), -t25 * t122 - t182 * t267, -t122 * t26 - t180 * t25 - t182 * t266 + t267 * t68, t267 * t140 + (qJD(3) * t122 - t182) * t238, t266 * t140 + (qJD(3) * t180 + t68) * t238, t140 * t238, -t56 * t180 + t160 * t26 + t198 * t68 + t266 * t65 + t289 * t140 + (-t9 + t253) * t238, t56 * t122 - t160 * t25 + t198 * t182 - t267 * t65 + t290 * t140 + (t10 - t252) * t238, -t180 * t4 + t26 * t66 + t270 * t68 + t266 * t23 + t272 * t140 + (t7 + t253) * t238, t1 * t180 + t122 * t2 + t181 * t25 + t182 * t272 - t26 * t81 - t266 * t8 - t267 * t7 + t273 * t68, -t122 * t4 + t25 * t66 - t270 * t182 + t267 * t23 + t273 * t140 + (-t8 + t252) * t238, t1 * t81 - t181 * t2 + t270 * t23 + t272 * t7 - t273 * t8 + t4 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120 * t119, -t119 ^ 2 + t120 ^ 2, t119 * t150 + t175, -t222 - t250, t155, -t88 * t120 - t210 + (-qJD(4) - t150) * t54, -t119 * t88 - t150 * t53 - t184, t274, t293, t13, t282, t155, -t11 * t140 - t291 + (-t120 * t68 + t140 * t228 + t278 * t155) * pkin(4) - t206, -t12 * t140 + t297 + (-t120 * t182 + t140 * t212 - t155 * t166) * pkin(4) + t207, -t275 - t31 * t68 + t197 * t140 + (pkin(5) - t159) * t155 - t206, -t153 * t26 - t159 * t25 + (t197 + t8) * t182 + (t255 + t7) * t68, t255 * t140 + t153 * t155 + t182 * t31 + t1 - t298, t1 * t153 + t159 * t2 + t197 * t7 - t23 * t31 - t255 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t274, t293, t13, t282, t155, t189 - t291, t190 + t297, -t39 * t68 + t189 + 0.2e1 * t201 - t275, pkin(5) * t25 - t26 * qJ(6) + (-t10 + t8) * t182 + (t7 - t257) * t68, t182 * t39 - 0.2e1 * t129 + 0.2e1 * t146 - t190 - t298, -t2 * pkin(5) + t1 * qJ(6) - t7 * t10 - t23 * t39 + t257 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155 + t274, t13, -t140 ^ 2 - t280, t140 * t8 + t2 + t275;];
tauc_reg  = t14;
