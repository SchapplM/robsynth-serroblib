% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauc_reg [6x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRRP5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:05:26
% EndTime: 2019-03-09 12:05:41
% DurationCPUTime: 5.67s
% Computational Cost: add. (9324->453), mult. (27854->629), div. (0->0), fcn. (22317->10), ass. (0->230)
t190 = sin(pkin(6));
t191 = cos(pkin(11));
t198 = cos(qJ(2));
t278 = t191 * t198;
t252 = t190 * t278;
t174 = qJD(1) * t252;
t189 = sin(pkin(11));
t195 = sin(qJ(2));
t268 = qJD(1) * t190;
t251 = t195 * t268;
t223 = t189 * t251 - t174;
t325 = qJD(4) + t223;
t192 = cos(pkin(6));
t312 = pkin(1) * t198;
t257 = t192 * t312;
t181 = qJD(1) * t257;
t310 = pkin(8) + qJ(3);
t247 = t195 * t310;
t234 = t190 * t247;
t144 = -qJD(1) * t234 + t181;
t313 = pkin(1) * t195;
t255 = t192 * t313;
t280 = t190 * t198;
t152 = t310 * t280 + t255;
t145 = t152 * qJD(1);
t279 = t191 * t145;
t101 = t144 * t189 + t279;
t194 = sin(qJ(4));
t197 = cos(qJ(4));
t331 = -t101 + t325 * (pkin(4) * t194 - pkin(10) * t197);
t316 = qJD(4) * t325;
t326 = t223 * t325 + t316;
t330 = t326 * t197;
t226 = t189 * t198 + t191 * t195;
t157 = t226 * t268;
t193 = sin(qJ(5));
t196 = cos(qJ(5));
t273 = t197 * t223;
t118 = t157 * t193 - t196 * t273;
t265 = qJD(4) * t197;
t329 = -t196 * t265 + t118;
t267 = qJD(1) * t192;
t237 = qJD(2) + t267;
t127 = t194 * t157 - t197 * t237;
t266 = qJD(2) * t190;
t250 = t195 * t266;
t232 = qJD(1) * t250;
t151 = qJD(2) * t174 - t189 * t232;
t274 = t197 * t151;
t201 = -qJD(4) * t127 + t274;
t328 = qJD(5) * t325 + t201;
t309 = -qJ(6) - pkin(10);
t327 = -qJ(6) * t127 + qJD(5) * t309;
t324 = pkin(2) * t250;
t117 = -t196 * t157 - t193 * t273;
t263 = qJD(5) * t196;
t229 = t194 * t263 - t117;
t129 = t197 * t157 + t194 * t237;
t87 = t196 * t129 + t193 * t325;
t323 = t229 * t87;
t85 = t129 * t193 - t196 * t325;
t322 = t325 * t85;
t261 = t194 * qJD(4);
t184 = t189 * pkin(2) + pkin(9);
t285 = t184 * t193;
t321 = t331 * t196 + t261 * t285;
t125 = qJD(5) + t127;
t264 = qJD(5) * t193;
t211 = -t194 * t264 - t329;
t320 = t211 * t125;
t186 = t190 ^ 2;
t217 = -pkin(8) * t280 - t255;
t319 = -t186 * t313 + t217 * t192;
t318 = t265 + t273;
t185 = -pkin(2) * t191 - pkin(3);
t169 = -pkin(4) * t197 - pkin(10) * t194 + t185;
t134 = t189 * t145;
t102 = t144 * t191 - t134;
t112 = pkin(2) * t251 + pkin(3) * t157 + pkin(9) * t223;
t292 = t197 * t102 + t194 * t112;
t45 = pkin(10) * t157 + t292;
t317 = -t169 * t263 - t331 * t193 + t196 * t45;
t315 = t87 ^ 2;
t199 = qJD(1) ^ 2;
t130 = qJD(2) * pkin(2) + t181 + (t192 * pkin(2) - t234) * qJD(1);
t80 = t189 * t130 + t279;
t72 = pkin(9) * t237 + t80;
t231 = (-pkin(2) * t198 - pkin(1)) * t190;
t224 = qJD(1) * t231;
t166 = qJD(3) + t224;
t97 = pkin(3) * t223 - pkin(9) * t157 + t166;
t43 = t194 * t97 + t197 * t72;
t33 = pkin(10) * t325 + t43;
t79 = t191 * t130 - t134;
t71 = -pkin(3) * t237 - t79;
t39 = t127 * pkin(4) - t129 * pkin(10) + t71;
t14 = -t193 * t33 + t196 * t39;
t9 = -qJ(6) * t87 + t14;
t8 = pkin(5) * t125 + t9;
t314 = t8 - t9;
t311 = pkin(5) * t193;
t42 = -t194 * t72 + t197 * t97;
t75 = pkin(4) * t129 + pkin(10) * t127;
t308 = t193 * t75 + t196 * t42;
t281 = t190 * t195;
t159 = t189 * t281 - t252;
t143 = (pkin(2) + t312) * t192 - t234;
t110 = t189 * t143 + t191 * t152;
t100 = pkin(9) * t192 + t110;
t160 = t226 * t190;
t116 = pkin(3) * t159 - pkin(9) * t160 + t231;
t293 = t197 * t100 + t194 * t116;
t50 = pkin(10) * t159 + t293;
t137 = t160 * t194 - t192 * t197;
t138 = t160 * t197 + t192 * t194;
t109 = t143 * t191 - t189 * t152;
t99 = -pkin(3) * t192 - t109;
t54 = pkin(4) * t137 - pkin(10) * t138 + t99;
t306 = t193 * t54 + t196 * t50;
t275 = t196 * t197;
t172 = t184 * t275;
t262 = qJD(6) * t196;
t287 = t223 * t194;
t290 = qJ(6) * t194;
t305 = pkin(5) * t287 + qJ(6) * t118 + t193 * t45 - t194 * t262 + (pkin(5) * t194 - qJ(6) * t275) * qJD(4) + (-t172 + (-t169 + t290) * t193) * qJD(5) + t321;
t276 = t194 * t196;
t304 = qJ(6) * t117 + (-qJ(6) * qJD(5) - qJD(4) * t184) * t276 + (-qJD(6) * t194 + (-qJ(6) * qJD(4) - qJD(5) * t184) * t197) * t193 - t317;
t303 = t125 * t85;
t302 = t223 * t87;
t215 = qJD(2) * t226;
t155 = t190 * t215;
t205 = qJD(1) * t155;
t177 = t181 * qJD(2);
t204 = (-qJD(2) * t247 + qJD(3) * t198) * t190;
t122 = qJD(1) * t204 + t177;
t132 = -t152 * qJD(2) - qJD(3) * t281;
t200 = qJD(1) * t132;
t62 = t191 * t122 + t189 * t200;
t176 = pkin(2) * t232;
t206 = qJD(2) * t157;
t98 = pkin(3) * t206 - t151 * pkin(9) + t176;
t241 = t194 * t62 - t197 * t98 + t97 * t261 + t72 * t265;
t17 = -pkin(4) * t205 + t241;
t301 = t17 * t193;
t300 = t17 * t196;
t277 = t194 * t151;
t84 = t129 * qJD(4) + t277;
t299 = t193 * t84;
t298 = t196 * t84;
t36 = t129 * t263 + t328 * t193 - t196 * t205;
t297 = t197 * t36;
t296 = t87 * t125;
t295 = t327 * t193 + t262 - t308;
t69 = t196 * t75;
t294 = -pkin(5) * t129 - t69 + t327 * t196 + (-qJD(6) + t42) * t193;
t289 = qJD(5) * t85;
t238 = t125 * t196;
t288 = t129 * t157;
t286 = t157 * t127;
t284 = t184 * t197;
t283 = t186 * t199;
t282 = t189 * t195;
t270 = t193 * t169 + t172;
t269 = t195 ^ 2 - t198 ^ 2;
t259 = qJD(4) - t174;
t258 = qJD(1) * qJD(2);
t254 = t87 * t265;
t253 = t198 * t283;
t249 = t190 * t192 * t199;
t246 = t186 * t258;
t245 = t184 + t311;
t244 = -t193 * t50 + t196 * t54;
t35 = t129 * t264 - t193 * t205 - t328 * t196;
t243 = t197 * t35 + t87 * t261;
t240 = -t194 * t102 + t112 * t197;
t239 = -t194 * t100 + t116 * t197;
t61 = t189 * t122 - t191 * t200;
t182 = qJD(2) * t257;
t131 = t182 + t204;
t73 = t131 * t189 - t191 * t132;
t236 = qJD(2) + 0.2e1 * t267;
t235 = -t36 * t276 + t329 * t85;
t233 = t198 * t246;
t15 = t193 * t39 + t196 * t33;
t10 = -qJ(6) * t85 + t15;
t228 = t10 * t196 - t193 * t8;
t227 = -t10 * t193 - t196 * t8;
t115 = t138 * t196 + t159 * t193;
t114 = t138 * t193 - t159 * t196;
t156 = (t278 - t282) * t266;
t113 = pkin(3) * t155 - pkin(9) * t156 + t324;
t74 = t131 * t191 + t132 * t189;
t225 = -t100 * t265 + t113 * t197 - t116 * t261 - t194 * t74;
t44 = -pkin(4) * t157 - t240;
t49 = -pkin(4) * t159 - t239;
t222 = -t125 * t263 - t299;
t221 = -t194 * t98 - t197 * t62 + t72 * t261 - t97 * t265;
t16 = pkin(10) * t205 - t221;
t28 = t84 * pkin(4) - pkin(10) * t201 + t61;
t5 = t196 * t16 + t193 * t28 + t39 * t263 - t33 * t264;
t216 = -t100 * t261 + t194 * t113 + t116 * t265 + t197 * t74;
t22 = pkin(10) * t155 + t216;
t107 = qJD(4) * t138 + t156 * t194;
t108 = -qJD(4) * t137 + t156 * t197;
t31 = pkin(4) * t107 - pkin(10) * t108 + t73;
t220 = t193 * t31 + t196 * t22 + t54 * t263 - t50 * t264;
t219 = (-t193 * t265 + t117) * t125;
t32 = -pkin(4) * t325 - t42;
t218 = -pkin(10) * t84 + t125 * t32;
t23 = -pkin(4) * t155 - t225;
t6 = -qJD(5) * t15 - t16 * t193 + t196 * t28;
t203 = -t306 * qJD(5) - t193 * t22 + t196 * t31;
t202 = -t326 * t194 + t197 * t205;
t7 = t36 * pkin(5) + t17;
t179 = t309 * t196;
t178 = t309 * t193;
t163 = t196 * t169;
t126 = -t193 * t290 + t270;
t121 = -qJ(6) * t276 + t163 + (-pkin(5) - t285) * t197;
t83 = t85 ^ 2;
t48 = -qJD(5) * t114 + t108 * t196 + t155 * t193;
t47 = qJD(5) * t115 + t108 * t193 - t155 * t196;
t25 = t85 * pkin(5) + qJD(6) + t32;
t18 = -qJ(6) * t114 + t306;
t11 = pkin(5) * t137 - qJ(6) * t115 + t244;
t4 = -qJ(6) * t47 - qJD(6) * t114 + t220;
t3 = pkin(5) * t107 - qJ(6) * t48 - qJD(6) * t115 + t203;
t2 = -qJ(6) * t36 - qJD(6) * t85 + t5;
t1 = pkin(5) * t84 + qJ(6) * t35 - qJD(6) * t87 + t6;
t12 = [0, 0, 0, 0.2e1 * t195 * t233, -0.2e1 * t269 * t246, t236 * t198 * t266, -t236 * t250, 0, t217 * qJD(2) ^ 2 + 0.2e1 * t319 * t258, -0.2e1 * pkin(1) * t233 - (-pkin(8) * t250 + t182) * t237 - (-pkin(8) * t232 + t177) * t192, -t109 * t151 - t110 * t206 - t80 * t155 - t79 * t156 + t73 * t157 - t62 * t159 + t61 * t160 - t223 * t74, -t61 * t109 + t62 * t110 - t79 * t73 + t80 * t74 + (t166 + t224) * t324, t129 * t108 + t138 * t201, -t129 * t107 - t108 * t127 - t137 * t201 - t138 * t84, t108 * t325 + t129 * t155 + t138 * t205 + t159 * t201, -t107 * t259 - t84 * t159 - t127 * t155 + (-t107 * t282 - t137 * t215) * t268, t259 * t155 + (t155 * t282 + t159 * t215) * t268, t71 * t107 + t73 * t127 + t61 * t137 + t42 * t155 - t159 * t241 + t205 * t239 + t225 * t325 + t99 * t84, t71 * t108 + t73 * t129 + t61 * t138 - t43 * t155 + t159 * t221 + t201 * t99 - t206 * t293 - t216 * t325, -t115 * t35 + t48 * t87, t114 * t35 - t115 * t36 - t47 * t87 - t48 * t85, t107 * t87 + t115 * t84 + t125 * t48 - t137 * t35, -t107 * t85 - t114 * t84 - t125 * t47 - t137 * t36, t107 * t125 + t137 * t84, t14 * t107 + t17 * t114 + t125 * t203 + t6 * t137 + t23 * t85 + t244 * t84 + t32 * t47 + t49 * t36, -t15 * t107 + t17 * t115 - t220 * t125 - t5 * t137 + t23 * t87 - t306 * t84 + t32 * t48 - t49 * t35, -t1 * t115 - t10 * t47 + t11 * t35 - t114 * t2 - t18 * t36 - t3 * t87 - t4 * t85 - t48 * t8, t2 * t18 + t10 * t4 + t1 * t11 + t8 * t3 + t7 * (pkin(5) * t114 + t49) + t25 * (pkin(5) * t47 + t23); 0, 0, 0, -t195 * t253, t269 * t283, -t198 * t249, t195 * t249, 0, -t319 * t199, pkin(1) * t253 + (-pkin(8) * t251 + t181) * t267 (t80 - t101) * t157 - (-t102 + t79) * t223 + (-t191 * t151 - t189 * t206) * pkin(2), t101 * t79 - t102 * t80 + (-t166 * t251 + t189 * t62 - t191 * t61) * pkin(2), -qJD(4) * t194 ^ 2 * t157 + ((qJD(4) * t237 + t151) * t194 + t325 * t129) * t197, -t194 * t84 + t197 * t201 + (-t261 - t287) * t129 - t318 * t127, t194 * t205 - t288 + t330, t202 + t286, -t325 * t157, -t316 * t284 + t185 * t84 - t61 * t197 - t240 * t325 - t42 * t157 - t101 * t127 + (-t184 * t205 + t325 * t71) * t194, -t101 * t129 + t43 * t157 + t185 * t201 - t206 * t284 + t325 * t292 + t318 * t71 + (t184 * t316 + t61) * t194, t211 * t87 - t276 * t35, -t323 + (-t254 + (t35 + t289) * t194) * t193 + t235 (t298 + t302) * t194 + t320 + t243, t297 + t219 + (t222 - t322) * t194, t125 * t194 * t325 - t197 * t84, -t32 * t117 + t163 * t84 - t44 * t85 + ((-qJD(5) * t169 + t45) * t193 + t321) * t125 + (t32 * t193 * qJD(4) - t6 + (qJD(4) * t85 + t222) * t184) * t197 + (t14 * t325 + t184 * t36 + t263 * t32 + t301) * t194, -t270 * t84 - t44 * t87 - t32 * t118 + t317 * t125 + (t184 * t125 * t264 + t5 + (t184 * t87 + t32 * t196) * qJD(4)) * t197 + (-t32 * t264 - t15 * t223 + t300 - t184 * t35 + (t184 * t238 - t15) * qJD(4)) * t194, t10 * t117 + t118 * t8 + t121 * t35 - t126 * t36 - t305 * t87 - t304 * t85 + t227 * t265 + (-qJD(5) * t228 - t1 * t196 - t193 * t2) * t194, t245 * t7 * t194 + t1 * t121 + t304 * t10 + t2 * t126 + t305 * t8 + (t229 * pkin(5) + t245 * t265 - t44) * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t157 ^ 2 - t223 ^ 2, t157 * t79 + t223 * t80 + t176, 0, 0, 0, 0, 0, t202 - t286, -t194 * t206 - t288 - t330, 0, 0, 0, 0, 0, -t297 + t219 + (t222 + t322) * t194 (-t298 + t302) * t194 - t320 + t243, t323 + (t254 + (-t35 + t289) * t194) * t193 + t235, -t10 * t118 + t117 * t8 + (qJD(4) * t228 - t7) * t197 + (qJD(5) * t227 - t1 * t193 + t196 * t2 + t25 * t325) * t194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129 * t127, -t127 ^ 2 + t129 ^ 2, t127 * t223 + t274, t129 * t223 - t277, t205, -t71 * t129 + t325 * t43 - t241, t71 * t127 + t325 * t42 + t221, -t193 * t35 + t238 * t87 (-t35 - t303) * t196 + (-t36 - t296) * t193, t125 * t238 - t129 * t87 + t299, -t125 ^ 2 * t193 + t129 * t85 + t298, -t125 * t129, -pkin(4) * t36 - t14 * t129 - t300 - t43 * t85 + (-pkin(10) * t263 - t69) * t125 + (t42 * t125 + t218) * t193, pkin(4) * t35 + t15 * t129 + t301 - t43 * t87 + (pkin(10) * t264 + t308) * t125 + t218 * t196, t178 * t35 + t179 * t36 - t294 * t87 - t295 * t85 + (-t125 * t8 + t2) * t196 + (-t10 * t125 - t1) * t193, -t2 * t179 + t1 * t178 + t7 * (-pkin(5) * t196 - pkin(4)) + t294 * t8 + (t125 * t311 - t43) * t25 + t295 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87 * t85, -t83 + t315, -t35 + t303, -t36 + t296, t84, t125 * t15 - t32 * t87 + t6, t125 * t14 + t32 * t85 - t5, pkin(5) * t35 - t314 * t85, t314 * t10 + (-t25 * t87 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83 - t315, t10 * t85 + t8 * t87 + t7;];
tauc_reg  = t12;
