% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% tauc_reg [6x35]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRRR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:14:33
% EndTime: 2019-03-09 07:14:51
% DurationCPUTime: 6.41s
% Computational Cost: add. (9107->436), mult. (23570->589), div. (0->0), fcn. (19062->10), ass. (0->221)
t223 = cos(pkin(11));
t231 = cos(qJ(3));
t290 = t231 * t223;
t212 = qJD(1) * t290;
t222 = sin(pkin(11));
t227 = sin(qJ(3));
t297 = t222 * t227;
t270 = qJD(1) * t297;
t179 = t212 - t270;
t329 = qJD(4) + qJD(5);
t354 = t179 - t329;
t224 = sin(qJ(6));
t228 = cos(qJ(6));
t192 = t222 * t231 + t223 * t227;
t180 = t192 * qJD(1);
t226 = sin(qJ(4));
t230 = cos(qJ(4));
t275 = t230 * qJD(3);
t157 = t180 * t226 - t275;
t159 = qJD(3) * t226 + t180 * t230;
t225 = sin(qJ(5));
t229 = cos(qJ(5));
t249 = t157 * t225 - t229 * t159;
t95 = t229 * t157 + t159 * t225;
t251 = t224 * t95 + t228 * t249;
t53 = t224 * t249 - t228 * t95;
t346 = t251 * t53;
t195 = t225 * t226 - t229 * t230;
t287 = t354 * t195;
t292 = t225 * t230;
t196 = t226 * t229 + t292;
t286 = t354 * t196;
t341 = t251 ^ 2 - t53 ^ 2;
t174 = qJD(4) - t179;
t214 = -pkin(2) * t223 - pkin(1);
t202 = t214 * qJD(1) + qJD(2);
t113 = -pkin(3) * t179 - pkin(8) * t180 + t202;
t325 = pkin(7) + qJ(2);
t203 = t325 * t222;
t193 = qJD(1) * t203;
t204 = t325 * t223;
t194 = qJD(1) * t204;
t144 = -t227 * t193 + t231 * t194;
t138 = qJD(3) * pkin(8) + t144;
t80 = t230 * t113 - t138 * t226;
t62 = -pkin(9) * t159 + t80;
t43 = pkin(4) * t174 + t62;
t81 = t113 * t226 + t138 * t230;
t63 = -pkin(9) * t157 + t81;
t59 = t229 * t63;
t25 = t225 * t43 + t59;
t349 = pkin(10) * t95;
t15 = t25 - t349;
t277 = qJD(6) * t224;
t13 = t15 * t277;
t331 = -t193 * t231 - t227 * t194;
t137 = -qJD(3) * pkin(3) - t331;
t88 = pkin(4) * t157 + t137;
t55 = pkin(5) * t95 + t88;
t340 = -t55 * t53 + t13;
t168 = qJD(5) + t174;
t164 = qJD(6) + t168;
t209 = qJD(3) * t212;
t169 = -qJD(3) * t270 + t209;
t281 = qJD(4) * t226;
t107 = qJD(4) * t275 + t230 * t169 - t180 * t281;
t108 = qJD(4) * t159 + t169 * t226;
t233 = qJD(5) * t249 - t107 * t225 - t229 * t108;
t276 = qJD(6) * t228;
t278 = qJD(5) * t229;
t279 = qJD(5) * t225;
t39 = t229 * t107 - t225 * t108 - t157 * t278 - t159 * t279;
t8 = t224 * t233 + t228 * t39 + t249 * t277 - t95 * t276;
t339 = -t164 * t53 + t8;
t182 = t192 * qJD(3);
t170 = qJD(1) * t182;
t122 = pkin(3) * t170 - pkin(8) * t169;
t118 = t230 * t122;
t191 = -t290 + t297;
t238 = t191 * qJD(2);
t99 = -qJD(1) * t238 + qJD(3) * t331;
t234 = -qJD(4) * t81 - t226 * t99 + t118;
t22 = pkin(4) * t170 - pkin(9) * t107 + t234;
t280 = qJD(4) * t230;
t241 = t113 * t280 + t226 * t122 - t138 * t281 + t230 * t99;
t28 = -pkin(9) * t108 + t241;
t266 = t229 * t22 - t225 * t28;
t236 = -qJD(5) * t25 + t266;
t2 = pkin(5) * t170 - pkin(10) * t39 + t236;
t261 = -t225 * t22 - t229 * t28 - t43 * t278 + t63 * t279;
t3 = pkin(10) * t233 - t261;
t272 = t228 * t2 - t224 * t3;
t353 = t55 * t251 + t272;
t235 = qJD(6) * t251 - t224 * t39 + t228 * t233;
t334 = -t164 * t251 + t235;
t326 = pkin(8) + pkin(9);
t271 = qJD(4) * t326;
t139 = pkin(3) * t180 - pkin(8) * t179;
t289 = t226 * t139 + t230 * t331;
t304 = t179 * t226;
t352 = -pkin(9) * t304 + t226 * t271 + t289;
t128 = t230 * t139;
t351 = pkin(4) * t180 - t226 * t331 + t128 + (-pkin(9) * t179 + t271) * t230;
t350 = t281 - t304;
t348 = pkin(10) * t249;
t147 = -t195 * t224 + t196 * t228;
t322 = qJD(6) * t147 + t287 * t224 - t286 * t228;
t347 = t249 * t95;
t344 = t192 * qJD(2);
t268 = t192 * t280;
t181 = t191 * qJD(3);
t303 = t181 * t226;
t343 = t268 - t303;
t342 = t249 ^ 2 - t95 ^ 2;
t338 = t168 * t95 + t39;
t337 = t88 * t95 + t261;
t336 = t88 * t249 + t236;
t57 = t225 * t63;
t24 = t229 * t43 - t57;
t14 = t24 + t348;
t12 = pkin(5) * t168 + t14;
t315 = t228 * t15;
t7 = t224 * t12 + t315;
t335 = -qJD(6) * t7 + t353;
t333 = -t168 * t249 + t233;
t332 = t351 * t229;
t130 = t196 * t192;
t255 = t350 * pkin(4) - t144;
t155 = t203 * t231 + t227 * t204;
t205 = t326 * t226;
t206 = t326 * t230;
t285 = -t225 * t205 + t229 * t206;
t330 = t205 * t278 + t206 * t279 + t351 * t225 + t352 * t229;
t146 = t228 * t195 + t196 * t224;
t323 = -qJD(6) * t146 + t286 * t224 + t287 * t228;
t328 = -t147 * t170 - t323 * t164;
t327 = -t287 * t168 - t170 * t196;
t324 = t229 * t62 - t57;
t141 = pkin(3) * t191 - pkin(8) * t192 + t214;
t133 = t230 * t141;
t156 = -t203 * t227 + t204 * t231;
t300 = t192 * t230;
t69 = pkin(4) * t191 - pkin(9) * t300 - t156 * t226 + t133;
t149 = t230 * t156;
t288 = t226 * t141 + t149;
t301 = t192 * t226;
t82 = -pkin(9) * t301 + t288;
t320 = t225 * t69 + t229 * t82;
t319 = t180 * t53;
t318 = t180 * t251;
t317 = t180 * t95;
t316 = t180 * t249;
t314 = -t286 * pkin(5) + t255;
t313 = t107 * t226;
t311 = t157 * t174;
t310 = t157 * t180;
t309 = t159 * t174;
t308 = t159 * t180;
t305 = t170 * t226;
t302 = t181 * t230;
t294 = t225 * t170;
t293 = t225 * t228;
t291 = t228 * t170;
t161 = t230 * t170;
t284 = t222 ^ 2 + t223 ^ 2;
t283 = qJD(3) * t227;
t282 = qJD(3) * t231;
t274 = qJD(1) * qJD(2);
t218 = -pkin(4) * t230 - pkin(3);
t269 = t192 * t281;
t267 = qJD(6) * t12 + t3;
t114 = -t155 * qJD(3) - t238;
t140 = pkin(3) * t182 + pkin(8) * t181;
t129 = t230 * t140;
t32 = pkin(9) * t302 + pkin(4) * t182 - t114 * t226 + t129 + (-t149 + (pkin(9) * t192 - t141) * t226) * qJD(4);
t240 = t230 * t114 + t226 * t140 + t141 * t280 - t156 * t281;
t35 = -pkin(9) * t343 + t240;
t264 = -t225 * t35 + t229 * t32;
t263 = -t225 * t62 - t59;
t262 = -t225 * t82 + t229 * t69;
t260 = t284 * qJD(1) ^ 2;
t258 = -t229 * t205 - t206 * t225;
t257 = t174 * t230;
t100 = t344 * qJD(1) - t193 * t283 + t194 * t282;
t115 = -t203 * t283 + t204 * t282 + t344;
t256 = -t146 * t170 - t322 * t164;
t124 = -pkin(10) * t196 + t258;
t254 = -t286 * pkin(10) - qJD(6) * t124 + t330;
t125 = -pkin(10) * t195 + t285;
t253 = pkin(5) * t180 + t287 * pkin(10) + t285 * qJD(5) + qJD(6) * t125 - t225 * t352 + t332;
t252 = t286 * t168 - t195 * t170;
t116 = pkin(4) * t301 + t155;
t131 = t195 * t192;
t84 = t228 * t130 - t131 * t224;
t85 = -t130 * t224 - t131 * t228;
t247 = 0.2e1 * t284 * t274;
t246 = -t350 * t174 + t161;
t87 = t343 * pkin(4) + t115;
t245 = t225 * t32 + t229 * t35 + t69 * t278 - t82 * t279;
t243 = -t269 - t302;
t68 = pkin(4) * t108 + t100;
t242 = -pkin(8) * t170 + t174 * t137;
t217 = pkin(4) * t229 + pkin(5);
t167 = pkin(5) * t195 + t218;
t145 = t170 * t191;
t86 = pkin(4) * t159 - pkin(5) * t249;
t83 = pkin(5) * t130 + t116;
t46 = -t181 * t292 - t225 * t269 - t279 * t301 + (t329 * t300 - t303) * t229;
t45 = -t329 * t130 + t195 * t181;
t33 = pkin(5) * t46 + t87;
t29 = -pkin(10) * t130 + t320;
t26 = pkin(5) * t191 + pkin(10) * t131 + t262;
t19 = -pkin(5) * t233 + t68;
t17 = t324 + t348;
t16 = t263 + t349;
t11 = qJD(6) * t85 + t224 * t45 + t228 * t46;
t10 = -qJD(6) * t84 - t224 * t46 + t228 * t45;
t6 = t228 * t12 - t15 * t224;
t5 = -pkin(10) * t46 + t245;
t4 = pkin(5) * t182 - pkin(10) * t45 - qJD(5) * t320 + t264;
t1 = [0, 0, 0, 0, 0, t247, qJ(2) * t247, t169 * t192 - t180 * t181, -t169 * t191 - t170 * t192 - t179 * t181 - t180 * t182, -t181 * qJD(3), -t182 * qJD(3), 0, -qJD(3) * t115 + t170 * t214 + t182 * t202, -qJD(3) * t114 + t169 * t214 - t181 * t202, t107 * t300 + t159 * t243 -(-t157 * t230 - t159 * t226) * t181 + (-t313 - t108 * t230 + (t157 * t226 - t159 * t230) * qJD(4)) * t192, t107 * t191 + t159 * t182 + t192 * t161 + t174 * t243, -t108 * t191 - t157 * t182 - t170 * t301 - t174 * t343, t174 * t182 + t145 (-t156 * t280 + t129) * t174 + t133 * t170 + (-t138 * t280 + t118) * t191 + t80 * t182 + t115 * t157 + t155 * t108 + t137 * t268 + ((-qJD(4) * t141 - t114) * t174 - t156 * t170 + (-qJD(4) * t113 - t99) * t191 + t100 * t192 - t137 * t181) * t226, t100 * t300 + t155 * t107 + t115 * t159 + t243 * t137 - t288 * t170 - t240 * t174 - t81 * t182 - t241 * t191, -t131 * t39 - t249 * t45, -t130 * t39 - t131 * t233 + t249 * t46 - t45 * t95, -t131 * t170 + t168 * t45 - t182 * t249 + t191 * t39, -t130 * t170 - t168 * t46 - t182 * t95 + t191 * t233, t168 * t182 + t145, t264 * t168 + t262 * t170 + t266 * t191 + t24 * t182 + t87 * t95 - t116 * t233 + t68 * t130 + t88 * t46 + (-t168 * t320 - t191 * t25) * qJD(5), t116 * t39 - t68 * t131 - t245 * t168 - t320 * t170 - t25 * t182 + t261 * t191 - t249 * t87 + t88 * t45, -t10 * t251 + t8 * t85, t10 * t53 + t11 * t251 + t235 * t85 - t8 * t84, t10 * t164 + t170 * t85 - t182 * t251 + t191 * t8, -t11 * t164 - t170 * t84 + t182 * t53 + t191 * t235, t164 * t182 + t145 (-t224 * t5 + t228 * t4) * t164 + (-t224 * t29 + t228 * t26) * t170 + t272 * t191 + t6 * t182 - t33 * t53 - t83 * t235 + t19 * t84 + t55 * t11 + ((-t224 * t26 - t228 * t29) * t164 - t7 * t191) * qJD(6), t55 * t10 + t13 * t191 - t7 * t182 + t19 * t85 - t33 * t251 + t83 * t8 + (-(-qJD(6) * t29 + t4) * t164 - t26 * t170 - t2 * t191) * t224 + (-(qJD(6) * t26 + t5) * t164 - t29 * t170 - t267 * t191) * t228; 0, 0, 0, 0, 0, -t260, -qJ(2) * t260, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t180, t209 + (t179 - t270) * qJD(3), 0, 0, 0, 0, 0, t246 - t310, -t174 ^ 2 * t230 - t305 - t308, 0, 0, 0, 0, 0, t252 - t317, t316 + t327, 0, 0, 0, 0, 0, t256 + t319, t318 + t328; 0, 0, 0, 0, 0, 0, 0, -t180 * t179, -t179 ^ 2 + t180 ^ 2, t209 + (-t179 - t270) * qJD(3), 0, 0, qJD(3) * t144 - t180 * t202 - t100, -t179 * t202 + t191 * t274, t159 * t257 + t313 (t107 - t311) * t230 + (-t108 - t309) * t226, t174 * t257 + t305 - t308, t246 + t310, -t174 * t180, -pkin(3) * t108 - t100 * t230 - t144 * t157 - t80 * t180 + (-pkin(8) * t280 - t128) * t174 + (t174 * t331 + t242) * t226, -pkin(3) * t107 + t100 * t226 - t144 * t159 + t81 * t180 + (pkin(8) * t281 + t289) * t174 + t242 * t230, t196 * t39 - t249 * t287, -t195 * t39 + t196 * t233 - t249 * t286 - t287 * t95, t316 - t327, t252 + t317, -t168 * t180, t258 * t170 - t218 * t233 + t68 * t195 - t24 * t180 + t255 * t95 - t286 * t88 + (-t206 * t278 + (qJD(5) * t205 + t352) * t225 - t332) * t168, t330 * t168 - t285 * t170 + t25 * t180 + t68 * t196 + t218 * t39 - t249 * t255 + t287 * t88, t147 * t8 - t251 * t323, -t146 * t8 + t147 * t235 + t251 * t322 + t323 * t53, t318 - t328, t256 - t319, -t164 * t180 (t124 * t228 - t125 * t224) * t170 - t167 * t235 + t19 * t146 - t6 * t180 + t322 * t55 - t314 * t53 + (t224 * t254 - t228 * t253) * t164 -(t124 * t224 + t125 * t228) * t170 + t167 * t8 + t19 * t147 + t7 * t180 + t323 * t55 - t314 * t251 + (t224 * t253 + t228 * t254) * t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159 * t157, -t157 ^ 2 + t159 ^ 2, t107 + t311, -t108 + t309, t170, -t137 * t159 + t174 * t81 + t234, t137 * t157 + t174 * t80 - t241, -t347, t342, t338, t333, t170, -t263 * t168 + (-t159 * t95 - t168 * t279 + t170 * t229) * pkin(4) + t336, t324 * t168 + (t159 * t249 - t168 * t278 - t294) * pkin(4) + t337, t346, t341, t339, t334, t170, t217 * t291 - (t16 * t228 - t17 * t224) * t164 + t86 * t53 + (-t224 * t294 + (-t224 * t229 - t293) * t164 * qJD(5)) * pkin(4) + ((-pkin(4) * t293 - t217 * t224) * t164 - t7) * qJD(6) + t353, t86 * t251 + (-t217 * t170 - t2 + (t16 - (-qJD(5) - qJD(6)) * t225 * pkin(4)) * t164) * t224 + (-pkin(4) * t294 + (-pkin(4) * t278 - qJD(6) * t217 + t17) * t164 - t267) * t228 + t340; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t347, t342, t338, t333, t170, t168 * t25 + t336, t168 * t24 + t337, t346, t341, t339, t334, t170 -(-t224 * t14 - t315) * t164 + (-t164 * t277 - t249 * t53 + t291) * pkin(5) + t335 (-t15 * t164 - t2) * t224 + (t14 * t164 - t267) * t228 + (-t164 * t276 - t224 * t170 - t249 * t251) * pkin(5) + t340; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t346, t341, t339, t334, t170, t164 * t7 + t335, t164 * t6 - t224 * t2 - t228 * t267 + t340;];
tauc_reg  = t1;
