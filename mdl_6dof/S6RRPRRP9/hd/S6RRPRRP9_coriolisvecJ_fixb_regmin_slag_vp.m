% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPRRP9
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
% tauc_reg [6x30]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRRP9_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP9_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:33:58
% EndTime: 2019-03-09 12:34:17
% DurationCPUTime: 6.69s
% Computational Cost: add. (10585->484), mult. (28071->684), div. (0->0), fcn. (22717->10), ass. (0->242)
t233 = cos(pkin(11));
t239 = cos(qJ(4));
t319 = t239 * t233;
t231 = sin(pkin(11));
t236 = sin(qJ(4));
t325 = t231 * t236;
t195 = -t319 + t325;
t232 = sin(pkin(6));
t240 = cos(qJ(2));
t322 = t232 * t240;
t249 = t195 * t322;
t161 = qJD(1) * t249;
t191 = t195 * qJD(4);
t354 = t191 - t161;
t196 = t231 * t239 + t233 * t236;
t250 = t196 * t322;
t314 = -qJD(1) * t250 + t196 * qJD(4);
t237 = sin(qJ(2));
t267 = pkin(2) * t237 - qJ(3) * t240;
t312 = qJD(1) * t232;
t183 = t267 * t312;
t310 = qJD(1) * t237;
t293 = t232 * t310;
t234 = cos(pkin(6));
t311 = qJD(1) * t234;
t299 = pkin(1) * t311;
t184 = -pkin(8) * t293 + t240 * t299;
t129 = t233 * t183 - t231 * t184;
t321 = t233 * t240;
t253 = t232 * (pkin(3) * t237 - pkin(9) * t321);
t107 = qJD(1) * t253 + t129;
t130 = t231 * t183 + t233 * t184;
t309 = qJD(1) * t240;
t221 = t232 * t309;
t273 = t231 * t221;
t114 = -pkin(9) * t273 + t130;
t348 = pkin(9) + qJ(3);
t208 = t348 * t231;
t209 = t348 * t233;
t264 = -t208 * t239 - t209 * t236;
t370 = -t195 * qJD(3) + t264 * qJD(4) - t236 * t107 - t239 * t114;
t185 = pkin(8) * t221 + t237 * t299;
t153 = pkin(3) * t273 + t185;
t369 = t314 * pkin(4) + pkin(10) * t354 - t153;
t364 = pkin(10) * t293 - t370;
t235 = sin(qJ(5));
t238 = cos(qJ(5));
t222 = qJD(2) + t311;
t172 = t222 * t233 - t231 * t293;
t173 = t222 * t231 + t233 * t293;
t265 = t172 * t236 + t173 * t239;
t301 = t221 - qJD(4);
t100 = -t235 * t301 + t238 * t265;
t120 = -t239 * t172 + t173 * t236;
t360 = qJD(5) + t120;
t330 = t360 * t235;
t368 = t100 * t330;
t300 = qJD(1) * qJD(2);
t288 = t232 * t300;
t271 = t240 * t288;
t263 = t231 * t271;
t306 = qJD(4) * t239;
t316 = t172 * t306 + t271 * t319;
t243 = (-qJD(4) * t173 - t263) * t236 + t316;
t270 = t237 * t288;
t242 = t235 * t243 - t238 * t270;
t41 = qJD(5) * t100 + t242;
t307 = qJD(4) * t236;
t162 = qJ(3) * t222 + t185;
t180 = (-pkin(2) * t240 - qJ(3) * t237 - pkin(1)) * t232;
t167 = qJD(1) * t180;
t108 = -t162 * t231 + t233 * t167;
t66 = -pkin(3) * t221 - pkin(9) * t173 + t108;
t165 = (t267 * qJD(2) - qJD(3) * t237) * t232;
t149 = qJD(1) * t165;
t298 = pkin(1) * qJD(2) * t234;
t274 = qJD(1) * t298;
t254 = -pkin(8) * t270 + t240 * t274;
t150 = qJD(3) * t222 + t254;
t103 = t233 * t149 - t231 * t150;
t248 = qJD(2) * t253;
t67 = qJD(1) * t248 + t103;
t109 = t233 * t162 + t231 * t167;
t80 = pkin(9) * t172 + t109;
t104 = t231 * t149 + t233 * t150;
t81 = -pkin(9) * t263 + t104;
t258 = -t236 * t67 - t239 * t81 - t66 * t306 + t80 * t307;
t14 = pkin(10) * t270 - t258;
t178 = pkin(8) * t271 + t237 * t274;
t145 = pkin(3) * t263 + t178;
t247 = qJD(2) * t250;
t246 = qJD(1) * t247;
t352 = t265 * qJD(4);
t79 = t246 + t352;
t30 = t79 * pkin(4) - pkin(10) * t243 + t145;
t304 = qJD(5) * t238;
t305 = qJD(5) * t235;
t35 = t236 * t66 + t239 * t80;
t32 = -t301 * pkin(10) + t35;
t155 = -pkin(2) * t222 + qJD(3) - t184;
t118 = -pkin(3) * t172 + t155;
t44 = pkin(4) * t120 - pkin(10) * t265 + t118;
t5 = t238 * t14 + t235 * t30 + t44 * t304 - t32 * t305;
t280 = t238 * t301;
t98 = t235 * t265 + t280;
t2 = -qJ(6) * t41 - qJD(6) * t98 + t5;
t16 = -t235 * t32 + t238 * t44;
t9 = -qJ(6) * t100 + t16;
t8 = pkin(5) * t360 + t9;
t367 = -t360 * t8 + t2;
t40 = qJD(5) * t280 - t235 * t270 - t238 * t243 + t265 * t305;
t17 = t235 * t44 + t238 * t32;
t6 = -t17 * qJD(5) - t235 * t14 + t238 * t30;
t1 = t79 * pkin(5) + t40 * qJ(6) - t100 * qJD(6) + t6;
t10 = -qJ(6) * t98 + t17;
t366 = t10 * t360 + t1;
t365 = t120 * t98;
t363 = t120 * t301;
t131 = -t161 * t235 - t238 * t293;
t362 = -t196 * t304 + t131;
t347 = -qJ(6) - pkin(10);
t361 = -qJ(6) * t120 + qJD(5) * t347;
t359 = t265 * t98;
t358 = t369 * t238;
t357 = t100 * t265;
t356 = t221 * t265;
t159 = -t208 * t236 + t209 * t239;
t355 = t196 * qJD(3) + t159 * qJD(4) + t239 * t107;
t227 = -pkin(3) * t233 - pkin(2);
t142 = pkin(4) * t195 - pkin(10) * t196 + t227;
t353 = t142 * t304 + t369 * t235 - t364 * t238;
t351 = t100 ^ 2;
t230 = t240 ^ 2;
t350 = t8 - t9;
t349 = pkin(1) * t237;
t132 = -t161 * t238 + t235 * t293;
t147 = t238 * t159;
t262 = qJ(6) * t191 - qJD(6) * t196;
t346 = qJ(6) * t132 - t147 * qJD(5) + t262 * t238 + t358 + ((qJ(6) * t196 - t142) * qJD(5) + t364) * t235 + t314 * pkin(5);
t345 = (-qJD(5) * t159 + t262) * t235 + t353 + t362 * qJ(6);
t34 = -t236 * t80 + t239 * t66;
t59 = pkin(4) * t265 + pkin(10) * t120;
t344 = t235 * t59 + t238 * t34;
t343 = -t235 * t41 - t98 * t304;
t179 = pkin(8) * t322 + (qJ(3) + t349) * t234;
t126 = t233 * t179 + t231 * t180;
t323 = t232 * t237;
t189 = t231 * t323 - t234 * t233;
t106 = -pkin(9) * t189 + t126;
t125 = -t179 * t231 + t233 * t180;
t190 = t231 * t234 + t233 * t323;
t88 = -pkin(3) * t322 - pkin(9) * t190 + t125;
t340 = t239 * t106 + t236 * t88;
t47 = -pkin(10) * t322 + t340;
t133 = t239 * t189 + t190 * t236;
t134 = -t189 * t236 + t190 * t239;
t182 = pkin(8) * t323 + (-pkin(1) * t240 - pkin(2)) * t234;
t140 = t189 * pkin(3) + t182;
t56 = t133 * pkin(4) - t134 * pkin(10) + t140;
t342 = t235 * t56 + t238 * t47;
t339 = t235 * t79;
t338 = t35 * t237;
t337 = t40 * t235;
t110 = t236 * t114;
t261 = pkin(4) * t293 - t110;
t336 = t261 + t355;
t335 = t238 * qJD(6) + t235 * t361 - t344;
t58 = t238 * t59;
t334 = -pkin(5) * t265 - t58 + t361 * t238 + (-qJD(6) + t34) * t235;
t332 = t108 * t237;
t331 = t109 * t237;
t329 = t191 * t235;
t328 = t196 * t235;
t327 = t196 * t238;
t228 = t232 ^ 2;
t326 = t228 * qJD(1) ^ 2;
t324 = t231 * t240;
t317 = t235 * t142 + t147;
t308 = qJD(2) * t237;
t292 = t232 * t308;
t260 = -pkin(8) * t292 + t240 * t298;
t171 = qJD(3) * t234 + t260;
t113 = t231 * t165 + t233 * t171;
t291 = qJD(2) * t322;
t186 = pkin(8) * t291 + t237 * t298;
t313 = t237 ^ 2 - t230;
t303 = qJD(2) - t222;
t297 = t240 * t326;
t296 = t236 * t324;
t295 = t235 * t322;
t272 = t231 * t291;
t154 = pkin(3) * t272 + t186;
t289 = t228 * t300;
t287 = -t235 * t47 + t238 * t56;
t286 = -t236 * t106 + t239 * t88;
t284 = -t236 * t81 + t239 * t67 - t80 * t306 - t66 * t307;
t283 = t131 + t329;
t282 = t191 * t238 + t132;
t112 = t233 * t165 - t231 * t171;
t281 = t360 ^ 2;
t278 = t360 * t238;
t277 = t301 * t232;
t276 = t222 + t311;
t275 = 0.2e1 * t289;
t46 = pkin(4) * t322 - t286;
t269 = -0.2e1 * pkin(1) * t289;
t102 = -pkin(9) * t272 + t113;
t85 = t112 + t248;
t266 = -t236 * t102 - t106 * t306 + t239 * t85 - t88 * t307;
t259 = -t173 * t307 + t316;
t115 = t134 * t235 + t238 * t322;
t255 = t239 * t102 - t106 * t307 + t236 * t85 + t88 * t306;
t21 = pkin(10) * t292 + t255;
t93 = -qJD(2) * t249 - t133 * qJD(4);
t94 = t134 * qJD(4) + t247;
t39 = pkin(4) * t94 - pkin(10) * t93 + t154;
t257 = t238 * t21 + t235 * t39 + t56 * t304 - t47 * t305;
t31 = t301 * pkin(4) - t34;
t256 = -pkin(10) * t79 + t31 * t360;
t252 = -t362 - t329;
t251 = -t196 * t305 - t282;
t15 = -pkin(4) * t270 - t284;
t22 = -pkin(4) * t292 - t266;
t245 = -qJ(3) * t308 + (-pkin(2) * qJD(2) + qJD(3) - t155) * t240;
t244 = -t342 * qJD(5) - t235 * t21 + t238 * t39;
t7 = pkin(5) * t41 + t15;
t214 = t347 * t238;
t213 = t347 * t235;
t139 = t238 * t142;
t116 = t134 * t238 - t295;
t95 = t98 ^ 2;
t73 = t238 * t79;
t70 = -qJ(6) * t328 + t317;
t60 = pkin(5) * t195 - qJ(6) * t327 - t159 * t235 + t139;
t49 = -qJD(5) * t295 + t134 * t304 + t235 * t93 - t238 * t292;
t48 = qJD(5) * t115 - t235 * t292 - t238 * t93;
t25 = t98 * pkin(5) + qJD(6) + t31;
t19 = -qJ(6) * t115 + t342;
t12 = pkin(5) * t133 - qJ(6) * t116 + t287;
t4 = -qJ(6) * t49 - qJD(6) * t115 + t257;
t3 = t94 * pkin(5) + t48 * qJ(6) - t116 * qJD(6) + t244;
t11 = [0, 0, 0, t237 * t240 * t275, -t313 * t275, t276 * t291, -t276 * t292, 0, -t178 * t234 - t186 * t222 + t237 * t269, -t222 * t260 - t234 * t254 + t240 * t269, -t186 * t172 + t178 * t189 + ((-qJD(1) * t112 - t103) * t240 + (t155 * t324 + t332 + (t125 * t237 + t182 * t324) * qJD(1)) * qJD(2)) * t232, t186 * t173 + t178 * t190 + ((qJD(1) * t113 + t104) * t240 + (t155 * t321 - t331 + (-t126 * t237 + t182 * t321) * qJD(1)) * qJD(2)) * t232, -t103 * t190 - t104 * t189 - t112 * t173 + t113 * t172 + (-t108 * t233 - t109 * t231 + (-t125 * t233 - t126 * t231) * qJD(1)) * t291, t103 * t125 + t104 * t126 + t108 * t112 + t109 * t113 + t155 * t186 + t178 * t182, t134 * t243 + t265 * t93, -t93 * t120 - t133 * t243 - t134 * t79 - t265 * t94, -t93 * t301 + (-t259 * t240 + (t265 * t237 + (t232 * t230 * t325 + t134 * t237) * qJD(1)) * qJD(2)) * t232, t94 * t301 + (t79 * t240 + (-qJD(1) * t133 - t120) * t308) * t232 (-t228 * t309 - t277) * t308, -t266 * t301 + t154 * t120 + t140 * t79 + t145 * t133 + t118 * t94 + (-t284 * t240 + (qJD(1) * t286 + t34) * t308) * t232, t255 * t301 + t154 * t265 + t140 * t259 + t145 * t134 + t118 * t93 + (-t258 * t240 + (-t338 + (-t140 * t296 - t237 * t340) * qJD(1)) * qJD(2)) * t232, -t100 * t48 - t116 * t40, -t100 * t49 + t115 * t40 - t116 * t41 + t48 * t98, t100 * t94 + t116 * t79 - t133 * t40 - t360 * t48, -t115 * t79 - t133 * t41 - t360 * t49 - t94 * t98, t133 * t79 + t360 * t94, t15 * t115 + t6 * t133 + t16 * t94 + t22 * t98 + t244 * t360 + t287 * t79 + t31 * t49 + t46 * t41, t22 * t100 + t15 * t116 - t5 * t133 - t17 * t94 - t257 * t360 - t31 * t48 - t342 * t79 - t46 * t40, -t1 * t116 - t10 * t49 - t100 * t3 - t115 * t2 + t12 * t40 - t19 * t41 - t4 * t98 + t48 * t8, t2 * t19 + t10 * t4 + t1 * t12 + t8 * t3 + t7 * (pkin(5) * t115 + t46) + t25 * (pkin(5) * t49 + t22); 0, 0, 0, -t237 * t297, t313 * t326, t303 * t221, -t303 * t293, 0, t185 * t222 + t326 * t349 - t178, pkin(1) * t297 + t184 * t222 - t254, t185 * t172 - t178 * t233 + (t129 * t240 + t231 * t245 - t332) * t312, -t185 * t173 + t178 * t231 + (-t130 * t240 + t233 * t245 + t331) * t312, t129 * t173 - t130 * t172 + (qJD(3) * t172 + t108 * t221 + t104) * t233 + (qJD(3) * t173 + t109 * t221 - t103) * t231, -t178 * pkin(2) - t108 * t129 - t109 * t130 - t155 * t185 + (-t108 * t231 + t109 * t233) * qJD(3) + (-t103 * t231 + t104 * t233) * qJ(3), t196 * t243 - t265 * t354, t120 * t354 - t195 * t243 - t196 * t79 - t265 * t314 (qJD(2) * t196 - t265) * t293 + t354 * t301 (-qJD(2) * t195 + t120) * t293 + t314 * t301, t277 * t310, t227 * t79 + t145 * t195 - t153 * t120 + t314 * t118 + (qJD(2) * t264 - t34) * t293 + (-t110 + t355) * t301, t227 * t259 + t145 * t196 - t153 * t265 - t354 * t118 + (t338 + (-t159 * t237 - t227 * t296) * qJD(2)) * t312 + t370 * t301, t100 * t251 - t327 * t40, t282 * t98 + t283 * t100 + (t337 - t238 * t41 + (-t100 * t238 + t235 * t98) * qJD(5)) * t196, t100 * t314 - t40 * t195 + t251 * t360 + t327 * t79, -t41 * t195 - t252 * t360 - t314 * t98 - t328 * t79, t79 * t195 + t314 * t360, t139 * t79 - t264 * t41 + t6 * t195 + t336 * t98 - t362 * t31 + t314 * t16 + (-t159 * t304 + t358) * t360 + (t15 * t196 - t159 * t79 - t31 * t191 + (-qJD(5) * t142 + t364) * t360) * t235, -t317 * t79 - t5 * t195 + t264 * t40 + t15 * t327 - t314 * t17 + (t159 * t305 - t353) * t360 + t336 * t100 + t251 * t31, t60 * t40 - t70 * t41 - t345 * t98 + t282 * t8 - t346 * t100 + t283 * t10 + (-t1 * t238 - t2 * t235 + (-t10 * t238 + t235 * t8) * qJD(5)) * t196, t2 * t70 + t1 * t60 + t7 * (pkin(5) * t328 - t264) + t346 * t8 + ((qJD(3) * t233 - qJD(4) * t208) * t236 + (qJD(3) * t231 + qJD(4) * t209 + t107) * t239 + t252 * pkin(5) + t261) * t25 + t345 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (qJD(2) * t231 - t173) * t221 (qJD(2) * t233 - t172) * t221, -t172 ^ 2 - t173 ^ 2, t108 * t173 - t109 * t172 + t178, 0, 0, 0, 0, 0, t246 + 0.2e1 * t352 - t356, t243 + t363, 0, 0, 0, 0, 0, -t235 * t281 - t359 + t73, -t238 * t281 - t339 - t357 (t40 - t365) * t238 + t368 + t343, t367 * t235 + t366 * t238 - t25 * t265; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t265 * t120, -t120 ^ 2 + t265 ^ 2, t243 - t363, -t196 * t271 - t356, t270, -t118 * t265 - t301 * t35 + t284, t118 * t120 - t301 * t34 + t258, t100 * t278 - t337 (-t40 - t365) * t238 - t368 + t343, t278 * t360 + t339 - t357, -t330 * t360 + t359 + t73, -t360 * t265, -pkin(4) * t41 - t16 * t265 - t15 * t238 - t35 * t98 + (-pkin(10) * t304 - t58) * t360 + (t34 * t360 + t256) * t235, pkin(4) * t40 - t35 * t100 + t17 * t265 + t15 * t235 + (pkin(10) * t305 + t344) * t360 + t256 * t238, -t334 * t100 + t213 * t40 + t214 * t41 - t366 * t235 + t367 * t238 - t335 * t98, -t2 * t214 + t1 * t213 + t7 * (-pkin(5) * t238 - pkin(4)) + t334 * t8 + (pkin(5) * t330 - t35) * t25 + t335 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100 * t98, -t95 + t351, t360 * t98 - t40, -t242 + (-qJD(5) + t360) * t100, t79, -t31 * t100 + t17 * t360 + t6, t16 * t360 + t31 * t98 - t5, pkin(5) * t40 - t350 * t98, t350 * t10 + (-t100 * t25 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95 - t351, t10 * t98 + t100 * t8 + t7;];
tauc_reg  = t11;
