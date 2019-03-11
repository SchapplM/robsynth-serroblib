% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauc_reg [6x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRPRP6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:01:09
% EndTime: 2019-03-09 17:01:28
% DurationCPUTime: 7.15s
% Computational Cost: add. (11127->466), mult. (29346->649), div. (0->0), fcn. (23150->10), ass. (0->238)
t235 = sin(pkin(6));
t243 = cos(qJ(2));
t313 = qJD(1) * t243;
t223 = t235 * t313;
t263 = t223 - qJD(3);
t239 = sin(qJ(3));
t240 = sin(qJ(2));
t242 = cos(qJ(3));
t314 = qJD(1) * t240;
t294 = t235 * t314;
t237 = cos(pkin(6));
t315 = qJD(1) * t237;
t302 = pkin(1) * t315;
t182 = -pkin(8) * t294 + t243 * t302;
t258 = (pkin(2) * t240 - pkin(9) * t243) * t235;
t183 = qJD(1) * t258;
t281 = -t239 * t182 + t242 * t183;
t349 = -qJ(4) - pkin(9);
t287 = qJD(3) * t349;
t316 = qJD(1) * t235;
t323 = t242 * t243;
t370 = (pkin(3) * t240 - qJ(4) * t323) * t316 + t281 + t239 * qJD(4) - t242 * t287;
t271 = t239 * t223;
t319 = t242 * t182 + t239 * t183;
t369 = qJ(4) * t271 + t242 * qJD(4) + t239 * t287 - t319;
t234 = sin(pkin(11));
t236 = cos(pkin(11));
t198 = t234 * t242 + t236 * t239;
t318 = t263 * t198;
t197 = t234 * t239 - t236 * t242;
t147 = t197 * t223;
t191 = t197 * qJD(3);
t368 = t191 - t147;
t336 = -t234 * t370 + t369 * t236;
t220 = t240 * t302;
t185 = pkin(8) * t223 + t220;
t310 = qJD(3) * t239;
t367 = -t185 + (-t271 + t310) * pkin(3);
t366 = -pkin(4) * t318 + t368 * pkin(10) + t367;
t365 = pkin(10) * t294 - t336;
t337 = t369 * t234 + t236 * t370;
t238 = sin(qJ(5));
t241 = cos(qJ(5));
t131 = -t147 * t238 - t241 * t294;
t307 = qJD(5) * t241;
t364 = -t198 * t307 + t131;
t303 = qJD(1) * qJD(2);
t289 = t235 * t303;
t270 = t243 * t289;
t176 = pkin(8) * t270 + qJD(2) * t220;
t273 = qJD(3) * t263;
t363 = -pkin(9) * t273 + t176;
t224 = qJD(2) + t315;
t165 = -t242 * t224 + t239 * t294;
t167 = t224 * t239 + t242 * t294;
t284 = -t236 * t165 - t167 * t234;
t354 = qJD(5) - t284;
t276 = t354 * t241;
t269 = qJD(3) * t294;
t309 = qJD(3) * t242;
t139 = t224 * t309 - t239 * t269 + t242 * t270;
t297 = t224 * t310 + t239 * t270 + t242 * t269;
t95 = t234 * t139 + t236 * t297;
t362 = -t238 * t95 - t354 * t276;
t360 = t366 * t241;
t261 = -t165 * t234 + t236 * t167;
t103 = -t238 * t263 + t241 * t261;
t277 = t354 * t238;
t359 = t103 * t277;
t357 = t165 * t263;
t356 = t167 * t263;
t338 = pkin(4) * t294 + t337;
t324 = t235 * t243;
t179 = pkin(8) * t324 + (pkin(1) * t240 + pkin(9)) * t237;
t180 = (-pkin(2) * t243 - pkin(9) * t240 - pkin(1)) * t235;
t320 = t242 * t179 + t239 * t180;
t295 = -pkin(3) * t242 - pkin(2);
t140 = pkin(4) * t197 - pkin(10) * t198 + t295;
t355 = t140 * t307 + t238 * t366 - t365 * t241;
t251 = t236 * t139 - t234 * t297;
t268 = t240 * t289;
t279 = t241 * t263;
t308 = qJD(5) * t238;
t51 = qJD(5) * t279 - t238 * t268 - t241 * t251 + t261 * t308;
t152 = pkin(9) * t224 + t185;
t159 = qJD(1) * t180;
t115 = t152 * t242 + t159 * t239;
t184 = qJD(2) * t258;
t174 = qJD(1) * t184;
t325 = t235 * t240;
t225 = pkin(8) * t325;
t350 = pkin(1) * t243;
t186 = (t237 * t350 - t225) * qJD(2);
t175 = qJD(1) * t186;
t283 = t242 * t174 - t239 * t175;
t246 = -qJD(3) * t115 + t283;
t43 = pkin(3) * t268 - t139 * qJ(4) - t167 * qJD(4) + t246;
t255 = t152 * t310 - t159 * t309 - t239 * t174 - t242 * t175;
t49 = -qJ(4) * t297 - t165 * qJD(4) - t255;
t18 = t234 * t43 + t236 * t49;
t15 = pkin(10) * t268 + t18;
t94 = -qJ(4) * t165 + t115;
t341 = t236 * t94;
t114 = -t239 * t152 + t242 * t159;
t93 = -t167 * qJ(4) + t114;
t83 = -pkin(3) * t263 + t93;
t42 = t234 * t83 + t341;
t36 = -pkin(10) * t263 + t42;
t151 = -pkin(2) * t224 - t182;
t120 = pkin(3) * t165 + qJD(4) + t151;
t55 = -pkin(4) * t284 - pkin(10) * t261 + t120;
t20 = t238 * t55 + t241 * t36;
t111 = pkin(3) * t297 + t176;
t33 = t95 * pkin(4) - pkin(10) * t251 + t111;
t6 = -qJD(5) * t20 - t238 * t15 + t241 * t33;
t1 = t95 * pkin(5) + t51 * qJ(6) - t103 * qJD(6) + t6;
t101 = t238 * t261 + t279;
t10 = -qJ(6) * t101 + t20;
t353 = t10 * t354 + t1;
t352 = t103 ^ 2;
t19 = -t238 * t36 + t241 * t55;
t9 = -qJ(6) * t103 + t19;
t8 = pkin(5) * t354 + t9;
t351 = t8 - t9;
t132 = -t147 * t241 + t238 * t294;
t216 = t349 * t239;
t217 = t349 * t242;
t150 = t216 * t234 - t217 * t236;
t142 = t241 * t150;
t260 = qJ(6) * t191 - qJD(6) * t198;
t348 = qJ(6) * t132 - t142 * qJD(5) + t260 * t241 + t360 + ((qJ(6) * t198 - t140) * qJD(5) + t365) * t238 - t318 * pkin(5);
t347 = (-qJD(5) * t150 + t260) * t238 + t355 + t364 * qJ(6);
t88 = t234 * t94;
t48 = t236 * t93 - t88;
t73 = pkin(3) * t167 + pkin(4) * t261 - pkin(10) * t284;
t346 = t238 * t73 + t241 * t48;
t247 = t238 * t251 - t241 * t268;
t52 = qJD(5) * t103 + t247;
t345 = -t101 * t307 - t238 * t52;
t192 = -t237 * t242 + t239 * t325;
t109 = -qJ(4) * t192 + t320;
t193 = t237 * t239 + t242 * t325;
t282 = -t179 * t239 + t242 * t180;
t98 = -pkin(3) * t324 - qJ(4) * t193 + t282;
t64 = t236 * t109 + t234 * t98;
t59 = -pkin(10) * t324 + t64;
t134 = t236 * t192 + t193 * t234;
t135 = -t192 * t234 + t193 * t236;
t178 = t225 + (-pkin(2) - t350) * t237;
t249 = t192 * pkin(3) + t178;
t79 = t134 * pkin(4) - t135 * pkin(10) + t249;
t344 = t238 * t79 + t241 * t59;
t292 = qJD(2) * t324;
t144 = -qJD(3) * t192 + t242 * t292;
t245 = -qJD(3) * t320 + t242 * t184 - t239 * t186;
t312 = qJD(2) * t240;
t293 = t235 * t312;
t60 = pkin(3) * t293 - t144 * qJ(4) - t193 * qJD(4) + t245;
t143 = qJD(3) * t193 + t239 * t292;
t254 = -t179 * t310 + t180 * t309 + t239 * t184 + t242 * t186;
t65 = -qJ(4) * t143 - qJD(4) * t192 + t254;
t27 = t234 * t60 + t236 * t65;
t342 = t150 * t95;
t339 = t51 * t238;
t229 = pkin(3) * t234 + pkin(10);
t322 = qJ(6) + t229;
t280 = qJD(5) * t322;
t330 = t284 * t238;
t335 = qJ(6) * t330 + t241 * qJD(6) - t238 * t280 - t346;
t72 = t241 * t73;
t334 = -pkin(5) * t261 - t72 + (qJ(6) * t284 - t280) * t241 + (-qJD(6) + t48) * t238;
t333 = t101 * t261;
t332 = t101 * t284;
t331 = t103 * t261;
t329 = t191 * t238;
t328 = t198 * t238;
t327 = t198 * t241;
t231 = t235 ^ 2;
t326 = t231 * qJD(1) ^ 2;
t321 = t238 * t140 + t142;
t187 = t237 * pkin(1) * t312 + pkin(8) * t292;
t317 = t240 ^ 2 - t243 ^ 2;
t311 = qJD(2) * t242;
t306 = t151 * qJD(3);
t305 = qJD(2) - t224;
t300 = t240 * t326;
t299 = t238 * t324;
t230 = -pkin(3) * t236 - pkin(4);
t290 = t231 * t303;
t47 = t234 * t93 + t341;
t17 = -t234 * t49 + t236 * t43;
t26 = -t234 * t65 + t236 * t60;
t41 = t236 * t83 - t88;
t288 = -t238 * t59 + t241 * t79;
t63 = -t234 * t109 + t236 * t98;
t286 = t131 + t329;
t285 = t191 * t241 + t132;
t149 = -t236 * t216 - t217 * t234;
t278 = t243 * t263;
t275 = t263 * t235;
t274 = t224 + t315;
t272 = 0.2e1 * t290;
t267 = t143 * pkin(3) + t187;
t58 = pkin(4) * t324 - t63;
t266 = -0.2e1 * pkin(1) * t290;
t5 = t241 * t15 + t238 * t33 + t55 * t307 - t308 * t36;
t2 = -qJ(6) * t52 - qJD(6) * t101 + t5;
t265 = -t354 * t8 + t2;
t262 = t143 * t234 - t144 * t236;
t259 = t241 * t95 + (-t308 + t330) * t354;
t117 = t135 * t238 + t241 * t324;
t24 = pkin(10) * t293 + t27;
t108 = t236 * t143 + t144 * t234;
t39 = t108 * pkin(4) + pkin(10) * t262 + t267;
t257 = t238 * t39 + t241 * t24 + t79 * t307 - t308 * t59;
t35 = pkin(4) * t263 - t41;
t256 = -t229 * t95 + t35 * t354;
t253 = -t364 - t329;
t252 = -t198 * t308 - t285;
t23 = -pkin(4) * t293 - t26;
t14 = -pkin(4) * t268 - t17;
t248 = -t344 * qJD(5) - t238 * t24 + t241 * t39;
t7 = pkin(5) * t52 + t14;
t195 = t322 * t241;
t194 = t322 * t238;
t138 = t241 * t140;
t118 = t135 * t241 - t299;
t100 = t101 ^ 2;
t86 = -qJ(6) * t328 + t321;
t81 = pkin(5) * t197 - qJ(6) * t327 - t150 * t238 + t138;
t68 = -qJD(5) * t299 + t135 * t307 - t238 * t262 - t241 * t293;
t67 = qJD(5) * t117 - t238 * t293 + t241 * t262;
t29 = t101 * pkin(5) + qJD(6) + t35;
t21 = -qJ(6) * t117 + t344;
t12 = pkin(5) * t134 - qJ(6) * t118 + t288;
t4 = -qJ(6) * t68 - qJD(6) * t117 + t257;
t3 = t108 * pkin(5) + t67 * qJ(6) - t118 * qJD(6) + t248;
t11 = [0, 0, 0, t240 * t243 * t272, -t317 * t272, t274 * t292, -t274 * t293, 0, -t176 * t237 - t187 * t224 + t240 * t266, -t175 * t237 - t186 * t224 + t243 * t266, t139 * t193 + t144 * t167, -t139 * t192 - t167 * t143 - t144 * t165 - t193 * t297, -t144 * t263 + (-t139 * t243 + (qJD(1) * t193 + t167) * t312) * t235, t143 * t263 + (t297 * t243 + (-qJD(1) * t192 - t165) * t312) * t235 (-t231 * t313 - t275) * t312, -t245 * t263 + t187 * t165 + t178 * t297 + t176 * t192 + t151 * t143 + (-t246 * t243 + (qJD(1) * t282 + t114) * t312) * t235, t254 * t263 + t187 * t167 + t178 * t139 + t176 * t193 + t151 * t144 + (-t255 * t243 + (-qJD(1) * t320 - t115) * t312) * t235, -t42 * t108 - t18 * t134 - t17 * t135 - t251 * t63 - t26 * t261 + t262 * t41 + t27 * t284 - t64 * t95, t111 * t249 + t120 * t267 + t17 * t63 + t18 * t64 + t41 * t26 + t42 * t27, -t103 * t67 - t118 * t51, t101 * t67 - t103 * t68 + t117 * t51 - t118 * t52, t103 * t108 + t118 * t95 - t134 * t51 - t354 * t67, -t101 * t108 - t117 * t95 - t134 * t52 - t354 * t68, t108 * t354 + t134 * t95, t23 * t101 + t19 * t108 + t14 * t117 + t6 * t134 + t248 * t354 + t288 * t95 + t35 * t68 + t58 * t52, t23 * t103 - t20 * t108 + t14 * t118 - t5 * t134 - t257 * t354 - t344 * t95 - t35 * t67 - t58 * t51, -t1 * t118 - t10 * t68 - t101 * t4 - t103 * t3 - t117 * t2 + t12 * t51 - t21 * t52 + t67 * t8, t2 * t21 + t10 * t4 + t1 * t12 + t8 * t3 + t7 * (pkin(5) * t117 + t58) + t29 * (pkin(5) * t68 + t23); 0, 0, 0, -t243 * t300, t317 * t326, t305 * t223, -t305 * t294, 0, pkin(1) * t300 + t185 * t224 - t176, pkin(8) * t268 + t182 * t224 + (-t237 * t303 + t326) * t350, t139 * t239 - t242 * t356 (t139 + t357) * t242 + (-t297 + t356) * t239, -t242 * t273 + (t242 * t278 + (qJD(2) * t239 - t167) * t240) * t316, t239 * t273 + (-t239 * t278 + (t165 + t311) * t240) * t316, t275 * t314, -pkin(2) * t297 + t239 * t306 + t281 * t263 - t185 * t165 - t363 * t242 + (-t114 * t240 + (-pkin(9) * t312 - t151 * t243) * t239) * t316, -pkin(2) * t139 + t242 * t306 - t319 * t263 - t185 * t167 + t363 * t239 + (-t151 * t323 + (-pkin(9) * t311 + t115) * t240) * t316, t149 * t251 - t17 * t198 - t18 * t197 + t337 * t261 + t336 * t284 + t318 * t42 + t368 * t41 - t342, t111 * t295 + t120 * t367 - t17 * t149 + t18 * t150 + t336 * t42 - t337 * t41, t103 * t252 - t327 * t51, t286 * t103 + t285 * t101 + (t339 - t241 * t52 + (t101 * t238 - t103 * t241) * qJD(5)) * t198, -t103 * t318 - t51 * t197 + t252 * t354 + t327 * t95, t101 * t318 - t52 * t197 - t253 * t354 - t328 * t95, t95 * t197 - t318 * t354, t138 * t95 + t149 * t52 + t6 * t197 - t364 * t35 - t318 * t19 + (-t150 * t307 + t360) * t354 + t338 * t101 + (t14 * t198 - t342 - t35 * t191 + (-qJD(5) * t140 + t365) * t354) * t238, -t321 * t95 - t5 * t197 - t149 * t51 + t14 * t327 + t318 * t20 + (t150 * t308 - t355) * t354 + t338 * t103 + t252 * t35, t81 * t51 - t86 * t52 + t285 * t8 - t348 * t103 - t347 * t101 + t286 * t10 + (-t1 * t241 - t2 * t238 + (-t10 * t241 + t238 * t8) * qJD(5)) * t198, t2 * t86 + t1 * t81 + t7 * (pkin(5) * t328 + t149) + t348 * t8 + (pkin(5) * t253 + t338) * t29 + t347 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t167 * t165, -t165 ^ 2 + t167 ^ 2, t139 - t357, -t297 - t356, t268, -t115 * t223 - t151 * t167 + t283, -t114 * t263 + t151 * t165 + t255 (-t234 * t95 - t236 * t251) * pkin(3) + (-t48 + t41) * t284 + (t42 - t47) * t261, t41 * t47 - t42 * t48 + (-t120 * t167 + t17 * t236 + t18 * t234) * pkin(3), t103 * t276 - t339 (-t51 + t332) * t241 - t359 + t345, -t331 - t362, t259 + t333, -t354 * t261, -t47 * t101 - t19 * t261 - t14 * t241 + t230 * t52 + (-t229 * t307 - t72) * t354 + (t354 * t48 + t256) * t238, -t47 * t103 + t20 * t261 + t14 * t238 - t230 * t51 + (t229 * t308 + t346) * t354 + t256 * t241, -t335 * t101 - t334 * t103 - t194 * t51 - t195 * t52 - t238 * t353 + t265 * t241, t2 * t195 - t1 * t194 + t7 * (-pkin(5) * t241 + t230) + t334 * t8 + (pkin(5) * t277 - t47) * t29 + t335 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t261 ^ 2 - t284 ^ 2, t261 * t41 - t284 * t42 + t111, 0, 0, 0, 0, 0, t259 - t333, -t331 + t362 (t51 + t332) * t241 + t359 + t345, t265 * t238 + t241 * t353 - t29 * t261; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103 * t101, -t100 + t352, t101 * t354 - t51, -t247 + (-qJD(5) + t354) * t103, t95, -t35 * t103 + t20 * t354 + t6, t101 * t35 + t19 * t354 - t5, pkin(5) * t51 - t351 * t101, t351 * t10 + (-t103 * t29 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100 - t352, t10 * t101 + t103 * t8 + t7;];
tauc_reg  = t11;
