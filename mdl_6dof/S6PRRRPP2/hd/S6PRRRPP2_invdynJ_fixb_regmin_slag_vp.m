% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRRPP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:52:44
% EndTime: 2019-03-08 22:52:57
% DurationCPUTime: 5.54s
% Computational Cost: add. (4471->535), mult. (10143->666), div. (0->0), fcn. (7781->10), ass. (0->255)
t203 = cos(qJ(4));
t200 = sin(qJ(4));
t301 = qJD(3) * t200;
t201 = sin(qJ(3));
t304 = qJD(2) * t201;
t150 = t203 * t304 + t301;
t204 = cos(qJ(3));
t303 = qJD(2) * t204;
t381 = qJD(4) - t303;
t325 = t150 * t381;
t289 = t201 * qJDD(2);
t297 = qJD(4) * t201;
t380 = qJD(2) * t297 - qJDD(3);
t62 = t200 * (qJD(3) * (qJD(4) + t303) + t289) + t380 * t203;
t387 = -t62 - t325;
t386 = -pkin(3) * t204 - pkin(2);
t320 = t200 * qJ(5);
t385 = pkin(3) + t320;
t157 = -pkin(9) * t201 + t386;
t316 = t203 * t204;
t185 = pkin(8) * t316;
t298 = qJD(4) * t200;
t199 = sin(pkin(6));
t306 = qJD(1) * t199;
t202 = sin(qJ(2));
t205 = cos(qJ(2));
t315 = t204 * t205;
t371 = t200 * t315 - t202 * t203;
t384 = qJD(4) * t185 + t157 * t298 - t371 * t306;
t104 = (t200 * t202 + t203 * t315) * t199;
t253 = pkin(3) * t201 - pkin(9) * t204;
t153 = t253 * qJD(3);
t296 = qJD(4) * t203;
t383 = -qJD(1) * t104 + t200 * t153 + t157 * t296;
t192 = t204 * qJDD(2);
t290 = qJD(2) * qJD(3);
t146 = t201 * t290 + qJDD(4) - t192;
t382 = t146 * qJ(5) + qJD(5) * t381;
t293 = t203 * qJD(3);
t148 = t200 * t304 - t293;
t379 = t62 * qJ(6) + t148 * qJD(6);
t300 = qJD(3) * t201;
t378 = qJ(5) * t300 - t204 * qJD(5) + t383;
t377 = -t203 * t153 + t384;
t154 = qJD(2) * pkin(8) + t202 * t306;
t336 = cos(pkin(6));
t262 = qJD(1) * t336;
t176 = t201 * t262;
t99 = t204 * t154 + t176;
t376 = t200 * qJD(5) + t99;
t327 = t148 * t381;
t272 = t204 * t290;
t61 = -qJD(4) * t293 + (-t272 - t289) * t203 + t380 * t200;
t375 = -t61 - t327;
t274 = t205 * t306;
t101 = qJD(2) * t157 - t274;
t89 = qJD(3) * pkin(9) + t99;
t29 = t203 * t101 - t200 * t89;
t313 = qJD(5) - t29;
t374 = 0.2e1 * t382;
t373 = -t201 * t154 + t204 * t262;
t145 = t150 ^ 2;
t372 = -t381 ^ 2 - t145;
t305 = qJD(2) * t199;
t279 = t202 * t305;
t323 = t199 * t202;
t131 = t201 * t323 - t204 * t336;
t278 = t205 * t305;
t80 = -qJD(3) * t131 + t204 * t278;
t322 = t199 * t204;
t132 = t201 * t336 + t202 * t322;
t321 = t199 * t205;
t83 = t132 * t203 - t200 * t321;
t13 = qJD(4) * t83 + t80 * t200 - t203 * t279;
t81 = qJD(3) * t132 + t201 * t278;
t82 = t132 * t200 + t203 * t321;
t370 = -t13 * t381 + t131 * t62 - t82 * t146 + t81 * t148;
t137 = t146 * pkin(4);
t369 = t137 - qJDD(5);
t14 = -qJD(4) * t82 + t200 * t279 + t80 * t203;
t368 = -t131 * t61 - t14 * t381 - t83 * t146 + t150 * t81;
t254 = qJD(3) * pkin(3) + t373;
t228 = t150 * qJ(5) + t254;
t31 = pkin(4) * t148 - t228;
t358 = pkin(9) * t146;
t367 = -t31 * t381 + t358;
t362 = pkin(4) + pkin(5);
t15 = -t148 * t362 + qJD(6) + t228;
t291 = qJD(1) * qJD(2);
t106 = qJDD(2) * pkin(8) + (qJDD(1) * t202 + t205 * t291) * t199;
t260 = qJDD(1) * t336;
t250 = t201 * t260;
t24 = qJDD(3) * pkin(9) + qJD(3) * t373 + t204 * t106 + t250;
t273 = t202 * t291;
t167 = t199 * t273;
t246 = -qJDD(1) * t321 + t167;
t46 = qJD(2) * t153 + qJDD(2) * t157 + t246;
t264 = t101 * t298 + t200 * t24 - t203 * t46 + t89 * t296;
t198 = sin(pkin(10));
t335 = cos(pkin(10));
t249 = t336 * t335;
t126 = t198 * t202 - t205 * t249;
t127 = t198 * t205 + t202 * t249;
t265 = t199 * t335;
t75 = t127 * t204 - t201 * t265;
t35 = -t126 * t203 + t200 * t75;
t266 = t198 * t336;
t128 = t202 * t335 + t205 * t266;
t129 = -t202 * t266 + t205 * t335;
t77 = t198 * t199 * t201 + t129 * t204;
t37 = -t128 * t203 + t200 * t77;
t219 = g(1) * t37 + g(2) * t35 + g(3) * t82 - t264;
t214 = t219 + t369;
t343 = t61 * qJ(6);
t366 = (qJD(6) + t15) * t150 + t214 - t343;
t365 = -t203 * t362 - t320;
t363 = t148 ^ 2;
t280 = -pkin(8) * t200 - pkin(4);
t292 = t203 * qJD(6);
t299 = qJD(3) * t204;
t361 = (-qJ(6) * t299 - t153) * t203 + (qJ(6) * t298 - t292 + (-pkin(5) + t280) * qJD(3)) * t201 + t384;
t318 = t201 * t203;
t360 = (-pkin(8) * qJD(3) + qJ(6) * qJD(4)) * t318 + (qJD(6) * t201 + (-pkin(8) * qJD(4) + qJ(6) * qJD(3)) * t204) * t200 + t378;
t357 = t146 * pkin(5);
t356 = pkin(9) - qJ(6);
t354 = (-t201 * t293 - t204 * t298) * pkin(8) + t378;
t352 = t280 * t300 + t377;
t30 = t200 * t101 + t203 * t89;
t351 = pkin(9) * qJD(4);
t350 = qJ(5) * t62;
t349 = qJD(2) * pkin(2);
t17 = qJ(6) * t148 + t30;
t173 = t381 * qJ(5);
t12 = t17 + t173;
t348 = t12 * t381;
t19 = t173 + t30;
t347 = t381 * t19;
t346 = t381 * t30;
t74 = -t127 * t201 - t204 * t265;
t345 = t203 * t74;
t76 = -t129 * t201 + t198 * t322;
t344 = t203 * t76;
t342 = t61 * t200;
t334 = qJ(5) * t203;
t238 = -t200 * t362 + t334;
t341 = t238 * t381 + t376;
t152 = t253 * qJD(2);
t337 = t200 * t152 + t203 * t373;
t42 = qJ(5) * t304 + t337;
t340 = -qJ(6) * t200 * t303 - t298 * t356 - t292 - t42;
t247 = pkin(4) * t200 - t334;
t339 = t247 * t381 - t376;
t163 = t356 * t203;
t86 = t200 * t373;
t263 = -t203 * t152 + t86;
t338 = qJD(4) * t163 - t200 * qJD(6) - (-qJ(6) * t316 - t201 * t362) * qJD(2) - t263;
t333 = qJ(6) * t201;
t331 = t126 * t201;
t330 = t128 * t201;
t329 = t131 * t203;
t328 = t148 * qJ(5);
t326 = t150 * t148;
t324 = t150 * t203;
t319 = t200 * t204;
t16 = t150 * qJ(6) + t29;
t314 = qJD(5) - t16;
t311 = qJDD(1) - g(3);
t308 = t200 * t157 + t185;
t195 = t201 ^ 2;
t307 = -t204 ^ 2 + t195;
t302 = qJD(3) * t150;
t295 = qJD(5) * t203;
t288 = pkin(4) * t345 + t385 * t74;
t287 = pkin(4) * t344 + t385 * t76;
t286 = t15 * t298;
t285 = t15 * t296;
t284 = t201 * t321;
t282 = -pkin(4) * t329 - t131 * t385;
t281 = -g(1) * t330 - g(2) * t331 + g(3) * t284;
t277 = t381 * t301;
t276 = t381 * t293;
t275 = t381 * t298;
t271 = t205 * t290;
t36 = t126 * t200 + t203 * t75;
t269 = -t35 * pkin(4) + qJ(5) * t36;
t38 = t128 * t200 + t203 * t77;
t268 = -t37 * pkin(4) + qJ(5) * t38;
t267 = -t82 * pkin(4) + qJ(5) * t83;
t184 = pkin(8) * t319;
t261 = t157 * t203 - t184;
t259 = -qJD(3) * t176 - t201 * t106 - t154 * t299 + t204 * t260;
t257 = t201 * t274;
t256 = t148 * t274;
t255 = t150 * t274;
t252 = g(1) * t128 + g(2) * t126;
t93 = -qJ(5) * t204 + t308;
t25 = -qJDD(3) * pkin(3) - t259;
t248 = pkin(4) * t203 + t320;
t18 = -pkin(4) * t381 + t313;
t244 = t18 * t203 - t19 * t200;
t208 = qJD(2) ^ 2;
t243 = qJDD(2) * t205 - t202 * t208;
t242 = pkin(8) + t247;
t5 = t264 - t369;
t207 = qJD(3) ^ 2;
t239 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t207 - t246;
t237 = t101 * t296 + t200 * t46 + t203 * t24 - t298 * t89;
t236 = t200 * t146 + t296 * t381;
t235 = t203 * t146 - t275;
t103 = t371 * t199;
t52 = -t126 * t319 - t127 * t203;
t54 = -t128 * t319 - t129 * t203;
t234 = g(1) * t54 + g(2) * t52 + g(3) * t103;
t53 = -t126 * t316 + t127 * t200;
t55 = -t128 * t316 + t129 * t200;
t233 = -g(1) * t55 - g(2) * t53 - g(3) * t104;
t232 = g(1) * t76 + g(2) * t74 - g(3) * t131;
t231 = g(1) * t77 + g(2) * t75 + g(3) * t132;
t229 = -pkin(8) + t238;
t6 = t62 * pkin(4) + t61 * qJ(5) - t150 * qJD(5) + t25;
t4 = -pkin(5) * t62 + qJDD(6) - t6;
t227 = -t232 + t4;
t226 = t199 * pkin(3) * t315 + pkin(2) * t321 + t104 * pkin(4) + pkin(8) * t323 + pkin(9) * t284 + qJ(5) * t103;
t224 = -t254 * t381 - t358;
t3 = t237 + t382;
t223 = t53 * pkin(4) + pkin(8) * t127 - pkin(9) * t331 + qJ(5) * t52 + t126 * t386;
t222 = t55 * pkin(4) + pkin(8) * t129 - pkin(9) * t330 + qJ(5) * t54 + t128 * t386;
t221 = -t146 + t326;
t220 = t13 * t150 - t14 * t148 - t61 * t82 - t62 * t83;
t218 = t351 * t381 + t232;
t217 = -t218 - t6;
t216 = t61 - t327;
t155 = -t274 - t349;
t215 = -pkin(8) * qJDD(3) + (t155 + t274 - t349) * qJD(3);
t213 = g(1) * t38 + g(2) * t36 + g(3) * t83 - t237;
t211 = t150 * t31 - t214;
t210 = t29 * t381 + t213;
t194 = t204 * pkin(4);
t162 = t356 * t200;
t156 = -pkin(3) - t248;
t144 = pkin(3) - t365;
t119 = t242 * t201;
t94 = t194 - t261;
t92 = t229 * t201;
t73 = pkin(4) * t150 + t328;
t68 = t200 * t333 + t93;
t63 = t204 * pkin(5) + t184 + t194 + (-t157 - t333) * t203;
t47 = -t150 * t362 - t328;
t45 = (qJD(4) * t248 - t295) * t201 + t242 * t299;
t44 = -pkin(4) * t304 + t263;
t22 = (qJD(4) * t365 + t295) * t201 + t229 * t299;
t7 = -t362 * t381 + t314;
t2 = t3 + t379;
t1 = -qJD(6) * t150 + t343 - t357 + t5;
t8 = [t311, 0, t243 * t199 (-qJDD(2) * t202 - t205 * t208) * t199, 0, 0, 0, 0, 0, -t81 * qJD(3) - t131 * qJDD(3) + (-t201 * t271 + t204 * t243) * t199, -t80 * qJD(3) - t132 * qJDD(3) + (-t201 * t243 - t204 * t271) * t199, 0, 0, 0, 0, 0, t370, t368, t370, t220, -t368, t13 * t18 + t131 * t6 + t14 * t19 + t3 * t83 + t31 * t81 + t5 * t82 - g(3), t370, -t368, -t220, t1 * t82 + t12 * t14 + t13 * t7 - t131 * t4 - t15 * t81 + t2 * t83 - g(3); 0, qJDD(2), t311 * t321 + t252, g(1) * t129 + g(2) * t127 - t311 * t323, qJDD(2) * t195 + 0.2e1 * t201 * t272, 0.2e1 * t192 * t201 - 0.2e1 * t290 * t307, qJDD(3) * t201 + t204 * t207, qJDD(3) * t204 - t201 * t207, 0, t215 * t201 + ((-g(3) * t205 + t273) * t199 + t239 + t252) * t204, t215 * t204 + (-t239 - t167) * t201 + t281, -t61 * t318 + (-t200 * t297 + t204 * t293) * t150 (-t148 * t203 - t150 * t200) * t299 + (t342 - t203 * t62 + (t148 * t200 - t324) * qJD(4)) * t201 (t61 + t276) * t204 + (t235 + t302) * t201 (t62 - t277) * t204 + (-qJD(3) * t148 - t236) * t201, -t146 * t204 + t300 * t381, t261 * t146 - t377 * t381 + ((pkin(8) * t148 - t200 * t254) * qJD(3) + t264) * t204 + (-t256 - t254 * t296 + t29 * qJD(3) + t25 * t200 + (t62 + t277) * pkin(8)) * t201 + t233, -t308 * t146 - t383 * t381 + (-t254 * t293 + (t275 + t302) * pkin(8) + t237) * t204 + (-t255 + t254 * t298 - t30 * qJD(3) + t25 * t203 + (-t61 + t276) * pkin(8)) * t201 + t234, t119 * t62 - t94 * t146 + t45 * t148 + (t301 * t31 + t5) * t204 - t352 * t381 + (-qJD(3) * t18 + t6 * t200 + t296 * t31 - t256) * t201 + t233, -t94 * t61 - t93 * t62 + t352 * t150 - t354 * t148 + t244 * t299 + (-t200 * t3 + t203 * t5 + (-t18 * t200 - t19 * t203) * qJD(4)) * t201 - t281, t119 * t61 + t93 * t146 - t45 * t150 + (-t293 * t31 - t3) * t204 + t354 * t381 + (qJD(3) * t19 - t6 * t203 + t298 * t31 + t255) * t201 - t234, t3 * t93 + t6 * t119 + t5 * t94 - g(1) * t222 - g(2) * t223 - g(3) * t226 + (t45 - t257) * t31 + t354 * t19 + t352 * t18, -t63 * t146 - t22 * t148 - t92 * t62 + (-t15 * t301 + t1) * t204 - t361 * t381 + (-qJD(3) * t7 - t4 * t200 - t256 - t285) * t201 + t233, t68 * t146 + t22 * t150 - t92 * t61 + (t15 * t293 - t2) * t204 + t360 * t381 + (qJD(3) * t12 + t4 * t203 + t255 - t286) * t201 - t234, t63 * t61 + t68 * t62 - t361 * t150 + t360 * t148 + (t12 * t200 - t203 * t7) * t299 + (-t1 * t203 + t2 * t200 + (t12 * t203 + t200 * t7) * qJD(4)) * t201 + t281, t2 * t68 + t1 * t63 + t4 * t92 - g(1) * (pkin(5) * t55 + qJ(6) * t330 + t222) - g(2) * (pkin(5) * t53 + qJ(6) * t331 + t223) - g(3) * (pkin(5) * t104 - qJ(6) * t284 + t226) + t361 * t7 + (t22 + t257) * t15 + t360 * t12; 0, 0, 0, 0, -t201 * t208 * t204, t307 * t208, t289, t192, qJDD(3), qJD(3) * t99 - t155 * t304 - t232 + t259, -t250 + (-qJD(2) * t155 - t106) * t204 + t231, t324 * t381 - t342, t200 * t387 + t375 * t203 (-t150 * t201 - t316 * t381) * qJD(2) + t236 (t148 * t201 + t319 * t381) * qJD(2) + t235, -t381 * t304, -t29 * t304 - pkin(3) * t62 - t99 * t148 + t86 * t381 + t224 * t200 + (-t25 - (t152 + t351) * t381 - t232) * t203, pkin(3) * t61 + t337 * t381 + t30 * t304 - t99 * t150 + t224 * t203 + (t218 + t25) * t200, t339 * t148 + t156 * t62 + t18 * t304 - t200 * t367 + t217 * t203 + t381 * t44, t42 * t148 - t44 * t150 + (t3 + t381 * t18 + (qJD(4) * t150 - t62) * pkin(9)) * t203 + (t5 - t347 + (qJD(4) * t148 - t61) * pkin(9)) * t200 - t231, -t339 * t150 + t156 * t61 - t19 * t304 + t217 * t200 + t203 * t367 - t381 * t42, t6 * t156 - t19 * t42 - t18 * t44 - g(1) * t287 - g(2) * t288 - g(3) * t282 + t339 * t31 + (qJD(4) * t244 + t5 * t200 + t3 * t203 - t231) * pkin(9), -t286 - t144 * t62 - t162 * t146 - t338 * t381 - t341 * t148 + (t15 * t319 + t201 * t7) * qJD(2) + t227 * t203, t285 - t144 * t61 + t163 * t146 + t340 * t381 + t341 * t150 + (-t12 * t201 - t15 * t316) * qJD(2) + t227 * t200, t162 * t61 + t163 * t62 - t338 * t150 + t340 * t148 + (-t381 * t7 - t2) * t203 + (-t1 + t348) * t200 + t231, t2 * t163 + t1 * t162 + t4 * t144 - g(1) * (pkin(5) * t344 + t356 * t77 + t287) - g(2) * (pkin(5) * t345 + t356 * t75 + t288) - g(3) * (-pkin(5) * t329 + t132 * t356 + t282) + t338 * t7 + t341 * t15 + t340 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t326, t145 - t363, -t216, -t62 + t325, t146, t150 * t254 + t219 + t346, -t148 * t254 + t210, -t148 * t73 + t137 - t211 + t346, pkin(4) * t61 - t350 + (t19 - t30) * t150 + (t18 - t313) * t148, -t148 * t31 + t150 * t73 - t210 + t374, -t5 * pkin(4) - g(1) * t268 - g(2) * t269 - g(3) * t267 + t3 * qJ(5) - t18 * t30 + t19 * t313 - t31 * t73, t148 * t47 + t17 * t381 + (pkin(5) + t362) * t146 + t366, t148 * t15 - t150 * t47 - t16 * t381 - t213 + t374 + t379, t350 - t362 * t61 + (-t12 + t17) * t150 + (-t7 + t314) * t148, t2 * qJ(5) - t1 * t362 - t7 * t17 - t15 * t47 - g(1) * (-pkin(5) * t37 + t268) - g(2) * (-pkin(5) * t35 + t269) - g(3) * (-pkin(5) * t82 + t267) + t314 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t221, -t216, t372, t211 - t347, t221, t372, t216, -t348 - t357 - t366; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t387, t375, -t145 - t363, -t12 * t148 + t150 * t7 + t227;];
tau_reg  = t8;
