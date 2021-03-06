% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6PRRPRR7
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRPRR7_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR7_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR7_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR7_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_invdynJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:34:19
% EndTime: 2019-03-08 22:34:34
% DurationCPUTime: 8.72s
% Computational Cost: add. (8184->669), mult. (17754->880), div. (0->0), fcn. (12808->14), ass. (0->324)
t232 = sin(qJ(3));
t358 = qJD(3) * t232;
t213 = pkin(3) * t358;
t236 = cos(qJ(3));
t397 = qJ(4) * t236;
t290 = pkin(9) * t232 - t397;
t355 = qJD(4) * t232;
t257 = qJD(3) * t290 - t355;
t117 = t213 + t257;
t239 = -pkin(3) - pkin(9);
t398 = qJ(4) * t232;
t319 = -pkin(2) - t398;
t152 = t236 * t239 + t319;
t356 = qJD(3) * t236;
t431 = pkin(4) + pkin(8);
t168 = t431 * t356;
t185 = t431 * t232;
t231 = sin(qJ(5));
t235 = cos(qJ(5));
t228 = sin(pkin(6));
t233 = sin(qJ(2));
t237 = cos(qJ(2));
t379 = t232 * t237;
t260 = (t231 * t379 + t233 * t235) * t228;
t353 = qJD(5) * t235;
t354 = qJD(5) * t231;
t403 = -qJD(1) * t260 + t235 * t117 - t152 * t354 + t231 * t168 + t185 * t353;
t376 = t235 * t237;
t261 = t228 * (-t231 * t233 + t232 * t376);
t460 = -qJD(1) * t261 - t117 * t231 + t235 * t168;
t359 = qJD(3) * t231;
t363 = qJD(2) * t236;
t159 = t235 * t363 + t359;
t327 = t231 * t363;
t357 = qJD(3) * t235;
t161 = -t327 + t357;
t230 = sin(qJ(6));
t234 = cos(qJ(6));
t288 = t159 * t230 - t234 * t161;
t75 = t234 * t159 + t161 * t230;
t425 = t75 * t288;
t164 = t231 * t185;
t381 = t231 * t232;
t279 = pkin(5) * t236 - pkin(10) * t381;
t318 = pkin(10) * t236 - t152;
t459 = t279 * qJD(3) + (t235 * t318 - t164) * qJD(5) + t460;
t352 = qJD(5) * t236;
t329 = t231 * t352;
t263 = t232 * t357 + t329;
t458 = -pkin(10) * t263 - t403;
t162 = t230 * t235 + t231 * t234;
t268 = t162 * t232;
t442 = qJD(5) + qJD(6);
t88 = t442 * t162;
t405 = qJD(2) * t268 + t88;
t423 = pkin(10) - t239;
t365 = qJD(2) * t232;
t214 = pkin(3) * t365;
t128 = qJD(2) * t290 + t214;
t369 = qJD(1) * t228;
t335 = t233 * t369;
t169 = qJD(2) * pkin(8) + t335;
t229 = cos(pkin(6));
t384 = t229 * t232;
t114 = qJD(1) * t384 + t169 * t236;
t212 = pkin(4) * t363;
t99 = t114 + t212;
t59 = -t128 * t231 + t235 * t99;
t457 = -qJD(2) * t279 + t423 * t354 - t59;
t172 = t423 * t235;
t330 = t235 * t365;
t60 = t235 * t128 + t231 * t99;
t456 = pkin(10) * t330 + qJD(5) * t172 + t60;
t349 = qJD(1) * qJD(2);
t324 = t237 * t349;
t360 = qJD(3) * t229;
t455 = qJDD(2) * pkin(8) + (qJDD(1) * t233 + t324) * t228 + qJD(1) * t360;
t204 = qJD(5) + t365;
t367 = qJD(1) * t237;
t334 = t228 * t367;
t102 = qJD(2) * t152 - t334;
t348 = qJD(2) * qJD(3);
t321 = t236 * t348;
t343 = t232 * qJDD(2);
t267 = t321 + t343;
t347 = qJDD(1) * t229;
t299 = t169 * t356 + t455 * t232 - t236 * t347;
t284 = qJDD(4) + t299;
t31 = t267 * pkin(4) + qJDD(3) * t239 + t284;
t325 = t233 * t349;
t385 = t228 * t237;
t289 = -qJDD(1) * t385 + t228 * t325;
t322 = t232 * t348;
t277 = pkin(3) * t322 + t289;
t41 = qJD(2) * t257 + qJDD(2) * t152 + t277;
t368 = qJD(1) * t236;
t195 = t229 * t368;
t443 = qJD(4) - t195;
t373 = (pkin(4) * qJD(2) + t169) * t232 + t443;
t73 = qJD(3) * t239 + t373;
t273 = t102 * t354 - t231 * t31 - t235 * t41 - t73 * t353;
t37 = -t102 * t231 + t235 * t73;
t454 = -t37 * t204 - t273;
t113 = t169 * t232 - t195;
t453 = -qJD(4) - t113;
t452 = t288 ^ 2 - t75 ^ 2;
t342 = t236 * qJDD(2);
t337 = qJD(3) * t353 + t231 * qJDD(3) + t235 * t342;
t250 = -qJD(2) * t263 + t337;
t350 = qJD(6) * t234;
t351 = qJD(6) * t230;
t276 = -t235 * qJDD(3) + (-t322 + t342) * t231;
t72 = qJD(5) * t159 + t276;
t19 = t159 * t350 + t161 * t351 + t230 * t250 + t234 * t72;
t191 = qJD(6) + t204;
t451 = t191 * t75 - t19;
t27 = -pkin(10) * t161 + t37;
t22 = pkin(5) * t204 + t27;
t38 = t102 * t235 + t231 * t73;
t28 = -pkin(10) * t159 + t38;
t158 = qJDD(5) + t267;
t9 = -qJD(5) * t38 - t231 * t41 + t235 * t31;
t6 = pkin(5) * t158 + pkin(10) * t72 + t9;
t7 = -pkin(10) * t250 - t273;
t1 = t234 * (qJD(6) * t22 + t7) + t230 * t6 - t28 * t351;
t400 = sin(pkin(11));
t309 = t400 * t233;
t401 = cos(pkin(11));
t310 = t401 * t237;
t138 = -t229 * t310 + t309;
t308 = t400 * t237;
t311 = t401 * t233;
t140 = t229 * t308 + t311;
t386 = t228 * t233;
t340 = t232 * t386;
t147 = -t229 * t236 + t340;
t227 = qJ(5) + qJ(6);
t217 = sin(t227);
t218 = cos(t227);
t224 = qJD(3) * qJ(4);
t107 = -t114 - t224;
t81 = -t107 + t212;
t62 = pkin(5) * t159 + t81;
t139 = t229 * t311 + t308;
t313 = t228 * t401;
t90 = t139 * t232 + t236 * t313;
t141 = -t229 * t309 + t310;
t312 = t228 * t400;
t92 = t141 * t232 - t236 * t312;
t450 = t62 * t75 - g(1) * (-t140 * t218 - t217 * t92) - g(2) * (-t138 * t218 - t217 * t90) - g(3) * (-t147 * t217 + t218 * t385) - t1;
t415 = t234 * t28;
t13 = t22 * t230 + t415;
t2 = -qJD(6) * t13 - t230 * t7 + t234 * t6;
t449 = t62 * t288 - g(1) * (-t140 * t217 + t218 * t92) - g(2) * (-t138 * t217 + t218 * t90) - g(3) * (t147 * t218 + t217 * t385) + t2;
t254 = qJD(6) * t288 + t230 * t72 - t234 * t250;
t448 = -t191 * t288 + t254;
t148 = t236 * t386 + t384;
t241 = qJD(2) ^ 2;
t280 = qJDD(2) * t237 - t233 * t241;
t320 = t237 * t348;
t331 = qJD(2) * t385;
t94 = -qJD(3) * t340 + (t331 + t360) * t236;
t447 = qJD(3) * t94 + qJDD(3) * t148 + t228 * (t232 * t280 + t236 * t320);
t95 = qJD(3) * t148 + t232 * t331;
t446 = -qJD(3) * t95 - qJDD(3) * t147 + t228 * (-t232 * t320 + t236 * t280);
t445 = t38 * t204 + t9;
t378 = t234 * t235;
t383 = t230 * t231;
t286 = -t378 + t383;
t129 = t286 * t236;
t106 = -qJD(3) * pkin(3) - t453;
t87 = t235 * t152 + t164;
t96 = t147 * t235 + t231 * t385;
t444 = -g(1) * (-t140 * t231 + t235 * t92) - g(2) * (-t138 * t231 + t235 * t90) - g(3) * t96;
t173 = -pkin(3) * t236 + t319;
t366 = qJD(2) * t173;
t116 = -t334 + t366;
t441 = t116 * t365 + qJDD(4);
t300 = t169 * t358 - t232 * t347 - t455 * t236;
t440 = (t113 * t236 - t114 * t232) * qJD(3) + t299 * t232 - t300 * t236;
t222 = qJDD(3) * qJ(4);
t223 = qJD(3) * qJD(4);
t46 = -t222 - t223 + t300;
t395 = qJDD(3) * pkin(3);
t47 = t284 - t395;
t439 = (t106 * t236 + t107 * t232) * qJD(3) + t47 * t232 - t46 * t236;
t438 = t19 * t286 + t288 * t405;
t406 = t230 * t354 + t231 * t351 - t234 * t330 + t365 * t383 - t378 * t442;
t437 = -t162 * t254 - t406 * t75;
t396 = qJDD(2) * pkin(2);
t119 = t289 - t396;
t293 = g(1) * t140 + g(2) * t138;
t240 = qJD(3) ^ 2;
t427 = pkin(8) * t240;
t436 = t228 * (-g(3) * t237 + t325) - t119 + t293 + t396 - t427;
t272 = -qJ(4) * t356 - t355;
t136 = t213 + t272;
t259 = g(3) * t385 - t293;
t346 = qJDD(2) * t173;
t61 = qJD(2) * t272 + t277 + t346;
t433 = qJD(2) * (-t136 + t335) - t259 - t346 - t427 - t61;
t417 = t230 * t28;
t12 = t22 * t234 - t417;
t270 = g(1) * t92 + g(2) * t90 + g(3) * t147;
t432 = t1 * t162 - t12 * t405 - t13 * t406 - t2 * t286 - t270;
t165 = t235 * t185;
t63 = pkin(5) * t232 + t231 * t318 + t165;
t377 = t235 * t236;
t70 = -pkin(10) * t377 + t87;
t25 = -t230 * t70 + t234 * t63;
t429 = qJD(6) * t25 + t459 * t230 - t458 * t234;
t26 = t230 * t63 + t234 * t70;
t428 = -qJD(6) * t26 + t458 * t230 + t459 * t234;
t426 = g(3) * t233;
t208 = pkin(5) * t235 + pkin(4);
t424 = pkin(8) + t208;
t171 = t423 * t231;
t100 = t171 * t230 - t172 * t234;
t422 = qJD(6) * t100 + t457 * t230 - t456 * t234;
t101 = -t171 * t234 - t172 * t230;
t421 = -qJD(6) * t101 + t456 * t230 + t457 * t234;
t420 = qJD(2) * pkin(2);
t416 = t231 * t37;
t414 = t235 * t72;
t407 = t72 * t231;
t404 = pkin(5) * t353 - (-qJD(2) * t208 - t169) * t232 + t443;
t402 = -qJD(5) * t87 + t460;
t399 = pkin(8) * qJDD(3);
t394 = t138 * t236;
t393 = t140 * t236;
t392 = t158 * t231;
t391 = t159 * t235;
t390 = t161 * t159;
t389 = t204 * t239;
t388 = t217 * t232;
t387 = t218 * t232;
t382 = t231 * t159;
t380 = t232 * t235;
t137 = t235 * t158;
t375 = t236 * t237;
t374 = t239 * t158;
t372 = qJDD(1) - g(3);
t371 = pkin(2) * t385 + pkin(8) * t386;
t186 = t431 * t236;
t225 = t232 ^ 2;
t226 = t236 ^ 2;
t370 = t225 - t226;
t364 = qJD(2) * t233;
t362 = qJD(3) * t159;
t361 = qJD(3) * t161;
t345 = qJDD(2) * t225;
t344 = qJDD(2) * t226;
t341 = t159 * t380;
t132 = t138 * pkin(2);
t339 = -pkin(3) * t394 - t138 * t398 - t132;
t133 = t140 * pkin(2);
t338 = -pkin(3) * t393 - t140 * t398 - t133;
t336 = qJ(4) * t363;
t333 = t236 * t367;
t332 = t228 * t364;
t328 = t235 * t352;
t206 = pkin(5) * t231 + qJ(4);
t84 = t90 * pkin(3);
t91 = t139 * t236 - t232 * t313;
t317 = qJ(4) * t91 - t84;
t85 = t92 * pkin(3);
t93 = t141 * t236 + t232 * t312;
t316 = qJ(4) * t93 - t85;
t306 = qJD(2) * t87 + t38;
t135 = t147 * pkin(3);
t304 = qJ(4) * t148 - t135;
t303 = -qJD(2) * t186 - t81;
t301 = t204 + t365;
t298 = t232 * t321;
t151 = qJDD(6) + t158;
t295 = -t286 * t151 - t405 * t191;
t294 = t337 * t231;
t292 = g(1) * t141 + g(2) * t139;
t97 = -t147 * t231 + t228 * t376;
t44 = t230 * t97 + t234 * t96;
t45 = t230 * t96 - t234 * t97;
t287 = -t161 * t235 + t382;
t285 = g(3) * (t371 + (pkin(3) * t375 + qJ(4) * t379) * t228);
t283 = t204 * t231;
t238 = -pkin(10) - pkin(9);
t278 = pkin(5) * t381 - t236 * t238;
t274 = -t204 * t353 - t392;
t271 = -t151 * t162 + t191 * t406;
t269 = -g(1) * t93 - g(2) * t91 - g(3) * t148;
t264 = pkin(4) * t342 - t46;
t32 = -pkin(4) * t322 + t264;
t262 = t269 + t32;
t258 = t235 * t322 - t337;
t255 = t270 - t299;
t170 = -t334 - t420;
t252 = -t399 + (t170 + t334 - t420) * qJD(3);
t251 = t399 + (-t116 - t334 - t366) * qJD(3);
t249 = qJD(3) * t114 + t255;
t248 = qJD(3) * t113 + t269 - t300;
t243 = (-t426 + (-t225 - t226) * t324) * t228 - t292 + (t345 + t344) * pkin(8);
t242 = (t147 * t232 + t148 * t236) * qJDD(2) + (t232 * t95 + t236 * t94 + (t147 * t236 - t148 * t232) * qJD(3)) * qJD(2);
t199 = t232 * t241 * t236;
t177 = t370 * t241;
t175 = qJDD(3) * t236 - t232 * t240;
t174 = qJDD(3) * t232 + t236 * t240;
t167 = t431 * t358;
t166 = t214 - t336;
t154 = -0.2e1 * t298 + t344;
t153 = 0.2e1 * t298 + t345;
t146 = pkin(5) * t377 + t186;
t130 = t162 * t236;
t115 = 0.2e1 * t232 * t342 - 0.2e1 * t348 * t370;
t103 = -pkin(5) * t329 - t358 * t424;
t86 = -t152 * t231 + t165;
t50 = t236 * t88 - t286 * t358;
t49 = qJD(3) * t268 + t129 * t442;
t36 = qJD(5) * t96 + t231 * t95 + t235 * t332;
t35 = qJD(5) * t97 - t231 * t332 + t235 * t95;
t18 = t337 * pkin(5) + (-pkin(4) * t358 - pkin(5) * t263) * qJD(2) + t264;
t15 = t234 * t27 - t417;
t14 = -t230 * t27 - t415;
t11 = -qJD(6) * t45 - t230 * t36 + t234 * t35;
t10 = qJD(6) * t44 + t230 * t35 + t234 * t36;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t372, 0, 0, 0, 0, 0, 0, t280 * t228 (-qJDD(2) * t233 - t237 * t241) * t228, 0, -g(3) + (t229 ^ 2 + (t233 ^ 2 + t237 ^ 2) * t228 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, t446, -t447, t242, t113 * t95 + t114 * t94 + t147 * t299 - t148 * t300 - g(3) + (-t119 * t237 + t170 * t364) * t228, 0, 0, 0, 0, 0, 0, t242, -t446, t447, t106 * t95 - t107 * t94 + t147 * t47 - t148 * t46 - g(3) + (t116 * t364 - t237 * t61) * t228, 0, 0, 0, 0, 0, 0, t148 * t250 + t96 * t158 + t94 * t159 + t35 * t204, -t148 * t72 + t158 * t97 + t161 * t94 - t204 * t36, -t36 * t159 - t35 * t161 + t250 * t97 + t96 * t72, t148 * t32 + t273 * t97 + t35 * t37 + t36 * t38 + t81 * t94 + t9 * t96 - g(3), 0, 0, 0, 0, 0, 0, t11 * t191 - t148 * t254 + t151 * t44 + t75 * t94, -t10 * t191 - t148 * t19 - t151 * t45 - t288 * t94, -t10 * t75 + t11 * t288 + t19 * t44 + t254 * t45, t1 * t45 + t10 * t13 + t11 * t12 + t148 * t18 + t2 * t44 + t62 * t94 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t372 * t385 + t293, -t372 * t386 + t292, 0, 0, t153, t115, t174, t154, t175, 0, t252 * t232 + t236 * t436, -t232 * t436 + t252 * t236, t243 + t440, -t119 * pkin(2) + g(1) * t133 + g(2) * t132 - g(3) * t371 + (-t170 * t233 + (-t113 * t232 - t114 * t236) * t237) * t369 + (-t292 + t440) * pkin(8), 0, -t174, -t175, t153, t115, t154, t243 + t439, t251 * t232 - t236 * t433, t232 * t433 + t251 * t236, t61 * t173 + t116 * t136 - g(1) * t338 - g(2) * t339 - t285 + (-t116 * t233 + (-t106 * t232 + t107 * t236) * t237) * t369 + (-t292 + t439) * pkin(8), t236 * t407 + (t231 * t358 - t328) * t161, -t287 * t358 + (-t258 * t231 + t414 + (t391 + (t161 - t327) * t231) * qJD(5)) * t236 (t204 * t359 - t72) * t232 + (t274 + t361) * t236, -t159 * t263 + t250 * t377 (t301 * t357 - t337) * t232 + (t301 * t354 - t137 - t362) * t236, t158 * t232 + t204 * t356, t86 * t158 - t167 * t159 + t186 * t337 - t292 * t235 + (t231 * t293 + t303 * t357 + t9) * t232 - g(3) * t260 + t402 * t204 + (t37 * qJD(3) - t159 * t334 + t32 * t235 + t303 * t354) * t236, -t87 * t158 - t167 * t161 - t186 * t72 + t292 * t231 + (t235 * t293 + t359 * t81 + t273) * t232 - g(3) * t261 - t403 * t204 + (-qJD(3) * t38 - t161 * t334 - t32 * t231 - t353 * t81) * t236, -t87 * t337 + t86 * t72 - t402 * t161 - t403 * t159 + (t235 * t306 - t416) * t358 + (t9 * t231 + t273 * t235 + (t231 * t306 + t235 * t37) * qJD(5) - t259) * t236, -t273 * t87 + t9 * t86 + t32 * t186 - t81 * t167 - g(1) * (-pkin(9) * t393 + t141 * t431 + t338) - g(2) * (-pkin(9) * t394 + t139 * t431 + t339) - t285 + t403 * t38 + t402 * t37 + (-t81 * t333 - g(3) * (pkin(4) * t233 + pkin(9) * t375)) * t228, t130 * t19 - t288 * t49, -t129 * t19 - t130 * t254 - t288 * t50 - t49 * t75, -t130 * t151 - t19 * t232 + t191 * t49 - t288 * t356, t129 * t254 - t50 * t75, t129 * t151 + t191 * t50 + t232 * t254 - t356 * t75, t151 * t232 + t191 * t356, t25 * t151 + t2 * t232 + t12 * t356 + t103 * t75 - t146 * t254 - t18 * t129 - t62 * t50 - g(1) * (-t140 * t388 + t141 * t218) - g(2) * (-t138 * t388 + t139 * t218) + t428 * t191 + (-t75 * t333 - g(3) * (t217 * t379 + t218 * t233)) * t228, -t26 * t151 - t1 * t232 - t13 * t356 - t103 * t288 - t146 * t19 - t18 * t130 + t62 * t49 - g(1) * (-t140 * t387 - t141 * t217) - g(2) * (-t138 * t387 - t139 * t217) - t429 * t191 + (t288 * t333 - g(3) * (-t217 * t233 + t218 * t379)) * t228, t1 * t129 - t12 * t49 + t13 * t50 + t130 * t2 + t19 * t25 - t236 * t259 + t254 * t26 + t288 * t428 - t429 * t75, t1 * t26 + t2 * t25 + t18 * t146 + t62 * t103 - g(1) * (-t140 * t278 + t141 * t424 + t338) - g(2) * (-t138 * t278 + t139 * t424 + t339) - t285 + t429 * t13 + t428 * t12 + (-t208 * t426 + (-g(3) * t278 - t368 * t62) * t237) * t228; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t199, t177, t343, t199, t342, qJDD(3), -t170 * t365 + t249, -t170 * t363 - t248, 0, 0, qJDD(3), -t343, -t342, -t199, t177, t199 (-pkin(3) * t232 + t397) * qJDD(2), -t166 * t363 - t249 - 0.2e1 * t395 + t441, 0.2e1 * t222 + 0.2e1 * t223 + (t116 * t236 + t166 * t232) * qJD(2) + t248, -t47 * pkin(3) - g(1) * t316 - g(2) * t317 - g(3) * t304 - t46 * qJ(4) - t106 * t114 + t453 * t107 - t116 * t166, -t161 * t283 - t414, -t337 * t235 + t407 + t287 * qJD(5) + (t231 * t328 + (t382 + (-t161 + t357) * t235) * t232) * qJD(2), -t204 * t354 + t137 + (-t161 * t236 - t204 * t381) * qJD(2), t294 + t159 * t353 + (-t231 * t263 + t341) * qJD(2) (t159 * t236 - t204 * t380) * qJD(2) + t274, -t204 * t363, qJ(4) * t337 - t59 * t204 - t37 * t363 + t373 * t159 + (t81 * qJD(5) + t374 + (t81 - t224) * t365) * t235 + ((-t336 - t389) * qJD(5) + t262) * t231, t38 * t363 - qJ(4) * t72 + t204 * t60 + t373 * t161 + (-t204 * t81 - t374) * t231 + (-qJD(5) * t389 + t262) * t235, t60 * t159 + t59 * t161 + (-t38 * t365 + t239 * t72 - t9 + (-t159 * t239 - t38) * qJD(5)) * t235 + (t239 * t258 + t273 + t37 * t365 + (t37 + (t161 + t327) * t239) * qJD(5)) * t231 + t270, t32 * qJ(4) - t38 * t60 - t37 * t59 - g(1) * (-pkin(9) * t92 + t316) - g(2) * (-pkin(9) * t90 + t317) - g(3) * (-pkin(9) * t147 + t304) + t373 * t81 + (-t273 * t231 + t9 * t235 + (t235 * t38 - t416) * qJD(5)) * t239, t438, t162 * t19 - t254 * t286 - t288 * t406 + t405 * t75, t288 * t363 + t295, t437, t363 * t75 + t271, -t191 * t363, t100 * t151 - t12 * t363 + t162 * t18 + t191 * t421 - t206 * t254 + t217 * t269 + t404 * t75 - t406 * t62, -t101 * t151 + t13 * t363 - t18 * t286 - t19 * t206 - t191 * t422 + t218 * t269 - t288 * t404 - t405 * t62, t100 * t19 + t101 * t254 + t288 * t421 - t422 * t75 - t432, t1 * t101 + t2 * t100 + t18 * t206 - g(1) * (t206 * t93 + t238 * t92 - t85) - g(2) * (t206 * t91 + t238 * t90 - t84) - g(3) * (t147 * t238 + t148 * t206 - t135) + t404 * t62 + t422 * t13 + t421 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t343, qJDD(3) + t199, -t225 * t241 - t240, qJD(3) * t107 - t255 - t395 + t441, 0, 0, 0, 0, 0, 0, -t204 * t283 + t137 - t362, -t204 ^ 2 * t235 - t361 - t392, -t294 + t414 + (t161 * t231 - t391) * qJD(5) + (-t341 + (t329 + (t161 + t357) * t232) * t231) * qJD(2), -qJD(3) * t81 + t454 * t231 + t445 * t235 - t270, 0, 0, 0, 0, 0, 0, -qJD(3) * t75 + t295, qJD(3) * t288 + t271, -t437 - t438, -qJD(3) * t62 + t432; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t390, -t159 ^ 2 + t161 ^ 2, -t276 + (-qJD(5) + t204) * t159, -t390, t161 * t204 - t250, t158, -t81 * t161 + t444 + t445, t81 * t159 - g(1) * (-t140 * t235 - t231 * t92) - g(2) * (-t138 * t235 - t231 * t90) - g(3) * t97 - t454, 0, 0, -t425, t452, t451, t425, t448, t151, -t14 * t191 + (t151 * t234 - t161 * t75 - t191 * t351) * pkin(5) + t449, t15 * t191 + (-t151 * t230 + t161 * t288 - t191 * t350) * pkin(5) + t450, -t12 * t75 - t13 * t288 - t14 * t288 + t15 * t75 + (t19 * t234 + t254 * t230 + (-t230 * t288 - t234 * t75) * qJD(6)) * pkin(5), -t12 * t14 - t13 * t15 + (t1 * t230 + t2 * t234 - t62 * t161 + (-t12 * t230 + t13 * t234) * qJD(6) + t444) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t425, t452, t451, t425, t448, t151, t13 * t191 + t449, t12 * t191 + t450, 0, 0;];
tau_reg  = t3;
