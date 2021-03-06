% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRR12_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR12_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:55:13
% EndTime: 2019-12-31 22:55:51
% DurationCPUTime: 15.57s
% Computational Cost: add. (23896->666), mult. (73105->963), div. (0->0), fcn. (59368->12), ass. (0->292)
t279 = cos(qJ(2));
t406 = cos(pkin(5));
t369 = pkin(1) * t406;
t265 = t279 * t369;
t257 = qJD(1) * t265;
t277 = sin(qJ(2));
t273 = sin(pkin(5));
t405 = cos(pkin(6));
t312 = t273 * (-pkin(9) * t405 - pkin(8));
t304 = t277 * t312;
t186 = qJD(1) * t304 + t257;
t264 = t277 * t369;
t284 = t279 * t312 - t264;
t187 = t284 * qJD(1);
t272 = sin(pkin(6));
t426 = pkin(9) * t272;
t303 = (pkin(2) * t277 - t279 * t426) * t273;
t214 = qJD(1) * t303;
t428 = cos(qJ(3));
t336 = t405 * t428;
t276 = sin(qJ(3));
t400 = t272 * t276;
t234 = pkin(2) * t336 - pkin(9) * t400;
t353 = t276 * t405;
t391 = t234 * qJD(3) - t428 * t186 - t187 * t353 - t214 * t400;
t295 = -t277 * t353 + t428 * t279;
t387 = qJD(1) * t273;
t206 = t295 * t387;
t360 = qJD(3) * t428;
t340 = t272 * t360;
t447 = t340 - t206;
t130 = -t187 * t272 + t405 * t214;
t293 = t276 * t279 + t277 * t336;
t205 = t293 * t387;
t446 = -pkin(3) * t205 + pkin(10) * t206 - t130 + (pkin(3) * t276 - t428 * pkin(10)) * t272 * qJD(3);
t363 = t277 * t387;
t343 = t272 * t363;
t445 = -pkin(10) * t343 + t391;
t275 = sin(qJ(4));
t278 = cos(qJ(4));
t232 = t275 * t400 - t278 * t405;
t393 = t232 * qJD(4) + t275 * t343 - t447 * t278;
t233 = t275 * t405 + t278 * t400;
t392 = qJD(4) * t233 + t447 * t275 + t278 * t343;
t385 = qJD(3) * t276;
t361 = t272 * t385;
t439 = t205 - t361;
t368 = t272 * t428;
t236 = pkin(2) * t353 + pkin(9) * t368;
t444 = t236 * qJD(3) - t276 * t186 + t187 * t336;
t342 = t276 * t363;
t348 = t406 * qJD(1);
t321 = t348 + qJD(2);
t305 = t272 * t321;
t316 = t279 * t336;
t306 = t273 * t316;
t438 = -qJD(1) * t306 - t428 * t305;
t154 = t342 + t438;
t298 = qJD(4) + t154;
t222 = pkin(10) * t405 + t236;
t223 = (-t428 * pkin(3) - pkin(10) * t276 - pkin(2)) * t272;
t382 = qJD(4) * t278;
t383 = qJD(4) * t275;
t411 = -t222 * t383 + t223 * t382 + t275 * t446 + t445 * t278;
t365 = t428 * t214;
t407 = -(-pkin(3) * t363 - t365) * t272 + t444;
t443 = t439 * pkin(11) - t411;
t442 = -t392 * pkin(4) - t393 * pkin(11) - t407;
t290 = qJD(4) * t298;
t379 = qJD(1) * qJD(2);
t356 = t273 * t379;
t339 = t277 * t356;
t322 = t272 * t339;
t398 = t273 * t279;
t389 = pkin(8) * t398 + t264;
t226 = t389 * qJD(1);
t351 = t279 * t405;
t337 = t273 * t351;
t149 = t226 + (qJD(1) * t337 + t305) * pkin(9);
t285 = pkin(2) * t406 + t304;
t153 = qJD(2) * pkin(2) + qJD(1) * t285 + t257;
t255 = qJD(2) * t257;
t296 = qJD(2) * t304;
t170 = qJD(1) * t296 + t255;
t189 = t284 * qJD(2);
t171 = qJD(1) * t189;
t211 = (-pkin(2) * t279 - t277 * t426 - pkin(1)) * t273;
t201 = qJD(1) * t211;
t215 = qJD(2) * t303;
t208 = qJD(1) * t215;
t335 = qJD(3) * t353;
t49 = -t149 * t360 - t153 * t335 - t276 * t170 + t171 * t336 - t201 * t361 + t208 * t368;
t47 = -pkin(3) * t322 - t49;
t441 = pkin(10) * t290 + t47;
t83 = t428 * t149 + (t153 * t405 + t201 * t272) * t276;
t440 = -t83 + t298 * (pkin(4) * t275 - pkin(11) * t278);
t115 = -t153 * t272 + t405 * t201;
t294 = t276 * t351 + t428 * t277;
t288 = t294 * t273;
t299 = t276 * t305;
t156 = qJD(1) * t288 + t299;
t71 = pkin(3) * t154 - pkin(10) * t156 + t115;
t372 = t272 * t398;
t380 = qJD(1) * t372 - qJD(3);
t286 = -t321 * t405 + t380;
t73 = -pkin(10) * t286 + t83;
t30 = -t275 * t73 + t278 * t71;
t315 = qJD(3) * t336;
t292 = t149 * t385 - t153 * t315 - t428 * t170 - t171 * t353 - t201 * t340 - t208 * t400;
t46 = pkin(10) * t322 - t292;
t318 = t428 * t356;
t123 = (qJD(2) * t405 + qJD(3)) * t342 - t279 * t318 + t438 * qJD(3);
t429 = qJD(2) * t293 + qJD(3) * t294;
t281 = t429 * t273;
t291 = qJD(3) * t299;
t124 = qJD(1) * t281 + t291;
t125 = -t171 * t272 + t405 * t208;
t61 = pkin(3) * t124 + pkin(10) * t123 + t125;
t309 = -t275 * t61 - t278 * t46 - t71 * t382 + t383 * t73;
t437 = -t298 * t30 - t309;
t269 = t273 ^ 2;
t436 = -0.2e1 * t269 * t379;
t200 = t278 * t286;
t120 = t156 * t275 + t200;
t119 = qJD(5) + t120;
t274 = sin(qJ(5));
t31 = t275 * t71 + t278 * t73;
t29 = pkin(11) * t298 + t31;
t122 = t278 * t156 - t275 * t286;
t82 = -t276 * t149 + t153 * t336 + t201 * t368;
t72 = pkin(3) * t286 - t82;
t36 = t120 * pkin(4) - t122 * pkin(11) + t72;
t427 = cos(qJ(5));
t14 = t274 * t36 + t427 * t29;
t68 = qJD(4) * t200 + t278 * t123 + t156 * t383 - t275 * t322;
t347 = -t275 * t123 - t278 * t322;
t384 = qJD(4) * t122;
t69 = t347 + t384;
t16 = pkin(4) * t69 + pkin(11) * t68 + t47;
t5 = pkin(11) * t124 - t309;
t2 = -qJD(5) * t14 + t427 * t16 - t274 * t5;
t434 = -t14 * t119 - t2;
t430 = t278 * t222 + t275 * t223;
t410 = -qJD(4) * t430 - t445 * t275 + t278 * t446;
t432 = t120 * t298;
t431 = t122 * t298;
t390 = t272 * t365 + t444;
t354 = t275 * t46 - t278 * t61;
t8 = -qJD(4) * t31 - t354;
t90 = t427 * t122 + t274 * t298;
t27 = qJD(5) * t90 - t427 * t124 - t274 * t68;
t185 = t265 + t285;
t127 = -t185 * t272 + t405 * t211;
t350 = t406 * t272;
t317 = t428 * t350;
t396 = t276 * t277;
t178 = t273 * t396 - t306 - t317;
t338 = t276 * t350;
t179 = t338 + t288;
t80 = pkin(3) * t178 - pkin(10) * t179 + t127;
t329 = t406 * t405;
t231 = -t329 + t372;
t172 = (t337 + t350) * pkin(9) + t389;
t97 = t428 * t172 + t185 * t353 + t211 * t400;
t85 = -pkin(10) * t231 + t97;
t420 = t275 * t80 + t278 * t85;
t386 = qJD(2) * t273;
t362 = t277 * t386;
t341 = t272 * t362;
t258 = qJD(2) * t265;
t188 = t258 + t296;
t62 = -t172 * t385 + t185 * t315 + t428 * t188 + t189 * t353 + t211 * t340 + t215 * t400;
t54 = pkin(10) * t341 + t62;
t128 = qJD(3) * t338 + t281;
t129 = qJD(3) * t317 + ((t316 - t396) * qJD(3) + t295 * qJD(2)) * t273;
t131 = -t189 * t272 + t405 * t215;
t66 = pkin(3) * t128 - pkin(10) * t129 + t131;
t12 = -qJD(4) * t420 - t275 * t54 + t278 * t66;
t310 = t274 * t29 - t427 * t36;
t1 = -qJD(5) * t310 + t274 * t16 + t427 * t5;
t280 = qJD(1) ^ 2;
t6 = -pkin(4) * t124 - t8;
t425 = t6 * t274;
t287 = t427 * t298;
t88 = t122 * t274 - t287;
t424 = t90 * t88;
t221 = -pkin(3) * t405 - t234;
t134 = t232 * pkin(4) - t233 * pkin(11) + t221;
t136 = -pkin(11) * t368 + t430;
t87 = t274 * t134 + t427 * t136;
t423 = t87 * qJD(5) - t443 * t274 + t442 * t427;
t86 = t427 * t134 - t274 * t136;
t422 = -t86 * qJD(5) + t442 * t274 + t443 * t427;
t421 = t439 * pkin(4) - t410;
t132 = t179 * t275 + t231 * t278;
t419 = t132 * t69;
t417 = t232 * t69;
t381 = qJD(5) * t274;
t26 = -qJD(5) * t287 + t122 * t381 - t274 * t124 + t427 * t68;
t416 = t26 * t274;
t415 = t274 * t69;
t414 = t274 * t88;
t413 = t278 * t69;
t412 = t90 * t119;
t254 = -pkin(4) * t278 - pkin(11) * t275 - pkin(3);
t358 = qJD(5) * t427;
t359 = qJD(4) * t427;
t109 = pkin(3) * t156 + pkin(10) * t154;
t57 = t275 * t109 + t278 * t82;
t45 = pkin(11) * t156 + t57;
t409 = -t427 * t45 + t254 * t358 + (-t275 * t359 - t278 * t381) * pkin(10) + t440 * t274;
t408 = t274 * t45 - t254 * t381 + (t274 * t383 - t278 * t358) * pkin(10) + t440 * t427;
t404 = t122 * t120;
t403 = t124 * t178;
t402 = t156 * t154;
t401 = t269 * t280;
t399 = t273 * t277;
t397 = t274 * t278;
t192 = t274 * t233 + t427 * t368;
t395 = t192 * qJD(5) + t439 * t274 + t393 * t427;
t344 = t274 * t368;
t394 = -qJD(5) * t344 + t233 * t358 - t393 * t274 + t439 * t427;
t388 = t277 ^ 2 - t279 ^ 2;
t378 = t1 * t427;
t377 = t6 * t427;
t376 = t26 * t427;
t375 = t27 * t427;
t374 = t427 * t69;
t373 = t279 * t401;
t367 = t278 * t427;
t366 = t428 * t124;
t364 = t428 * t215;
t355 = t279 * t379;
t346 = t298 * t275;
t345 = t277 * t373;
t334 = t273 * t280 * t406;
t333 = pkin(1) * t436;
t106 = -t154 * t397 - t427 * t156;
t330 = -t274 * t382 + t106;
t39 = -t275 * t85 + t278 * t80;
t56 = t109 * t278 - t275 * t82;
t133 = t179 * t278 - t231 * t275;
t139 = -t275 * t222 + t278 * t223;
t323 = t269 * t277 * t355;
t320 = 0.2e1 * t348 + qJD(2);
t107 = -t154 * t367 + t274 * t156;
t313 = t278 * t359 - t107;
t38 = pkin(11) * t178 + t420;
t96 = -t276 * t172 + t185 * t336 + t211 * t368;
t84 = t231 * pkin(3) - t96;
t50 = t132 * pkin(4) - t133 * pkin(11) + t84;
t17 = -t274 * t38 + t427 * t50;
t311 = -t14 * t274 + t310 * t427;
t18 = t274 * t50 + t427 * t38;
t11 = t275 * t66 + t278 * t54 + t80 * t382 - t383 * t85;
t103 = t133 * t427 + t178 * t274;
t308 = -t133 * t274 + t178 * t427;
t307 = -pkin(10) * t124 + t298 * t72;
t218 = -pkin(8) * t339 + t255;
t301 = t119 * t381 - t374;
t300 = t427 * t120 + t358;
t230 = t389 * qJD(2);
t289 = -t172 * t360 - t185 * t335 - t276 * t188 + t189 * t336 - t211 * t361;
t283 = qJD(3) * t286;
t282 = t154 * t298 + t290;
t55 = (-pkin(3) * t362 - t364) * t272 - t289;
t268 = t272 ^ 2;
t235 = -pkin(8) * t399 + t265;
t229 = -pkin(8) * t362 + t258;
t224 = -pkin(8) * t363 + t257;
t219 = qJD(1) * t230;
t217 = pkin(10) * t367 + t274 * t254;
t216 = -pkin(10) * t397 + t427 * t254;
t193 = t427 * t233 - t344;
t135 = pkin(4) * t368 - t139;
t76 = -qJD(4) * t132 + t129 * t278 + t275 * t341;
t75 = qJD(4) * t133 + t129 * t275 - t278 * t341;
t74 = pkin(4) * t122 + pkin(11) * t120;
t63 = t272 * t364 + t289;
t44 = -pkin(4) * t156 - t56;
t37 = -pkin(4) * t178 - t39;
t35 = t308 * qJD(5) + t128 * t274 + t76 * t427;
t34 = t103 * qJD(5) - t128 * t427 + t76 * t274;
t28 = -pkin(4) * t298 - t30;
t23 = t274 * t74 + t427 * t30;
t22 = -t274 * t30 + t427 * t74;
t19 = t75 * pkin(4) - t76 * pkin(11) + t55;
t10 = -pkin(4) * t128 - t12;
t9 = pkin(11) * t128 + t11;
t4 = -t18 * qJD(5) + t427 * t19 - t274 * t9;
t3 = t17 * qJD(5) + t274 * t19 + t427 * t9;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t323, t388 * t436, t320 * t279 * t386, -0.2e1 * t323, -t320 * t362, 0, -t219 * t406 - t230 * t321 + t277 * t333, -t218 * t406 - t229 * t321 + t279 * t333, (t218 * t279 + t219 * t277 + (-t224 * t279 - t226 * t277) * qJD(2) + (t229 * t279 + t230 * t277 + (-t235 * t279 - t277 * t389) * qJD(2)) * qJD(1)) * t273, t218 * t389 - t219 * t235 - t224 * t230 + t226 * t229, -t123 * t179 + t129 * t156, t123 * t178 - t124 * t179 - t128 * t156 - t129 * t154, -t129 * t286 + t123 * t231 + (qJD(1) * t179 + t156) * t341, t128 * t154 + t403, t128 * t286 + t124 * t231 + (-qJD(1) * t178 - t154) * t341, (-qJD(2) * t286 - t231 * t379) * t272 * t399, -t63 * t286 - t49 * t231 + t131 * t154 + t127 * t124 + t125 * t178 + t115 * t128 + (qJD(1) * t96 + t82) * t341, t62 * t286 - t292 * t231 + t131 * t156 - t127 * t123 + t125 * t179 + t115 * t129 + (-qJD(1) * t97 - t83) * t341, t123 * t96 - t124 * t97 - t128 * t83 - t129 * t82 - t154 * t62 - t156 * t63 + t178 * t292 - t179 * t49, t115 * t131 + t125 * t127 - t292 * t97 + t49 * t96 + t62 * t83 + t63 * t82, t122 * t76 - t133 * t68, -t120 * t76 - t122 * t75 + t132 * t68 - t133 * t69, t122 * t128 + t133 * t124 - t68 * t178 + t298 * t76, t120 * t75 + t419, -t120 * t128 - t132 * t124 - t69 * t178 - t298 * t75, t128 * t298 + t403, t12 * t298 + t55 * t120 + t39 * t124 + t30 * t128 + t47 * t132 + t8 * t178 + t84 * t69 + t72 * t75, -t11 * t298 + t55 * t122 - t124 * t420 - t31 * t128 + t47 * t133 + t178 * t309 - t84 * t68 + t72 * t76, -t11 * t120 - t12 * t122 + t132 * t309 - t133 * t8 - t30 * t76 - t31 * t75 + t39 * t68 - t420 * t69, t11 * t31 + t12 * t30 - t309 * t420 + t39 * t8 + t47 * t84 + t55 * t72, -t103 * t26 + t35 * t90, -t103 * t27 - t26 * t308 - t34 * t90 - t35 * t88, t103 * t69 + t119 * t35 - t132 * t26 + t75 * t90, -t27 * t308 + t34 * t88, -t119 * t34 - t132 * t27 + t308 * t69 - t75 * t88, t119 * t75 + t419, t10 * t88 + t119 * t4 + t132 * t2 + t17 * t69 + t27 * t37 + t28 * t34 - t308 * t6 - t310 * t75, -t1 * t132 + t10 * t90 + t103 * t6 - t119 * t3 - t14 * t75 - t18 * t69 - t26 * t37 + t28 * t35, t1 * t308 - t103 * t2 - t14 * t34 + t17 * t26 - t18 * t27 - t3 * t88 + t310 * t35 - t4 * t90, t1 * t18 + t10 * t28 + t14 * t3 + t17 * t2 - t310 * t4 + t37 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t345, t388 * t401, -t279 * t334, t345, t277 * t334, 0, -t273 * pkin(8) * t355 + t226 * t321 + (-qJD(2) * t348 + t401) * t277 * pkin(1), pkin(1) * t373 + t224 * t321 - t218, 0, 0, -t123 * t400 + t447 * t156, t206 * t154 + t156 * t205 + (-t428 * t123 - t124 * t276 + (-t428 * t154 - t156 * t276) * qJD(3)) * t272, t268 * t276 * t339 - t123 * t405 + t206 * t286 + (-t156 * t363 - t428 * t283) * t272, -t154 * t439 - t272 * t366, t268 * t277 * t318 - t124 * t405 - t205 * t286 + (t154 * t363 + t276 * t283) * t272, -(qJD(1) * t329 - t380) * t343, t49 * t405 - t130 * t154 - t115 * t205 + (t115 * t385 - t125 * t428 - pkin(2) * t124 + (qJD(2) * t234 - t82) * t363) * t272 + t390 * t286, t292 * t405 - t130 * t156 - t115 * t206 + (t115 * t360 + pkin(2) * t123 + t125 * t276 + (-qJD(2) * t236 + t83) * t363) * t272 + t391 * t286, t234 * t123 - t236 * t124 + t83 * t205 + t82 * t206 + t390 * t156 - t391 * t154 + (-t428 * t292 - t276 * t49 + (-t276 * t83 - t428 * t82) * qJD(3)) * t272, -pkin(2) * t125 * t272 - t115 * t130 + t234 * t49 - t236 * t292 - t390 * t82 + t391 * t83, -t122 * t393 - t233 * t68, t120 * t393 - t122 * t392 + t232 * t68 - t233 * t69, -t122 * t439 + t233 * t124 - t298 * t393 + t368 * t68, t120 * t392 + t417, t120 * t439 - t232 * t124 - t298 * t392 + t368 * t69, -t298 * t205 + (t298 * t385 - t366) * t272, t139 * t124 + t221 * t69 + t47 * t232 - t30 * t205 + t392 * t72 + (t30 * t385 - t8 * t428) * t272 + t407 * t120 + t410 * t298, -t430 * t124 - t221 * t68 + t47 * t233 + t31 * t205 - t393 * t72 + (-t309 * t428 - t31 * t385) * t272 + t407 * t122 - t411 * t298, -t120 * t411 - t122 * t410 + t139 * t68 + t232 * t309 - t233 * t8 + t30 * t393 - t31 * t392 - t430 * t69, t139 * t8 + t221 * t47 + t30 * t410 - t309 * t430 + t31 * t411 + t407 * t72, -t193 * t26 - t395 * t90, t192 * t26 - t193 * t27 - t394 * t90 + t395 * t88, -t119 * t395 + t193 * t69 - t232 * t26 + t392 * t90, t192 * t27 + t394 * t88, -t119 * t394 - t192 * t69 - t232 * t27 - t392 * t88, t119 * t392 + t417, -t423 * t119 + t135 * t27 + t192 * t6 + t2 * t232 + t394 * t28 - t310 * t392 + t421 * t88 + t69 * t86, -t1 * t232 + t422 * t119 - t135 * t26 - t392 * t14 + t193 * t6 - t395 * t28 + t421 * t90 - t69 * t87, -t1 * t192 - t394 * t14 - t193 * t2 + t26 * t86 - t27 * t87 - t310 * t395 + t422 * t88 + t423 * t90, t1 * t87 + t135 * t6 - t422 * t14 + t2 * t86 + t421 * t28 + t310 * t423; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t402, -t154 ^ 2 + t156 ^ 2, -t154 * t286 - t123, -t402, -t156 * t286 - t387 * t429 - t291, t322, -t115 * t156 - t286 * t83 + t49, t115 * t154 - t286 * t82 + t292, 0, 0, -t275 * t68 + t278 * t431, (-t68 - t432) * t278 + (-t69 - t431) * t275, -t122 * t156 + t275 * t124 + t278 * t282, t120 * t346 - t413, t120 * t156 + t278 * t124 - t275 * t282, -t298 * t156, -pkin(3) * t69 - t83 * t120 - t30 * t156 + t307 * t275 - t441 * t278 - t56 * t298, pkin(3) * t68 - t83 * t122 + t31 * t156 + t441 * t275 + t307 * t278 + t57 * t298, t120 * t57 + t122 * t56 + ((-t69 + t384) * pkin(10) + t437) * t278 + (-t8 - t298 * t31 + (qJD(4) * t120 - t68) * pkin(10)) * t275, -pkin(3) * t47 - t30 * t56 - t31 * t57 - t72 * t83 + (-t275 * t8 - t278 * t309 + (-t275 * t31 - t278 * t30) * qJD(4)) * pkin(10), -t275 * t376 + (-t275 * t381 + t313) * t90, t90 * t106 + t107 * t88 + (-t274 * t90 - t427 * t88) * t382 + (-t375 + t416 + (-t427 * t90 + t414) * qJD(5)) * t275, t26 * t278 + t313 * t119 + (t298 * t90 - t301) * t275, t27 * t274 * t275 + (t275 * t358 - t330) * t88, t27 * t278 + t330 * t119 + (-t119 * t358 - t298 * t88 - t415) * t275, t119 * t346 - t413, -t28 * t106 + t216 * t69 - t44 * t88 + t408 * t119 + (-t2 + (pkin(10) * t88 + t274 * t28) * qJD(4)) * t278 + (pkin(10) * t27 + t28 * t358 - t298 * t310 + t425) * t275, -t28 * t107 - t217 * t69 - t44 * t90 - t409 * t119 + (t1 + (pkin(10) * t90 + t427 * t28) * qJD(4)) * t278 + (-pkin(10) * t26 - t14 * t298 - t28 * t381 + t377) * t275, t14 * t106 - t310 * t107 + t216 * t26 - t217 * t27 - t408 * t90 - t409 * t88 + t311 * t382 + (-t427 * t2 - t1 * t274 + (-t427 * t14 - t274 * t310) * qJD(5)) * t275, t1 * t217 + t2 * t216 - t28 * t44 + t409 * t14 - t408 * t310 + (t275 * t6 + t28 * t382) * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t404, -t120 ^ 2 + t122 ^ 2, -t68 + t432, -t404, t122 * t154 - t347, t124, -t72 * t122 + t154 * t31 - t354, t72 * t120 - t437, 0, 0, t300 * t90 - t416, -t376 - t300 * t88 + (-t27 - t412) * t274, t119 * t300 - t90 * t122 + t415, t119 * t414 - t375, -t119 ^ 2 * t274 + t88 * t122 + t374, -t119 * t122, -t377 - pkin(4) * t27 + t310 * t122 - t31 * t88 + (-pkin(11) * t358 - t22) * t119 + (-pkin(11) * t69 + t119 * t28) * t274, pkin(4) * t26 + pkin(11) * t301 + t23 * t119 + t14 * t122 + t28 * t300 - t31 * t90 + t425, t378 + t22 * t90 + t23 * t88 + t300 * t310 + (t358 * t90 - t375) * pkin(11) + ((qJD(5) * t88 - t26) * pkin(11) + t434) * t274, -t6 * pkin(4) + t310 * t22 - t14 * t23 - t28 * t31 + (qJD(5) * t311 - t2 * t274 + t378) * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t424, -t88 ^ 2 + t90 ^ 2, t119 * t88 - t26, -t424, -t27 + t412, t69, -t28 * t90 - t434, -t119 * t310 + t28 * t88 - t1, 0, 0;];
tauc_reg = t7;
