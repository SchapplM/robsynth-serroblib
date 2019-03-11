% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RRRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRRPR10_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR10_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:08:32
% EndTime: 2019-03-09 23:09:04
% DurationCPUTime: 13.68s
% Computational Cost: add. (23161->716), mult. (59569->928), div. (0->0), fcn. (47003->10), ass. (0->322)
t269 = cos(qJ(2));
t267 = sin(qJ(2));
t262 = sin(pkin(6));
t399 = qJD(1) * t262;
t372 = t267 * t399;
t263 = cos(pkin(6));
t398 = qJD(1) * t263;
t382 = pkin(1) * t398;
t211 = -pkin(8) * t372 + t269 * t382;
t309 = (pkin(2) * t267 - pkin(9) * t269) * t262;
t212 = qJD(1) * t309;
t266 = sin(qJ(3));
t447 = cos(qJ(3));
t154 = -t266 * t211 + t447 * t212;
t383 = t447 * pkin(9);
t241 = pkin(10) * t447 + t383;
t375 = t269 * t447;
t500 = (pkin(3) * t267 - pkin(10) * t375) * t399 + t154 + qJD(3) * t241;
t155 = t447 * t211 + t266 * t212;
t240 = (-pkin(10) - pkin(9)) * t266;
t397 = qJD(1) * t269;
t371 = t262 * t397;
t339 = t266 * t371;
t499 = -pkin(10) * t339 - qJD(3) * t240 + t155;
t265 = sin(qJ(4));
t446 = cos(qJ(4));
t226 = t265 * t266 - t446 * t447;
t334 = t447 * t399;
t315 = t269 * t334;
t457 = qJD(3) + qJD(4);
t402 = t226 * t457 - t265 * t339 + t315 * t446;
t227 = t265 * t447 + t266 * t446;
t173 = t457 * t227;
t180 = t227 * t371;
t401 = t173 - t180;
t214 = pkin(8) * t371 + t267 * t382;
t393 = qJD(3) * t266;
t331 = -t214 + (-t339 + t393) * pkin(3);
t344 = qJD(2) + t398;
t301 = t447 * t344;
t336 = t266 * t372;
t195 = t301 - t336;
t196 = t266 * t344 + t267 * t334;
t144 = -t446 * t195 + t196 * t265;
t239 = -qJD(3) + t371;
t234 = -qJD(4) + t239;
t394 = qJD(2) * t267;
t370 = t262 * t394;
t242 = qJD(1) * t370;
t264 = sin(qJ(6));
t268 = cos(qJ(6));
t389 = qJD(6) * t268;
t390 = qJD(6) * t264;
t232 = qJD(3) * t336;
t292 = qJD(2) * t315 - t232;
t276 = -qJD(3) * t301 - t292;
t306 = t265 * t195 + t196 * t446;
t385 = qJD(1) * qJD(2);
t362 = t269 * t385;
t333 = t262 * t362;
t377 = t196 * qJD(3) + t266 * t333;
t84 = qJD(4) * t306 - t265 * t276 + t446 * t377;
t46 = -t144 * t389 - t234 * t390 - t268 * t242 - t264 * t84;
t121 = -t268 * t144 - t234 * t264;
t458 = qJD(6) + t306;
t481 = t458 * t121;
t498 = t46 - t481;
t487 = t234 * t306 + t84;
t184 = t265 * t240 + t241 * t446;
t427 = qJD(4) * t184 - t499 * t265 + t500 * t446;
t482 = t268 * t458;
t365 = qJD(4) * t446;
t392 = qJD(4) * t265;
t83 = -t195 * t365 + t196 * t392 + t265 * t377 + t446 * t276;
t293 = t264 * t83 - t458 * t482;
t497 = -t121 * t144 + t293;
t123 = t144 * t264 - t234 * t268;
t483 = t264 * t458;
t308 = -t268 * t83 - t458 * t483;
t496 = t123 * t144 + t308;
t495 = -t402 * qJ(5) + qJD(5) * t227 - t331;
t43 = t268 * t46;
t494 = -t123 * t483 - t43;
t391 = qJD(6) * t123;
t47 = t242 * t264 - t268 * t84 + t391;
t436 = t264 * t47;
t493 = t268 * t481 + t436;
t492 = -t268 * t47 - t123 * t482 + (t46 + t481) * t264;
t445 = t144 * pkin(5);
t449 = pkin(4) + pkin(11);
t491 = t401 * t449 - t495;
t178 = -pkin(2) * t344 - t211;
t149 = -t195 * pkin(3) + t178;
t277 = -qJ(5) * t306 + t149;
t68 = t144 * pkin(4) + t277;
t490 = t144 * t68;
t179 = pkin(9) * t344 + t214;
t209 = (-pkin(2) * t269 - pkin(9) * t267 - pkin(1)) * t262;
t191 = qJD(1) * t209;
t134 = t179 * t447 + t266 * t191;
t116 = t195 * pkin(10) + t134;
t114 = t446 * t116;
t133 = -t179 * t266 + t447 * t191;
t115 = -pkin(10) * t196 + t133;
t291 = -pkin(3) * t239 + t115;
t285 = t265 * t291;
t54 = t114 + t285;
t51 = t234 * qJ(5) - t54;
t39 = -t51 - t445;
t489 = t39 * t458;
t410 = t262 * t267;
t342 = t449 * t410;
t488 = -pkin(5) * t402 + qJD(1) * t342 + t427;
t429 = -t240 * t365 + t241 * t392 + t500 * t265 + t499 * t446;
t424 = qJ(5) * t144;
t484 = t144 * t306;
t480 = t458 * t144;
t105 = t446 * t291;
t407 = t265 * t116;
t53 = -t105 + t407;
t405 = -qJD(5) - t53;
t450 = t306 ^ 2;
t478 = -t144 ^ 2 + t450;
t467 = pkin(5) * t306;
t404 = t467 - t405;
t37 = t234 * t449 + t404;
t48 = t144 * t449 + t277;
t19 = t264 * t37 + t268 * t48;
t236 = qJ(5) * t242;
t213 = qJD(2) * t309;
t203 = qJD(1) * t213;
t249 = pkin(8) * t410;
t408 = t263 * t269;
t221 = pkin(1) * t408 - t249;
t215 = t221 * qJD(2);
t204 = qJD(1) * t215;
t357 = t447 * t203 - t266 * t204;
t67 = -t292 * pkin(10) + pkin(3) * t242 + (-pkin(10) * t301 - t134) * qJD(3) + t357;
t366 = qJD(3) * t447;
t96 = -t179 * t393 + t191 * t366 + t266 * t203 + t447 * t204;
t73 = -pkin(10) * t377 + t96;
t349 = -qJD(4) * t105 + t116 * t392 - t265 * t67 - t446 * t73;
t328 = qJD(5) * t234 + t349;
t12 = -t236 + t328;
t7 = -pkin(5) * t84 - t12;
t477 = -t144 * t19 + t7 * t268;
t252 = t263 * t267 * pkin(1);
t409 = t262 * t269;
t222 = pkin(8) * t409 + t252;
t216 = t222 * qJD(2);
t205 = qJD(1) * t216;
t132 = pkin(3) * t377 + t205;
t273 = t83 * qJ(5) - qJD(5) * t306 + t132;
t13 = t449 * t84 + t273;
t16 = -qJD(4) * t285 - t116 * t365 - t265 * t73 + t446 * t67;
t319 = qJD(2) * t342;
t9 = -pkin(5) * t83 - qJD(1) * t319 - t16;
t2 = -qJD(6) * t19 - t13 * t264 + t268 * t9;
t476 = t19 * t458 + t2;
t475 = t123 * t458 - t47;
t474 = t144 * t234 + t83;
t473 = t144 * t149 + t349;
t321 = t264 * t48 - t268 * t37;
t472 = -t144 * t321 + t7 * t264 + t39 * t389;
t376 = t267 * t447;
t217 = t262 * t376 + t263 * t266;
t379 = t266 * t410;
t297 = -t263 * t447 + t379;
t282 = t446 * t297;
t164 = t217 * t265 + t282;
t288 = t217 * qJD(3);
t369 = qJD(2) * t409;
t279 = t266 * t369 + t288;
t289 = t265 * t297;
t364 = t447 * qJD(2);
t453 = qJD(3) * t297 - t364 * t409;
t99 = -qJD(4) * t289 + t217 * t365 - t265 * t453 + t279 * t446;
t471 = -t234 * t99 + t262 * ((qJD(1) * t164 + t144) * t394 - t269 * t84);
t165 = t217 * t446 - t289;
t98 = qJD(4) * t282 + t217 * t392 + t265 * t279 + t446 * t453;
t470 = t234 * t98 + t262 * ((qJD(1) * t165 + t306) * t394 + t269 * t83);
t259 = t262 ^ 2;
t469 = -0.2e1 * t259 * t385;
t468 = pkin(4) * t306;
t466 = t306 * t39;
t441 = t306 * t68;
t430 = qJ(5) * t372 + t429;
t61 = t115 * t446 - t407;
t426 = -pkin(3) * t365 - qJD(5) + t61;
t58 = t83 * t165;
t465 = -t306 * t98 - t58;
t74 = t83 * t227;
t464 = -t402 * t306 - t74;
t97 = -qJD(3) * t134 + t357;
t463 = t134 * t239 - t97;
t461 = t149 * t306;
t459 = t306 * t449;
t208 = pkin(9) * t263 + t222;
t153 = t447 * t208 + t266 * t209;
t1 = -qJD(6) * t321 + t13 * t268 + t264 * t9;
t367 = t306 * t321 + t1;
t456 = -t306 * t19 - t2;
t152 = -t266 * t208 + t447 * t209;
t120 = -pkin(3) * t409 - pkin(10) * t217 + t152;
t128 = -pkin(10) * t297 + t153;
t356 = t447 * t213 - t266 * t215;
t85 = pkin(3) * t370 + pkin(10) * t453 - t208 * t366 - t209 * t393 + t356;
t106 = -t208 * t393 + t209 * t366 + t266 * t213 + t447 * t215;
t91 = -pkin(10) * t279 + t106;
t302 = -t120 * t365 + t128 * t392 - t265 * t85 - t446 * t91;
t20 = -t262 * (qJ(5) * t394 - qJD(5) * t269) + t302;
t452 = t234 * t402 + (qJD(2) * t227 - t306) * t372;
t451 = t234 * t401 - (qJD(2) * t226 - t144) * t372;
t271 = qJD(1) ^ 2;
t444 = t269 * pkin(1);
t257 = -pkin(3) * t447 - pkin(2);
t299 = -t227 * qJ(5) + t257;
t150 = t226 * t449 + t299;
t183 = -t446 * t240 + t241 * t265;
t159 = pkin(5) * t227 + t183;
t103 = -t150 * t264 + t159 * t268;
t443 = qJD(6) * t103 + t488 * t264 + t491 * t268;
t104 = t150 * t268 + t159 * t264;
t442 = -qJD(6) * t104 - t491 * t264 + t488 * t268;
t256 = -pkin(3) * t446 - pkin(4);
t253 = -pkin(11) + t256;
t437 = t253 * t83;
t432 = t449 * t83;
t431 = -t401 * pkin(5) - t430;
t428 = pkin(4) * t372 + t427;
t425 = -t426 + t467;
t423 = t123 * t121;
t417 = t196 * t195;
t416 = t196 * t239;
t415 = t226 * t264;
t414 = t226 * t268;
t412 = t239 * t266;
t411 = t259 * t271;
t403 = -t401 * pkin(4) + t495;
t76 = t265 * t120 + t446 * t128;
t400 = t267 ^ 2 - t269 ^ 2;
t396 = qJD(2) * t183;
t395 = qJD(2) * t184;
t388 = qJD(6) * t449;
t381 = pkin(3) * t392;
t380 = t96 * t447;
t378 = t268 * t409;
t374 = t447 * t178;
t373 = t259 * t397;
t360 = -t47 + t391;
t60 = t115 * t265 + t114;
t157 = -t268 * t180 + t264 * t372;
t359 = t173 * t268 + t157;
t158 = t180 * t264 + t268 * t372;
t358 = -t173 * t264 + t158;
t343 = qJD(2) + 0.2e1 * t398;
t341 = t269 * t267 * t411;
t338 = t234 * t372;
t332 = -t60 + t381;
t330 = pkin(1) * t469;
t75 = t120 * t446 - t265 * t128;
t326 = pkin(3) * t196 + t424;
t325 = t144 * t99 + t164 * t84;
t324 = t19 * t264 - t268 * t321;
t323 = t19 * t268 + t264 * t321;
t322 = -t183 * t83 - t184 * t84;
t70 = pkin(4) * t409 - t75;
t49 = t165 * pkin(5) + pkin(11) * t409 + t70;
t166 = pkin(3) * t379 + t249 + (t257 - t444) * t263;
t275 = -t165 * qJ(5) + t166;
t63 = t164 * t449 + t275;
t30 = -t264 * t63 + t268 * t49;
t31 = t264 * t49 + t268 * t63;
t320 = pkin(4) * t242;
t314 = t259 * t267 * t362;
t69 = qJ(5) * t409 - t76;
t311 = t377 * t447;
t140 = t164 * t268 + t264 * t409;
t307 = -t234 * t54 + t16;
t304 = t120 * t392 + t128 * t365 + t265 * t91 - t446 * t85;
t300 = t144 * t401 + t226 * t84;
t296 = t226 * t390 - t359;
t295 = t226 * t389 - t358;
t294 = pkin(1) * (-t263 * t385 + t411);
t14 = -t320 - t16;
t284 = t315 - t366;
t283 = t144 * t98 + t164 * t83 - t165 * t84 - t306 * t99;
t281 = qJD(6) * t323 + t1 * t264 + t2 * t268;
t280 = t144 * t402 + t226 * t83 - t227 * t84 - t306 * t401;
t34 = pkin(3) * t279 + t99 * pkin(4) + pkin(8) * t369 + t98 * qJ(5) + qJD(2) * t252 - t165 * qJD(5);
t254 = pkin(3) * t265 + qJ(5);
t207 = t249 + (-pkin(2) - t444) * t263;
t171 = (-t234 * t262 - t373) * t394;
t167 = t226 * pkin(4) + t299;
t160 = -t226 * pkin(5) + t184;
t151 = pkin(3) * t288 + (t252 + (pkin(3) * t266 + pkin(8)) * t409) * qJD(2);
t141 = -t264 * t164 + t378;
t107 = -qJD(3) * t153 + t356;
t100 = t424 + t468;
t92 = t164 * pkin(4) + t275;
t86 = t326 + t468;
t66 = t424 + t459;
t57 = qJD(6) * t140 + t264 * t99 + t268 * t370;
t56 = -qJD(6) * t378 - t268 * t99 + (qJD(6) * t164 + t370) * t264;
t55 = t326 + t459;
t52 = -pkin(5) * t164 - t69;
t50 = pkin(4) * t234 - t405;
t44 = t60 - t445;
t41 = t54 - t445;
t29 = t84 * pkin(4) + t273;
t28 = t264 * t41 + t268 * t66;
t27 = -t264 * t66 + t268 * t41;
t24 = t264 * t44 + t268 * t55;
t23 = -t264 * t55 + t268 * t44;
t22 = -pkin(4) * t370 + t304;
t21 = t99 * pkin(11) + t34;
t17 = t321 * t390;
t11 = -pkin(5) * t99 - t20;
t10 = -t98 * pkin(5) + t304 - t319;
t4 = -qJD(6) * t31 + t10 * t268 - t21 * t264;
t3 = qJD(6) * t30 + t10 * t264 + t21 * t268;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t314, t400 * t469, t343 * t369, -0.2e1 * t314, -t343 * t370, 0, -t205 * t263 - t216 * t344 + t267 * t330, -t204 * t263 - t215 * t344 + t269 * t330 (t204 * t269 + t205 * t267 + (-t211 * t269 - t214 * t267) * qJD(2) + (t215 * t269 + t216 * t267 + (-t221 * t269 - t222 * t267) * qJD(2)) * qJD(1)) * t262, t204 * t222 - t205 * t221 - t211 * t216 + t214 * t215, -t196 * t453 - t217 * t276, -t195 * t453 - t196 * t279 - t217 * t377 + t276 * t297, t196 * t370 + t217 * t242 + t239 * t453 + t276 * t409, -t195 * t279 + t297 * t377, t195 * t370 + t239 * t279 - t242 * t297 + t377 * t409 (-t239 * t262 - t373) * t394, -t107 * t239 + t133 * t370 + t152 * t242 + t178 * t279 - t216 * t195 + t205 * t297 + t207 * t377 - t97 * t409, t106 * t239 - t134 * t370 - t153 * t242 - t178 * t453 + t216 * t196 + t205 * t217 - t207 * t276 + t409 * t96, t106 * t195 - t107 * t196 + t133 * t453 - t134 * t279 + t152 * t276 - t153 * t377 - t97 * t217 - t96 * t297, t106 * t134 + t107 * t133 + t152 * t97 + t153 * t96 + t178 * t216 + t205 * t207, t465, t283, t470, t325, -t471, t171, t132 * t164 + t144 * t151 + t149 * t99 + t166 * t84 + t234 * t304 + (-t16 * t269 + (qJD(1) * t75 - t53) * t394) * t262, t132 * t165 + t306 * t151 - t149 * t98 - t166 * t83 - t234 * t302 + (-t349 * t269 + (-qJD(1) * t76 - t54) * t394) * t262, t144 * t302 - t16 * t165 + t164 * t349 + t304 * t306 - t53 * t98 - t54 * t99 + t75 * t83 - t76 * t84, t132 * t166 + t149 * t151 + t16 * t75 - t302 * t54 + t304 * t53 - t349 * t76, t171, -t470, t471, t465, t283, t325, t12 * t164 + t14 * t165 + t144 * t20 + t22 * t306 - t50 * t98 + t51 * t99 + t69 * t84 - t70 * t83, -t144 * t34 - t164 * t29 - t22 * t234 - t68 * t99 - t84 * t92 + (-t14 * t269 + (qJD(1) * t70 + t50) * t394) * t262, -t306 * t34 - t165 * t29 + t20 * t234 + t68 * t98 + t83 * t92 + (t12 * t269 + (-qJD(1) * t69 - t51) * t394) * t262, t12 * t69 + t14 * t70 + t20 * t51 + t22 * t50 + t29 * t92 + t34 * t68, t123 * t57 + t141 * t46, -t121 * t57 - t123 * t56 - t140 * t46 + t141 * t47, -t123 * t98 + t141 * t83 - t165 * t46 + t458 * t57, t121 * t56 - t140 * t47, t121 * t98 - t140 * t83 - t165 * t47 - t458 * t56, -t458 * t98 - t58, t11 * t121 - t140 * t7 + t165 * t2 - t30 * t83 + t321 * t98 + t39 * t56 + t4 * t458 + t47 * t52, -t1 * t165 + t11 * t123 - t141 * t7 + t19 * t98 - t3 * t458 + t31 * t83 + t39 * t57 - t46 * t52, t1 * t140 - t121 * t3 - t123 * t4 + t141 * t2 - t19 * t56 + t30 * t46 - t31 * t47 + t321 * t57, t1 * t31 + t11 * t39 + t19 * t3 + t2 * t30 - t321 * t4 + t52 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t341, t400 * t411, -t262 * t271 * t408, t341, t344 * t372 - t242, 0, -pkin(8) * t333 + t214 * t344 + t267 * t294, pkin(8) * t242 + t211 * t344 + t269 * t294, 0, 0, -t196 * t284 - t266 * t276, -t276 * t447 - t284 * t195 + (-t377 + t416) * t266, -t239 * t366 + (t239 * t375 + (qJD(2) * t266 - t196) * t267) * t399, t195 * t412 - t311, t239 * t393 + (-t269 * t412 + (t364 - t195) * t267) * t399, t239 * t372, -pkin(2) * t377 - t205 * t447 + t154 * t239 + t214 * t195 + (t178 * t266 + t239 * t383) * qJD(3) + (-t133 * t267 + (-pkin(9) * t394 - t178 * t269) * t266) * t399, pkin(2) * t232 - t155 * t239 - t214 * t196 + t205 * t266 + (-pkin(2) * t301 - pkin(9) * t412 + t374) * qJD(3) + (-t269 * t374 + t134 * t267 + (-pkin(2) * t375 - pkin(9) * t376) * qJD(2)) * t399, t380 + t154 * t196 - t155 * t195 + t284 * t133 + (t196 * t366 - t311) * pkin(9) + (((-t195 + t301) * qJD(3) + t292) * pkin(9) + t463) * t266, -t205 * pkin(2) - t133 * t154 - t134 * t155 - t178 * t214 + (t380 - t97 * t266 + (-t133 * t447 - t134 * t266) * qJD(3)) * pkin(9), t464, t280, t452, t300, t451, t338, t132 * t226 + t257 * t84 + t427 * t234 + t401 * t149 + t331 * t144 + (t53 - t396) * t372, t132 * t227 - t257 * t83 - t429 * t234 - t402 * t149 + t331 * t306 + (t54 - t395) * t372, t144 * t429 - t16 * t227 + t226 * t349 + t306 * t427 - t401 * t54 - t402 * t53 + t322, t132 * t257 + t149 * t331 - t16 * t183 - t184 * t349 + t427 * t53 - t429 * t54, t338, -t452, -t451, t464, t280, t300, t12 * t226 + t14 * t227 + t144 * t430 + t306 * t428 + t401 * t51 - t402 * t50 + t322, -t167 * t84 - t226 * t29 - t401 * t68 - t428 * t234 + t403 * t144 + (-t50 + t396) * t372, t167 * t83 - t227 * t29 + t402 * t68 + t430 * t234 + t403 * t306 + (t51 + t395) * t372, -t12 * t184 + t14 * t183 + t167 * t29 - t403 * t68 + t428 * t50 + t430 * t51, t123 * t295 - t415 * t46, t359 * t123 + t358 * t121 + (-t436 - t43 + (-t121 * t268 - t123 * t264) * qJD(6)) * t226, -t123 * t402 - t227 * t46 + t295 * t458 - t415 * t83, t121 * t296 - t414 * t47, t121 * t402 - t227 * t47 - t296 * t458 - t414 * t83, -t402 * t458 - t74, -t103 * t83 + t121 * t431 + t160 * t47 + t2 * t227 + t296 * t39 + t321 * t402 - t414 * t7 + t442 * t458, -t1 * t227 + t104 * t83 + t123 * t431 - t160 * t46 + t19 * t402 + t295 * t39 + t415 * t7 - t443 * t458, t103 * t46 - t104 * t47 + t157 * t19 - t158 * t321 + t323 * t173 - t442 * t123 - t443 * t121 + (-qJD(6) * t324 + t1 * t268 - t2 * t264) * t226, t1 * t104 + t103 * t2 + t160 * t7 + t19 * t443 - t321 * t442 + t39 * t431; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t417, -t195 ^ 2 + t196 ^ 2, t195 * t239 - t276, t417, -t377 - t416, t242, -t178 * t196 - t463, -t133 * t239 - t178 * t195 - t96, 0, 0, t484, t478, -t474, -t484, -t487, t242, -t461 - t60 * t234 + (-t144 * t196 + t234 * t392 + t242 * t446) * pkin(3) + t16, -t61 * t234 + (-t196 * t306 + t234 * t365 - t242 * t265) * pkin(3) + t473, t54 * t306 + t61 * t144 + t144 * t53 - t60 * t306 + (t446 * t83 - t265 * t84 + (-t144 * t446 + t265 * t306) * qJD(4)) * pkin(3), -t53 * t60 - t54 * t61 + (t446 * t16 - t149 * t196 - t349 * t265 + (t265 * t53 + t446 * t54) * qJD(4)) * pkin(3), t242, t474, t487, t484, t478, -t484, -t254 * t84 - t256 * t83 + (t332 - t51) * t306 + (t426 + t50) * t144, t441 + t144 * t86 - t332 * t234 + (-pkin(4) + t256) * t242 - t16, t234 * t426 + t242 * t254 + t306 * t86 - t12 - t490, -t12 * t254 + t14 * t256 + t332 * t50 + t426 * t51 - t68 * t86, t494, t492, t496, t493, t497, t480, t254 * t47 + (-t437 + t466) * t268 + t425 * t121 + (-t253 * t390 + t268 * t381 - t23) * t458 + t472, -t254 * t46 + (-t253 * t389 + t24) * t458 + t425 * t123 + (-t381 * t458 + t437 - t489) * t264 + t477, t121 * t24 + t123 * t23 - t17 + (-t121 * t381 + t253 * t360 - t367) * t264 + (-t123 * t381 + t253 * t46 + (-t121 * t253 - t19) * qJD(6) + t456) * t268, -t19 * t24 + t23 * t321 + t253 * t281 + t254 * t7 + t324 * t381 + t39 * t425; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t484, t478, -t474, -t484, -t487, t242, t307 - t461, t234 * t53 + t473, 0, 0, t242, t474, t487, t484, t478, -t484, pkin(4) * t83 - qJ(5) * t84 + (-t51 - t54) * t306 + (t50 + t405) * t144, t100 * t144 - t307 - 0.2e1 * t320 + t441, t100 * t306 + t234 * t405 + 0.2e1 * t236 - t328 - t490, -pkin(4) * t14 - qJ(5) * t12 - t100 * t68 + t405 * t51 - t50 * t54, t494, t492, t496, t493, t497, t480, qJ(5) * t47 + (t432 + t466) * t268 + (t264 * t388 - t27) * t458 + t404 * t121 + t472, -qJ(5) * t46 + (t268 * t388 + t28) * t458 + t404 * t123 + (-t432 - t489) * t264 + t477, t121 * t28 + t123 * t27 - t17 + (-t360 * t449 - t367) * t264 + (-t449 * t46 + (t121 * t449 - t19) * qJD(6) + t456) * t268, qJ(5) * t7 - t19 * t28 + t27 * t321 - t281 * t449 + t39 * t404; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t474, t242 - t484, -t234 ^ 2 - t450, -t234 * t51 + t14 + t441, 0, 0, 0, 0, 0, 0, t121 * t234 + t308, t123 * t234 + t293, t475 * t264 + t498 * t268, t234 * t39 + t367 * t264 + t268 * t476 + t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t423, -t121 ^ 2 + t123 ^ 2, -t498, -t423, t475, -t83, -t123 * t39 + t476, t121 * t39 - t321 * t458 - t1, 0, 0;];
tauc_reg  = t5;
