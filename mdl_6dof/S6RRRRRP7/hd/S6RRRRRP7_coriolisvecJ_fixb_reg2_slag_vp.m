% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RRRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRRRP7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:43:10
% EndTime: 2019-03-10 01:43:47
% DurationCPUTime: 16.05s
% Computational Cost: add. (29089->748), mult. (74670->990), div. (0->0), fcn. (59121->10), ass. (0->324)
t333 = sin(pkin(6));
t338 = sin(qJ(2));
t427 = qJD(1) * t338;
t403 = t333 * t427;
t334 = cos(pkin(6));
t340 = cos(qJ(2));
t487 = pkin(1) * t340;
t416 = t334 * t487;
t271 = -pkin(8) * t403 + qJD(1) * t416;
t358 = (pkin(2) * t338 - pkin(9) * t340) * t333;
t272 = qJD(1) * t358;
t337 = sin(qJ(3));
t489 = cos(qJ(3));
t211 = -t337 * t271 + t489 * t272;
t414 = t489 * pkin(9);
t311 = pkin(10) * t489 + t414;
t405 = t340 * t489;
t429 = qJD(1) * t333;
t539 = (pkin(3) * t338 - pkin(10) * t405) * t429 + t211 + qJD(3) * t311;
t212 = t489 * t271 + t337 * t272;
t310 = (-pkin(10) - pkin(9)) * t337;
t426 = qJD(1) * t340;
t402 = t333 * t426;
t377 = t337 * t402;
t538 = -pkin(10) * t377 - qJD(3) * t310 + t212;
t336 = sin(qJ(4));
t488 = cos(qJ(4));
t352 = -t336 * t337 + t488 * t489;
t417 = qJD(3) + qJD(4);
t230 = t352 * t417;
t241 = t352 * t402;
t537 = t230 - t241;
t289 = t336 * t489 + t337 * t488;
t432 = (-t402 + t417) * t289;
t339 = cos(qJ(5));
t420 = qJD(5) * t339;
t428 = qJD(1) * t334;
t317 = qJD(2) + t428;
t406 = t333 * t489;
t374 = qJD(1) * t406;
t364 = t338 * t374;
t256 = t337 * t317 + t364;
t376 = t337 * t403;
t255 = t317 * t489 - t376;
t502 = t488 * t255;
t203 = t256 * t336 - t502;
t517 = t203 * t339;
t536 = t420 + t517;
t335 = sin(qJ(5));
t421 = qJD(5) * t335;
t518 = t203 * t335;
t535 = t421 + t518;
t396 = qJD(4) * t488;
t423 = qJD(4) * t336;
t435 = -t310 * t396 + t311 * t423 + t336 * t539 + t538 * t488;
t441 = t333 * t340;
t431 = t334 * t338 * pkin(1) + pkin(8) * t441;
t274 = t431 * qJD(1);
t424 = qJD(3) * t337;
t371 = -t274 + (-t377 + t424) * pkin(3);
t307 = -qJD(3) + t402;
t360 = -qJD(4) + t307;
t351 = t339 * t360;
t357 = t336 * t255 + t256 * t488;
t177 = t335 * t357 + t351;
t179 = -t335 * t360 + t339 * t357;
t516 = qJD(5) + t203;
t530 = t516 * t335;
t425 = qJD(2) * t338;
t401 = t333 * t425;
t312 = qJD(1) * t401;
t365 = t340 * t374;
t397 = qJD(3) * t489;
t223 = -qJD(2) * t365 + qJD(3) * t376 - t317 * t397;
t418 = qJD(1) * qJD(2);
t393 = t340 * t418;
t373 = t333 * t393;
t408 = qJD(3) * t364 + t317 * t424 + t337 * t373;
t346 = t223 * t488 + t256 * t423 + t336 * t408;
t343 = -t255 * t396 + t346;
t87 = qJD(5) * t351 - t335 * t312 + t339 * t343 + t357 * t421;
t422 = qJD(5) * t179;
t88 = -t339 * t312 - t335 * t343 + t422;
t534 = -t177 * t536 - t179 * t530 - t335 * t88 - t87 * t339;
t239 = pkin(9) * t317 + t274;
t269 = (-pkin(2) * t340 - pkin(9) * t338 - pkin(1)) * t333;
t251 = qJD(1) * t269;
t194 = t239 * t489 + t337 * t251;
t273 = qJD(2) * t358;
t263 = qJD(1) * t273;
t442 = t333 * t338;
t318 = pkin(8) * t442;
t283 = -t318 + t416;
t275 = t283 * qJD(2);
t264 = qJD(1) * t275;
t145 = -qJD(3) * t194 + t489 * t263 - t337 * t264;
t111 = pkin(3) * t312 + t223 * pkin(10) + t145;
t144 = -t239 * t424 + t251 * t397 + t337 * t263 + t489 * t264;
t118 = -pkin(10) * t408 + t144;
t193 = -t239 * t337 + t489 * t251;
t165 = -pkin(10) * t256 + t193;
t154 = -pkin(3) * t307 + t165;
t166 = t255 * pkin(10) + t194;
t378 = -t336 * t111 - t488 * t118 - t154 * t396 + t166 * t423;
t28 = pkin(11) * t312 - t378;
t238 = -pkin(2) * t317 - t271;
t207 = -pkin(3) * t255 + t238;
t112 = pkin(4) * t203 - pkin(11) * t357 + t207;
t164 = t488 * t166;
t100 = t336 * t154 + t164;
t97 = -pkin(11) * t360 + t100;
t53 = t112 * t335 + t339 * t97;
t385 = t336 * t223 - t488 * t408;
t128 = qJD(4) * t357 - t385;
t276 = t431 * qJD(2);
t265 = qJD(1) * t276;
t192 = pkin(3) * t408 + t265;
t61 = t128 * pkin(4) + pkin(11) * t343 + t192;
t8 = -qJD(5) * t53 - t28 * t335 + t339 * t61;
t344 = qJ(6) * t87 + t8;
t486 = pkin(5) * t128;
t1 = -qJD(6) * t179 + t344 + t486;
t388 = -t112 * t420 - t339 * t28 - t335 * t61 + t97 * t421;
t359 = qJ(6) * t88 + t388;
t3 = -qJD(6) * t177 - t359;
t52 = t339 * t112 - t335 * t97;
t33 = -qJ(6) * t179 + t52;
t32 = pkin(5) * t516 + t33;
t473 = t32 * t339;
t34 = -qJ(6) * t177 + t53;
t476 = t516 * t34;
t533 = t3 * t339 + t335 * (-t1 - t476) - t516 * t473;
t124 = t339 * t128;
t532 = t177 * t357 - t516 * t530 + t124;
t529 = -pkin(11) * t403 - t435;
t528 = -t432 * pkin(4) + pkin(11) * t537 - t371;
t214 = t241 * t335 - t339 * t403;
t527 = -t289 * t420 + t214;
t526 = t535 * pkin(5);
t31 = t488 * t111 - t336 * t118 - t154 * t423 - t166 * t396;
t29 = -pkin(4) * t312 - t31;
t20 = pkin(5) * t88 + t29;
t389 = t177 * pkin(5) + qJD(6);
t163 = t336 * t166;
t99 = t154 * t488 - t163;
t96 = pkin(4) * t360 - t99;
t78 = t389 + t96;
t524 = t20 * t335 + t34 * t357 + t536 * t78;
t523 = -t20 * t339 - t32 * t357 + t535 * t78;
t522 = -qJ(6) * t518 + t339 * qJD(6);
t84 = t87 * t335;
t520 = t179 * t536 - t84;
t122 = t335 * t128;
t500 = -t420 * t516 - t122;
t519 = -t179 * t357 + t516 * t517 - t500;
t450 = t203 * t357;
t244 = t336 * t310 + t311 * t488;
t434 = -qJD(4) * t244 + t538 * t336 - t488 * t539;
t449 = t230 * t335;
t354 = t449 - t527;
t515 = -t203 ^ 2 + t357 ^ 2;
t150 = pkin(4) * t357 + pkin(11) * t203;
t514 = t207 * t203 + t378;
t329 = t339 * qJ(6);
t513 = -pkin(5) * t357 - t203 * t329;
t512 = -t203 * t360 - t343;
t330 = t333 ^ 2;
t510 = -0.2e1 * t330 * t418;
t507 = -t516 * t53 - t8;
t327 = -pkin(3) * t489 - pkin(2);
t224 = -pkin(4) * t352 - t289 * pkin(11) + t327;
t506 = t224 * t420 - t335 * t528 + t339 * t529;
t380 = pkin(3) * t396;
t106 = t165 * t488 - t163;
t132 = pkin(3) * t256 + t150;
t65 = t339 * t106 + t335 * t132;
t505 = -t339 * t380 + t65;
t504 = t335 * t529 + t339 * t528;
t503 = t516 * t357;
t105 = t165 * t336 + t164;
t370 = pkin(3) * t423 - t105;
t268 = pkin(9) * t334 + t431;
t209 = -t268 * t337 + t489 * t269;
t280 = t334 * t337 + t338 * t406;
t174 = -pkin(3) * t441 - pkin(10) * t280 + t209;
t210 = t489 * t268 + t337 * t269;
t279 = -t334 * t489 + t337 * t442;
t187 = -pkin(10) * t279 + t210;
t120 = t336 * t174 + t488 * t187;
t115 = -pkin(11) * t441 + t120;
t217 = t279 * t488 + t280 * t336;
t218 = -t336 * t279 + t280 * t488;
t267 = t318 + (-pkin(2) - t487) * t334;
t222 = pkin(3) * t279 + t267;
t141 = pkin(4) * t217 - pkin(11) * t218 + t222;
t70 = t339 * t115 + t335 * t141;
t436 = pkin(4) * t403 - t434;
t501 = t194 * t307 - t145;
t235 = t339 * t244;
t185 = t335 * t224 + t235;
t499 = -t335 * t52 + t339 * t53;
t498 = -t207 * t357 + t31;
t495 = -t29 * t339 - t52 * t357 + t96 * t421;
t494 = t29 * t335 + t53 * t357 + t96 * t420;
t491 = -t307 * t357 + t385;
t490 = t179 ^ 2;
t485 = t339 * pkin(5);
t6 = t388 * t339;
t484 = -qJ(6) - pkin(11);
t481 = t32 - t33;
t215 = t241 * t339 + t335 * t403;
t362 = -qJ(6) * t230 - qJD(6) * t289;
t480 = t362 * t339 + (-t235 + (qJ(6) * t289 - t224) * t335) * qJD(5) + qJ(6) * t215 - t504 + t432 * pkin(5);
t479 = (-qJD(5) * t244 + t362) * t335 + t506 + t527 * qJ(6);
t478 = t185 * qJD(5) + t504;
t477 = t244 * t421 - t506;
t475 = t516 * t52;
t86 = t88 * t339;
t324 = pkin(3) * t336 + pkin(11);
t439 = -qJ(6) - t324;
t384 = qJD(5) * t439;
t469 = t335 * t384 - t505 + t522;
t64 = -t106 * t335 + t339 * t132;
t468 = (-t380 - qJD(6)) * t335 + t339 * t384 - t64 + t513;
t390 = qJD(5) * t484;
t68 = t335 * t150 + t339 * t99;
t467 = t335 * t390 + t522 - t68;
t67 = t339 * t150 - t335 * t99;
t466 = -t335 * qJD(6) + t339 * t390 + t513 - t67;
t465 = t370 + t526;
t464 = t128 * t217;
t463 = t128 * t352;
t462 = t177 * t516;
t461 = t177 * t335;
t460 = t179 * t177;
t459 = t179 * t516;
t448 = t256 * t255;
t447 = t256 * t307;
t446 = t289 * t335;
t445 = t289 * t339;
t444 = t307 * t337;
t443 = t330 * qJD(1) ^ 2;
t438 = t354 * pkin(5) + t436;
t430 = t338 ^ 2 - t340 ^ 2;
t410 = t335 * t441;
t407 = t144 * t489;
t404 = t489 * t238;
t400 = qJD(2) * t441;
t398 = t330 * t426;
t395 = t489 * qJD(2);
t69 = -t115 * t335 + t339 * t141;
t386 = -t230 * t339 + t215;
t184 = t339 * t224 - t244 * t335;
t243 = -t488 * t310 + t311 * t336;
t381 = t317 + t428;
t379 = t338 * t340 * t443;
t325 = -pkin(3) * t488 - pkin(4);
t372 = -t100 + t526;
t369 = pkin(1) * t510;
t119 = t174 * t488 - t336 * t187;
t368 = t335 * t53 + t339 * t52;
t366 = -t128 * t324 + t203 * t96;
t363 = t330 * t338 * t393;
t361 = t408 * t489;
t114 = pkin(4) * t441 - t119;
t199 = t218 * t335 + t339 * t441;
t156 = -t210 * qJD(3) + t489 * t273 - t337 * t275;
t233 = -qJD(3) * t279 + t395 * t441;
t131 = pkin(3) * t401 - t233 * pkin(10) + t156;
t155 = -t268 * t424 + t269 * t397 + t337 * t273 + t489 * t275;
t232 = qJD(3) * t280 + t337 * t400;
t138 = -pkin(10) * t232 + t155;
t47 = t336 * t131 + t488 * t138 + t174 * t396 - t187 * t423;
t44 = pkin(11) * t401 + t47;
t148 = t336 * t232 - t233 * t488 + t279 * t396 + t280 * t423;
t149 = qJD(4) * t218 + t232 * t488 + t336 * t233;
t208 = pkin(3) * t232 + t276;
t77 = pkin(4) * t149 + pkin(11) * t148 + t208;
t11 = -t115 * t421 + t141 * t420 + t335 * t77 + t339 * t44;
t48 = t131 * t488 - t336 * t138 - t174 * t423 - t187 * t396;
t353 = -t289 * t421 - t386;
t350 = t360 * t333;
t349 = pkin(1) * (-t334 * t418 + t443);
t347 = t365 - t397;
t12 = -t70 * qJD(5) - t335 * t44 + t339 * t77;
t45 = -pkin(4) * t401 - t48;
t345 = -qJD(5) * t368 - t8 * t335 - t6;
t326 = -pkin(4) - t485;
t309 = pkin(11) * t339 + t329;
t308 = t484 * t335;
t306 = t325 - t485;
t287 = t324 * t339 + t329;
t286 = t439 * t335;
t213 = pkin(5) * t446 + t243;
t200 = t218 * t339 - t410;
t176 = t177 ^ 2;
t162 = -qJ(6) * t446 + t185;
t151 = -pkin(5) * t352 - t289 * t329 + t184;
t103 = -qJD(5) * t410 - t148 * t335 + t218 * t420 - t339 * t401;
t102 = qJD(5) * t199 + t339 * t148 - t335 * t401;
t92 = t199 * pkin(5) + t114;
t89 = -t176 + t490;
t62 = t432 * t516 - t463;
t58 = -qJ(6) * t199 + t70;
t57 = t149 * t516 + t464;
t56 = t459 - t88;
t55 = -t87 + t462;
t50 = pkin(5) * t217 - qJ(6) * t200 + t69;
t38 = t461 * t516 - t86;
t37 = t177 * t530 - t86;
t25 = t103 * t177 + t199 * t88;
t24 = -t102 * t179 - t200 * t87;
t23 = t177 * t354 + t446 * t88;
t22 = t179 * t353 - t445 * t87;
t21 = t103 * pkin(5) + t45;
t18 = -t122 * t289 - t177 * t432 + t352 * t88 - t354 * t516;
t17 = t124 * t289 + t179 * t432 + t352 * t87 + t353 * t516;
t14 = -t103 * t516 - t128 * t199 - t149 * t177 - t217 * t88;
t13 = -t102 * t516 + t128 * t200 + t149 * t179 - t217 * t87;
t10 = t102 * t177 - t103 * t179 + t199 * t87 - t200 * t88;
t9 = (t214 - t449) * t179 + t386 * t177 + (t84 - t86 + (-t179 * t339 + t461) * qJD(5)) * t289;
t5 = -qJ(6) * t103 - qJD(6) * t199 + t11;
t4 = pkin(5) * t149 + qJ(6) * t102 - qJD(6) * t200 + t12;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t363, t430 * t510, t381 * t400, -0.2e1 * t363, -t381 * t401, 0, -t265 * t334 - t276 * t317 + t338 * t369, -t264 * t334 - t275 * t317 + t340 * t369 (t264 * t340 + t265 * t338 + (-t271 * t340 - t274 * t338) * qJD(2) + (t275 * t340 + t276 * t338 + (-t283 * t340 - t338 * t431) * qJD(2)) * qJD(1)) * t333, t264 * t431 - t265 * t283 - t271 * t276 + t274 * t275, -t223 * t280 + t233 * t256, t223 * t279 - t256 * t232 + t233 * t255 - t280 * t408, -t233 * t307 + (t223 * t340 + (qJD(1) * t280 + t256) * t425) * t333, -t255 * t232 + t279 * t408, t232 * t307 + (t408 * t340 + (-qJD(1) * t279 + t255) * t425) * t333 (-t307 * t333 - t398) * t425, -t156 * t307 - t276 * t255 + t267 * t408 + t265 * t279 + t238 * t232 + (-t145 * t340 + (qJD(1) * t209 + t193) * t425) * t333, t155 * t307 - t223 * t267 + t233 * t238 + t256 * t276 + t265 * t280 + (t144 * t340 + (-qJD(1) * t210 - t194) * t425) * t333, -t144 * t279 - t145 * t280 + t155 * t255 - t156 * t256 - t193 * t233 - t194 * t232 + t209 * t223 - t210 * t408, t144 * t210 + t145 * t209 + t155 * t194 + t156 * t193 + t238 * t276 + t265 * t267, -t148 * t357 - t218 * t343, -t218 * t128 + t148 * t203 - t149 * t357 + t217 * t343, -t148 * t417 + (t343 * t340 + t357 * t425 + (t148 * t340 + t218 * t425) * qJD(1)) * t333, t149 * t203 + t464, -t149 * t417 + (-t203 * t425 + t128 * t340 + (t149 * t340 - t217 * t425) * qJD(1)) * t333 (-t350 - t398) * t425, t48 * t417 + t208 * t203 + t222 * t128 + t192 * t217 + t207 * t149 + (t99 * t425 - t31 * t340 + (t119 * t425 - t340 * t48) * qJD(1)) * t333, -t47 * t417 + t208 * t357 - t222 * t343 + t192 * t218 - t207 * t148 + (-t100 * t425 - t378 * t340 + (-t120 * t425 + t340 * t47) * qJD(1)) * t333, -t100 * t149 + t119 * t343 - t120 * t128 + t99 * t148 - t47 * t203 + t217 * t378 - t31 * t218 - t357 * t48, t100 * t47 + t119 * t31 - t120 * t378 + t192 * t222 + t207 * t208 + t48 * t99, t24, t10, t13, t25, t14, t57, t103 * t96 + t114 * t88 + t12 * t516 + t128 * t69 + t149 * t52 + t177 * t45 + t199 * t29 + t217 * t8, -t102 * t96 - t11 * t516 - t114 * t87 - t128 * t70 - t149 * t53 + t179 * t45 + t200 * t29 + t217 * t388, t102 * t52 - t103 * t53 - t11 * t177 - t12 * t179 + t199 * t388 - t200 * t8 + t69 * t87 - t70 * t88, t11 * t53 + t114 * t29 + t12 * t52 - t388 * t70 + t45 * t96 + t69 * t8, t24, t10, t13, t25, t14, t57, t1 * t217 + t103 * t78 + t128 * t50 + t149 * t32 + t177 * t21 + t199 * t20 + t4 * t516 + t88 * t92, -t102 * t78 - t128 * t58 - t149 * t34 + t179 * t21 + t20 * t200 - t217 * t3 - t5 * t516 - t87 * t92, -t1 * t200 + t102 * t32 - t103 * t34 - t177 * t5 - t179 * t4 - t199 * t3 + t50 * t87 - t58 * t88, t1 * t50 + t20 * t92 + t21 * t78 + t3 * t58 + t32 * t4 + t34 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t379, t430 * t443 (qJD(2) - t317) * t402, t379, t317 * t403 - t312, 0, -pkin(8) * t373 + t274 * t317 + t338 * t349, pkin(8) * t312 + t271 * t317 + t340 * t349, 0, 0, -t223 * t337 - t256 * t347, -t223 * t489 - t347 * t255 + (-t408 + t447) * t337, -t307 * t397 + (t307 * t405 + (qJD(2) * t337 - t256) * t338) * t429, t255 * t444 - t361, t307 * t424 + (-t340 * t444 + (t395 - t255) * t338) * t429, t307 * t403, -pkin(2) * t408 - t265 * t489 + t211 * t307 + t274 * t255 + (t238 * t337 + t307 * t414) * qJD(3) + (-t193 * t338 + (-pkin(9) * t425 - t238 * t340) * t337) * t429, pkin(2) * t223 - t212 * t307 - t274 * t256 + t265 * t337 + (-pkin(9) * t444 + t404) * qJD(3) + (-t340 * t404 + (-pkin(9) * t395 + t194) * t338) * t429, t407 + t211 * t256 - t212 * t255 + t347 * t193 + (t256 * t397 - t361) * pkin(9) + ((-qJD(3) * t255 - t223) * pkin(9) + t501) * t337, -t265 * pkin(2) - t193 * t211 - t194 * t212 - t238 * t274 + (t407 - t145 * t337 + (-t193 * t489 - t194 * t337) * qJD(3)) * pkin(9), -t289 * t343 + t357 * t537, -t289 * t128 - t203 * t537 - t343 * t352 - t357 * t432 (-t537 * t340 + (qJD(2) * t289 - t357) * t338) * t429 + t537 * t417, t203 * t432 - t463 (t432 * t340 + (qJD(2) * t352 + t203) * t338) * t429 - t432 * t417, t350 * t427, t327 * t128 - t192 * t352 + t432 * t207 + t371 * t203 + (-t434 * t340 + (-qJD(2) * t243 - t99) * t338) * t429 + t434 * t417, -t327 * t343 + t192 * t289 + t537 * t207 + t371 * t357 + (-t435 * t340 + (-qJD(2) * t244 + t100) * t338) * t429 + t435 * t417, -t100 * t432 - t244 * t128 + t203 * t435 - t243 * t343 - t31 * t289 - t352 * t378 - t357 * t434 - t537 * t99, -t100 * t435 + t192 * t327 + t207 * t371 - t243 * t31 - t244 * t378 + t434 * t99, t22, t9, t17, t23, t18, t62, t128 * t184 + t177 * t436 + t243 * t88 + t29 * t446 - t352 * t8 + t354 * t96 + t432 * t52 - t478 * t516, -t128 * t185 + t179 * t436 - t243 * t87 + t29 * t445 - t352 * t388 + t353 * t96 - t432 * t53 + t477 * t516, t184 * t87 - t185 * t88 + t214 * t53 + t215 * t52 - t368 * t230 + t478 * t179 + t477 * t177 + (-qJD(5) * t499 + t335 * t388 - t339 * t8) * t289, t184 * t8 - t185 * t388 + t243 * t29 + t436 * t96 - t477 * t53 - t478 * t52, t22, t9, t17, t23, t18, t62, -t1 * t352 + t128 * t151 + t177 * t438 + t20 * t446 + t213 * t88 + t32 * t432 + t354 * t78 + t480 * t516, -t128 * t162 + t179 * t438 + t20 * t445 - t213 * t87 + t3 * t352 - t34 * t432 + t353 * t78 - t479 * t516, t151 * t87 - t162 * t88 + t214 * t34 + t215 * t32 - (t335 * t34 + t473) * t230 - t480 * t179 - t479 * t177 + (-t1 * t339 - t3 * t335 + (t32 * t335 - t339 * t34) * qJD(5)) * t289, t1 * t151 + t162 * t3 + t20 * t213 + t32 * t480 + t34 * t479 + t438 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t448, -t255 ^ 2 + t256 ^ 2, t255 * t307 - t223, t448, -t408 - t447, t312, -t238 * t256 - t501, -t193 * t307 - t238 * t255 - t144, 0, 0, t450, t515, t512, -t450, t491, t312, -t105 * t360 + (-t256 * t203 + t312 * t488 + t360 * t423) * pkin(3) + t498, -t106 * t360 + (-t256 * t357 - t312 * t336 + t360 * t396) * pkin(3) + t514, t100 * t357 - t105 * t357 + t106 * t203 - t99 * t203 + (-t336 * t128 + t488 * t346 + (t336 * t357 + (-t203 - t502) * t488) * qJD(4)) * pkin(3), -t100 * t106 + t99 * t105 + (t488 * t31 - t207 * t256 - t378 * t336 + (t100 * t488 - t336 * t99) * qJD(4)) * pkin(3), t520, t534, t519, t38, t532, -t503, t325 * t88 + t366 * t335 + t370 * t177 + (-t324 * t420 - t335 * t380 - t64) * t516 + t495, -t325 * t87 + t366 * t339 + t370 * t179 + (t324 * t421 + t505) * t516 + t494, t65 * t177 + t64 * t179 - t6 + (-t177 * t380 - t203 * t52 - t324 * t88 + (t179 * t324 - t52) * qJD(5)) * t339 + (t179 * t380 - t203 * t53 - t324 * t87 - t8 + (t177 * t324 - t53) * qJD(5)) * t335, -t96 * t105 + t29 * t325 - t52 * t64 - t53 * t65 + (t336 * t96 + t488 * t499) * qJD(4) * pkin(3) + t345 * t324, t520, t534, t519, t38, t532, -t503, t128 * t286 + t177 * t465 + t306 * t88 + t468 * t516 + t523, -t128 * t287 + t179 * t465 - t306 * t87 - t469 * t516 + t524, -t469 * t177 - t468 * t179 + t286 * t87 - t287 * t88 + t533, t1 * t286 + t20 * t306 + t287 * t3 + t32 * t468 + t34 * t469 + t465 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t450, t515, t512, -t450, t491, t312, -t100 * t360 + t498, -t360 * t99 + t514, 0, 0, t520, t534, t519, t37, t532, -t503, -pkin(4) * t88 + pkin(11) * t500 - t100 * t177 - t516 * t67 + t518 * t96 + t495, t96 * t517 + pkin(4) * t87 - t100 * t179 + t516 * t68 + (t421 * t516 - t124) * pkin(11) + t494, t177 * t68 + t179 * t67 - t6 + (-t475 + (-t88 + t422) * pkin(11)) * t339 + ((qJD(5) * t177 - t87) * pkin(11) + t507) * t335, -pkin(4) * t29 + pkin(11) * t345 - t100 * t96 - t52 * t67 - t53 * t68, t520, t534, t519, t37, t532, -t503, t128 * t308 + t177 * t372 + t326 * t88 + t466 * t516 + t523, -t128 * t309 + t179 * t372 - t326 * t87 - t467 * t516 + t524, -t467 * t177 - t466 * t179 + t308 * t87 - t309 * t88 + t533, t1 * t308 + t20 * t326 + t3 * t309 + t32 * t466 + t34 * t467 + t372 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t460, t89, t55, -t460, t56, t128, -t179 * t96 - t507, t177 * t96 + t388 + t475, 0, 0, t460, t89, t55, -t460, t56, t128, 0.2e1 * t486 + t476 + (-t389 - t78) * t179 + t344, -pkin(5) * t490 + t516 * t33 + (qJD(6) + t78) * t177 + t359, pkin(5) * t87 - t177 * t481, t481 * t34 + (-t179 * t78 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88 + t459, -t87 - t462, -t176 - t490, t177 * t34 + t179 * t32 + t20;];
tauc_reg  = t2;
