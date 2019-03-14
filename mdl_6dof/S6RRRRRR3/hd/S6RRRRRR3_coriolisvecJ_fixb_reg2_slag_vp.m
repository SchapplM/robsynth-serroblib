% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRRRR3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:42:39
% EndTime: 2019-03-10 03:43:10
% DurationCPUTime: 14.04s
% Computational Cost: add. (34749->689), mult. (80747->904), div. (0->0), fcn. (60502->10), ass. (0->336)
t376 = cos(qJ(2));
t507 = cos(qJ(3));
t424 = t507 * t376;
t352 = qJD(1) * t424;
t372 = sin(qJ(3));
t373 = sin(qJ(2));
t442 = qJD(1) * t373;
t422 = t372 * t442;
t311 = -t352 + t422;
t330 = t372 * t376 + t373 * t507;
t313 = qJD(1) * t330;
t253 = pkin(3) * t313 + pkin(9) * t311;
t431 = pkin(2) * t442;
t231 = t253 + t431;
t509 = -pkin(8) - pkin(7);
t347 = t509 * t376;
t337 = qJD(1) * t347;
t315 = t372 * t337;
t345 = t509 * t373;
t335 = qJD(1) * t345;
t266 = t335 * t507 + t315;
t371 = sin(qJ(4));
t375 = cos(qJ(4));
t175 = t375 * t231 - t266 * t371;
t419 = qJD(3) * t507;
t407 = pkin(2) * t419;
t548 = -t371 * t407 - t175;
t176 = t371 * t231 + t375 * t266;
t547 = -t375 * t407 + t176;
t468 = t311 * t375;
t404 = t313 * pkin(4) + pkin(10) * t468;
t359 = pkin(2) * t372 + pkin(9);
t500 = -pkin(10) - t359;
t414 = qJD(4) * t500;
t546 = t375 * t414 - t404 + t548;
t494 = qJD(2) * pkin(2);
t320 = t335 + t494;
t257 = t320 * t507 + t315;
t179 = t375 * t253 - t257 * t371;
t508 = -pkin(10) - pkin(9);
t425 = qJD(4) * t508;
t545 = t375 * t425 - t179 - t404;
t469 = t311 * t371;
t432 = pkin(10) * t469;
t544 = -t371 * t414 + t432 + t547;
t180 = t371 * t253 + t375 * t257;
t543 = -t371 * t425 + t180 + t432;
t374 = cos(qJ(5));
t433 = qJD(4) + qJD(5);
t437 = qJD(5) * t374;
t439 = qJD(4) * t375;
t370 = sin(qJ(5));
t462 = t370 * t371;
t515 = t374 * t375 - t462;
t447 = -t515 * t311 - t374 * t439 - t375 * t437 + t433 * t462;
t461 = t370 * t375;
t329 = t371 * t374 + t461;
t271 = t433 * t329;
t446 = t329 * t311 + t271;
t434 = qJD(2) + qJD(3);
t408 = t375 * t434;
t279 = t313 * t371 - t408;
t281 = t375 * t313 + t371 * t434;
t200 = t279 * t374 + t281 * t370;
t369 = sin(qJ(6));
t395 = -t279 * t370 + t374 * t281;
t506 = cos(qJ(6));
t118 = t506 * t200 + t369 * t395;
t526 = -t369 * t200 + t395 * t506;
t486 = t118 * t526;
t532 = -t118 ^ 2 + t526 ^ 2;
t537 = pkin(11) * t200;
t363 = -pkin(2) * t376 - pkin(1);
t343 = qJD(1) * t363;
t227 = t311 * pkin(3) - t313 * pkin(9) + t343;
t318 = t507 * t337;
t258 = t372 * t320 - t318;
t235 = pkin(9) * t434 + t258;
t163 = t375 * t227 - t235 * t371;
t131 = -pkin(10) * t281 + t163;
t305 = qJD(4) + t311;
t114 = pkin(4) * t305 + t131;
t164 = t227 * t371 + t235 * t375;
t132 = -pkin(10) * t279 + t164;
t128 = t374 * t132;
t73 = t114 * t370 + t128;
t57 = t73 - t537;
t493 = t369 * t57;
t296 = qJD(5) + t305;
t521 = pkin(11) * t395;
t126 = t370 * t132;
t72 = t374 * t114 - t126;
t56 = t72 - t521;
t52 = pkin(5) * t296 + t56;
t24 = t506 * t52 - t493;
t428 = t506 * t57;
t25 = t369 * t52 + t428;
t542 = -t24 * t118 + t25 * t526;
t293 = qJD(6) + t296;
t457 = t372 * t373;
t400 = t434 * t457;
t444 = t434 * t352;
t249 = qJD(1) * t400 - t444;
t440 = qJD(4) * t371;
t189 = -qJD(4) * t408 + t375 * t249 + t313 * t440;
t190 = qJD(4) * t281 - t371 * t249;
t438 = qJD(5) * t370;
t394 = -t370 * t189 + t190 * t374 - t279 * t438 + t281 * t437;
t418 = qJD(6) * t506;
t436 = qJD(6) * t369;
t95 = t374 * t189 + t370 * t190 + t279 * t437 + t281 * t438;
t32 = t200 * t418 + t369 * t394 + t395 * t436 + t506 * t95;
t530 = t118 * t293 - t32;
t234 = -pkin(3) * t434 - t257;
t191 = t279 * pkin(4) + t234;
t115 = t200 * pkin(5) + t191;
t273 = t434 * t330;
t250 = t273 * qJD(1);
t435 = qJD(1) * qJD(2);
t417 = t373 * t435;
t173 = pkin(2) * t417 + pkin(3) * t250 + pkin(9) * t249;
t426 = qJD(2) * t509;
t406 = qJD(1) * t426;
t321 = t373 * t406;
t322 = t376 * t406;
t441 = qJD(3) * t372;
t181 = t320 * t419 + t321 * t507 + t372 * t322 + t337 * t441;
t79 = -qJD(4) * t164 + t375 * t173 - t371 * t181;
t50 = t250 * pkin(4) + t189 * pkin(10) + t79;
t78 = t371 * t173 + t375 * t181 + t227 * t439 - t235 * t440;
t55 = -pkin(10) * t190 + t78;
t409 = -t114 * t437 + t132 * t438 - t370 * t50 - t374 * t55;
t10 = -pkin(11) * t394 - t409;
t16 = -qJD(5) * t73 - t370 * t55 + t374 * t50;
t9 = t250 * pkin(5) + t95 * pkin(11) + t16;
t382 = -t10 * t506 - t369 * t9 - t418 * t52 + t57 * t436;
t528 = t115 * t118 + t382;
t324 = t500 * t371;
t366 = t375 * pkin(10);
t325 = t359 * t375 + t366;
t492 = t324 * t437 - t325 * t438 + t370 * t546 - t544 * t374;
t252 = t370 * t324 + t374 * t325;
t491 = -qJD(5) * t252 + t544 * t370 + t374 * t546;
t344 = t508 * t371;
t346 = pkin(9) * t375 + t366;
t490 = t344 * t437 - t346 * t438 + t370 * t545 - t374 * t543;
t284 = t370 * t344 + t374 * t346;
t489 = -qJD(5) * t284 + t370 * t543 + t374 * t545;
t539 = t446 * pkin(11);
t538 = t313 * pkin(5) - pkin(11) * t447;
t4 = -qJD(6) * t25 - t369 * t10 + t506 * t9;
t513 = -t115 * t526 + t4;
t33 = qJD(6) * t526 - t369 * t95 + t506 * t394;
t511 = t293 * t526 - t33;
t536 = -t539 + t492;
t535 = t538 - t491;
t534 = -t539 + t490;
t533 = t538 - t489;
t480 = t200 * t395;
t451 = t329 * t436 + t369 * t446 - t418 * t515 + t447 * t506;
t263 = t329 * t506 + t369 * t515;
t450 = qJD(6) * t263 - t369 * t447 + t446 * t506;
t265 = t335 * t372 - t318;
t403 = pkin(2) * t441 - t265;
t272 = -qJD(2) * t424 - t376 * t419 + t400;
t479 = t272 * t371;
t391 = t330 * t439 - t479;
t531 = -t200 ^ 2 + t395 ^ 2;
t529 = t200 * t296 - t95;
t527 = t191 * t200 + t409;
t430 = t373 * t494;
t523 = 0.2e1 * t430;
t522 = -0.2e1 * t435;
t328 = -t424 + t457;
t256 = pkin(3) * t328 - pkin(9) * t330 + t363;
t285 = t372 * t345 - t347 * t507;
t192 = t375 * t256 - t285 * t371;
t465 = t330 * t375;
t151 = pkin(4) * t328 - pkin(10) * t465 + t192;
t275 = t375 * t285;
t193 = t371 * t256 + t275;
t466 = t330 * t371;
t170 = -pkin(10) * t466 + t193;
t99 = t370 * t151 + t374 * t170;
t364 = pkin(4) * t440;
t518 = pkin(5) * t446 + t364;
t517 = t507 * t345 + t372 * t347;
t294 = pkin(4) * t469;
t516 = t294 + t403;
t514 = -t163 * t371 + t164 * t375;
t512 = -t191 * t395 + t16;
t510 = t296 * t395 - t394;
t505 = pkin(11) * t329;
t502 = t515 * pkin(5);
t501 = t375 * pkin(4);
t251 = t374 * t324 - t325 * t370;
t212 = t251 - t505;
t323 = t515 * pkin(11);
t213 = t323 + t252;
t141 = t369 * t212 + t213 * t506;
t499 = qJD(6) * t141 + t369 * t536 + t506 * t535;
t140 = t212 * t506 - t369 * t213;
t498 = -qJD(6) * t140 + t369 * t535 - t506 * t536;
t282 = t374 * t344 - t346 * t370;
t232 = t282 - t505;
t233 = t323 + t284;
t167 = t369 * t232 + t233 * t506;
t497 = qJD(6) * t167 + t369 * t534 + t506 * t533;
t166 = t232 * t506 - t369 * t233;
t496 = -qJD(6) * t166 + t369 * t533 - t506 * t534;
t495 = pkin(2) * qJD(3);
t76 = t78 * t375;
t360 = pkin(4) * t374 + pkin(5);
t463 = t369 * t370;
t81 = -t131 * t370 - t128;
t58 = t81 + t537;
t82 = t374 * t131 - t126;
t59 = t82 - t521;
t488 = -t369 * t58 - t506 * t59 + t360 * t418 + (-t370 * t436 + (t374 * t506 - t463) * qJD(5)) * pkin(4);
t423 = t506 * t370;
t487 = t369 * t59 - t506 * t58 - t360 * t436 + (-t370 * t418 + (-t369 * t374 - t423) * qJD(5)) * pkin(4);
t298 = t507 * t322;
t411 = t372 * t321 - t298;
t182 = qJD(3) * t258 + t411;
t483 = t182 * t517;
t482 = t189 * t371;
t481 = t190 * t375;
t214 = t250 * t328;
t478 = t272 * t375;
t477 = t279 * t305;
t476 = t279 * t371;
t475 = t281 * t279;
t474 = t281 * t305;
t473 = t281 * t375;
t472 = t293 * t313;
t471 = t296 * t313;
t470 = t305 * t313;
t467 = t313 * t311;
t464 = t343 * t313;
t459 = t371 * t250;
t455 = t375 * t250;
t378 = qJD(1) ^ 2;
t454 = t376 * t378;
t377 = qJD(2) ^ 2;
t453 = t377 * t373;
t452 = t377 * t376;
t210 = -t294 + t258;
t449 = -t210 + t518;
t448 = t516 + t518;
t445 = t364 + t516;
t443 = t373 ^ 2 - t376 ^ 2;
t429 = -t163 * t468 - t164 * t469 + t76;
t427 = t373 * t454;
t362 = -pkin(3) - t501;
t421 = t330 * t440;
t413 = pkin(1) * t522;
t98 = t374 * t151 - t170 * t370;
t188 = pkin(3) * t273 + pkin(9) * t272 + t430;
t336 = t373 * t426;
t338 = t376 * t426;
t206 = qJD(3) * t517 + t507 * t336 + t372 * t338;
t412 = t375 * t188 - t371 * t206;
t410 = t305 * t375;
t361 = -pkin(2) * t507 - pkin(3);
t405 = t376 * t417;
t402 = -t210 + t364;
t401 = t164 * t313 + t182 * t371 + t234 * t439;
t229 = pkin(4) * t466 - t517;
t397 = t163 * t375 + t164 * t371;
t396 = t234 * t311 - t250 * t359;
t242 = t515 * t330;
t80 = pkin(5) * t328 - pkin(11) * t242 + t98;
t241 = t329 * t330;
t87 = -pkin(11) * t241 + t99;
t36 = -t369 * t87 + t506 * t80;
t37 = t369 * t80 + t506 * t87;
t393 = -t163 * t313 - t182 * t375 + t234 * t440;
t172 = -t369 * t241 + t242 * t506;
t390 = -t421 - t478;
t342 = t361 - t501;
t65 = pkin(10) * t478 + t273 * pkin(4) + (-t275 + (pkin(10) * t330 - t256) * t371) * qJD(4) + t412;
t96 = t371 * t188 + t375 * t206 + t256 * t439 - t285 * t440;
t77 = -pkin(10) * t391 + t96;
t22 = t151 * t437 - t170 * t438 + t370 * t65 + t374 * t77;
t262 = t329 * t369 - t506 * t515;
t109 = t190 * pkin(4) + t182;
t51 = pkin(5) * t394 + t109;
t389 = t115 * t450 - t24 * t313 + t51 * t262;
t388 = -t115 * t451 + t25 * t313 + t51 * t263;
t387 = -t109 * t515 + t191 * t446 - t72 * t313;
t386 = t109 * t329 - t191 * t447 + t73 * t313;
t385 = t24 * t451 - t25 * t450 + t382 * t262 - t4 * t263;
t384 = -t16 * t329 - t409 * t515 - t446 * t73 + t447 * t72;
t383 = -qJD(4) * t397 - t79 * t371;
t23 = -qJD(5) * t99 - t370 * t77 + t374 * t65;
t207 = qJD(3) * t285 + t372 * t336 - t507 * t338;
t380 = t383 + t76;
t137 = pkin(4) * t391 + t207;
t379 = t343 * t311 - t181;
t308 = pkin(4) * t423 + t369 * t360;
t307 = -pkin(4) * t463 + t360 * t506;
t290 = t362 - t502;
t288 = t342 - t502;
t215 = -t311 ^ 2 + t313 ^ 2;
t208 = t444 + (t311 - t422) * t434;
t174 = pkin(5) * t241 + t229;
t171 = t241 * t506 + t242 * t369;
t162 = pkin(4) * t281 + pkin(5) * t395;
t113 = -t281 * t313 + t305 * t410 + t459;
t112 = -t305 ^ 2 * t371 + t279 * t313 + t455;
t108 = t305 * t476 - t481;
t107 = t281 * t410 - t482;
t106 = -t272 * t461 - t370 * t421 - t438 * t466 + (t433 * t465 - t479) * t374;
t105 = t271 * t330 + t272 * t515;
t97 = -qJD(4) * t193 + t412;
t84 = t200 * t313 + t250 * t515 - t296 * t446;
t83 = t329 * t250 - t296 * t447 - t313 * t395;
t71 = t106 * pkin(5) + t137;
t60 = (-t189 - t477) * t375 + (-t190 - t474) * t371;
t43 = qJD(6) * t172 - t369 * t105 + t106 * t506;
t42 = t105 * t506 + t369 * t106 + t241 * t418 + t242 * t436;
t41 = t200 * t446 - t394 * t515;
t40 = -t95 * t329 - t395 * t447;
t39 = t118 * t313 - t262 * t250 - t293 * t450;
t38 = t263 * t250 - t293 * t451 - t313 * t526;
t27 = t506 * t56 - t493;
t26 = -t369 * t56 - t428;
t19 = t200 * t447 - t329 * t394 - t395 * t446 - t515 * t95;
t18 = -pkin(11) * t106 + t22;
t17 = t273 * pkin(5) + t105 * pkin(11) + t23;
t12 = t118 * t450 + t33 * t262;
t11 = -t32 * t263 - t451 * t526;
t7 = t118 * t451 + t32 * t262 - t263 * t33 - t450 * t526;
t6 = -qJD(6) * t37 + t17 * t506 - t369 * t18;
t5 = qJD(6) * t36 + t369 * t17 + t18 * t506;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t405, t443 * t522, t452, -0.2e1 * t405, -t453, 0, -pkin(7) * t452 + t373 * t413, pkin(7) * t453 + t376 * t413, 0, 0, -t249 * t330 - t272 * t313, t249 * t328 - t250 * t330 + t272 * t311 - t273 * t313, -t272 * t434, t273 * t311 + t214, -t273 * t434, 0, t363 * t250 + t343 * t273 - t207 * t434 + (qJD(1) * t328 + t311) * t430, -t206 * t434 - t363 * t249 - t343 * t272 + t313 * t523, -t181 * t328 + t182 * t330 - t206 * t311 + t207 * t313 + t249 * t517 - t250 * t285 + t257 * t272 - t258 * t273, t181 * t285 + t258 * t206 - t257 * t207 + t343 * t523 - t483, -t189 * t465 + t281 * t390 (t279 * t375 + t281 * t371) * t272 + (t482 - t481 + (-t473 + t476) * qJD(4)) * t330, -t189 * t328 + t281 * t273 + t305 * t390 + t330 * t455, t190 * t466 + t279 * t391, -t190 * t328 - t279 * t273 - t305 * t391 - t330 * t459, t273 * t305 + t214, t163 * t273 + t182 * t466 - t190 * t517 + t192 * t250 + t207 * t279 + t234 * t391 + t97 * t305 + t79 * t328, -t164 * t273 + t182 * t465 + t189 * t517 - t193 * t250 + t207 * t281 + t234 * t390 - t96 * t305 - t78 * t328, t192 * t189 - t193 * t190 - t96 * t279 - t97 * t281 + t397 * t272 + (-qJD(4) * t514 - t371 * t78 - t375 * t79) * t330, t163 * t97 + t164 * t96 + t192 * t79 + t193 * t78 + t207 * t234 - t483, -t105 * t395 - t242 * t95, t105 * t200 - t106 * t395 + t95 * t241 - t242 * t394, -t105 * t296 + t242 * t250 + t273 * t395 - t328 * t95, t200 * t106 + t241 * t394, -t106 * t296 - t200 * t273 - t241 * t250 - t328 * t394, t273 * t296 + t214, t191 * t106 + t109 * t241 + t137 * t200 + t16 * t328 + t229 * t394 + t23 * t296 + t98 * t250 + t72 * t273, -t105 * t191 + t109 * t242 + t137 * t395 - t22 * t296 - t229 * t95 - t250 * t99 - t273 * t73 + t328 * t409, t72 * t105 - t73 * t106 - t16 * t242 - t22 * t200 - t23 * t395 + t241 * t409 - t394 * t99 + t98 * t95, t109 * t229 + t137 * t191 + t16 * t98 + t22 * t73 + t23 * t72 - t409 * t99, -t172 * t32 - t42 * t526, t118 * t42 + t171 * t32 - t172 * t33 - t43 * t526, t172 * t250 + t273 * t526 - t293 * t42 - t32 * t328, t118 * t43 + t171 * t33, -t118 * t273 - t171 * t250 - t293 * t43 - t328 * t33, t273 * t293 + t214, t115 * t43 + t118 * t71 + t171 * t51 + t174 * t33 + t24 * t273 + t250 * t36 + t293 * t6 + t328 * t4, -t115 * t42 + t172 * t51 - t174 * t32 - t25 * t273 - t250 * t37 - t293 * t5 + t328 * t382 + t526 * t71, -t118 * t5 + t171 * t382 - t172 * t4 + t24 * t42 - t25 * t43 + t32 * t36 - t33 * t37 - t526 * t6, t115 * t71 + t174 * t51 + t24 * t6 + t25 * t5 + t36 * t4 - t37 * t382; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t427, t443 * t378, 0, t427, 0, 0, t378 * pkin(1) * t373, pkin(1) * t454, 0, 0, t467, t215, t208, -t467, 0, 0, t337 * t419 + t298 - t311 * t431 - t464 + t265 * t434 + (-qJD(3) * t320 - t434 * t495 - t321) * t372, t266 * t434 + (-t313 * t442 - t419 * t434) * pkin(2) + t379 (t258 - t265) * t313 + (-t257 + t266) * t311 + (t507 * t249 - t250 * t372 + (-t311 * t507 + t313 * t372) * qJD(3)) * pkin(2), t257 * t265 - t258 * t266 + (-t343 * t442 - t507 * t182 + t181 * t372 + (-t257 * t372 + t258 * t507) * qJD(3)) * pkin(2), t107, t60, t113, t108, t112, -t470, t361 * t190 + t396 * t371 + t403 * t279 + (-t359 * t439 + t548) * t305 + t393, -t361 * t189 + t396 * t375 + t403 * t281 + (t359 * t440 + t547) * t305 + t401, t175 * t281 + t176 * t279 + (-t279 * t407 - t190 * t359 + (t281 * t359 - t163) * qJD(4)) * t375 + (t281 * t407 - t189 * t359 - t79 + (t279 * t359 - t164) * qJD(4)) * t371 + t429, -t163 * t175 - t164 * t176 + t182 * t361 - t234 * t265 + (t234 * t372 + t507 * t514) * t495 + t380 * t359, t40, t19, t83, t41, t84, -t471, t200 * t445 + t251 * t250 + t296 * t491 + t342 * t394 + t387, -t252 * t250 - t296 * t492 - t342 * t95 + t395 * t445 + t386, -t200 * t492 + t251 * t95 - t252 * t394 - t395 * t491 + t384, t109 * t342 + t16 * t251 + t191 * t445 - t252 * t409 + t491 * t72 + t492 * t73, t11, t7, t38, t12, t39, -t472, t118 * t448 + t140 * t250 + t288 * t33 - t293 * t499 + t389, -t141 * t250 - t288 * t32 + t293 * t498 + t448 * t526 + t388, t118 * t498 + t140 * t32 - t141 * t33 + t499 * t526 + t385, t115 * t448 + t4 * t140 - t141 * t382 - t24 * t499 - t25 * t498 + t51 * t288; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t467, t215, t208, -t467, 0, 0, t258 * qJD(2) - t411 - t464, t257 * t434 + t379, 0, 0, t107, t60, t113, t108, t112, -t470, t234 * t469 - pkin(3) * t190 - t179 * t305 - t258 * t279 + (-t305 * t439 - t459) * pkin(9) + t393, t234 * t468 + pkin(3) * t189 + t180 * t305 - t258 * t281 + (t305 * t440 - t455) * pkin(9) + t401, t179 * t281 + t180 * t279 + (-t482 - t481 + (t473 + t476) * qJD(4)) * pkin(9) + t383 + t429, -t182 * pkin(3) + pkin(9) * t380 - t163 * t179 - t164 * t180 - t234 * t258, t40, t19, t83, t41, t84, -t471, t200 * t402 + t282 * t250 + t296 * t489 + t362 * t394 + t387, -t284 * t250 - t296 * t490 - t362 * t95 + t395 * t402 + t386, -t200 * t490 + t282 * t95 - t284 * t394 - t395 * t489 + t384, t109 * t362 + t16 * t282 + t191 * t402 - t284 * t409 + t489 * t72 + t490 * t73, t11, t7, t38, t12, t39, -t472, t118 * t449 + t166 * t250 + t290 * t33 - t293 * t497 + t389, -t167 * t250 - t290 * t32 + t293 * t496 + t449 * t526 + t388, t118 * t496 + t166 * t32 - t167 * t33 + t497 * t526 + t385, t115 * t449 + t4 * t166 - t167 * t382 - t24 * t497 - t25 * t496 + t51 * t290; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t475, -t279 ^ 2 + t281 ^ 2, -t189 + t477, -t475, -t190 + t474, t250, t164 * t305 - t234 * t281 + t79, t163 * t305 + t234 * t279 - t78, 0, 0, t480, t531, t529, -t480, t510, t250, -t81 * t296 + (-t200 * t281 + t250 * t374 - t296 * t438) * pkin(4) + t512, t82 * t296 + (-t250 * t370 - t281 * t395 - t296 * t437) * pkin(4) + t527, t73 * t395 + t82 * t200 - t72 * t200 + t81 * t395 + (-t370 * t394 + t374 * t95 + (-t200 * t374 + t370 * t395) * qJD(5)) * pkin(4), -t72 * t81 - t73 * t82 + (-t409 * t370 + t16 * t374 - t191 * t281 + (-t370 * t72 + t374 * t73) * qJD(5)) * pkin(4), t486, t532, t530, -t486, t511, t250, -t162 * t118 + t307 * t250 + t293 * t487 + t513, -t162 * t526 - t308 * t250 - t293 * t488 + t528, -t118 * t488 + t307 * t32 - t308 * t33 - t487 * t526 + t542, -t115 * t162 + t24 * t487 + t25 * t488 + t4 * t307 - t308 * t382; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t480, t531, t529, -t480, t510, t250, t73 * t296 + t512, t296 * t72 + t527, 0, 0, t486, t532, t530, -t486, t511, t250, -t26 * t293 + (-t118 * t395 + t250 * t506 - t293 * t436) * pkin(5) + t513, t27 * t293 + (-t250 * t369 - t293 * t418 - t395 * t526) * pkin(5) + t528, t27 * t118 + t26 * t526 + (t506 * t32 - t33 * t369 + (-t118 * t506 + t369 * t526) * qJD(6)) * pkin(5) + t542, -t24 * t26 - t25 * t27 + (t506 * t4 - t115 * t395 - t382 * t369 + (-t24 * t369 + t25 * t506) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t486, t532, t530, -t486, t511, t250, t25 * t293 + t513, t24 * t293 + t528, 0, 0;];
tauc_reg  = t1;