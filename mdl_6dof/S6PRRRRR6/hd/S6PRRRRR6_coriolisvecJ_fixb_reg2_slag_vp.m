% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6PRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRRRR6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:25:05
% EndTime: 2019-03-09 01:26:00
% DurationCPUTime: 26.13s
% Computational Cost: add. (34209->795), mult. (101775->1154), div. (0->0), fcn. (86435->16), ass. (0->360)
t335 = cos(qJ(3));
t502 = cos(pkin(7));
t428 = t335 * t502;
t318 = pkin(2) * t428;
t333 = sin(qJ(2));
t336 = cos(qJ(2));
t332 = sin(qJ(3));
t430 = t332 * t502;
t360 = -t333 * t430 + t335 * t336;
t327 = sin(pkin(6));
t474 = qJD(1) * t327;
t560 = -qJD(3) * t318 + t360 * t474;
t326 = sin(pkin(7));
t501 = cos(pkin(8));
t375 = t326 * (-pkin(11) * t501 - pkin(10));
t366 = t332 * t375;
t559 = qJD(3) * t366 - t560;
t317 = pkin(2) * t430;
t231 = (t335 * t375 - t317) * qJD(3);
t358 = -t332 * t336 - t333 * t428;
t246 = t358 * t474;
t558 = -t231 + t246;
t443 = t333 * t474;
t472 = qJD(2) * t326;
t287 = pkin(10) * t472 + t443;
t432 = t326 * t501;
t415 = pkin(11) * t432;
t523 = qJD(2) * pkin(2);
t299 = t336 * t474 + t523;
t425 = t502 * t299;
t328 = cos(pkin(6));
t473 = qJD(1) * t328;
t541 = t326 * t473 + t425;
t476 = t541 * t335;
t175 = (-qJD(2) * t415 - t287) * t332 + t476;
t484 = t335 * t287;
t361 = -t332 * t425 - t484;
t427 = t335 * t501;
t402 = qJD(2) * t427;
t176 = (-pkin(11) * t402 - t332 * t473) * t326 + t361;
t325 = sin(pkin(8));
t527 = pkin(11) * t325;
t380 = pkin(3) * t332 - t335 * t527;
t258 = t380 * t472;
t529 = cos(qJ(4));
t403 = t501 * t529;
t331 = sin(qJ(4));
t492 = t325 * t331;
t279 = pkin(3) * t403 - pkin(11) * t492;
t431 = t331 * t501;
t508 = t279 * qJD(4) - t529 * t175 - t176 * t431 - t258 * t492;
t490 = t326 * t335;
t282 = pkin(10) * t490 + t317;
t424 = t502 * t325;
t218 = (t326 * t427 + t424) * pkin(11) + t282;
t229 = pkin(3) * t502 + t318 + t366;
t379 = -pkin(3) * t335 - t332 * t527;
t254 = (-pkin(2) + t379) * t326;
t364 = t380 * qJD(3);
t259 = t326 * t364;
t378 = qJD(4) * t403;
t437 = qJD(4) * t529;
t406 = t325 * t437;
t489 = t327 * t333;
t454 = t326 * t489;
t416 = t325 * t454;
t469 = qJD(4) * t331;
t504 = -t218 * t469 + t229 * t378 + t231 * t431 + t254 * t406 + t259 * t492 - (qJD(1) * t416 + t246 * t501) * t331 + t559 * t529;
t383 = t432 * t489;
t481 = -qJD(1) * t383 + t501 * t259 + t325 * t558;
t247 = (-t332 * t431 + t335 * t529) * t472;
t557 = t406 - t247;
t470 = qJD(3) * t326;
t440 = t332 * t470;
t410 = t325 * t440;
t556 = -pkin(12) * t410 - t504;
t385 = qJD(3) * t501 + qJD(4);
t491 = t326 * t332;
t453 = t331 * t491;
t351 = t385 * t453;
t407 = t529 * t470;
t368 = t403 * t490;
t533 = -t529 * t424 - t368;
t173 = qJD(4) * t533 - t335 * t407 + t351;
t353 = t331 * t335 + t332 * t403;
t354 = t331 * t427 + t332 * t529;
t532 = qJD(3) * t353 + qJD(4) * t354;
t338 = t532 * t326;
t404 = t331 * t424;
t174 = qJD(4) * t404 + t338;
t555 = pkin(4) * t174 + pkin(12) * t173 + t481;
t131 = -t176 * t325 + t501 * t258;
t245 = t353 * t472;
t554 = -pkin(4) * t245 + pkin(12) * t247 - t131 + (pkin(4) * t331 - pkin(12) * t529) * t325 * qJD(4);
t442 = t332 * t472;
t413 = t325 * t442;
t553 = -pkin(12) * t413 + t508;
t371 = t529 * t416;
t401 = qJD(4) * t431;
t438 = t325 * t469;
t552 = qJD(1) * t371 + t218 * t437 + t229 * t401 + t254 * t438 + t331 * t559 + t558 * t403;
t330 = sin(qJ(5));
t334 = cos(qJ(5));
t277 = t330 * t492 - t334 * t501;
t480 = qJD(5) * t277 + t330 * t413 - t334 * t557;
t278 = t330 * t501 + t334 * t492;
t479 = qJD(5) * t278 + t330 * t557 + t334 * t413;
t544 = t245 - t438;
t449 = t325 * t529;
t281 = pkin(3) * t431 + pkin(11) * t449;
t551 = t281 * qJD(4) - t331 * t175 + t176 * t403;
t411 = t331 * t442;
t421 = t502 * qJD(2);
t387 = t421 + qJD(3);
t367 = t325 * t387;
t543 = -qJD(2) * t368 - t529 * t367;
t204 = t411 + t543;
t355 = qJD(5) + t204;
t264 = pkin(12) * t501 + t281;
t265 = (-pkin(4) * t529 - pkin(12) * t331 - pkin(3)) * t325;
t466 = qJD(5) * t334;
t467 = qJD(5) * t330;
t514 = -t264 * t467 + t265 * t466 + t330 * t554 + t334 * t553;
t169 = -t229 * t325 + t501 * t254;
t222 = t453 + t533;
t345 = t354 * t326;
t225 = t404 + t345;
t103 = pkin(4) * t222 - pkin(12) * t225 + t169;
t116 = t529 * t218 + t229 * t431 + t254 * t492;
t395 = t501 * t502;
t456 = t325 * t490;
t273 = -t395 + t456;
t108 = -pkin(12) * t273 + t116;
t512 = t103 * t466 - t108 * t467 + t555 * t330 - t334 * t556;
t446 = t529 * t258;
t507 = -(-pkin(4) * t442 - t446) * t325 + t551;
t445 = t529 * t259;
t505 = (-pkin(4) * t440 - t445) * t325 + t552;
t550 = -pkin(13) * t174 - t512;
t549 = -pkin(13) * t544 + t514;
t548 = -pkin(5) * t479 - pkin(13) * t480 - t507;
t180 = t225 * t334 - t273 * t330;
t100 = qJD(5) * t180 - t173 * t330 - t334 * t410;
t99 = t334 * t173 + t225 * t467 + t273 * t466 - t330 * t410;
t547 = t100 * pkin(5) + t99 * pkin(13) + t505;
t422 = qJD(1) * t502;
t381 = t422 * t489;
t399 = qJD(3) * t425;
t439 = t335 * t470;
t408 = t328 * t439;
t488 = t327 * t336;
t441 = qJD(2) * t488;
t451 = t335 * t399 + (t335 * t441 + t408) * qJD(1);
t471 = qJD(3) * t287;
t129 = (-t471 + (-qJD(3) * t415 - t381) * qJD(2)) * t332 + t451;
t349 = t358 * qJD(2);
t344 = t327 * t349;
t130 = qJD(1) * t344 + qJD(3) * t176;
t190 = t332 * t541 + t484;
t155 = (t326 * t402 + t367) * pkin(11) + t190;
t161 = pkin(3) * t387 + t175;
t311 = t328 * t422;
t194 = t311 + (qJD(2) * t379 - t299) * t326;
t214 = (t364 + t443) * t472;
t36 = -t331 * t129 + t130 * t403 - t155 * t437 - t161 * t401 - t194 * t438 + t214 * t449;
t405 = qJD(2) * t440;
t388 = t325 * t405;
t34 = -pkin(4) * t388 - t36;
t348 = qJD(5) * t355;
t546 = pkin(12) * t348 + t34;
t81 = t529 * t155 + (t161 * t501 + t194 * t325) * t331;
t545 = -t81 + t355 * (pkin(5) * t330 - pkin(13) * t334);
t464 = qJD(2) * t456 - qJD(4);
t342 = -t387 * t501 + t464;
t69 = -pkin(12) * t342 + t81;
t109 = -t161 * t325 + t501 * t194;
t356 = t331 * t367;
t206 = qJD(2) * t345 + t356;
t79 = pkin(4) * t204 - pkin(12) * t206 + t109;
t37 = -t330 * t69 + t334 * t79;
t352 = -t529 * t129 - t130 * t431 + t155 * t469 - t161 * t378 - t194 * t406 - t214 * t492;
t33 = pkin(12) * t388 - t352;
t384 = qJD(2) * t407;
t450 = qJD(4) * t543 - t335 * t384;
t159 = qJD(2) * t351 + t450;
t350 = qJD(4) * t356;
t160 = qJD(2) * t338 + t350;
t98 = -t130 * t325 + t501 * t214;
t59 = pkin(4) * t160 + pkin(12) * t159 + t98;
t370 = -t334 * t33 - t330 * t59 - t79 * t466 + t467 * t69;
t542 = -t355 * t37 - t370;
t38 = t330 * t79 + t334 * t69;
t31 = pkin(13) * t355 + t38;
t329 = sin(qJ(6));
t241 = t334 * t342;
t156 = t206 * t330 + t241;
t158 = t334 * t206 - t330 * t342;
t80 = -t331 * t155 + t161 * t403 + t194 * t449;
t68 = pkin(4) * t342 - t80;
t39 = t156 * pkin(5) - t158 * pkin(13) + t68;
t528 = cos(qJ(6));
t12 = t31 * t528 + t329 * t39;
t154 = qJD(6) + t156;
t88 = qJD(5) * t241 + t334 * t159 + t206 * t467 - t330 * t388;
t420 = -t330 * t159 - t334 * t388;
t468 = qJD(5) * t158;
t89 = t420 + t468;
t14 = pkin(5) * t89 + pkin(13) * t88 + t34;
t5 = pkin(13) * t160 - t370;
t2 = -qJD(6) * t12 + t528 * t14 - t329 * t5;
t539 = -t12 * t154 - t2;
t534 = t334 * t264 + t330 * t265;
t513 = -qJD(5) * t534 - t330 * t553 + t334 * t554;
t535 = t330 * t103 + t334 * t108;
t511 = -qJD(5) * t535 + t330 * t556 + t555 * t334;
t506 = t325 * t446 + t551;
t503 = t325 * t445 - t552;
t537 = t156 * t355;
t536 = t158 * t355;
t271 = qJD(3) * t282;
t478 = -t246 - t271;
t477 = pkin(10) * t440 + t560;
t433 = t330 * t33 - t334 * t59;
t8 = -qJD(5) * t38 - t433;
t114 = t158 * t528 + t329 * t355;
t41 = qJD(6) * t114 - t528 * t160 - t329 * t88;
t373 = t329 * t31 - t39 * t528;
t1 = -qJD(6) * t373 + t329 * t14 + t5 * t528;
t337 = qJD(2) ^ 2;
t53 = pkin(13) * t222 + t535;
t115 = -t331 * t218 + t229 * t403 + t254 * t449;
t107 = t273 * pkin(4) - t115;
t179 = t225 * t330 + t334 * t273;
t67 = t179 * pkin(5) - t180 * pkin(13) + t107;
t25 = -t329 * t53 + t528 * t67;
t531 = qJD(6) * t25 + t547 * t329 - t528 * t550;
t26 = t329 * t67 + t528 * t53;
t530 = -qJD(6) * t26 + t329 * t550 + t547 * t528;
t6 = -pkin(5) * t160 - t8;
t526 = t6 * t329;
t263 = -pkin(4) * t501 - t279;
t181 = t277 * pkin(5) - t278 * pkin(13) + t263;
t183 = -pkin(13) * t449 + t534;
t111 = t329 * t181 + t183 * t528;
t525 = qJD(6) * t111 + t329 * t549 + t528 * t548;
t110 = t181 * t528 - t329 * t183;
t524 = -qJD(6) * t110 + t329 * t548 - t528 * t549;
t521 = t179 * t89;
t520 = t277 * t89;
t519 = t329 * t89;
t518 = t334 * t89;
t343 = t528 * t355;
t465 = qJD(6) * t329;
t40 = -qJD(6) * t343 + t158 * t465 - t329 * t160 + t528 * t88;
t517 = t40 * t329;
t516 = -pkin(5) * t174 - t511;
t515 = pkin(5) * t544 - t513;
t308 = -pkin(5) * t334 - pkin(13) * t330 - pkin(4);
t435 = qJD(6) * t528;
t436 = qJD(5) * t528;
t136 = pkin(4) * t206 + pkin(12) * t204;
t55 = t330 * t136 + t334 * t80;
t45 = pkin(13) * t206 + t55;
t510 = -t45 * t528 + t308 * t435 + (-t330 * t436 - t334 * t465) * pkin(12) + t545 * t329;
t509 = t329 * t45 - t308 * t465 + (t329 * t467 - t334 * t435) * pkin(12) + t545 * t528;
t112 = t158 * t329 - t343;
t500 = t112 * t329;
t499 = t114 * t112;
t498 = t114 * t154;
t497 = t158 * t156;
t496 = t160 * t222;
t495 = t204 * t206;
t250 = -t299 * t326 + t311;
t494 = t250 * t326;
t322 = t326 ^ 2;
t493 = t322 * t337;
t487 = t327 * t337;
t486 = t329 * t334;
t485 = t332 * t287;
t234 = t329 * t278 + t449 * t528;
t483 = qJD(6) * t234 + t329 * t544 + t480 * t528;
t414 = t329 * t449;
t482 = -qJD(6) * t414 + t278 * t435 - t329 * t480 + t528 * t544;
t475 = t332 ^ 2 - t335 ^ 2;
t463 = qJD(2) * qJD(3);
t462 = t1 * t528;
t461 = t6 * t528;
t460 = t322 * t523;
t459 = t40 * t528;
t458 = t41 * t528;
t457 = t528 * t89;
t452 = t333 * t487;
t448 = t334 * t528;
t447 = t529 * t160;
t434 = t335 * t463;
t426 = t336 * t502;
t419 = t355 * t330;
t418 = t322 * t452;
t417 = t332 * t335 * t493;
t412 = qJD(2) * t454;
t409 = t328 * t440;
t400 = t326 * t337 * t502;
t132 = -t204 * t486 - t206 * t528;
t396 = -t329 * t466 + t132;
t54 = t136 * t334 - t330 * t80;
t62 = t103 * t334 - t108 * t330;
t357 = -t332 * t333 + t335 * t426;
t223 = t327 * t357 + t328 * t490;
t359 = t332 * t426 + t333 * t335;
t224 = t327 * t359 + t328 * t491;
t274 = -t326 * t488 + t328 * t502;
t122 = t224 * t529 + (t223 * t501 + t274 * t325) * t331;
t178 = -t223 * t325 + t274 * t501;
t91 = t122 * t334 + t178 * t330;
t90 = t122 * t330 - t178 * t334;
t191 = -t330 * t264 + t334 * t265;
t389 = t322 * t332 * t434;
t386 = 0.2e1 * t421 + qJD(3);
t133 = -t204 * t448 + t329 * t206;
t376 = t334 * t436 - t133;
t374 = -t12 * t329 + t373 * t528;
t372 = -t460 + t494;
t340 = t223 * t403 - t224 * t331 + t274 * t449;
t60 = -t91 * t329 - t340 * t528;
t61 = -t329 * t340 + t528 * t91;
t126 = t180 * t528 + t329 * t222;
t369 = -pkin(12) * t160 + t355 * t68;
t363 = t154 * t465 - t457;
t362 = t156 * t528 + t435;
t347 = -qJD(2) * t381 - t471;
t341 = qJD(4) * t342;
t339 = t204 * t355 + t348;
t321 = t325 ^ 2;
t280 = -pkin(10) * t491 + t318;
t261 = pkin(12) * t448 + t329 * t308;
t260 = -pkin(12) * t486 + t308 * t528;
t235 = t278 * t528 - t414;
t189 = t476 - t485;
t182 = pkin(5) * t449 - t191;
t172 = t408 + (qJD(2) * t360 + qJD(3) * t357) * t327;
t171 = -t409 + (-qJD(3) * t359 + t349) * t327;
t140 = qJD(2) * t383 - t171 * t325;
t139 = t361 * qJD(3) + (t344 - t409) * qJD(1);
t138 = t332 * t347 + t451;
t125 = t180 * t329 - t222 * t528;
t95 = pkin(5) * t158 + pkin(13) * t156;
t66 = t172 * t529 + (t171 * t501 + t325 * t412) * t331 + t340 * qJD(4);
t65 = -qJD(2) * t371 + qJD(4) * t122 - t171 * t403 + t172 * t331;
t52 = -pkin(5) * t222 - t62;
t51 = qJD(6) * t126 - t174 * t528 - t329 * t99;
t50 = -t329 * t174 + t180 * t465 - t222 * t435 + t528 * t99;
t44 = -pkin(5) * t206 - t54;
t30 = -pkin(5) * t355 - t37;
t28 = -qJD(5) * t90 + t140 * t330 + t334 * t66;
t27 = qJD(5) * t91 - t140 * t334 + t330 * t66;
t24 = t329 * t95 + t37 * t528;
t23 = -t329 * t37 + t528 * t95;
t10 = qJD(6) * t60 + t28 * t528 + t65 * t329;
t9 = -qJD(6) * t61 - t28 * t329 + t528 * t65;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t452, -t336 * t487, 0, 0, 0, 0, 0, 0, 0, 0, t171 * t387 + t274 * t405 - t335 * t418, t274 * t326 * t434 - t172 * t387 + t332 * t418 (-t171 * t332 + t172 * t335 + (-t223 * t335 - t224 * t332) * qJD(3)) * t472, t138 * t224 + t139 * t223 + t171 * t189 + t172 * t190 + (qJD(1) * t274 + t250) * t412, 0, 0, 0, 0, 0, 0, t140 * t204 + t178 * t160 + t340 * t388 + t342 * t65, -t122 * t388 + t140 * t206 - t178 * t159 + t342 * t66, -t122 * t160 + t159 * t340 - t204 * t66 + t206 * t65, t109 * t140 - t122 * t352 + t178 * t98 + t340 * t36 - t65 * t80 + t66 * t81, 0, 0, 0, 0, 0, 0, t65 * t156 - t90 * t160 - t27 * t355 - t340 * t89, t65 * t158 - t91 * t160 - t28 * t355 + t340 * t88, -t156 * t28 + t158 * t27 - t88 * t90 - t89 * t91, -t27 * t37 + t28 * t38 - t34 * t340 - t370 * t91 + t65 * t68 - t8 * t90, 0, 0, 0, 0, 0, 0, t112 * t27 + t154 * t9 + t41 * t90 + t60 * t89, -t10 * t154 + t114 * t27 - t40 * t90 - t61 * t89, -t10 * t112 - t114 * t9 + t40 * t60 - t41 * t61, t1 * t61 + t10 * t12 + t2 * t60 + t27 * t30 - t373 * t9 + t6 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t389, -0.2e1 * t475 * t322 * t463, t386 * t439, -0.2e1 * t389, -t386 * t440, 0 (t332 * t372 + t478) * qJD(3) + (qJD(2) * t478 + t139) * t502 (t335 * t372 + t477) * qJD(3) + (qJD(2) * t477 - t138) * t502 (t138 * t335 - t139 * t332 + (-t189 * t335 - t190 * t332) * qJD(3) + ((-qJD(3) * t280 - t477) * t335 + (-t478 - t271) * t332) * qJD(2)) * t326, t138 * t282 + t139 * t280 - t477 * t190 + t478 * t189 + (-t460 - t494) * t443, -t159 * t225 - t173 * t206, t159 * t222 - t160 * t225 + t173 * t204 - t174 * t206, t273 * t159 + t342 * t173 + (qJD(2) * t225 + t206) * t410, t174 * t204 + t496, t273 * t160 + t342 * t174 + (-qJD(2) * t222 - t204) * t410 (-qJD(3) * t342 - t273 * t463) * t325 * t491, t169 * t160 + t98 * t222 + t109 * t174 - t36 * t273 + (qJD(2) * t115 + t80) * t410 + t481 * t204 - t503 * t342, -t169 * t159 + t98 * t225 - t109 * t173 - t352 * t273 + (-qJD(2) * t116 - t81) * t410 + t481 * t206 + t504 * t342, t115 * t159 - t116 * t160 + t173 * t80 - t174 * t81 - t204 * t504 - t206 * t503 + t222 * t352 - t225 * t36, t109 * t481 + t115 * t36 - t116 * t352 + t169 * t98 + t503 * t80 + t504 * t81, -t158 * t99 - t180 * t88, -t100 * t158 + t156 * t99 + t179 * t88 - t180 * t89, t158 * t174 + t180 * t160 - t88 * t222 - t355 * t99, t100 * t156 + t521, -t100 * t355 - t156 * t174 - t179 * t160 - t89 * t222, t174 * t355 + t496, t68 * t100 + t107 * t89 + t156 * t505 + t62 * t160 + t37 * t174 + t34 * t179 + t8 * t222 + t355 * t511, -t107 * t88 + t158 * t505 - t160 * t535 - t38 * t174 + t34 * t180 + t222 * t370 - t355 * t512 - t68 * t99, -t100 * t38 - t156 * t512 - t158 * t511 + t179 * t370 - t180 * t8 + t37 * t99 - t535 * t89 + t62 * t88, t107 * t34 + t37 * t511 - t370 * t535 + t38 * t512 + t505 * t68 + t62 * t8, -t114 * t50 - t126 * t40, t112 * t50 - t114 * t51 + t125 * t40 - t126 * t41, t100 * t114 + t126 * t89 - t154 * t50 - t179 * t40, t112 * t51 + t125 * t41, -t100 * t112 - t125 * t89 - t154 * t51 - t179 * t41, t100 * t154 + t521, -t100 * t373 + t112 * t516 + t125 * t6 + t154 * t530 + t179 * t2 + t25 * t89 + t30 * t51 + t41 * t52, -t1 * t179 - t100 * t12 + t114 * t516 + t126 * t6 - t154 * t531 - t26 * t89 - t30 * t50 - t40 * t52, -t1 * t125 - t112 * t531 - t114 * t530 - t12 * t51 - t126 * t2 + t25 * t40 - t26 * t41 - t373 * t50, t1 * t26 + t12 * t531 + t2 * t25 + t30 * t516 - t373 * t530 + t52 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t417, t475 * t493, -t335 * t400, t417, t332 * t400, 0, t190 * t387 + t347 * t335 + (-t399 - t250 * t472 + (-t328 * t470 - t441) * qJD(1)) * t332 (t189 + t485) * qJD(3) + (t189 * t502 - t250 * t490 + t332 * t381) * qJD(2) - t451, 0, 0, -t159 * t492 + t206 * t557, t204 * t247 + t245 * t206 + (-t529 * t159 - t160 * t331 + (-t204 * t529 - t206 * t331) * qJD(4)) * t325, -t501 * t159 + t321 * t331 * t405 + t342 * t247 + (-t206 * t442 - t341 * t529) * t325, -t204 * t544 - t447 * t325, -t501 * t160 + t321 * t332 * t384 - t342 * t245 + (t204 * t442 + t331 * t341) * t325 -(qJD(2) * t395 - t464) * t413, t36 * t501 - t131 * t204 - t109 * t245 + (t109 * t469 - t98 * t529 - pkin(3) * t160 + (qJD(3) * t279 - t80) * t442) * t325 + t506 * t342, t352 * t501 - t131 * t206 - t109 * t247 + (t109 * t437 + pkin(3) * t159 + t98 * t331 + (-qJD(3) * t281 + t81) * t442) * t325 + t508 * t342, t279 * t159 - t281 * t160 + t81 * t245 + t80 * t247 + t506 * t206 - t508 * t204 + (-t529 * t352 - t331 * t36 + (-t331 * t81 - t529 * t80) * qJD(4)) * t325, -pkin(3) * t325 * t98 - t109 * t131 + t279 * t36 - t281 * t352 - t506 * t80 + t508 * t81, -t158 * t480 - t278 * t88, t156 * t480 - t158 * t479 + t277 * t88 - t278 * t89, -t158 * t544 + t278 * t160 - t355 * t480 + t449 * t88, t156 * t479 + t520, t156 * t544 - t277 * t160 - t355 * t479 + t449 * t89, -t355 * t245 + (t355 * t469 - t447) * t325, t191 * t160 + t263 * t89 + t34 * t277 - t37 * t245 + t479 * t68 + (t37 * t469 - t529 * t8) * t325 + t507 * t156 + t513 * t355, -t534 * t160 - t263 * t88 + t34 * t278 + t38 * t245 - t480 * t68 + (-t370 * t529 - t38 * t469) * t325 + t507 * t158 - t514 * t355, -t156 * t514 - t158 * t513 + t191 * t88 + t277 * t370 - t278 * t8 + t37 * t480 - t38 * t479 - t534 * t89, t191 * t8 + t263 * t34 + t37 * t513 - t370 * t534 + t38 * t514 + t507 * t68, -t114 * t483 - t235 * t40, t112 * t483 - t114 * t482 + t234 * t40 - t235 * t41, t114 * t479 - t154 * t483 + t235 * t89 - t277 * t40, t112 * t482 + t234 * t41, -t112 * t479 - t154 * t482 - t234 * t89 - t277 * t41, t154 * t479 + t520, t110 * t89 + t112 * t515 - t154 * t525 + t182 * t41 + t2 * t277 + t234 * t6 + t30 * t482 - t373 * t479, -t1 * t277 - t111 * t89 + t114 * t515 - t12 * t479 + t154 * t524 - t182 * t40 + t235 * t6 - t30 * t483, -t1 * t234 + t110 * t40 - t111 * t41 + t112 * t524 + t114 * t525 - t12 * t482 - t2 * t235 - t373 * t483, t1 * t111 + t110 * t2 - t12 * t524 + t182 * t6 + t30 * t515 + t373 * t525; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t495, -t204 ^ 2 + t206 ^ 2, -t204 * t342 - t385 * t411 - t450, -t495, -t342 * t206 - t472 * t532 - t350, t388, -t109 * t206 - t342 * t81 + t36, t109 * t204 - t342 * t80 + t352, 0, 0, -t330 * t88 + t334 * t536 (-t88 - t537) * t334 + (-t89 - t536) * t330, -t158 * t206 + t330 * t160 + t334 * t339, t156 * t419 - t518, t156 * t206 + t334 * t160 - t330 * t339, -t355 * t206, -pkin(4) * t89 - t81 * t156 - t37 * t206 + t369 * t330 - t334 * t546 - t54 * t355, pkin(4) * t88 - t81 * t158 + t38 * t206 + t330 * t546 + t369 * t334 + t55 * t355, t156 * t55 + t158 * t54 + ((-t89 + t468) * pkin(12) + t542) * t334 + (-t8 - t355 * t38 + (qJD(5) * t156 - t88) * pkin(12)) * t330, -pkin(4) * t34 - t37 * t54 - t38 * t55 - t68 * t81 + (-t330 * t8 - t334 * t370 + (-t330 * t38 - t334 * t37) * qJD(5)) * pkin(12), -t330 * t459 + (-t330 * t465 + t376) * t114, t133 * t112 + t114 * t132 + (-t112 * t528 - t114 * t329) * t466 + (-t458 + t517 + (-t114 * t528 + t500) * qJD(6)) * t330, t40 * t334 + t376 * t154 + (t114 * t355 - t363) * t330, t41 * t329 * t330 + (t330 * t435 - t396) * t112, t41 * t334 + t396 * t154 + (-t112 * t355 - t154 * t435 - t519) * t330, t154 * t419 - t518, -t44 * t112 - t30 * t132 + t260 * t89 + t509 * t154 + (-t2 + (pkin(12) * t112 + t30 * t329) * qJD(5)) * t334 + (pkin(12) * t41 + t30 * t435 - t355 * t373 + t526) * t330, -t44 * t114 - t30 * t133 - t261 * t89 - t510 * t154 + (t1 + (pkin(12) * t114 + t30 * t528) * qJD(5)) * t334 + (-pkin(12) * t40 - t12 * t355 - t30 * t465 + t461) * t330, -t373 * t133 + t12 * t132 + t260 * t40 - t261 * t41 - t509 * t114 - t510 * t112 + t374 * t466 + (-t528 * t2 - t1 * t329 + (-t12 * t528 - t329 * t373) * qJD(6)) * t330, t1 * t261 + t2 * t260 - t30 * t44 + t510 * t12 - t509 * t373 + (t30 * t466 + t330 * t6) * pkin(12); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t497, -t156 ^ 2 + t158 ^ 2, -t88 + t537, -t497, t158 * t204 - t420, t160, -t68 * t158 + t204 * t38 - t433, t68 * t156 - t542, 0, 0, t114 * t362 - t517, -t459 - t362 * t112 + (-t41 - t498) * t329, -t114 * t158 + t154 * t362 + t519, t154 * t500 - t458, -t154 ^ 2 * t329 + t112 * t158 + t457, -t154 * t158, -t461 - pkin(5) * t41 + t373 * t158 - t38 * t112 + (-pkin(13) * t435 - t23) * t154 + (-pkin(13) * t89 + t154 * t30) * t329, pkin(5) * t40 + pkin(13) * t363 - t38 * t114 + t12 * t158 + t24 * t154 + t30 * t362 + t526, t462 + t24 * t112 + t23 * t114 + t362 * t373 + (t114 * t435 - t458) * pkin(13) + ((qJD(6) * t112 - t40) * pkin(13) + t539) * t329, -t6 * pkin(5) + t373 * t23 - t12 * t24 - t30 * t38 + (qJD(6) * t374 - t2 * t329 + t462) * pkin(13); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t499, -t112 ^ 2 + t114 ^ 2, t112 * t154 - t40, -t499, -t41 + t498, t89, -t30 * t114 - t539, t30 * t112 - t154 * t373 - t1, 0, 0;];
tauc_reg  = t3;
