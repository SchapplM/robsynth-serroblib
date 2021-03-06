% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RRRPRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRPRP6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP6_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_invdynJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:01:17
% EndTime: 2019-03-09 17:01:55
% DurationCPUTime: 21.77s
% Computational Cost: add. (26417->888), mult. (65357->1134), div. (0->0), fcn. (52720->14), ass. (0->393)
t359 = sin(pkin(6));
t369 = cos(qJ(2));
t489 = qJD(1) * t369;
t335 = t359 * t489;
t409 = t335 - qJD(3);
t365 = sin(qJ(2));
t491 = qJD(1) * t359;
t458 = t365 * t491;
t360 = cos(pkin(6));
t490 = qJD(1) * t360;
t474 = pkin(1) * t490;
t278 = -pkin(8) * t458 + t369 * t474;
t419 = pkin(2) * t365 - pkin(9) * t369;
t279 = t419 * t491;
t364 = sin(qJ(3));
t368 = cos(qJ(3));
t200 = -t278 * t364 + t368 * t279;
t362 = -qJ(4) - pkin(9);
t447 = qJD(3) * t362;
t605 = (-qJ(4) * t368 * t369 + pkin(3) * t365) * t491 + t200 + qJD(4) * t364 - t368 * t447;
t201 = t368 * t278 + t364 * t279;
t427 = t364 * t335;
t604 = qJ(4) * t427 + qJD(4) * t368 + t364 * t447 - t201;
t358 = sin(pkin(11));
t541 = cos(pkin(11));
t445 = t541 * t364;
t299 = t358 * t368 + t445;
t496 = t409 * t299;
t513 = t358 * t364;
t393 = t368 * t541 - t513;
t235 = t393 * t335;
t286 = t393 * qJD(3);
t495 = t286 - t235;
t499 = -t358 * t605 + t604 * t541;
t509 = t359 * t369;
t570 = pkin(1) * t365;
t494 = pkin(8) * t509 + t360 * t570;
t281 = t494 * qJD(1);
t486 = qJD(3) * t364;
t420 = -t281 + (-t427 + t486) * pkin(3);
t480 = qJD(1) * qJD(2);
t451 = t369 * t480;
t478 = qJDD(1) * t365;
t603 = t451 + t478;
t363 = sin(qJ(5));
t367 = cos(qJ(5));
t336 = qJD(2) + t490;
t240 = pkin(9) * t336 + t281;
t404 = -pkin(2) * t369 - pkin(9) * t365 - pkin(1);
t272 = t404 * t359;
t245 = qJD(1) * t272;
t163 = -t364 * t240 + t368 * t245;
t260 = t336 * t364 + t368 * t458;
t140 = -t260 * qJ(4) + t163;
t129 = -pkin(3) * t409 + t140;
t164 = t240 * t368 + t245 * t364;
t426 = t364 * t458;
t259 = t336 * t368 - t426;
t141 = qJ(4) * t259 + t164;
t446 = t541 * t141;
t85 = t358 * t129 + t446;
t80 = -pkin(10) * t409 + t85;
t239 = -pkin(2) * t336 - t278;
t177 = -pkin(3) * t259 + qJD(4) + t239;
t179 = -t541 * t259 + t260 * t358;
t394 = t358 * t259 + t260 * t541;
t96 = pkin(4) * t179 - pkin(10) * t394 + t177;
t42 = -t363 * t80 + t367 * t96;
t591 = qJD(5) + t179;
t602 = t42 * t591;
t43 = t363 * t96 + t367 * t80;
t601 = t43 * t591;
t433 = t367 * t409;
t149 = t363 * t394 + t433;
t30 = -qJ(6) * t149 + t43;
t600 = t591 * t30;
t599 = pkin(10) * t458 - t499;
t598 = -t496 * pkin(4) - t495 * pkin(10) + t420;
t205 = t235 * t363 - t367 * t458;
t440 = -t286 * t363 + t205;
t512 = t359 * t365;
t287 = -t360 * t368 + t364 * t512;
t597 = pkin(3) * t287;
t596 = t179 * t394;
t500 = t604 * t358 + t541 * t605;
t483 = qJD(5) * t367;
t455 = t299 * t483;
t383 = t455 - t440;
t479 = qJDD(1) * t360;
t425 = qJDD(2) + t479;
t430 = qJD(2) * t474;
t472 = pkin(1) * t479;
t594 = t603 * t359;
t428 = t594 * pkin(8) + t365 * t430 - t369 * t472;
t187 = -pkin(2) * t425 + t428;
t571 = cos(qJ(1));
t460 = t571 * t369;
t366 = sin(qJ(1));
t505 = t365 * t366;
t289 = -t360 * t460 + t505;
t461 = t571 * t365;
t504 = t366 * t369;
t291 = t360 * t504 + t461;
t577 = -g(1) * t291 - g(2) * t289 + g(3) * t509;
t375 = t187 + t577;
t431 = qJD(3) * t409;
t595 = -pkin(9) * t431 + t375;
t292 = -t360 * t505 + t460;
t510 = t359 * t368;
t231 = -t292 * t364 + t366 * t510;
t485 = qJD(3) * t368;
t165 = qJD(3) * t426 - t336 * t485 - t364 * t425 - t594 * t368;
t487 = qJD(2) * t369;
t456 = t364 * t487;
t166 = t336 * t486 + t359 * (qJD(1) * (t365 * t485 + t456) + t364 * t478) - t368 * t425;
t443 = t358 * t165 - t541 * t166;
t120 = qJDD(5) - t443;
t434 = t367 * t591;
t593 = -t120 * t363 - t434 * t591;
t452 = t365 * t480;
t424 = t359 * t452;
t477 = qJDD(1) * t369;
t334 = t359 * t477;
t464 = pkin(8) * t334 + t365 * t472 + t369 * t430;
t212 = -pkin(8) * t424 + t464;
t186 = pkin(9) * t425 + t212;
t399 = t419 * qJD(2);
t194 = (qJD(1) * t399 + qJDD(1) * t404) * t359;
t390 = -t368 * t186 - t364 * t194 + t240 * t486 - t245 * t485;
t592 = t409 * t163 - t390;
t337 = pkin(8) * t512;
t569 = pkin(1) * t369;
t293 = t360 * t569 - t337;
t282 = qJD(2) * t293;
t290 = t360 * t461 + t504;
t355 = qJ(3) + pkin(11);
t351 = sin(t355);
t352 = cos(t355);
t463 = t359 * t571;
t220 = t290 * t351 + t352 * t463;
t511 = t359 * t366;
t224 = t292 * t351 - t352 * t511;
t268 = t351 * t512 - t360 * t352;
t387 = g(1) * t224 + g(2) * t220 + g(3) * t268;
t221 = t290 * t352 - t351 * t463;
t590 = t221 * t363 - t289 * t367;
t525 = t289 * t363;
t589 = t221 * t367 + t525;
t354 = t359 ^ 2;
t588 = 0.2e1 * t354;
t558 = g(3) * t359;
t567 = pkin(3) * t358;
t347 = pkin(10) + t567;
t502 = qJ(6) + t347;
t437 = qJD(5) * t502;
t568 = pkin(3) * t260;
t114 = pkin(4) * t394 + pkin(10) * t179 + t568;
t136 = t358 * t141;
t89 = t140 * t541 - t136;
t50 = t363 * t114 + t367 * t89;
t529 = t179 * t363;
t586 = -qJ(6) * t529 + qJD(6) * t367 - t363 * t437 - t50;
t49 = t367 * t114 - t363 * t89;
t528 = t179 * t367;
t542 = -pkin(5) * t394 - qJ(6) * t528 - qJD(6) * t363 - t367 * t437 - t49;
t585 = t599 * t363 + t598 * t367;
t350 = pkin(3) * t368 + pkin(2);
t217 = -pkin(4) * t393 - pkin(10) * t299 - t350;
t584 = t217 * t483 + t598 * t363 - t599 * t367;
t484 = qJD(5) * t363;
t583 = t347 * t484 + t50;
t581 = t179 * t409;
t580 = t259 * t409;
t579 = t260 * t409;
t374 = t577 * t351;
t501 = pkin(4) * t458 + t500;
t323 = t362 * t368;
t238 = -t323 * t541 + t362 * t513;
t227 = t367 * t238;
t146 = t363 * t217 + t227;
t271 = pkin(9) * t360 + t494;
t189 = t368 * t271 + t364 * t272;
t283 = t494 * qJD(2);
t225 = t292 * t352 + t351 * t511;
t171 = -t225 * t363 + t291 * t367;
t269 = t351 * t360 + t352 * t512;
t503 = t367 * t369;
t467 = t359 * t503;
t576 = -g(3) * (-t269 * t363 - t467) + g(2) * t590 - g(1) * t171;
t151 = -t363 * t409 + t367 * t394;
t476 = qJDD(3) - t334;
t385 = t424 + t476;
t442 = -t364 * t186 + t368 * t194;
t95 = -qJD(3) * t164 + t442;
t64 = pkin(3) * t385 + t165 * qJ(4) - t260 * qJD(4) + t95;
t67 = -qJ(4) * t166 + qJD(4) * t259 - t390;
t24 = t358 * t64 + t541 * t67;
t20 = pkin(10) * t385 + t24;
t125 = t166 * pkin(3) + qJDD(4) + t187;
t395 = t165 * t541 + t358 * t166;
t44 = -pkin(4) * t443 + pkin(10) * t395 + t125;
t4 = -qJD(5) * t43 - t20 * t363 + t367 * t44;
t74 = qJD(5) * t433 - t363 * t385 + t367 * t395 + t394 * t484;
t549 = qJ(6) * t74;
t566 = pkin(5) * t120;
t1 = -qJD(6) * t151 + t4 + t549 + t566;
t575 = t1 + t600;
t574 = -pkin(9) * t385 - t409 * t239;
t75 = qJD(5) * t151 - t363 * t395 - t367 * t385;
t572 = t151 ^ 2;
t302 = t350 * t509;
t559 = g(3) * t302;
t557 = t149 * pkin(5);
t556 = t367 * pkin(5);
t29 = -qJ(6) * t151 + t42;
t25 = pkin(5) * t591 + t29;
t555 = -t29 + t25;
t554 = t358 * t67 - t541 * t64;
t206 = t235 * t367 + t363 * t458;
t406 = -qJ(6) * t286 - qJD(6) * t299;
t553 = qJ(6) * t206 + t406 * t367 + (-t227 + (qJ(6) * t299 - t217) * t363) * qJD(5) + t585 - t496 * pkin(5);
t552 = (-qJD(5) * t238 + t406) * t363 + t584 + (t205 - t455) * qJ(6);
t551 = -t238 * t484 + t584;
t550 = -t146 * qJD(5) + t585;
t548 = qJ(6) * t75;
t547 = t363 * t74;
t546 = t367 * t75;
t280 = t359 * t399;
t131 = -t189 * qJD(3) + t368 * t280 - t282 * t364;
t229 = -qJD(3) * t287 + t487 * t510;
t288 = t360 * t364 + t365 * t510;
t488 = qJD(2) * t365;
t457 = t359 * t488;
t101 = pkin(3) * t457 - qJ(4) * t229 - qJD(4) * t288 + t131;
t130 = -t271 * t486 + t272 * t485 + t364 * t280 + t368 * t282;
t228 = qJD(3) * t288 + t359 * t456;
t106 = -qJ(4) * t228 - qJD(4) * t287 + t130;
t55 = t358 * t101 + t541 * t106;
t545 = t383 * pkin(5) + t501;
t544 = -t149 * t483 - t363 * t75;
t188 = -t271 * t364 + t368 * t272;
t144 = -pkin(3) * t509 - qJ(4) * t288 + t188;
t158 = -qJ(4) * t287 + t189;
t105 = t358 * t144 + t541 * t158;
t100 = -pkin(10) * t509 + t105;
t208 = t287 * t541 + t288 * t358;
t209 = -t358 * t287 + t288 * t541;
t270 = t337 + (-pkin(2) - t569) * t360;
t216 = t270 + t597;
t124 = pkin(4) * t208 - pkin(10) * t209 + t216;
t59 = t367 * t100 + t363 * t124;
t539 = t149 * t591;
t538 = t149 * t394;
t537 = t149 * t179;
t536 = t151 * t149;
t535 = t151 * t591;
t534 = t151 * t394;
t533 = t151 * t363;
t532 = t591 * t394;
t531 = t394 ^ 2;
t530 = t179 ^ 2;
t527 = t260 * t259;
t523 = t290 * t363;
t522 = t290 * t364;
t521 = t291 * t363;
t520 = t292 * t363;
t518 = t299 * t363;
t517 = t299 * t367;
t516 = t352 * t363;
t515 = t352 * t367;
t514 = t354 * qJD(1) ^ 2;
t507 = t362 * t365;
t506 = t363 * t369;
t115 = t367 * t120;
t498 = -t289 * t350 - t290 * t362;
t497 = -t291 * t350 - t292 * t362;
t493 = t571 * pkin(1) + pkin(8) * t511;
t356 = t365 ^ 2;
t357 = t369 ^ 2;
t492 = t356 - t357;
t482 = qJD(2) - t336;
t471 = t369 * t514;
t469 = t364 * t511;
t325 = t359 * t506;
t465 = t387 * t363;
t462 = t368 * t571;
t459 = t541 * pkin(3);
t453 = pkin(1) * t588;
t448 = -pkin(1) * t366 + pkin(8) * t463;
t3 = t367 * t20 + t363 * t44 + t96 * t483 - t80 * t484;
t444 = qJD(6) + t557;
t58 = -t100 * t363 + t367 * t124;
t88 = t140 * t358 + t446;
t439 = -t286 * t367 + t206;
t145 = t367 * t217 - t238 * t363;
t326 = t364 * t463;
t438 = t290 * t368 - t326;
t237 = -t323 * t358 - t362 * t445;
t436 = t394 * t409;
t435 = t363 * t591;
t432 = t369 * t409;
t429 = t365 * t471;
t422 = t365 * t451;
t421 = t231 * pkin(3);
t348 = -t459 - pkin(4);
t418 = pkin(4) * t352 + pkin(10) * t351;
t417 = -g(1) * t590 - g(2) * t171;
t172 = t225 * t367 + t521;
t416 = g(1) * t589 - g(2) * t172;
t415 = -g(1) * t220 + g(2) * t224;
t413 = g(1) * t289 - g(2) * t291;
t412 = g(1) * t292 + g(2) * t290;
t2 = -qJD(6) * t149 + t3 - t548;
t410 = -t25 * t591 + t2;
t54 = t101 * t541 - t358 * t106;
t84 = t129 * t541 - t136;
t104 = t144 * t541 - t358 * t158;
t408 = pkin(3) * t469 - t291 * t362 + t292 * t350 + t493;
t349 = pkin(4) + t556;
t361 = -qJ(6) - pkin(10);
t407 = t349 * t352 - t351 * t361;
t402 = g(1) * t571 + g(2) * t366;
t401 = t115 + (-t484 - t529) * t591;
t400 = t409 * t458;
t99 = pkin(4) * t509 - t104;
t174 = t209 * t363 + t467;
t397 = pkin(3) * t326 + t289 * t362 - t290 * t350 + t448;
t52 = pkin(10) * t457 + t55;
t156 = t228 * t541 + t229 * t358;
t157 = -t358 * t228 + t229 * t541;
t183 = pkin(3) * t228 + t283;
t83 = pkin(4) * t156 - pkin(10) * t157 + t183;
t13 = -t100 * t484 + t124 * t483 + t363 * t83 + t367 * t52;
t79 = pkin(4) * t409 - t84;
t391 = -t120 * t347 + t591 * t79;
t389 = -g(1) * (t291 * t516 + t292 * t367) - g(2) * (t289 * t516 + t290 * t367) - (-t352 * t506 + t365 * t367) * t558;
t388 = -g(1) * (-t291 * t515 + t520) - g(2) * (-t289 * t515 + t523) - (t352 * t503 + t363 * t365) * t558;
t386 = g(1) * t225 + g(2) * t221 + g(3) * t269;
t384 = -t359 * t462 - t522;
t382 = -t299 * t484 - t439;
t19 = -pkin(4) * t385 + t554;
t17 = t75 * pkin(5) + qJDD(6) + t19;
t379 = -t17 + t387;
t377 = -g(3) * t512 - t412;
t376 = t384 * pkin(3);
t51 = -pkin(4) * t457 - t54;
t372 = g(1) * t172 + g(2) * t589 - g(3) * (-t269 * t367 + t325) - t3;
t14 = -qJD(5) * t59 - t363 * t52 + t367 * t83;
t371 = t4 + t576;
t311 = t348 - t556;
t296 = t502 * t367;
t295 = t502 * t363;
t232 = t292 * t368 + t469;
t195 = (-t476 * t369 + (-t335 - t409) * t488) * t359;
t193 = pkin(5) * t518 + t237;
t175 = t209 * t367 - t325;
t148 = t149 ^ 2;
t134 = -qJ(6) * t518 + t146;
t127 = -pkin(5) * t393 - qJ(6) * t517 + t145;
t109 = -qJD(5) * t325 + t157 * t363 + t209 * t483 - t367 * t457;
t108 = qJD(5) * t174 - t367 * t157 - t363 * t457;
t77 = t174 * pkin(5) + t99;
t76 = -t148 + t572;
t71 = -pkin(5) * t529 + t88;
t63 = t444 + t79;
t57 = t120 * t208 + t156 * t591;
t56 = -t120 * t393 - t496 * t591;
t47 = -qJ(6) * t174 + t59;
t46 = t535 - t75;
t45 = -t74 + t539;
t38 = pkin(5) * t208 - qJ(6) * t175 + t58;
t37 = -t534 + t593;
t36 = -t534 - t593;
t35 = t401 + t538;
t34 = t401 - t538;
t32 = t149 * t435 - t546;
t31 = t151 * t434 - t547;
t28 = t109 * pkin(5) + t51;
t27 = t109 * t149 + t174 * t75;
t26 = -t108 * t151 - t175 * t74;
t22 = t149 * t383 + t518 * t75;
t21 = t151 * t382 - t517 * t74;
t16 = -t109 * t591 - t120 * t174 - t149 * t156 - t208 * t75;
t15 = -t108 * t591 + t120 * t175 + t151 * t156 - t208 * t74;
t12 = -t120 * t518 + t149 * t496 - t383 * t591 + t393 * t75;
t11 = t115 * t299 - t151 * t496 + t382 * t591 + t393 * t74;
t10 = (-t74 - t537) * t367 - t591 * t533 + t544;
t9 = (t74 - t537) * t367 + t151 * t435 + t544;
t8 = t108 * t149 - t109 * t151 + t174 * t74 - t175 * t75;
t7 = -qJ(6) * t109 - qJD(6) * t174 + t13;
t6 = pkin(5) * t156 + qJ(6) * t108 - qJD(6) * t175 + t14;
t5 = t440 * t151 + t439 * t149 + (t547 - t546 + (t149 * t363 - t151 * t367) * qJD(5)) * t299;
t18 = [0, 0, 0, 0, 0, qJDD(1), g(1) * t366 - g(2) * t571, t402, 0, 0 (qJDD(1) * t356 + 0.2e1 * t422) * t354 (t365 * t477 - t480 * t492) * t588 (t336 * t487 + t365 * qJDD(2) + (t451 + 0.2e1 * t478) * t360) * t359 (qJDD(1) * t357 - 0.2e1 * t422) * t354 (-t336 * t488 + t369 * qJDD(2) + (-t452 + 0.2e1 * t477) * t360) * t359, t425 * t360, -t283 * t336 + t293 * t425 - t428 * t360 + g(1) * t290 - g(2) * t292 + (-t452 + t477) * t453, -t212 * t360 - t282 * t336 - t494 * t425 - t603 * t453 - t413 ((-qJD(2) * t278 + qJDD(1) * t494 + t212) * t369 + (-qJD(2) * t281 - qJDD(1) * t293 + t428) * t365 - t402) * t359, t354 * qJDD(1) * pkin(1) ^ 2 - g(1) * t448 - g(2) * t493 + t212 * t494 - t278 * t283 + t281 * t282 - t293 * t428, -t165 * t288 + t229 * t260, t165 * t287 - t166 * t288 - t228 * t260 + t229 * t259, -t229 * t409 + t288 * t476 + (t165 * t369 + (qJD(1) * t288 + t260) * t488) * t359, t166 * t287 - t228 * t259, t228 * t409 - t287 * t476 + (t166 * t369 + (-qJD(1) * t287 + t259) * t488) * t359, t195, -t131 * t409 + t188 * t476 - t283 * t259 + t270 * t166 + t187 * t287 + t239 * t228 + g(1) * t438 - g(2) * t232 + (-t95 * t369 + (qJD(1) * t188 + t163) * t488) * t359, t130 * t409 - t189 * t476 + t283 * t260 - t270 * t165 + t187 * t288 + t239 * t229 - g(1) * t522 - g(2) * t231 + (-g(1) * t462 - t390 * t369 + (-qJD(1) * t189 - t164) * t488) * t359, t130 * t259 - t131 * t260 - t163 * t229 - t164 * t228 + t165 * t188 - t166 * t189 + t287 * t390 - t288 * t95 + t413, -t390 * t189 + t164 * t130 + t95 * t188 + t163 * t131 + t187 * t270 + t239 * t283 - g(1) * (-pkin(2) * t290 - pkin(9) * t289 + t448) - g(2) * (pkin(2) * t292 + pkin(9) * t291 + t493) t157 * t394 - t209 * t395, -t156 * t394 - t157 * t179 + t208 * t395 + t209 * t443, -t157 * t409 + t209 * t476 + (t395 * t369 + (qJD(1) * t209 + t394) * t488) * t359, t156 * t179 - t208 * t443, t156 * t409 - t208 * t476 + (-t443 * t369 + (-qJD(1) * t208 - t179) * t488) * t359, t195, -t54 * t409 + t104 * t476 + t183 * t179 - t216 * t443 + t125 * t208 + t177 * t156 + g(1) * t221 - g(2) * t225 + (t554 * t369 + (qJD(1) * t104 + t84) * t488) * t359, t55 * t409 - t105 * t476 + t183 * t394 - t216 * t395 + t125 * t209 + t177 * t157 + (t24 * t369 + (-qJD(1) * t105 - t85) * t488) * t359 + t415, t104 * t395 + t105 * t443 - t85 * t156 - t84 * t157 - t55 * t179 - t24 * t208 + t209 * t554 - t394 * t54 + t413, -g(1) * t397 - g(2) * t408 - t104 * t554 + t24 * t105 + t125 * t216 + t177 * t183 + t84 * t54 + t85 * t55, t26, t8, t15, t27, t16, t57, t109 * t79 + t120 * t58 + t14 * t591 + t149 * t51 + t156 * t42 + t174 * t19 + t208 * t4 + t75 * t99 + t416, -t108 * t79 - t120 * t59 - t13 * t591 + t151 * t51 - t156 * t43 + t175 * t19 - t208 * t3 - t74 * t99 + t417, t108 * t42 - t109 * t43 - t13 * t149 - t14 * t151 - t174 * t3 - t175 * t4 + t58 * t74 - t59 * t75 - t415, t3 * t59 + t43 * t13 + t4 * t58 + t42 * t14 + t19 * t99 + t79 * t51 - g(1) * (-pkin(4) * t221 - pkin(10) * t220 + t397) - g(2) * (pkin(4) * t225 + pkin(10) * t224 + t408) t26, t8, t15, t27, t16, t57, t1 * t208 + t109 * t63 + t120 * t38 + t149 * t28 + t156 * t25 + t17 * t174 + t591 * t6 + t75 * t77 + t416, -t108 * t63 - t120 * t47 + t151 * t28 - t156 * t30 + t17 * t175 - t2 * t208 - t591 * t7 - t74 * t77 + t417, -t1 * t175 + t108 * t25 - t109 * t30 - t149 * t7 - t151 * t6 - t174 * t2 + t38 * t74 - t47 * t75 - t415, t2 * t47 + t30 * t7 + t1 * t38 + t25 * t6 + t17 * t77 + t63 * t28 - g(1) * (-pkin(5) * t525 + t220 * t361 - t221 * t349 + t397) - g(2) * (pkin(5) * t521 - t224 * t361 + t225 * t349 + t408); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t429, t492 * t514 (t482 * t489 + t478) * t359, t429, -t458 * t482 + t334, t425, t281 * t336 + t514 * t570 - t428 - t577, pkin(1) * t471 + t278 * t336 + (pkin(8) * t480 + g(3)) * t512 + t412 - t464, 0, 0, -t165 * t364 - t368 * t579 (-t165 - t580) * t368 + (-t166 + t579) * t364, -t368 * t431 + t364 * t476 + (t368 * t432 + (qJD(2) * t364 - t260) * t365) * t491, -t166 * t368 + t364 * t580, t364 * t431 + t368 * t476 + (-t364 * t432 + (qJD(2) * t368 - t259) * t365) * t491, t400, -pkin(2) * t166 - t163 * t458 + t200 * t409 + t281 * t259 + t574 * t364 - t595 * t368, pkin(2) * t165 + t164 * t458 - t201 * t409 - t281 * t260 + t595 * t364 + t574 * t368, t200 * t260 - t201 * t259 + ((qJD(3) * t260 - t166) * pkin(9) + t592) * t368 + (-t95 + t409 * t164 + (-qJD(3) * t259 - t165) * pkin(9)) * t364 + t377, -t163 * t200 - t164 * t201 - t239 * t281 - t375 * pkin(2) + (-t95 * t364 - t390 * t368 + (-t163 * t368 - t164 * t364) * qJD(3) + t377) * pkin(9), -t299 * t395 + t394 * t495, -t179 * t495 + t299 * t443 - t393 * t395 + t394 * t496, t299 * t385 - t394 * t458 - t409 * t495, -t179 * t496 + t393 * t443, t179 * t458 + t385 * t393 - t409 * t496, t400, -t125 * t393 - t496 * t177 + t420 * t179 - t237 * t385 + t350 * t443 - t352 * t577 + t500 * t409 - t84 * t458, t125 * t299 + t177 * t495 - t238 * t385 + t350 * t395 + t394 * t420 + t499 * t409 + t458 * t85 + t374, -t179 * t499 - t237 * t395 + t238 * t443 + t24 * t393 + t299 * t554 + t394 * t500 - t495 * t84 + t496 * t85 + t377, t24 * t238 + t554 * t237 - t125 * t350 - g(1) * t497 - g(2) * t498 - g(3) * (-t359 * t507 + t302) + t499 * t85 - t500 * t84 + t420 * t177, t21, t5, t11, t22, t12, t56, t120 * t145 + t149 * t501 + t19 * t518 + t237 * t75 + t383 * t79 - t393 * t4 - t42 * t496 + t550 * t591 + t388, -t120 * t146 + t151 * t501 + t19 * t517 - t237 * t74 + t3 * t393 + t382 * t79 + t43 * t496 - t551 * t591 + t389, t145 * t74 - t146 * t75 + t205 * t43 + t206 * t42 + (-t363 * t43 - t367 * t42) * t286 - t550 * t151 - t551 * t149 - t374 + (-t3 * t363 - t367 * t4 + (t363 * t42 - t367 * t43) * qJD(5)) * t299, t3 * t146 + t4 * t145 + t19 * t237 - g(1) * (-t291 * t418 + t497) - g(2) * (-t289 * t418 + t498) - t559 + t501 * t79 + t551 * t43 + t550 * t42 - (t369 * t418 - t507) * t558, t21, t5, t11, t22, t12, t56, -t1 * t393 + t120 * t127 + t149 * t545 + t17 * t518 + t193 * t75 - t25 * t496 + t383 * t63 + t553 * t591 + t388, -t120 * t134 + t151 * t545 + t17 * t517 - t193 * t74 + t2 * t393 + t30 * t496 + t382 * t63 - t552 * t591 + t389, t127 * t74 - t134 * t75 + t440 * t30 + t439 * t25 - t553 * t151 - t552 * t149 - t374 + (-t1 * t367 - t2 * t363 + (t25 * t363 - t30 * t367) * qJD(5)) * t299, t2 * t134 + t1 * t127 + t17 * t193 - g(1) * (pkin(5) * t520 - t291 * t407 + t497) - g(2) * (pkin(5) * t523 - t289 * t407 + t498) - t559 + t545 * t63 + t552 * t30 + t553 * t25 - (t407 * t369 + (pkin(5) * t363 - t362) * t365) * t558; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t527, -t259 ^ 2 + t260 ^ 2, -t165 + t580, t527, -t166 - t579, t385, -g(1) * t231 - g(2) * t384 + g(3) * t287 - t164 * t335 - t239 * t260 + t442, g(1) * t232 + g(2) * t438 + g(3) * t288 - t239 * t259 - t592, 0, 0, t596, -t530 + t531, -t395 - t581, -t596, -t436 + t443, t385, -t88 * t409 - t177 * t394 + (-t260 * t179 + t385 * t541) * pkin(3) + t387 - t554, -t89 * t409 + t177 * t179 + (-t260 * t394 - t358 * t385) * pkin(3) + t386 - t24 (t358 * t443 + t395 * t541) * pkin(3) + (t85 - t88) * t394 + (t89 - t84) * t179, -g(1) * t421 - g(2) * t376 + g(3) * t597 - t177 * t568 + t24 * t567 - t459 * t554 + t84 * t88 - t85 * t89, t31, t10, t36, t32, t35, -t532, -t149 * t88 - t591 * t49 - t394 * t42 + t348 * t75 + t391 * t363 + (-qJD(5) * t347 * t591 - t19 + t387) * t367, -t151 * t88 + t19 * t363 - t348 * t74 + t391 * t367 + t394 * t43 + t583 * t591 - t465, t149 * t50 + t151 * t49 + (-t179 * t42 - t347 * t75 + t3 + (t151 * t347 - t42) * qJD(5)) * t367 + (-t179 * t43 - t347 * t74 - t4 + (t149 * t347 - t43) * qJD(5)) * t363 - t386, t19 * t348 - t42 * t49 - t79 * t88 - g(1) * (-pkin(4) * t224 + pkin(10) * t225 + t421) - g(2) * (-t220 * pkin(4) + t221 * pkin(10) + t376) - g(3) * (-pkin(4) * t268 + pkin(10) * t269 - t597) - t583 * t43 + (t3 * t367 - t4 * t363 - t42 * t483) * t347, t31, t10, t36, t32, t35, -t532, -t120 * t295 - t149 * t71 - t394 * t25 + t311 * t75 + t542 * t591 + (t179 * t63 + (t63 + t557) * qJD(5)) * t363 + t379 * t367, t63 * t528 - t120 * t296 - t151 * t71 + t17 * t363 + t394 * t30 - t311 * t74 - t586 * t591 + (pkin(5) * t533 + t367 * t63) * qJD(5) - t465, -t149 * t586 - t542 * t151 - t295 * t74 - t296 * t75 - t363 * t575 + t410 * t367 - t386, t2 * t296 - t1 * t295 + t17 * t311 - g(1) * (-t224 * t349 - t225 * t361 + t421) - g(2) * (-t220 * t349 - t221 * t361 + t376) - g(3) * (-t268 * t349 - t269 * t361 - t597) + (pkin(5) * t484 - t71) * t63 + t586 * t30 + t542 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t436 - t443, -t395 + t581, -t530 - t531, t179 * t85 + t394 * t84 + t125 + t577, 0, 0, 0, 0, 0, 0, t34, t37, t9, -t394 * t79 + (t4 + t601) * t367 + (t3 - t602) * t363 + t577, 0, 0, 0, 0, 0, 0, t34, t37, t9, t410 * t363 + t367 * t575 - t394 * t63 + t577; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t536, t76, t45, -t536, t46, t120, -t151 * t79 + t371 + t601, t149 * t79 + t372 + t602, 0, 0, t536, t76, t45, -t536, t46, t120, 0.2e1 * t566 + t549 + t600 + (-t444 - t63) * t151 + t371, -pkin(5) * t572 + t548 + t591 * t29 + (qJD(6) + t63) * t149 + t372, pkin(5) * t74 - t149 * t555, t555 * t30 + (-t63 * t151 + t1 + t576) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75 + t535, -t74 - t539, -t148 - t572, t30 * t149 + t25 * t151 - t379;];
tau_reg  = t18;
