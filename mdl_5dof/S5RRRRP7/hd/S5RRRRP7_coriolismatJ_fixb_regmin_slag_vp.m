% Calculate minimal parameter regressor of coriolis matrix for
% S5RRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x28]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRRRP7_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP7_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:57:47
% EndTime: 2019-12-31 21:58:06
% DurationCPUTime: 8.97s
% Computational Cost: add. (8171->548), mult. (16182->686), div. (0->0), fcn. (16696->6), ass. (0->433)
t348 = cos(qJ(2));
t577 = sin(qJ(3));
t334 = t577 * t348;
t346 = sin(qJ(2));
t578 = cos(qJ(3));
t335 = t578 * t346;
t306 = -t335 - t334;
t575 = pkin(4) * t306;
t345 = sin(qJ(4));
t595 = -pkin(7) - pkin(6);
t316 = t595 * t346;
t317 = t595 * t348;
t407 = -t578 * t316 - t577 * t317;
t616 = t407 * t345;
t623 = t575 / 0.2e1 - t616 / 0.2e1;
t347 = cos(qJ(4));
t529 = t306 * t347;
t530 = t306 * t345;
t410 = -pkin(4) * t530 + qJ(5) * t529;
t127 = t407 + t410;
t547 = t127 * t345;
t121 = t547 / 0.2e1;
t313 = pkin(4) * t345 - qJ(5) * t347;
t526 = t313 * t345;
t439 = -t526 / 0.2e1;
t596 = pkin(4) / 0.2e1;
t622 = t121 + (t439 + t596) * t306 + t623;
t472 = qJD(2) + qJD(3);
t615 = t407 * t347;
t621 = -t615 / 0.2e1;
t193 = t616 / 0.2e1;
t295 = t306 * qJ(5);
t572 = t306 * pkin(3);
t304 = t577 * t346 - t578 * t348;
t573 = t304 * pkin(8);
t209 = -t572 + t573;
t576 = pkin(2) * t346;
t181 = t209 + t576;
t166 = t345 * t181;
t511 = t166 / 0.2e1 + t621;
t619 = t511 - t295;
t522 = t345 * qJ(5);
t570 = t347 * pkin(4);
t403 = -t522 - t570;
t312 = -pkin(3) + t403;
t470 = t578 * pkin(2);
t294 = -t470 + t312;
t275 = t294 * t345;
t302 = t312 * t345;
t503 = t275 / 0.2e1 + t302 / 0.2e1;
t525 = t313 * t347;
t618 = t525 - t503;
t617 = t575 - t616;
t601 = t472 * t306;
t405 = t304 * t601;
t174 = t345 * t304;
t177 = t347 * t304;
t408 = t577 * t316 - t578 * t317;
t360 = -pkin(4) * t174 + qJ(5) * t177 + t408;
t605 = t360 * t345;
t339 = -pkin(2) * t348 - pkin(1);
t571 = t306 * pkin(8);
t411 = t304 * pkin(3) + t571;
t370 = t339 + t411;
t603 = t408 * t345;
t111 = -t347 * t370 + t603;
t574 = t304 * pkin(4);
t91 = t111 - t574;
t614 = (t91 - t605) * t306;
t613 = (t111 - t603) * t306;
t602 = t408 * t347;
t112 = t345 * t370 + t602;
t612 = (t112 - t602) * t306;
t611 = t472 * t407;
t604 = t360 * t347;
t533 = t304 * qJ(5);
t90 = t112 + t533;
t610 = t127 * t177 + (-t90 + t604) * t306;
t609 = t472 * t408;
t608 = t127 * t360;
t588 = t294 / 0.2e1;
t388 = t470 / 0.2e1 + t588;
t586 = -t312 / 0.2e1;
t607 = t304 * (t586 + t388);
t469 = t577 * pkin(2);
t337 = t469 + pkin(8);
t461 = t337 / 0.2e1 - pkin(8) / 0.2e1;
t606 = t306 * t461;
t550 = t360 * t312;
t206 = t472 * t304;
t303 = t306 ^ 2;
t598 = t304 ^ 2;
t471 = t303 - t598;
t343 = t345 ^ 2;
t344 = t347 ^ 2;
t173 = (-t343 / 0.2e1 + t344 / 0.2e1) * t306;
t519 = t345 * t347;
t452 = qJD(1) * t519;
t255 = t303 * t452;
t599 = t173 * t472 + t255;
t326 = t344 - t343;
t416 = t306 * t452;
t155 = t326 * t472 + 0.2e1 * t416;
t597 = pkin(3) / 0.2e1;
t594 = -qJ(5) / 0.2e1;
t517 = t347 * t209;
t101 = -t517 + t617;
t593 = t101 / 0.2e1;
t172 = t403 * t306;
t592 = -t172 / 0.2e1;
t591 = -t181 / 0.2e1;
t590 = -t209 / 0.2e1;
t589 = t615 / 0.2e1;
t587 = -t306 / 0.2e1;
t585 = t312 / 0.2e1;
t584 = t313 / 0.2e1;
t338 = -t470 - pkin(3);
t583 = -t338 / 0.2e1;
t582 = -t345 / 0.2e1;
t581 = t345 / 0.2e1;
t580 = -t347 / 0.2e1;
t579 = t347 / 0.2e1;
t569 = t304 * t90;
t504 = -t615 + t166;
t94 = -t295 + t504;
t518 = t347 * t181;
t95 = -t518 + t617;
t7 = t90 * t94 + t91 * t95 + t608;
t568 = t7 * qJD(1);
t185 = t345 * t209;
t505 = t185 - t615;
t100 = -t295 + t505;
t8 = t100 * t90 + t101 * t91 + t608;
t567 = t8 * qJD(1);
t9 = -t111 * t90 + t112 * t91 + t127 * t172;
t566 = t9 * qJD(1);
t565 = t94 * t347;
t564 = t95 * t345;
t18 = t111 * t580 + t90 * t582 + t112 * t581 + t91 * t579 + (t522 / 0.2e1 + t570 / 0.2e1) * t304;
t563 = t18 * qJD(4);
t459 = -t90 / 0.2e1 + t112 / 0.2e1;
t460 = -t111 / 0.2e1 + t91 / 0.2e1;
t17 = (t574 / 0.2e1 - t460) * t347 + (t533 / 0.2e1 - t459) * t345;
t562 = -t17 * qJD(4) + t177 * qJD(5);
t551 = t112 * t304;
t41 = -t551 + (-t127 * t347 - t172 * t345) * t306;
t561 = qJD(1) * t41;
t552 = t111 * t304;
t42 = -t552 + (t172 * t347 - t547) * t306;
t560 = qJD(1) * t42;
t51 = t127 * t529 + t569;
t559 = qJD(1) * t51;
t61 = t407 * t530 + t552;
t558 = qJD(1) * t61;
t62 = -t407 * t529 - t551;
t557 = qJD(1) * t62;
t68 = t91 * t177;
t10 = t68 + t95 * t529 + (-t306 * t94 - t569) * t345;
t556 = t10 * qJD(1);
t555 = t100 * t347;
t554 = t101 * t345;
t11 = t68 + t101 * t529 + (-t100 * t306 - t569) * t345;
t553 = t11 * qJD(1);
t13 = t112 * t529 + (-t347 * t90 + (t111 - t91) * t345) * t306;
t544 = t13 * qJD(1);
t14 = t17 * qJD(1);
t19 = t94 * t304 + t610;
t543 = t19 * qJD(1);
t20 = t614 + (-t95 - t547) * t304;
t542 = t20 * qJD(1);
t21 = t100 * t304 + t610;
t541 = t21 * qJD(1);
t22 = t614 + (-t101 - t547) * t304;
t540 = t22 * qJD(1);
t536 = t294 * t304;
t535 = t294 * t313;
t534 = t294 * t347;
t532 = t304 * t338;
t531 = t306 * t337;
t528 = t312 * t347;
t527 = t313 * t312;
t524 = t337 * t304;
t523 = t338 * t306;
t37 = t518 * t304 + t613;
t515 = t37 * qJD(1);
t38 = t612 + (-t504 - t615) * t304;
t514 = t38 * qJD(1);
t39 = t209 * t177 + t613;
t513 = t39 * qJD(1);
t40 = t612 + (-t505 - t615) * t304;
t512 = t40 * qJD(1);
t509 = t518 / 0.2e1 + t193;
t508 = t185 / 0.2e1 + t621;
t507 = -t185 / 0.2e1 + t589;
t180 = t517 / 0.2e1;
t506 = t180 + t193;
t501 = t343 * t306 / 0.2e1 + t344 * t587;
t488 = qJD(3) * t347;
t491 = qJD(2) * t347;
t500 = (t488 + t491) * t345;
t428 = t577 * qJD(3);
t419 = pkin(2) * t428;
t325 = t345 * t419;
t342 = t343 * qJD(5);
t499 = t342 - t325;
t164 = t304 * t576 - t306 * t339;
t498 = qJD(1) * t164;
t165 = -t304 * t339 - t306 * t576;
t497 = qJD(1) * t165;
t496 = qJD(1) * t304;
t495 = qJD(1) * t306;
t494 = qJD(1) * t339;
t493 = qJD(1) * t348;
t492 = qJD(2) * t345;
t490 = qJD(3) * t339;
t489 = qJD(3) * t345;
t487 = qJD(4) * t111;
t486 = qJD(4) * t304;
t485 = qJD(4) * t345;
t484 = qJD(4) * t347;
t483 = qJD(5) * t345;
t142 = t471 * t345;
t482 = t142 * qJD(1);
t143 = t471 * t347;
t481 = t143 * qJD(1);
t480 = t471 * qJD(1);
t479 = t174 * qJD(1);
t171 = t177 * qJD(1);
t184 = t326 * t303;
t478 = t184 * qJD(1);
t297 = t335 / 0.2e1 + t334 / 0.2e1;
t477 = t297 * qJD(1);
t293 = t304 * qJD(5);
t327 = -t346 ^ 2 + t348 ^ 2;
t476 = t327 * qJD(1);
t475 = t346 * qJD(2);
t474 = t347 * qJD(5);
t473 = t348 * qJD(2);
t468 = pkin(1) * t346 * qJD(1);
t467 = pkin(1) * t493;
t466 = pkin(8) * t485;
t465 = pkin(8) * t484;
t464 = -t575 / 0.2e1;
t463 = pkin(8) * t581;
t462 = t597 + t583;
t122 = -t547 / 0.2e1;
t435 = t177 / 0.2e1;
t441 = t529 / 0.2e1;
t458 = t294 * t441 + t337 * t435 + t122;
t457 = pkin(8) * t435 + t312 * t441 + t122;
t456 = t345 * t578;
t455 = t577 * t306;
t454 = t304 * t494;
t453 = t306 * t494;
t451 = t306 * t486;
t450 = t337 * t485;
t449 = t337 * t484;
t217 = t304 * t495;
t333 = t345 * t484;
t448 = t306 * t483;
t330 = t345 * t474;
t447 = t346 * t493;
t446 = t127 * t584;
t443 = t531 / 0.2e1;
t442 = -t529 / 0.2e1;
t440 = -t528 / 0.2e1;
t438 = -t523 / 0.2e1;
t436 = -t177 / 0.2e1;
t432 = t585 + t588;
t431 = t578 * qJD(2);
t430 = t578 * qJD(3);
t429 = t577 * qJD(2);
t367 = (-t469 / 0.2e1 + t461) * t306;
t28 = (t367 - t607) * t345;
t418 = pkin(2) * t429;
t324 = t347 * t418;
t427 = t28 * qJD(1) - t324;
t404 = pkin(2) * t347 * t455;
t251 = t404 / 0.2e1;
t30 = t251 + (-t606 + t607) * t347;
t323 = t345 * t418;
t426 = -qJD(1) * t30 + t323;
t250 = -t404 / 0.2e1;
t417 = -t470 / 0.2e1;
t389 = t417 + t583;
t366 = (-pkin(3) / 0.2e1 + t389) * t304;
t50 = t250 + (t366 + t606) * t347;
t425 = -qJD(1) * t50 - t323;
t47 = (t367 + t366) * t345;
t424 = t47 * qJD(1) - t324;
t422 = t472 * t347;
t421 = t306 * t463;
t420 = qJD(4) + t496;
t415 = -t456 / 0.2e1;
t414 = t456 / 0.2e1;
t413 = t578 * t579;
t412 = t592 - t573 / 0.2e1;
t409 = t592 - t524 / 0.2e1;
t139 = t297 + t501;
t115 = qJD(1) * t139 + t500;
t134 = qJD(1) * t173 - t500;
t406 = t345 * t422;
t402 = -t304 * t312 + t571;
t401 = t564 + t565;
t400 = t347 * t419;
t371 = (t343 + t344) * t578;
t144 = (t577 * t294 + t371 * t337) * pkin(2);
t380 = t555 / 0.2e1 + t554 / 0.2e1;
t349 = t380 * t337 + (t90 * t413 + t127 * t577 / 0.2e1 + t91 * t414) * pkin(2) + t360 * t588;
t383 = -t565 / 0.2e1 - t564 / 0.2e1;
t2 = -t550 / 0.2e1 + t383 * pkin(8) + t349;
t399 = t2 * qJD(1) + t144 * qJD(2);
t398 = t554 + t555;
t397 = -t531 + t536;
t396 = t531 - t532;
t191 = t526 + t534;
t381 = (t306 * t584 - t127 / 0.2e1) * t347;
t352 = t381 + (t294 * t587 + t409) * t345;
t34 = t352 - t619;
t395 = -qJD(1) * t34 + qJD(2) * t191;
t192 = -t275 + t525;
t189 = t294 * t442;
t25 = t189 + (t591 + t409) * t347 + t622;
t394 = -qJD(1) * t25 + qJD(2) * t192;
t24 = (t100 / 0.2e1 - t94 / 0.2e1) * t347 + (t593 - t95 / 0.2e1) * t345;
t296 = t371 * pkin(2);
t393 = -qJD(1) * t24 - qJD(2) * t296;
t392 = t306 * t420;
t208 = t303 * t344 + t598;
t391 = -qJD(1) * t208 - t486;
t390 = qJD(4) * t313 - t483;
t387 = t94 * t594 + t95 * t596;
t386 = t306 * t406;
t385 = t601 * t519;
t384 = pkin(4) * t593 + t100 * t594;
t358 = t345 * t459 + t347 * t460;
t350 = t172 * t588 + t358 * t337 + t446;
t4 = t350 + t387;
t382 = t4 * qJD(1) + qJD(2) * t535;
t44 = t464 + t458 + t509;
t379 = -qJD(1) * t44 + t294 * t492;
t234 = t337 * t436;
t54 = t234 + (t438 + t591) * t347;
t378 = -qJD(1) * t54 - t338 * t492;
t355 = (t524 / 0.2e1 + t523 / 0.2e1) * t345 + t589;
t52 = t355 + t511;
t377 = -qJD(1) * t52 - t338 * t491;
t376 = t347 * t392;
t146 = -qJD(4) * t297 + t217;
t375 = t330 - t400;
t373 = t306 * t439 + t121 + 0.2e1 * t464;
t372 = -0.2e1 * t386;
t320 = pkin(2) * t415;
t130 = t320 + t618;
t219 = -t302 + t525;
t215 = t306 * t440;
t31 = t215 + (t590 + t412) * t347 + t622;
t369 = -qJD(1) * t31 + qJD(2) * t130 + qJD(3) * t219;
t321 = pkin(2) * t413;
t131 = t347 * t432 + t321 + t526;
t218 = t526 + t528;
t353 = t381 + (t306 * t586 + t412) * t345;
t36 = t295 + t353 + t507;
t368 = -qJD(1) * t36 + qJD(2) * t131 + qJD(3) * t218;
t242 = t345 * t462 + t320;
t266 = pkin(8) * t436;
t58 = t266 + (t572 / 0.2e1 + t590) * t347;
t364 = pkin(3) * t489 - qJD(1) * t58 + qJD(2) * t242;
t322 = t347 * t417;
t243 = t347 * t462 + t322;
t357 = (t573 / 0.2e1 - t572 / 0.2e1) * t345 + t589;
t56 = t357 + t508;
t363 = pkin(3) * t488 - qJD(1) * t56 + qJD(2) * t243;
t356 = (pkin(4) * t415 + qJ(5) * t413) * pkin(2);
t117 = -t313 * t432 + t356;
t351 = t358 * pkin(8) + t172 * t585 + t446;
t6 = t351 + t384;
t362 = t6 * qJD(1) - t117 * qJD(2) + qJD(3) * t527;
t319 = pkin(2) * t414;
t151 = t319 + t503;
t46 = t464 + t457 + t506;
t361 = -qJD(1) * t46 + qJD(2) * t151 + t312 * t489;
t359 = qJD(4) * t403 + t474;
t354 = (-t578 * t304 / 0.2e1 - t455 / 0.2e1) * pkin(2) + t443;
t315 = t326 * qJD(4);
t284 = t296 * qJD(3);
t256 = t306 * t330;
t254 = t420 * qJ(5);
t245 = pkin(3) * t580 + t338 * t579 + t322;
t244 = pkin(3) * t582 + t338 * t581 + t320;
t205 = t343 * t472 - t416;
t182 = t472 * t297;
t167 = t173 * qJD(4);
t152 = t319 - t503;
t150 = t171 + t484;
t149 = -t479 - t485;
t140 = -t297 + t501;
t138 = t140 * qJD(5);
t137 = t139 * qJD(5);
t133 = t320 - t618;
t132 = -t534 / 0.2e1 - t526 + t440 + t321;
t129 = 0.2e1 * t345 * t376;
t118 = t527 / 0.2e1 + t535 / 0.2e1 + t356;
t113 = -t217 * t344 - t167;
t110 = t112 * qJD(4);
t86 = qJD(4) * t177 - t481;
t85 = -qJD(4) * t174 + t482;
t67 = -t177 * t472 - t217 * t345;
t66 = -t167 + (t344 * t495 - t406) * t304;
t65 = -t345 * t601 + t481;
t64 = -t306 * t422 - t482;
t63 = t345 * t392;
t60 = 0.2e1 * (qJD(4) - t496) * t306 * t519 - t326 * t206;
t59 = pkin(3) * t441 + t180 + 0.2e1 * t193 + t266;
t57 = t357 + t507;
t55 = t347 * t438 + t193 + t234 + t509;
t53 = t355 - t511;
t49 = t250 + t603 + pkin(8) * t441 + pkin(3) * t435 + (t304 * t389 + t443) * t347;
t48 = -t602 + t421 + t174 * t597 + (-t532 / 0.2e1 + t354) * t345;
t45 = -t517 / 0.2e1 + t457 + t623;
t43 = -t518 / 0.2e1 + t458 + t623;
t35 = -t295 + t353 + t508;
t33 = t352 + t619;
t32 = t347 * t412 + t215 + t373 + t506;
t29 = t251 - t605 + pkin(8) * t442 + t312 * t435 + (-t531 / 0.2e1 + t388 * t304) * t347;
t27 = -t604 + t174 * t586 + t421 + (-t536 / 0.2e1 + t354) * t345;
t26 = t347 * t409 + t189 + t373 + t509;
t23 = t380 - t383;
t12 = t359 - t14;
t5 = t351 - t384;
t3 = t350 - t387;
t1 = pkin(8) * t565 / 0.2e1 + t550 / 0.2e1 + t95 * t463 + t349;
t15 = [0, 0, 0, t346 * t473, t327 * qJD(2), 0, 0, 0, -pkin(1) * t475, -pkin(1) * t473, t405, -t472 * t471, 0, 0, 0, qJD(2) * t164 - t306 * t490, qJD(2) * t165 - t304 * t490, -t303 * t333 + t344 * t405, -t184 * qJD(4) + t304 * t372, t143 * t472 + t345 * t451, -t142 * t472 + t347 * t451, -t405, qJD(2) * t37 + qJD(3) * t39 + qJD(4) * t62, qJD(2) * t38 + qJD(3) * t40 + qJD(4) * t61, qJD(2) * t20 + qJD(3) * t22 + qJD(4) * t41 - t303 * t330, -qJD(2) * t10 - qJD(3) * t11 - qJD(4) * t13 + t304 * t448, qJD(2) * t19 + qJD(3) * t21 + qJD(4) * t42 + qJD(5) * t208, qJD(2) * t7 + qJD(3) * t8 + qJD(4) * t9 + qJD(5) * t51; 0, 0, 0, t447, t476, t473, -t475, 0, -pkin(6) * t473 - t468, pkin(6) * t475 - t467, t217, -t480, -t206, t601, 0, t498 - t609, t497 + t611, t66, t60, t65, t64, -t146, t515 + (t345 * t396 - t602) * qJD(2) + t48 * qJD(3) + t55 * qJD(4), t514 + (t347 * t396 + t603) * qJD(2) + t49 * qJD(3) + t53 * qJD(4), t542 + (-t345 * t397 - t604) * qJD(2) + t27 * qJD(3) + t26 * qJD(4) + t138, qJD(2) * t401 + t23 * qJD(3) - t556 + t563, t543 + (t347 * t397 - t605) * qJD(2) + t29 * qJD(3) + t33 * qJD(4) - t256, t568 + (t294 * t360 + t337 * t401) * qJD(2) + t1 * qJD(3) + t3 * qJD(4) + t43 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t217, -t480, -t206, t601, 0, -t453 - t609, -t454 + t611, t66, t60, t65, t64, -t146, t513 + t48 * qJD(2) + (t345 * t411 - t602) * qJD(3) + t59 * qJD(4), t512 + t49 * qJD(2) + (t347 * t411 + t603) * qJD(3) + t57 * qJD(4), t540 + t27 * qJD(2) + (t345 * t402 - t604) * qJD(3) + t32 * qJD(4) + t138, t23 * qJD(2) + qJD(3) * t398 - t553 + t563, t541 + t29 * qJD(2) + (-t347 * t402 - t605) * qJD(3) + t35 * qJD(4) - t256, t567 + t1 * qJD(2) + (pkin(8) * t398 + t550) * qJD(3) + t5 * qJD(4) + t45 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t599, 0.2e1 * t385 - t478, t63, t376, t182, qJD(2) * t55 + qJD(3) * t59 - t110 + t557, qJD(2) * t53 + qJD(3) * t57 + t487 + t558, qJD(2) * t26 + qJD(3) * t32 - t110 + t561, qJD(4) * t410 + t18 * t472 + t448 - t544, qJD(2) * t33 + qJD(3) * t35 + t293 - t487 + t560, t566 + t3 * qJD(2) + t5 * qJD(3) + (-pkin(4) * t112 - qJ(5) * t111) * qJD(4) + t90 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140 * t472 - t255, t63, -t386 - t391, qJD(2) * t43 + qJD(3) * t45 + qJD(4) * t90 + t559; 0, 0, 0, -t447, -t476, 0, 0, 0, t468, t467, -t217, t480, 0, 0, 0, -t498, -t497, t113, t129, t86, t85, t146, qJD(3) * t47 + qJD(4) * t54 - t515, qJD(3) * t50 + qJD(4) * t52 - t514, qJD(3) * t28 + qJD(4) * t25 + t137 - t542, qJD(3) * t24 + t556 + t562, qJD(3) * t30 + qJD(4) * t34 - t256 - t543, qJD(3) * t2 + qJD(4) * t4 + qJD(5) * t44 - t568; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t419, -pkin(2) * t430, t333, t315, 0, 0, 0, t338 * t485 - t400, t338 * t484 + t325, -t192 * qJD(4) + t375, t284, -qJD(4) * t191 + t499, t144 * qJD(3) + t294 * t390; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t429 - t428) * pkin(2), (-t431 - t430) * pkin(2), t333, t315, 0, 0, 0, t244 * qJD(4) - t400 + t424, qJD(4) * t245 + t325 - t425, t133 * qJD(4) + t375 + t427, t284 - t393, qJD(4) * t132 - t426 + t499, (pkin(8) * t371 + t312 * t577) * pkin(2) * qJD(3) + t118 * qJD(4) + t152 * qJD(5) + t399; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134, t155, t150, t149, -t477, qJD(3) * t244 - t378 - t449, qJD(3) * t245 - t377 + t450, qJD(3) * t133 - t394 - t449, t12, qJD(3) * t132 - t395 - t450, t118 * qJD(3) + t337 * t359 + t382; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, t150, t205, qJD(3) * t152 - t379 + t449; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t217, t480, 0, 0, 0, t453, t454, t113, t129, t86, t85, t146, -qJD(2) * t47 + qJD(4) * t58 - t513, -qJD(2) * t50 + qJD(4) * t56 - t512, -qJD(2) * t28 + qJD(4) * t31 + t137 - t540, -qJD(2) * t24 + t553 + t562, -qJD(2) * t30 + qJD(4) * t36 - t256 - t541, -qJD(2) * t2 + qJD(4) * t6 + qJD(5) * t46 - t567; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t418, pkin(2) * t431, t333, t315, 0, 0, 0, -qJD(4) * t242 - t424, -qJD(4) * t243 + t425, -qJD(4) * t130 + t330 - t427, t393, -qJD(4) * t131 + t342 + t426, -qJD(4) * t117 - qJD(5) * t151 - t399; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t333, t315, 0, 0, 0, -pkin(3) * t485, -pkin(3) * t484, -qJD(4) * t219 + t330, 0, -qJD(4) * t218 + t342, t390 * t312; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134, t155, t150, t149, -t477, -t364 - t465, -t363 + t466, -t369 - t465, t12, -t368 - t466, pkin(8) * t359 + t362; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, t150, t205, -t361 + t465; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t599, t372 + t478, t67, t174 * t472 - t217 * t347, t182, -qJD(2) * t54 - qJD(3) * t58 - t557, -qJD(2) * t52 - qJD(3) * t56 - t558, -qJD(2) * t25 - qJD(3) * t31 - t561, t17 * t472 + t544, -qJD(2) * t34 - qJD(3) * t36 + t293 - t560, qJ(5) * t293 - qJD(2) * t4 - qJD(3) * t6 - t566; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, -t155, -t171, t479, t477, qJD(3) * t242 + t378, qJD(3) * t243 + t377, qJD(3) * t130 + t394, t14, qJD(3) * t131 + t395, qJD(3) * t117 - t382; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, -t155, -t171, t479, t477, t364, t363, t369, t14, t368, -t362; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), qJ(5) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t420, t254; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t139 * t472 + t255, t67, t385 + t391, -qJ(5) * t486 - qJD(2) * t44 - qJD(3) * t46 - t559; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t115, -t171, -t205, qJD(3) * t151 + t379; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t115, -t171, -t205, t361; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t420, -t254; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t15;
