% Calculate inertial parameters regressor of coriolis matrix for
% S6RPRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPRRPP2_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP2_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:33:28
% EndTime: 2019-03-09 04:33:44
% DurationCPUTime: 10.39s
% Computational Cost: add. (5594->636), mult. (10954->730), div. (0->0), fcn. (9228->6), ass. (0->428)
t356 = sin(qJ(4));
t345 = t356 * qJ(5);
t358 = cos(qJ(4));
t397 = t358 * pkin(4) + t345;
t460 = t358 * qJD(5);
t593 = -qJD(4) * t397 + t460;
t352 = t356 ^ 2;
t354 = t358 ^ 2;
t309 = t354 + t352;
t579 = pkin(4) + pkin(5);
t592 = -t358 * t579 - t345;
t591 = 0.2e1 * t356;
t590 = pkin(8) - qJ(6);
t329 = -cos(pkin(9)) * pkin(1) - pkin(2);
t359 = cos(qJ(3));
t357 = sin(qJ(3));
t557 = t357 * pkin(8);
t402 = -t359 * pkin(3) - t557;
t364 = t402 + t329;
t220 = t358 * t364;
t328 = sin(pkin(9)) * pkin(1) + pkin(7);
t492 = t359 * t328;
t449 = t356 * t492;
t145 = -t220 + t449;
t501 = t357 * t358;
t448 = qJ(6) * t501;
t110 = -t145 + t448;
t348 = t359 * pkin(4);
t80 = -t448 - t220 + t348 + (t328 * t356 + pkin(5)) * t359;
t589 = t110 + t80;
t250 = pkin(3) - t592;
t504 = t357 * t250;
t427 = -t504 / 0.2e1;
t494 = t359 * qJ(6);
t588 = t427 + t494 / 0.2e1;
t463 = t357 * qJD(3);
t318 = t356 * t463;
t355 = t359 ^ 2;
t353 = t357 ^ 2;
t497 = t358 * t353;
t272 = t358 * t355 - t497;
t468 = t272 * qJD(1);
t587 = t468 - t318;
t457 = t359 * qJD(4);
t442 = t358 * t457;
t197 = t468 - t442;
t461 = t358 * qJD(3);
t432 = t357 * t461;
t507 = t356 * t353;
t270 = t356 * t355 - t507;
t469 = t270 * qJD(1);
t195 = t469 + t432;
t443 = t356 * t457;
t586 = t469 - t443;
t555 = t359 * pkin(8);
t559 = t357 * pkin(3);
t293 = -t555 + t559;
t273 = t356 * t293;
t187 = -t328 * t501 + t273;
t527 = t187 * t357;
t502 = t357 * t328;
t275 = t356 * t502;
t498 = t358 * t293;
t186 = t275 + t498;
t528 = t186 * t357;
t276 = t358 * t492;
t146 = t356 * t364 + t276;
t533 = t146 * t359;
t534 = t145 * t359;
t28 = (t527 / 0.2e1 + t533 / 0.2e1) * t358 + (-t528 / 0.2e1 + t534 / 0.2e1) * t356 + (t353 / 0.2e1 - t355 / 0.2e1) * t328;
t500 = t357 * t359;
t193 = (-0.1e1 + t309) * t500;
t473 = t193 * qJD(2);
t585 = -t28 * qJD(1) - t473;
t458 = t359 * qJD(3);
t459 = t359 * qJD(1);
t324 = -qJD(4) + t459;
t310 = t354 - t352;
t104 = t501 * t591 * (qJD(4) + t459) - t310 * t458;
t382 = t324 * t357;
t218 = t358 * t382;
t181 = t218 * t591;
t582 = pkin(4) / 0.2e1;
t581 = pkin(5) / 0.2e1;
t506 = t356 * t357;
t111 = qJ(6) * t506 + t146;
t495 = t359 * qJ(5);
t90 = t111 - t495;
t580 = -t90 / 0.2e1;
t499 = t358 * qJ(5);
t323 = t357 * t499;
t483 = -pkin(4) * t506 + t323;
t143 = (-t356 * pkin(5) - t328) * t357 + t483;
t578 = -t143 / 0.2e1;
t189 = -t483 + t502;
t577 = t189 / 0.2e1;
t202 = t592 * t357;
t576 = -t202 / 0.2e1;
t575 = t202 / 0.2e1;
t244 = t397 * t357;
t574 = -t244 / 0.2e1;
t573 = -t275 / 0.2e1;
t571 = -t293 / 0.2e1;
t347 = t357 * pkin(4);
t570 = -t347 / 0.2e1;
t569 = t352 / 0.2e1;
t568 = -t354 / 0.2e1;
t567 = -t356 / 0.2e1;
t566 = t356 / 0.2e1;
t565 = -t357 / 0.2e1;
t564 = t357 / 0.2e1;
t563 = -t358 / 0.2e1;
t562 = t358 / 0.2e1;
t561 = t359 / 0.2e1;
t560 = t579 / 0.2e1;
t558 = t357 * pkin(5);
t346 = t357 * qJ(5);
t150 = t187 + t346;
t113 = t356 * t494 + t150;
t258 = -t356 * t579 + t499;
t144 = (-t328 + t258) * t359;
t551 = t90 * t358;
t95 = -t558 - t275 - t347 + (-t293 - t494) * t358;
t7 = (t551 / 0.2e1 + t80 * t566 + t144 / 0.2e1) * t359 + (t113 * t562 + t566 * t95 + t578) * t357;
t554 = t7 * qJD(1);
t553 = t80 * t358;
t552 = t90 * t356;
t550 = t90 * t359;
t134 = t146 - t495;
t135 = t145 + t348;
t415 = -t145 / 0.2e1 + t135 / 0.2e1;
t362 = t415 * t358 + (-t134 / 0.2e1 + t146 / 0.2e1) * t356;
t18 = t357 * t362 + t359 * t574;
t549 = qJD(1) * t18;
t44 = t110 * t359 + t143 * t506 - t202 * t501;
t548 = qJD(1) * t44;
t114 = t143 * t501;
t45 = t111 * t359 - t202 * t506 - t114;
t547 = qJD(1) * t45;
t46 = (t552 - t553) * t357;
t546 = qJD(1) * t46;
t522 = t244 * t356;
t525 = t189 * t358;
t47 = t533 + (t522 + t525) * t357;
t545 = qJD(1) * t47;
t48 = -t189 * t506 + t244 * t501 - t534;
t544 = qJD(1) * t48;
t51 = -t114 + t550;
t543 = qJD(1) * t51;
t537 = t134 * t359;
t60 = t189 * t501 + t537;
t542 = qJD(1) * t60;
t77 = -t328 * t507 - t534;
t541 = qJD(1) * t77;
t78 = -t328 * t497 - t533;
t540 = qJD(1) * t78;
t421 = t495 / 0.2e1;
t297 = t356 * t421;
t446 = t90 / 0.2e1 - t111 / 0.2e1;
t405 = t446 * t356;
t416 = t359 * t560;
t447 = t110 / 0.2e1 + t80 / 0.2e1;
t10 = t297 - t405 + (t416 + t447) * t358;
t539 = t10 * qJD(1);
t13 = t202 * t561 + (t358 * t447 - t405) * t357;
t538 = t13 * qJD(1);
t290 = pkin(4) * t356 - t499;
t190 = (t290 + t328) * t359;
t153 = -t186 - t347;
t530 = t153 * t356;
t531 = t150 * t358;
t14 = (t134 * t562 - t190 / 0.2e1 + t135 * t566) * t359 + (t531 / 0.2e1 + t577 + t530 / 0.2e1) * t357;
t536 = t14 * qJD(1);
t535 = t143 * t356;
t15 = (-t95 * t357 - t80 * t359) * t358 + (t113 * t357 + t550) * t356;
t532 = t15 * qJD(1);
t16 = ((-t111 + t90) * t358 + t589 * t356) * t357;
t529 = t16 * qJD(1);
t526 = t189 * t356;
t20 = -t80 * t357 + t95 * t359 + (-t143 * t359 - t144 * t357) * t356;
t524 = t20 * qJD(1);
t422 = -t495 / 0.2e1;
t298 = t356 * t422;
t373 = t134 * t566 + t146 * t567;
t23 = t298 + (-t348 / 0.2e1 - t415) * t358 + t373;
t523 = t23 * qJD(1);
t521 = t250 * t358;
t520 = t258 * t356;
t496 = t358 * t359;
t27 = t113 * t359 - t143 * t496 - t144 * t501 - t90 * t357;
t519 = t27 * qJD(1);
t281 = -pkin(3) - t397;
t517 = t281 * t356;
t289 = t590 * t356;
t516 = t289 * t359;
t29 = ((t134 - t146) * t358 + (t135 - t145) * t356) * t357;
t515 = t29 * qJD(1);
t291 = t590 * t358;
t514 = t291 * t357;
t513 = t291 * t359;
t32 = -t153 * t501 - t135 * t496 + (t150 * t357 + t537) * t356;
t512 = t32 * qJD(1);
t33 = (t528 - t534) * t358 + (t527 + t533) * t356;
t511 = t33 * qJD(1);
t387 = t189 * t359 + t190 * t357;
t34 = -t134 * t357 + t150 * t359 + t358 * t387;
t510 = t34 * qJD(1);
t35 = -t135 * t357 + t153 * t359 + t356 * t387;
t509 = t35 * qJD(1);
t508 = t352 * t353;
t505 = t356 * t359;
t503 = t357 * t281;
t493 = t359 * t290;
t491 = t579 * t357;
t414 = t573 + t570;
t363 = -t535 / 0.2e1 + t513 / 0.2e1 + t414;
t420 = -t494 / 0.2e1;
t371 = t427 + t420;
t452 = -t558 / 0.2e1;
t42 = t452 + (t571 + t371) * t358 + t363;
t489 = t42 * qJD(1);
t52 = t145 * t357 + (t186 - 0.2e1 * t275) * t359;
t488 = t52 * qJD(1);
t53 = t187 * t359 + (-t146 + 0.2e1 * t276) * t357;
t487 = t53 * qJD(1);
t477 = qJD(5) * t359;
t486 = t356 * t477 + t473;
t257 = t498 / 0.2e1;
t485 = t257 + t275 / 0.2e1;
t296 = t358 * t421;
t429 = -t505 / 0.2e1;
t484 = pkin(4) * t429 + t296;
t482 = t309 * t555;
t481 = qJD(4) * t145;
t480 = qJD(4) * t291;
t479 = qJD(4) * t356;
t343 = qJD(4) * t358;
t478 = qJD(5) * t356;
t476 = qJD(6) * t356;
t475 = qJD(6) * t358;
t474 = qJD(6) * t359;
t413 = t569 + t568;
t254 = t413 * t357;
t472 = t254 * qJD(4);
t336 = t354 * t353;
t267 = t336 + t508;
t471 = t267 * qJD(1);
t247 = t309 * t458;
t467 = t309 * qJD(3);
t292 = t310 * qJD(4);
t311 = t355 - t353;
t466 = t311 * qJD(1);
t465 = t356 * qJD(3);
t464 = t357 * qJD(1);
t462 = t357 * qJD(4);
t456 = t582 + t581;
t455 = pkin(8) * t479;
t454 = pkin(8) * t343;
t453 = pkin(3) * t566;
t451 = t555 / 0.2e1;
t450 = t581 + t560;
t255 = t273 / 0.2e1;
t425 = -t501 / 0.2e1;
t445 = t328 * t425 + t255 + t346;
t444 = qJ(5) * t477;
t441 = t359 * t460;
t440 = t357 * t475;
t439 = t329 * t464;
t438 = t329 * t459;
t317 = t356 * t343;
t315 = t356 * t461;
t314 = t356 * t460;
t437 = t356 * t459;
t436 = t356 * t474;
t435 = t357 * t458;
t434 = t357 * t478;
t433 = t357 * t459;
t431 = t358 * t464;
t319 = t357 * t460;
t430 = t358 * t459;
t428 = t505 / 0.2e1;
t426 = t503 / 0.2e1;
t424 = t501 / 0.2e1;
t423 = -t496 / 0.2e1;
t418 = t356 * t560;
t417 = t491 / 0.2e1;
t223 = (-0.1e1 / 0.2e1 + t413) * t357;
t163 = t223 * qJD(1) - t315;
t191 = t254 * qJD(1) - t315;
t284 = t356 * qJD(1) * t497;
t171 = qJD(3) * t254 + t284;
t410 = t359 * t434;
t409 = t353 * t317;
t408 = t356 * t431;
t407 = t356 * t432;
t231 = t289 * t425;
t406 = t582 + t450;
t403 = t514 / 0.2e1 - t80 / 0.2e1;
t401 = 0.2e1 * t407;
t400 = t571 + t420;
t399 = t258 * t565 + t578;
t398 = t290 * t564 + t577;
t396 = -t281 * t359 + t557;
t9 = t113 * t90 + t143 * t144 + t80 * t95;
t395 = t9 * qJD(1) + t7 * qJD(2);
t394 = t289 * t424 + t231;
t12 = t110 * t90 + t111 * t80 + t143 * t202;
t393 = t12 * qJD(1) + t13 * qJD(2);
t17 = t134 * t150 + t135 * t153 + t189 * t190;
t392 = t17 * qJD(1) + t14 * qJD(2);
t19 = -t134 * t145 + t135 * t146 + t189 * t244;
t391 = t19 * qJD(1) + t18 * qJD(2);
t36 = t328 ^ 2 * t500 - t145 * t186 + t146 * t187;
t390 = t36 * qJD(1) + t28 * qJD(2);
t389 = t530 + t531;
t388 = -t186 * t356 + t187 * t358;
t174 = t289 * t356 + t291 * t358;
t128 = t520 + t521;
t30 = -t516 / 0.2e1 + t399 * t358 + (t576 - t371) * t356 + t445;
t386 = -qJD(1) * t30 + qJD(3) * t128;
t129 = -t250 * t356 + t258 * t358;
t21 = (t575 + t400) * t358 + (-t520 / 0.2e1 - t521 / 0.2e1 - t450) * t357 + t363;
t385 = qJD(1) * t21 + qJD(3) * t129;
t166 = t281 * t358 + t290 * t356;
t256 = -t273 / 0.2e1;
t301 = pkin(8) * t428;
t39 = t301 - t522 / 0.2e1 - t525 / 0.2e1 - t346 + t256 + (t517 / 0.2e1 + (-t290 / 0.2e1 + t328 / 0.2e1) * t358) * t357;
t384 = -qJD(1) * t39 + qJD(3) * t166;
t167 = t290 * t358 - t517;
t378 = t426 + t451;
t368 = t574 + t378;
t374 = t398 * t356;
t41 = -t347 + t573 + t374 + (t571 + t368) * t358;
t383 = -qJD(1) * t41 + qJD(3) * t167;
t142 = (t358 * t406 + t345) * t357;
t201 = -t356 * t406 + t499;
t381 = qJD(1) * t142 - qJD(3) * t201;
t380 = t451 - t559 / 0.2e1;
t379 = -t150 * qJ(5) / 0.2e1 + t153 * t582;
t377 = t113 * qJ(5) / 0.2e1 - t95 * t579 / 0.2e1;
t302 = pkin(8) * t429;
t154 = t357 * t453 + t255 + t302;
t376 = pkin(3) * t461 - qJD(1) * t154;
t155 = (t571 + t380) * t358;
t375 = pkin(3) * t465 - qJD(1) * t155;
t372 = t250 * t576 + t258 * t578;
t49 = t526 / 0.2e1 + (t571 + t378) * t358 + t414;
t370 = qJD(1) * t49 + t281 * t465;
t265 = t431 + t465;
t268 = -t336 + t508;
t183 = -t268 * qJD(1) + t401;
t207 = -t310 * qJD(3) + 0.2e1 * t408;
t184 = t270 * qJD(3) + t357 * t442;
t185 = -t272 * qJD(3) + t357 * t443;
t369 = t570 - t513 / 0.2e1 + t452;
t3 = t289 * t446 - t291 * t447 + t372 + t377;
t61 = t250 * t258;
t63 = (t258 / 0.2e1 - t499 / 0.2e1 + t418) * t359 + t394;
t367 = -t3 * qJD(1) + t63 * qJD(2) + t61 * qJD(3);
t313 = t354 * t564;
t224 = t313 + (-0.1e1 / 0.2e1 + t569) * t357;
t25 = t231 + t492 / 0.2e1 + (t580 + t422) * t358 + (t359 * t456 + t403) * t356;
t366 = qJD(1) * t25 - qJD(2) * t224 - qJD(3) * t174;
t160 = -t268 * qJD(4) + t359 * t401;
t151 = t493 / 0.2e1 + t484;
t361 = t362 * pkin(8) + t290 * t577 + t244 * t281 / 0.2e1;
t6 = t361 + t379;
t365 = t281 * t290 * qJD(3) + t6 * qJD(1) - t151 * qJD(2);
t283 = t355 + t336;
t161 = t283 * qJD(1) + t407 - t457;
t350 = qJ(5) * qJD(5);
t342 = t352 * qJD(5);
t334 = qJ(5) * t457;
t333 = -t464 / 0.2e1;
t332 = t464 / 0.2e1;
t331 = t463 / 0.2e1;
t320 = t358 * t458;
t312 = t352 * t565;
t285 = t356 * t319;
t282 = t324 * qJ(5);
t277 = t283 * qJD(5);
t266 = -t343 + t430;
t264 = t324 * t356;
t263 = t356 * t464 - t461;
t260 = (t459 - qJD(4) / 0.2e1) * t357;
t242 = -t318 + t442;
t241 = t356 * t458 + t358 * t462;
t240 = -t432 - t443;
t239 = t356 * t462 - t320;
t238 = t352 * qJD(3) + t408;
t227 = -t356 * t456 + t418;
t226 = t312 + (t568 - 0.1e1 / 0.2e1) * t357;
t225 = t312 + t313 + t565;
t217 = t265 * t359;
t216 = t356 * t382;
t215 = t356 * t433 - t320;
t213 = t225 * qJD(5);
t212 = t223 * qJD(5);
t206 = t354 * t435 - t409;
t205 = t352 * t435 + t409;
t188 = t358 * t417 + t425 * t579;
t178 = t193 * qJD(3);
t170 = -t354 * t433 - t472;
t169 = -t352 * t433 + t472;
t152 = -t493 / 0.2e1 + t484;
t148 = qJD(3) * t225 - t284;
t147 = qJD(3) * t223 + t284;
t138 = t146 * qJD(4);
t118 = -t472 + (t354 * t464 + t315) * t359;
t117 = t472 + (t352 * t464 - t315) * t359;
t94 = t358 * t380 + t257 + t275;
t93 = t302 + t256 + (t328 * t358 + t453) * t357;
t62 = t258 * t561 - t428 * t579 + t296 + t394;
t50 = pkin(8) * t423 + t281 * t425 - t526 / 0.2e1 - t498 / 0.2e1 + t414;
t43 = t250 * t424 + t535 / 0.2e1 + t573 + t400 * t358 + t369;
t40 = t358 * t368 + t347 + t374 + t485;
t38 = t301 - t398 * t358 + (t426 + t574) * t356 + t445;
t31 = t516 / 0.2e1 + t258 * t424 + t202 * t566 + t143 * t562 + t445 + t588 * t356;
t26 = -t551 / 0.2e1 + t231 - t492 / 0.2e1 + pkin(5) * t429 + t403 * t356 + t484;
t24 = pkin(4) * t423 + t135 * t562 + t145 * t563 + t298 - t373;
t22 = t399 * t356 - t369 + t417 + t485 + (t575 + t588) * t358;
t11 = t110 * t563 + t552 / 0.2e1 + t111 * t567 - t553 / 0.2e1 + t297 + t358 * t416;
t8 = qJD(3) * t28;
t5 = t361 - t379;
t4 = -t372 + t377 + t589 * t291 / 0.2e1 + (t580 + t111 / 0.2e1) * t289;
t2 = qJD(3) * t14 + qJD(4) * t18;
t1 = qJD(3) * t7 + qJD(4) * t13;
t37 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t435, t311 * qJD(3), 0, -t435, 0, 0, t329 * t463, t329 * t458, 0, 0, t206, -t160, t185, t205, t184, -t435, -qJD(3) * t52 - qJD(4) * t78, qJD(3) * t53 + qJD(4) * t77, -qJD(3) * t33, qJD(3) * t36, t206, t185, t160, -t435, -t184, t205, qJD(3) * t35 + qJD(4) * t47 - t314 * t353, -qJD(3) * t32 - qJD(4) * t29 + t410, -qJD(3) * t34 - qJD(4) * t48 + t277, qJD(3) * t17 + qJD(4) * t19 - qJD(5) * t60, t206, t160, -t185, t205, t184, -t435, t20 * qJD(3) + t45 * qJD(4) + (-t353 * t478 - t357 * t474) * t358, -qJD(3) * t27 - qJD(4) * t44 - t357 * t436 + t277, qJD(3) * t15 + qJD(4) * t16 + qJD(6) * t267 - t410, qJD(3) * t9 + qJD(4) * t12 - qJD(5) * t51 + qJD(6) * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t433, t466, t458, -t433, -t463, 0, -t328 * t458 + t439, t328 * t463 + t438, 0, 0, t118, -t104, -t587, t117, t195, -t260, -t488 + (t356 * t402 - t276) * qJD(3) + t94 * qJD(4), t487 + (t358 * t402 + t449) * qJD(3) + t93 * qJD(4), qJD(3) * t388 - t511 (-pkin(3) * t492 + pkin(8) * t388) * qJD(3) + t390, t118, -t587, t104, -t260, -t195, t117, t509 + (-t190 * t358 - t356 * t396) * qJD(3) + t40 * qJD(4) + t213, qJD(3) * t389 + t24 * qJD(4) - t512, -t510 + (-t190 * t356 + t358 * t396) * qJD(3) + t38 * qJD(4) + t285 (pkin(8) * t389 + t190 * t281) * qJD(3) + t5 * qJD(4) + t50 * qJD(5) + t392, t118, t104, t587, t117, t195, -t260, t524 + (t144 * t358 - t250 * t505 - t289 * t357) * qJD(3) + t22 * qJD(4) + t213 - t436, -t519 + (t144 * t356 + t250 * t496 + t514) * qJD(3) + t31 * qJD(4) + t285 + t358 * t474, t532 + ((-t113 - t516) * t358 + (-t95 + t513) * t356) * qJD(3) + t11 * qJD(4) (t113 * t291 + t144 * t250 + t289 * t95) * qJD(3) + t4 * qJD(4) + t43 * qJD(5) + t26 * qJD(6) + t395; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t171, -t183, t216, t171, t218, t331, qJD(3) * t94 - t138 - t540, qJD(3) * t93 + t481 + t541, 0, 0, -t171, t216, t183, t331, -t218, t171, qJD(3) * t40 - t138 + t545, t24 * qJD(3) - qJD(4) * t483 - t434 - t515, qJD(3) * t38 - t477 - t481 - t544, t5 * qJD(3) + (-pkin(4) * t146 - qJ(5) * t145) * qJD(4) + t134 * qJD(5) + t391, -t171, t183, -t216, t171, t218, t331, qJD(3) * t22 - qJD(4) * t111 + t547, qJD(3) * t31 + qJD(4) * t110 - t477 - t548, t529 + t11 * qJD(3) + (-t356 * t491 + t323) * qJD(4) + t434, t4 * qJD(3) + (t110 * qJ(5) - t111 * t579) * qJD(4) + t90 * qJD(5) + t188 * qJD(6) + t393; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, t216, t161, qJD(3) * t50 + qJD(4) * t134 - t542, 0, 0, 0, 0, 0, 0, t148, t161, -t216, qJD(3) * t43 + qJD(4) * t90 - t543; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t217, -t215, t471, qJD(3) * t26 + qJD(4) * t188 + t546; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t178, 0, 0, 0, 0, 0, 0, 0, 0, 0, t178, 0, 0, 0, 0, 0, 0, 0, 0, 0, t178; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t463, -t458, 0, 0, 0, 0, 0, 0, 0, 0, t240, -t242, t247 (t482 - t559) * qJD(3) - t585, 0, 0, 0, 0, 0, 0, t240, t247, t242, t536 + (t482 + t503) * qJD(3) + t152 * qJD(4) + t486, 0, 0, 0, 0, 0, 0, t240, t242, -t247, t554 + (t174 * t359 - t504) * qJD(3) + t62 * qJD(4) + t226 * qJD(6) + t486; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t241, t239, 0, 0, 0, 0, 0, 0, 0, 0, -t241, 0, -t239, qJD(3) * t152 - qJD(4) * t244 + t319 + t549, 0, 0, 0, 0, 0, 0, -t241, -t239, 0, t62 * qJD(3) + t462 * t592 + t319 + t538; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t241, 0, 0, 0, 0, 0, 0, 0, 0, 0, t241; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t226 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t433, -t466, 0, t433, 0, 0, -t439, -t438, 0, 0, t170, t181, t197, t169, -t586, t260, qJD(4) * t155 + t488, qJD(4) * t154 - t487, t511, -t390, t170, t197, -t181, t260, t586, t169, qJD(4) * t41 - t212 - t509, -qJD(4) * t23 - t441 + t512, qJD(4) * t39 + t285 + t510, qJD(4) * t6 - qJD(5) * t49 - t392, t170, -t181, -t197, t169, -t586, t260, qJD(4) * t21 - t212 - t524, -qJD(4) * t30 + t285 + t519, -qJD(4) * t10 + t441 - t532, -qJD(4) * t3 - qJD(5) * t42 + qJD(6) * t25 - t395; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t585, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t151 - t473 - t536, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4) * t63 - qJD(6) * t224 - t473 - t554; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t317, t292, 0, -t317, 0, 0, -pkin(3) * t479, -pkin(3) * t343, 0, 0, t317, 0, -t292, 0, 0, -t317, -qJD(4) * t167 + t314, 0, -qJD(4) * t166 + t342 (qJD(4) * t290 - t478) * t281, t317, -t292, 0, -t317, 0, 0, qJD(4) * t129 + t314, qJD(4) * t128 + t342, t309 * qJD(6), qJD(4) * t61 - qJD(6) * t174 + t250 * t478; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t191, -t207, -t266, t191, t264, t333, -t375 - t454, -t376 + t455, 0, 0, -t191, -t266, t207, t333, -t264, t191, -t383 - t454, -t523 + t593, -t384 - t455, pkin(8) * t593 + t365, -t191, t207, t266, t191, t264, t333, t385 - t480, -qJD(4) * t289 + t386, -qJD(4) * t592 - t460 - t539 (-t289 * qJ(5) - t291 * t579) * qJD(4) + t291 * qJD(5) + t227 * qJD(6) + t367; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t163, -t266, t238, -t370 + t454, 0, 0, 0, 0, 0, 0, -t163, t238, t266, t250 * t465 + t480 - t489; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t467, qJD(4) * t227 + t366; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t171, t183, -t215, -t171, -t217, t331, -qJD(3) * t155 + t540, -qJD(3) * t154 - t541, 0, 0, t171, -t215, -t183, t331, t217, -t171, -qJD(3) * t41 - t545, qJD(3) * t23 + t515, -qJD(3) * t39 - t477 + t544, -qJD(3) * t6 - t391 - t444, t171, -t183, t215, -t171, -t217, t331, -qJD(3) * t21 + t440 - t547, qJD(3) * t30 + t357 * t476 - t477 + t548, qJD(3) * t10 - t529, qJD(3) * t3 + qJD(6) * t142 - t393 - t444; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t151 - t549, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t63 - t538; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t191, t207, t430, -t191, -t437, t332, t375, t376, 0, 0, t191, t430, -t207, t332, t437, -t191, t383, t523, t384, -t365, t191, -t207, -t430, -t191, -t437, t332, -t385 + t476, -t386 - t475, t539, -qJD(6) * t201 - t367; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), t350, 0, 0, 0, 0, 0, 0, 0, qJD(5), 0, t350; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t324, -t282, 0, 0, 0, 0, 0, 0, 0, -t324, 0, -t282; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t265, t263, 0, t381; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, -t215, -t161, qJD(3) * t49 + t334 + t542, 0, 0, 0, 0, 0, 0, t147, -t161, t215, qJD(3) * t42 + t334 - t440 + t543; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t163, t430, -t238, t370, 0, 0, 0, 0, 0, 0, t163, -t238, -t430, t489 + (-qJD(3) * t250 - qJD(6)) * t356; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t324, t282, 0, 0, 0, 0, 0, 0, 0, t324, 0, t282; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t265; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t218, t216, -t471, -qJD(3) * t25 - qJD(4) * t142 + t319 - t546; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t224 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t479, t343, -t467, qJD(4) * t201 - t366 + t478; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t265, -t263, 0, -t381; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t265; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t37;