% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPP5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPP5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RRPRPP5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:06:25
% EndTime: 2019-03-09 10:06:36
% DurationCPUTime: 6.48s
% Computational Cost: add. (3228->464), mult. (6987->580), div. (0->0), fcn. (3806->4), ass. (0->197)
t490 = MDP(23) - MDP(28);
t431 = cos(qJ(4));
t429 = sin(qJ(4));
t507 = qJD(2) * t429;
t432 = cos(qJ(2));
t508 = qJD(1) * t432;
t374 = t431 * t508 + t507;
t430 = sin(qJ(2));
t491 = qJD(1) * qJD(2);
t478 = t430 * t491;
t332 = t374 * qJD(4) - t429 * t478;
t509 = qJD(1) * t430;
t549 = -t509 - qJD(4);
t528 = t374 * t549;
t552 = t332 + t528;
t557 = t552 * t490;
t556 = t332 - t528;
t484 = t431 * t509;
t501 = qJD(4) * t431;
t555 = t484 + t501;
t554 = MDP(22) + MDP(26);
t553 = MDP(24) + MDP(27);
t540 = pkin(3) + pkin(7);
t479 = t429 * t508;
t333 = qJD(2) * t501 - qJD(4) * t479 - t431 * t478;
t505 = qJD(2) * t431;
t376 = -t479 + t505;
t527 = t376 * t549;
t551 = -t333 + t527;
t477 = t432 * t491;
t550 = qJ(5) * t477 - qJD(5) * t549;
t415 = pkin(7) * t509;
t495 = pkin(3) * t509 + qJD(3) + t415;
t548 = -0.2e1 * t491;
t427 = t430 ^ 2;
t428 = t432 ^ 2;
t547 = MDP(5) * (t427 - t428);
t546 = t333 * qJ(6) + t374 * qJD(6);
t545 = t555 * t549;
t544 = 0.2e1 * t550;
t542 = qJD(5) * t431 - t495;
t434 = -pkin(2) - pkin(8);
t476 = -qJ(3) * t430 - pkin(1);
t368 = t432 * t434 + t476;
t346 = t368 * qJD(1);
t348 = qJD(2) * t434 + t495;
t308 = -t429 * t346 + t431 * t348;
t493 = qJD(5) - t308;
t534 = qJ(5) * t431;
t539 = pkin(4) + pkin(5);
t541 = t429 * t539 - t534;
t535 = qJ(5) * t429;
t450 = t431 * t539 + t535;
t373 = t376 ^ 2;
t405 = t549 ^ 2;
t538 = qJD(2) * pkin(2);
t537 = qJ(5) * t333;
t536 = qJ(5) * t374;
t425 = qJD(2) * qJD(3);
t386 = pkin(7) * t478 - t425;
t353 = -pkin(3) * t478 - t386;
t292 = t333 * pkin(4) + t332 * qJ(5) - t376 * qJD(5) + t353;
t533 = t292 * t429;
t532 = t292 * t431;
t531 = t332 * t431;
t530 = t353 * t429;
t529 = t353 * t431;
t526 = t376 * t432;
t525 = t549 * t434;
t524 = t429 * t430;
t523 = t429 * t432;
t522 = t430 * t431;
t435 = qJD(2) ^ 2;
t521 = t430 * t435;
t520 = t432 * t435;
t436 = qJD(1) ^ 2;
t519 = t432 * t436;
t518 = qJ(6) + t434;
t471 = qJD(4) * t518;
t420 = pkin(2) * t509;
t463 = pkin(8) * t430 - qJ(3) * t432;
t356 = qJD(1) * t463 + t420;
t417 = pkin(7) * t508;
t382 = pkin(3) * t508 + t417;
t473 = -t429 * t356 + t382 * t431;
t497 = qJD(6) * t431;
t517 = (-qJ(6) * t524 - t432 * t539) * qJD(1) - t473 - t429 * t471 + t497;
t513 = t431 * t356 + t429 * t382;
t314 = qJ(5) * t508 + t513;
t498 = qJD(6) * t429;
t516 = qJ(6) * t484 + t431 * t471 - t314 + t498;
t515 = -t450 * t549 - t542;
t465 = pkin(4) * t431 + t535;
t514 = t465 * t549 + t542;
t309 = t431 * t346 + t429 * t348;
t395 = t540 * t430;
t512 = t431 * t368 + t429 * t395;
t396 = t540 * t432;
t392 = -pkin(2) * t432 + t476;
t361 = qJD(1) * t392;
t426 = qJD(2) * qJ(3);
t506 = qJD(2) * t430;
t504 = qJD(2) * t432;
t503 = qJD(3) * t430;
t502 = qJD(4) * t429;
t500 = qJD(4) * t432;
t499 = qJD(5) * t429;
t496 = t432 * MDP(19);
t300 = qJ(6) * t376 + t308;
t494 = qJD(5) - t300;
t360 = t426 + t382;
t456 = qJ(5) * t376 - t360;
t302 = -t374 * t539 + qJD(6) + t456;
t492 = qJD(6) + t302;
t489 = t429 * t525;
t488 = t431 * t525;
t487 = t430 * t519;
t411 = pkin(2) * t478;
t439 = qJD(2) * t463 - t503;
t325 = qJD(1) * t439 + t411;
t409 = pkin(7) * t477;
t367 = pkin(3) * t477 + t409;
t486 = -t431 * t325 - t348 * t501 - t429 * t367;
t318 = t430 * qJ(5) + t512;
t485 = t429 * t509;
t483 = t429 * t506;
t482 = t434 * t504;
t384 = t549 * t502;
t481 = t429 * t500;
t480 = t431 * t500;
t475 = pkin(1) * t548;
t474 = qJD(3) - t538;
t472 = -t429 * t368 + t395 * t431;
t400 = t431 * t477;
t470 = -qJD(2) * t374 + t400;
t398 = t429 * t477;
t469 = -qJD(2) * t376 - t398;
t467 = -t429 * t325 - t346 * t501 - t348 * t502 + t431 * t367;
t301 = qJ(6) * t374 + t309;
t289 = -pkin(4) * t477 - t467;
t403 = t549 * qJ(5);
t305 = -t403 + t309;
t466 = -t305 * t509 + t289;
t464 = -pkin(4) * t429 + t534;
t443 = t346 * t502 + t486;
t285 = -t443 + t550;
t283 = t285 + t546;
t295 = t539 * t549 + t494;
t298 = t301 - t403;
t462 = t283 * t429 + t555 * t298 + (t485 + t502) * t295;
t303 = pkin(4) * t549 + t493;
t461 = t303 * t429 + t305 * t431;
t460 = -qJD(1) * t428 - t430 * t549;
t459 = -0.2e1 * qJD(2) * t361;
t419 = pkin(2) * t506;
t340 = t419 + t439;
t383 = t540 * t504;
t455 = -t429 * t340 - t368 * t501 + t383 * t431 - t395 * t502;
t447 = -qJ(3) * t504 - t503;
t343 = qJD(1) * t447 + t411;
t359 = t419 + t447;
t452 = pkin(7) * t435 + qJD(1) * t359 + t343;
t286 = -pkin(5) * t333 - t292;
t449 = -t286 * t429 - t302 * t501;
t448 = -t286 * t431 + t302 * t502;
t445 = qJ(6) * t332 - t467;
t442 = t431 * t340 - t368 * t502 + t429 * t383 + t395 * t501;
t291 = qJ(5) * t504 + t430 * qJD(5) + t442;
t441 = -t308 * t549 + t443;
t438 = -t477 * t539 + t445;
t390 = t415 + t474;
t393 = -t417 - t426;
t437 = -t386 * t432 + (t390 * t432 + (t393 + t417) * t430) * qJD(2);
t391 = qJ(3) - t464;
t389 = t518 * t431;
t388 = t518 * t429;
t387 = t434 * t400;
t381 = t540 * t506;
t379 = -qJ(3) * t508 + t420;
t364 = -qJ(3) - t541;
t347 = t361 * t509;
t339 = t432 * t465 + t396;
t323 = pkin(4) * t376 + t536;
t322 = -t432 * t450 - t396;
t319 = -pkin(4) * t430 - t472;
t316 = -pkin(4) * t508 - t473;
t315 = qJ(6) * t431 * t432 + t318;
t313 = -t376 * t539 - t536;
t312 = pkin(4) * t374 - t456;
t310 = qJ(6) * t523 - t430 * t539 - t472;
t306 = (qJD(4) * t464 + t499) * t432 + (-t465 - t540) * t506;
t299 = (qJD(4) * t541 - t499) * t432 + (t450 + t540) * t506;
t294 = -pkin(4) * t504 - t455;
t288 = t432 * t497 + (-t505 * t430 - t481) * qJ(6) + t291;
t287 = -qJ(6) * t483 + (qJ(6) * t501 - qJD(2) * t539 + t498) * t432 - t455;
t284 = -qJD(6) * t376 + t438;
t1 = [t547 * t548 + (-t430 * t452 + t432 * t459) * MDP(13) + (t430 * t459 + t432 * t452) * MDP(12) + (-t455 * t549 - t381 * t374 + t396 * t333 + (-t360 * t505 + t467) * t430 + (-t360 * t502 + t529 + (qJD(1) * t472 + t308) * qJD(2)) * t432) * MDP(20) + (t442 * t549 - t381 * t376 - t396 * t332 + ((qJD(2) * t360 + qJD(4) * t346) * t429 + t486) * t430 + (-t360 * t501 - t530 + (-qJD(1) * t512 - t309) * qJD(2)) * t432) * MDP(21) + (t294 * t549 + t306 * t374 + t333 * t339 + (-t312 * t505 - t289) * t430 + (-t312 * t502 + t532 + (-qJD(1) * t319 - t303) * qJD(2)) * t432) * MDP(22) + (-t291 * t549 - t306 * t376 + t332 * t339 + (-t312 * t507 + t285) * t430 + (t312 * t501 + t533 + (qJD(1) * t318 + t305) * qJD(2)) * t432) * MDP(24) + (t549 * t480 - t332 * t430 + (t429 * t460 + t526) * qJD(2)) * MDP(17) + (t287 * t549 - t299 * t374 - t322 * t333 + (t302 * t505 - t284) * t430 + ((-qJD(1) * t310 - t295) * qJD(2) + t448) * t432) * MDP(26) + (-t288 * t549 + t299 * t376 - t322 * t332 + (t302 * t507 + t283) * t430 + ((qJD(1) * t315 + t298) * qJD(2) + t449) * t432) * MDP(27) + ((-t374 * t429 + t376 * t431) * t506 + (t531 + t333 * t429 + (t374 * t431 + t376 * t429) * qJD(4)) * t432) * MDP(16) + (-t549 * t481 - t333 * t430 + (-t374 * t432 + t431 * t460) * qJD(2)) * MDP(18) + (t332 * t523 + (-t480 + t483) * t376) * MDP(15) + (-pkin(7) * t520 + t430 * t475) * MDP(9) - MDP(7) * t521 + (pkin(7) * t521 + t432 * t475) * MDP(10) + (-t287 * t376 + t288 * t374 + t310 * t332 + t315 * t333 + (-t295 * t429 - t298 * t431) * t506 + (t283 * t431 + t284 * t429 + (t295 * t431 - t298 * t429) * qJD(4)) * t432) * MDP(28) + (-t291 * t374 + t294 * t376 - t318 * t333 - t319 * t332 + t461 * t506 + (-t285 * t431 - t289 * t429 + (-t303 * t431 + t305 * t429) * qJD(4)) * t432) * MDP(23) + (t285 * t318 + t289 * t319 + t291 * t305 + t292 * t339 + t294 * t303 + t306 * t312) * MDP(25) + (t283 * t315 + t284 * t310 + t286 * t322 + t287 * t295 + t288 * t298 + t299 * t302) * MDP(29) + t437 * MDP(11) + (pkin(7) * t437 + t343 * t392 + t359 * t361) * MDP(14) + (-t549 + t509) * qJD(2) * t496 + 0.2e1 * t430 * MDP(4) * t477 + MDP(6) * t520; -MDP(4) * t487 + t436 * t547 + (-t379 * t508 + t347) * MDP(12) + 0.2e1 * t425 * MDP(13) + (-qJ(3) * t386 - qJD(3) * t393 - t361 * t379) * MDP(14) + (t429 * t527 - t531) * MDP(15) + (t556 * t429 + t551 * t431) * MDP(16) + (t384 + t400) * MDP(17) + (t374 * t508 - t398 + t545) * MDP(18) + (t387 + qJ(3) * t333 + t530 + t473 * t549 + t495 * t374 + (t360 * t431 + t489) * qJD(4)) * MDP(20) + (-qJ(3) * t332 + t529 - t513 * t549 + t495 * t376 + (-t360 * t429 + t488) * qJD(4)) * MDP(21) + (t533 - t316 * t549 + t333 * t391 + t387 - t514 * t374 + (t312 * t431 + t489) * qJD(4)) * MDP(22) + (t314 * t374 - t316 * t376 + (t332 * t434 + (-t374 * t434 - t305) * qJD(4) + t466) * t431 + (-t303 * t509 - t333 * t434 - t285 + (t376 * t434 - t303) * qJD(4)) * t429) * MDP(23) + (-t532 + t314 * t549 + t332 * t391 + t514 * t376 + (t312 * t429 - t488) * qJD(4)) * MDP(24) + (t292 * t391 - t303 * t316 - t305 * t314 - t514 * t312 + (qJD(4) * t461 + t285 * t429 - t289 * t431) * t434) * MDP(25) + (-t333 * t364 + t515 * t374 - t517 * t549 + t449) * MDP(26) + (-t332 * t364 - t515 * t376 - t516 * t549 - t448) * MDP(27) + (-t284 * t431 - t332 * t389 + t333 * t388 + t374 * t516 + t376 * t517 + t462) * MDP(28) + (t283 * t388 - t284 * t389 + t286 * t364 - t295 * t517 + t298 * t516 - t302 * t515) * MDP(29) + (((-t393 - t426) * t430 + (-t390 + t474) * t432) * MDP(11) + (t361 * t432 + t379 * t430) * MDP(13) + (-t393 * t430 + (-t390 - t538) * t432) * pkin(7) * MDP(14) + (t524 * t549 - t526) * MDP(17) + t549 * t496 + (-t308 * t432 + t360 * t522) * MDP(20) + (t309 * t432 + (-t360 * t430 - t482) * t429) * MDP(21) + (t303 * t432 + t312 * t522) * MDP(22) + (-t305 * t432 + (t312 * t430 + t482) * t429) * MDP(24) + (-t302 * t522 + (qJD(2) * t389 + t295) * t432) * MDP(26) + (-t302 * t524 + (qJD(2) * t388 - t298) * t432) * MDP(27)) * qJD(1) + (MDP(9) * t430 * t436 + MDP(10) * t519) * pkin(1); MDP(12) * t487 + (-t427 * t436 - t435) * MDP(13) + (qJD(2) * t393 + t347 + t409) * MDP(14) + t470 * MDP(20) + t469 * MDP(21) - t312 * qJD(2) * MDP(25) + (qJD(2) * t302 + t462) * MDP(29) + ((qJD(4) * t305 - t466) * MDP(25) - t284 * MDP(29) - MDP(21) * t405 + t557) * t431 + (-MDP(20) * t405 + t285 * MDP(25) - t490 * t333 - (t303 * MDP(25) + t376 * t490) * t549) * t429 + t554 * (t485 * t549 + t384 + t470) + t553 * (-t469 - t545); -t552 * MDP(17) - t333 * MDP(18) + t441 * MDP(21) + (pkin(4) * t332 - t537) * MDP(23) + (-t441 + t544) * MDP(24) + (-pkin(4) * t289 + qJ(5) * t285 - t303 * t309 + t305 * t493 - t312 * t323) * MDP(25) + (-t301 * t549 - t445) * MDP(26) + (t300 * t549 - t443 + t544 + t546) * MDP(27) + (-t332 * t539 + t537) * MDP(28) + (qJ(5) * t283 - t284 * t539 - t295 * t301 + t298 * t494 - t302 * t313) * MDP(29) + (0.2e1 * pkin(4) * MDP(22) + 0.2e1 * t539 * MDP(26) + MDP(19)) * t477 + (-t549 * MDP(18) - t360 * MDP(20) - t312 * MDP(22) + (t305 - t309) * MDP(23) + t323 * MDP(24) + t492 * MDP(26) - t313 * MDP(27) + (-t298 + t301) * MDP(28) + MDP(16) * t376) * t376 + (t376 * MDP(15) + t360 * MDP(21) - t323 * MDP(22) + (t303 - t493) * MDP(23) - t312 * MDP(24) + t313 * MDP(26) + t302 * MDP(27) + (-t295 + t494) * MDP(28) - MDP(16) * t374) * t374 + (MDP(20) + MDP(22)) * (-t309 * t549 + t467); (t305 * t549 + t312 * t376 + t289) * MDP(25) + (t298 * t549 - t492 * t376 + t438) * MDP(29) + t554 * (t374 * t376 - t477) + t553 * (-t405 - t373) - t557; t551 * MDP(26) - t556 * MDP(27) + (-t374 ^ 2 - t373) * MDP(28) + (t295 * t376 - t298 * t374 + t286) * MDP(29);];
tauc  = t1;
