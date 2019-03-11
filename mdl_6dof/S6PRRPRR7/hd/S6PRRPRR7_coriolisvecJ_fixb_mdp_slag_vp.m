% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRPRR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:34:21
% EndTime: 2019-03-08 22:34:30
% DurationCPUTime: 5.22s
% Computational Cost: add. (2419->434), mult. (5798->607), div. (0->0), fcn. (4052->10), ass. (0->210)
t566 = pkin(4) + pkin(8);
t446 = sin(qJ(2));
t441 = sin(pkin(6));
t531 = qJD(1) * t441;
t503 = t446 * t531;
t404 = qJD(2) * pkin(8) + t503;
t445 = sin(qJ(3));
t449 = cos(qJ(3));
t442 = cos(pkin(6));
t530 = qJD(1) * t442;
t364 = t449 * t404 + t445 * t530;
t526 = qJD(2) * t449;
t349 = pkin(4) * t526 + t364;
t438 = qJD(3) * qJ(4);
t338 = t438 + t349;
t448 = cos(qJ(5));
t444 = sin(qJ(5));
t524 = qJD(3) * t444;
t392 = t448 * t526 + t524;
t320 = pkin(5) * t392 + t338;
t447 = cos(qJ(6));
t496 = t444 * t526;
t522 = qJD(3) * t448;
t394 = -t496 + t522;
t443 = sin(qJ(6));
t555 = t394 * t443;
t334 = t447 * t392 + t555;
t580 = t320 * t334;
t527 = qJD(2) * t445;
t430 = qJD(5) + t527;
t421 = qJD(6) + t430;
t579 = t334 * t421;
t472 = t392 * t443 - t447 * t394;
t578 = t421 * t472;
t450 = cos(qJ(2));
t528 = qJD(2) * t441;
t495 = qJD(1) * t528;
t480 = t450 * t495;
t525 = qJD(3) * t442;
t577 = qJD(1) * t525 + t480;
t363 = t445 * t404 - t449 * t530;
t576 = qJD(4) + t363;
t569 = qJD(5) + qJD(6);
t512 = qJD(2) * qJD(3);
t493 = t445 * t512;
t359 = -qJD(5) * t392 + t444 * t493;
t451 = -pkin(3) - pkin(9);
t514 = pkin(4) * t527 + t576;
t332 = qJD(3) * t451 + t514;
t491 = -qJ(4) * t445 - pkin(2);
t388 = t449 * t451 + t491;
t502 = t450 * t531;
t350 = qJD(2) * t388 - t502;
t312 = t332 * t444 + t350 * t448;
t521 = qJD(3) * t449;
t329 = t404 * t521 + t577 * t445;
t492 = t449 * t512;
t321 = pkin(4) * t492 + t329;
t476 = pkin(9) * t445 - qJ(4) * t449;
t520 = qJD(4) * t445;
t456 = qJD(3) * t476 - t520;
t535 = pkin(3) * t493 + t446 * t495;
t330 = qJD(2) * t456 + t535;
t485 = t448 * t321 - t330 * t444;
t455 = -qJD(5) * t312 + t485;
t296 = pkin(5) * t492 - pkin(10) * t359 + t455;
t518 = qJD(5) * t448;
t360 = qJD(3) * t518 - qJD(5) * t496 - t448 * t493;
t507 = -t444 * t321 - t448 * t330 - t332 * t518;
t519 = qJD(5) * t444;
t460 = -t350 * t519 - t507;
t297 = -pkin(10) * t360 + t460;
t487 = t447 * t296 - t443 * t297;
t575 = t320 * t472 + t487;
t574 = MDP(27) * t492 + (-t334 ^ 2 + t472 ^ 2) * MDP(24) - t334 * MDP(23) * t472;
t439 = t445 ^ 2;
t440 = t449 ^ 2;
t573 = (t439 - t440) * MDP(6);
t523 = qJD(3) * t445;
t433 = pkin(3) * t523;
t366 = t433 + t456;
t401 = t566 * t521;
t544 = t448 * t450;
t572 = -(-t444 * t446 + t445 * t544) * t531 - t366 * t444 + t448 * t401;
t411 = t566 * t445;
t548 = t444 * t450;
t571 = (t445 * t548 + t446 * t448) * t531 - t448 * t366 + t388 * t519 - t444 * t401 - t411 * t518;
t397 = t444 * t411;
t536 = t448 * t388 + t397;
t357 = -t438 - t364;
t356 = -qJD(3) * pkin(3) + t576;
t462 = -qJ(4) * t521 - t520;
t341 = qJD(2) * t462 + t535;
t376 = t433 + t462;
t452 = qJD(3) ^ 2;
t568 = qJD(2) * (-t376 + t503) - pkin(8) * t452 - t341;
t483 = t359 * t443 + t447 * t360;
t304 = -qJD(6) * t472 + t483;
t567 = t449 * t569;
t565 = pkin(10) - t451;
t564 = qJD(2) * pkin(2);
t311 = t448 * t332 - t350 * t444;
t306 = -pkin(10) * t394 + t311;
t301 = pkin(5) * t430 + t306;
t563 = t301 * t447;
t307 = -pkin(10) * t392 + t312;
t562 = t307 * t447;
t437 = qJD(3) * qJD(4);
t505 = -t404 * t523 + t577 * t449;
t323 = -t437 - t505;
t317 = -pkin(4) * t493 - t323;
t561 = t317 * t444;
t560 = t317 * t448;
t559 = t359 * t448;
t408 = -pkin(3) * t449 + t491;
t529 = qJD(2) * t408;
t365 = -t502 + t529;
t558 = t365 * t446;
t557 = t392 * t430;
t556 = t394 * t430;
t554 = t394 * t449;
t553 = t430 * t448;
t552 = t430 * t451;
t551 = t441 * t446;
t550 = t443 * t444;
t549 = t444 * t445;
t547 = t445 * t452;
t546 = t447 * t448;
t545 = t448 * t449;
t543 = t449 * t452;
t467 = pkin(5) * t449 - pkin(10) * t549;
t489 = pkin(10) * t449 - t388;
t542 = t467 * qJD(3) + (t448 * t489 - t397) * qJD(5) + t572;
t517 = qJD(5) * t449;
t498 = t444 * t517;
t541 = -(t445 * t522 + t498) * pkin(10) + t571;
t395 = t443 * t448 + t444 * t447;
t342 = t569 * t395;
t471 = -t546 + t550;
t540 = -t342 * t421 - t471 * t492;
t461 = t395 * t445;
t369 = qJD(2) * t461;
t539 = -t342 - t369;
t499 = t448 * t527;
t516 = qJD(6) * t443;
t538 = -t443 * t519 - t444 * t516 + t447 * t499 - t527 * t550 + t569 * t546;
t434 = pkin(3) * t527;
t371 = qJD(2) * t476 + t434;
t537 = t444 * t349 + t448 * t371;
t504 = -pkin(5) * t448 - pkin(4);
t534 = pkin(5) * t518 - t504 * t527 + t576;
t412 = t566 * t449;
t515 = qJD(6) * t447;
t511 = -MDP(10) + MDP(13);
t510 = MDP(11) - MDP(14);
t509 = t445 * t551;
t453 = qJD(2) ^ 2;
t508 = t445 * t449 * t453;
t506 = t447 * t359 - t443 * t360 - t392 * t515;
t501 = t446 * t528;
t500 = t450 * t528;
t497 = t448 * t517;
t407 = t565 * t448;
t305 = t307 * t516;
t486 = t443 * t296 - t305;
t484 = t448 * t349 - t371 * t444;
t481 = qJD(6) * t301 + t297;
t417 = t445 * t492;
t406 = t565 * t444;
t478 = qJD(2) * t467 - qJD(6) * t406 - t565 * t519 + t484;
t477 = pkin(10) * t499 + t569 * t407 + t537;
t299 = t301 * t443 + t562;
t398 = t448 * t411;
t322 = pkin(5) * t445 + t444 * t489 + t398;
t328 = -pkin(10) * t545 + t536;
t475 = t322 * t443 + t328 * t447;
t382 = -t442 * t449 + t509;
t346 = t382 * t448 + t441 * t548;
t464 = -t382 * t444 + t441 * t544;
t474 = t346 * t447 + t443 * t464;
t473 = t346 * t443 - t447 * t464;
t470 = -qJD(2) * t440 + t430 * t445;
t469 = t430 * t444;
t466 = qJD(3) * t363 + t505;
t465 = qJD(3) * t364 - t329;
t383 = t442 * t445 + t449 * t551;
t463 = t338 * t445 + t451 * t521;
t303 = -t394 * t516 + t506;
t405 = -t502 - t564;
t458 = qJD(3) * (t405 + t502 - t564);
t457 = qJD(3) * (-t365 - t502 - t529);
t454 = -t323 * t449 + t329 * t445 + (t356 * t449 + t357 * t445) * qJD(3);
t431 = pkin(5) * t444 + qJ(4);
t418 = t448 * t492;
t400 = t566 * t523;
t399 = -qJ(4) * t526 + t434;
t381 = pkin(5) * t545 + t412;
t373 = t395 * t449;
t372 = t471 * t449;
t352 = t365 * t527;
t351 = -pkin(5) * t498 + (-pkin(8) + t504) * t523;
t345 = qJD(3) * t383 + t445 * t500;
t344 = -qJD(3) * t509 + (t500 + t525) * t449;
t316 = t395 * t567 - t471 * t523;
t315 = qJD(3) * t461 + t471 * t567;
t310 = pkin(5) * t360 + t317;
t309 = qJD(5) * t346 + t345 * t444 + t448 * t501;
t308 = qJD(5) * t464 + t345 * t448 - t444 * t501;
t298 = -t307 * t443 + t563;
t1 = [(-t323 * t383 + t329 * t382 - t344 * t357 + t345 * t356) * MDP(15) + (t308 * t430 + t344 * t392 + t360 * t383) * MDP(21) + (-t309 * t430 + t344 * t394 + t359 * t383) * MDP(22) + ((-qJD(6) * t473 + t308 * t447 - t309 * t443) * t421 + t344 * t334 + t383 * t304) * MDP(28) + (-(qJD(6) * t474 + t308 * t443 + t309 * t447) * t421 - t344 * t472 + t383 * t303) * MDP(29) + (t344 * t449 + t345 * t445) * MDP(12) * qJD(2) + (t511 * t345 - t510 * t344 + (-t383 * t445 * MDP(12) + (t382 * MDP(12) + t346 * MDP(21) + MDP(22) * t464 + MDP(28) * t474 - MDP(29) * t473) * t449) * qJD(2)) * qJD(3) + (-t341 * t450 * MDP(15) + (MDP(15) * t558 + (t445 * t511 - t449 * t510) * t450 * qJD(3)) * qJD(2) + (-t450 * MDP(4) + (t445 * t510 + t449 * t511 - MDP(3)) * t446) * t453) * t441; MDP(7) * t543 - MDP(8) * t547 - 0.2e1 * t512 * t573 + (-pkin(8) * t543 + t445 * t458) * MDP(10) + (pkin(8) * t547 + t449 * t458) * MDP(11) + ((-t439 - t440) * t480 + t454) * MDP(12) + (t445 * t457 - t568 * t449) * MDP(13) + (t568 * t445 + t449 * t457) * MDP(14) + (t341 * t408 + t365 * t376 + (-t558 + (-t356 * t445 + t357 * t449) * t450) * t531 + t454 * pkin(8)) * MDP(15) + (-t359 * t444 * t449 + (t444 * t523 - t497) * t394) * MDP(16) + ((-t392 * t444 + t394 * t448) * t523 + (-t559 + t360 * t444 + (t392 * t448 + t394 * t444) * qJD(5)) * t449) * MDP(17) + (-t430 * t497 + t359 * t445 + (t444 * t470 + t554) * qJD(3)) * MDP(18) + (t430 * t498 - t360 * t445 + (-t392 * t449 + t448 * t470) * qJD(3)) * MDP(19) + (t430 * t521 + t417) * MDP(20) + (t412 * t360 - t400 * t392 + (-t338 * t522 + t485) * t445 + t572 * t430 + (-t312 * t445 - t536 * t430) * qJD(5) + (-t392 * t502 - t338 * t519 + t560 + ((-t388 * t444 + t398) * qJD(2) + t311) * qJD(3)) * t449) * MDP(21) + (t412 * t359 - t400 * t394 + ((qJD(3) * t338 + qJD(5) * t350) * t444 + t507) * t445 + t571 * t430 + (-t394 * t502 - t338 * t518 - t561 + (-qJD(2) * t536 - t312) * qJD(3)) * t449) * MDP(22) + (-t303 * t373 - t315 * t472) * MDP(23) + (t303 * t372 + t304 * t373 - t315 * t334 - t316 * t472) * MDP(24) + (t303 * t445 + t315 * t421 + (-qJD(2) * t373 - t472) * t521) * MDP(25) + (-t304 * t445 + t316 * t421 + (qJD(2) * t372 - t334) * t521) * MDP(26) + (t421 * t521 + t417) * MDP(27) + (t487 * t445 + t351 * t334 + t381 * t304 - t310 * t372 - t320 * t316 + (t443 * t541 + t447 * t542) * t421 + (-t299 * t445 - t421 * t475) * qJD(6) + (-t334 * t502 + ((t322 * t447 - t328 * t443) * qJD(2) + t298) * qJD(3)) * t449) * MDP(28) + (-(t447 * t481 + t486) * t445 - t351 * t472 + t381 * t303 - t310 * t373 + t320 * t315 + ((-qJD(6) * t322 + t541) * t447 + (qJD(6) * t328 - t542) * t443) * t421 + (t472 * t502 + (-qJD(2) * t475 - t299) * qJD(3)) * t449) * MDP(29) + 0.2e1 * MDP(5) * t417; -MDP(5) * t508 + t453 * t573 + (-t405 * t527 + t465) * MDP(10) - t466 * MDP(11) + (t352 - t465) * MDP(13) + (0.2e1 * t437 + (t365 * t449 + t399 * t445) * qJD(2) + t466) * MDP(14) + (-pkin(3) * t329 - qJ(4) * t323 - t356 * t364 - t357 * t576 - t365 * t399) * MDP(15) + (-t394 * t469 + t559) * MDP(16) + ((-t360 - t556) * t448 + (-t359 + t557) * t444) * MDP(17) + (-t430 * t519 + t418 + (-t430 * t549 - t554) * qJD(2)) * MDP(18) + (-t430 * t518 + (-t445 * t553 + (t392 - t524) * t449) * qJD(2)) * MDP(19) + (qJ(4) * t360 + t561 - t484 * t430 + t514 * t392 + (t338 * t448 - t444 * t552) * qJD(5) + (-t311 * t449 + t448 * t463) * qJD(2)) * MDP(21) + (qJ(4) * t359 + t560 + t537 * t430 + t514 * t394 + (-t338 * t444 - t448 * t552) * qJD(5) + (t312 * t449 - t444 * t463) * qJD(2)) * MDP(22) + (-t303 * t471 - t472 * t539) * MDP(23) + (-t303 * t395 + t304 * t471 - t334 * t539 + t472 * t538) * MDP(24) + (-t369 * t421 + t540) * MDP(25) - t538 * t421 * MDP(26) + (t431 * t304 + t310 * t395 + (t443 * t477 - t447 * t478) * t421 + t534 * t334 + t538 * t320) * MDP(28) + (t431 * t303 - t310 * t471 + (t443 * t478 + t447 * t477) * t421 - t534 * t472 + t539 * t320) * MDP(29) + (-t405 * MDP(11) - t399 * MDP(13) - t430 * MDP(20) + t472 * MDP(25) + (-qJD(3) * t395 + t334) * MDP(26) - t421 * MDP(27) + ((t406 * t443 - t407 * t447) * qJD(3) - t298) * MDP(28) + (-(-t406 * t447 - t407 * t443) * qJD(3) + t299) * MDP(29)) * t526; MDP(13) * t508 + (-t439 * t453 - t452) * MDP(14) + (t352 + t329) * MDP(15) + t418 * MDP(21) + t540 * MDP(28) + (-t369 * MDP(28) - MDP(29) * t538) * t421 + (-MDP(21) * t469 - MDP(22) * t553) * t430 + (t357 * MDP(15) - t392 * MDP(21) + (-t394 - t496) * MDP(22) - t334 * MDP(28) + (-t395 * t526 + t472) * MDP(29)) * qJD(3); t394 * t392 * MDP(16) + (-t392 ^ 2 + t394 ^ 2) * MDP(17) + (t359 + t557) * MDP(18) + (-t360 + t556) * MDP(19) + MDP(20) * t492 + (t312 * t430 - t338 * t394 + t455) * MDP(21) + (t311 * t430 + t338 * t392 - t460) * MDP(22) + (t303 + t579) * MDP(25) + (-t304 - t578) * MDP(26) + (-(-t306 * t443 - t562) * t421 - t299 * qJD(6) + (-t334 * t394 - t421 * t516 + t447 * t492) * pkin(5) + t575) * MDP(28) + (t580 + t305 + (-t307 * t421 - t296) * t443 + (t306 * t421 - t481) * t447 + (t394 * t472 - t421 * t515 - t443 * t492) * pkin(5)) * MDP(29) + t574; (t506 + t579) * MDP(25) + (-t483 - t578) * MDP(26) + (t299 * t421 + t575) * MDP(28) + (-t447 * t297 + t298 * t421 - t486 + t580) * MDP(29) + (-MDP(25) * t555 + MDP(26) * t472 - MDP(28) * t299 - MDP(29) * t563) * qJD(6) + t574;];
tauc  = t1;
