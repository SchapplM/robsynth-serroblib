% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRPRR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:07:58
% EndTime: 2019-03-08 22:08:10
% DurationCPUTime: 8.46s
% Computational Cost: add. (5851->493), mult. (17694->718), div. (0->0), fcn. (14949->14), ass. (0->221)
t474 = cos(pkin(7));
t482 = cos(qJ(3));
t483 = cos(qJ(2));
t575 = t482 * t483;
t478 = sin(qJ(3));
t479 = sin(qJ(2));
t580 = t478 * t479;
t500 = -t474 * t580 + t575;
t557 = qJD(3) * t482;
t538 = t474 * t557;
t472 = sin(pkin(6));
t565 = qJD(1) * t472;
t615 = pkin(2) * t538 - t500 * t565;
t471 = sin(pkin(7));
t603 = pkin(9) + qJ(4);
t534 = t603 * t478;
t555 = qJD(4) * t482;
t614 = (-qJD(3) * t534 + t555) * t471 + t615;
t578 = t479 * t482;
t579 = t478 * t483;
t502 = -t474 * t578 - t579;
t431 = t502 * t565;
t584 = t474 * t478;
t587 = t471 * t482;
t436 = pkin(2) * t584 + t587 * t603;
t556 = qJD(4) * t478;
t613 = -qJD(3) * t436 - t471 * t556 - t431;
t470 = sin(pkin(13));
t473 = cos(pkin(13));
t570 = t613 * t470 + t614 * t473;
t504 = t470 * t482 + t473 * t478;
t444 = t504 * t471;
t439 = qJD(3) * t444;
t559 = qJD(3) * t471;
t585 = t473 * t482;
t440 = (-t470 * t478 + t585) * t559;
t536 = t479 * t565;
t558 = qJD(3) * t478;
t540 = t471 * t558;
t612 = pkin(3) * t540 + pkin(4) * t439 - pkin(10) * t440 - t471 * t536;
t467 = t471 ^ 2;
t611 = (t478 * t482 * MDP(5) - (t478 ^ 2 - t482 ^ 2) * MDP(6)) * t467;
t563 = qJD(2) * t471;
t541 = t478 * t563;
t546 = t471 * t585;
t438 = qJD(2) * t546 - t470 * t541;
t434 = qJD(5) - t438;
t441 = t504 * t563;
t481 = cos(qJ(5));
t561 = qJD(2) * t474;
t524 = qJD(3) + t561;
t459 = t481 * t524;
t477 = sin(qJ(5));
t406 = t441 * t477 - t459;
t405 = qJD(6) + t406;
t571 = t614 * t470 - t613 * t473;
t609 = -t477 * t570 + t612 * t481;
t423 = (pkin(2) * t482 + pkin(3)) * t474 - t471 * t534;
t390 = t470 * t423 + t473 * t436;
t377 = pkin(10) * t474 + t390;
t588 = t471 * t478;
t443 = t470 * t588 - t546;
t543 = -pkin(3) * t482 - pkin(2);
t400 = pkin(4) * t443 - pkin(10) * t444 + t471 * t543;
t551 = qJD(5) * t481;
t553 = qJD(5) * t477;
t608 = t377 * t553 - t400 * t551 - t612 * t477 - t481 * t570;
t607 = t481 * t377 + t477 * t400;
t535 = t483 * t565;
t602 = qJD(2) * pkin(2);
t456 = t535 + t602;
t475 = cos(pkin(6));
t564 = qJD(1) * t475;
t542 = t471 * t564;
t499 = -t456 * t474 - t542;
t449 = pkin(9) * t563 + t536;
t576 = t482 * t449;
t606 = (t478 * t499 - t576) * qJD(3) + ((-qJ(4) * t557 - t556) * t471 + t431) * qJD(2);
t429 = qJD(2) * t439;
t457 = t482 * t542;
t583 = t474 * t482;
t567 = t456 * t583 + t457;
t398 = (-qJ(4) * t563 - t449) * t478 + t567;
t375 = pkin(3) * t524 + t398;
t560 = qJD(2) * t482;
t399 = t456 * t584 + t576 + (qJ(4) * t560 + t478 * t564) * t471;
t586 = t473 * t399;
t342 = t470 * t375 + t586;
t340 = pkin(10) * t524 + t342;
t463 = t474 * t564;
t412 = qJD(4) + t463 + (-pkin(3) * t560 - t456) * t471;
t354 = -pkin(4) * t438 - pkin(10) * t441 + t412;
t324 = t340 * t481 + t354 * t477;
t562 = qJD(2) * t472;
t533 = qJD(1) * t562;
t518 = t483 * t533;
t544 = -qJD(3) * t457 - t456 * t538 - t482 * t518;
t350 = -t449 * t558 + (-t536 * t584 + (-qJ(4) * t558 + t555) * t471) * qJD(2) - t544;
t330 = t473 * t350 + t470 * t606;
t548 = qJD(2) * qJD(3);
t532 = t471 * t548;
t516 = t482 * t532;
t517 = t478 * t532;
t430 = -t470 * t517 + t473 * t516;
t519 = t479 * t533;
t435 = pkin(3) * t517 + t471 * t519;
t361 = pkin(4) * t429 - pkin(10) * t430 + t435;
t529 = t330 * t477 - t481 * t361;
t604 = -qJD(5) * t324 - t529;
t312 = -pkin(5) * t429 - t604;
t492 = -t481 * t441 - t477 * t524;
t605 = t405 * (-pkin(5) * t492 + pkin(11) * t405) + t312;
t370 = -qJD(5) * t492 + t477 * t430;
t484 = qJD(2) ^ 2;
t369 = qJD(5) * t459 + t481 * t430 - t441 * t553;
t476 = sin(qJ(6));
t480 = cos(qJ(6));
t549 = qJD(6) * t480;
t545 = t480 * t369 + t476 * t429 + t434 * t549;
t550 = qJD(6) * t476;
t333 = t492 * t550 + t545;
t600 = t333 * t476;
t592 = t492 * t476;
t371 = -t480 * t434 - t592;
t599 = t371 * t405;
t373 = t434 * t476 - t480 * t492;
t598 = t373 * t405;
t597 = t373 * t438;
t526 = t405 * t480;
t596 = t406 * t434;
t595 = t406 * t441;
t594 = t492 * t434;
t593 = t492 * t441;
t433 = -t456 * t471 + t463;
t591 = t433 * t471;
t590 = t434 * t477;
t589 = t438 * t481;
t381 = t470 * t399;
t582 = t476 * t370;
t577 = t480 * t370;
t574 = -pkin(5) * t439 + qJD(5) * t607 - t609;
t347 = t398 * t470 + t586;
t573 = t347 - t434 * (pkin(5) * t477 - pkin(11) * t481);
t348 = t398 * t473 - t381;
t394 = pkin(3) * t541 + pkin(4) * t441 - pkin(10) * t438;
t572 = t481 * t348 + t477 * t394;
t568 = t481 * t429 + t438 * t590;
t554 = qJD(5) * t476;
t552 = qJD(5) * t480;
t466 = -pkin(3) * t473 - pkin(4);
t539 = t471 * t557;
t496 = t481 * t330 - t340 * t553 + t354 * t551 + t477 * t361;
t311 = pkin(11) * t429 + t496;
t329 = t350 * t470 - t473 * t606;
t316 = pkin(5) * t370 - pkin(11) * t369 + t329;
t530 = -t311 * t476 + t480 * t316;
t528 = t369 * t476 - t480 * t429;
t341 = t473 * t375 - t381;
t389 = t423 * t473 - t470 * t436;
t527 = t434 * t481;
t448 = -pkin(5) * t481 - pkin(11) * t477 + t466;
t525 = pkin(11) * t441 - qJD(6) * t448 + t572;
t520 = t471 * t479 * t562;
t514 = -t449 + t536;
t376 = -pkin(4) * t474 - t389;
t415 = t444 * t477 - t481 * t474;
t416 = t444 * t481 + t474 * t477;
t345 = pkin(5) * t415 - pkin(11) * t416 + t376;
t513 = -pkin(11) * t439 - qJD(6) * t345 + t608;
t338 = pkin(11) * t443 + t607;
t387 = -qJD(5) * t415 + t440 * t481;
t388 = qJD(5) * t416 + t440 * t477;
t512 = -pkin(5) * t388 + pkin(11) * t387 + qJD(6) * t338 - t571;
t511 = t311 * t480 + t316 * t476;
t320 = pkin(11) * t434 + t324;
t339 = -pkin(4) * t524 - t341;
t325 = t406 * pkin(5) + pkin(11) * t492 + t339;
t314 = t320 * t480 + t325 * t476;
t510 = t320 * t476 - t325 * t480;
t323 = -t340 * t477 + t354 * t481;
t501 = t474 * t579 + t578;
t413 = t472 * t501 + t475 * t588;
t503 = t474 * t575 - t580;
t487 = t472 * t503 + t475 * t587;
t367 = t473 * t413 + t470 * t487;
t447 = -t471 * t472 * t483 + t474 * t475;
t352 = t367 * t481 + t447 * t477;
t366 = t413 * t470 - t473 * t487;
t509 = t352 * t480 + t366 * t476;
t508 = -t352 * t476 + t366 * t480;
t351 = t367 * t477 - t447 * t481;
t505 = -t377 * t477 + t400 * t481;
t397 = t416 * t480 + t443 * t476;
t396 = t416 * t476 - t480 * t443;
t498 = -t405 * t549 - t582;
t497 = t405 * t550 - t577;
t465 = pkin(3) * t470 + pkin(10);
t494 = t339 * t434 - t465 * t429;
t402 = t441 * t476 + t480 * t589;
t493 = -t477 * t550 + t480 * t551 - t402;
t491 = -qJD(3) * t449 - t474 * t519;
t490 = qJD(5) * t371 + t498;
t319 = -pkin(5) * t434 - t323;
t489 = -pkin(11) * t370 + (t319 + t323) * t405;
t488 = t474 * t499 + t591;
t401 = -t480 * t441 + t476 * t589;
t386 = t475 * t539 + (qJD(2) * t500 + qJD(3) * t503) * t472;
t385 = -t475 * t540 + (qJD(2) * t502 - qJD(3) * t501) * t472;
t364 = t373 * t553;
t344 = t385 * t470 + t386 * t473;
t343 = -t473 * t385 + t386 * t470;
t337 = -pkin(5) * t443 - t505;
t336 = qJD(6) * t397 + t387 * t476 - t480 * t439;
t335 = -qJD(6) * t396 + t387 * t480 + t439 * t476;
t334 = qJD(6) * t373 + t528;
t326 = -pkin(5) * t441 + t348 * t477 - t394 * t481;
t322 = qJD(5) * t352 + t344 * t477 - t481 * t520;
t321 = -qJD(5) * t351 + t344 * t481 + t477 * t520;
t310 = -qJD(6) * t314 + t530;
t309 = -qJD(6) * t510 + t511;
t1 = [(t385 * t524 + t447 * t517) * MDP(10) + (-t386 * t524 + t447 * t516) * MDP(11) + (t343 * t441 + t344 * t438 + t366 * t430 - t367 * t429) * MDP(12) + (t329 * t366 + t330 * t367 - t341 * t343 + t342 * t344 + t435 * t447) * MDP(13) + (-t322 * t434 + t343 * t406 - t351 * t429 + t366 * t370) * MDP(19) + (-t321 * t434 - t343 * t492 - t352 * t429 + t366 * t369) * MDP(20) + ((-qJD(6) * t509 - t321 * t476 + t343 * t480) * t405 + t508 * t370 + t322 * t371 + t351 * t334) * MDP(26) + (-(qJD(6) * t508 + t321 * t480 + t343 * t476) * t405 - t509 * t370 + t322 * t373 + t351 * t333) * MDP(27) + (-t484 * t483 * MDP(4) + (MDP(13) * t412 * t563 + (-MDP(3) + (-MDP(10) * t482 + MDP(11) * t478) * t467) * t484) * t479) * t472; (-t431 * t524 + (-pkin(9) * t524 * t559 + t474 * t491) * t482 + (-t474 * t518 + t488 * qJD(3) + (-qJD(3) * t474 + (-t474 ^ 2 - t467) * qJD(2)) * qJD(3) * pkin(2)) * t478) * MDP(10) + (-(t478 * t491 - t544) * t474 + (-t467 * t602 + t591) * t557 + (pkin(9) * t540 - t615) * t524) * MDP(11) + (t329 * t444 - t330 * t443 - t341 * t440 - t342 * t439 - t389 * t430 - t390 * t429 + t438 * t570 + t441 * t571) * MDP(12) + (-t329 * t389 + t330 * t390 + t570 * t342 - t571 * t341 + (t435 * t543 + (pkin(3) * t558 - t536) * t412) * t471) * MDP(13) + (t369 * t416 - t387 * t492) * MDP(14) + (-t369 * t415 - t370 * t416 - t387 * t406 + t388 * t492) * MDP(15) + (t369 * t443 + t387 * t434 + t416 * t429 - t439 * t492) * MDP(16) + (-t370 * t443 - t388 * t434 - t406 * t439 - t415 * t429) * MDP(17) + (t429 * t443 + t434 * t439) * MDP(18) + (t505 * t429 - t529 * t443 + t323 * t439 + t376 * t370 + t329 * t415 + t339 * t388 + t609 * t434 + t571 * t406 + (-t324 * t443 - t434 * t607) * qJD(5)) * MDP(19) + (-t324 * t439 + t329 * t416 + t339 * t387 + t376 * t369 - t429 * t607 + t434 * t608 - t496 * t443 - t492 * t571) * MDP(20) + (t333 * t397 + t335 * t373) * MDP(21) + (-t333 * t396 - t334 * t397 - t335 * t371 - t336 * t373) * MDP(22) + (t333 * t415 + t335 * t405 + t370 * t397 + t373 * t388) * MDP(23) + (-t334 * t415 - t336 * t405 - t370 * t396 - t371 * t388) * MDP(24) + (t370 * t415 + t388 * t405) * MDP(25) + ((-t338 * t476 + t345 * t480) * t370 + t310 * t415 - t510 * t388 + t337 * t334 + t312 * t396 + t319 * t336 + (t476 * t513 - t480 * t512) * t405 + t574 * t371) * MDP(26) + (-(t338 * t480 + t345 * t476) * t370 - t309 * t415 - t314 * t388 + t337 * t333 + t312 * t397 + t319 * t335 + (t476 * t512 + t480 * t513) * t405 + t574 * t373) * MDP(27) + 0.2e1 * t611 * t548 + (MDP(7) * t539 - MDP(8) * t540) * (qJD(3) + 0.2e1 * t561); (-t514 * t583 + (-t488 - t535) * t478) * qJD(2) * MDP(10) + ((-t433 * t587 + (t478 * t514 + t567) * t474) * qJD(2) + t567 * qJD(3) + t544) * MDP(11) + ((t342 - t347) * t441 + (t341 - t348) * t438 + (-t429 * t470 - t430 * t473) * pkin(3)) * MDP(12) + (t341 * t347 - t342 * t348 + (-t329 * t473 + t330 * t470 - t412 * t541) * pkin(3)) * MDP(13) + (t369 * t477 - t492 * t527) * MDP(14) + ((t369 - t596) * t481 + (-t370 + t594) * t477) * MDP(15) + (t477 * t429 + t434 * t527 + t593) * MDP(16) + (-t434 * t553 + t568 + t595) * MDP(17) - t434 * t441 * MDP(18) + (-t323 * t441 - t347 * t406 + t466 * t370 + (-t329 + (-qJD(5) * t465 - t394) * t434) * t481 + (t348 * t434 + t494) * t477) * MDP(19) + (t324 * t441 + t329 * t477 + t347 * t492 + t466 * t369 + (t465 * t553 + t572) * t434 + t494 * t481) * MDP(20) + (t333 * t477 * t480 + t373 * t493) * MDP(21) + (t371 * t402 + t373 * t401 + (-t371 * t480 - t373 * t476) * t551 + (-t600 - t334 * t480 + (t371 * t476 - t373 * t480) * qJD(6)) * t477) * MDP(22) + (-t333 * t481 + t364 + (t577 - t597) * t477 + t493 * t405) * MDP(23) + (t334 * t481 + (-t476 * t551 + t401) * t405 + (-t371 * t434 + t498) * t477) * MDP(24) + (-t370 * t481 + t405 * t590) * MDP(25) + (t448 * t577 - t319 * t401 - t326 * t371 + (t476 * t525 - t480 * t573) * t405 + (t319 * t554 + t465 * t490 - t310) * t481 + (t319 * t549 + t312 * t476 + t510 * t438 + t465 * t334 + (t405 * t465 * t476 - t510) * qJD(5)) * t477) * MDP(26) + (-t448 * t582 - t319 * t402 - t326 * t373 + (t476 * t573 + t480 * t525) * t405 + (t319 * t552 + t309 + (qJD(5) * t373 + t497) * t465) * t481 + (-t319 * t550 + t312 * t480 + t314 * t438 + t465 * t333 + (t465 * t526 - t314) * qJD(5)) * t477) * MDP(27) + ((-MDP(7) * t482 + MDP(8) * t478) * t471 * t474 - t611) * t484; (-t438 ^ 2 - t441 ^ 2) * MDP(12) + (t341 * t441 - t342 * t438 + t435) * MDP(13) + (t568 - t595) * MDP(19) + MDP(20) * t593 + t401 * t405 * MDP(26) + (t402 * t405 + t364) * MDP(27) + ((-t405 * t554 - t334) * MDP(26) + (-t405 * t552 - t333) * MDP(27) - t434 ^ 2 * MDP(20)) * t481 + (-qJD(5) * t434 * MDP(19) - t429 * MDP(20) + (-t371 * t438 + t490) * MDP(26) + (t497 - t597) * MDP(27)) * t477; -t406 ^ 2 * MDP(15) + (t369 + t596) * MDP(16) + (-t370 - t594) * MDP(17) + t429 * MDP(18) + (t324 * t434 + t604) * MDP(19) + (t323 * t434 + t339 * t406 - t496) * MDP(20) + (t373 * t526 + t600) * MDP(21) + ((t333 - t599) * t480 + (-t334 - t598) * t476) * MDP(22) + (t405 * t526 + t582) * MDP(23) + (-t405 ^ 2 * t476 + t577) * MDP(24) + (-pkin(5) * t334 - t324 * t371 + t489 * t476 - t480 * t605) * MDP(26) + (-pkin(5) * t333 - t324 * t373 + t476 * t605 + t489 * t480) * MDP(27) - (MDP(14) * t406 - MDP(15) * t492 - t339 * MDP(19) - MDP(23) * t373 + t371 * MDP(24) - t405 * MDP(25) + MDP(26) * t510 + MDP(27) * t314) * t492; t373 * t371 * MDP(21) + (-t371 ^ 2 + t373 ^ 2) * MDP(22) + (t545 + t599) * MDP(23) + (-t528 + t598) * MDP(24) + t370 * MDP(25) + (t314 * t405 - t319 * t373 + t530) * MDP(26) + (t319 * t371 - t405 * t510 - t511) * MDP(27) + (MDP(23) * t592 - MDP(24) * t373 - MDP(26) * t314 + MDP(27) * t510) * qJD(6);];
tauc  = t1;
