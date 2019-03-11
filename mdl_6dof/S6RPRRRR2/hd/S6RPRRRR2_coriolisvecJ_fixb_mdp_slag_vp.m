% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRRR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:58:41
% EndTime: 2019-03-09 06:58:50
% DurationCPUTime: 5.34s
% Computational Cost: add. (5501->418), mult. (12539->561), div. (0->0), fcn. (9208->10), ass. (0->195)
t602 = qJD(5) + qJD(6);
t497 = sin(qJ(4));
t501 = cos(qJ(3));
t600 = cos(qJ(4));
t544 = qJD(1) * t600;
t498 = sin(qJ(3));
t565 = qJD(1) * t498;
t619 = -t497 * t565 + t501 * t544;
t620 = -t619 + t602;
t481 = sin(pkin(11)) * pkin(1) + pkin(7);
t598 = pkin(8) + t481;
t539 = t598 * qJD(1);
t443 = qJD(2) * t498 + t501 * t539;
t434 = t497 * t443;
t442 = t501 * qJD(2) - t539 * t498;
t391 = t442 * t600 - t434;
t543 = qJD(4) * t600;
t614 = -pkin(3) * t543 + t391;
t581 = t497 * t501;
t455 = -qJD(1) * t581 - t498 * t544;
t418 = -pkin(4) * t455 - pkin(9) * t619;
t405 = pkin(3) * t565 + t418;
t496 = sin(qJ(5));
t500 = cos(qJ(5));
t618 = -t500 * t405 + t614 * t496;
t490 = qJD(3) + qJD(4);
t437 = -t455 * t496 - t500 * t490;
t499 = cos(qJ(6));
t495 = sin(qJ(6));
t517 = t455 * t500 - t490 * t496;
t588 = t517 * t495;
t380 = t499 * t437 - t588;
t452 = qJD(5) - t619;
t447 = qJD(6) + t452;
t617 = t380 * t447;
t518 = t437 * t495 + t499 * t517;
t616 = t447 * t518;
t463 = t495 * t500 + t496 * t499;
t615 = t620 * t463;
t461 = t495 * t496 - t499 * t500;
t571 = t620 * t461;
t435 = t600 * t443;
t596 = qJD(3) * pkin(3);
t436 = t442 + t596;
t386 = t497 * t436 + t435;
t377 = pkin(9) * t490 + t386;
t482 = -cos(pkin(11)) * pkin(1) - pkin(2);
t468 = -pkin(3) * t501 + t482;
t456 = t468 * qJD(1);
t396 = -pkin(4) * t619 + pkin(9) * t455 + t456;
t351 = t377 * t500 + t396 * t496;
t336 = -pkin(10) * t437 + t351;
t560 = qJD(6) * t495;
t333 = t336 * t560;
t385 = t436 * t600 - t434;
t376 = -t490 * pkin(4) - t385;
t354 = t437 * pkin(5) + t376;
t613 = t354 * t380 + t333;
t416 = t619 * t490;
t561 = qJD(5) * t500;
t562 = qJD(5) * t496;
t372 = t500 * t416 + t455 * t562 + t490 * t561;
t464 = t498 * t600 + t581;
t427 = t490 * t464;
t417 = t427 * qJD(1);
t432 = t442 * qJD(3);
t433 = t443 * qJD(3);
t563 = qJD(4) * t497;
t344 = t432 * t600 - t497 * t433 + t436 * t543 - t443 * t563;
t556 = qJD(1) * qJD(3);
t542 = t498 * t556;
t359 = pkin(3) * t542 + pkin(4) * t417 - pkin(9) * t416;
t357 = t500 * t359;
t506 = -qJD(5) * t351 - t344 * t496 + t357;
t313 = pkin(5) * t417 - pkin(10) * t372 + t506;
t373 = -qJD(5) * t517 + t416 * t496;
t510 = t500 * t344 + t496 * t359 - t377 * t562 + t396 * t561;
t315 = -pkin(10) * t373 + t510;
t538 = t499 * t313 - t495 * t315;
t612 = t354 * t518 + t538;
t611 = t417 * MDP(30) + (-t380 ^ 2 + t518 ^ 2) * MDP(27) - t380 * MDP(26) * t518;
t610 = MDP(5) * t501;
t609 = MDP(6) * (t498 ^ 2 - t501 ^ 2);
t411 = t463 * t464;
t390 = t497 * t442 + t435;
t531 = pkin(3) * t563 - t390;
t587 = t619 * t496;
t608 = (t562 - t587) * pkin(5);
t607 = t496 * t405 + t614 * t500;
t457 = t598 * t498;
t458 = t598 * t501;
t606 = -t600 * t457 - t497 * t458;
t605 = -t498 * MDP(10) - t501 * MDP(11);
t514 = -t497 * t498 + t501 * t600;
t426 = t490 * t514;
t579 = t500 * t426;
t604 = -t464 * t562 + t579;
t603 = qJD(1) * t464;
t537 = t372 * t495 + t499 * t373;
t327 = -qJD(6) * t518 + t537;
t601 = -pkin(9) - pkin(10);
t599 = t500 * pkin(5);
t485 = pkin(3) * t497 + pkin(9);
t597 = -pkin(10) - t485;
t350 = -t377 * t496 + t500 * t396;
t335 = pkin(10) * t517 + t350;
t330 = pkin(5) * t452 + t335;
t595 = t330 * t499;
t594 = t336 * t499;
t593 = t372 * t496;
t591 = t417 * t500;
t590 = t437 * t452;
t589 = t517 * t452;
t586 = t619 * t500;
t585 = t464 * t496;
t584 = t464 * t500;
t583 = t496 * t417;
t582 = t496 * t426;
t502 = qJD(3) ^ 2;
t580 = t498 * t502;
t415 = -t497 * t457 + t458 * t600;
t406 = t500 * t415;
t578 = t501 * t502;
t559 = qJD(6) * t499;
t549 = t499 * t372 - t495 * t373 - t437 * t559;
t326 = t517 * t560 + t549;
t577 = -t326 * t514 - t427 * t518;
t332 = -t560 * t585 + (t584 * t602 + t582) * t499 + t604 * t495;
t576 = -t332 * t447 - t411 * t417;
t575 = -t372 * t514 - t427 * t517;
t574 = t500 * t385 + t496 * t418;
t408 = -pkin(4) * t514 - pkin(9) * t464 + t468;
t570 = t496 * t408 + t406;
t569 = t608 + t531;
t471 = qJD(1) * t482;
t555 = pkin(10) * t587;
t554 = t498 * t596;
t551 = t464 * t583;
t550 = t417 * t584;
t548 = qJD(5) * t601;
t375 = t376 * t561;
t541 = qJD(3) * t598;
t540 = qJD(5) * t597;
t536 = -t385 * t496 + t500 * t418;
t535 = t452 * t500;
t534 = qJD(6) * t330 + t315;
t345 = t497 * t432 + t600 * t433 + t436 * t563 + t443 * t543;
t486 = -pkin(3) * t600 - pkin(4);
t532 = -t455 * pkin(5) - pkin(10) * t586;
t530 = -t386 + t608;
t529 = t345 * t496 - t351 * t455 + t375;
t489 = t500 * pkin(10);
t460 = t485 * t500 + t489;
t528 = qJD(6) * t460 - t500 * t540 + t532 - t618;
t473 = pkin(9) * t500 + t489;
t527 = qJD(6) * t473 - t500 * t548 + t532 + t536;
t459 = t597 * t496;
t526 = -qJD(6) * t459 - t496 * t540 - t555 + t607;
t472 = t601 * t496;
t525 = -qJD(6) * t472 - t496 * t548 - t555 + t574;
t522 = t327 * t514 - t380 * t427;
t318 = t330 * t495 + t594;
t331 = -t411 * t602 - t461 * t426;
t412 = t461 * t464;
t521 = -t331 * t447 + t412 * t417;
t520 = t373 * t514 - t427 * t437;
t519 = -t376 * t619 - t417 * t485;
t516 = 0.2e1 * qJD(3) * t471;
t515 = -t345 * t500 + t350 * t455 + t376 * t562;
t513 = t464 * t561 + t582;
t511 = t455 * t456 - t345;
t449 = t498 * t541;
t450 = t501 * t541;
t361 = qJD(4) * t606 - t600 * t449 - t497 * t450;
t371 = pkin(4) * t427 - pkin(9) * t426 + t554;
t509 = t500 * t361 + t496 * t371 + t408 * t561 - t415 * t562;
t317 = -t336 * t495 + t595;
t328 = pkin(5) * t373 + t345;
t508 = t317 * t455 + t328 * t461 + t615 * t354;
t507 = -t318 * t455 + t328 * t463 - t571 * t354;
t505 = -t456 * t619 - t344;
t362 = qJD(4) * t415 - t497 * t449 + t450 * t600;
t504 = (-t326 * t461 - t327 * t463 + t380 * t571 + t518 * t615) * MDP(27) + (t326 * t463 + t518 * t571) * MDP(26) + ((t372 - t590) * t500 + (-t373 + t589) * t496) * MDP(20) + (t417 * t463 - t447 * t571 - t455 * t518) * MDP(28) + (-t380 * t455 - t417 * t461 - t447 * t615) * MDP(29) + (-t517 * t535 + t593) * MDP(19) + (-t452 ^ 2 * t496 - t437 * t455 + t591) * MDP(22) + (t452 * t535 - t455 * t517 + t583) * MDP(21) + t416 * MDP(14) + (t455 ^ 2 - t619 ^ 2) * MDP(13) + (MDP(12) * t619 + MDP(23) * t452 + MDP(30) * t447) * t455 + (-t619 * MDP(14) + (-t455 - t603) * MDP(15)) * t490;
t487 = -pkin(4) - t599;
t469 = t486 - t599;
t404 = t500 * t408;
t397 = t417 * t514;
t393 = pkin(5) * t585 - t606;
t367 = t500 * t371;
t353 = -pkin(10) * t585 + t570;
t352 = -pkin(5) * t514 - pkin(10) * t584 - t415 * t496 + t404;
t341 = pkin(5) * t513 + t362;
t320 = -pkin(10) * t513 + t509;
t316 = -pkin(10) * t579 + pkin(5) * t427 - t361 * t496 + t367 + (-t406 + (pkin(10) * t464 - t408) * t496) * qJD(5);
t1 = [((-t415 * t561 + t367) * t452 + t404 * t417 - (-t377 * t561 + t357) * t514 + t350 * t427 + t362 * t437 - t606 * t373 + t464 * t375 + ((-qJD(5) * t408 - t361) * t452 - t415 * t417 - (-qJD(5) * t396 - t344) * t514 + t345 * t464 + t376 * t426) * t496) * MDP(24) + (t345 * t584 - t351 * t427 - t362 * t517 - t372 * t606 + t376 * t604 - t417 * t570 - t452 * t509 + t510 * t514) * MDP(25) + (t452 * t604 + t550 + t575) * MDP(21) - 0.2e1 * t556 * t609 + 0.2e1 * t542 * t610 + (t416 * t468 + t426 * t456 + (-t455 + t603) * t554) * MDP(18) + ((-t437 * t500 + t496 * t517) * t426 + (-t593 - t373 * t500 + (t437 * t496 + t500 * t517) * qJD(5)) * t464) * MDP(20) + (-t318 * t427 + t393 * t326 - t328 * t412 + t354 * t331 - t333 * t514 - t341 * t518 + (-(-qJD(6) * t353 + t316) * t447 - t352 * t417 + t313 * t514) * t495 + (-(qJD(6) * t352 + t320) * t447 - t353 * t417 + t534 * t514) * t499) * MDP(32) + ((t316 * t499 - t320 * t495) * t447 + (t352 * t499 - t353 * t495) * t417 - t538 * t514 + t317 * t427 + t341 * t380 + t393 * t327 + t328 * t411 + t354 * t332 + ((-t352 * t495 - t353 * t499) * t447 + t318 * t514) * qJD(6)) * MDP(31) + (t427 * t452 - t397) * MDP(23) + (t427 * t447 - t397) * MDP(30) + (-t326 * t411 + t327 * t412 - t331 * t380 + t332 * t518) * MDP(27) + (-t326 * t412 - t331 * t518) * MDP(26) + (t372 * t584 - t517 * t604) * MDP(19) + (-t481 * t578 + t498 * t516) * MDP(10) - MDP(8) * t580 + (t481 * t580 + t501 * t516) * MDP(11) + (t522 + t576) * MDP(29) + (-t521 + t577) * MDP(28) + (-t452 * t513 + t520 - t551) * MDP(22) + (MDP(14) * t426 - MDP(15) * t427 - MDP(17) * t362 - MDP(18) * t361) * t490 + (t417 * t468 + t427 * t456 + (-qJD(1) * t514 - t619) * t554) * MDP(17) + (t416 * t514 - t417 * t464 + t426 * t619 + t427 * t455) * MDP(13) + (t416 * t464 - t426 * t455) * MDP(12) + MDP(7) * t578; (-t520 - t551) * MDP(24) + (-t550 + t575) * MDP(25) + (-t522 + t576) * MDP(31) + (t521 + t577) * MDP(32) + t605 * t502 + (-t427 * MDP(17) - t426 * MDP(18)) * t490 + (-MDP(24) * t513 - MDP(25) * t604) * t452; (t390 * t490 + (-t490 * t563 + t565 * t619) * pkin(3) + t511) * MDP(17) + t504 + (-(t459 * t495 + t460 * t499) * t417 + t469 * t326 + (t495 * t528 + t499 * t526) * t447 - t569 * t518 + t507) * MDP(32) + (t486 * t373 + t519 * t496 + t531 * t437 + (-t485 * t561 + t618) * t452 + t515) * MDP(24) + (t486 * t372 + t519 * t500 - t531 * t517 + (t485 * t562 + t607) * t452 + t529) * MDP(25) + ((t459 * t499 - t460 * t495) * t417 + t469 * t327 + (t495 * t526 - t499 * t528) * t447 + t569 * t380 + t508) * MDP(31) + (t391 * t490 + (t455 * t565 - t490 * t543) * pkin(3) + t505) * MDP(18) + (t605 * t471 + (-t498 * t610 + t609) * qJD(1)) * qJD(1); t504 + (t385 * t490 + t505) * MDP(18) + (t386 * t490 + t511) * MDP(17) + ((t472 * t499 - t473 * t495) * t417 + t487 * t327 + (t495 * t525 - t499 * t527) * t447 + t530 * t380 + t508) * MDP(31) + (-(t472 * t495 + t473 * t499) * t417 + t487 * t326 + (t495 * t527 + t499 * t525) * t447 - t530 * t518 + t507) * MDP(32) + (-pkin(4) * t373 - t536 * t452 - t386 * t437 - t376 * t587 + (-t452 * t561 - t583) * pkin(9) + t515) * MDP(24) + (-pkin(4) * t372 + t574 * t452 + t386 * t517 - t376 * t586 + (t452 * t562 - t591) * pkin(9) + t529) * MDP(25); -t517 * t437 * MDP(19) + (-t437 ^ 2 + t517 ^ 2) * MDP(20) + (t372 + t590) * MDP(21) + (-t373 - t589) * MDP(22) + t417 * MDP(23) + (t351 * t452 + t376 * t517 + t506) * MDP(24) + (t350 * t452 + t376 * t437 - t510) * MDP(25) + (t326 + t617) * MDP(28) + (-t327 - t616) * MDP(29) + (-(-t335 * t495 - t594) * t447 - t318 * qJD(6) + (t380 * t517 + t499 * t417 - t447 * t560) * pkin(5) + t612) * MDP(31) + ((-t336 * t447 - t313) * t495 + (t335 * t447 - t534) * t499 + (-t495 * t417 - t447 * t559 - t517 * t518) * pkin(5) + t613) * MDP(32) + t611; (t549 + t617) * MDP(28) + (-t537 - t616) * MDP(29) + (t318 * t447 + t612) * MDP(31) + (-t495 * t313 - t499 * t315 + t317 * t447 + t613) * MDP(32) + (MDP(28) * t588 + MDP(29) * t518 - MDP(31) * t318 - MDP(32) * t595) * qJD(6) + t611;];
tauc  = t1;
