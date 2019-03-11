% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RRPRRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:41:24
% EndTime: 2019-03-09 11:41:34
% DurationCPUTime: 5.62s
% Computational Cost: add. (7643->409), mult. (19605->532), div. (0->0), fcn. (14954->8), ass. (0->200)
t506 = cos(qJ(5));
t501 = sin(pkin(10));
t502 = cos(pkin(10));
t505 = sin(qJ(2));
t508 = cos(qJ(2));
t478 = t501 * t508 + t502 * t505;
t468 = t478 * qJD(2);
t456 = qJD(1) * t468;
t546 = qJD(1) * qJD(2);
t538 = t508 * t546;
t539 = t505 * t546;
t457 = -t501 * t539 + t502 * t538;
t477 = -t501 * t505 + t502 * t508;
t467 = t477 * qJD(1);
t469 = t478 * qJD(1);
t504 = sin(qJ(4));
t507 = cos(qJ(4));
t553 = qJD(4) * t507;
t554 = qJD(4) * t504;
t515 = -t504 * t456 + t507 * t457 + t467 * t553 - t469 * t554;
t522 = t467 * t504 + t507 * t469;
t545 = qJD(2) + qJD(4);
t525 = t506 * t545;
t503 = sin(qJ(5));
t552 = qJD(5) * t503;
t350 = -qJD(5) * t525 - t506 * t515 + t522 * t552;
t348 = t350 * t503;
t349 = t506 * t350;
t405 = t503 * t522 - t525;
t408 = t503 * t545 + t506 * t522;
t422 = t507 * t467 - t469 * t504;
t531 = t507 * t456 + t504 * t457;
t549 = t522 * qJD(2);
t551 = qJD(5) * t506;
t599 = t522 * qJD(4);
t378 = t531 + t599;
t375 = t503 * t378;
t616 = qJD(5) - t422;
t563 = t551 * t616 + t375;
t513 = t503 * t515;
t351 = qJD(5) * t408 + t513;
t567 = -t503 * t351 - t405 * t551;
t583 = t408 * t522;
t624 = t616 * t503;
t606 = t408 * t624;
t619 = t422 * t506;
t628 = (t549 - t531) * MDP(16) - t422 ^ 2 * MDP(14) + (-t422 * MDP(13) + MDP(14) * t522 - t616 * MDP(24)) * t522 + (-t348 + (t551 - t619) * t408) * MDP(20) + (-t616 * t619 + t563 - t583) * MDP(22) + (t405 * t619 - t349 + t567 - t606) * MDP(21);
t627 = pkin(5) * t624;
t589 = -qJ(3) - pkin(7);
t537 = qJD(2) * t589;
t462 = qJD(3) * t508 + t505 * t537;
t444 = t462 * qJD(1);
t463 = -qJD(3) * t505 + t508 * t537;
t445 = t463 * qJD(1);
t404 = -t444 * t501 + t502 * t445;
t387 = -pkin(8) * t457 + t404;
t407 = t502 * t444 + t501 * t445;
t388 = -pkin(8) * t456 + t407;
t488 = t589 * t508;
t483 = qJD(1) * t488;
t472 = t501 * t483;
t486 = t589 * t505;
t482 = qJD(1) * t486;
t587 = qJD(2) * pkin(2);
t476 = t482 + t587;
t424 = t502 * t476 + t472;
t590 = pkin(8) * t469;
t399 = qJD(2) * pkin(3) + t424 - t590;
t574 = t502 * t483;
t425 = t501 * t476 - t574;
t591 = pkin(8) * t467;
t402 = t425 + t591;
t329 = t507 * (qJD(4) * t399 + t388) + t504 * t387 - t402 * t554;
t492 = pkin(2) * t539;
t431 = t456 * pkin(3) + t492;
t342 = t378 * pkin(4) - pkin(9) * t515 + t431;
t358 = t504 * t399 + t507 * t402;
t355 = pkin(9) * t545 + t358;
t542 = -pkin(2) * t508 - pkin(1);
t524 = t542 * qJD(1);
t484 = qJD(3) + t524;
t430 = -pkin(3) * t467 + t484;
t362 = -pkin(4) * t422 - pkin(9) * t522 + t430;
t514 = t506 * t329 + t503 * t342 - t355 * t552 + t362 * t551;
t318 = -qJ(6) * t351 - qJD(6) * t405 + t514;
t335 = t355 * t506 + t362 * t503;
t341 = t506 * t342;
t511 = -qJD(5) * t335 - t329 * t503 + t341;
t316 = pkin(5) * t378 + qJ(6) * t350 - qJD(6) * t408 + t511;
t326 = -qJ(6) * t405 + t335;
t598 = t326 * t616 + t316;
t626 = t318 * t506 - t503 * t598;
t428 = -t482 * t501 + t574;
t410 = t428 - t591;
t429 = t502 * t482 + t472;
t411 = t429 - t590;
t494 = pkin(2) * t502 + pkin(3);
t593 = pkin(2) * t501;
t603 = t494 * t507 - t504 * t593;
t604 = -t603 * qJD(4) + t410 * t504 + t411 * t507;
t620 = t422 * t503;
t623 = qJ(6) * t620 + t506 * qJD(6);
t618 = t422 * t545;
t382 = pkin(4) * t522 - pkin(9) * t422;
t615 = -t430 * t422 - t329;
t498 = t506 * qJ(6);
t614 = -pkin(5) * t522 + t422 * t498;
t611 = -0.2e1 * t546;
t608 = MDP(4) * t505;
t607 = MDP(5) * (t505 ^ 2 - t508 ^ 2);
t584 = t405 * t522;
t432 = t502 * t486 + t488 * t501;
t416 = -pkin(8) * t478 + t432;
t433 = t501 * t486 - t502 * t488;
t417 = pkin(8) * t477 + t433;
t373 = -t416 * t507 + t504 * t417;
t518 = t494 * t504 + t507 * t593;
t559 = t518 * qJD(4) + t507 * t410 - t504 * t411;
t377 = t506 * t378;
t602 = t552 * t616 - t377;
t555 = qJD(1) * t505;
t439 = pkin(2) * t555 + pkin(3) * t469;
t367 = t382 + t439;
t600 = t503 * t367 + t604 * t506;
t330 = -t507 * t387 + t504 * t388 + t399 * t554 + t402 * t553;
t597 = -t430 * t522 - t330;
t334 = -t355 * t503 + t506 * t362;
t357 = t507 * t399 - t504 * t402;
t354 = -pkin(4) * t545 - t357;
t585 = t330 * t506;
t596 = -t334 * t522 + t354 * t552 - t585;
t353 = t354 * t551;
t595 = t330 * t503 + t335 * t522 + t353;
t594 = t408 ^ 2;
t592 = pkin(5) * t506;
t588 = -qJ(6) - pkin(9);
t325 = -qJ(6) * t408 + t334;
t322 = pkin(5) * t616 + t325;
t586 = t322 * t506;
t427 = t477 * t504 + t478 * t507;
t577 = t427 * t503;
t509 = qJD(2) ^ 2;
t573 = t505 * t509;
t374 = t416 * t504 + t417 * t507;
t371 = t506 * t374;
t572 = t508 * t509;
t510 = qJD(1) ^ 2;
t571 = t508 * t510;
t461 = pkin(9) + t518;
t570 = -qJ(6) - t461;
t569 = t322 - t325;
t566 = t506 * t357 + t503 * t382;
t447 = -pkin(3) * t477 + t542;
t521 = t507 * t477 - t478 * t504;
t372 = -pkin(4) * t521 - pkin(9) * t427 + t447;
t564 = t503 * t372 + t371;
t530 = qJD(5) * t570;
t562 = t503 * t530 - t600 + t623;
t364 = t506 * t367;
t561 = t506 * t530 - t364 + (-qJD(6) + t604) * t503 + t614;
t415 = t502 * t462 + t501 * t463;
t536 = qJD(5) * t588;
t558 = t503 * t536 - t566 + t623;
t533 = -t357 * t503 + t506 * t382;
t557 = -qJD(6) * t503 + t506 * t536 - t533 + t614;
t496 = t505 * t587;
t414 = -t462 * t501 + t502 * t463;
t471 = t477 * qJD(2);
t390 = -pkin(8) * t471 + t414;
t391 = -pkin(8) * t468 + t415;
t338 = -qJD(4) * t373 + t390 * t504 + t391 * t507;
t383 = qJD(4) * t521 - t468 * t504 + t471 * t507;
t384 = qJD(4) * t427 + t507 * t468 + t471 * t504;
t440 = pkin(3) * t468 + t496;
t345 = pkin(4) * t384 - pkin(9) * t383 + t440;
t544 = t506 * t338 + t503 * t345 + t372 * t551;
t540 = t427 * t551;
t535 = pkin(1) * t611;
t528 = t616 * t506;
t460 = -pkin(4) - t603;
t523 = -t354 * t422 - t378 * t461;
t520 = -qJ(6) * t383 - qJD(6) * t427;
t519 = t616 * t620 - t602;
t517 = t383 * t503 + t540;
t516 = t383 * t506 - t427 * t552;
t321 = pkin(5) * t351 + t330;
t339 = qJD(4) * t374 - t390 * t507 + t391 * t504;
t487 = pkin(9) * t506 + t498;
t485 = t588 * t503;
t435 = t461 * t506 + t498;
t434 = t570 * t503;
t403 = t405 ^ 2;
t370 = t506 * t372;
t346 = t405 * pkin(5) + qJD(6) + t354;
t344 = t506 * t345;
t336 = -qJ(6) * t577 + t564;
t332 = -pkin(5) * t521 - t374 * t503 - t427 * t498 + t370;
t320 = -qJ(6) * t540 + (-qJD(5) * t374 + t520) * t503 + t544;
t319 = pkin(5) * t384 - t338 * t503 + t344 + t520 * t506 + (-t371 + (qJ(6) * t427 - t372) * t503) * qJD(5);
t1 = [0.2e1 * t538 * t608 + t607 * t611 + MDP(6) * t572 - MDP(7) * t573 + (-pkin(7) * t572 + t505 * t535) * MDP(9) + (pkin(7) * t573 + t508 * t535) * MDP(10) + (-t404 * t478 + t407 * t477 - t414 * t469 + t415 * t467 - t424 * t471 - t425 * t468 - t432 * t457 - t433 * t456) * MDP(11) + (t404 * t432 + t407 * t433 + t424 * t414 + t425 * t415 + (t484 + t524) * t496) * MDP(12) + (t383 * t522 + t427 * t515) * MDP(13) + (-t427 * t378 + t383 * t422 - t384 * t522 + t515 * t521) * MDP(14) + (t447 * t378 + t430 * t384 - t422 * t440 - t431 * t521) * MDP(18) + (t430 * t383 + t431 * t427 + t440 * t522 + t447 * t515) * MDP(19) + (-t349 * t427 + t408 * t516) * MDP(20) + ((-t405 * t506 - t408 * t503) * t383 + (t348 - t351 * t506 + (t405 * t503 - t408 * t506) * qJD(5)) * t427) * MDP(21) + (t350 * t521 + t377 * t427 + t384 * t408 + t516 * t616) * MDP(22) + (t351 * t521 - t375 * t427 - t384 * t405 - t517 * t616) * MDP(23) + (-t378 * t521 + t384 * t616) * MDP(24) + ((-t374 * t551 + t344) * t616 + t370 * t378 - (-t355 * t551 + t341) * t521 + t334 * t384 + t339 * t405 + t373 * t351 + t427 * t353 + ((-qJD(5) * t372 - t338) * t616 - t374 * t378 - (-qJD(5) * t362 - t329) * t521 + t330 * t427 + t354 * t383) * t503) * MDP(25) + (-(-t374 * t552 + t544) * t616 - t564 * t378 + t514 * t521 - t335 * t384 + t339 * t408 - t373 * t350 + t427 * t585 + t516 * t354) * MDP(26) + (-t319 * t408 - t320 * t405 + t332 * t350 - t336 * t351 + (-t326 * t503 - t586) * t383 + (-t316 * t506 - t318 * t503 + (t322 * t503 - t326 * t506) * qJD(5)) * t427) * MDP(27) + (t318 * t336 + t326 * t320 + t316 * t332 + t322 * t319 + t321 * (pkin(5) * t577 + t373) + t346 * (pkin(5) * t517 + t339)) * MDP(28) + (t383 * MDP(15) - t384 * MDP(16) - t339 * MDP(18) - t338 * MDP(19)) * t545; (t422 * t439 + t597) * MDP(18) + ((t425 + t428) * t469 + (t424 - t429) * t467 + (-t456 * t501 - t457 * t502) * pkin(2)) * MDP(11) + (-t439 * t522 + t615) * MDP(19) + (MDP(9) * t505 * t510 + MDP(10) * t571) * pkin(1) - t571 * t608 + t510 * t607 + (-t460 * t350 + t523 * t506 + t559 * t408 + (t461 * t552 + t600) * t616 + t595) * MDP(26) + (t460 * t351 + t523 * t503 + t559 * t405 + (-t461 * t551 + t503 * t604 - t364) * t616 + t596) * MDP(25) + (-t422 * MDP(15) - MDP(18) * t559 + MDP(19) * t604) * t545 + t515 * MDP(15) + (-t322 * t528 + t350 * t434 - t351 * t435 - t562 * t405 - t561 * t408 + t626) * MDP(27) + (-t424 * t428 - t425 * t429 + (t404 * t502 + t407 * t501 - t484 * t555) * pkin(2)) * MDP(12) + (t519 + t584) * MDP(23) + (t318 * t435 + t316 * t434 + t321 * (t460 - t592) + (t559 + t627) * t346 + t562 * t326 + t561 * t322) * MDP(28) + t628; (-t467 ^ 2 - t469 ^ 2) * MDP(11) + (t424 * t469 - t425 * t467 + t492) * MDP(12) + (t531 + t549 + 0.2e1 * t599) * MDP(18) + (t515 + t618) * MDP(19) + (t519 - t584) * MDP(25) + (-t528 * t616 - t375 - t583) * MDP(26) + ((t405 * t422 + t350) * t506 + t606 + t567) * MDP(27) + (-t346 * t522 + t598 * t506 + (-t322 * t616 + t318) * t503) * MDP(28); (t515 - t618) * MDP(15) + (t358 * t545 + t597) * MDP(18) + (t357 * t545 + t615) * MDP(19) + (-t616 * t624 + t377 + t584) * MDP(23) + (-pkin(4) * t351 - pkin(9) * t563 - t354 * t620 - t358 * t405 - t533 * t616 + t596) * MDP(25) + (pkin(4) * t350 + pkin(9) * t602 - t354 * t619 - t358 * t408 + t566 * t616 + t595) * MDP(26) + (t350 * t485 - t351 * t487 - t558 * t405 - t557 * t408 - t616 * t586 + t626) * MDP(27) + (t318 * t487 + t316 * t485 + t321 * (-pkin(4) - t592) + (-t358 + t627) * t346 + t558 * t326 + t557 * t322) * MDP(28) + t628; t408 * t405 * MDP(20) + (-t403 + t594) * MDP(21) + (t405 * t616 - t350) * MDP(22) + (-t513 + (-qJD(5) + t616) * t408) * MDP(23) + t378 * MDP(24) + (t335 * t616 - t354 * t408 + t511) * MDP(25) + (t334 * t616 + t354 * t405 - t514) * MDP(26) + (pkin(5) * t350 - t405 * t569) * MDP(27) + (t569 * t326 + (-t346 * t408 + t316) * pkin(5)) * MDP(28); (-t403 - t594) * MDP(27) + (t322 * t408 + t326 * t405 + t321) * MDP(28);];
tauc  = t1;
