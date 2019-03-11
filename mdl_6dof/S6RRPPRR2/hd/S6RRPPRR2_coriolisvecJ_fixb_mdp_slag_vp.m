% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPPRR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:52:49
% EndTime: 2019-03-09 08:53:04
% DurationCPUTime: 8.56s
% Computational Cost: add. (7201->482), mult. (18635->651), div. (0->0), fcn. (14776->10), ass. (0->207)
t536 = cos(qJ(2));
t609 = cos(pkin(10));
t564 = t609 * t536;
t518 = qJD(1) * t564;
t529 = sin(pkin(10));
t533 = sin(qJ(2));
t578 = qJD(1) * t533;
t492 = t529 * t578 - t518;
t487 = qJD(5) + t492;
t508 = t529 * t536 + t533 * t609;
t495 = t508 * qJD(1);
t528 = sin(pkin(11));
t530 = cos(pkin(11));
t463 = qJD(2) * t528 + t495 * t530;
t464 = t530 * qJD(2) - t495 * t528;
t532 = sin(qJ(5));
t535 = cos(qJ(5));
t409 = t463 * t532 - t464 * t535;
t534 = cos(qJ(6));
t531 = sin(qJ(6));
t620 = -t463 * t535 - t464 * t532;
t602 = t620 * t531;
t355 = t534 * t409 - t602;
t479 = qJD(6) + t487;
t624 = t355 * t479;
t628 = t409 * t487;
t548 = t409 * t531 + t534 * t620;
t623 = t479 * t548;
t591 = t535 * t530;
t507 = t528 * t532 - t591;
t582 = t487 * t507;
t509 = t528 * t535 + t530 * t532;
t499 = t509 * qJD(5);
t581 = t509 * t492 + t499;
t627 = t487 * t620;
t435 = pkin(2) * t578 + pkin(3) * t495 + qJ(4) * t492;
t611 = -qJ(3) - pkin(7);
t516 = t611 * t536;
t512 = qJD(1) * t516;
t500 = t529 * t512;
t515 = t611 * t533;
t511 = qJD(1) * t515;
t454 = t511 * t609 + t500;
t391 = t528 * t435 + t530 * t454;
t596 = t492 * t528;
t375 = pkin(8) * t596 + t391;
t625 = -qJD(4) * t530 + t375;
t450 = -t507 * t531 + t509 * t534;
t585 = qJD(6) * t450 - t582 * t531 + t581 * t534;
t569 = -pkin(2) * t536 - pkin(1);
t553 = t569 * qJD(1);
t514 = qJD(3) + t553;
t426 = pkin(3) * t492 - qJ(4) * t495 + t514;
t610 = qJD(2) * pkin(2);
t505 = t511 + t610;
t565 = t609 * t512;
t447 = t529 * t505 - t565;
t443 = qJD(2) * qJ(4) + t447;
t380 = t530 * t426 - t443 * t528;
t348 = pkin(4) * t492 - pkin(8) * t463 + t380;
t381 = t528 * t426 + t530 * t443;
t362 = pkin(8) * t464 + t381;
t333 = t348 * t532 + t362 * t535;
t327 = -pkin(9) * t409 + t333;
t574 = qJD(6) * t531;
t325 = t327 * t574;
t446 = t505 * t609 + t500;
t438 = -qJD(2) * pkin(3) + qJD(4) - t446;
t402 = -pkin(4) * t464 + t438;
t353 = t409 * pkin(5) + t402;
t622 = t353 * t355 + t325;
t572 = qJD(1) * qJD(2);
t568 = t533 * t572;
t481 = qJD(2) * t518 - t529 * t568;
t575 = qJD(5) * t535;
t598 = t481 * t528;
t366 = t464 * t575 + t481 * t591 + (-qJD(5) * t463 - t598) * t532;
t494 = t508 * qJD(2);
t480 = qJD(1) * t494;
t520 = pkin(2) * t568;
t401 = pkin(3) * t480 - qJ(4) * t481 - qJD(4) * t495 + t520;
t566 = qJD(2) * t611;
t488 = qJD(3) * t536 + t533 * t566;
t470 = t488 * qJD(1);
t489 = -qJD(3) * t533 + t536 * t566;
t471 = t489 * qJD(1);
t422 = t609 * t470 + t529 * t471;
t416 = qJD(2) * qJD(4) + t422;
t350 = t530 * t401 - t416 * t528;
t597 = t481 * t530;
t339 = pkin(4) * t480 - pkin(8) * t597 + t350;
t351 = t528 * t401 + t530 * t416;
t341 = -pkin(8) * t598 + t351;
t561 = t535 * t339 - t341 * t532;
t539 = -t333 * qJD(5) + t561;
t316 = pkin(5) * t480 - pkin(9) * t366 + t539;
t367 = -qJD(5) * t620 + t481 * t509;
t576 = qJD(5) * t532;
t543 = t532 * t339 + t535 * t341 + t348 * t575 - t362 * t576;
t317 = -pkin(9) * t367 + t543;
t562 = t534 * t316 - t531 * t317;
t621 = t353 * t548 + t562;
t619 = t480 * MDP(28) + (-t355 ^ 2 + t548 ^ 2) * MDP(25) - t355 * MDP(24) * t548;
t618 = -0.2e1 * t572;
t617 = MDP(5) * (t533 ^ 2 - t536 ^ 2);
t544 = -t529 * t533 + t564;
t445 = -pkin(3) * t544 - qJ(4) * t508 + t569;
t459 = t529 * t515 - t516 * t609;
t397 = t530 * t445 - t459 * t528;
t593 = t508 * t530;
t374 = -pkin(4) * t544 - pkin(8) * t593 + t397;
t398 = t528 * t445 + t530 * t459;
t594 = t508 * t528;
t384 = -pkin(8) * t594 + t398;
t587 = t532 * t374 + t535 * t384;
t521 = pkin(2) * t529 + qJ(4);
t612 = pkin(8) + t521;
t503 = t612 * t528;
t504 = t612 * t530;
t583 = -t532 * t503 + t535 * t504;
t449 = t534 * t507 + t509 * t531;
t586 = -qJD(6) * t449 - t531 * t581 - t534 * t582;
t616 = -t450 * t480 - t479 * t586;
t615 = -t480 * t509 + t487 * t582;
t559 = t531 * t366 + t534 * t367;
t323 = -qJD(6) * t548 + t559;
t390 = t530 * t435 - t454 * t528;
t613 = pkin(8) * t530;
t361 = pkin(4) * t495 + t492 * t613 + t390;
t545 = qJD(4) * t528 + qJD(5) * t504;
t614 = t503 * t575 + t625 * t535 + (t361 + t545) * t532;
t332 = t535 * t348 - t362 * t532;
t326 = pkin(9) * t620 + t332;
t324 = pkin(5) * t487 + t326;
t608 = t324 * t534;
t607 = t327 * t534;
t606 = t355 * t495;
t605 = t548 * t495;
t604 = t409 * t495;
t603 = t620 * t495;
t421 = t470 * t529 - t609 * t471;
t458 = -t609 * t515 - t516 * t529;
t601 = t421 * t458;
t497 = t544 * qJD(2);
t595 = t497 * t528;
t537 = qJD(2) ^ 2;
t592 = t533 * t537;
t590 = t536 * t537;
t538 = qJD(1) ^ 2;
t589 = t536 * t538;
t571 = t533 * t610;
t415 = pkin(3) * t494 - qJ(4) * t497 - qJD(4) * t508 + t571;
t437 = t488 * t609 + t529 * t489;
t370 = t528 * t415 + t530 * t437;
t573 = qJD(6) * t534;
t570 = t534 * t366 - t531 * t367 - t409 * t573;
t453 = t511 * t529 - t565;
t418 = -pkin(4) * t596 + t453;
t567 = pkin(5) * t581 - t418;
t563 = pkin(1) * t618;
t369 = t530 * t415 - t437 * t528;
t345 = pkin(4) * t494 - t497 * t613 + t369;
t352 = -pkin(8) * t595 + t370;
t560 = t535 * t345 - t352 * t532;
t558 = t535 * t374 - t384 * t532;
t436 = t488 * t529 - t609 * t489;
t557 = -t535 * t503 - t504 * t532;
t556 = qJD(6) * t324 + t317;
t524 = -pkin(2) * t609 - pkin(3);
t555 = -t449 * t480 - t479 * t585;
t554 = -t507 * t480 - t487 * t581;
t396 = pkin(4) * t598 + t421;
t404 = pkin(4) * t595 + t436;
t432 = pkin(4) * t594 + t458;
t360 = t535 * t361;
t420 = -pkin(9) * t507 + t583;
t552 = pkin(5) * t495 - pkin(9) * t582 + t509 * qJD(4) + qJD(5) * t583 + qJD(6) * t420 - t375 * t532 + t360;
t419 = -pkin(9) * t509 + t557;
t551 = pkin(9) * t581 - qJD(6) * t419 + t614;
t319 = t324 * t531 + t607;
t549 = -t380 * t528 + t381 * t530;
t547 = t421 * t508 + t458 * t481;
t439 = t509 * t508;
t440 = t507 * t508;
t387 = t534 * t439 - t440 * t531;
t388 = -t439 * t531 - t440 * t534;
t513 = -t530 * pkin(4) + t524;
t542 = t532 * t345 + t535 * t352 + t374 * t575 - t384 * t576;
t322 = t574 * t620 + t570;
t541 = t438 * t497 + t547;
t540 = -t480 * t521 + t481 * t524 + (-qJD(4) + t438) * t492;
t490 = t492 ^ 2;
t460 = t507 * pkin(5) + t513;
t451 = t480 * t544;
t389 = pkin(5) * t439 + t432;
t386 = t497 * t509 + t575 * t593 - t576 * t594;
t385 = -t497 * t507 - t499 * t508;
t342 = pkin(5) * t386 + t404;
t336 = pkin(5) * t367 + t396;
t335 = -pkin(9) * t439 + t587;
t334 = -pkin(5) * t544 + pkin(9) * t440 + t558;
t331 = qJD(6) * t388 + t385 * t531 + t534 * t386;
t330 = -qJD(6) * t387 + t385 * t534 - t386 * t531;
t321 = -pkin(9) * t386 + t542;
t320 = pkin(5) * t494 - pkin(9) * t385 - qJD(5) * t587 + t560;
t318 = -t327 * t531 + t608;
t1 = [(-t366 * t439 + t367 * t440 - t385 * t409 + t386 * t620) * MDP(18) + (-t366 * t440 - t385 * t620) * MDP(17) + (-t366 * t544 + t385 * t487 - t440 * t480 - t494 * t620) * MDP(19) + (-t333 * t494 + t432 * t366 + t402 * t385 - t396 * t440 - t404 * t620 - t480 * t587 - t487 * t542 + t543 * t544) * MDP(23) + (t560 * t487 + t558 * t480 - t561 * t544 + t332 * t494 + t404 * t409 + t432 * t367 + t396 * t439 + t402 * t386 + (t333 * t544 - t487 * t587) * qJD(5)) * MDP(22) + (t479 * t494 - t451) * MDP(28) + (t322 * t388 - t330 * t548) * MDP(24) + (-t322 * t387 - t323 * t388 - t330 * t355 + t331 * t548) * MDP(25) + (-t369 * t463 + t370 * t464 + (-t350 * t508 - t380 * t497 - t397 * t481) * t530 + (-t351 * t508 - t381 * t497 - t398 * t481) * t528) * MDP(15) + (t422 * t544 + t436 * t495 - t437 * t492 - t446 * t497 - t447 * t494 - t459 * t480 + t547) * MDP(11) + (t351 * t544 - t370 * t492 - t381 * t494 - t398 * t480 + t436 * t463 + t530 * t541) * MDP(14) + (t367 * t544 - t386 * t487 - t409 * t494 - t439 * t480) * MDP(20) + (t323 * t544 - t331 * t479 - t355 * t494 - t387 * t480) * MDP(27) + ((t320 * t534 - t321 * t531) * t479 + (t334 * t534 - t335 * t531) * t480 - t562 * t544 + t318 * t494 + t342 * t355 + t389 * t323 + t336 * t387 + t353 * t331 + ((-t334 * t531 - t335 * t534) * t479 + t319 * t544) * qJD(6)) * MDP(29) + (-t322 * t544 + t330 * t479 + t388 * t480 - t494 * t548) * MDP(26) + (-t319 * t494 + t389 * t322 - t325 * t544 + t353 * t330 + t336 * t388 - t342 * t548 + (-(-qJD(6) * t335 + t320) * t479 - t334 * t480 + t316 * t544) * t531 + (-(qJD(6) * t334 + t321) * t479 - t335 * t480 + t556 * t544) * t534) * MDP(30) + (-t350 * t544 + t369 * t492 + t380 * t494 + t397 * t480 - t436 * t464 + t528 * t541) * MDP(13) + (t487 * t494 - t451) * MDP(21) + 0.2e1 * t536 * MDP(4) * t568 + t617 * t618 + MDP(6) * t590 + (-pkin(7) * t590 + t533 * t563) * MDP(9) - MDP(7) * t592 + (pkin(7) * t592 + t536 * t563) * MDP(10) + (t350 * t397 + t351 * t398 + t369 * t380 + t370 * t381 + t436 * t438 + t601) * MDP(16) + (t601 + t422 * t459 - t446 * t436 + t447 * t437 + (t514 + t553) * t571) * MDP(12); (-t380 * t390 - t381 * t391 + t421 * t524 - t438 * t453 + (-t350 * t528 + t351 * t530) * t521 + t549 * qJD(4)) * MDP(16) + (t366 * t509 + t582 * t620) * MDP(17) + (-t366 * t507 - t367 * t509 + t409 * t582 + t581 * t620) * MDP(18) + (t603 - t615) * MDP(19) + (t554 + t604) * MDP(20) + (t557 * t480 + t513 * t367 + t396 * t507 - t418 * t409 + (-t360 - t545 * t535 + (qJD(5) * t503 + t625) * t532) * t487 + t581 * t402) * MDP(22) + (t513 * t366 + t396 * t509 - t582 * t402 + t418 * t620 - t583 * t480 + t487 * t614) * MDP(23) + (t322 * t450 - t548 * t586) * MDP(24) + (-t322 * t449 - t323 * t450 - t355 * t586 + t548 * t585) * MDP(25) + (t605 - t616) * MDP(26) + (t555 + t606) * MDP(27) + ((t419 * t534 - t420 * t531) * t480 + t460 * t323 + t336 * t449 + (t531 * t551 - t534 * t552) * t479 + t567 * t355 + t585 * t353) * MDP(29) + (-(t419 * t531 + t420 * t534) * t480 + t460 * t322 + t336 * t450 + (t531 * t552 + t534 * t551) * t479 - t567 * t548 + t586 * t353) * MDP(30) - t533 * MDP(4) * t589 + t538 * t617 + ((-t446 + t454) * t492 + (-t480 * t529 - t481 * t609) * pkin(2)) * MDP(11) + (t446 * t453 - t447 * t454 + (-t421 * t609 + t422 * t529 - t514 * t578) * pkin(2)) * MDP(12) + (-t390 * t492 - t421 * t530 + t453 * t464 + t528 * t540) * MDP(13) + (t391 * t492 + t421 * t528 - t453 * t463 + t530 * t540) * MDP(14) + (t390 * t463 - t391 * t464 + (qJD(4) * t464 - t380 * t492 + t351) * t530 + (qJD(4) * t463 - t381 * t492 - t350) * t528) * MDP(15) - (t332 * MDP(22) - t333 * MDP(23) + t318 * MDP(29) - t319 * MDP(30) + (-t447 + t453) * MDP(11) + t380 * MDP(13) - t381 * MDP(14) + t487 * MDP(21) + t479 * MDP(28)) * t495 + (MDP(9) * t533 * t538 + MDP(10) * t589) * pkin(1); (-t495 ^ 2 - t490) * MDP(11) + (t446 * t495 + t520) * MDP(12) + (t464 * t495 + t480 * t530) * MDP(13) + (-t463 * t495 - t480 * t528 - t490 * t530) * MDP(14) + (t350 * t530 + t351 * t528 - t438 * t495) * MDP(16) + (t554 - t604) * MDP(22) + (t603 + t615) * MDP(23) + (t555 - t606) * MDP(29) + (t605 + t616) * MDP(30) + (-t528 ^ 2 - t530 ^ 2) * MDP(15) * t481 + (t447 * MDP(12) + (t463 * t528 + t464 * t530) * MDP(15) + t549 * MDP(16) - MDP(13) * t596) * t492; (t463 * t492 + t598) * MDP(13) + (t464 * t492 + t597) * MDP(14) + (-t463 ^ 2 - t464 ^ 2) * MDP(15) + (t380 * t463 - t381 * t464 + t421) * MDP(16) + (t367 - t627) * MDP(22) + (t366 - t628) * MDP(23) + (t323 - t623) * MDP(29) + (t322 - t624) * MDP(30); -t620 * t409 * MDP(17) + (-t409 ^ 2 + t620 ^ 2) * MDP(18) + (t366 + t628) * MDP(19) + (-t367 - t627) * MDP(20) + t480 * MDP(21) + (t333 * t487 + t402 * t620 + t539) * MDP(22) + (t332 * t487 + t402 * t409 - t543) * MDP(23) + (t322 + t624) * MDP(26) + (-t323 - t623) * MDP(27) + (-(-t326 * t531 - t607) * t479 - t319 * qJD(6) + (t355 * t620 - t479 * t574 + t480 * t534) * pkin(5) + t621) * MDP(29) + ((-t327 * t479 - t316) * t531 + (t326 * t479 - t556) * t534 + (-t479 * t573 - t480 * t531 - t548 * t620) * pkin(5) + t622) * MDP(30) + t619; (t570 + t624) * MDP(26) + (-t559 - t623) * MDP(27) + (t319 * t479 + t621) * MDP(29) + (-t531 * t316 - t534 * t317 + t318 * t479 + t622) * MDP(30) + (MDP(26) * t602 + MDP(27) * t548 - MDP(29) * t319 - MDP(30) * t608) * qJD(6) + t619;];
tauc  = t1;
