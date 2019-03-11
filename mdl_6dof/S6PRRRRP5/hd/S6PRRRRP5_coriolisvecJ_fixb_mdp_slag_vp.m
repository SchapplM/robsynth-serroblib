% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRRRP5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRRRP5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:26:12
% EndTime: 2019-03-09 00:26:28
% DurationCPUTime: 10.20s
% Computational Cost: add. (6418->522), mult. (17229->737), div. (0->0), fcn. (13950->12), ass. (0->226)
t508 = sin(pkin(7));
t514 = sin(qJ(3));
t634 = t508 * t514;
t501 = pkin(9) * t634;
t510 = cos(pkin(7));
t518 = cos(qJ(3));
t519 = cos(qJ(2));
t622 = t518 * t519;
t515 = sin(qJ(2));
t627 = t514 * t515;
t535 = -t510 * t627 + t622;
t509 = sin(pkin(6));
t602 = qJD(1) * t509;
t630 = t510 * t518;
t611 = t535 * t602 - (pkin(2) * t630 - t501) * qJD(3);
t549 = pkin(3) * t514 - pkin(10) * t518;
t534 = t549 * qJD(3);
t581 = t515 * t602;
t664 = (t534 - t581) * t508;
t597 = qJD(2) * t518;
t500 = t508 * t597;
t545 = t500 - qJD(4);
t633 = t508 * t518;
t586 = pkin(9) * t633;
t463 = t586 + (pkin(2) * t514 + pkin(10)) * t510;
t550 = -pkin(3) * t518 - pkin(10) * t514;
t464 = (-pkin(2) + t550) * t508;
t513 = sin(qJ(4));
t517 = cos(qJ(4));
t592 = qJD(4) * t517;
t594 = qJD(4) * t513;
t663 = -t463 * t594 + t464 * t592 + t664 * t513 - t611 * t517;
t625 = t515 * t518;
t626 = t514 * t519;
t537 = t510 * t625 + t626;
t631 = t510 * t514;
t609 = -t537 * t602 + (pkin(2) * t631 + t586) * qJD(3);
t596 = qJD(3) * t514;
t579 = t508 * t596;
t662 = -pkin(11) * t579 - t663;
t475 = -t517 * t510 + t513 * t634;
t595 = qJD(3) * t518;
t578 = t508 * t595;
t438 = -qJD(4) * t475 + t517 * t578;
t476 = t510 * t513 + t517 * t634;
t439 = qJD(4) * t476 + t513 * t578;
t661 = pkin(4) * t439 - pkin(11) * t438 + t609;
t511 = cos(pkin(6));
t601 = qJD(1) * t511;
t582 = t508 * t601;
t489 = t514 * t582;
t600 = qJD(2) * t508;
t481 = pkin(9) * t600 + t581;
t645 = qJD(2) * pkin(2);
t488 = t519 * t602 + t645;
t652 = t518 * t481 + t488 * t631;
t417 = t489 + t652;
t660 = t417 + t545 * (pkin(4) * t513 - pkin(11) * t517);
t598 = qJD(2) * t510;
t565 = qJD(3) + t598;
t580 = t514 * t600;
t524 = -t513 * t580 + t517 * t565;
t587 = qJD(2) * qJD(3);
t573 = t508 * t587;
t552 = t518 * t573;
t542 = t517 * t552;
t521 = t524 * qJD(4) + t542;
t659 = -qJD(5) * t545 + t521;
t505 = t508 ^ 2;
t658 = (t514 * t518 * MDP(5) - (t514 ^ 2 - t518 ^ 2) * MDP(6)) * t505;
t647 = -qJ(6) - pkin(11);
t657 = qJ(6) * t524 + qJD(5) * t647;
t512 = sin(qJ(5));
t516 = cos(qJ(5));
t462 = t501 + (-pkin(2) * t518 - pkin(3)) * t510;
t408 = pkin(4) * t475 - pkin(11) * t476 + t462;
t610 = t517 * t463 + t513 * t464;
t410 = -pkin(11) * t633 + t610;
t613 = t512 * t408 + t516 * t410;
t655 = -qJD(5) * t613 + t512 * t662 + t661 * t516;
t590 = qJD(5) * t516;
t591 = qJD(5) * t512;
t654 = -t408 * t590 + t410 * t591 - t661 * t512 + t516 * t662;
t585 = pkin(10) * t594;
t653 = t512 * t585 - t516 * t660;
t651 = t463 * t592 + t464 * t594 - t611 * t513 - t664 * t517;
t470 = t514 * t481;
t635 = t488 * t510;
t416 = t518 * (t582 + t635) - t470;
t465 = t549 * t600;
t612 = t517 * t416 + t513 * t465;
t382 = pkin(11) * t580 + t612;
t495 = -pkin(4) * t517 - pkin(11) * t513 - pkin(3);
t650 = t516 * t382 - t495 * t590 + t512 * t660;
t457 = t513 * t565 + t517 * t580;
t428 = t516 * t457 - t512 * t545;
t649 = t428 ^ 2;
t520 = qJD(2) ^ 2;
t648 = pkin(5) * t512;
t646 = pkin(10) * qJD(4);
t642 = qJ(6) * t513;
t553 = t514 * t573;
t526 = t535 * qJD(2);
t555 = t511 * t578;
t384 = (t488 * t630 - t470) * qJD(3) + (t509 * t526 + t555) * qJD(1);
t405 = pkin(10) * t565 + t417;
t499 = t510 * t601;
t423 = t499 + (qJD(2) * t550 - t488) * t508;
t433 = (t534 + t581) * t600;
t561 = -t513 * t384 - t405 * t592 - t423 * t594 + t517 * t433;
t346 = -pkin(4) * t553 - t561;
t641 = t346 * t512;
t640 = t346 * t516;
t374 = t457 * t591 - t512 * t553 - t516 * t659;
t639 = t374 * t512;
t426 = t457 * t512 + t516 * t545;
t451 = qJD(5) - t524;
t638 = t426 * t451;
t637 = t428 * t451;
t543 = t513 * t552;
t430 = qJD(4) * t457 + t543;
t636 = t430 * t516;
t632 = t509 * t520;
t629 = t512 * t430;
t628 = t513 * t516;
t624 = t516 * t517;
t623 = t516 * t518;
t440 = t476 * t512 + t508 * t623;
t387 = qJD(5) * t440 - t516 * t438 - t512 * t579;
t584 = t512 * t633;
t441 = t476 * t516 - t584;
t621 = pkin(5) * t439 + qJ(6) * t387 - qJD(6) * t441 + t655;
t388 = -qJD(5) * t584 + t438 * t512 + t476 * t590 - t516 * t579;
t620 = -qJ(6) * t388 - qJD(6) * t440 - t654;
t367 = t517 * t405 + t513 * t423;
t362 = -pkin(11) * t545 + t367;
t404 = -pkin(3) * t565 - t416;
t369 = -pkin(4) * t524 - t457 * pkin(11) + t404;
t347 = -t362 * t512 + t516 * t369;
t342 = -qJ(6) * t428 + t347;
t341 = pkin(5) * t451 + t342;
t619 = t341 - t342;
t366 = -t513 * t405 + t517 * t423;
t415 = pkin(4) * t457 - pkin(11) * t524;
t618 = t516 * t366 + t512 * t415;
t617 = -pkin(4) * t579 + t651;
t449 = (t512 * t514 + t517 * t623) * t600;
t503 = pkin(10) * t624;
t558 = t513 * t500;
t589 = qJD(6) * t516;
t615 = -pkin(5) * t558 + qJ(6) * t449 + t382 * t512 - t513 * t589 + (pkin(5) * t513 - qJ(6) * t624) * qJD(4) + (-t503 + (-t495 + t642) * t512) * qJD(5) + t653;
t557 = t517 * t500;
t448 = t512 * t557 - t516 * t580;
t614 = qJ(6) * t448 + (-qJ(6) * qJD(5) - t646) * t628 + (-qJD(6) * t513 + (-pkin(10) * qJD(5) - qJ(6) * qJD(4)) * t517) * t512 - t650;
t608 = t512 * t657 + t589 - t618;
t412 = t516 * t415;
t607 = -pkin(5) * t457 - t412 + t657 * t516 + (-qJD(6) + t366) * t512;
t604 = t512 * t495 + t503;
t599 = qJD(2) * t509;
t593 = qJD(4) * t516;
t588 = t404 * qJD(4);
t583 = t515 * t632;
t577 = t451 * t591;
t575 = qJD(3) * t635;
t574 = qJD(1) * t599;
t571 = t516 * t408 - t410 * t512;
t570 = -t513 * t463 + t464 * t517;
t569 = t518 * t545;
t568 = t451 * t516;
t567 = qJD(4) * t545;
t562 = t505 * t583;
t559 = t508 * t515 * t599;
t556 = t511 * t579;
t551 = MDP(16) * t580;
t409 = pkin(4) * t633 - t570;
t547 = -t512 * t592 + t448;
t546 = t516 * t592 - t449;
t348 = t362 * t516 + t369 * t512;
t536 = t510 * t626 + t625;
t437 = t509 * t536 + t511 * t634;
t474 = -t508 * t509 * t519 + t510 * t511;
t407 = t437 * t517 + t474 * t513;
t538 = t510 * t622 - t627;
t436 = -t509 * t538 - t511 * t633;
t373 = t407 * t516 + t436 * t512;
t372 = -t407 * t512 + t436 * t516;
t406 = t437 * t513 - t474 * t517;
t544 = t510 * t515 * t574;
t413 = t513 * t416;
t541 = pkin(4) * t580 - t413;
t450 = -t488 * t508 + t499;
t539 = t450 * t508 - t505 * t645;
t532 = -t451 * t590 - t629;
t361 = pkin(4) * t545 - t366;
t531 = -pkin(11) * t430 + t361 * t451;
t530 = -t517 * t384 + t405 * t594 - t423 * t592 - t513 * t433;
t345 = pkin(11) * t553 - t530;
t357 = t430 * pkin(4) - pkin(11) * t521 + qJD(3) * t489 + t481 * t595 + t514 * t575 + t518 * t544 + t574 * t626;
t334 = t516 * t345 + t512 * t357 - t362 * t591 + t369 * t590;
t527 = t537 * qJD(2);
t525 = qJD(3) * t481 + t544;
t335 = -qJD(5) * t348 - t345 * t512 + t516 * t357;
t375 = t457 * t590 + t512 * t659 - t516 * t553;
t338 = pkin(5) * t375 + t346;
t522 = -t450 * t600 - t575 + (-qJD(3) * t508 * t511 - t519 * t599) * qJD(1);
t497 = t647 * t516;
t496 = t647 * t512;
t483 = t516 * t495;
t442 = -t512 * t642 + t604;
t435 = -qJ(6) * t628 + t483 + (-pkin(10) * t512 - pkin(5)) * t517;
t425 = t426 ^ 2;
t400 = t555 + (qJD(3) * t538 + t526) * t509;
t399 = t556 + (qJD(3) * t536 + t527) * t509;
t385 = t652 * qJD(3) + (t509 * t527 + t556) * qJD(1);
t381 = -t465 * t517 - t541;
t360 = -qJD(4) * t406 + t400 * t517 + t513 * t559;
t359 = qJD(4) * t407 + t400 * t513 - t517 * t559;
t358 = -qJ(6) * t440 + t613;
t354 = t426 * pkin(5) + qJD(6) + t361;
t353 = pkin(5) * t475 - qJ(6) * t441 + t571;
t343 = -qJ(6) * t426 + t348;
t340 = qJD(5) * t372 + t360 * t516 + t399 * t512;
t339 = -qJD(5) * t373 - t360 * t512 + t399 * t516;
t333 = -qJ(6) * t375 - qJD(6) * t426 + t334;
t332 = pkin(5) * t430 + qJ(6) * t374 - qJD(6) * t428 + t335;
t1 = [-MDP(3) * t583 - t519 * MDP(4) * t632 + (-t399 * t565 + t474 * t553 - t518 * t562) * MDP(10) + (-t400 * t565 + t474 * t552 + t514 * t562) * MDP(11) + (t359 * t545 - t399 * t524 - t406 * t553 + t436 * t430) * MDP(17) + (t360 * t545 + t399 * t457 - t407 * t553 + t436 * t521) * MDP(18) + (t339 * t451 + t359 * t426 + t372 * t430 + t375 * t406) * MDP(24) + (-t340 * t451 + t359 * t428 - t373 * t430 - t374 * t406) * MDP(25) + (-t339 * t428 - t340 * t426 + t372 * t374 - t373 * t375) * MDP(26) + (t332 * t372 + t333 * t373 + t338 * t406 + t339 * t341 + t340 * t343 + t354 * t359) * MDP(27); ((-qJD(2) * t609 - t385) * t510 + (t514 * t539 - t609) * qJD(3)) * MDP(10) + ((qJD(2) * t611 - t384) * t510 + (t518 * t539 + t611) * qJD(3)) * MDP(11) + (t457 * t438 + t476 * t521) * MDP(12) + (-t476 * t430 + t438 * t524 - t457 * t439 - t475 * t521) * MDP(13) + (-t438 * t545 + t457 * t579 + t476 * t553 - t521 * t633) * MDP(14) + (t439 * t545 + (t430 * t518 + (-qJD(2) * t475 + t524) * t596) * t508) * MDP(15) + (-t505 * t597 - t508 * t545) * MDP(16) * t596 + (t462 * t430 + t385 * t475 + t404 * t439 - t609 * t524 + (-t561 * t518 + (qJD(2) * t570 + t366) * t596) * t508 + t651 * t545) * MDP(17) + (-t367 * t579 + t385 * t476 + t404 * t438 + t609 * t457 + t462 * t521 - t530 * t633 + t545 * t663 - t553 * t610) * MDP(18) + (-t374 * t441 - t387 * t428) * MDP(19) + (t374 * t440 - t375 * t441 + t387 * t426 - t388 * t428) * MDP(20) + (-t374 * t475 - t387 * t451 + t428 * t439 + t430 * t441) * MDP(21) + (-t375 * t475 - t388 * t451 - t426 * t439 - t430 * t440) * MDP(22) + (t430 * t475 + t439 * t451) * MDP(23) + (t335 * t475 + t346 * t440 + t347 * t439 + t361 * t388 + t409 * t375 + t617 * t426 + t571 * t430 + t451 * t655) * MDP(24) + (-t334 * t475 + t346 * t441 - t348 * t439 - t361 * t387 - t409 * t374 + t617 * t428 - t613 * t430 + t451 * t654) * MDP(25) + (-t332 * t441 - t333 * t440 + t341 * t387 - t343 * t388 + t353 * t374 - t358 * t375 - t426 * t620 - t428 * t621) * MDP(26) + (t333 * t358 + t332 * t353 + t338 * (pkin(5) * t440 + t409) + (pkin(5) * t388 + t617) * t354 + t620 * t343 + t621 * t341) * MDP(27) + 0.2e1 * t658 * t587 + (MDP(7) * t578 - MDP(8) * t579) * (qJD(3) + 0.2e1 * t598); (t417 * t565 + t514 * t522 - t518 * t525) * MDP(10) + (t416 * t565 + t514 * t525 + t518 * t522) * MDP(11) + (-qJD(4) * t513 ^ 2 * t580 + ((qJD(4) * t565 + t552) * t513 - t545 * t457) * t517) * MDP(12) + (-t513 * t430 + t517 * t521 + (t558 - t594) * t457 - (t557 - t592) * t524) * MDP(13) + (-t517 * t567 + (t517 * t569 + (qJD(3) * t513 - t457) * t514) * t600) * MDP(14) + (t513 * t567 + (-t513 * t569 + (qJD(3) * t517 - t524) * t514) * t600) * MDP(15) + t545 * t551 + (-pkin(3) * t430 + t513 * t588 - t413 * t545 + t417 * t524 + (pkin(10) * t567 + t465 * t545 - t385) * t517 + (-t366 * t514 + (-pkin(10) * t596 - t404 * t518) * t513) * t600) * MDP(17) + (-pkin(3) * t521 + t367 * t580 + t385 * t513 - t404 * t557 - t417 * t457 + (-t585 - t612) * t545 + (-pkin(10) * t553 + t588) * t517) * MDP(18) + (-t374 * t628 + (-t513 * t591 + t546) * t428) * MDP(19) + (t426 * t449 + t428 * t448 + (-t426 * t516 - t428 * t512) * t592 + (t639 - t375 * t516 + (t426 * t512 - t428 * t516) * qJD(5)) * t513) * MDP(20) + (t374 * t517 + t546 * t451 + (-t428 * t545 - t577 + t636) * t513) * MDP(21) + (t375 * t517 + t547 * t451 + (t426 * t545 + t532) * t513) * MDP(22) + (-t451 * t513 * t545 - t430 * t517) * MDP(23) + (-t361 * t448 - t381 * t426 + t483 * t430 + ((-qJD(5) * t495 + t382) * t512 + t653) * t451 + (t361 * t512 * qJD(4) - t335 + (qJD(4) * t426 + t532) * pkin(10)) * t517 + (pkin(10) * t375 - t347 * t545 + t361 * t590 + t641) * t513) * MDP(24) + (-t604 * t430 - t381 * t428 - t361 * t449 + t650 * t451 + (t361 * t593 + t334 + (qJD(4) * t428 + t577) * pkin(10)) * t517 + (-t361 * t591 + t640 + t545 * t348 + (t451 * t593 - t374) * pkin(10)) * t513) * MDP(25) + (t341 * t449 + t343 * t448 + t374 * t435 - t375 * t442 - t615 * t428 - t614 * t426 + (-t341 * t516 - t343 * t512) * t592 + (-t332 * t516 - t333 * t512 + (t341 * t512 - t343 * t516) * qJD(5)) * t513) * MDP(26) + (t332 * t435 + t333 * t442 + t338 * (pkin(10) + t648) * t513 + ((t465 + t646) * t517 + (t513 * t590 - t547) * pkin(5) + t541) * t354 + t614 * t343 + t615 * t341) * MDP(27) + ((-MDP(7) * t518 + MDP(8) * t514) * t508 * t510 - t658) * t520; -t524 ^ 2 * MDP(13) + (t500 * t524 + t542) * MDP(14) - t543 * MDP(15) + qJD(3) * t551 + (-t367 * t545 + t561) * MDP(17) + (-t366 * t545 - t404 * t524 + t530) * MDP(18) + (t428 * t568 - t639) * MDP(19) + ((-t374 - t638) * t516 + (-t375 - t637) * t512) * MDP(20) + (t451 * t568 + t629) * MDP(21) + (-t451 ^ 2 * t512 + t636) * MDP(22) + (-pkin(4) * t375 - t640 - t367 * t426 + (-pkin(11) * t590 - t412) * t451 + (t366 * t451 + t531) * t512) * MDP(24) + (pkin(4) * t374 + t641 - t367 * t428 + (pkin(11) * t591 + t618) * t451 + t531 * t516) * MDP(25) + (t374 * t496 + t375 * t497 - t607 * t428 - t608 * t426 + (-t341 * t451 + t333) * t516 + (-t343 * t451 - t332) * t512) * MDP(26) + (-t333 * t497 + t332 * t496 + t338 * (-pkin(5) * t516 - pkin(4)) + (t451 * t648 - t367) * t354 + t608 * t343 + t607 * t341) * MDP(27) + (-MDP(12) * t524 + MDP(13) * t457 - t500 * MDP(15) - t404 * MDP(17) - t428 * MDP(21) + t426 * MDP(22) - t451 * MDP(23) - t347 * MDP(24) + t348 * MDP(25)) * t457; t428 * t426 * MDP(19) + (-t425 + t649) * MDP(20) + (-t374 + t638) * MDP(21) + (-t375 + t637) * MDP(22) + t430 * MDP(23) + (t348 * t451 - t361 * t428 + t335) * MDP(24) + (t347 * t451 + t361 * t426 - t334) * MDP(25) + (pkin(5) * t374 - t426 * t619) * MDP(26) + (t619 * t343 + (-t354 * t428 + t332) * pkin(5)) * MDP(27); (-t425 - t649) * MDP(26) + (t341 * t428 + t343 * t426 + t338) * MDP(27);];
tauc  = t1;
