% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRRR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [14x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRPRRR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:56:10
% EndTime: 2019-03-08 20:56:25
% DurationCPUTime: 8.67s
% Computational Cost: add. (8698->522), mult. (28074->779), div. (0->0), fcn. (25004->16), ass. (0->238)
t502 = sin(pkin(14));
t505 = sin(pkin(6));
t506 = cos(pkin(14));
t516 = cos(qJ(2));
t512 = sin(qJ(2));
t643 = cos(pkin(7));
t583 = t512 * t643;
t521 = t505 * (-t502 * t516 - t506 * t583);
t462 = qJD(1) * t521;
t503 = sin(pkin(8));
t507 = cos(pkin(8));
t504 = sin(pkin(7));
t610 = qJD(1) * t505;
t594 = t512 * t610;
t571 = t504 * t594;
t608 = qJD(3) * t504;
t633 = t502 * t503;
t542 = t462 * t503 - t507 * t571 + t608 * t633;
t515 = cos(qJ(4));
t581 = t643 * t503;
t564 = t515 * t581;
t624 = t507 * t515;
t598 = t506 * t624;
t661 = t504 * t598 + t564;
t522 = t505 * (-t502 * t583 + t506 * t516);
t466 = qJD(1) * t522;
t511 = sin(qJ(4));
t525 = t462 * t507 + t503 * t571;
t625 = t507 * t511;
t526 = (-t502 * t625 + t506 * t515) * t504;
t628 = t504 * t507;
t523 = (t506 * t628 + t581) * pkin(10);
t588 = t502 * t643;
t629 = t504 * t506;
t612 = pkin(2) * t588 + qJ(3) * t629;
t451 = t523 + t612;
t586 = t506 * t643;
t497 = pkin(2) * t586;
t596 = t643 * pkin(3);
t456 = t596 + t497 + (-pkin(10) * t507 - qJ(3)) * t504 * t502;
t543 = -pkin(3) * t506 - pkin(10) * t633;
t469 = (-pkin(2) + t543) * t504;
t545 = t456 * t507 + t469 * t503;
t647 = -t511 * t451 + t515 * t545;
t660 = qJD(3) * t526 + t647 * qJD(4) - t466 * t515 - t511 * t525;
t632 = t502 * t511;
t446 = (t564 + (t598 - t632) * t504) * qJD(4);
t455 = t504 * (t502 * t515 + t506 * t625) + t511 * t581;
t447 = t455 * qJD(4);
t659 = pkin(4) * t447 - pkin(11) * t446 + t542;
t609 = qJD(2) * t504;
t482 = qJ(3) * t609 + t594;
t593 = t516 * t610;
t487 = qJD(2) * pkin(2) + t593;
t508 = cos(pkin(6));
t627 = t504 * t508;
t595 = qJD(1) * t627;
t425 = t506 * t482 + t487 * t588 + t502 * t595;
t411 = qJD(2) * t523 + t425;
t424 = -t502 * t482 + t487 * t586 + t506 * t595;
t412 = (-pkin(10) * t502 * t628 + t596) * qJD(2) + t424;
t584 = t508 * t643;
t599 = qJD(1) * t584 + qJD(3);
t426 = (qJD(2) * t543 - t487) * t504 + t599;
t548 = t412 * t507 + t426 * t503;
t651 = t511 * t411 - t548 * t515;
t467 = qJD(2) * t526;
t606 = qJD(4) * t515;
t658 = -t503 * t606 + t467;
t592 = t502 * t609;
t657 = qJD(2) * t661 - t511 * t592;
t442 = qJD(5) - t657;
t445 = t455 * qJD(2);
t514 = cos(qJ(5));
t591 = t503 * t609;
t488 = t506 * t591;
t579 = qJD(2) * t643;
t540 = t507 * t579 - t488;
t529 = -qJD(4) - t540;
t468 = t514 * t529;
t510 = sin(qJ(5));
t414 = t445 * t510 + t468;
t413 = qJD(6) + t414;
t500 = t504 ^ 2;
t656 = t500 * (t502 ^ 2 + t506 ^ 2);
t518 = t451 * t515 + t511 * t545;
t527 = (t502 * t624 + t506 * t511) * t504;
t616 = qJD(3) * t527 + qJD(4) * t518 - t466 * t511 + t515 * t525;
t654 = -t510 * t660 + t659 * t514;
t417 = -t456 * t503 + t507 * t469;
t454 = t504 * t632 - t661;
t381 = pkin(4) * t454 - pkin(11) * t455 + t417;
t585 = t507 * t643;
t475 = t503 * t629 - t585;
t386 = -pkin(11) * t475 + t518;
t602 = qJD(5) * t514;
t604 = qJD(5) * t510;
t653 = -t381 * t602 + t386 * t604 - t659 * t510 - t514 * t660;
t652 = t510 * t381 + t514 * t386;
t650 = t447 * MDP(12);
t438 = qJD(2) * t447;
t410 = t515 * t411;
t631 = t503 * t511;
t366 = t412 * t625 + t426 * t631 + t410;
t362 = -pkin(11) * t529 + t366;
t387 = -t412 * t503 + t507 * t426;
t367 = -pkin(4) * t657 - pkin(11) * t445 + t387;
t343 = t362 * t514 + t367 * t510;
t477 = (t593 + t608) * qJD(2);
t626 = t505 * t512;
t590 = qJD(2) * t626;
t566 = qJD(1) * t590;
t532 = t643 * t566;
t448 = -t502 * t477 - t506 * t532;
t449 = t506 * t477 - t502 * t532;
t544 = t504 * t566;
t538 = t503 * t544;
t646 = (t448 * t507 + t538) * t511 + t515 * t449;
t351 = -t651 * qJD(4) + t646;
t421 = -t448 * t503 + t507 * t544;
t437 = t657 * qJD(4);
t374 = pkin(4) * t438 - pkin(11) * t437 + t421;
t577 = t351 * t510 - t514 * t374;
t644 = -qJD(5) * t343 - t577;
t333 = -pkin(5) * t438 - t644;
t416 = t514 * t445 - t510 * t529;
t649 = t413 * (pkin(5) * t416 + pkin(12) * t413) + t333;
t648 = t511 * t548 + t410;
t582 = t516 * t643;
t453 = t506 * t626 + (t505 * t582 + t627) * t502;
t452 = t506 * t627 + (-t502 * t512 + t506 * t582) * t505;
t476 = -t505 * t516 * t504 + t584;
t546 = t452 * t507 + t476 * t503;
t645 = -t453 * t511 + t515 * t546;
t389 = -qJD(5) * t468 + t514 * t437 - t445 * t604;
t509 = sin(qJ(6));
t513 = cos(qJ(6));
t600 = qJD(6) * t513;
t597 = t513 * t389 + t509 * t438 + t442 * t600;
t601 = qJD(6) * t509;
t356 = -t416 * t601 + t597;
t642 = t356 * t509;
t637 = t416 * t509;
t391 = -t513 * t442 + t637;
t641 = t391 * t413;
t393 = t416 * t513 + t442 * t509;
t640 = t393 * t413;
t639 = t414 * t442;
t638 = t416 * t442;
t636 = t657 * t514;
t635 = t448 * t502;
t630 = t503 * t515;
t622 = t510 * t437;
t390 = qJD(5) * t416 + t622;
t623 = t509 * t390;
t621 = t513 * t390;
t619 = -pkin(5) * t447 + qJD(5) * t652 - t654;
t618 = t366 - t442 * (pkin(5) * t510 - pkin(12) * t514);
t407 = pkin(4) * t445 - pkin(11) * t657;
t617 = t510 * t407 - t514 * t651;
t478 = -t507 * t514 + t510 * t631;
t570 = t502 * t591;
t614 = qJD(5) * t478 + t510 * t570 + t514 * t658;
t479 = t507 * t510 + t514 * t631;
t613 = qJD(5) * t479 - t510 * t658 + t514 * t570;
t607 = qJD(4) * t511;
t605 = qJD(5) * t509;
t603 = qJD(5) * t513;
t587 = t504 * t643;
t531 = t514 * t351 - t362 * t604 + t367 * t602 + t510 * t374;
t332 = pkin(12) * t438 + t531;
t560 = -t448 * t624 + t511 * t449 - t515 * t538;
t352 = t648 * qJD(4) + t560;
t339 = pkin(5) * t390 - pkin(12) * t389 + t352;
t578 = -t332 * t509 + t513 * t339;
t576 = t389 * t509 - t513 * t438;
t575 = t442 * t514;
t574 = t413 * t513;
t492 = -pkin(5) * t514 - pkin(12) * t510 - pkin(4);
t573 = pkin(12) * t445 - qJD(6) * t492 + t617;
t568 = t504 * t590;
t562 = qJD(3) * t587;
t463 = qJD(2) * t527;
t559 = t503 * t607 - t463;
t403 = t445 * t509 + t513 * t636;
t558 = t513 * t602 - t403;
t385 = pkin(4) * t475 - t647;
t419 = t455 * t510 + t514 * t475;
t420 = t455 * t514 - t475 * t510;
t358 = pkin(5) * t419 - pkin(12) * t420 + t385;
t557 = -pkin(12) * t447 - qJD(6) * t358 + t653;
t354 = pkin(12) * t454 + t652;
t394 = -qJD(5) * t419 + t446 * t514;
t395 = qJD(5) * t420 + t446 * t510;
t556 = -pkin(5) * t395 + pkin(12) * t394 + qJD(6) * t354 - t616;
t555 = t332 * t513 + t339 * t509;
t341 = pkin(12) * t442 + t343;
t361 = pkin(4) * t529 + t651;
t346 = t414 * pkin(5) - t416 * pkin(12) + t361;
t335 = t341 * t513 + t346 * t509;
t554 = t341 * t509 - t346 * t513;
t342 = -t362 * t510 + t367 * t514;
t397 = t453 * t515 + t511 * t546;
t418 = -t452 * t503 + t476 * t507;
t371 = t397 * t514 + t418 * t510;
t553 = t371 * t513 - t509 * t645;
t552 = -t371 * t509 - t513 * t645;
t550 = t381 * t514 - t386 * t510;
t370 = t397 * t510 - t418 * t514;
t399 = t420 * t513 + t454 * t509;
t398 = t420 * t509 - t513 * t454;
t547 = -t424 * t502 + t425 * t506;
t541 = qJD(6) * t630 + t614;
t537 = -t413 * t600 - t623;
t536 = -t413 * t601 + t621;
t535 = -pkin(11) * t438 + t361 * t442;
t528 = -qJD(6) * t479 + t559;
t464 = qJD(2) * t521;
t524 = t464 * t507 + t503 * t568;
t340 = -pkin(5) * t442 - t342;
t520 = -pkin(12) * t390 + (t340 + t342) * t413;
t517 = qJD(2) ^ 2;
t465 = qJD(2) * t522;
t461 = -t487 * t504 + t599;
t431 = -t464 * t503 + t507 * t568;
t402 = -t513 * t445 + t509 * t636;
t369 = t645 * qJD(4) + t465 * t515 + t524 * t511;
t368 = qJD(4) * t397 + t465 * t511 - t515 * t524;
t360 = qJD(6) * t399 + t394 * t509 - t513 * t447;
t359 = -qJD(6) * t398 + t394 * t513 + t447 * t509;
t357 = qJD(6) * t393 + t576;
t353 = -pkin(5) * t454 - t550;
t348 = -pkin(5) * t445 - t407 * t514 - t510 * t651;
t345 = -qJD(5) * t370 + t369 * t514 + t431 * t510;
t344 = qJD(5) * t371 + t369 * t510 - t431 * t514;
t331 = -qJD(6) * t335 + t578;
t330 = -qJD(6) * t554 + t555;
t1 = [(t424 * t464 + t425 * t465 + t448 * t452 + t449 * t453) * MDP(8) + (t368 * t529 + t418 * t438 - t431 * t657) * MDP(14) + (t369 * t529 + t418 * t437 + t431 * t445) * MDP(15) + (-t344 * t442 + t368 * t414 - t370 * t438 - t390 * t645) * MDP(21) + (-t345 * t442 + t368 * t416 - t371 * t438 - t389 * t645) * MDP(22) + ((-qJD(6) * t553 - t345 * t509 + t368 * t513) * t413 + t552 * t390 + t344 * t391 + t370 * t357) * MDP(28) + (-(qJD(6) * t552 + t345 * t513 + t368 * t509) * t413 - t553 * t390 + t344 * t393 + t370 * t356) * MDP(29) + (-MDP(4) * t516 + (-MDP(3) + (-MDP(5) * t506 + MDP(6) * t502) * t500) * t512) * t517 * t505 + ((-t464 * t502 + t465 * t506) * MDP(7) + (qJD(1) * t476 + t461) * MDP(8) * t626) * t609 + (MDP(5) * t464 - MDP(6) * t465) * t579; (t448 * t643 + (-t462 * t643 - t502 * t562) * qJD(2)) * MDP(5) + (-t449 * t643 + (t466 * t643 - t506 * t562) * qJD(2)) * MDP(6) + ((t449 * t506 - t635) * t504 + ((t462 * t502 - t466 * t506) * t504 + qJD(3) * t656) * qJD(2)) * MDP(7) + (t449 * t612 + t448 * t497 - t500 * pkin(2) * t566 - t425 * t466 - t424 * t462 + (-qJ(3) * t635 + qJD(3) * t547 - t461 * t594) * t504) * MDP(8) + (t437 * t455 + t445 * t446) * MDP(9) + (-t437 * t454 - t438 * t455 - t445 * t447 + t446 * t657) * MDP(10) + (t387 * t447 + t417 * t438 + t421 * t454 - t542 * t657) * MDP(14) + (t387 * t446 + t417 * t437 + t421 * t455 + t445 * t542) * MDP(15) + (t389 * t420 + t394 * t416) * MDP(16) + (-t389 * t419 - t390 * t420 - t394 * t414 - t395 * t416) * MDP(17) + (t389 * t454 + t394 * t442 + t416 * t447 + t420 * t438) * MDP(18) + (-t390 * t454 - t395 * t442 - t414 * t447 - t419 * t438) * MDP(19) + (t438 * t454 + t442 * t447) * MDP(20) + (t550 * t438 - t577 * t454 + t342 * t447 + t385 * t390 + t352 * t419 + t361 * t395 + t654 * t442 + t616 * t414 + (-t343 * t454 - t442 * t652) * qJD(5)) * MDP(21) + (-t343 * t447 + t352 * t420 + t361 * t394 + t385 * t389 + t616 * t416 - t438 * t652 + t442 * t653 - t531 * t454) * MDP(22) + (t356 * t399 + t359 * t393) * MDP(23) + (-t356 * t398 - t357 * t399 - t359 * t391 - t360 * t393) * MDP(24) + (t356 * t419 + t359 * t413 + t390 * t399 + t393 * t395) * MDP(25) + (-t357 * t419 - t360 * t413 - t390 * t398 - t391 * t395) * MDP(26) + (t390 * t419 + t395 * t413) * MDP(27) + ((-t354 * t509 + t358 * t513) * t390 + t331 * t419 - t554 * t395 + t353 * t357 + t333 * t398 + t340 * t360 + (t509 * t557 - t513 * t556) * t413 + t619 * t391) * MDP(28) + (-(t354 * t513 + t358 * t509) * t390 - t330 * t419 - t335 * t395 + t353 * t356 + t333 * t399 + t340 * t359 + (t509 * t556 + t513 * t557) * t413 + t619 * t393) * MDP(29) + (-t437 * MDP(11) + MDP(12) * t438 + MDP(14) * t352 + MDP(15) * t351) * t475 + (-t446 * MDP(11) + t616 * MDP(14) + MDP(15) * t660 + t650) * t529; (-t547 + t594) * MDP(8) * t609 + (t507 * t438 - t463 * t529 + (t529 * t607 + t592 * t657) * t503) * MDP(14) + (t507 * t437 - t467 * t529 + (-t445 * t592 + t529 * t606) * t503) * MDP(15) + (-t390 * t630 + t414 * t559 - t438 * t478 - t442 * t613) * MDP(21) + (-t389 * t630 + t416 * t559 - t438 * t479 + t442 * t614) * MDP(22) + ((-t479 * t509 - t513 * t630) * t390 + t478 * t357 + (t509 * t541 + t513 * t528) * t413 + t613 * t391) * MDP(28) + (-(t479 * t513 - t509 * t630) * t390 + t478 * t356 + (-t509 * t528 + t513 * t541) * t413 + t613 * t393) * MDP(29) + (-MDP(7) * t656 + (MDP(5) * t502 + MDP(6) * t506) * t587) * t517; -t657 ^ 2 * MDP(10) + (t529 * t657 + t437) * MDP(11) - qJD(2) * t650 + (t366 * t540 + (t366 - t648) * qJD(4) - t560) * MDP(14) + (-t387 * t657 - t540 * t651 - t646) * MDP(15) + (t389 * t510 + t416 * t575) * MDP(16) + ((t389 - t639) * t514 + (-t390 - t638) * t510) * MDP(17) + (t510 * t438 + t442 * t575) * MDP(18) + (-t442 ^ 2 * t510 + t514 * t438) * MDP(19) + (-pkin(4) * t390 - t366 * t414 + (-t352 + (-pkin(11) * qJD(5) - t407) * t442) * t514 + (-t442 * t651 + t535) * t510) * MDP(21) + (-pkin(4) * t389 + t352 * t510 - t366 * t416 + (pkin(11) * t604 + t617) * t442 + t535 * t514) * MDP(22) + (t356 * t510 * t513 + (-t510 * t601 + t558) * t393) * MDP(23) + (t391 * t403 + t393 * t402 + (-t391 * t513 - t393 * t509) * t602 + (-t642 - t357 * t513 + (t391 * t509 - t393 * t513) * qJD(6)) * t510) * MDP(24) + (-t356 * t514 + t558 * t413 + (t393 * t442 + t536) * t510) * MDP(25) + (t357 * t514 + (-t509 * t602 + t402) * t413 + (-t391 * t442 + t537) * t510) * MDP(26) + (t413 * t442 * t510 - t390 * t514) * MDP(27) + (t492 * t621 - t340 * t402 - t348 * t391 + (t509 * t573 - t513 * t618) * t413 + (t340 * t605 - t331 + (qJD(5) * t391 + t537) * pkin(11)) * t514 + (t340 * t600 + t333 * t509 - t442 * t554 + (t413 * t605 + t357) * pkin(11)) * t510) * MDP(28) + (-t492 * t623 - t340 * t403 - t348 * t393 + (t509 * t618 + t513 * t573) * t413 + (t340 * t603 + t330 + (qJD(5) * t393 - t536) * pkin(11)) * t514 + (-t340 * t601 + t333 * t513 - t442 * t335 + (t413 * t603 + t356) * pkin(11)) * t510) * MDP(29) + (-t657 * MDP(9) + (qJD(2) * t585 + qJD(4) - t488) * MDP(12) - t387 * MDP(14) - t416 * MDP(18) + t414 * MDP(19) - t442 * MDP(20) - t342 * MDP(21) + t343 * MDP(22) + MDP(10) * t445) * t445; -t414 ^ 2 * MDP(17) + (t389 + t639) * MDP(18) + (-t622 + t638) * MDP(19) + t438 * MDP(20) + (t343 * t442 + t644) * MDP(21) + (t342 * t442 + t361 * t414 - t531) * MDP(22) + (t393 * t574 + t642) * MDP(23) + ((t356 - t641) * t513 + (-t357 - t640) * t509) * MDP(24) + (t413 * t574 + t623) * MDP(25) + (-t413 ^ 2 * t509 + t621) * MDP(26) + (-pkin(5) * t357 - t343 * t391 + t520 * t509 - t513 * t649) * MDP(28) + (-pkin(5) * t356 - t343 * t393 + t509 * t649 + t520 * t513) * MDP(29) + (MDP(16) * t414 + MDP(17) * t416 - MDP(19) * qJD(5) - MDP(21) * t361 - MDP(25) * t393 + MDP(26) * t391 - MDP(27) * t413 + MDP(28) * t554 + MDP(29) * t335) * t416; t393 * t391 * MDP(23) + (-t391 ^ 2 + t393 ^ 2) * MDP(24) + (t597 + t641) * MDP(25) + (-t576 + t640) * MDP(26) + t390 * MDP(27) + (t335 * t413 - t340 * t393 + t578) * MDP(28) + (t340 * t391 - t413 * t554 - t555) * MDP(29) + (-MDP(25) * t637 - MDP(26) * t393 - MDP(28) * t335 + MDP(29) * t554) * qJD(6);];
tauc  = t1;
