% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRRPRR4
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
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRRPRR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR4_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:10:25
% EndTime: 2019-03-08 22:10:44
% DurationCPUTime: 10.57s
% Computational Cost: add. (13267->589), mult. (29469->767), div. (0->0), fcn. (30826->10), ass. (0->333)
t387 = sin(qJ(3));
t390 = cos(qJ(3));
t383 = sin(pkin(6));
t388 = sin(qJ(2));
t538 = t383 * t388;
t566 = cos(pkin(6));
t303 = t387 * t566 + t390 * t538;
t386 = sin(qJ(5));
t462 = t387 * t538 - t566 * t390;
t600 = cos(qJ(5));
t201 = t303 * t600 + t386 * t462;
t385 = sin(qJ(6));
t389 = cos(qJ(6));
t391 = cos(qJ(2));
t537 = t383 * t391;
t141 = t201 * t389 + t385 * t537;
t562 = t141 * t389;
t140 = -t201 * t385 + t389 * t537;
t563 = t140 * t385;
t450 = t562 - t563;
t636 = t303 * t386 - t462 * t600;
t709 = (t201 - t450) * m(7) * t636;
t554 = t709 * qJD(1);
t569 = t389 * mrSges(7,1);
t590 = mrSges(7,2) * t385;
t349 = -t569 + t590;
t710 = t349 - mrSges(6,1);
t493 = t600 * t387;
t514 = t390 * t537;
t275 = t386 * t514 - t493 * t537;
t332 = t386 * t387 + t390 * t600;
t276 = t332 * t537;
t392 = -pkin(3) - pkin(4);
t344 = -t386 * qJ(4) + t392 * t600;
t336 = pkin(5) - t344;
t345 = qJ(4) * t600 + t386 * t392;
t337 = -pkin(10) + t345;
t215 = t276 * t389 - t385 * t538;
t555 = t215 * t389;
t214 = -t276 * t385 - t389 * t538;
t556 = t214 * t385;
t448 = t555 - t556;
t510 = t276 * mrSges(6,2) / 0.2e1;
t625 = m(7) / 0.2e1;
t627 = m(6) / 0.2e1;
t638 = (t349 / 0.2e1 - mrSges(6,1) / 0.2e1) * t275;
t646 = (-t555 / 0.2e1 + t556 / 0.2e1) * mrSges(7,3) - t638;
t708 = (-t275 * t344 + t276 * t345) * t627 + (t275 * t336 + t337 * t448) * t625 + t510 + t646;
t333 = t386 * t390 - t493;
t240 = mrSges(6,1) * t332 - mrSges(6,2) * t333;
t649 = -t390 * pkin(3) - t387 * qJ(4);
t346 = -pkin(2) + t649;
t310 = t390 * pkin(4) - t346;
t705 = m(6) * t310 + t240;
t678 = pkin(8) - pkin(9);
t359 = t678 * t390;
t476 = t678 * t387;
t635 = t386 * t359 - t476 * t600;
t530 = t386 * t635;
t267 = t600 * t359 + t386 * t476;
t686 = t600 * t267;
t703 = t686 + t530;
t531 = t386 * t636;
t687 = t600 * t201;
t702 = t687 + t531;
t701 = t710 * t201;
t700 = t635 * mrSges(6,2) + t710 * t267;
t459 = mrSges(6,1) * t333 + mrSges(6,2) * t332;
t487 = t537 / 0.2e1;
t568 = t389 * mrSges(7,2);
t571 = t385 * mrSges(7,1);
t351 = t568 + t571;
t221 = t351 * t333;
t613 = -t221 / 0.2e1;
t699 = t201 * t613 - t459 * t487;
t587 = Ifges(7,4) * t389;
t356 = -Ifges(7,2) * t385 + t587;
t148 = Ifges(7,6) * t333 + t332 * t356;
t698 = t148 / 0.4e1;
t588 = Ifges(7,4) * t385;
t358 = Ifges(7,1) * t389 - t588;
t151 = Ifges(7,5) * t333 + t332 * t358;
t697 = t151 / 0.4e1;
t696 = pkin(5) * t267;
t544 = t332 * t385;
t228 = mrSges(7,2) * t333 + mrSges(7,3) * t544;
t668 = t389 * t228;
t695 = pkin(10) * t668;
t379 = t385 ^ 2;
t381 = t389 ^ 2;
t464 = mrSges(7,3) * (t381 / 0.2e1 + t379 / 0.2e1);
t692 = t201 * t336;
t673 = t201 * t635;
t690 = t267 * t336;
t688 = t267 * t635;
t330 = t386 * t349;
t519 = t379 + t381;
t434 = t519 * t600;
t416 = -t386 * mrSges(6,1) - t600 * mrSges(6,2) + mrSges(7,3) * t434;
t685 = t330 + t416;
t570 = t385 * mrSges(7,3);
t508 = -t570 / 0.2e1;
t684 = t570 / 0.2e1 + t508;
t172 = pkin(5) * t332 + pkin(10) * t333 + t310;
t103 = t172 * t389 - t267 * t385;
t104 = t172 * t385 + t267 * t389;
t543 = t332 * t389;
t231 = -mrSges(7,1) * t333 + mrSges(7,3) * t543;
t653 = t351 * t332;
t682 = t103 * t231 + t104 * t228 - t267 * t221 - t635 * t653;
t665 = t519 * t636;
t681 = pkin(5) * t201 + pkin(10) * t665;
t604 = t386 / 0.2e1;
t677 = -t231 / 0.2e1;
t680 = t385 * t386 * t677 + t668 * t604;
t602 = -t389 / 0.2e1;
t605 = t385 / 0.2e1;
t371 = Ifges(7,6) * t385;
t585 = Ifges(7,5) * t389;
t647 = t371 - t585;
t659 = -t332 / 0.2e1;
t679 = -t148 * t605 - t151 * t602 + (Ifges(7,3) + Ifges(6,2)) * t659 + (-t647 / 0.2e1 - Ifges(6,4)) * t333;
t612 = t228 / 0.2e1;
t676 = t332 / 0.2e1;
t675 = t333 / 0.2e1;
t672 = t385 * pkin(10) * t231;
t670 = t345 * t635;
t669 = t385 * t635;
t667 = t389 * t635;
t440 = t568 / 0.2e1 + t571 / 0.2e1;
t431 = t440 * t636;
t664 = t387 ^ 2 + t390 ^ 2;
t662 = 0.2e1 * t332;
t661 = Ifges(6,4) * t676;
t658 = mrSges(4,1) + mrSges(5,1);
t657 = Ifges(4,4) - Ifges(5,5);
t655 = mrSges(7,3) * t519;
t444 = t585 / 0.2e1 - t371 / 0.2e1;
t652 = (-t647 / 0.4e1 + t444) * t332;
t651 = t444 * t332;
t350 = -t390 * mrSges(5,1) - t387 * mrSges(5,3);
t650 = -t390 * mrSges(4,1) + t387 * mrSges(4,2) + t350;
t367 = t390 * qJ(4);
t648 = -pkin(3) * t387 + t367;
t522 = t684 * t141;
t645 = t684 * t104;
t243 = -pkin(5) * t333 + pkin(10) * t332;
t117 = t243 * t389 + t669;
t118 = t243 * t385 - t667;
t451 = -t117 * t385 + t118 * t389;
t616 = t636 / 0.2e1;
t617 = -t636 / 0.2e1;
t644 = t617 + t616;
t504 = t385 * t600;
t467 = -t504 / 0.2e1;
t503 = t389 * t600;
t623 = -mrSges(7,2) / 0.2e1;
t418 = mrSges(7,1) * t467 + t503 * t623;
t497 = t600 * t351;
t409 = t497 / 0.2e1 + t418;
t643 = qJD(6) * t409;
t642 = t409 * qJD(4);
t410 = -t497 / 0.2e1 + t418;
t641 = t410 * qJD(6);
t639 = t332 * t647 / 0.4e1 + Ifges(7,3) * t675;
t153 = t332 * Ifges(7,5) - t333 * t358;
t528 = t389 * t153;
t150 = t332 * Ifges(7,6) - t333 * t356;
t536 = t385 * t150;
t621 = Ifges(7,3) / 0.2e1;
t637 = -(t621 - Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t333 + Ifges(6,1) * t675 + t661 + t536 / 0.2e1 - t528 / 0.2e1;
t355 = Ifges(7,2) * t389 + t588;
t357 = Ifges(7,1) * t385 + t587;
t607 = t357 / 0.4e1;
t634 = t385 * (t356 / 0.4e1 + t607) + (t355 / 0.4e1 - t358 / 0.4e1) * t389;
t480 = t356 / 0.2e1 + t357 / 0.2e1;
t606 = -t358 / 0.2e1;
t482 = t355 / 0.2e1 + t606;
t633 = -t482 * t385 + t480 * t389;
t222 = t333 * t355;
t223 = t333 * t357;
t631 = t389 * t222 / 0.4e1 + t528 / 0.4e1 + t385 * t223 / 0.4e1 - t536 / 0.4e1 + t635 * t351 / 0.2e1;
t630 = 2 * qJD(3);
t629 = m(5) / 0.2e1;
t628 = -m(6) / 0.2e1;
t626 = -m(7) / 0.2e1;
t624 = -mrSges(7,1) / 0.2e1;
t622 = mrSges(7,2) / 0.2e1;
t620 = t117 / 0.2e1;
t619 = -t118 / 0.2e1;
t618 = t141 / 0.2e1;
t218 = t349 * t333;
t615 = -t218 / 0.2e1;
t614 = t653 / 0.2e1;
t609 = t336 / 0.2e1;
t608 = -t344 / 0.2e1;
t603 = -t387 / 0.2e1;
t601 = t389 / 0.2e1;
t599 = m(7) * (-t600 + t434) * t386;
t598 = m(7) * (-pkin(5) * t386 + pkin(10) * t434);
t596 = pkin(5) * t351;
t421 = -t140 * t504 + t141 * t503;
t34 = (-t687 + (0.1e1 - t519) * t531 + t421) * t625;
t594 = t34 * qJD(5);
t593 = mrSges(4,1) * t387;
t592 = mrSges(5,1) * t387;
t591 = mrSges(4,2) * t390;
t589 = mrSges(5,3) * t390;
t584 = Ifges(6,6) * t333;
t580 = t636 * mrSges(6,2);
t573 = t332 * mrSges(6,3);
t567 = t389 * mrSges(7,3);
t559 = t636 * t275;
t558 = t636 * t351;
t540 = t383 ^ 2 * t388;
t26 = m(7) * (t140 * t214 + t141 * t215 + t559) + m(6) * (t201 * t276 - t391 * t540 + t559) + (m(5) + m(4)) * (t303 * t514 + (t383 * t387 * t462 - t540) * t391);
t553 = t26 * qJD(1);
t550 = t635 * t275;
t353 = Ifges(7,5) * t385 + Ifges(7,6) * t389;
t542 = t333 * t353;
t541 = t336 * t351;
t229 = -mrSges(7,2) * t332 + t333 * t570;
t534 = t385 * t229;
t232 = mrSges(7,1) * t332 + t333 * t567;
t533 = t385 * t232;
t524 = t389 * t229;
t523 = t389 * t232;
t521 = t664 * pkin(8) * t537;
t516 = -t600 / 0.2e1;
t515 = qJD(4) * t599;
t506 = t567 / 0.2e1;
t505 = t345 * t600;
t498 = t600 * t275;
t488 = -t537 / 0.2e1;
t486 = t524 / 0.2e1;
t477 = m(5) * t346 + t350;
t473 = t519 * t344;
t468 = t387 * t488;
t329 = t387 * t392 + t367;
t179 = -t243 + t329;
t105 = t179 * t389 - t669;
t106 = t179 * t385 + t667;
t449 = -t267 * t636 - t673;
t453 = -t103 * t385 + t104 * t389;
t393 = (t329 * t537 + t449 + t673) * t627 + (t105 * t140 + t106 * t141 + t449) * t625 + t140 * t677 - t228 * t618 + t533 * t617 + t644 * t573 + (-m(5) * t648 - t589 + t591 + t592 + t593) * t488 + (t267 * t627 + t453 * t625 + t486 + t614) * t636 - t699;
t435 = t648 * t629;
t2 = t435 * t537 + t468 * t658 + t487 * t589 + t488 * t591 - t393 + t708;
t5 = t105 * t232 + t106 * t229 - t477 * t648 + m(7) * (t103 * t105 + t104 * t106 - t688) + (t310 * mrSges(6,1) - t679) * t333 + (t310 * mrSges(6,2) + (-Ifges(6,4) / 0.2e1 + t444) * t332 - t637) * t332 + (-pkin(2) * mrSges(4,2) - t346 * mrSges(5,3) + t657 * t390) * t390 + (-pkin(2) * mrSges(4,1) + t346 * mrSges(5,1) + (Ifges(4,1) + Ifges(5,1) - Ifges(4,2) - Ifges(5,3)) * t390 - t657 * t387) * t387 - t682 + t705 * t329;
t458 = -t2 * qJD(1) + t5 * qJD(2);
t457 = t558 / 0.2e1 + t522;
t13 = t103 * t229 - t104 * t232 + t635 * t218 + (t353 * t676 + t150 * t601 + t223 * t602 + t453 * mrSges(7,3) + (t222 + t153) * t605) * t333;
t420 = -t140 * t229 / 0.2e1 + t232 * t618 + t636 * t615;
t441 = t214 * mrSges(7,1) / 0.2e1 + t215 * t623;
t16 = (-t562 / 0.2e1 + t563 / 0.2e1) * t333 * mrSges(7,3) + t420 + t441;
t456 = -t16 * qJD(1) + t13 * qJD(2);
t41 = (t534 + t523 + m(7) * (t103 * t389 + t104 * t385) - t477 + t705) * t387;
t400 = (t386 * t276 - t498) * t627 + (t386 * t448 - t498) * t625;
t413 = m(7) * (t140 * t389 + t141 * t385) * t603 + m(6) * t468;
t43 = t400 + t413;
t455 = -qJD(1) * t43 + qJD(2) * t41;
t454 = t34 * qJD(4) + t554;
t452 = -t105 * t385 + t106 * t389;
t447 = t670 + t690;
t438 = -t534 / 0.2e1 - t523 / 0.2e1;
t404 = (t333 * t464 + t438) * t386 + t218 * t516;
t39 = (t590 / 0.2e1 - t569 / 0.2e1) * t387 + t404;
t446 = t39 * qJD(2) - qJD(3) * t410;
t445 = t598 / 0.2e1 + t330 / 0.2e1;
t443 = t105 * t624 + t106 * t622;
t442 = t117 * t624 + t118 * t622;
t439 = t345 * t613 - t609 * t653;
t437 = t355 * t605 + t357 * t602;
t436 = -Ifges(6,5) + t437;
t411 = (pkin(10) * t452 + t696) * t626 + pkin(5) * t614;
t425 = t232 * t608 + t337 * t677 + t697;
t426 = t698 + t337 * t612 + t344 * t229 / 0.2e1;
t4 = (pkin(10) * t612 - t148 / 0.4e1 + t426) * t389 + (pkin(10) * t677 - t151 / 0.4e1 + t425) * t385 + (t337 * t451 + t344 * t453 + t447) * t625 + ((t619 - t106 / 0.2e1) * t389 + (t620 + t105 / 0.2e1) * t385) * mrSges(7,3) + t411 + t439;
t417 = -t344 * mrSges(6,2) + mrSges(7,3) * t473 + t710 * t345;
t49 = -m(7) * (t336 * t345 + t337 * t473) + t417;
t405 = m(7) * (t692 + t450 * t344 + (-t337 * t519 + t345) * t636);
t419 = t681 * t625;
t9 = t644 * mrSges(6,2) + t419 - t405 / 0.2e1;
t433 = -t9 * qJD(1) + t4 * qJD(2) - t49 * qJD(3);
t412 = -t221 * t604 + t232 * t467 + t486 * t600;
t422 = -t103 * t504 + t104 * t503;
t21 = (-t686 + (t635 + t451) * t386 + t422) * t625 - t653 * t516 + t412 + t680;
t394 = (-t524 / 0.2e1 + t533 / 0.2e1 + (t267 - t453) * t625 - t653 / 0.2e1) * t636 + (t117 * t140 + t118 * t141 + t673) * t625 + t140 * t231 / 0.2e1 + t141 * t612 + t699;
t408 = (-pkin(5) * t275 + pkin(10) * t448) * t626 + t510;
t7 = t394 + t408 + t646;
t8 = t118 * t229 + t117 * t232 + m(7) * (t103 * t117 + t104 * t118 + t688) - t310 * t459 + t679 * t333 + (t661 - t651 + t637) * t332 + t682;
t432 = t7 * qJD(1) + t8 * qJD(2) + t21 * qJD(4);
t401 = t218 * t609 + t438 * t337 - t631 + t645;
t406 = t337 * t464 - t634;
t11 = -t652 + (-Ifges(7,3) / 0.2e1 + t406) * t333 + t401 + t443;
t31 = -t431 + t457;
t92 = -t541 + t633;
t428 = -t31 * qJD(1) + t11 * qJD(2) + t92 * qJD(3);
t395 = t703 * t627 + (t422 + t530) * t625 + t412;
t396 = t703 * t628 + (t386 * t452 + t686) * t626 + t600 * t614 + t680;
t19 = t395 + t396;
t398 = t702 * t628 + (t421 + t531) * t626;
t402 = t702 * t627 + (t386 * t665 + t687) * t625;
t28 = t398 + t402;
t423 = t434 * t337;
t66 = -m(7) * (t336 * t386 + t423) - m(6) * (-t344 * t386 + t505) - m(5) * qJ(4) - mrSges(5,3) + t685;
t427 = t28 * qJD(1) - t19 * qJD(2) + t66 * qJD(3);
t397 = -t330 / 0.2e1 + (-t505 + t423 + (t336 + t473) * t386) * t625;
t47 = -t397 + t416 + t445;
t415 = t34 * qJD(1) + t21 * qJD(2) - t47 * qJD(3) + t515;
t111 = t596 - t633;
t403 = pkin(5) * t615 + t438 * pkin(10) + t631 + t645;
t407 = pkin(10) * t464 + t634;
t14 = t652 + (t621 + t407) * t333 + t403 + t442;
t35 = (-t351 / 0.2e1 + t440) * t636 + t522;
t56 = (-pkin(5) / 0.2e1 - t336 / 0.2e1) * t351 + (mrSges(7,2) * t608 + t480) * t389 + (mrSges(7,1) * t608 - t482) * t385;
t414 = t35 * qJD(1) - t14 * qJD(2) + t56 * qJD(3) + t111 * qJD(5) + t642;
t57 = t596 / 0.2e1 + t356 * t602 + t541 / 0.2e1 + t385 * t606 - t440 * t344 + t437;
t54 = t397 + t445;
t42 = m(5) * t387 * t537 + t400 - t413;
t40 = t590 * t603 + t387 * t569 / 0.2e1 + t404;
t36 = t457 + t431;
t32 = -t558 / 0.2e1 - t431 + t522;
t27 = m(5) * t303 - t398 + t402;
t20 = t21 * qJD(5);
t18 = (m(5) * pkin(8) + mrSges(5,2)) * t390 + t395 - t396 + (-t386 * t333 + t516 * t662) * mrSges(6,3);
t17 = -t420 + t441 + (t140 * t508 + t141 * t506) * t333;
t15 = t407 * t333 - Ifges(7,5) * t543 / 0.2e1 + Ifges(7,6) * t544 / 0.2e1 + t403 - t442 - t639;
t12 = t406 * t333 + t401 - t443 + t639 + t651;
t10 = t405 / 0.2e1 - t580 / 0.2e1 + t419 + mrSges(6,2) * t617 + t636 * t464 + t616 * t655 - t701;
t6 = t214 * t508 + t215 * t506 + t394 - t408 + t638;
t3 = t447 * t625 - t411 + (t353 / 0.4e1 - Ifges(6,6) / 0.2e1) * t333 + t106 * t506 - t695 / 0.2e1 + t672 / 0.2e1 + t105 * t508 + t542 / 0.4e1 - t584 / 0.2e1 + Ifges(6,5) * t332 + t439 + (t607 * t662 + t698 + (t104 * t344 + t118 * t337) * t625 + mrSges(7,3) * t619 + t426) * t389 + (t355 * t659 + t697 + (-t103 * t344 - t117 * t337) * t625 + mrSges(7,3) * t620 + t425) * t385 - t700;
t1 = t393 + (t435 - t591 / 0.2e1 + t589 / 0.2e1 - t593 / 0.2e1 - t592 / 0.2e1) * t537 + t708;
t22 = [t26 * qJD(2) + (-qJD(3) + qJD(5)) * t709, t1 * qJD(3) + t42 * qJD(4) + t6 * qJD(5) + t17 * qJD(6) + t553 + (t214 * t232 + t215 * t229 - t276 * t573 + ((-mrSges(3,1) - t240 + t650) * t388 + (-mrSges(3,2) + (mrSges(5,2) + mrSges(4,3)) * t664) * t391) * t383 + 0.2e1 * (t103 * t214 + t104 * t215 + t550) * t625 + 0.2e1 * (t267 * t276 - t310 * t538 + t550) * t627 + 0.2e1 * (t346 * t538 + t521) * t629 + m(4) * (-pkin(2) * t538 + t521) + (-t333 * mrSges(6,3) - t221) * t275) * qJD(2), -t554 + t1 * qJD(2) + t27 * qJD(4) + t10 * qJD(5) + t32 * qJD(6) + ((t337 * t665 - t692) * t625 + (t201 * t344 + t345 * t636) * t627 + (-pkin(3) * t303 - qJ(4) * t462) * t629) * t630 + ((mrSges(6,2) - t655) * t636 + (mrSges(4,2) - mrSges(5,3)) * t462 - t658 * t303 + t701) * qJD(3), qJD(2) * t42 + qJD(3) * t27 + t594, t6 * qJD(2) + t10 * qJD(3) + (-m(7) * t681 - mrSges(7,3) * t665 + t580 + t701) * qJD(5) + t36 * qJD(6) + t454, t17 * qJD(2) + t32 * qJD(3) + t36 * qJD(5) + (-mrSges(7,1) * t141 - mrSges(7,2) * t140) * qJD(6); -qJD(3) * t2 - qJD(4) * t43 + qJD(5) * t7 - qJD(6) * t16 - t553, qJD(3) * t5 + qJD(4) * t41 + qJD(5) * t8 + qJD(6) * t13, t18 * qJD(4) + t3 * qJD(5) + t12 * qJD(6) + ((t337 * t452 - t690) * t625 + (t267 * t344 + t670) * t627) * t630 + t458 + (t336 * t653 + (-t337 * t228 - t106 * mrSges(7,3) - t148 / 0.2e1) * t389 + (t337 * t231 + t105 * mrSges(7,3) - t151 / 0.2e1) * t385 + (-pkin(3) * mrSges(5,2) + Ifges(5,4) + Ifges(4,5)) * t390 + (-qJ(4) * mrSges(5,2) - Ifges(4,6) + Ifges(5,6)) * t387 + (-t345 * mrSges(6,3) - t353 / 0.2e1 + Ifges(6,6)) * t333 + (-t344 * mrSges(6,3) + t436) * t332 + (m(5) * t649 + t650) * pkin(8) + t700) * qJD(3), qJD(3) * t18 + qJD(6) * t40 + t20 + t455, t3 * qJD(3) + t15 * qJD(6) + t432 + (t695 - t672 + m(7) * (pkin(10) * t451 - t696) - t151 * t605 - t148 * t601 + pkin(5) * t653 - t542 / 0.2e1 + t584 + t436 * t332 + t451 * mrSges(7,3) + t700) * qJD(5), t12 * qJD(3) + t40 * qJD(4) + t15 * qJD(5) + (-mrSges(7,1) * t104 - mrSges(7,2) * t103 + t542) * qJD(6) + t456; qJD(2) * t2 - qJD(4) * t28 - qJD(5) * t9 - qJD(6) * t31 + t554, qJD(4) * t19 + qJD(5) * t4 + qJD(6) * t11 - t458, -qJD(4) * t66 - qJD(5) * t49 + qJD(6) * t92, t54 * qJD(5) - t427 + t515 + t643, t54 * qJD(4) + (m(7) * (-pkin(5) * t345 + pkin(10) * t473) + t417) * qJD(5) + t57 * qJD(6) + t433, t642 + t57 * qJD(5) + (t337 * t349 + t647) * qJD(6) + t428; qJD(2) * t43 + qJD(3) * t28 + t594, -qJD(3) * t19 + qJD(6) * t39 + t20 - t455, -qJD(5) * t47 + t427 - t641, qJD(5) * t599 (t598 + t685) * qJD(5) + t641 + t415, qJD(5) * t410 + qJD(6) * t330 + t446; -qJD(2) * t7 + qJD(3) * t9 - qJD(6) * t35 - t454, -qJD(3) * t4 + qJD(6) * t14 - t432, qJD(4) * t47 - qJD(6) * t56 - t433, -t415 - t643, -t111 * qJD(6) (pkin(10) * t349 - t647) * qJD(6) - t414; t16 * qJD(2) + t31 * qJD(3) + t35 * qJD(5), -qJD(3) * t11 - qJD(4) * t39 - qJD(5) * t14 - t456, qJD(4) * t410 + qJD(5) * t56 - t428, qJD(5) * t409 - t446, t414, 0;];
Cq  = t22;
