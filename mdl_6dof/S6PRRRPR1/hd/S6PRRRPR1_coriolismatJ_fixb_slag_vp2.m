% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRRRPR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_coriolismatJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:00:26
% EndTime: 2019-03-08 23:00:49
% DurationCPUTime: 13.01s
% Computational Cost: add. (25781->552), mult. (56513->753), div. (0->0), fcn. (67157->12), ass. (0->328)
t355 = sin(qJ(4));
t356 = sin(qJ(3));
t359 = cos(qJ(4));
t360 = cos(qJ(3));
t325 = -t355 * t356 + t359 * t360;
t326 = -t355 * t360 - t359 * t356;
t351 = sin(pkin(12));
t549 = cos(pkin(12));
t399 = t351 * t325 - t326 * t549;
t653 = t399 * mrSges(6,3);
t342 = -pkin(3) * t360 - pkin(2);
t298 = -pkin(4) * t325 + t342;
t722 = m(6) * t298;
t353 = cos(pkin(6));
t352 = sin(pkin(6));
t357 = sin(qJ(2));
t505 = t352 * t357;
t312 = t353 * t360 - t356 * t505;
t313 = t353 * t356 + t360 * t505;
t262 = t312 * t355 + t313 * t359;
t441 = t359 * t312 - t313 * t355;
t175 = t262 * t351 - t441 * t549;
t658 = t262 * t549 + t351 * t441;
t610 = -pkin(9) - pkin(8);
t484 = t610 * t356;
t645 = t610 * t360;
t287 = -t355 * t484 + t359 * t645;
t385 = t325 * qJ(5) - t287;
t657 = t355 * t645 + t359 * t484;
t671 = t326 * qJ(5) + t657;
t704 = -t351 * t385 + t549 * t671;
t675 = t704 * t658;
t703 = t351 * t671 + t385 * t549;
t721 = t175 * t703 - t675;
t677 = t657 * mrSges(5,2);
t693 = t287 * mrSges(5,1);
t716 = t704 * mrSges(6,2);
t720 = -t677 / 0.2e1 + t693 / 0.2e1 - t716 / 0.2e1;
t354 = sin(qJ(6));
t358 = cos(qJ(6));
t432 = mrSges(7,1) * t358 - mrSges(7,2) * t354;
t710 = t703 * t432;
t717 = t703 * mrSges(6,1);
t719 = -t710 / 0.2e1 - t717 / 0.2e1;
t276 = t549 * t325 + t326 * t351;
t343 = Ifges(7,4) * t358;
t584 = Ifges(7,2) * t354;
t643 = t343 - t584;
t120 = Ifges(7,6) * t399 + t276 * t643;
t587 = Ifges(7,4) * t354;
t430 = Ifges(7,1) * t358 - t587;
t121 = Ifges(7,5) * t399 + t276 * t430;
t427 = Ifges(7,5) * t354 + Ifges(7,6) * t358;
t596 = t358 / 0.2e1;
t597 = t354 / 0.2e1;
t604 = t399 / 0.2e1;
t332 = Ifges(7,2) * t358 + t587;
t333 = Ifges(7,1) * t354 + t343;
t695 = -t354 / 0.2e1;
t660 = t332 * t695 + t333 * t596;
t387 = Ifges(5,5) * t325 + Ifges(5,6) * t326 - Ifges(6,6) * t399 + t120 * t596 + t121 * t597 + t427 * t604 + (Ifges(6,5) + t660) * t276;
t718 = t387 - t677 + t693 - t710 - t716 - t717;
t341 = pkin(3) * t359 + pkin(4);
t506 = t351 * t355;
t306 = -pkin(3) * t506 + t341 * t549;
t303 = -pkin(5) - t306;
t715 = t303 * t703;
t447 = t549 * t355;
t307 = pkin(3) * t447 + t351 * t341;
t714 = t307 * t704;
t470 = t549 * pkin(4);
t338 = -t470 - pkin(5);
t713 = t338 * t703;
t712 = t354 * t704;
t711 = t358 * t704;
t542 = t703 * t704;
t709 = qJD(3) + qJD(4);
t708 = t351 * t704 - t549 * t703;
t555 = t358 * mrSges(7,2);
t558 = t354 * mrSges(7,1);
t431 = t555 + t558;
t196 = t431 * t276;
t197 = t431 * t399;
t557 = t354 * mrSges(7,3);
t200 = -mrSges(7,2) * t399 - t276 * t557;
t202 = -mrSges(7,3) * t276 * t358 + mrSges(7,1) * t399;
t582 = Ifges(7,6) * t354;
t585 = Ifges(7,5) * t358;
t428 = -t582 + t585;
t271 = t399 * mrSges(6,1);
t446 = t276 * mrSges(6,2) + t271;
t349 = t358 ^ 2;
t654 = t349 / 0.2e1;
t664 = t276 * mrSges(6,3);
t161 = -pkin(5) * t276 - pkin(10) * t399 + t298;
t74 = t161 * t358 - t354 * t703;
t75 = t161 * t354 + t358 * t703;
t707 = (t342 * mrSges(5,2) + Ifges(5,4) * t325 + (-Ifges(5,1) + Ifges(5,2)) * t326) * t325 + t75 * t200 + t74 * t202 + t298 * t446 - (t196 + t664) * t704 + ((Ifges(6,1) + Ifges(7,1) * t654 + (-t343 + t584 / 0.2e1) * t354) * t276 - Ifges(6,4) * t399 + t428 * t604 + t120 * t695 + t121 * t596) * t399 - (-mrSges(6,3) * t704 + (Ifges(6,2) + Ifges(7,3)) * t399 + (-Ifges(6,4) + t428) * t276) * t276 + t703 * t197;
t559 = t349 * mrSges(7,3);
t347 = t354 ^ 2;
t560 = t347 * mrSges(7,3);
t639 = t560 / 0.2e1 + t559 / 0.2e1;
t706 = t639 * t276;
t406 = t175 * t431;
t410 = -t555 / 0.2e1 - t558 / 0.2e1;
t688 = t410 * t175;
t699 = -t406 / 0.2e1 - t688;
t705 = qJD(1) * t699;
t361 = cos(qJ(2));
t504 = t352 * t361;
t135 = -t354 * t504 + t358 * t658;
t473 = t557 / 0.2e1;
t474 = -t557 / 0.2e1;
t659 = t473 + t474;
t697 = -t688 + t406 / 0.2e1 + t659 * t135;
t545 = t135 * t358;
t134 = -t354 * t658 - t358 * t504;
t547 = t134 * t354;
t416 = t545 - t547;
t23 = m(7) * (-t416 + t658) * t175;
t700 = t23 * qJD(1);
t702 = t697 * qJD(6) + t700;
t701 = -qJD(6) * t699 - t700;
t345 = t356 * pkin(3);
t694 = m(5) * t345;
t692 = t175 * t276;
t691 = t175 * t307;
t315 = (t359 * t549 - t506) * pkin(3);
t686 = -t306 + t315;
t404 = t338 * t431;
t405 = t303 * t431;
t437 = t430 * t597 + t643 * t596 + t660;
t684 = t437 + t405 / 0.2e1 + t404 / 0.2e1;
t652 = t441 * mrSges(5,2);
t665 = t262 * mrSges(5,1);
t672 = t658 * t432;
t676 = t658 * mrSges(6,1);
t683 = -t652 - t665 - t672 - t676;
t682 = -t672 / 0.2e1 - t676 / 0.2e1;
t594 = pkin(4) * t326;
t302 = t345 - t594;
t424 = -t354 * t74 + t358 * t75;
t608 = -t197 / 0.2e1;
t619 = -m(7) / 0.2e1;
t621 = -m(6) / 0.2e1;
t518 = t399 * t354;
t201 = mrSges(7,2) * t276 - mrSges(7,3) * t518;
t500 = t358 * t201;
t517 = t399 * t358;
t203 = -mrSges(7,1) * t276 - mrSges(7,3) * t517;
t503 = t354 * t203;
t641 = t500 / 0.2e1 - t503 / 0.2e1;
t179 = pkin(5) * t399 - pkin(10) * t276 - t594;
t166 = t179 + t345;
t76 = t166 * t358 - t712;
t77 = t166 * t354 + t711;
t681 = t641 * t175 + (t134 * t76 + t135 * t77 - t175 * t424 + t721) * t619 + t658 * t608 + (t694 / 0.2e1 - t302 * t621) * t504;
t304 = pkin(10) + t307;
t496 = t347 + t349;
t444 = t496 * t304;
t680 = -t175 * t444 + t303 * t658;
t679 = -t202 / 0.2e1;
t678 = t262 / 0.2e1;
t662 = t432 * t399;
t462 = t662 / 0.2e1;
t673 = t658 * t399;
t668 = -t652 / 0.2e1 - t665 / 0.2e1;
t586 = Ifges(7,5) * t276;
t583 = Ifges(7,6) * t276;
t561 = t326 * mrSges(5,3);
t663 = t561 * t678;
t655 = -t175 * t196 / 0.2e1 - t135 * t200 / 0.2e1 + t134 * t679;
t435 = (t654 + t347 / 0.2e1) * mrSges(7,3);
t651 = -mrSges(6,1) - t432;
t617 = m(6) * pkin(4);
t492 = t617 / 0.2e1;
t648 = t351 * t492;
t642 = t639 * t399;
t603 = t304 / 0.2e1;
t81 = t179 * t354 + t711;
t638 = t200 * t603 + t81 * mrSges(7,3) / 0.2e1;
t613 = -mrSges(7,2) / 0.2e1;
t615 = mrSges(7,1) / 0.2e1;
t80 = t179 * t358 - t712;
t637 = t81 * t613 + t80 * t615;
t636 = t77 * t613 + t76 * t615;
t422 = -t80 * t354 + t81 * t358;
t635 = -t653 - t197;
t634 = -mrSges(4,1) * t360 + mrSges(4,2) * t356;
t472 = -t432 / 0.2e1 - mrSges(6,1) / 0.2e1;
t633 = -t662 / 0.2e1 + t706;
t631 = (-mrSges(5,1) * t355 - mrSges(5,2) * t359) * pkin(3);
t630 = m(7) * t338 - t549 * t617;
t292 = t326 * t504;
t293 = t325 * t504;
t219 = t351 * t292 + t293 * t549;
t209 = t219 * t358 + t354 * t505;
t524 = t209 * t358;
t208 = -t219 * t354 + t358 * t505;
t525 = t208 * t354;
t614 = -mrSges(6,2) / 0.2e1;
t629 = t219 * t614 + t292 * mrSges(5,1) / 0.2e1 - t293 * mrSges(5,2) / 0.2e1 + (-t525 / 0.2e1 + t524 / 0.2e1) * mrSges(7,3);
t618 = m(7) / 0.2e1;
t628 = t338 * t618 - t492 * t549;
t593 = pkin(4) * t351;
t337 = pkin(10) + t593;
t443 = t496 * t337;
t627 = m(7) * t443 + t351 * t617 + t559 + t560;
t510 = t307 * t399;
t511 = t306 * t276;
t514 = t303 * t196;
t314 = (t351 * t359 + t447) * pkin(3);
t538 = t704 * t314;
t601 = -t315 / 0.2e1;
t602 = t314 / 0.2e1;
t620 = m(6) / 0.2e1;
t626 = (t686 * t703 - t538 + t714) * t620 + (t304 * t422 + t315 * t424 - t538 + t715) * t618 + t514 / 0.2e1 + t197 * t602 + (-t276 * t601 + t399 * t602 - t510 / 0.2e1 - t511 / 0.2e1) * mrSges(6,3) + t720;
t402 = t361 * (-mrSges(5,1) * t326 + mrSges(5,2) * t325);
t420 = t361 * t446;
t625 = (t402 / 0.2e1 + t420 / 0.2e1) * t352 + t655;
t611 = mrSges(6,3) / 0.2e1;
t612 = -mrSges(6,3) / 0.2e1;
t624 = t673 * t612 + t692 * t611 - t655 + t663 - (t402 + t420) * t352 / 0.2e1;
t623 = 2 * qJD(3);
t622 = m(5) / 0.2e1;
t609 = t704 / 0.2e1;
t607 = t201 / 0.2e1;
t606 = -t203 / 0.2e1;
t605 = -t262 / 0.2e1;
t599 = t337 / 0.2e1;
t575 = t175 * mrSges(6,2);
t553 = t76 * t354;
t552 = t77 * t358;
t218 = -t292 * t549 + t293 * t351;
t540 = t704 * t218;
t531 = t175 * t218;
t530 = t175 * t399;
t528 = t175 * t314;
t377 = (t276 * t307 - t306 * t399) * t620 + (t276 * t444 + t303 * t399) * t618;
t408 = t200 * t597 + t202 * t596;
t380 = t302 * t620 + (t354 * t77 + t358 * t76) * t618 + t408;
t25 = t271 + t462 + (mrSges(6,2) - t435) * t276 - t377 + t380;
t521 = t25 * qJD(2);
t507 = t352 ^ 2 * t357;
t482 = t361 * t507;
t26 = m(7) * (t134 * t208 + t135 * t209 + t531) + m(6) * (t219 * t658 - t482 + t531) + m(5) * (t262 * t293 + t292 * t441 - t482) + m(4) * (-t507 + (-t312 * t356 + t313 * t360) * t352) * t361;
t520 = t26 * qJD(1);
t373 = (t276 * t443 + t338 * t399) * t618 + (t276 * t351 - t399 * t549) * t492;
t493 = -t617 / 0.2e1;
t379 = (t354 * t81 + t358 * t80) * t618 + t326 * t493 + t408;
t28 = -t373 + t379 + t446 + t462 - t706;
t516 = t28 * qJD(2);
t515 = t303 * t662;
t513 = t304 * t354;
t512 = t304 * t358;
t509 = t338 * t662;
t508 = t338 * t196;
t501 = t354 * t337;
t498 = t358 * t337;
t395 = t410 * t276;
t39 = t395 - t641;
t497 = t39 * qJD(2);
t495 = t356 ^ 2 + t360 ^ 2;
t494 = pkin(3) * t622;
t491 = mrSges(7,3) * t552;
t481 = t202 * t501;
t480 = Ifges(7,3) * t604;
t459 = -t513 / 0.2e1;
t450 = -t498 / 0.2e1;
t442 = t653 * t593;
t434 = mrSges(7,3) * t496 - mrSges(6,2);
t433 = -t342 * mrSges(5,1) - Ifges(5,4) * t326;
t413 = t524 - t525;
t365 = t472 * t218 + (-t218 * t306 + t219 * t307) * t620 + (t218 * t303 + t304 * t413) * t618 + (t292 * t359 + t293 * t355) * t494 + t629;
t2 = (t605 + t678) * t561 + t365 + t625 + t681;
t207 = -mrSges(6,1) * t276 + mrSges(6,2) * t399;
t279 = -mrSges(5,1) * t325 - mrSges(5,2) * t326;
t3 = t77 * t201 + t76 * t203 + t342 * t694 + m(7) * (t74 * t76 + t75 * t77 - t542) + (-mrSges(4,2) * pkin(2) + Ifges(4,4) * t360 + (Ifges(4,1) - Ifges(4,2)) * t356) * t360 + (-pkin(2) * mrSges(4,1) - Ifges(4,4) * t356 + pkin(3) * t279) * t356 + t433 * t326 + t707 + (t207 + t722) * t302;
t426 = -t2 * qJD(1) + t3 * qJD(2);
t5 = t81 * t201 + t80 * t203 + (-pkin(4) * t207 + t433) * t326 - t594 * t722 + m(7) * (t74 * t80 + t75 * t81 - t542) + t707;
t363 = (t504 * t594 + t675 + t721) * t621 + (t134 * t80 + t135 * t81 + t721) * t619 + (t608 - t653 / 0.2e1) * t658 - (t703 * t621 + t424 * t619 - t664 / 0.2e1 - t641) * t175;
t367 = t219 * t648 + t337 * t413 * t618 + (t628 + t472) * t218 + t629;
t7 = mrSges(6,3) * (t673 / 0.2e1 - t692 / 0.2e1) + t363 + t367 + t625;
t425 = -t7 * qJD(1) + t5 * qJD(2);
t423 = t552 - t553;
t419 = t470 * t664;
t11 = -t704 * t662 + t74 * t201 - t75 * t203 + ((-mrSges(7,3) * t75 - Ifges(7,4) * t517 + t583) * t358 + (t74 * mrSges(7,3) + t586 + (t587 + (-Ifges(7,1) + Ifges(7,2)) * t358) * t399) * t354) * t399;
t375 = (-t545 / 0.2e1 + t547 / 0.2e1) * t399 * mrSges(7,3) + t134 * t607 + t135 * t606 + t175 * t462;
t411 = t208 * t615 + t209 * t613;
t16 = t375 - t411;
t418 = t16 * qJD(1) + t11 * qJD(2);
t20 = (-(-m(7) - m(6)) * t704 + t635) * t399 + (-m(6) * t703 - m(7) * t424 - t500 + t503 - t664) * t276;
t374 = (t276 * t658 + t530) * t620 + (t276 * t416 + t530) * t618;
t383 = (t208 * t358 + t209 * t354) * t618 + t505 * t620;
t31 = t374 - t383;
t417 = -qJD(1) * t31 + qJD(2) * t20;
t412 = t399 * t435;
t394 = t410 * t315;
t368 = t628 * t658 - (t443 * t618 + t614 + t639 + t648) * t175 + t668 + t682;
t369 = (t658 * t686 + t528 - t691) * t620 + (t315 * t416 + t528 + t680) * t618 + t668;
t10 = t575 / 0.2e1 - t368 + t369 - t639 * t175 + t682;
t52 = t651 * t314 + t631 + t434 * t315 + m(7) * (t303 * t314 + t315 * t444) + m(6) * (-t306 * t314 + t307 * t315);
t362 = (t337 * t423 + t713) * t619 - t508 / 0.2e1 + t708 * t493 + t481 / 0.2e1 + t200 * t450 + t76 * t473 - t491 / 0.2e1 + t442 / 0.2e1 + t419 / 0.2e1 - t719 - t720;
t8 = t202 * t459 + t641 * t315 + t638 * t358 + t80 * t474 + t362 + t626 + t719;
t391 = t10 * qJD(1) + t8 * qJD(2) + t52 * qJD(3);
t372 = t586 + (0.5e1 / 0.4e1 * t587 + t332 / 0.4e1 + (-Ifges(7,1) / 0.2e1 + Ifges(7,2) / 0.4e1) * t358) * t399 + mrSges(7,2) * t609;
t376 = -t583 + (t343 / 0.4e1 + t333 / 0.4e1 + (Ifges(7,1) / 0.4e1 - Ifges(7,2) / 0.2e1) * t354) * t399 + mrSges(7,1) * t609;
t12 = t480 - t515 / 0.2e1 + t304 * t412 + (t201 * t603 + t376) * t354 + (t203 * t603 + t372) * t358 + t636;
t186 = t405 + t437;
t390 = -t12 * qJD(2) + t186 * qJD(3) - t705;
t103 = t394 - t684;
t14 = t480 - t509 / 0.2e1 + t337 * t412 + (t201 * t599 + t376) * t354 + (t203 * t599 + t372) * t358 + t637;
t206 = t404 + t437;
t384 = t14 * qJD(2) + t103 * qJD(3) - t206 * qJD(4) + t705;
t371 = -t354 * (t399 * t643 - t583) / 0.4e1 - t704 * t431 / 0.2e1 - t332 * t517 / 0.2e1 + t480 + (0.2e1 * t430 * t399 - t586) * t358 / 0.4e1 + t659 * t75 - (0.2e1 * t333 + t643) * t518 / 0.4e1 + (-t428 / 0.4e1 - t582 / 0.2e1 + t585 / 0.2e1) * t276;
t104 = t394 + t684;
t40 = t395 + t641;
t30 = t374 + t383;
t29 = t373 + t379 + t633;
t27 = t377 + t380 + t633;
t17 = t375 + t411;
t15 = t371 + t509 / 0.2e1 + t203 * t450 - t201 * t501 / 0.2e1 - t642 * t337 + t637;
t13 = t371 + t515 / 0.2e1 + t512 * t606 + t201 * t459 - t642 * t304 + t636;
t9 = -(t614 + t435) * t175 + t472 * t658 + t368 + t369;
t6 = t561 * t605 - t363 + t367 + t624;
t4 = -t362 + (t203 * t601 + t304 * t679 - t80 * mrSges(7,3) / 0.2e1) * t354 + (t315 * t607 + t638) * t358 + t472 * t703 + t387 + t626;
t1 = -t663 + t692 * t612 + t673 * t611 + t365 + (-t356 * mrSges(4,1) - t360 * mrSges(4,2)) * t504 + t624 - t681;
t18 = [t26 * qJD(2) + t23 * t709, t1 * qJD(3) + t6 * qJD(4) + t30 * qJD(5) + t17 * qJD(6) + t520 + (t209 * t201 + t208 * t203 + t219 * t664 + t292 * t561 + t293 * t325 * mrSges(5,3) + 0.2e1 * (t208 * t74 + t209 * t75 - t540) * t618 + 0.2e1 * (t219 * t703 + t298 * t505 - t540) * t620 + 0.2e1 * (-t287 * t293 + t292 * t657 + t342 * t505) * t622 + ((mrSges(4,3) * t495 - mrSges(3,2)) * t361 + (-mrSges(3,1) + t207 + t279 + t634) * t357 + m(4) * (pkin(8) * t361 * t495 - pkin(2) * t357)) * t352 - t635 * t218) * qJD(2), t1 * qJD(2) + t9 * qJD(4) + (t680 * t618 + (-t306 * t658 - t691) * t620 + (-t262 * t359 + t355 * t441) * t494) * t623 + (-t313 * mrSges(4,1) - t312 * mrSges(4,2) - t175 * t434 + t683) * qJD(3) + t702, t6 * qJD(2) + t9 * qJD(3) + (-t175 * t627 + t630 * t658 + t575 + t683) * qJD(4) + t702, qJD(2) * t30, t17 * qJD(2) + (-mrSges(7,1) * t135 - mrSges(7,2) * t134) * qJD(6) + t709 * t697; -qJD(3) * t2 - qJD(4) * t7 + qJD(5) * t31 + qJD(6) * t16 - t520, qJD(3) * t3 + qJD(4) * t5 - qJD(5) * t20 + qJD(6) * t11, t4 * qJD(4) + t27 * qJD(5) + t13 * qJD(6) + ((t304 * t423 + t715) * t618 + (-t306 * t703 + t714) * t620) * t623 + t426 + (-mrSges(7,3) * t553 + Ifges(4,5) * t360 - Ifges(4,6) * t356 + t200 * t512 - t202 * t513 + t491 + t514 + (m(5) * (t287 * t359 + t355 * t657) + (-t325 * t359 + t326 * t355) * mrSges(5,3)) * pkin(3) + t634 * pkin(8) + (-t510 - t511) * mrSges(6,3) + t718) * qJD(3), t4 * qJD(3) + (t508 + t200 * t498 - t481 + t708 * t617 + m(7) * (t337 * t422 + t713) - t442 - t419 + t422 * mrSges(7,3) + t718) * qJD(4) + t29 * qJD(5) + t15 * qJD(6) + t425, qJD(3) * t27 + qJD(4) * t29 + qJD(6) * t40 - t417, t13 * qJD(3) + t15 * qJD(4) + t40 * qJD(5) + (-mrSges(7,1) * t75 - mrSges(7,2) * t74 - t399 * t427) * qJD(6) + t418; qJD(2) * t2 + qJD(4) * t10 + t701, qJD(4) * t8 - qJD(5) * t25 - qJD(6) * t12 - t426, qJD(4) * t52 + qJD(6) * t186, t104 * qJD(6) + t391 + (t631 + (t630 + t651) * t314 + (-mrSges(6,2) + t627) * t315) * qJD(4), -t521, t104 * qJD(4) + (-t304 * t432 + t428) * qJD(6) + t390; qJD(2) * t7 - qJD(3) * t10 + t701, -qJD(3) * t8 - qJD(5) * t28 - qJD(6) * t14 - t425, -qJD(6) * t103 - t391, t206 * qJD(6), -t516 (-t337 * t432 + t428) * qJD(6) - t384; -qJD(2) * t31, qJD(3) * t25 + qJD(4) * t28 - qJD(6) * t39 + t417, t521, t516, 0, -qJD(6) * t431 - t497; -t16 * qJD(2) + t699 * t709, qJD(3) * t12 + qJD(4) * t14 + qJD(5) * t39 - t418, qJD(4) * t103 - t390, t384, t497, 0;];
Cq  = t18;
