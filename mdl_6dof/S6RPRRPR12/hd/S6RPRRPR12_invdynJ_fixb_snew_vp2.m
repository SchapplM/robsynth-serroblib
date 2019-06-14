% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
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
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 00:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRPR12_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR12_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR12_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR12_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR12_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR12_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR12_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 00:42:49
% EndTime: 2019-05-06 00:43:00
% DurationCPUTime: 10.31s
% Computational Cost: add. (136915->339), mult. (427055->449), div. (0->0), fcn. (362228->14), ass. (0->163)
t653 = Ifges(5,1) + Ifges(6,2);
t642 = Ifges(5,4) + Ifges(6,6);
t641 = Ifges(5,5) - Ifges(6,4);
t652 = -Ifges(5,2) - Ifges(6,3);
t640 = Ifges(5,6) - Ifges(6,5);
t651 = Ifges(5,3) + Ifges(6,1);
t579 = sin(pkin(12));
t582 = cos(pkin(12));
t587 = sin(qJ(3));
t583 = cos(pkin(7));
t590 = cos(qJ(3));
t627 = t583 * t590;
t650 = -t579 * t587 + t582 * t627;
t581 = sin(pkin(6));
t584 = cos(pkin(6));
t628 = t583 * t587;
t580 = sin(pkin(7));
t634 = t580 * t587;
t597 = t584 * t634 + (t579 * t590 + t582 * t628) * t581;
t553 = t597 * qJD(1);
t633 = t580 * t590;
t619 = t584 * t633;
t539 = -t553 * qJD(3) + (t650 * t581 + t619) * qJDD(1);
t629 = t581 * t583;
t566 = (t580 * t584 + t582 * t629) * qJD(1) * pkin(9);
t592 = qJD(1) ^ 2;
t588 = sin(qJ(1));
t591 = cos(qJ(1));
t613 = -t591 * g(1) - t588 * g(2);
t632 = t581 * qJ(2);
t570 = -t592 * pkin(1) + qJDD(1) * t632 + t613;
t643 = pkin(9) * t580;
t609 = -t582 * pkin(2) - t579 * t643;
t622 = qJD(1) * t581;
t638 = pkin(9) * qJDD(1);
t604 = qJD(1) * t609 * t622 + t583 * t638;
t617 = qJD(2) * t622;
t630 = t581 * t582;
t616 = t588 * g(1) - t591 * g(2);
t569 = qJDD(1) * pkin(1) + t592 * t632 + t616;
t636 = t569 * t584;
t610 = -g(3) * t630 - 0.2e1 * t579 * t617 + t582 * t636;
t515 = (pkin(2) * qJDD(1) + qJD(1) * t566) * t584 + (-t604 * t581 - t570) * t579 + t610;
t644 = pkin(9) * t579;
t571 = (t584 * pkin(2) - t629 * t644) * qJD(1);
t647 = 0.2e1 * t582;
t618 = t582 * t570 + t579 * t636 + t617 * t647;
t516 = (-qJD(1) * t571 + t580 * t638) * t584 + (-g(3) * t579 + t604 * t582) * t581 + t618;
t615 = -t584 * g(3) + qJDD(2);
t525 = (-t569 + t609 * qJDD(1) + (-t566 * t582 + t571 * t579) * qJD(1)) * t581 + t615;
t478 = -t587 * t516 + (t515 * t583 + t525 * t580) * t590;
t605 = -t580 * t630 + t584 * t583;
t567 = t605 * qJD(1) + qJD(3);
t586 = sin(qJ(4));
t645 = cos(qJ(4));
t546 = t586 * t553 - t645 * t567;
t552 = qJD(1) * t619 + t650 * t622;
t551 = qJD(4) - t552;
t526 = t546 * mrSges(6,1) - t551 * mrSges(6,3);
t536 = qJDD(4) - t539;
t479 = t515 * t628 + t590 * t516 + t525 * t634;
t538 = -t552 * pkin(3) - t553 * pkin(10);
t563 = t567 ^ 2;
t564 = t605 * qJDD(1) + qJDD(3);
t475 = -t563 * pkin(3) + t564 * pkin(10) + t552 * t538 + t479;
t490 = -t580 * t515 + t583 * t525;
t540 = t552 * qJD(3) + t597 * qJDD(1);
t477 = (-t552 * t567 - t540) * pkin(10) + (t553 * t567 - t539) * pkin(3) + t490;
t470 = -t586 * t475 + t645 * t477;
t547 = t645 * t553 + t586 * t567;
t517 = t546 * pkin(4) - t547 * qJ(5);
t550 = t551 ^ 2;
t468 = -t536 * pkin(4) - t550 * qJ(5) + t547 * t517 + qJDD(5) - t470;
t510 = -t546 * qJD(4) + t645 * t540 + t586 * t564;
t637 = t546 * t551;
t463 = (t546 * t547 - t536) * pkin(11) + (t510 + t637) * pkin(5) + t468;
t509 = t547 * qJD(4) + t586 * t540 - t645 * t564;
t530 = t547 * pkin(5) - t551 * pkin(11);
t545 = t546 ^ 2;
t474 = -t564 * pkin(3) - t563 * pkin(10) + t553 * t538 - t478;
t646 = -2 * qJD(5);
t594 = (-t510 + t637) * qJ(5) + t474 + (t551 * pkin(4) + t646) * t547;
t466 = -t545 * pkin(5) - t547 * t530 + (pkin(4) + pkin(11)) * t509 + t594;
t585 = sin(qJ(6));
t589 = cos(qJ(6));
t461 = t589 * t463 - t585 * t466;
t523 = t589 * t546 - t585 * t551;
t485 = t523 * qJD(6) + t585 * t509 + t589 * t536;
t524 = t585 * t546 + t589 * t551;
t491 = -t523 * mrSges(7,1) + t524 * mrSges(7,2);
t542 = qJD(6) + t547;
t494 = -t542 * mrSges(7,2) + t523 * mrSges(7,3);
t506 = qJDD(6) + t510;
t458 = m(7) * t461 + t506 * mrSges(7,1) - t485 * mrSges(7,3) - t524 * t491 + t542 * t494;
t462 = t585 * t463 + t589 * t466;
t484 = -t524 * qJD(6) + t589 * t509 - t585 * t536;
t495 = t542 * mrSges(7,1) - t524 * mrSges(7,3);
t459 = m(7) * t462 - t506 * mrSges(7,2) + t484 * mrSges(7,3) + t523 * t491 - t542 * t495;
t449 = t589 * t458 + t585 * t459;
t519 = -t546 * mrSges(6,2) - t547 * mrSges(6,3);
t600 = -m(6) * t468 - t510 * mrSges(6,1) - t547 * t519 - t449;
t447 = t536 * mrSges(6,2) + t551 * t526 - t600;
t471 = t645 * t475 + t586 * t477;
t599 = -t550 * pkin(4) + t536 * qJ(5) - t546 * t517 + t471;
t465 = -t509 * pkin(5) - t545 * pkin(11) + ((2 * qJD(5)) + t530) * t551 + t599;
t486 = Ifges(7,5) * t524 + Ifges(7,6) * t523 + Ifges(7,3) * t542;
t488 = Ifges(7,1) * t524 + Ifges(7,4) * t523 + Ifges(7,5) * t542;
t450 = -mrSges(7,1) * t465 + mrSges(7,3) * t462 + Ifges(7,4) * t485 + Ifges(7,2) * t484 + Ifges(7,6) * t506 - t524 * t486 + t542 * t488;
t487 = Ifges(7,4) * t524 + Ifges(7,2) * t523 + Ifges(7,6) * t542;
t451 = mrSges(7,2) * t465 - mrSges(7,3) * t461 + Ifges(7,1) * t485 + Ifges(7,4) * t484 + Ifges(7,5) * t506 + t523 * t486 - t542 * t487;
t467 = t551 * t646 - t599;
t527 = t547 * mrSges(6,1) + t551 * mrSges(6,2);
t601 = -m(7) * t465 + t484 * mrSges(7,1) - t485 * mrSges(7,2) + t523 * t494 - t524 * t495;
t596 = -m(6) * t467 + t536 * mrSges(6,3) + t551 * t527 - t601;
t623 = -t642 * t546 + t653 * t547 + t641 * t551;
t624 = t652 * t546 + t642 * t547 + t640 * t551;
t648 = -t640 * t509 + t641 * t510 + t651 * t536 + t623 * t546 + t624 * t547 + mrSges(5,1) * t470 - mrSges(5,2) * t471 + mrSges(6,2) * t468 - mrSges(6,3) * t467 - pkin(4) * t447 - pkin(11) * t449 + qJ(5) * (-t509 * mrSges(6,1) - t546 * t519 + t596) - t585 * t450 + t589 * t451;
t631 = t581 * t579;
t518 = t546 * mrSges(5,1) + t547 * mrSges(5,2);
t528 = -t551 * mrSges(5,2) - t546 * mrSges(5,3);
t446 = m(5) * t470 - t510 * mrSges(5,3) - t547 * t518 + (-t526 + t528) * t551 + (mrSges(5,1) - mrSges(6,2)) * t536 + t600;
t529 = t551 * mrSges(5,1) - t547 * mrSges(5,3);
t454 = m(5) * t471 - t536 * mrSges(5,2) - t551 * t529 + (-t518 - t519) * t546 + (-mrSges(5,3) - mrSges(6,1)) * t509 + t596;
t441 = t645 * t446 + t586 * t454;
t626 = -t585 * t458 + t589 * t459;
t625 = t640 * t546 - t641 * t547 - t651 * t551;
t537 = -t552 * mrSges(4,1) + t553 * mrSges(4,2);
t549 = t567 * mrSges(4,1) - t553 * mrSges(4,3);
t614 = -t586 * t446 + t645 * t454;
t438 = m(4) * t479 - t564 * mrSges(4,2) + t539 * mrSges(4,3) + t552 * t537 - t567 * t549 + t614;
t548 = -t567 * mrSges(4,2) + t552 * mrSges(4,3);
t440 = m(4) * t490 - t539 * mrSges(4,1) + t540 * mrSges(4,2) - t552 * t548 + t553 * t549 + t441;
t469 = t509 * pkin(4) + t594;
t603 = -m(6) * t469 + t509 * mrSges(6,2) + t546 * t526 - t626;
t595 = -m(5) * t474 - t509 * mrSges(5,1) - t546 * t528 + (t527 - t529) * t547 + (-mrSges(5,2) + mrSges(6,3)) * t510 + t603;
t444 = m(4) * t478 + t564 * mrSges(4,1) - t540 * mrSges(4,3) - t553 * t537 + t567 * t548 + t595;
t429 = t438 * t634 + t583 * t440 + t444 * t633;
t434 = t590 * t438 - t587 * t444;
t430 = t438 * t628 - t580 * t440 + t444 * t627;
t612 = -t582 * mrSges(3,1) + t579 * mrSges(3,2);
t608 = t584 * mrSges(3,1) - mrSges(3,3) * t631;
t607 = -t584 * mrSges(3,2) + mrSges(3,3) * t630;
t448 = -t510 * mrSges(6,3) - t547 * t527 - t603;
t431 = -mrSges(5,1) * t474 - mrSges(6,1) * t467 + mrSges(6,2) * t469 + mrSges(5,3) * t471 - pkin(4) * t448 - pkin(5) * t601 - pkin(11) * t626 - t589 * t450 - t585 * t451 + t652 * t509 + t642 * t510 + t640 * t536 + t625 * t547 + t623 * t551;
t598 = mrSges(7,1) * t461 - mrSges(7,2) * t462 + Ifges(7,5) * t485 + Ifges(7,6) * t484 + Ifges(7,3) * t506 + t524 * t487 - t523 * t488;
t432 = mrSges(6,1) * t468 + mrSges(5,2) * t474 - mrSges(5,3) * t470 - mrSges(6,3) * t469 + pkin(5) * t449 - qJ(5) * t448 - t642 * t509 + t653 * t510 + t641 * t536 + t625 * t546 - t624 * t551 + t598;
t532 = Ifges(4,5) * t553 + Ifges(4,6) * t552 + Ifges(4,3) * t567;
t533 = Ifges(4,4) * t553 + Ifges(4,2) * t552 + Ifges(4,6) * t567;
t425 = mrSges(4,2) * t490 - mrSges(4,3) * t478 + Ifges(4,1) * t540 + Ifges(4,4) * t539 + Ifges(4,5) * t564 - pkin(10) * t441 - t586 * t431 + t645 * t432 + t552 * t532 - t567 * t533;
t534 = Ifges(4,1) * t553 + Ifges(4,4) * t552 + Ifges(4,5) * t567;
t426 = -mrSges(4,1) * t490 + mrSges(4,3) * t479 + Ifges(4,4) * t540 + Ifges(4,2) * t539 + Ifges(4,6) * t564 - pkin(3) * t441 - t553 * t532 + t567 * t534 - t648;
t602 = pkin(9) * t434 + t425 * t587 + t426 * t590;
t573 = t607 * qJD(1);
t572 = t608 * qJD(1);
t568 = t612 * t622;
t554 = -t581 * t569 + t615;
t544 = -g(3) * t631 + t618;
t543 = -t579 * t570 + t610;
t433 = m(3) * t544 + t607 * qJDD(1) + (t568 * t630 - t572 * t584) * qJD(1) + t434;
t428 = m(3) * t554 + (t612 * qJDD(1) + (t572 * t579 - t573 * t582) * qJD(1)) * t581 + t429;
t427 = m(3) * t543 + t608 * qJDD(1) + (-t568 * t631 + t573 * t584) * qJD(1) + t430;
t424 = mrSges(4,1) * t478 - mrSges(4,2) * t479 + Ifges(4,5) * t540 + Ifges(4,6) * t539 + Ifges(4,3) * t564 + pkin(3) * t595 + pkin(10) * t614 + t645 * t431 + t586 * t432 + t553 * t533 - t552 * t534;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t616 - mrSges(2,2) * t613 + (mrSges(3,1) * t543 - mrSges(3,2) * t544 + pkin(2) * t430 + t583 * t424 + pkin(1) * (t427 * t582 + t433 * t579) + Ifges(3,3) * t584 * qJDD(1) + t602 * t580) * t584 + (t579 * (mrSges(3,2) * t554 - mrSges(3,3) * t543 + t590 * t425 - t587 * t426 - t429 * t643) + t582 * (-mrSges(3,1) * t554 + mrSges(3,3) * t544 - pkin(2) * t429 - t580 * t424) - pkin(1) * t428 + qJ(2) * (-t579 * t427 + t582 * t433) + (-t430 * t644 + t582 * t602) * t583 + ((Ifges(3,2) * t582 ^ 2 + (Ifges(3,1) * t579 + Ifges(3,4) * t647) * t579) * t581 + 0.2e1 * t584 * (Ifges(3,5) * t579 + Ifges(3,6) * t582)) * qJDD(1)) * t581; t428; t424; t648; t447; t598;];
tauJ  = t1;
