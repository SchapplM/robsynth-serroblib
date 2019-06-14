% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-05-05 21:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRPRR13_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR13_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR13_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR13_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR13_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR13_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR13_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 20:54:45
% EndTime: 2019-05-05 20:54:55
% DurationCPUTime: 9.70s
% Computational Cost: add. (120738->337), mult. (382581->446), div. (0->0), fcn. (321699->14), ass. (0->163)
t676 = Ifges(4,1) + Ifges(5,2);
t667 = Ifges(4,5) - Ifges(5,4);
t675 = -Ifges(4,2) - Ifges(5,3);
t666 = Ifges(4,6) - Ifges(5,5);
t665 = -Ifges(5,6) - Ifges(4,4);
t674 = Ifges(4,3) + Ifges(5,1);
t603 = sin(pkin(6));
t652 = qJDD(1) * t603;
t604 = cos(pkin(12));
t602 = sin(pkin(7));
t664 = cos(pkin(6));
t640 = t664 * t602;
t663 = cos(pkin(7));
t643 = t603 * t663;
t584 = (t604 * t643 + t640) * qJD(1) * pkin(9);
t612 = qJD(1) ^ 2;
t608 = sin(qJ(1));
t611 = cos(qJ(1));
t634 = -g(1) * t611 - g(2) * t608;
t589 = -pkin(1) * t612 + qJ(2) * t652 + t634;
t601 = sin(pkin(12));
t669 = pkin(9) * t602;
t627 = -t604 * pkin(2) - t601 * t669;
t647 = t663 * pkin(9);
t654 = qJD(1) * t603;
t621 = qJD(1) * t627 * t654 + qJDD(1) * t647;
t645 = t608 * g(1) - g(2) * t611;
t588 = qJ(2) * t603 * t612 + qJDD(1) * pkin(1) + t645;
t642 = t604 * t664;
t646 = qJD(2) * t654;
t659 = t603 * t604;
t629 = -g(3) * t659 + t588 * t642 - 0.2e1 * t601 * t646;
t636 = qJDD(1) * t664;
t639 = qJD(1) * t664;
t522 = pkin(2) * t636 + t584 * t639 + (-t621 * t603 - t589) * t601 + t629;
t590 = (-pkin(9) * t601 * t643 + t664 * pkin(2)) * qJD(1);
t644 = t601 * t664;
t671 = 0.2e1 * t604;
t649 = t588 * t644 + t604 * t589 + t646 * t671;
t523 = t636 * t669 - t590 * t639 + (-g(3) * t601 + t621 * t604) * t603 + t649;
t631 = -t664 * g(3) + qJDD(2);
t533 = (-t588 + t627 * qJDD(1) + (-t584 * t604 + t590 * t601) * qJD(1)) * t603 + t631;
t607 = sin(qJ(3));
t641 = t607 * t663;
t661 = t602 * t607;
t670 = cos(qJ(3));
t496 = t522 * t641 + t670 * t523 + t533 * t661;
t635 = t663 * t670;
t660 = t603 * t601;
t672 = t607 * t660 - t635 * t659 - t670 * t640;
t569 = t672 * qJD(1);
t614 = t607 * t640 + (t670 * t601 + t604 * t641) * t603;
t570 = t614 * qJD(1);
t549 = pkin(3) * t569 - qJ(4) * t570;
t673 = -t602 * t659 + t664 * t663;
t585 = -t673 * qJD(1) - qJD(3);
t581 = t585 ^ 2;
t582 = t673 * qJDD(1) + qJDD(3);
t491 = pkin(3) * t581 - t582 * qJ(4) + 0.2e1 * qJD(4) * t585 + t569 * t549 - t496;
t668 = mrSges(4,1) - mrSges(5,2);
t662 = t569 * t585;
t648 = t602 * t670;
t495 = t522 * t635 - t607 * t523 + t533 * t648;
t493 = -t582 * pkin(3) - t581 * qJ(4) + t570 * t549 + qJDD(4) - t495;
t553 = -t569 * qJD(3) + t614 * qJDD(1);
t486 = (t569 * t570 - t582) * pkin(10) + (t553 - t662) * pkin(4) + t493;
t552 = qJD(3) * t570 + t672 * qJDD(1);
t563 = pkin(4) * t570 + pkin(10) * t585;
t568 = t569 ^ 2;
t504 = -t522 * t602 + t663 * t533;
t617 = (-t553 - t662) * qJ(4) + t504 + (-t585 * pkin(3) - 0.2e1 * qJD(4)) * t570;
t490 = -pkin(4) * t568 - t563 * t570 + (pkin(3) + pkin(10)) * t552 + t617;
t606 = sin(qJ(5));
t610 = cos(qJ(5));
t483 = t606 * t486 + t610 * t490;
t557 = t569 * t610 + t585 * t606;
t558 = t569 * t606 - t585 * t610;
t525 = -pkin(5) * t557 - pkin(11) * t558;
t548 = qJDD(5) + t553;
t567 = qJD(5) + t570;
t566 = t567 ^ 2;
t480 = -pkin(5) * t566 + pkin(11) * t548 + t525 * t557 + t483;
t488 = -pkin(4) * t552 - pkin(10) * t568 - t585 * t563 - t491;
t517 = -qJD(5) * t558 + t552 * t610 - t582 * t606;
t518 = qJD(5) * t557 + t552 * t606 + t582 * t610;
t484 = (-t557 * t567 - t518) * pkin(11) + (t558 * t567 - t517) * pkin(5) + t488;
t605 = sin(qJ(6));
t609 = cos(qJ(6));
t476 = -t480 * t605 + t484 * t609;
t531 = -t558 * t605 + t567 * t609;
t499 = qJD(6) * t531 + t518 * t609 + t548 * t605;
t532 = t558 * t609 + t567 * t605;
t506 = -mrSges(7,1) * t531 + mrSges(7,2) * t532;
t554 = qJD(6) - t557;
t507 = -mrSges(7,2) * t554 + mrSges(7,3) * t531;
t515 = qJDD(6) - t517;
t474 = m(7) * t476 + mrSges(7,1) * t515 - mrSges(7,3) * t499 - t506 * t532 + t507 * t554;
t477 = t480 * t609 + t484 * t605;
t498 = -qJD(6) * t532 - t518 * t605 + t548 * t609;
t508 = mrSges(7,1) * t554 - mrSges(7,3) * t532;
t475 = m(7) * t477 - mrSges(7,2) * t515 + mrSges(7,3) * t498 + t506 * t531 - t508 * t554;
t465 = t609 * t474 + t605 * t475;
t658 = t569 * t666 - t570 * t667 + t585 * t674;
t657 = t569 * t675 - t570 * t665 - t585 * t666;
t656 = -t665 * t569 - t570 * t676 + t667 * t585;
t560 = mrSges(5,1) * t569 + mrSges(5,3) * t585;
t655 = mrSges(4,2) * t585 - mrSges(4,3) * t569 - t560;
t550 = mrSges(4,1) * t569 + mrSges(4,2) * t570;
t524 = -mrSges(6,1) * t557 + mrSges(6,2) * t558;
t535 = mrSges(6,1) * t567 - mrSges(6,3) * t558;
t637 = -t474 * t605 + t609 * t475;
t463 = m(6) * t483 - mrSges(6,2) * t548 + mrSges(6,3) * t517 + t524 * t557 - t535 * t567 + t637;
t482 = t486 * t610 - t490 * t606;
t534 = -mrSges(6,2) * t567 + mrSges(6,3) * t557;
t479 = -pkin(5) * t548 - pkin(11) * t566 + t525 * t558 - t482;
t620 = -m(7) * t479 + t498 * mrSges(7,1) - mrSges(7,2) * t499 + t531 * t507 - t508 * t532;
t470 = m(6) * t482 + mrSges(6,1) * t548 - mrSges(6,3) * t518 - t524 * t558 + t534 * t567 + t620;
t457 = t463 * t606 + t470 * t610;
t551 = -mrSges(5,2) * t569 - mrSges(5,3) * t570;
t619 = -m(5) * t493 - t553 * mrSges(5,1) - t570 * t551 - t457;
t452 = m(4) * t495 - mrSges(4,3) * t553 - t550 * t570 + t668 * t582 - t655 * t585 + t619;
t562 = -mrSges(4,1) * t585 - mrSges(4,3) * t570;
t494 = pkin(3) * t552 + t617;
t561 = mrSges(5,1) * t570 - mrSges(5,2) * t585;
t638 = t610 * t463 - t470 * t606;
t625 = m(5) * t494 - t553 * mrSges(5,3) - t570 * t561 + t638;
t454 = m(4) * t504 + mrSges(4,2) * t553 + t668 * t552 + t562 * t570 + t655 * t569 + t625;
t618 = -m(6) * t488 + mrSges(6,1) * t517 - t518 * mrSges(6,2) + t534 * t557 - t558 * t535 - t465;
t615 = -m(5) * t491 + t582 * mrSges(5,3) - t585 * t561 - t618;
t461 = (-mrSges(4,3) - mrSges(5,1)) * t552 + (-t550 - t551) * t569 + m(4) * t496 - mrSges(4,2) * t582 + t562 * t585 + t615;
t444 = t452 * t648 + t663 * t454 + t461 * t661;
t447 = -t607 * t452 + t670 * t461;
t445 = t452 * t635 - t602 * t454 + t461 * t641;
t632 = -t604 * mrSges(3,1) + t601 * mrSges(3,2);
t624 = mrSges(3,1) * t664 - mrSges(3,3) * t660;
t623 = -mrSges(3,2) * t664 + mrSges(3,3) * t659;
t500 = Ifges(7,5) * t532 + Ifges(7,6) * t531 + Ifges(7,3) * t554;
t502 = Ifges(7,1) * t532 + Ifges(7,4) * t531 + Ifges(7,5) * t554;
t468 = -mrSges(7,1) * t479 + mrSges(7,3) * t477 + Ifges(7,4) * t499 + Ifges(7,2) * t498 + Ifges(7,6) * t515 - t500 * t532 + t502 * t554;
t501 = Ifges(7,4) * t532 + Ifges(7,2) * t531 + Ifges(7,6) * t554;
t469 = mrSges(7,2) * t479 - mrSges(7,3) * t476 + Ifges(7,1) * t499 + Ifges(7,4) * t498 + Ifges(7,5) * t515 + t500 * t531 - t501 * t554;
t510 = Ifges(6,4) * t558 + Ifges(6,2) * t557 + Ifges(6,6) * t567;
t511 = Ifges(6,1) * t558 + Ifges(6,4) * t557 + Ifges(6,5) * t567;
t616 = mrSges(6,1) * t482 - mrSges(6,2) * t483 + Ifges(6,5) * t518 + Ifges(6,6) * t517 + Ifges(6,3) * t548 + pkin(5) * t620 + pkin(11) * t637 + t609 * t468 + t605 * t469 + t558 * t510 - t557 * t511;
t613 = mrSges(7,1) * t476 - mrSges(7,2) * t477 + Ifges(7,5) * t499 + Ifges(7,6) * t498 + Ifges(7,3) * t515 + t501 * t532 - t502 * t531;
t592 = t623 * qJD(1);
t591 = t624 * qJD(1);
t587 = t632 * t654;
t571 = -t603 * t588 + t631;
t556 = -g(3) * t660 + t649;
t555 = -t601 * t589 + t629;
t509 = Ifges(6,5) * t558 + Ifges(6,6) * t557 + Ifges(6,3) * t567;
t456 = mrSges(5,2) * t582 - t560 * t585 - t619;
t455 = -mrSges(5,2) * t552 - t560 * t569 + t625;
t449 = -mrSges(6,1) * t488 + mrSges(6,3) * t483 + Ifges(6,4) * t518 + Ifges(6,2) * t517 + Ifges(6,6) * t548 - pkin(5) * t465 - t509 * t558 + t511 * t567 - t613;
t448 = mrSges(6,2) * t488 - mrSges(6,3) * t482 + Ifges(6,1) * t518 + Ifges(6,4) * t517 + Ifges(6,5) * t548 - pkin(11) * t465 - t468 * t605 + t469 * t609 + t509 * t557 - t510 * t567;
t446 = m(3) * t556 + t623 * qJDD(1) + (t587 * t659 - t664 * t591) * qJD(1) + t447;
t443 = m(3) * t571 + (t632 * qJDD(1) + (t591 * t601 - t592 * t604) * qJD(1)) * t603 + t444;
t442 = m(3) * t555 + t624 * qJDD(1) + (-t587 * t660 + t664 * t592) * qJD(1) + t445;
t441 = mrSges(5,1) * t493 + mrSges(4,2) * t504 - mrSges(4,3) * t495 - mrSges(5,3) * t494 + pkin(4) * t457 - qJ(4) * t455 + t665 * t552 + t553 * t676 + t658 * t569 + t667 * t582 + t657 * t585 + t616;
t440 = -mrSges(4,1) * t504 - mrSges(5,1) * t491 + mrSges(5,2) * t494 + mrSges(4,3) * t496 - pkin(3) * t455 - pkin(4) * t618 - pkin(10) * t638 - t606 * t448 - t610 * t449 + t552 * t675 - t665 * t553 + t658 * t570 + t666 * t582 + t656 * t585;
t439 = mrSges(4,1) * t495 - mrSges(4,2) * t496 + mrSges(5,2) * t493 - mrSges(5,3) * t491 + t610 * t448 - t606 * t449 - pkin(10) * t457 - pkin(3) * t456 + qJ(4) * t615 + t674 * t582 + t657 * t570 + (-qJ(4) * t551 - t656) * t569 + t667 * t553 + (-qJ(4) * mrSges(5,1) - t666) * t552;
t1 = [mrSges(2,1) * t645 - mrSges(2,2) * t634 + t664 * (mrSges(3,1) * t555 - mrSges(3,2) * t556 + Ifges(3,3) * t636 + pkin(2) * t445 + t663 * t439 + t440 * t648 + t441 * t661 + t447 * t669) + pkin(1) * (t442 * t642 + t446 * t644) + Ifges(2,3) * qJDD(1) + (t601 * (mrSges(3,2) * t571 - mrSges(3,3) * t555 - t607 * t440 + t670 * t441 - t444 * t669 - t445 * t647) + t604 * (-mrSges(3,1) * t571 + mrSges(3,3) * t556 - pkin(2) * t444 - t602 * t439 + t440 * t635 + t441 * t641 + t447 * t647) - pkin(1) * t443 + qJ(2) * (-t442 * t601 + t446 * t604) + 0.2e1 * (t601 * Ifges(3,5) + t604 * Ifges(3,6)) * t636 + (Ifges(3,2) * t604 ^ 2 + (Ifges(3,1) * t601 + Ifges(3,4) * t671) * t601) * t652) * t603; t443; t439; t456; t616; t613;];
tauJ  = t1;
