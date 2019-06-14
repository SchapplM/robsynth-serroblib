% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPPRR5
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-05-06 10:44
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPPRR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 10:38:48
% EndTime: 2019-05-06 10:38:55
% DurationCPUTime: 3.63s
% Computational Cost: add. (22764->328), mult. (52400->402), div. (0->0), fcn. (35515->10), ass. (0->142)
t665 = Ifges(3,1) + Ifges(4,1) + Ifges(5,1);
t642 = Ifges(3,4) - Ifges(4,5) - Ifges(5,4);
t664 = Ifges(5,5) - Ifges(3,5) - Ifges(4,4);
t640 = Ifges(3,6) + Ifges(5,6) - Ifges(4,6);
t663 = Ifges(4,3) + Ifges(5,2) + Ifges(3,2);
t662 = -Ifges(5,3) - Ifges(3,3) - Ifges(4,2);
t611 = qJD(1) ^ 2;
t606 = sin(qJ(1));
t610 = cos(qJ(1));
t626 = -g(1) * t610 - g(2) * t606;
t601 = sin(pkin(6));
t643 = qJDD(1) * t601;
t568 = -pkin(1) * t611 + pkin(8) * t643 + t626;
t605 = sin(qJ(2));
t609 = cos(qJ(2));
t649 = t601 * t609;
t631 = t606 * g(1) - g(2) * t610;
t656 = pkin(8) * t601;
t567 = qJDD(1) * pkin(1) + t611 * t656 + t631;
t602 = cos(pkin(6));
t652 = t567 * t602;
t515 = -g(3) * t649 - t605 * t568 + t609 * t652;
t645 = qJD(1) * t601;
t569 = (-t609 * pkin(2) - t605 * qJ(3)) * t645;
t593 = qJDD(1) * t602 + qJDD(2);
t634 = t605 * t645;
t594 = qJD(1) * t602 + qJD(2);
t660 = t594 ^ 2;
t498 = -t593 * pkin(2) - qJ(3) * t660 + t569 * t634 + qJDD(3) - t515;
t644 = qJD(1) * t609;
t633 = t601 * t644;
t566 = mrSges(4,2) * t633 + mrSges(4,3) * t594;
t661 = -m(4) * t498 + t593 * mrSges(4,1) + t594 * t566;
t659 = -2 * qJD(4);
t657 = pkin(2) * t594;
t655 = -mrSges(5,2) - mrSges(4,3);
t654 = mrSges(3,3) + mrSges(4,2);
t653 = qJ(4) * t609;
t651 = t601 ^ 2 * t611;
t650 = t601 * t605;
t574 = (qJD(2) * t644 + qJDD(1) * t605) * t601;
t575 = -qJD(2) * t634 + t609 * t643;
t560 = -pkin(3) * t594 - qJ(4) * t634;
t532 = -t602 * g(3) - t601 * t567;
t627 = t594 * t633;
t625 = -t575 * pkin(2) + t532 + (-t574 - t627) * qJ(3);
t628 = 0.2e1 * t634;
t638 = t609 ^ 2 * t651;
t615 = -qJ(4) * t638 + qJD(3) * t628 + t560 * t634 + qJDD(4) - t625;
t478 = -pkin(9) * t574 + (pkin(3) + pkin(4)) * t575 + (-pkin(9) * t609 + (-pkin(2) - pkin(4)) * t605) * t594 * t645 + t615;
t573 = (t609 * pkin(4) - t605 * pkin(9)) * t645;
t648 = t609 * t568 + t605 * t652;
t622 = -pkin(2) * t660 + t593 * qJ(3) + 0.2e1 * qJD(3) * t594 + t569 * t633 + t648;
t613 = -pkin(3) * t638 - t575 * qJ(4) + t594 * t560 + t633 * t659 + t622;
t480 = -pkin(4) * t660 - pkin(9) * t593 + (-g(3) * t605 - t573 * t644) * t601 + t613;
t604 = sin(qJ(5));
t608 = cos(qJ(5));
t475 = t604 * t478 + t608 * t480;
t546 = -t594 * t608 - t604 * t634;
t547 = -t594 * t604 + t608 * t634;
t518 = -pkin(5) * t546 - pkin(10) * t547;
t558 = qJDD(5) + t575;
t580 = qJD(5) + t633;
t578 = t580 ^ 2;
t472 = -pkin(5) * t578 + pkin(10) * t558 + t518 * t546 + t475;
t618 = -t574 * qJ(4) + t498 + (-t605 * t609 * t651 - t593) * pkin(3);
t482 = pkin(4) * t593 - pkin(9) * t660 - qJ(4) * t627 + qJD(4) * t628 + t573 * t634 - t618;
t513 = -qJD(5) * t547 - t574 * t604 - t593 * t608;
t514 = qJD(5) * t546 + t574 * t608 - t593 * t604;
t476 = t482 + (t547 * t580 - t513) * pkin(5) + (-t546 * t580 - t514) * pkin(10);
t603 = sin(qJ(6));
t607 = cos(qJ(6));
t469 = -t472 * t603 + t476 * t607;
t519 = -t547 * t603 + t580 * t607;
t489 = qJD(6) * t519 + t514 * t607 + t558 * t603;
t520 = t547 * t607 + t580 * t603;
t499 = -mrSges(7,1) * t519 + mrSges(7,2) * t520;
t544 = qJD(6) - t546;
t501 = -mrSges(7,2) * t544 + mrSges(7,3) * t519;
t510 = qJDD(6) - t513;
t466 = m(7) * t469 + mrSges(7,1) * t510 - mrSges(7,3) * t489 - t499 * t520 + t501 * t544;
t470 = t472 * t607 + t476 * t603;
t488 = -qJD(6) * t520 - t514 * t603 + t558 * t607;
t502 = mrSges(7,1) * t544 - mrSges(7,3) * t520;
t467 = m(7) * t470 - mrSges(7,2) * t510 + mrSges(7,3) * t488 + t499 * t519 - t502 * t544;
t457 = t607 * t466 + t603 * t467;
t561 = -mrSges(5,1) * t594 - mrSges(5,3) * t634;
t647 = -mrSges(3,1) * t594 + mrSges(3,3) * t634 + t561;
t564 = mrSges(5,2) * t594 - mrSges(5,3) * t633;
t646 = -mrSges(3,2) * t594 + mrSges(3,3) * t633 + t564;
t639 = g(3) * t650;
t637 = (t664 * t605 - t640 * t609) * t645 + t662 * t594;
t636 = (-t642 * t605 - t663 * t609) * t645 - t640 * t594;
t635 = (-t665 * t605 - t642 * t609) * t645 + t664 * t594;
t570 = (-t609 * mrSges(4,1) - t605 * mrSges(4,3)) * t645;
t632 = t570 * t645;
t517 = -mrSges(6,1) * t546 + mrSges(6,2) * t547;
t522 = mrSges(6,1) * t580 - mrSges(6,3) * t547;
t629 = -t466 * t603 + t607 * t467;
t455 = m(6) * t475 - mrSges(6,2) * t558 + mrSges(6,3) * t513 + t517 * t546 - t522 * t580 + t629;
t474 = t478 * t608 - t480 * t604;
t521 = -mrSges(6,2) * t580 + mrSges(6,3) * t546;
t471 = -pkin(5) * t558 - pkin(10) * t578 + t518 * t547 - t474;
t620 = -m(7) * t471 + t488 * mrSges(7,1) - mrSges(7,2) * t489 + t519 * t501 - t502 * t520;
t462 = m(6) * t474 + mrSges(6,1) * t558 - mrSges(6,3) * t514 - t517 * t547 + t521 * t580 + t620;
t630 = t608 * t455 - t462 * t604;
t450 = t604 * t455 + t608 * t462;
t484 = t613 - t639;
t624 = m(5) * t484 - t575 * mrSges(5,3) + t630;
t623 = -m(6) * t482 + t513 * mrSges(6,1) - t514 * mrSges(6,2) + t546 * t521 - t547 * t522 - t457;
t486 = pkin(3) * t575 - t634 * t657 + t615;
t621 = m(5) * t486 + t575 * mrSges(5,1) + t450;
t492 = t622 - t639;
t563 = -mrSges(4,1) * t594 + mrSges(4,2) * t634;
t619 = m(4) * t492 + t593 * mrSges(4,3) + t594 * t563 + t609 * t632 + t624;
t493 = (-0.2e1 * qJD(3) + t657) * t634 + t625;
t617 = m(4) * t493 - t575 * mrSges(4,1) - t621;
t485 = (t594 * t653 + t605 * t659) * t645 + t618;
t571 = (t609 * mrSges(5,1) + t605 * mrSges(5,2)) * t645;
t616 = -m(5) * t485 + t574 * mrSges(5,3) + t571 * t634 - t623;
t494 = Ifges(7,5) * t520 + Ifges(7,6) * t519 + Ifges(7,3) * t544;
t496 = Ifges(7,1) * t520 + Ifges(7,4) * t519 + Ifges(7,5) * t544;
t460 = -mrSges(7,1) * t471 + mrSges(7,3) * t470 + Ifges(7,4) * t489 + Ifges(7,2) * t488 + Ifges(7,6) * t510 - t494 * t520 + t496 * t544;
t495 = Ifges(7,4) * t520 + Ifges(7,2) * t519 + Ifges(7,6) * t544;
t461 = mrSges(7,2) * t471 - mrSges(7,3) * t469 + Ifges(7,1) * t489 + Ifges(7,4) * t488 + Ifges(7,5) * t510 + t494 * t519 - t495 * t544;
t504 = Ifges(6,4) * t547 + Ifges(6,2) * t546 + Ifges(6,6) * t580;
t505 = Ifges(6,1) * t547 + Ifges(6,4) * t546 + Ifges(6,5) * t580;
t614 = mrSges(6,1) * t474 - mrSges(6,2) * t475 + Ifges(6,5) * t514 + Ifges(6,6) * t513 + Ifges(6,3) * t558 + pkin(5) * t620 + pkin(10) * t629 + t607 * t460 + t603 * t461 + t547 * t504 - t546 * t505;
t453 = -t593 * mrSges(5,1) - t594 * t564 - t616;
t612 = mrSges(7,1) * t469 - mrSges(7,2) * t470 + Ifges(7,5) * t489 + Ifges(7,6) * t488 + Ifges(7,3) * t510 + t495 * t520 - t496 * t519;
t572 = (-t609 * mrSges(3,1) + t605 * mrSges(3,2)) * t645;
t516 = -t639 + t648;
t503 = Ifges(6,5) * t547 + Ifges(6,6) * t546 + Ifges(6,3) * t580;
t452 = t574 * mrSges(4,2) + t605 * t632 + t453 - t661;
t451 = t616 + t646 * t594 + (mrSges(3,1) + mrSges(5,1)) * t593 - t654 * t574 + m(3) * t515 + (-t570 - t572) * t634 + t661;
t449 = t574 * mrSges(5,2) + (t561 * t605 + t564 * t609) * t645 + t621;
t448 = t655 * t574 + ((-t564 - t566) * t609 + (-t561 - t563) * t605) * t645 + t617;
t447 = m(3) * t516 + t647 * t594 + (-mrSges(3,2) + mrSges(5,2)) * t593 + t654 * t575 + (-t571 + t572) * t633 + t619;
t446 = -mrSges(6,1) * t482 + mrSges(6,3) * t475 + Ifges(6,4) * t514 + Ifges(6,2) * t513 + Ifges(6,6) * t558 - pkin(5) * t457 - t503 * t547 + t505 * t580 - t612;
t445 = mrSges(6,2) * t482 - mrSges(6,3) * t474 + Ifges(6,1) * t514 + Ifges(6,4) * t513 + Ifges(6,5) * t558 - pkin(10) * t457 - t460 * t603 + t461 * t607 + t503 * t546 - t504 * t580;
t444 = -t608 * t446 + qJ(3) * (t561 * t594 + t619) - pkin(9) * t630 - t604 * t445 - pkin(4) * t623 + mrSges(3,1) * t515 - mrSges(3,2) * t516 + mrSges(4,3) * t492 - mrSges(4,1) * t498 + mrSges(5,2) * t484 - mrSges(5,1) * t485 - pkin(3) * t453 - pkin(2) * t452 + (qJ(3) * mrSges(5,2) - t662) * t593 + (qJ(3) * mrSges(4,2) + t640) * t575 - t664 * t574 + (-t636 * t605 + (-qJ(3) * t571 + t635) * t609) * t645;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t631 - mrSges(2,2) * t626 + (mrSges(3,2) * t532 + mrSges(4,2) * t498 + mrSges(5,2) * t486 - mrSges(3,3) * t515 - mrSges(4,3) * t493 - mrSges(5,3) * t485 - pkin(9) * t450 - qJ(3) * t448 - qJ(4) * t453 + t608 * t445 - t604 * t446 + t636 * t594 - t664 * t593 + t642 * t575 + t665 * t574 - t637 * t633) * t650 + ((-qJ(4) * mrSges(5,2) + t640) * t593 + t663 * t575 + t642 * t574 - qJ(4) * t624 + t614 + (-qJ(4) * t561 - t635) * t594 - mrSges(3,1) * t532 + mrSges(3,3) * t516 + mrSges(4,2) * t492 - mrSges(4,1) * t493 - mrSges(5,3) * t484 + mrSges(5,1) * t486 + pkin(4) * t450 + pkin(3) * t449 - pkin(2) * t448 + (t571 * t653 + t605 * t637) * t645) * t649 + t602 * t444 + pkin(1) * ((t447 * t605 + t451 * t609) * t602 + (-m(3) * t532 + t575 * mrSges(3,1) + (-mrSges(3,2) - t655) * t574 + ((t566 + t646) * t609 + (t563 + t647) * t605) * t645 - t617) * t601) + (t447 * t609 - t451 * t605) * t656; t444; t452; t449; t614; t612;];
tauJ  = t1;
