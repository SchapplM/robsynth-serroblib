% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRRP8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 01:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRRP8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP8_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP8_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP8_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP8_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP8_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:46:39
% EndTime: 2019-05-06 01:46:48
% DurationCPUTime: 4.90s
% Computational Cost: add. (52717->315), mult. (102995->375), div. (0->0), fcn. (67158->8), ass. (0->122)
t679 = Ifges(6,1) + Ifges(7,1);
t671 = Ifges(6,4) - Ifges(7,5);
t678 = -Ifges(6,5) - Ifges(7,4);
t677 = Ifges(6,2) + Ifges(7,3);
t667 = Ifges(6,6) - Ifges(7,6);
t676 = -Ifges(6,3) - Ifges(7,2);
t637 = sin(qJ(1));
t640 = cos(qJ(1));
t620 = -t640 * g(1) - t637 * g(2);
t648 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t620;
t635 = sin(qJ(4));
t636 = sin(qJ(3));
t638 = cos(qJ(4));
t639 = cos(qJ(3));
t607 = (t635 * t639 + t636 * t638) * qJD(1);
t675 = -pkin(1) - pkin(7);
t674 = cos(qJ(5));
t673 = mrSges(2,1) - mrSges(3,2);
t672 = -mrSges(6,3) - mrSges(7,2);
t669 = Ifges(2,5) - Ifges(3,4);
t668 = (-Ifges(2,6) + Ifges(3,5));
t619 = t637 * g(1) - t640 * g(2);
t641 = qJD(1) ^ 2;
t647 = -t641 * qJ(2) + qJDD(2) - t619;
t599 = t675 * qJDD(1) + t647;
t589 = t636 * g(3) + t639 * t599;
t659 = qJD(1) * qJD(3);
t656 = t636 * t659;
t615 = qJDD(1) * t639 - t656;
t565 = (-t615 - t656) * pkin(8) + (-t636 * t639 * t641 + qJDD(3)) * pkin(3) + t589;
t590 = -g(3) * t639 + t636 * t599;
t614 = -qJDD(1) * t636 - t639 * t659;
t660 = qJD(1) * t639;
t618 = qJD(3) * pkin(3) - pkin(8) * t660;
t631 = t636 ^ 2;
t566 = -pkin(3) * t631 * t641 + pkin(8) * t614 - qJD(3) * t618 + t590;
t547 = t635 * t565 + t638 * t566;
t608 = (-t635 * t636 + t638 * t639) * qJD(1);
t577 = -qJD(4) * t608 + t614 * t638 - t615 * t635;
t587 = mrSges(5,1) * t607 + mrSges(5,2) * t608;
t628 = qJD(3) + qJD(4);
t597 = mrSges(5,1) * t628 - mrSges(5,3) * t608;
t627 = qJDD(3) + qJDD(4);
t571 = -t614 * pkin(3) + t618 * t660 + (-pkin(8) * t631 + t675) * t641 + t648;
t578 = -qJD(4) * t607 + t614 * t635 + t615 * t638;
t542 = (t607 * t628 - t578) * pkin(9) + (t608 * t628 - t577) * pkin(4) + t571;
t588 = pkin(4) * t607 - pkin(9) * t608;
t626 = t628 ^ 2;
t545 = -pkin(4) * t626 + pkin(9) * t627 - t588 * t607 + t547;
t634 = sin(qJ(5));
t540 = t634 * t542 + t674 * t545;
t592 = t674 * t608 + t634 * t628;
t550 = t592 * qJD(5) + t634 * t578 - t674 * t627;
t576 = qJDD(5) - t577;
t603 = qJD(5) + t607;
t581 = mrSges(6,1) * t603 - mrSges(6,3) * t592;
t591 = t634 * t608 - t674 * t628;
t567 = pkin(5) * t591 - qJ(6) * t592;
t601 = t603 ^ 2;
t536 = -pkin(5) * t601 + qJ(6) * t576 + 0.2e1 * qJD(6) * t603 - t567 * t591 + t540;
t582 = -mrSges(7,1) * t603 + mrSges(7,2) * t592;
t657 = m(7) * t536 + t576 * mrSges(7,3) + t603 * t582;
t568 = mrSges(7,1) * t591 - mrSges(7,3) * t592;
t662 = -mrSges(6,1) * t591 - mrSges(6,2) * t592 - t568;
t531 = m(6) * t540 - t576 * mrSges(6,2) + t672 * t550 - t603 * t581 + t662 * t591 + t657;
t539 = t674 * t542 - t634 * t545;
t551 = -t591 * qJD(5) + t674 * t578 + t634 * t627;
t580 = -mrSges(6,2) * t603 - mrSges(6,3) * t591;
t537 = -t576 * pkin(5) - t601 * qJ(6) + t592 * t567 + qJDD(6) - t539;
t579 = -mrSges(7,2) * t591 + mrSges(7,3) * t603;
t651 = -m(7) * t537 + t576 * mrSges(7,1) + t603 * t579;
t533 = m(6) * t539 + t576 * mrSges(6,1) + t672 * t551 + t603 * t580 + t662 * t592 + t651;
t652 = t674 * t531 - t533 * t634;
t523 = m(5) * t547 - mrSges(5,2) * t627 + mrSges(5,3) * t577 - t587 * t607 - t597 * t628 + t652;
t546 = t638 * t565 - t635 * t566;
t596 = -mrSges(5,2) * t628 - mrSges(5,3) * t607;
t544 = -t627 * pkin(4) - t626 * pkin(9) + t608 * t588 - t546;
t538 = -0.2e1 * qJD(6) * t592 + (t591 * t603 - t551) * qJ(6) + (t592 * t603 + t550) * pkin(5) + t544;
t534 = m(7) * t538 + mrSges(7,1) * t550 - t551 * mrSges(7,3) + t579 * t591 - t592 * t582;
t643 = -m(6) * t544 - t550 * mrSges(6,1) - mrSges(6,2) * t551 - t591 * t580 - t581 * t592 - t534;
t528 = m(5) * t546 + mrSges(5,1) * t627 - mrSges(5,3) * t578 - t587 * t608 + t596 * t628 + t643;
t517 = t635 * t523 + t638 * t528;
t613 = (mrSges(4,1) * t636 + mrSges(4,2) * t639) * qJD(1);
t661 = qJD(1) * t636;
t616 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t661;
t515 = m(4) * t589 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t615 + qJD(3) * t616 - t613 * t660 + t517;
t617 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t660;
t653 = t638 * t523 - t528 * t635;
t516 = m(4) * t590 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t614 - qJD(3) * t617 - t613 * t661 + t653;
t512 = t639 * t515 + t636 * t516;
t602 = -qJDD(1) * pkin(1) + t647;
t646 = -m(3) * t602 + (t641 * mrSges(3,3)) - t512;
t509 = m(2) * t619 - (t641 * mrSges(2,2)) + t673 * qJDD(1) + t646;
t600 = t641 * pkin(1) - t648;
t598 = t675 * t641 + t648;
t526 = t634 * t531 + t674 * t533;
t645 = m(5) * t571 - mrSges(5,1) * t577 + t578 * mrSges(5,2) + t596 * t607 + t608 * t597 + t526;
t644 = -m(4) * t598 + mrSges(4,1) * t614 - t615 * mrSges(4,2) - t616 * t661 - t617 * t660 - t645;
t642 = -m(3) * t600 + (t641 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t644;
t520 = m(2) * t620 - (mrSges(2,1) * t641) - qJDD(1) * mrSges(2,2) + t642;
t666 = t640 * t509 + t637 * t520;
t665 = t677 * t591 - t671 * t592 - t667 * t603;
t664 = t667 * t591 + t678 * t592 + t676 * t603;
t663 = -t671 * t591 + t679 * t592 - t678 * t603;
t655 = -t509 * t637 + t640 * t520;
t654 = -t636 * t515 + t639 * t516;
t606 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t639 - Ifges(4,4) * t636) * qJD(1);
t605 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t639 - Ifges(4,2) * t636) * qJD(1);
t604 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t639 - Ifges(4,6) * t636) * qJD(1);
t585 = Ifges(5,1) * t608 - Ifges(5,4) * t607 + Ifges(5,5) * t628;
t584 = Ifges(5,4) * t608 - Ifges(5,2) * t607 + Ifges(5,6) * t628;
t583 = Ifges(5,5) * t608 - Ifges(5,6) * t607 + Ifges(5,3) * t628;
t525 = mrSges(6,2) * t544 + mrSges(7,2) * t537 - mrSges(6,3) * t539 - mrSges(7,3) * t538 - qJ(6) * t534 - t671 * t550 + t679 * t551 - t576 * t678 + t664 * t591 + t665 * t603;
t524 = -mrSges(6,1) * t544 - mrSges(7,1) * t538 + mrSges(7,2) * t536 + mrSges(6,3) * t540 - pkin(5) * t534 - t677 * t550 + t671 * t551 + t667 * t576 + t664 * t592 + t663 * t603;
t513 = Ifges(5,4) * t578 + Ifges(5,2) * t577 + Ifges(5,6) * t627 - t608 * t583 + t628 * t585 - mrSges(5,1) * t571 + mrSges(5,3) * t547 - mrSges(6,1) * t539 + mrSges(6,2) * t540 + mrSges(7,1) * t537 - mrSges(7,3) * t536 - pkin(5) * t651 - qJ(6) * t657 - pkin(4) * t526 + (pkin(5) * t568 + t665) * t592 + (qJ(6) * t568 - t663) * t591 + t676 * t576 + (mrSges(7,2) * pkin(5) + t678) * t551 + (mrSges(7,2) * qJ(6) + t667) * t550;
t511 = -m(3) * g(3) + t654;
t510 = mrSges(5,2) * t571 - mrSges(5,3) * t546 + Ifges(5,1) * t578 + Ifges(5,4) * t577 + Ifges(5,5) * t627 - pkin(9) * t526 - t634 * t524 + t674 * t525 - t607 * t583 - t628 * t584;
t507 = mrSges(4,2) * t598 - mrSges(4,3) * t589 + Ifges(4,1) * t615 + Ifges(4,4) * t614 + Ifges(4,5) * qJDD(3) - pkin(8) * t517 - qJD(3) * t605 + t510 * t638 - t513 * t635 - t604 * t661;
t506 = -mrSges(4,1) * t598 + mrSges(4,3) * t590 + Ifges(4,4) * t615 + Ifges(4,2) * t614 + Ifges(4,6) * qJDD(3) - pkin(3) * t645 + pkin(8) * t653 + qJD(3) * t606 + t635 * t510 + t638 * t513 - t604 * t660;
t505 = t674 * t524 + pkin(2) * t512 + t669 * qJDD(1) + (t668 * t641) - qJ(2) * t511 + pkin(9) * t652 + pkin(4) * t643 + Ifges(4,3) * qJDD(3) + (t605 * t639 + t606 * t636) * qJD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + pkin(3) * t517 + t634 * t525 - mrSges(2,3) * t619 + Ifges(5,3) * t627 + Ifges(4,6) * t614 + Ifges(4,5) * t615 + mrSges(3,1) * t602 + t607 * t585 + t608 * t584 + mrSges(4,1) * t589 - mrSges(4,2) * t590 + Ifges(5,6) * t577 + Ifges(5,5) * t578 + mrSges(5,1) * t546 - mrSges(5,2) * t547;
t504 = -mrSges(3,1) * t600 + mrSges(2,3) * t620 - pkin(1) * t511 - pkin(2) * t644 - pkin(7) * t654 + t673 * g(3) - t668 * qJDD(1) - t639 * t506 - t636 * t507 + t669 * t641;
t1 = [-m(1) * g(1) + t655; -m(1) * g(2) + t666; (-m(1) - m(2) - m(3)) * g(3) + t654; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t666 - t637 * t504 + t640 * t505; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t655 + t640 * t504 + t637 * t505; pkin(1) * t646 + qJ(2) * t642 - mrSges(2,2) * t620 + t639 * t507 - t636 * t506 - pkin(7) * t512 + mrSges(2,1) * t619 + mrSges(3,2) * t602 - mrSges(3,3) * t600 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-mrSges(3,2) * pkin(1) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;
