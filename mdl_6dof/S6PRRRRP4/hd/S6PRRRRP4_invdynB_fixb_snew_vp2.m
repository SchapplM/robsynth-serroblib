% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRRRP4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-05-05 09:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRRRP4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:51:32
% EndTime: 2019-05-05 09:51:48
% DurationCPUTime: 10.54s
% Computational Cost: add. (173493->317), mult. (334039->395), div. (0->0), fcn. (236614->12), ass. (0->133)
t698 = Ifges(6,1) + Ifges(7,1);
t693 = Ifges(6,4) - Ifges(7,5);
t692 = Ifges(7,4) + Ifges(6,5);
t697 = Ifges(6,2) + Ifges(7,3);
t691 = Ifges(6,6) - Ifges(7,6);
t696 = -Ifges(6,3) - Ifges(7,2);
t655 = sin(pkin(11));
t657 = cos(pkin(11));
t645 = g(1) * t655 - g(2) * t657;
t646 = -g(1) * t657 - g(2) * t655;
t654 = -g(3) + qJDD(1);
t656 = sin(pkin(6));
t658 = cos(pkin(6));
t662 = sin(qJ(2));
t665 = cos(qJ(2));
t605 = -t662 * t646 + (t645 * t658 + t654 * t656) * t665;
t695 = cos(qJ(5));
t694 = -mrSges(6,3) - mrSges(7,2);
t667 = qJD(2) ^ 2;
t688 = t658 * t662;
t689 = t656 * t662;
t606 = t645 * t688 + t665 * t646 + t654 * t689;
t600 = -pkin(2) * t667 + qJDD(2) * pkin(8) + t606;
t623 = -t645 * t656 + t654 * t658;
t661 = sin(qJ(3));
t664 = cos(qJ(3));
t591 = t664 * t600 + t661 * t623;
t642 = (-pkin(3) * t664 - pkin(9) * t661) * qJD(2);
t666 = qJD(3) ^ 2;
t681 = qJD(2) * t664;
t577 = -pkin(3) * t666 + qJDD(3) * pkin(9) + t642 * t681 + t591;
t599 = -qJDD(2) * pkin(2) - t667 * pkin(8) - t605;
t680 = qJD(2) * qJD(3);
t678 = t664 * t680;
t643 = qJDD(2) * t661 + t678;
t653 = t661 * t680;
t644 = qJDD(2) * t664 - t653;
t580 = (-t643 - t678) * pkin(9) + (-t644 + t653) * pkin(3) + t599;
t660 = sin(qJ(4));
t663 = cos(qJ(4));
t567 = -t660 * t577 + t663 * t580;
t682 = qJD(2) * t661;
t639 = qJD(3) * t663 - t660 * t682;
t615 = qJD(4) * t639 + qJDD(3) * t660 + t643 * t663;
t636 = qJDD(4) - t644;
t640 = qJD(3) * t660 + t663 * t682;
t652 = qJD(4) - t681;
t564 = (t639 * t652 - t615) * pkin(10) + (t639 * t640 + t636) * pkin(4) + t567;
t568 = t663 * t577 + t660 * t580;
t614 = -qJD(4) * t640 + qJDD(3) * t663 - t643 * t660;
t622 = pkin(4) * t652 - pkin(10) * t640;
t635 = t639 ^ 2;
t566 = -pkin(4) * t635 + pkin(10) * t614 - t622 * t652 + t568;
t659 = sin(qJ(5));
t560 = t659 * t564 + t566 * t695;
t617 = t659 * t639 + t640 * t695;
t573 = qJD(5) * t617 - t614 * t695 + t615 * t659;
t651 = qJD(5) + t652;
t603 = mrSges(6,1) * t651 - mrSges(6,3) * t617;
t616 = -t639 * t695 + t640 * t659;
t632 = qJDD(5) + t636;
t592 = pkin(5) * t616 - qJ(6) * t617;
t650 = t651 ^ 2;
t557 = -pkin(5) * t650 + qJ(6) * t632 + 0.2e1 * qJD(6) * t651 - t592 * t616 + t560;
t604 = -mrSges(7,1) * t651 + mrSges(7,2) * t617;
t679 = m(7) * t557 + t632 * mrSges(7,3) + t651 * t604;
t593 = mrSges(7,1) * t616 - mrSges(7,3) * t617;
t683 = -mrSges(6,1) * t616 - mrSges(6,2) * t617 - t593;
t552 = m(6) * t560 - t632 * mrSges(6,2) + t573 * t694 - t651 * t603 + t616 * t683 + t679;
t559 = t564 * t695 - t659 * t566;
t574 = -t616 * qJD(5) + t659 * t614 + t615 * t695;
t602 = -mrSges(6,2) * t651 - mrSges(6,3) * t616;
t558 = -t632 * pkin(5) - t650 * qJ(6) + t617 * t592 + qJDD(6) - t559;
t601 = -mrSges(7,2) * t616 + mrSges(7,3) * t651;
t673 = -m(7) * t558 + t632 * mrSges(7,1) + t651 * t601;
t554 = m(6) * t559 + t632 * mrSges(6,1) + t574 * t694 + t651 * t602 + t617 * t683 + t673;
t547 = t659 * t552 + t554 * t695;
t618 = -mrSges(5,1) * t639 + mrSges(5,2) * t640;
t620 = -mrSges(5,2) * t652 + mrSges(5,3) * t639;
t543 = m(5) * t567 + mrSges(5,1) * t636 - mrSges(5,3) * t615 - t618 * t640 + t620 * t652 + t547;
t621 = mrSges(5,1) * t652 - mrSges(5,3) * t640;
t674 = t552 * t695 - t554 * t659;
t544 = m(5) * t568 - mrSges(5,2) * t636 + mrSges(5,3) * t614 + t618 * t639 - t621 * t652 + t674;
t541 = t543 * t663 + t544 * t660;
t647 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t682;
t648 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t681;
t669 = -m(4) * t599 + t644 * mrSges(4,1) - mrSges(4,2) * t643 - t647 * t682 + t648 * t681 - t541;
t537 = m(3) * t605 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t667 + t669;
t690 = t537 * t665;
t641 = (-mrSges(4,1) * t664 + mrSges(4,2) * t661) * qJD(2);
t675 = -t543 * t660 + t663 * t544;
t540 = m(4) * t591 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t644 - qJD(3) * t647 + t641 * t681 + t675;
t590 = -t661 * t600 + t664 * t623;
t576 = -qJDD(3) * pkin(3) - t666 * pkin(9) + t642 * t682 - t590;
t569 = -t614 * pkin(4) - t635 * pkin(10) + t640 * t622 + t576;
t562 = -0.2e1 * qJD(6) * t617 + (t616 * t651 - t574) * qJ(6) + (t617 * t651 + t573) * pkin(5) + t569;
t555 = m(7) * t562 + t573 * mrSges(7,1) - t574 * mrSges(7,3) + t616 * t601 - t617 * t604;
t670 = m(6) * t569 + t573 * mrSges(6,1) + t574 * mrSges(6,2) + t616 * t602 + t617 * t603 + t555;
t668 = -m(5) * t576 + t614 * mrSges(5,1) - t615 * mrSges(5,2) + t639 * t620 - t640 * t621 - t670;
t549 = m(4) * t590 + qJDD(3) * mrSges(4,1) - t643 * mrSges(4,3) + qJD(3) * t648 - t641 * t682 + t668;
t676 = t664 * t540 - t549 * t661;
t531 = m(3) * t606 - mrSges(3,1) * t667 - qJDD(2) * mrSges(3,2) + t676;
t534 = t661 * t540 + t664 * t549;
t533 = m(3) * t623 + t534;
t520 = t531 * t688 - t533 * t656 + t658 * t690;
t518 = m(2) * t645 + t520;
t524 = t665 * t531 - t537 * t662;
t523 = m(2) * t646 + t524;
t687 = t657 * t518 + t655 * t523;
t686 = t697 * t616 - t693 * t617 - t691 * t651;
t685 = t691 * t616 - t692 * t617 + t696 * t651;
t684 = -t693 * t616 + t698 * t617 + t692 * t651;
t519 = t531 * t689 + t658 * t533 + t656 * t690;
t677 = -t518 * t655 + t657 * t523;
t545 = -mrSges(6,1) * t569 - mrSges(7,1) * t562 + mrSges(7,2) * t557 + mrSges(6,3) * t560 - pkin(5) * t555 - t697 * t573 + t693 * t574 + t685 * t617 + t691 * t632 + t684 * t651;
t546 = mrSges(6,2) * t569 + mrSges(7,2) * t558 - mrSges(6,3) * t559 - mrSges(7,3) * t562 - qJ(6) * t555 - t693 * t573 + t698 * t574 + t685 * t616 + t692 * t632 + t686 * t651;
t608 = Ifges(5,5) * t640 + Ifges(5,6) * t639 + Ifges(5,3) * t652;
t610 = Ifges(5,1) * t640 + Ifges(5,4) * t639 + Ifges(5,5) * t652;
t526 = -mrSges(5,1) * t576 + mrSges(5,3) * t568 + Ifges(5,4) * t615 + Ifges(5,2) * t614 + Ifges(5,6) * t636 - pkin(4) * t670 + pkin(10) * t674 + t545 * t695 + t659 * t546 - t640 * t608 + t652 * t610;
t609 = Ifges(5,4) * t640 + Ifges(5,2) * t639 + Ifges(5,6) * t652;
t527 = mrSges(5,2) * t576 - mrSges(5,3) * t567 + Ifges(5,1) * t615 + Ifges(5,4) * t614 + Ifges(5,5) * t636 - pkin(10) * t547 - t659 * t545 + t546 * t695 + t639 * t608 - t652 * t609;
t629 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t661 + Ifges(4,6) * t664) * qJD(2);
t630 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t661 + Ifges(4,2) * t664) * qJD(2);
t516 = mrSges(4,2) * t599 - mrSges(4,3) * t590 + Ifges(4,1) * t643 + Ifges(4,4) * t644 + Ifges(4,5) * qJDD(3) - pkin(9) * t541 - qJD(3) * t630 - t526 * t660 + t527 * t663 + t629 * t681;
t631 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t661 + Ifges(4,4) * t664) * qJD(2);
t525 = Ifges(4,6) * qJDD(3) + t696 * t632 + Ifges(4,2) * t644 - Ifges(5,3) * t636 + t639 * t610 - t640 * t609 + Ifges(4,4) * t643 + qJD(3) * t631 - Ifges(5,6) * t614 - Ifges(5,5) * t615 + mrSges(4,3) * t591 - mrSges(4,1) * t599 - mrSges(5,1) * t567 + mrSges(5,2) * t568 + mrSges(6,2) * t560 - mrSges(7,3) * t557 + mrSges(7,1) * t558 - mrSges(6,1) * t559 - pkin(4) * t547 - pkin(3) * t541 - pkin(5) * t673 - qJ(6) * t679 - t629 * t682 + (qJ(6) * t593 - t684) * t616 + (pkin(5) * t593 + t686) * t617 + (mrSges(7,2) * qJ(6) + t691) * t573 + (mrSges(7,2) * pkin(5) - t692) * t574;
t514 = mrSges(3,2) * t623 - mrSges(3,3) * t605 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t667 - pkin(8) * t534 + t516 * t664 - t525 * t661;
t515 = Ifges(3,6) * qJDD(2) + t667 * Ifges(3,5) - mrSges(3,1) * t623 + mrSges(3,3) * t606 - Ifges(4,5) * t643 - Ifges(4,6) * t644 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t590 + mrSges(4,2) * t591 - t660 * t527 - t663 * t526 - pkin(3) * t668 - pkin(9) * t675 - pkin(2) * t534 + (-t630 * t661 + t631 * t664) * qJD(2);
t671 = pkin(7) * t524 + t514 * t662 + t515 * t665;
t513 = mrSges(3,1) * t605 - mrSges(3,2) * t606 + Ifges(3,3) * qJDD(2) + pkin(2) * t669 + pkin(8) * t676 + t661 * t516 + t664 * t525;
t512 = mrSges(2,2) * t654 - mrSges(2,3) * t645 + t665 * t514 - t662 * t515 + (-t519 * t656 - t520 * t658) * pkin(7);
t511 = -mrSges(2,1) * t654 + mrSges(2,3) * t646 - pkin(1) * t519 - t656 * t513 + t658 * t671;
t1 = [-m(1) * g(1) + t677; -m(1) * g(2) + t687; -m(1) * g(3) + m(2) * t654 + t519; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t687 - t655 * t511 + t657 * t512; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t677 + t657 * t511 + t655 * t512; -mrSges(1,1) * g(2) + mrSges(2,1) * t645 + mrSges(1,2) * g(1) - mrSges(2,2) * t646 + pkin(1) * t520 + t658 * t513 + t656 * t671;];
tauB  = t1;
