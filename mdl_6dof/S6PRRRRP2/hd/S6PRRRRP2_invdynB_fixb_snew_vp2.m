% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRRRP2
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
% Datum: 2019-05-05 09:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRRRP2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:35:16
% EndTime: 2019-05-05 09:35:28
% DurationCPUTime: 10.18s
% Computational Cost: add. (164949->318), mult. (322222->400), div. (0->0), fcn. (229338->12), ass. (0->133)
t710 = Ifges(6,1) + Ifges(7,1);
t704 = Ifges(6,4) - Ifges(7,5);
t709 = -Ifges(6,5) - Ifges(7,4);
t708 = Ifges(6,2) + Ifges(7,3);
t702 = Ifges(6,6) - Ifges(7,6);
t707 = -Ifges(6,3) - Ifges(7,2);
t670 = sin(qJ(4));
t671 = sin(qJ(3));
t673 = cos(qJ(4));
t674 = cos(qJ(3));
t642 = (t670 * t671 - t673 * t674) * qJD(2);
t665 = sin(pkin(11));
t667 = cos(pkin(11));
t652 = g(1) * t665 - g(2) * t667;
t653 = -g(1) * t667 - g(2) * t665;
t664 = -g(3) + qJDD(1);
t666 = sin(pkin(6));
t668 = cos(pkin(6));
t672 = sin(qJ(2));
t675 = cos(qJ(2));
t624 = -t672 * t653 + (t652 * t668 + t664 * t666) * t675;
t706 = cos(qJ(5));
t705 = -mrSges(6,3) - mrSges(7,2);
t676 = qJD(2) ^ 2;
t680 = -qJDD(2) * pkin(2) - t624;
t619 = -pkin(8) * t676 + t680;
t691 = qJD(2) * qJD(3);
t689 = t674 * t691;
t650 = qJDD(2) * t671 + t689;
t651 = qJDD(2) * t674 - t671 * t691;
t693 = qJD(2) * t671;
t654 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t693;
t692 = qJD(2) * t674;
t655 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t692;
t699 = t668 * t672;
t700 = t666 * t672;
t625 = t652 * t699 + t675 * t653 + t664 * t700;
t620 = -pkin(2) * t676 + qJDD(2) * pkin(8) + t625;
t635 = -t652 * t666 + t664 * t668;
t600 = -t620 * t671 + t674 * t635;
t585 = (-t650 + t689) * pkin(9) + (t671 * t674 * t676 + qJDD(3)) * pkin(3) + t600;
t601 = t674 * t620 + t671 * t635;
t657 = qJD(3) * pkin(3) - pkin(9) * t693;
t663 = t674 ^ 2;
t586 = -pkin(3) * t663 * t676 + pkin(9) * t651 - qJD(3) * t657 + t601;
t581 = t670 * t585 + t673 * t586;
t643 = (t670 * t674 + t671 * t673) * qJD(2);
t628 = pkin(4) * t642 - pkin(10) * t643;
t662 = qJD(3) + qJD(4);
t660 = t662 ^ 2;
t661 = qJDD(3) + qJDD(4);
t577 = -pkin(4) * t660 + pkin(10) * t661 - t628 * t642 + t581;
t591 = -pkin(3) * t651 + t657 * t693 + (-pkin(9) * t663 - pkin(8)) * t676 + t680;
t613 = -qJD(4) * t643 - t650 * t670 + t651 * t673;
t614 = -qJD(4) * t642 + t650 * t673 + t651 * t670;
t579 = (t642 * t662 - t614) * pkin(10) + (t643 * t662 - t613) * pkin(4) + t591;
t669 = sin(qJ(5));
t574 = t706 * t577 + t669 * t579;
t630 = t643 * t706 + t669 * t662;
t589 = qJD(5) * t630 + t614 * t669 - t661 * t706;
t611 = qJDD(5) - t613;
t637 = qJD(5) + t642;
t617 = mrSges(6,1) * t637 - mrSges(6,3) * t630;
t629 = t643 * t669 - t662 * t706;
t604 = pkin(5) * t629 - qJ(6) * t630;
t636 = t637 ^ 2;
t570 = -pkin(5) * t636 + qJ(6) * t611 + 0.2e1 * qJD(6) * t637 - t604 * t629 + t574;
t618 = -mrSges(7,1) * t637 + mrSges(7,2) * t630;
t690 = m(7) * t570 + t611 * mrSges(7,3) + t637 * t618;
t605 = mrSges(7,1) * t629 - mrSges(7,3) * t630;
t694 = -mrSges(6,1) * t629 - mrSges(6,2) * t630 - t605;
t565 = m(6) * t574 - mrSges(6,2) * t611 + t589 * t705 - t617 * t637 + t629 * t694 + t690;
t573 = -t669 * t577 + t579 * t706;
t590 = -t629 * qJD(5) + t614 * t706 + t669 * t661;
t616 = -mrSges(6,2) * t637 - mrSges(6,3) * t629;
t571 = -t611 * pkin(5) - t636 * qJ(6) + t630 * t604 + qJDD(6) - t573;
t615 = -mrSges(7,2) * t629 + mrSges(7,3) * t637;
t684 = -m(7) * t571 + t611 * mrSges(7,1) + t637 * t615;
t567 = m(6) * t573 + mrSges(6,1) * t611 + t590 * t705 + t616 * t637 + t630 * t694 + t684;
t560 = t669 * t565 + t706 * t567;
t633 = -mrSges(5,2) * t662 - mrSges(5,3) * t642;
t634 = mrSges(5,1) * t662 - mrSges(5,3) * t643;
t679 = m(5) * t591 - t613 * mrSges(5,1) + mrSges(5,2) * t614 + t642 * t633 + t634 * t643 + t560;
t677 = -m(4) * t619 + t651 * mrSges(4,1) - mrSges(4,2) * t650 - t654 * t693 + t655 * t692 - t679;
t554 = m(3) * t624 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t676 + t677;
t701 = t554 * t675;
t627 = mrSges(5,1) * t642 + mrSges(5,2) * t643;
t685 = t706 * t565 - t567 * t669;
t557 = m(5) * t581 - mrSges(5,2) * t661 + mrSges(5,3) * t613 - t627 * t642 - t634 * t662 + t685;
t580 = t585 * t673 - t670 * t586;
t576 = -pkin(4) * t661 - pkin(10) * t660 + t643 * t628 - t580;
t572 = -0.2e1 * qJD(6) * t630 + (t629 * t637 - t590) * qJ(6) + (t630 * t637 + t589) * pkin(5) + t576;
t568 = m(7) * t572 + mrSges(7,1) * t589 - t590 * mrSges(7,3) + t615 * t629 - t630 * t618;
t678 = -m(6) * t576 - t589 * mrSges(6,1) - mrSges(6,2) * t590 - t629 * t616 - t617 * t630 - t568;
t562 = m(5) * t580 + mrSges(5,1) * t661 - mrSges(5,3) * t614 - t627 * t643 + t633 * t662 + t678;
t551 = t670 * t557 + t673 * t562;
t649 = (-mrSges(4,1) * t674 + mrSges(4,2) * t671) * qJD(2);
t549 = m(4) * t600 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t650 + qJD(3) * t655 - t649 * t693 + t551;
t686 = t673 * t557 - t562 * t670;
t550 = m(4) * t601 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t651 - qJD(3) * t654 + t649 * t692 + t686;
t687 = -t549 * t671 + t674 * t550;
t540 = m(3) * t625 - mrSges(3,1) * t676 - qJDD(2) * mrSges(3,2) + t687;
t543 = t674 * t549 + t671 * t550;
t542 = m(3) * t635 + t543;
t531 = t540 * t699 - t542 * t666 + t668 * t701;
t529 = m(2) * t652 + t531;
t536 = t675 * t540 - t554 * t672;
t535 = m(2) * t653 + t536;
t698 = t667 * t529 + t665 * t535;
t697 = t629 * t708 - t630 * t704 - t637 * t702;
t696 = t629 * t702 + t630 * t709 + t637 * t707;
t695 = -t704 * t629 + t630 * t710 - t709 * t637;
t530 = t540 * t700 + t668 * t542 + t666 * t701;
t688 = -t529 * t665 + t667 * t535;
t558 = -mrSges(6,1) * t576 - mrSges(7,1) * t572 + mrSges(7,2) * t570 + mrSges(6,3) * t574 - pkin(5) * t568 - t589 * t708 + t704 * t590 + t702 * t611 + t696 * t630 + t695 * t637;
t559 = mrSges(6,2) * t576 + mrSges(7,2) * t571 - mrSges(6,3) * t573 - mrSges(7,3) * t572 - qJ(6) * t568 - t704 * t589 + t590 * t710 - t611 * t709 + t696 * t629 + t697 * t637;
t621 = Ifges(5,5) * t643 - Ifges(5,6) * t642 + Ifges(5,3) * t662;
t622 = Ifges(5,4) * t643 - Ifges(5,2) * t642 + Ifges(5,6) * t662;
t544 = mrSges(5,2) * t591 - mrSges(5,3) * t580 + Ifges(5,1) * t614 + Ifges(5,4) * t613 + Ifges(5,5) * t661 - pkin(10) * t560 - t669 * t558 + t559 * t706 - t642 * t621 - t662 * t622;
t623 = Ifges(5,1) * t643 - Ifges(5,4) * t642 + Ifges(5,5) * t662;
t545 = Ifges(5,4) * t614 + Ifges(5,2) * t613 + Ifges(5,6) * t661 - t643 * t621 + t662 * t623 - mrSges(5,1) * t591 + mrSges(5,3) * t581 - mrSges(6,1) * t573 + mrSges(6,2) * t574 + mrSges(7,1) * t571 - mrSges(7,3) * t570 - pkin(5) * t684 - qJ(6) * t690 - pkin(4) * t560 + (pkin(5) * t605 + t697) * t630 + (qJ(6) * t605 - t695) * t629 + t707 * t611 + (pkin(5) * mrSges(7,2) + t709) * t590 + (mrSges(7,2) * qJ(6) + t702) * t589;
t639 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t671 + Ifges(4,6) * t674) * qJD(2);
t641 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t671 + Ifges(4,4) * t674) * qJD(2);
t527 = -mrSges(4,1) * t619 + mrSges(4,3) * t601 + Ifges(4,4) * t650 + Ifges(4,2) * t651 + Ifges(4,6) * qJDD(3) - pkin(3) * t679 + pkin(9) * t686 + qJD(3) * t641 + t670 * t544 + t673 * t545 - t639 * t693;
t640 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t671 + Ifges(4,2) * t674) * qJD(2);
t532 = mrSges(4,2) * t619 - mrSges(4,3) * t600 + Ifges(4,1) * t650 + Ifges(4,4) * t651 + Ifges(4,5) * qJDD(3) - pkin(9) * t551 - qJD(3) * t640 + t544 * t673 - t545 * t670 + t639 * t692;
t525 = mrSges(3,2) * t635 - mrSges(3,3) * t624 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t676 - pkin(8) * t543 - t527 * t671 + t532 * t674;
t526 = Ifges(3,6) * qJDD(2) - pkin(2) * t543 + mrSges(3,3) * t625 - mrSges(3,1) * t635 - pkin(3) * t551 + mrSges(4,2) * t601 - Ifges(4,5) * t650 - Ifges(4,6) * t651 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t600 - pkin(10) * t685 - Ifges(5,5) * t614 - Ifges(5,6) * t613 - Ifges(5,3) * t661 - mrSges(5,1) * t580 + mrSges(5,2) * t581 - t669 * t559 - t706 * t558 - pkin(4) * t678 - t643 * t622 - t642 * t623 + t676 * Ifges(3,5) + (-t671 * t640 + t674 * t641) * qJD(2);
t681 = pkin(7) * t536 + t525 * t672 + t526 * t675;
t524 = mrSges(3,1) * t624 - mrSges(3,2) * t625 + Ifges(3,3) * qJDD(2) + pkin(2) * t677 + pkin(8) * t687 + t674 * t527 + t671 * t532;
t523 = mrSges(2,2) * t664 - mrSges(2,3) * t652 + t525 * t675 - t526 * t672 + (-t530 * t666 - t531 * t668) * pkin(7);
t522 = -mrSges(2,1) * t664 + mrSges(2,3) * t653 - pkin(1) * t530 - t524 * t666 + t668 * t681;
t1 = [-m(1) * g(1) + t688; -m(1) * g(2) + t698; -m(1) * g(3) + m(2) * t664 + t530; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t698 - t665 * t522 + t667 * t523; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t688 + t667 * t522 + t665 * t523; -mrSges(1,1) * g(2) + mrSges(2,1) * t652 + mrSges(1,2) * g(1) - mrSges(2,2) * t653 + pkin(1) * t531 + t524 * t668 + t666 * t681;];
tauB  = t1;
