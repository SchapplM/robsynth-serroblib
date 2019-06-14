% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-05-05 22:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRPR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR6_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR6_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR6_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR6_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR6_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:49:23
% EndTime: 2019-05-05 22:49:49
% DurationCPUTime: 24.98s
% Computational Cost: add. (378312->365), mult. (912610->459), div. (0->0), fcn. (703334->12), ass. (0->151)
t733 = qJD(1) ^ 2;
t723 = cos(pkin(10));
t761 = pkin(2) * t723;
t721 = sin(pkin(10));
t760 = mrSges(3,2) * t721;
t719 = t723 ^ 2;
t759 = t719 * t733;
t727 = sin(qJ(1));
t731 = cos(qJ(1));
t709 = -g(1) * t731 - g(2) * t727;
t705 = -pkin(1) * t733 + qJDD(1) * qJ(2) + t709;
t754 = qJD(1) * qJD(2);
t751 = -t723 * g(3) - 0.2e1 * t721 * t754;
t677 = (-pkin(7) * qJDD(1) + t733 * t761 - t705) * t721 + t751;
t692 = -g(3) * t721 + (t705 + 0.2e1 * t754) * t723;
t752 = qJDD(1) * t723;
t678 = -pkin(2) * t759 + pkin(7) * t752 + t692;
t726 = sin(qJ(3));
t730 = cos(qJ(3));
t652 = t726 * t677 + t730 * t678;
t756 = qJD(1) * t723;
t757 = qJD(1) * t721;
t703 = -t726 * t757 + t730 * t756;
t739 = t721 * t730 + t723 * t726;
t704 = t739 * qJD(1);
t684 = -mrSges(4,1) * t703 + mrSges(4,2) * t704;
t700 = t704 * qJD(3);
t753 = qJDD(1) * t721;
t689 = -t726 * t753 + t730 * t752 - t700;
t697 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t704;
t687 = -pkin(3) * t703 - pkin(8) * t704;
t732 = qJD(3) ^ 2;
t635 = -pkin(3) * t732 + qJDD(3) * pkin(8) + t687 * t703 + t652;
t718 = t721 ^ 2;
t708 = t727 * g(1) - t731 * g(2);
t744 = qJDD(2) - t708;
t688 = (-pkin(1) - t761) * qJDD(1) + (-qJ(2) + (-t718 - t719) * pkin(7)) * t733 + t744;
t755 = t703 * qJD(3);
t690 = qJDD(1) * t739 + t755;
t643 = (-t690 - t755) * pkin(8) + (-t689 + t700) * pkin(3) + t688;
t725 = sin(qJ(4));
t729 = cos(qJ(4));
t625 = -t725 * t635 + t729 * t643;
t694 = qJD(3) * t729 - t704 * t725;
t663 = qJD(4) * t694 + qJDD(3) * t725 + t690 * t729;
t686 = qJDD(4) - t689;
t695 = qJD(3) * t725 + t704 * t729;
t701 = qJD(4) - t703;
t617 = (t694 * t701 - t663) * qJ(5) + (t694 * t695 + t686) * pkin(4) + t625;
t626 = t729 * t635 + t725 * t643;
t662 = -qJD(4) * t695 + qJDD(3) * t729 - t690 * t725;
t673 = pkin(4) * t701 - qJ(5) * t695;
t693 = t694 ^ 2;
t619 = -pkin(4) * t693 + qJ(5) * t662 - t673 * t701 + t626;
t720 = sin(pkin(11));
t722 = cos(pkin(11));
t668 = t694 * t720 + t695 * t722;
t611 = -0.2e1 * qJD(5) * t668 + t722 * t617 - t720 * t619;
t638 = t662 * t720 + t663 * t722;
t667 = t694 * t722 - t695 * t720;
t609 = (t667 * t701 - t638) * pkin(9) + (t667 * t668 + t686) * pkin(5) + t611;
t612 = 0.2e1 * qJD(5) * t667 + t720 * t617 + t722 * t619;
t637 = t662 * t722 - t663 * t720;
t655 = pkin(5) * t701 - pkin(9) * t668;
t666 = t667 ^ 2;
t610 = -pkin(5) * t666 + pkin(9) * t637 - t655 * t701 + t612;
t724 = sin(qJ(6));
t728 = cos(qJ(6));
t607 = t609 * t728 - t610 * t724;
t648 = t667 * t728 - t668 * t724;
t623 = qJD(6) * t648 + t637 * t724 + t638 * t728;
t649 = t667 * t724 + t668 * t728;
t632 = -mrSges(7,1) * t648 + mrSges(7,2) * t649;
t699 = qJD(6) + t701;
t641 = -mrSges(7,2) * t699 + mrSges(7,3) * t648;
t683 = qJDD(6) + t686;
t603 = m(7) * t607 + mrSges(7,1) * t683 - mrSges(7,3) * t623 - t632 * t649 + t641 * t699;
t608 = t609 * t724 + t610 * t728;
t622 = -qJD(6) * t649 + t637 * t728 - t638 * t724;
t642 = mrSges(7,1) * t699 - mrSges(7,3) * t649;
t604 = m(7) * t608 - mrSges(7,2) * t683 + mrSges(7,3) * t622 + t632 * t648 - t642 * t699;
t597 = t728 * t603 + t724 * t604;
t650 = -mrSges(6,1) * t667 + mrSges(6,2) * t668;
t653 = -mrSges(6,2) * t701 + mrSges(6,3) * t667;
t595 = m(6) * t611 + mrSges(6,1) * t686 - mrSges(6,3) * t638 - t650 * t668 + t653 * t701 + t597;
t654 = mrSges(6,1) * t701 - mrSges(6,3) * t668;
t745 = -t603 * t724 + t728 * t604;
t596 = m(6) * t612 - mrSges(6,2) * t686 + mrSges(6,3) * t637 + t650 * t667 - t654 * t701 + t745;
t591 = t722 * t595 + t720 * t596;
t669 = -mrSges(5,1) * t694 + mrSges(5,2) * t695;
t672 = -mrSges(5,2) * t701 + mrSges(5,3) * t694;
t589 = m(5) * t625 + mrSges(5,1) * t686 - mrSges(5,3) * t663 - t669 * t695 + t672 * t701 + t591;
t674 = mrSges(5,1) * t701 - mrSges(5,3) * t695;
t746 = -t595 * t720 + t722 * t596;
t590 = m(5) * t626 - mrSges(5,2) * t686 + mrSges(5,3) * t662 + t669 * t694 - t674 * t701 + t746;
t747 = -t589 * t725 + t729 * t590;
t582 = m(4) * t652 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t689 - qJD(3) * t697 + t684 * t703 + t747;
t651 = t677 * t730 - t726 * t678;
t696 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t703;
t634 = -qJDD(3) * pkin(3) - pkin(8) * t732 + t704 * t687 - t651;
t624 = -pkin(4) * t662 - qJ(5) * t693 + t695 * t673 + qJDD(5) + t634;
t614 = -pkin(5) * t637 - pkin(9) * t666 + t655 * t668 + t624;
t740 = m(7) * t614 - t622 * mrSges(7,1) + t623 * mrSges(7,2) - t648 * t641 + t649 * t642;
t736 = m(6) * t624 - t637 * mrSges(6,1) + mrSges(6,2) * t638 - t667 * t653 + t654 * t668 + t740;
t734 = -m(5) * t634 + t662 * mrSges(5,1) - mrSges(5,2) * t663 + t694 * t672 - t674 * t695 - t736;
t606 = m(4) * t651 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t690 + qJD(3) * t696 - t684 * t704 + t734;
t577 = t726 * t582 + t730 * t606;
t691 = -t721 * t705 + t751;
t738 = mrSges(3,3) * qJDD(1) + t733 * (-mrSges(3,1) * t723 + t760);
t575 = m(3) * t691 - t721 * t738 + t577;
t748 = t730 * t582 - t726 * t606;
t576 = m(3) * t692 + t723 * t738 + t748;
t749 = -t575 * t721 + t723 * t576;
t567 = m(2) * t709 - mrSges(2,1) * t733 - qJDD(1) * mrSges(2,2) + t749;
t702 = -qJDD(1) * pkin(1) - t733 * qJ(2) + t744;
t583 = t729 * t589 + t725 * t590;
t737 = m(4) * t688 - t689 * mrSges(4,1) + t690 * mrSges(4,2) - t703 * t696 + t704 * t697 + t583;
t735 = -m(3) * t702 + mrSges(3,1) * t752 - t737 + (t718 * t733 + t759) * mrSges(3,3);
t579 = (mrSges(2,1) - t760) * qJDD(1) + t735 - t733 * mrSges(2,2) + m(2) * t708;
t758 = t727 * t567 + t731 * t579;
t568 = t723 * t575 + t721 * t576;
t750 = t731 * t567 - t579 * t727;
t743 = Ifges(3,1) * t721 + Ifges(3,4) * t723;
t742 = Ifges(3,4) * t721 + Ifges(3,2) * t723;
t741 = Ifges(3,5) * t721 + Ifges(3,6) * t723;
t707 = t741 * qJD(1);
t681 = Ifges(4,1) * t704 + Ifges(4,4) * t703 + Ifges(4,5) * qJD(3);
t680 = Ifges(4,4) * t704 + Ifges(4,2) * t703 + Ifges(4,6) * qJD(3);
t679 = Ifges(4,5) * t704 + Ifges(4,6) * t703 + Ifges(4,3) * qJD(3);
t658 = Ifges(5,1) * t695 + Ifges(5,4) * t694 + Ifges(5,5) * t701;
t657 = Ifges(5,4) * t695 + Ifges(5,2) * t694 + Ifges(5,6) * t701;
t656 = Ifges(5,5) * t695 + Ifges(5,6) * t694 + Ifges(5,3) * t701;
t646 = Ifges(6,1) * t668 + Ifges(6,4) * t667 + Ifges(6,5) * t701;
t645 = Ifges(6,4) * t668 + Ifges(6,2) * t667 + Ifges(6,6) * t701;
t644 = Ifges(6,5) * t668 + Ifges(6,6) * t667 + Ifges(6,3) * t701;
t629 = Ifges(7,1) * t649 + Ifges(7,4) * t648 + Ifges(7,5) * t699;
t628 = Ifges(7,4) * t649 + Ifges(7,2) * t648 + Ifges(7,6) * t699;
t627 = Ifges(7,5) * t649 + Ifges(7,6) * t648 + Ifges(7,3) * t699;
t599 = mrSges(7,2) * t614 - mrSges(7,3) * t607 + Ifges(7,1) * t623 + Ifges(7,4) * t622 + Ifges(7,5) * t683 + t627 * t648 - t628 * t699;
t598 = -mrSges(7,1) * t614 + mrSges(7,3) * t608 + Ifges(7,4) * t623 + Ifges(7,2) * t622 + Ifges(7,6) * t683 - t627 * t649 + t629 * t699;
t585 = mrSges(6,2) * t624 - mrSges(6,3) * t611 + Ifges(6,1) * t638 + Ifges(6,4) * t637 + Ifges(6,5) * t686 - pkin(9) * t597 - t598 * t724 + t599 * t728 + t644 * t667 - t645 * t701;
t584 = -mrSges(6,1) * t624 + mrSges(6,3) * t612 + Ifges(6,4) * t638 + Ifges(6,2) * t637 + Ifges(6,6) * t686 - pkin(5) * t740 + pkin(9) * t745 + t728 * t598 + t724 * t599 - t668 * t644 + t701 * t646;
t571 = mrSges(5,2) * t634 - mrSges(5,3) * t625 + Ifges(5,1) * t663 + Ifges(5,4) * t662 + Ifges(5,5) * t686 - qJ(5) * t591 - t584 * t720 + t585 * t722 + t656 * t694 - t657 * t701;
t570 = -mrSges(5,1) * t634 + mrSges(5,3) * t626 + Ifges(5,4) * t663 + Ifges(5,2) * t662 + Ifges(5,6) * t686 - pkin(4) * t736 + qJ(5) * t746 + t722 * t584 + t720 * t585 - t695 * t656 + t701 * t658;
t569 = (-Ifges(5,3) - Ifges(6,3)) * t686 - t704 * t679 - t695 * t657 + Ifges(4,4) * t690 + t694 * t658 - Ifges(7,3) * t683 - mrSges(4,1) * t688 + Ifges(4,2) * t689 + qJD(3) * t681 - Ifges(5,6) * t662 - Ifges(5,5) * t663 + t667 * t646 - t668 * t645 + t648 * t629 - t649 * t628 + mrSges(4,3) * t652 - Ifges(6,6) * t637 - Ifges(6,5) * t638 - Ifges(7,6) * t622 - Ifges(7,5) * t623 - mrSges(5,1) * t625 + mrSges(5,2) * t626 - mrSges(6,1) * t611 + mrSges(6,2) * t612 + mrSges(7,2) * t608 - mrSges(7,1) * t607 - pkin(5) * t597 - pkin(3) * t583 + Ifges(4,6) * qJDD(3) - pkin(4) * t591;
t564 = mrSges(4,2) * t688 - mrSges(4,3) * t651 + Ifges(4,1) * t690 + Ifges(4,4) * t689 + Ifges(4,5) * qJDD(3) - pkin(8) * t583 - qJD(3) * t680 - t570 * t725 + t571 * t729 + t679 * t703;
t563 = mrSges(3,2) * t702 - mrSges(3,3) * t691 - pkin(7) * t577 + qJDD(1) * t743 + t730 * t564 - t726 * t569 + t707 * t756;
t562 = mrSges(2,1) * g(3) - pkin(1) * t568 + mrSges(2,3) * t709 - pkin(2) * t577 - mrSges(3,1) * t691 + mrSges(3,2) * t692 - t725 * t571 - t729 * t570 - pkin(3) * t734 - pkin(8) * t747 - mrSges(4,1) * t651 + mrSges(4,2) * t652 - Ifges(4,5) * t690 - Ifges(4,6) * t689 - Ifges(4,3) * qJDD(3) - t704 * t680 + t703 * t681 + (Ifges(2,6) - t741) * qJDD(1) + (-t721 * t742 + t723 * t743 + Ifges(2,5)) * t733;
t561 = -mrSges(3,1) * t702 + mrSges(3,3) * t692 - pkin(2) * t737 + pkin(7) * t748 + qJDD(1) * t742 + t726 * t564 + t730 * t569 - t707 * t757;
t560 = -mrSges(2,2) * g(3) - mrSges(2,3) * t708 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t733 - qJ(2) * t568 - t561 * t721 + t563 * t723;
t1 = [-m(1) * g(1) + t750; -m(1) * g(2) + t758; (-m(1) - m(2)) * g(3) + t568; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t758 + t731 * t560 - t727 * t562; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t750 + t727 * t560 + t731 * t562; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t708 - mrSges(2,2) * t709 + t721 * t563 + t723 * t561 + pkin(1) * (-mrSges(3,2) * t753 + t735) + qJ(2) * t749;];
tauB  = t1;
