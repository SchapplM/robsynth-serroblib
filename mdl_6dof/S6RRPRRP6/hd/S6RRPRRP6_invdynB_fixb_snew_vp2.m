% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 18:11
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRRP6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP6_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP6_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP6_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP6_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP6_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:59:48
% EndTime: 2019-05-06 18:00:12
% DurationCPUTime: 22.01s
% Computational Cost: add. (329686->373), mult. (861033->474), div. (0->0), fcn. (676069->12), ass. (0->153)
t795 = -2 * qJD(3);
t794 = Ifges(6,1) + Ifges(7,1);
t787 = Ifges(6,4) - Ifges(7,5);
t793 = -Ifges(6,5) - Ifges(7,4);
t792 = Ifges(6,2) + Ifges(7,3);
t785 = Ifges(6,6) - Ifges(7,6);
t791 = -Ifges(6,3) - Ifges(7,2);
t744 = sin(pkin(11));
t746 = cos(pkin(11));
t750 = sin(qJ(2));
t753 = cos(qJ(2));
t745 = sin(pkin(6));
t774 = qJD(1) * t745;
t725 = (t744 * t750 - t746 * t753) * t774;
t772 = qJD(1) * qJD(2);
t734 = (qJDD(1) * t750 + t753 * t772) * t745;
t747 = cos(pkin(6));
t739 = t747 * qJDD(1) + qJDD(2);
t740 = t747 * qJD(1) + qJD(2);
t751 = sin(qJ(1));
t754 = cos(qJ(1));
t736 = t751 * g(1) - t754 * g(2);
t755 = qJD(1) ^ 2;
t789 = pkin(8) * t745;
t731 = qJDD(1) * pkin(1) + t755 * t789 + t736;
t737 = -t754 * g(1) - t751 * g(2);
t732 = -t755 * pkin(1) + qJDD(1) * t789 + t737;
t780 = t747 * t753;
t762 = t731 * t780 - t750 * t732;
t784 = t745 ^ 2 * t755;
t669 = t739 * pkin(2) - t734 * qJ(3) + (pkin(2) * t750 * t784 + (qJ(3) * qJD(1) * t740 - g(3)) * t745) * t753 + t762;
t781 = t747 * t750;
t783 = t745 * t750;
t699 = -g(3) * t783 + t731 * t781 + t753 * t732;
t769 = t750 * t774;
t728 = t740 * pkin(2) - qJ(3) * t769;
t735 = (qJDD(1) * t753 - t750 * t772) * t745;
t771 = t753 ^ 2 * t784;
t674 = -pkin(2) * t771 + t735 * qJ(3) - t740 * t728 + t699;
t726 = (t744 * t753 + t746 * t750) * t774;
t644 = t746 * t669 - t744 * t674 + t726 * t795;
t790 = cos(qJ(5));
t788 = -mrSges(6,3) - mrSges(7,2);
t782 = t745 * t753;
t645 = t744 * t669 + t746 * t674 + t725 * t795;
t700 = t725 * mrSges(4,1) + t726 * mrSges(4,2);
t705 = -t744 * t734 + t746 * t735;
t712 = t740 * mrSges(4,1) - t726 * mrSges(4,3);
t702 = t725 * pkin(3) - t726 * pkin(9);
t738 = t740 ^ 2;
t643 = -t738 * pkin(3) + t739 * pkin(9) - t725 * t702 + t645;
t716 = -t747 * g(3) - t745 * t731;
t686 = -t735 * pkin(2) - qJ(3) * t771 + t728 * t769 + qJDD(3) + t716;
t706 = t746 * t734 + t744 * t735;
t647 = (t725 * t740 - t706) * pkin(9) + (t726 * t740 - t705) * pkin(3) + t686;
t749 = sin(qJ(4));
t752 = cos(qJ(4));
t639 = t752 * t643 + t749 * t647;
t710 = t752 * t726 + t749 * t740;
t683 = -t710 * qJD(4) - t749 * t706 + t752 * t739;
t709 = -t749 * t726 + t752 * t740;
t687 = -t709 * mrSges(5,1) + t710 * mrSges(5,2);
t724 = qJD(4) + t725;
t693 = t724 * mrSges(5,1) - t710 * mrSges(5,3);
t704 = qJDD(4) - t705;
t688 = -t709 * pkin(4) - t710 * pkin(10);
t723 = t724 ^ 2;
t635 = -t723 * pkin(4) + t704 * pkin(10) + t709 * t688 + t639;
t642 = -t739 * pkin(3) - t738 * pkin(9) + t726 * t702 - t644;
t684 = t709 * qJD(4) + t752 * t706 + t749 * t739;
t637 = (-t709 * t724 - t684) * pkin(10) + (t710 * t724 - t683) * pkin(4) + t642;
t748 = sin(qJ(5));
t632 = t790 * t635 + t748 * t637;
t691 = t790 * t710 + t748 * t724;
t650 = t691 * qJD(5) + t748 * t684 - t790 * t704;
t708 = qJD(5) - t709;
t672 = t708 * mrSges(6,1) - t691 * mrSges(6,3);
t682 = qJDD(5) - t683;
t690 = t748 * t710 - t790 * t724;
t662 = t690 * pkin(5) - t691 * qJ(6);
t707 = t708 ^ 2;
t628 = -t707 * pkin(5) + t682 * qJ(6) + 0.2e1 * qJD(6) * t708 - t690 * t662 + t632;
t673 = -t708 * mrSges(7,1) + t691 * mrSges(7,2);
t770 = m(7) * t628 + t682 * mrSges(7,3) + t708 * t673;
t663 = t690 * mrSges(7,1) - t691 * mrSges(7,3);
t775 = -t690 * mrSges(6,1) - t691 * mrSges(6,2) - t663;
t624 = m(6) * t632 - t682 * mrSges(6,2) + t788 * t650 - t708 * t672 + t775 * t690 + t770;
t631 = -t748 * t635 + t790 * t637;
t651 = -t690 * qJD(5) + t790 * t684 + t748 * t704;
t671 = -t708 * mrSges(6,2) - t690 * mrSges(6,3);
t629 = -t682 * pkin(5) - t707 * qJ(6) + t691 * t662 + qJDD(6) - t631;
t670 = -t690 * mrSges(7,2) + t708 * mrSges(7,3);
t761 = -m(7) * t629 + t682 * mrSges(7,1) + t708 * t670;
t625 = m(6) * t631 + t682 * mrSges(6,1) + t788 * t651 + t708 * t671 + t775 * t691 + t761;
t764 = t790 * t624 - t748 * t625;
t617 = m(5) * t639 - t704 * mrSges(5,2) + t683 * mrSges(5,3) + t709 * t687 - t724 * t693 + t764;
t638 = -t749 * t643 + t752 * t647;
t692 = -t724 * mrSges(5,2) + t709 * mrSges(5,3);
t634 = -t704 * pkin(4) - t723 * pkin(10) + t710 * t688 - t638;
t630 = -0.2e1 * qJD(6) * t691 + (t690 * t708 - t651) * qJ(6) + (t691 * t708 + t650) * pkin(5) + t634;
t626 = m(7) * t630 + t650 * mrSges(7,1) - t651 * mrSges(7,3) + t690 * t670 - t691 * t673;
t756 = -m(6) * t634 - t650 * mrSges(6,1) - t651 * mrSges(6,2) - t690 * t671 - t691 * t672 - t626;
t622 = m(5) * t638 + t704 * mrSges(5,1) - t684 * mrSges(5,3) - t710 * t687 + t724 * t692 + t756;
t765 = t752 * t617 - t749 * t622;
t609 = m(4) * t645 - t739 * mrSges(4,2) + t705 * mrSges(4,3) - t725 * t700 - t740 * t712 + t765;
t711 = -t740 * mrSges(4,2) - t725 * mrSges(4,3);
t620 = t748 * t624 + t790 * t625;
t757 = -m(5) * t642 + t683 * mrSges(5,1) - t684 * mrSges(5,2) + t709 * t692 - t710 * t693 - t620;
t614 = m(4) * t644 + t739 * mrSges(4,1) - t706 * mrSges(4,3) - t726 * t700 + t740 * t711 + t757;
t605 = t744 * t609 + t746 * t614;
t698 = -g(3) * t782 + t762;
t768 = t753 * t774;
t730 = -t740 * mrSges(3,2) + mrSges(3,3) * t768;
t733 = (-mrSges(3,1) * t753 + mrSges(3,2) * t750) * t774;
t603 = m(3) * t698 + t739 * mrSges(3,1) - t734 * mrSges(3,3) + t740 * t730 - t733 * t769 + t605;
t729 = t740 * mrSges(3,1) - mrSges(3,3) * t769;
t766 = t746 * t609 - t744 * t614;
t604 = m(3) * t699 - t739 * mrSges(3,2) + t735 * mrSges(3,3) - t740 * t729 + t733 * t768 + t766;
t612 = t749 * t617 + t752 * t622;
t758 = m(4) * t686 - t705 * mrSges(4,1) + t706 * mrSges(4,2) + t725 * t711 + t726 * t712 + t612;
t611 = m(3) * t716 - t735 * mrSges(3,1) + t734 * mrSges(3,2) + (t729 * t750 - t730 * t753) * t774 + t758;
t591 = t603 * t780 + t604 * t781 - t745 * t611;
t589 = m(2) * t736 + qJDD(1) * mrSges(2,1) - t755 * mrSges(2,2) + t591;
t595 = -t750 * t603 + t753 * t604;
t594 = m(2) * t737 - t755 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t595;
t779 = t754 * t589 + t751 * t594;
t778 = t792 * t690 - t787 * t691 - t785 * t708;
t777 = t785 * t690 + t793 * t691 + t791 * t708;
t776 = -t787 * t690 + t794 * t691 - t793 * t708;
t590 = t603 * t782 + t604 * t783 + t747 * t611;
t767 = -t751 * t589 + t754 * t594;
t618 = -mrSges(6,1) * t634 - mrSges(7,1) * t630 + mrSges(7,2) * t628 + mrSges(6,3) * t632 - pkin(5) * t626 - t792 * t650 + t787 * t651 + t785 * t682 + t777 * t691 + t776 * t708;
t619 = mrSges(6,2) * t634 + mrSges(7,2) * t629 - mrSges(6,3) * t631 - mrSges(7,3) * t630 - qJ(6) * t626 - t787 * t650 + t794 * t651 - t682 * t793 + t777 * t690 + t778 * t708;
t675 = Ifges(5,5) * t710 + Ifges(5,6) * t709 + Ifges(5,3) * t724;
t676 = Ifges(5,4) * t710 + Ifges(5,2) * t709 + Ifges(5,6) * t724;
t597 = mrSges(5,2) * t642 - mrSges(5,3) * t638 + Ifges(5,1) * t684 + Ifges(5,4) * t683 + Ifges(5,5) * t704 - pkin(10) * t620 - t748 * t618 + t790 * t619 + t709 * t675 - t724 * t676;
t677 = Ifges(5,1) * t710 + Ifges(5,4) * t709 + Ifges(5,5) * t724;
t606 = Ifges(5,4) * t684 + Ifges(5,2) * t683 + Ifges(5,6) * t704 - t710 * t675 + t724 * t677 - mrSges(5,1) * t642 + mrSges(5,3) * t639 - mrSges(6,1) * t631 + mrSges(6,2) * t632 + mrSges(7,1) * t629 - mrSges(7,3) * t628 - pkin(5) * t761 - qJ(6) * t770 - pkin(4) * t620 + (pkin(5) * t663 + t778) * t691 + (qJ(6) * t663 - t776) * t690 + t791 * t682 + (pkin(5) * mrSges(7,2) + t793) * t651 + (qJ(6) * mrSges(7,2) + t785) * t650;
t694 = Ifges(4,5) * t726 - Ifges(4,6) * t725 + Ifges(4,3) * t740;
t695 = Ifges(4,4) * t726 - Ifges(4,2) * t725 + Ifges(4,6) * t740;
t587 = mrSges(4,2) * t686 - mrSges(4,3) * t644 + Ifges(4,1) * t706 + Ifges(4,4) * t705 + Ifges(4,5) * t739 - pkin(9) * t612 + t752 * t597 - t749 * t606 - t725 * t694 - t740 * t695;
t696 = Ifges(4,1) * t726 - Ifges(4,4) * t725 + Ifges(4,5) * t740;
t596 = Ifges(4,4) * t706 + Ifges(4,2) * t705 + Ifges(4,6) * t739 - t726 * t694 + t740 * t696 - mrSges(4,1) * t686 + mrSges(4,3) * t645 - Ifges(5,5) * t684 - Ifges(5,6) * t683 - Ifges(5,3) * t704 - t710 * t676 + t709 * t677 - mrSges(5,1) * t638 + mrSges(5,2) * t639 - t748 * t619 - t790 * t618 - pkin(4) * t756 - pkin(10) * t764 - pkin(3) * t612;
t713 = Ifges(3,3) * t740 + (Ifges(3,5) * t750 + Ifges(3,6) * t753) * t774;
t715 = Ifges(3,5) * t740 + (Ifges(3,1) * t750 + Ifges(3,4) * t753) * t774;
t584 = -mrSges(3,1) * t716 + mrSges(3,3) * t699 + Ifges(3,4) * t734 + Ifges(3,2) * t735 + Ifges(3,6) * t739 - pkin(2) * t758 + qJ(3) * t766 + t744 * t587 + t746 * t596 - t713 * t769 + t740 * t715;
t714 = Ifges(3,6) * t740 + (Ifges(3,4) * t750 + Ifges(3,2) * t753) * t774;
t585 = mrSges(3,2) * t716 - mrSges(3,3) * t698 + Ifges(3,1) * t734 + Ifges(3,4) * t735 + Ifges(3,5) * t739 - qJ(3) * t605 + t746 * t587 - t744 * t596 + t713 * t768 - t740 * t714;
t759 = pkin(8) * t595 + t584 * t753 + t585 * t750;
t586 = Ifges(3,5) * t734 + Ifges(3,6) * t735 + mrSges(3,1) * t698 - mrSges(3,2) * t699 + Ifges(4,5) * t706 + Ifges(4,6) * t705 + t726 * t695 + t725 * t696 + mrSges(4,1) * t644 - mrSges(4,2) * t645 + t749 * t597 + t752 * t606 + pkin(3) * t757 + pkin(9) * t765 + pkin(2) * t605 + (Ifges(3,3) + Ifges(4,3)) * t739 + (t714 * t750 - t715 * t753) * t774;
t583 = -mrSges(2,2) * g(3) - mrSges(2,3) * t736 + Ifges(2,5) * qJDD(1) - t755 * Ifges(2,6) - t750 * t584 + t753 * t585 + (-t590 * t745 - t591 * t747) * pkin(8);
t582 = mrSges(2,1) * g(3) + mrSges(2,3) * t737 + t755 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t590 - t745 * t586 + t759 * t747;
t1 = [-m(1) * g(1) + t767; -m(1) * g(2) + t779; (-m(1) - m(2)) * g(3) + t590; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t779 - t751 * t582 + t754 * t583; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t767 + t754 * t582 + t751 * t583; -mrSges(1,1) * g(2) + mrSges(2,1) * t736 + mrSges(1,2) * g(1) - mrSges(2,2) * t737 + Ifges(2,3) * qJDD(1) + pkin(1) * t591 + t747 * t586 + t759 * t745;];
tauB  = t1;
