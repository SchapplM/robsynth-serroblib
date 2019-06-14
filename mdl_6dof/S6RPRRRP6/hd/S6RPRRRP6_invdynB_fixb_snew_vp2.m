% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 01:38
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRRP6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP6_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP6_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP6_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP6_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP6_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:34:56
% EndTime: 2019-05-06 01:35:09
% DurationCPUTime: 11.29s
% Computational Cost: add. (163822->344), mult. (387036->416), div. (0->0), fcn. (293454->10), ass. (0->143)
t775 = Ifges(6,1) + Ifges(7,1);
t771 = Ifges(6,4) + Ifges(7,4);
t770 = Ifges(6,5) + Ifges(7,5);
t774 = Ifges(6,2) + Ifges(7,2);
t769 = -Ifges(6,6) - Ifges(7,6);
t773 = -Ifges(6,3) - Ifges(7,3);
t737 = qJD(1) ^ 2;
t727 = cos(pkin(10));
t772 = pkin(2) * t727;
t726 = sin(pkin(10));
t768 = mrSges(3,2) * t726;
t725 = t727 ^ 2;
t767 = t725 * t737;
t731 = sin(qJ(1));
t735 = cos(qJ(1));
t715 = -g(1) * t735 - g(2) * t731;
t711 = -pkin(1) * t737 + qJDD(1) * qJ(2) + t715;
t759 = qJD(1) * qJD(2);
t754 = -t727 * g(3) - 0.2e1 * t726 * t759;
t682 = (-pkin(7) * qJDD(1) + t737 * t772 - t711) * t726 + t754;
t698 = -g(3) * t726 + (t711 + 0.2e1 * t759) * t727;
t757 = qJDD(1) * t727;
t683 = -pkin(2) * t767 + pkin(7) * t757 + t698;
t730 = sin(qJ(3));
t734 = cos(qJ(3));
t657 = t730 * t682 + t734 * t683;
t761 = qJD(1) * t727;
t762 = qJD(1) * t726;
t709 = -t730 * t762 + t734 * t761;
t743 = t726 * t734 + t727 * t730;
t710 = t743 * qJD(1);
t690 = -mrSges(4,1) * t709 + mrSges(4,2) * t710;
t706 = t710 * qJD(3);
t758 = qJDD(1) * t726;
t695 = -t730 * t758 + t734 * t757 - t706;
t703 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t710;
t693 = -pkin(3) * t709 - pkin(8) * t710;
t736 = qJD(3) ^ 2;
t640 = -pkin(3) * t736 + qJDD(3) * pkin(8) + t693 * t709 + t657;
t724 = t726 ^ 2;
t714 = t731 * g(1) - t735 * g(2);
t748 = qJDD(2) - t714;
t694 = (-pkin(1) - t772) * qJDD(1) + (-qJ(2) + (-t724 - t725) * pkin(7)) * t737 + t748;
t760 = t709 * qJD(3);
t696 = qJDD(1) * t743 + t760;
t643 = (-t696 - t760) * pkin(8) + (-t695 + t706) * pkin(3) + t694;
t729 = sin(qJ(4));
t733 = cos(qJ(4));
t629 = -t729 * t640 + t733 * t643;
t700 = qJD(3) * t733 - t710 * t729;
t670 = qJD(4) * t700 + qJDD(3) * t729 + t696 * t733;
t692 = qJDD(4) - t695;
t701 = qJD(3) * t729 + t710 * t733;
t707 = qJD(4) - t709;
t625 = (t700 * t707 - t670) * pkin(9) + (t700 * t701 + t692) * pkin(4) + t629;
t630 = t733 * t640 + t729 * t643;
t669 = -qJD(4) * t701 + qJDD(3) * t733 - t696 * t729;
t681 = pkin(4) * t707 - pkin(9) * t701;
t699 = t700 ^ 2;
t627 = -pkin(4) * t699 + pkin(9) * t669 - t681 * t707 + t630;
t728 = sin(qJ(5));
t732 = cos(qJ(5));
t619 = t732 * t625 - t728 * t627;
t672 = t700 * t732 - t701 * t728;
t636 = qJD(5) * t672 + t669 * t728 + t670 * t732;
t673 = t700 * t728 + t701 * t732;
t653 = -mrSges(7,1) * t672 + mrSges(7,2) * t673;
t654 = -mrSges(6,1) * t672 + mrSges(6,2) * t673;
t705 = qJD(5) + t707;
t659 = -mrSges(6,2) * t705 + mrSges(6,3) * t672;
t689 = qJDD(5) + t692;
t616 = -0.2e1 * qJD(6) * t673 + (t672 * t705 - t636) * qJ(6) + (t672 * t673 + t689) * pkin(5) + t619;
t658 = -mrSges(7,2) * t705 + mrSges(7,3) * t672;
t756 = m(7) * t616 + t689 * mrSges(7,1) + t705 * t658;
t608 = m(6) * t619 + t689 * mrSges(6,1) + t705 * t659 + (-t653 - t654) * t673 + (-mrSges(6,3) - mrSges(7,3)) * t636 + t756;
t620 = t728 * t625 + t732 * t627;
t635 = -qJD(5) * t673 + t669 * t732 - t670 * t728;
t661 = mrSges(7,1) * t705 - mrSges(7,3) * t673;
t662 = mrSges(6,1) * t705 - mrSges(6,3) * t673;
t660 = pkin(5) * t705 - qJ(6) * t673;
t671 = t672 ^ 2;
t618 = -pkin(5) * t671 + qJ(6) * t635 + 0.2e1 * qJD(6) * t672 - t660 * t705 + t620;
t755 = m(7) * t618 + t635 * mrSges(7,3) + t672 * t653;
t611 = m(6) * t620 + t635 * mrSges(6,3) + t672 * t654 + (-t661 - t662) * t705 + (-mrSges(6,2) - mrSges(7,2)) * t689 + t755;
t606 = t732 * t608 + t728 * t611;
t674 = -mrSges(5,1) * t700 + mrSges(5,2) * t701;
t677 = -mrSges(5,2) * t707 + mrSges(5,3) * t700;
t603 = m(5) * t629 + mrSges(5,1) * t692 - mrSges(5,3) * t670 - t674 * t701 + t677 * t707 + t606;
t678 = mrSges(5,1) * t707 - mrSges(5,3) * t701;
t749 = -t608 * t728 + t732 * t611;
t604 = m(5) * t630 - mrSges(5,2) * t692 + mrSges(5,3) * t669 + t674 * t700 - t678 * t707 + t749;
t750 = -t603 * t729 + t733 * t604;
t597 = m(4) * t657 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t695 - qJD(3) * t703 + t690 * t709 + t750;
t656 = t682 * t734 - t730 * t683;
t702 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t709;
t639 = -qJDD(3) * pkin(3) - pkin(8) * t736 + t710 * t693 - t656;
t628 = -pkin(4) * t669 - pkin(9) * t699 + t701 * t681 + t639;
t622 = -pkin(5) * t635 - qJ(6) * t671 + t660 * t673 + qJDD(6) + t628;
t744 = m(7) * t622 - t635 * mrSges(7,1) + t636 * mrSges(7,2) - t672 * t658 + t673 * t661;
t740 = m(6) * t628 - t635 * mrSges(6,1) + mrSges(6,2) * t636 - t672 * t659 + t662 * t673 + t744;
t738 = -m(5) * t639 + t669 * mrSges(5,1) - mrSges(5,2) * t670 + t700 * t677 - t678 * t701 - t740;
t613 = m(4) * t656 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t696 + qJD(3) * t702 - t690 * t710 + t738;
t592 = t730 * t597 + t734 * t613;
t697 = -t726 * t711 + t754;
t742 = mrSges(3,3) * qJDD(1) + t737 * (-mrSges(3,1) * t727 + t768);
t590 = m(3) * t697 - t726 * t742 + t592;
t751 = t734 * t597 - t730 * t613;
t591 = m(3) * t698 + t727 * t742 + t751;
t752 = -t590 * t726 + t727 * t591;
t582 = m(2) * t715 - mrSges(2,1) * t737 - qJDD(1) * mrSges(2,2) + t752;
t708 = -qJDD(1) * pkin(1) - t737 * qJ(2) + t748;
t598 = t733 * t603 + t729 * t604;
t741 = m(4) * t694 - t695 * mrSges(4,1) + t696 * mrSges(4,2) - t709 * t702 + t710 * t703 + t598;
t739 = -m(3) * t708 + mrSges(3,1) * t757 - t741 + (t724 * t737 + t767) * mrSges(3,3);
t594 = t739 - t737 * mrSges(2,2) + m(2) * t714 + (mrSges(2,1) - t768) * qJDD(1);
t766 = t731 * t582 + t735 * t594;
t583 = t727 * t590 + t726 * t591;
t765 = t769 * t672 - t770 * t673 + t773 * t705;
t764 = -t774 * t672 - t771 * t673 + t769 * t705;
t763 = t771 * t672 + t775 * t673 + t770 * t705;
t753 = t735 * t582 - t594 * t731;
t747 = Ifges(3,1) * t726 + Ifges(3,4) * t727;
t746 = Ifges(3,4) * t726 + Ifges(3,2) * t727;
t745 = Ifges(3,5) * t726 + Ifges(3,6) * t727;
t713 = t745 * qJD(1);
t686 = Ifges(4,1) * t710 + Ifges(4,4) * t709 + Ifges(4,5) * qJD(3);
t685 = Ifges(4,4) * t710 + Ifges(4,2) * t709 + Ifges(4,6) * qJD(3);
t684 = Ifges(4,5) * t710 + Ifges(4,6) * t709 + Ifges(4,3) * qJD(3);
t665 = Ifges(5,1) * t701 + Ifges(5,4) * t700 + Ifges(5,5) * t707;
t664 = Ifges(5,4) * t701 + Ifges(5,2) * t700 + Ifges(5,6) * t707;
t663 = Ifges(5,5) * t701 + Ifges(5,6) * t700 + Ifges(5,3) * t707;
t614 = -t636 * mrSges(7,3) - t673 * t653 + t756;
t605 = mrSges(6,2) * t628 + mrSges(7,2) * t622 - mrSges(6,3) * t619 - mrSges(7,3) * t616 - qJ(6) * t614 + t771 * t635 + t775 * t636 - t765 * t672 + t770 * t689 + t764 * t705;
t599 = -mrSges(6,1) * t628 + mrSges(6,3) * t620 - mrSges(7,1) * t622 + mrSges(7,3) * t618 - pkin(5) * t744 + qJ(6) * t755 + (-qJ(6) * t661 + t763) * t705 + (-mrSges(7,2) * qJ(6) - t769) * t689 + t765 * t673 + t771 * t636 + t774 * t635;
t586 = mrSges(5,2) * t639 - mrSges(5,3) * t629 + Ifges(5,1) * t670 + Ifges(5,4) * t669 + Ifges(5,5) * t692 - pkin(9) * t606 - t599 * t728 + t605 * t732 + t663 * t700 - t664 * t707;
t585 = -mrSges(5,1) * t639 + mrSges(5,3) * t630 + Ifges(5,4) * t670 + Ifges(5,2) * t669 + Ifges(5,6) * t692 - pkin(4) * t740 + pkin(9) * t749 + t732 * t599 + t728 * t605 - t701 * t663 + t707 * t665;
t584 = t773 * t689 + Ifges(4,6) * qJDD(3) + t769 * t635 - t770 * t636 - t710 * t684 + t700 * t665 - t701 * t664 - Ifges(5,3) * t692 - mrSges(4,1) * t694 + Ifges(4,2) * t695 + Ifges(4,4) * t696 + qJD(3) * t686 - Ifges(5,6) * t669 - Ifges(5,5) * t670 + mrSges(4,3) * t657 - mrSges(5,1) * t629 + mrSges(5,2) * t630 - mrSges(6,1) * t619 + mrSges(6,2) * t620 + mrSges(7,2) * t618 - mrSges(7,1) * t616 - pkin(5) * t614 - pkin(4) * t606 - pkin(3) * t598 + t763 * t672 + t764 * t673;
t579 = mrSges(4,2) * t694 - mrSges(4,3) * t656 + Ifges(4,1) * t696 + Ifges(4,4) * t695 + Ifges(4,5) * qJDD(3) - pkin(8) * t598 - qJD(3) * t685 - t585 * t729 + t586 * t733 + t684 * t709;
t578 = mrSges(3,2) * t708 - mrSges(3,3) * t697 - pkin(7) * t592 + qJDD(1) * t747 + t734 * t579 - t730 * t584 + t713 * t761;
t577 = mrSges(2,1) * g(3) - pkin(1) * t583 + mrSges(2,3) * t715 - pkin(2) * t592 - mrSges(3,1) * t697 + mrSges(3,2) * t698 - t729 * t586 - t733 * t585 - pkin(3) * t738 - pkin(8) * t750 - Ifges(4,5) * t696 - Ifges(4,6) * t695 - Ifges(4,3) * qJDD(3) - t710 * t685 + t709 * t686 - mrSges(4,1) * t656 + mrSges(4,2) * t657 + (Ifges(2,6) - t745) * qJDD(1) + (-t726 * t746 + t727 * t747 + Ifges(2,5)) * t737;
t576 = -mrSges(3,1) * t708 + mrSges(3,3) * t698 - pkin(2) * t741 + pkin(7) * t751 + qJDD(1) * t746 + t730 * t579 + t734 * t584 - t713 * t762;
t575 = -mrSges(2,2) * g(3) - mrSges(2,3) * t714 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t737 - qJ(2) * t583 - t576 * t726 + t578 * t727;
t1 = [-m(1) * g(1) + t753; -m(1) * g(2) + t766; (-m(1) - m(2)) * g(3) + t583; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t766 + t735 * t575 - t731 * t577; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t753 + t731 * t575 + t735 * t577; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t714 - mrSges(2,2) * t715 + t726 * t578 + t727 * t576 + pkin(1) * (-mrSges(3,2) * t758 + t739) + qJ(2) * t752;];
tauB  = t1;
