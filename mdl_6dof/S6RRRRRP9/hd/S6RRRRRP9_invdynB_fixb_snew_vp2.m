% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRRP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 06:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRRP9_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP9_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP9_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP9_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP9_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP9_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP9_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 05:59:37
% EndTime: 2019-05-08 06:00:26
% DurationCPUTime: 28.72s
% Computational Cost: add. (476571->377), mult. (1013792->474), div. (0->0), fcn. (812500->12), ass. (0->155)
t809 = Ifges(6,1) + Ifges(7,1);
t804 = Ifges(6,4) + Ifges(7,4);
t803 = Ifges(6,5) + Ifges(7,5);
t808 = Ifges(6,2) + Ifges(7,2);
t802 = -Ifges(6,6) - Ifges(7,6);
t807 = -Ifges(6,3) - Ifges(7,3);
t763 = sin(pkin(6));
t768 = sin(qJ(2));
t773 = cos(qJ(2));
t790 = qJD(1) * qJD(2);
t752 = (-qJDD(1) * t773 + t768 * t790) * t763;
t806 = pkin(8) * t763;
t764 = cos(pkin(6));
t805 = t764 * g(3);
t801 = t763 * t768;
t800 = t763 * t773;
t799 = t764 * t768;
t798 = t764 * t773;
t769 = sin(qJ(1));
t774 = cos(qJ(1));
t756 = t769 * g(1) - t774 * g(2);
t775 = qJD(1) ^ 2;
t747 = qJDD(1) * pkin(1) + t775 * t806 + t756;
t757 = -t774 * g(1) - t769 * g(2);
t748 = -t775 * pkin(1) + qJDD(1) * t806 + t757;
t793 = t747 * t799 + t773 * t748;
t722 = -g(3) * t801 + t793;
t760 = t764 * qJD(1) + qJD(2);
t792 = qJD(1) * t763;
t787 = t768 * t792;
t745 = t760 * mrSges(3,1) - mrSges(3,3) * t787;
t749 = (-mrSges(3,1) * t773 + mrSges(3,2) * t768) * t792;
t759 = t764 * qJDD(1) + qJDD(2);
t750 = (-pkin(2) * t773 - pkin(9) * t768) * t792;
t758 = t760 ^ 2;
t791 = qJD(1) * t773;
t701 = -t758 * pkin(2) + t759 * pkin(9) + (-g(3) * t768 + t750 * t791) * t763 + t793;
t751 = (qJDD(1) * t768 + t773 * t790) * t763;
t702 = t752 * pkin(2) - t751 * pkin(9) - t805 + (-t747 + (pkin(2) * t768 - pkin(9) * t773) * t760 * qJD(1)) * t763;
t767 = sin(qJ(3));
t772 = cos(qJ(3));
t671 = t772 * t701 + t767 * t702;
t740 = t767 * t760 + t772 * t787;
t719 = -t740 * qJD(3) - t767 * t751 + t772 * t759;
t739 = t772 * t760 - t767 * t787;
t723 = -t739 * mrSges(4,1) + t740 * mrSges(4,2);
t786 = t763 * t791;
t755 = qJD(3) - t786;
t729 = t755 * mrSges(4,1) - t740 * mrSges(4,3);
t744 = qJDD(3) + t752;
t724 = -t739 * pkin(3) - t740 * pkin(10);
t753 = t755 ^ 2;
t665 = -t753 * pkin(3) + t744 * pkin(10) + t739 * t724 + t671;
t721 = -g(3) * t800 + t747 * t798 - t768 * t748;
t700 = -t759 * pkin(2) - t758 * pkin(9) + t750 * t787 - t721;
t720 = t739 * qJD(3) + t772 * t751 + t767 * t759;
t668 = (-t739 * t755 - t720) * pkin(10) + (t740 * t755 - t719) * pkin(3) + t700;
t766 = sin(qJ(4));
t771 = cos(qJ(4));
t654 = -t766 * t665 + t771 * t668;
t726 = -t766 * t740 + t771 * t755;
t686 = t726 * qJD(4) + t771 * t720 + t766 * t744;
t717 = qJDD(4) - t719;
t727 = t771 * t740 + t766 * t755;
t738 = qJD(4) - t739;
t651 = (t726 * t738 - t686) * pkin(11) + (t726 * t727 + t717) * pkin(4) + t654;
t655 = t771 * t665 + t766 * t668;
t685 = -t727 * qJD(4) - t766 * t720 + t771 * t744;
t710 = t738 * pkin(4) - t727 * pkin(11);
t725 = t726 ^ 2;
t653 = -t725 * pkin(4) + t685 * pkin(11) - t738 * t710 + t655;
t765 = sin(qJ(5));
t770 = cos(qJ(5));
t645 = t770 * t651 - t765 * t653;
t704 = t770 * t726 - t765 * t727;
t662 = t704 * qJD(5) + t765 * t685 + t770 * t686;
t705 = t765 * t726 + t770 * t727;
t681 = -t704 * mrSges(7,1) + t705 * mrSges(7,2);
t682 = -t704 * mrSges(6,1) + t705 * mrSges(6,2);
t736 = qJD(5) + t738;
t688 = -t736 * mrSges(6,2) + t704 * mrSges(6,3);
t712 = qJDD(5) + t717;
t642 = -0.2e1 * qJD(6) * t705 + (t704 * t736 - t662) * qJ(6) + (t704 * t705 + t712) * pkin(5) + t645;
t687 = -t736 * mrSges(7,2) + t704 * mrSges(7,3);
t789 = m(7) * t642 + t712 * mrSges(7,1) + t736 * t687;
t634 = m(6) * t645 + t712 * mrSges(6,1) + t736 * t688 + (-t681 - t682) * t705 + (-mrSges(6,3) - mrSges(7,3)) * t662 + t789;
t646 = t765 * t651 + t770 * t653;
t661 = -t705 * qJD(5) + t770 * t685 - t765 * t686;
t690 = t736 * mrSges(7,1) - t705 * mrSges(7,3);
t691 = t736 * mrSges(6,1) - t705 * mrSges(6,3);
t689 = t736 * pkin(5) - t705 * qJ(6);
t703 = t704 ^ 2;
t644 = -t703 * pkin(5) + t661 * qJ(6) + 0.2e1 * qJD(6) * t704 - t736 * t689 + t646;
t788 = m(7) * t644 + t661 * mrSges(7,3) + t704 * t681;
t637 = m(6) * t646 + t661 * mrSges(6,3) + t704 * t682 + (-t690 - t691) * t736 + (-mrSges(6,2) - mrSges(7,2)) * t712 + t788;
t632 = t770 * t634 + t765 * t637;
t706 = -t726 * mrSges(5,1) + t727 * mrSges(5,2);
t708 = -t738 * mrSges(5,2) + t726 * mrSges(5,3);
t629 = m(5) * t654 + t717 * mrSges(5,1) - t686 * mrSges(5,3) - t727 * t706 + t738 * t708 + t632;
t709 = t738 * mrSges(5,1) - t727 * mrSges(5,3);
t782 = -t765 * t634 + t770 * t637;
t630 = m(5) * t655 - t717 * mrSges(5,2) + t685 * mrSges(5,3) + t726 * t706 - t738 * t709 + t782;
t783 = -t766 * t629 + t771 * t630;
t625 = m(4) * t671 - t744 * mrSges(4,2) + t719 * mrSges(4,3) + t739 * t723 - t755 * t729 + t783;
t670 = -t767 * t701 + t772 * t702;
t728 = -t755 * mrSges(4,2) + t739 * mrSges(4,3);
t664 = -t744 * pkin(3) - t753 * pkin(10) + t740 * t724 - t670;
t656 = -t685 * pkin(4) - t725 * pkin(11) + t727 * t710 + t664;
t648 = -t661 * pkin(5) - t703 * qJ(6) + t705 * t689 + qJDD(6) + t656;
t781 = m(7) * t648 - t661 * mrSges(7,1) + t662 * mrSges(7,2) - t704 * t687 + t705 * t690;
t778 = m(6) * t656 - t661 * mrSges(6,1) + t662 * mrSges(6,2) - t704 * t688 + t705 * t691 + t781;
t776 = -m(5) * t664 + t685 * mrSges(5,1) - t686 * mrSges(5,2) + t726 * t708 - t727 * t709 - t778;
t639 = m(4) * t670 + t744 * mrSges(4,1) - t720 * mrSges(4,3) - t740 * t723 + t755 * t728 + t776;
t784 = t772 * t625 - t767 * t639;
t616 = m(3) * t722 - t759 * mrSges(3,2) - t752 * mrSges(3,3) - t760 * t745 + t749 * t786 + t784;
t619 = t767 * t625 + t772 * t639;
t733 = -t763 * t747 - t805;
t746 = -t760 * mrSges(3,2) + mrSges(3,3) * t786;
t618 = m(3) * t733 + t752 * mrSges(3,1) + t751 * mrSges(3,2) + (t745 * t768 - t746 * t773) * t792 + t619;
t626 = t771 * t629 + t766 * t630;
t777 = -m(4) * t700 + t719 * mrSges(4,1) - t720 * mrSges(4,2) + t739 * t728 - t740 * t729 - t626;
t622 = m(3) * t721 + t759 * mrSges(3,1) - t751 * mrSges(3,3) + t760 * t746 - t749 * t787 + t777;
t605 = t616 * t799 - t763 * t618 + t622 * t798;
t603 = m(2) * t756 + qJDD(1) * mrSges(2,1) - t775 * mrSges(2,2) + t605;
t609 = t773 * t616 - t768 * t622;
t608 = m(2) * t757 - t775 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t609;
t797 = t774 * t603 + t769 * t608;
t796 = t704 * t802 - t705 * t803 + t736 * t807;
t795 = -t704 * t808 - t705 * t804 + t736 * t802;
t794 = t804 * t704 + t705 * t809 + t803 * t736;
t604 = t616 * t801 + t764 * t618 + t622 * t800;
t785 = -t769 * t603 + t774 * t608;
t627 = -mrSges(6,1) * t656 + mrSges(6,3) * t646 - mrSges(7,1) * t648 + mrSges(7,3) * t644 - pkin(5) * t781 + qJ(6) * t788 + (-qJ(6) * t690 + t794) * t736 + (-qJ(6) * mrSges(7,2) - t802) * t712 + t796 * t705 + t804 * t662 + t808 * t661;
t640 = -t662 * mrSges(7,3) - t705 * t681 + t789;
t631 = mrSges(6,2) * t656 + mrSges(7,2) * t648 - mrSges(6,3) * t645 - mrSges(7,3) * t642 - qJ(6) * t640 + t804 * t661 + t662 * t809 - t796 * t704 + t803 * t712 + t795 * t736;
t692 = Ifges(5,5) * t727 + Ifges(5,6) * t726 + Ifges(5,3) * t738;
t694 = Ifges(5,1) * t727 + Ifges(5,4) * t726 + Ifges(5,5) * t738;
t611 = -mrSges(5,1) * t664 + mrSges(5,3) * t655 + Ifges(5,4) * t686 + Ifges(5,2) * t685 + Ifges(5,6) * t717 - pkin(4) * t778 + pkin(11) * t782 + t770 * t627 + t765 * t631 - t727 * t692 + t738 * t694;
t693 = Ifges(5,4) * t727 + Ifges(5,2) * t726 + Ifges(5,6) * t738;
t612 = mrSges(5,2) * t664 - mrSges(5,3) * t654 + Ifges(5,1) * t686 + Ifges(5,4) * t685 + Ifges(5,5) * t717 - pkin(11) * t632 - t765 * t627 + t770 * t631 + t726 * t692 - t738 * t693;
t713 = Ifges(4,5) * t740 + Ifges(4,6) * t739 + Ifges(4,3) * t755;
t714 = Ifges(4,4) * t740 + Ifges(4,2) * t739 + Ifges(4,6) * t755;
t601 = mrSges(4,2) * t700 - mrSges(4,3) * t670 + Ifges(4,1) * t720 + Ifges(4,4) * t719 + Ifges(4,5) * t744 - pkin(10) * t626 - t766 * t611 + t771 * t612 + t739 * t713 - t755 * t714;
t715 = Ifges(4,1) * t740 + Ifges(4,4) * t739 + Ifges(4,5) * t755;
t610 = t755 * t715 + Ifges(4,6) * t744 - t740 * t713 + t726 * t694 - t727 * t693 - Ifges(5,3) * t717 + Ifges(4,2) * t719 + Ifges(4,4) * t720 - mrSges(4,1) * t700 - Ifges(5,6) * t685 - Ifges(5,5) * t686 + mrSges(4,3) * t671 + mrSges(5,2) * t655 - mrSges(5,1) * t654 + mrSges(6,2) * t646 - mrSges(6,1) * t645 + mrSges(7,2) * t644 - mrSges(7,1) * t642 - pkin(5) * t640 - pkin(4) * t632 + t807 * t712 - pkin(3) * t626 + t802 * t661 - t803 * t662 + t794 * t704 + t795 * t705;
t730 = Ifges(3,3) * t760 + (Ifges(3,5) * t768 + Ifges(3,6) * t773) * t792;
t731 = Ifges(3,6) * t760 + (Ifges(3,4) * t768 + Ifges(3,2) * t773) * t792;
t599 = mrSges(3,2) * t733 - mrSges(3,3) * t721 + Ifges(3,1) * t751 - Ifges(3,4) * t752 + Ifges(3,5) * t759 - pkin(9) * t619 + t772 * t601 - t767 * t610 + t730 * t786 - t760 * t731;
t732 = Ifges(3,5) * t760 + (Ifges(3,1) * t768 + Ifges(3,4) * t773) * t792;
t600 = Ifges(3,4) * t751 - Ifges(3,2) * t752 + Ifges(3,6) * t759 - t730 * t787 + t760 * t732 - mrSges(3,1) * t733 + mrSges(3,3) * t722 - Ifges(4,5) * t720 - Ifges(4,6) * t719 - Ifges(4,3) * t744 - t740 * t714 + t739 * t715 - mrSges(4,1) * t670 + mrSges(4,2) * t671 - t766 * t612 - t771 * t611 - pkin(3) * t776 - pkin(10) * t783 - pkin(2) * t619;
t779 = pkin(8) * t609 + t599 * t768 + t600 * t773;
t598 = Ifges(3,5) * t751 - Ifges(3,6) * t752 + Ifges(3,3) * t759 + mrSges(3,1) * t721 - mrSges(3,2) * t722 + t767 * t601 + t772 * t610 + pkin(2) * t777 + pkin(9) * t784 + (t731 * t768 - t732 * t773) * t792;
t597 = -mrSges(2,2) * g(3) - mrSges(2,3) * t756 + Ifges(2,5) * qJDD(1) - t775 * Ifges(2,6) + t773 * t599 - t768 * t600 + (-t604 * t763 - t605 * t764) * pkin(8);
t596 = mrSges(2,1) * g(3) + mrSges(2,3) * t757 + t775 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t604 - t763 * t598 + t779 * t764;
t1 = [-m(1) * g(1) + t785; -m(1) * g(2) + t797; (-m(1) - m(2)) * g(3) + t604; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t797 - t769 * t596 + t774 * t597; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t785 + t774 * t596 + t769 * t597; -mrSges(1,1) * g(2) + mrSges(2,1) * t756 + mrSges(1,2) * g(1) - mrSges(2,2) * t757 + Ifges(2,3) * qJDD(1) + pkin(1) * t605 + t764 * t598 + t779 * t763;];
tauB  = t1;
