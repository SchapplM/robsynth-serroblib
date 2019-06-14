% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRRP7
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
% Datum: 2019-05-08 05:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRRP7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP7_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP7_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP7_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP7_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP7_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 05:22:49
% EndTime: 2019-05-08 05:23:32
% DurationCPUTime: 27.46s
% Computational Cost: add. (450038->377), mult. (956380->473), div. (0->0), fcn. (771872->12), ass. (0->155)
t800 = Ifges(6,1) + Ifges(7,1);
t794 = Ifges(6,4) + Ifges(7,4);
t793 = Ifges(6,5) + Ifges(7,5);
t799 = Ifges(6,2) + Ifges(7,2);
t798 = Ifges(6,6) + Ifges(7,6);
t797 = Ifges(6,3) + Ifges(7,3);
t754 = cos(pkin(6));
t796 = t754 * g(3);
t795 = -mrSges(6,2) - mrSges(7,2);
t753 = sin(pkin(6));
t758 = sin(qJ(2));
t791 = t753 * t758;
t763 = cos(qJ(2));
t790 = t753 * t763;
t789 = t754 * t758;
t788 = t754 * t763;
t759 = sin(qJ(1));
t764 = cos(qJ(1));
t745 = t759 * g(1) - g(2) * t764;
t765 = qJD(1) ^ 2;
t736 = pkin(8) * t753 * t765 + qJDD(1) * pkin(1) + t745;
t746 = -g(1) * t764 - g(2) * t759;
t779 = qJDD(1) * t753;
t737 = -pkin(1) * t765 + pkin(8) * t779 + t746;
t782 = t736 * t789 + t763 * t737;
t713 = -g(3) * t791 + t782;
t750 = qJD(1) * t754 + qJD(2);
t781 = qJD(1) * t753;
t776 = t758 * t781;
t734 = mrSges(3,1) * t750 - mrSges(3,3) * t776;
t738 = (-mrSges(3,1) * t763 + mrSges(3,2) * t758) * t781;
t741 = -qJD(2) * t776 + t763 * t779;
t749 = qJDD(1) * t754 + qJDD(2);
t739 = (-pkin(2) * t763 - pkin(9) * t758) * t781;
t748 = t750 ^ 2;
t780 = qJD(1) * t763;
t697 = -t748 * pkin(2) + t749 * pkin(9) + (-g(3) * t758 + t739 * t780) * t753 + t782;
t740 = (qJD(2) * t780 + qJDD(1) * t758) * t753;
t698 = -t741 * pkin(2) - t740 * pkin(9) - t796 + (-t736 + (pkin(2) * t758 - pkin(9) * t763) * t750 * qJD(1)) * t753;
t757 = sin(qJ(3));
t762 = cos(qJ(3));
t663 = -t757 * t697 + t762 * t698;
t728 = t750 * t762 - t757 * t776;
t711 = qJD(3) * t728 + t740 * t762 + t749 * t757;
t729 = t750 * t757 + t762 * t776;
t733 = qJDD(3) - t741;
t775 = t753 * t780;
t744 = qJD(3) - t775;
t652 = (t728 * t744 - t711) * pkin(10) + (t728 * t729 + t733) * pkin(3) + t663;
t664 = t762 * t697 + t757 * t698;
t710 = -qJD(3) * t729 - t740 * t757 + t749 * t762;
t720 = pkin(3) * t744 - pkin(10) * t729;
t727 = t728 ^ 2;
t655 = -pkin(3) * t727 + pkin(10) * t710 - t720 * t744 + t664;
t756 = sin(qJ(4));
t761 = cos(qJ(4));
t650 = t756 * t652 + t761 * t655;
t716 = t728 * t756 + t729 * t761;
t677 = -qJD(4) * t716 + t710 * t761 - t711 * t756;
t715 = t728 * t761 - t729 * t756;
t691 = -mrSges(5,1) * t715 + mrSges(5,2) * t716;
t743 = qJD(4) + t744;
t703 = mrSges(5,1) * t743 - mrSges(5,3) * t716;
t732 = qJDD(4) + t733;
t692 = -pkin(4) * t715 - pkin(11) * t716;
t742 = t743 ^ 2;
t645 = -pkin(4) * t742 + pkin(11) * t732 + t692 * t715 + t650;
t712 = -g(3) * t790 + t736 * t788 - t758 * t737;
t696 = -t749 * pkin(2) - t748 * pkin(9) + t739 * t776 - t712;
t661 = -t710 * pkin(3) - t727 * pkin(10) + t729 * t720 + t696;
t678 = qJD(4) * t715 + t710 * t756 + t711 * t761;
t648 = (-t715 * t743 - t678) * pkin(11) + (t716 * t743 - t677) * pkin(4) + t661;
t755 = sin(qJ(5));
t760 = cos(qJ(5));
t640 = -t755 * t645 + t760 * t648;
t700 = -t716 * t755 + t743 * t760;
t660 = qJD(5) * t700 + t678 * t760 + t732 * t755;
t676 = qJDD(5) - t677;
t701 = t716 * t760 + t743 * t755;
t680 = -mrSges(7,1) * t700 + mrSges(7,2) * t701;
t681 = -mrSges(6,1) * t700 + mrSges(6,2) * t701;
t714 = qJD(5) - t715;
t683 = -mrSges(6,2) * t714 + mrSges(6,3) * t700;
t637 = -0.2e1 * qJD(6) * t701 + (t700 * t714 - t660) * qJ(6) + (t700 * t701 + t676) * pkin(5) + t640;
t682 = -mrSges(7,2) * t714 + mrSges(7,3) * t700;
t778 = m(7) * t637 + t676 * mrSges(7,1) + t714 * t682;
t629 = m(6) * t640 + t676 * mrSges(6,1) + t714 * t683 + (-t680 - t681) * t701 + (-mrSges(6,3) - mrSges(7,3)) * t660 + t778;
t641 = t760 * t645 + t755 * t648;
t659 = -qJD(5) * t701 - t678 * t755 + t732 * t760;
t684 = pkin(5) * t714 - qJ(6) * t701;
t699 = t700 ^ 2;
t639 = -pkin(5) * t699 + qJ(6) * t659 + 0.2e1 * qJD(6) * t700 - t684 * t714 + t641;
t777 = m(7) * t639 + t659 * mrSges(7,3) + t700 * t680;
t685 = mrSges(7,1) * t714 - mrSges(7,3) * t701;
t783 = -mrSges(6,1) * t714 + mrSges(6,3) * t701 - t685;
t632 = m(6) * t641 + t659 * mrSges(6,3) + t676 * t795 + t700 * t681 + t714 * t783 + t777;
t771 = -t629 * t755 + t760 * t632;
t625 = m(5) * t650 - mrSges(5,2) * t732 + mrSges(5,3) * t677 + t691 * t715 - t703 * t743 + t771;
t649 = t652 * t761 - t756 * t655;
t702 = -mrSges(5,2) * t743 + mrSges(5,3) * t715;
t644 = -pkin(4) * t732 - pkin(11) * t742 + t716 * t692 - t649;
t642 = -pkin(5) * t659 - qJ(6) * t699 + t684 * t701 + qJDD(6) + t644;
t770 = m(7) * t642 - t659 * mrSges(7,1) - t700 * t682;
t767 = -m(6) * t644 + t659 * mrSges(6,1) + t660 * t795 + t700 * t683 + t701 * t783 - t770;
t634 = m(5) * t649 + t732 * mrSges(5,1) - t678 * mrSges(5,3) - t716 * t691 + t743 * t702 + t767;
t618 = t756 * t625 + t761 * t634;
t717 = -mrSges(4,1) * t728 + mrSges(4,2) * t729;
t718 = -mrSges(4,2) * t744 + mrSges(4,3) * t728;
t616 = m(4) * t663 + mrSges(4,1) * t733 - mrSges(4,3) * t711 - t717 * t729 + t718 * t744 + t618;
t719 = mrSges(4,1) * t744 - mrSges(4,3) * t729;
t772 = t761 * t625 - t634 * t756;
t617 = m(4) * t664 - mrSges(4,2) * t733 + mrSges(4,3) * t710 + t717 * t728 - t719 * t744 + t772;
t773 = -t616 * t757 + t762 * t617;
t608 = m(3) * t713 - mrSges(3,2) * t749 + mrSges(3,3) * t741 - t734 * t750 + t738 * t775 + t773;
t611 = t762 * t616 + t757 * t617;
t724 = -t753 * t736 - t796;
t735 = -mrSges(3,2) * t750 + mrSges(3,3) * t775;
t610 = m(3) * t724 - t741 * mrSges(3,1) + t740 * mrSges(3,2) + (t734 * t758 - t735 * t763) * t781 + t611;
t627 = t760 * t629 + t755 * t632;
t768 = m(5) * t661 - t677 * mrSges(5,1) + mrSges(5,2) * t678 - t715 * t702 + t703 * t716 + t627;
t766 = -m(4) * t696 + t710 * mrSges(4,1) - mrSges(4,2) * t711 + t728 * t718 - t719 * t729 - t768;
t622 = m(3) * t712 + mrSges(3,1) * t749 - mrSges(3,3) * t740 + t735 * t750 - t738 * t776 + t766;
t598 = t608 * t789 - t610 * t753 + t622 * t788;
t596 = m(2) * t745 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t765 + t598;
t603 = t763 * t608 - t622 * t758;
t602 = m(2) * t746 - mrSges(2,1) * t765 - qJDD(1) * mrSges(2,2) + t603;
t787 = t764 * t596 + t759 * t602;
t786 = t798 * t700 + t793 * t701 + t797 * t714;
t785 = -t799 * t700 - t794 * t701 - t798 * t714;
t784 = t794 * t700 + t800 * t701 + t793 * t714;
t597 = t608 * t791 + t754 * t610 + t622 * t790;
t774 = -t596 * t759 + t764 * t602;
t619 = -mrSges(6,1) * t644 + mrSges(6,3) * t641 - mrSges(7,1) * t642 + mrSges(7,3) * t639 - pkin(5) * t770 + qJ(6) * t777 + (-qJ(6) * t685 + t784) * t714 + (-pkin(5) * t685 - t786) * t701 + (-mrSges(7,2) * qJ(6) + t798) * t676 + (-mrSges(7,2) * pkin(5) + t794) * t660 + t799 * t659;
t635 = -t660 * mrSges(7,3) - t701 * t680 + t778;
t626 = mrSges(6,2) * t644 + mrSges(7,2) * t642 - mrSges(6,3) * t640 - mrSges(7,3) * t637 - qJ(6) * t635 + t794 * t659 + t800 * t660 + t793 * t676 + t786 * t700 + t785 * t714;
t687 = Ifges(5,5) * t716 + Ifges(5,6) * t715 + Ifges(5,3) * t743;
t688 = Ifges(5,4) * t716 + Ifges(5,2) * t715 + Ifges(5,6) * t743;
t604 = mrSges(5,2) * t661 - mrSges(5,3) * t649 + Ifges(5,1) * t678 + Ifges(5,4) * t677 + Ifges(5,5) * t732 - pkin(11) * t627 - t619 * t755 + t626 * t760 + t687 * t715 - t688 * t743;
t689 = Ifges(5,1) * t716 + Ifges(5,4) * t715 + Ifges(5,5) * t743;
t612 = -mrSges(5,1) * t661 - mrSges(6,1) * t640 - mrSges(7,1) * t637 + mrSges(6,2) * t641 + mrSges(7,2) * t639 + mrSges(5,3) * t650 + Ifges(5,4) * t678 + Ifges(5,2) * t677 + Ifges(5,6) * t732 - pkin(4) * t627 - pkin(5) * t635 - t716 * t687 + t743 * t689 + t785 * t701 + t784 * t700 - t797 * t676 - t793 * t660 - t798 * t659;
t704 = Ifges(4,5) * t729 + Ifges(4,6) * t728 + Ifges(4,3) * t744;
t706 = Ifges(4,1) * t729 + Ifges(4,4) * t728 + Ifges(4,5) * t744;
t594 = -mrSges(4,1) * t696 + mrSges(4,3) * t664 + Ifges(4,4) * t711 + Ifges(4,2) * t710 + Ifges(4,6) * t733 - pkin(3) * t768 + pkin(10) * t772 + t756 * t604 + t761 * t612 - t729 * t704 + t744 * t706;
t705 = Ifges(4,4) * t729 + Ifges(4,2) * t728 + Ifges(4,6) * t744;
t599 = mrSges(4,2) * t696 - mrSges(4,3) * t663 + Ifges(4,1) * t711 + Ifges(4,4) * t710 + Ifges(4,5) * t733 - pkin(10) * t618 + t604 * t761 - t612 * t756 + t704 * t728 - t705 * t744;
t721 = Ifges(3,3) * t750 + (Ifges(3,5) * t758 + Ifges(3,6) * t763) * t781;
t722 = Ifges(3,6) * t750 + (Ifges(3,4) * t758 + Ifges(3,2) * t763) * t781;
t592 = mrSges(3,2) * t724 - mrSges(3,3) * t712 + Ifges(3,1) * t740 + Ifges(3,4) * t741 + Ifges(3,5) * t749 - pkin(9) * t611 - t594 * t757 + t599 * t762 + t721 * t775 - t722 * t750;
t723 = Ifges(3,5) * t750 + (Ifges(3,1) * t758 + Ifges(3,4) * t763) * t781;
t593 = -t721 * t776 - pkin(11) * t771 - pkin(4) * t767 - t760 * t619 + Ifges(3,6) * t749 + t750 * t723 - t755 * t626 - Ifges(5,3) * t732 - Ifges(4,3) * t733 + Ifges(3,4) * t740 + Ifges(3,2) * t741 - mrSges(3,1) * t724 + t728 * t706 - t729 * t705 + t715 * t689 - t716 * t688 - Ifges(4,5) * t711 + mrSges(3,3) * t713 - Ifges(4,6) * t710 - Ifges(5,6) * t677 - Ifges(5,5) * t678 - mrSges(4,1) * t663 + mrSges(4,2) * t664 - mrSges(5,1) * t649 + mrSges(5,2) * t650 - pkin(3) * t618 - pkin(2) * t611;
t769 = pkin(8) * t603 + t592 * t758 + t593 * t763;
t591 = Ifges(3,5) * t740 + Ifges(3,6) * t741 + Ifges(3,3) * t749 + mrSges(3,1) * t712 - mrSges(3,2) * t713 + t757 * t599 + t762 * t594 + pkin(2) * t766 + pkin(9) * t773 + (t722 * t758 - t723 * t763) * t781;
t590 = -mrSges(2,2) * g(3) - mrSges(2,3) * t745 + Ifges(2,5) * qJDD(1) - t765 * Ifges(2,6) + t763 * t592 - t758 * t593 + (-t597 * t753 - t598 * t754) * pkin(8);
t589 = mrSges(2,1) * g(3) + mrSges(2,3) * t746 + t765 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t597 - t753 * t591 + t754 * t769;
t1 = [-m(1) * g(1) + t774; -m(1) * g(2) + t787; (-m(1) - m(2)) * g(3) + t597; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t787 - t759 * t589 + t764 * t590; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t774 + t764 * t589 + t759 * t590; -mrSges(1,1) * g(2) + mrSges(2,1) * t745 + mrSges(1,2) * g(1) - mrSges(2,2) * t746 + Ifges(2,3) * qJDD(1) + pkin(1) * t598 + t754 * t591 + t753 * t769;];
tauB  = t1;
