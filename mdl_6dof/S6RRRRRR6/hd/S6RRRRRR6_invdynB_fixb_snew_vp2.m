% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 11:29
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRRR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR6_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR6_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR6_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR6_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR6_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 11:11:08
% EndTime: 2019-05-08 11:12:30
% DurationCPUTime: 59.09s
% Computational Cost: add. (994885->399), mult. (2102764->515), div. (0->0), fcn. (1726862->14), ass. (0->164)
t772 = cos(pkin(6));
t807 = t772 * g(3);
t771 = sin(pkin(6));
t777 = sin(qJ(2));
t806 = t771 * t777;
t783 = cos(qJ(2));
t805 = t771 * t783;
t804 = t772 * t777;
t803 = t772 * t783;
t778 = sin(qJ(1));
t784 = cos(qJ(1));
t763 = t778 * g(1) - g(2) * t784;
t785 = qJD(1) ^ 2;
t754 = pkin(8) * t771 * t785 + qJDD(1) * pkin(1) + t763;
t764 = -g(1) * t784 - g(2) * t778;
t798 = qJDD(1) * t771;
t755 = -pkin(1) * t785 + pkin(8) * t798 + t764;
t801 = t754 * t804 + t783 * t755;
t728 = -g(3) * t806 + t801;
t768 = qJD(1) * t772 + qJD(2);
t800 = qJD(1) * t771;
t797 = t777 * t800;
t752 = mrSges(3,1) * t768 - mrSges(3,3) * t797;
t756 = (-mrSges(3,1) * t783 + mrSges(3,2) * t777) * t800;
t759 = -qJD(2) * t797 + t783 * t798;
t767 = qJDD(1) * t772 + qJDD(2);
t757 = (-pkin(2) * t783 - pkin(9) * t777) * t800;
t766 = t768 ^ 2;
t799 = qJD(1) * t783;
t711 = -t766 * pkin(2) + t767 * pkin(9) + (-g(3) * t777 + t757 * t799) * t771 + t801;
t758 = (qJD(2) * t799 + qJDD(1) * t777) * t771;
t712 = -t759 * pkin(2) - t758 * pkin(9) - t807 + (-t754 + (pkin(2) * t777 - pkin(9) * t783) * t768 * qJD(1)) * t771;
t776 = sin(qJ(3));
t782 = cos(qJ(3));
t683 = -t776 * t711 + t782 * t712;
t746 = t768 * t782 - t776 * t797;
t726 = qJD(3) * t746 + t758 * t782 + t767 * t776;
t747 = t768 * t776 + t782 * t797;
t751 = qJDD(3) - t759;
t796 = t771 * t799;
t762 = qJD(3) - t796;
t668 = (t746 * t762 - t726) * pkin(10) + (t746 * t747 + t751) * pkin(3) + t683;
t684 = t782 * t711 + t776 * t712;
t725 = -qJD(3) * t747 - t758 * t776 + t767 * t782;
t737 = pkin(3) * t762 - pkin(10) * t747;
t745 = t746 ^ 2;
t675 = -pkin(3) * t745 + pkin(10) * t725 - t737 * t762 + t684;
t775 = sin(qJ(4));
t781 = cos(qJ(4));
t666 = t775 * t668 + t781 * t675;
t733 = t746 * t775 + t747 * t781;
t692 = -t733 * qJD(4) + t725 * t781 - t775 * t726;
t732 = t746 * t781 - t775 * t747;
t705 = -mrSges(5,1) * t732 + mrSges(5,2) * t733;
t761 = qJD(4) + t762;
t717 = mrSges(5,1) * t761 - mrSges(5,3) * t733;
t750 = qJDD(4) + t751;
t706 = -pkin(4) * t732 - pkin(11) * t733;
t760 = t761 ^ 2;
t658 = -pkin(4) * t760 + pkin(11) * t750 + t706 * t732 + t666;
t727 = -g(3) * t805 + t754 * t803 - t777 * t755;
t710 = -t767 * pkin(2) - t766 * pkin(9) + t757 * t797 - t727;
t680 = -t725 * pkin(3) - t745 * pkin(10) + t747 * t737 + t710;
t693 = qJD(4) * t732 + t725 * t775 + t726 * t781;
t661 = (-t732 * t761 - t693) * pkin(11) + (t733 * t761 - t692) * pkin(4) + t680;
t774 = sin(qJ(5));
t780 = cos(qJ(5));
t653 = -t774 * t658 + t780 * t661;
t714 = -t733 * t774 + t761 * t780;
t679 = qJD(5) * t714 + t693 * t780 + t750 * t774;
t691 = qJDD(5) - t692;
t715 = t733 * t780 + t761 * t774;
t731 = qJD(5) - t732;
t651 = (t714 * t731 - t679) * pkin(12) + (t714 * t715 + t691) * pkin(5) + t653;
t654 = t780 * t658 + t774 * t661;
t678 = -qJD(5) * t715 - t693 * t774 + t750 * t780;
t700 = pkin(5) * t731 - pkin(12) * t715;
t713 = t714 ^ 2;
t652 = -pkin(5) * t713 + pkin(12) * t678 - t700 * t731 + t654;
t773 = sin(qJ(6));
t779 = cos(qJ(6));
t649 = t651 * t779 - t652 * t773;
t695 = t714 * t779 - t715 * t773;
t664 = qJD(6) * t695 + t678 * t773 + t679 * t779;
t696 = t714 * t773 + t715 * t779;
t676 = -mrSges(7,1) * t695 + mrSges(7,2) * t696;
t729 = qJD(6) + t731;
t681 = -mrSges(7,2) * t729 + mrSges(7,3) * t695;
t689 = qJDD(6) + t691;
t647 = m(7) * t649 + mrSges(7,1) * t689 - mrSges(7,3) * t664 - t676 * t696 + t681 * t729;
t650 = t651 * t773 + t652 * t779;
t663 = -qJD(6) * t696 + t678 * t779 - t679 * t773;
t682 = mrSges(7,1) * t729 - mrSges(7,3) * t696;
t648 = m(7) * t650 - mrSges(7,2) * t689 + mrSges(7,3) * t663 + t676 * t695 - t682 * t729;
t639 = t779 * t647 + t773 * t648;
t697 = -mrSges(6,1) * t714 + mrSges(6,2) * t715;
t698 = -mrSges(6,2) * t731 + mrSges(6,3) * t714;
t637 = m(6) * t653 + mrSges(6,1) * t691 - mrSges(6,3) * t679 - t697 * t715 + t698 * t731 + t639;
t699 = mrSges(6,1) * t731 - mrSges(6,3) * t715;
t791 = -t647 * t773 + t779 * t648;
t638 = m(6) * t654 - mrSges(6,2) * t691 + mrSges(6,3) * t678 + t697 * t714 - t699 * t731 + t791;
t792 = -t637 * t774 + t780 * t638;
t632 = m(5) * t666 - mrSges(5,2) * t750 + mrSges(5,3) * t692 + t705 * t732 - t717 * t761 + t792;
t665 = t668 * t781 - t775 * t675;
t716 = -mrSges(5,2) * t761 + mrSges(5,3) * t732;
t657 = -pkin(4) * t750 - pkin(11) * t760 + t733 * t706 - t665;
t655 = -pkin(5) * t678 - pkin(12) * t713 + t700 * t715 + t657;
t789 = m(7) * t655 - t663 * mrSges(7,1) + mrSges(7,2) * t664 - t695 * t681 + t682 * t696;
t787 = -m(6) * t657 + t678 * mrSges(6,1) - mrSges(6,2) * t679 + t714 * t698 - t699 * t715 - t789;
t643 = m(5) * t665 + mrSges(5,1) * t750 - mrSges(5,3) * t693 - t705 * t733 + t716 * t761 + t787;
t624 = t775 * t632 + t781 * t643;
t734 = -mrSges(4,1) * t746 + mrSges(4,2) * t747;
t735 = -mrSges(4,2) * t762 + mrSges(4,3) * t746;
t622 = m(4) * t683 + mrSges(4,1) * t751 - mrSges(4,3) * t726 - t734 * t747 + t735 * t762 + t624;
t736 = mrSges(4,1) * t762 - mrSges(4,3) * t747;
t793 = t781 * t632 - t643 * t775;
t623 = m(4) * t684 - mrSges(4,2) * t751 + mrSges(4,3) * t725 + t734 * t746 - t736 * t762 + t793;
t794 = -t622 * t776 + t782 * t623;
t614 = m(3) * t728 - mrSges(3,2) * t767 + mrSges(3,3) * t759 - t752 * t768 + t756 * t796 + t794;
t617 = t782 * t622 + t776 * t623;
t741 = -t771 * t754 - t807;
t753 = -mrSges(3,2) * t768 + mrSges(3,3) * t796;
t616 = m(3) * t741 - t759 * mrSges(3,1) + t758 * mrSges(3,2) + (t752 * t777 - t753 * t783) * t800 + t617;
t633 = t780 * t637 + t774 * t638;
t788 = m(5) * t680 - t692 * mrSges(5,1) + mrSges(5,2) * t693 - t732 * t716 + t717 * t733 + t633;
t786 = -m(4) * t710 + t725 * mrSges(4,1) - mrSges(4,2) * t726 + t746 * t735 - t736 * t747 - t788;
t629 = m(3) * t727 + mrSges(3,1) * t767 - mrSges(3,3) * t758 + t753 * t768 - t756 * t797 + t786;
t605 = t614 * t804 - t616 * t771 + t629 * t803;
t603 = m(2) * t763 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t785 + t605;
t609 = t783 * t614 - t629 * t777;
t608 = m(2) * t764 - mrSges(2,1) * t785 - qJDD(1) * mrSges(2,2) + t609;
t802 = t784 * t603 + t778 * t608;
t604 = t614 * t806 + t772 * t616 + t629 * t805;
t795 = -t603 * t778 + t784 * t608;
t669 = Ifges(7,5) * t696 + Ifges(7,6) * t695 + Ifges(7,3) * t729;
t671 = Ifges(7,1) * t696 + Ifges(7,4) * t695 + Ifges(7,5) * t729;
t640 = -mrSges(7,1) * t655 + mrSges(7,3) * t650 + Ifges(7,4) * t664 + Ifges(7,2) * t663 + Ifges(7,6) * t689 - t669 * t696 + t671 * t729;
t670 = Ifges(7,4) * t696 + Ifges(7,2) * t695 + Ifges(7,6) * t729;
t641 = mrSges(7,2) * t655 - mrSges(7,3) * t649 + Ifges(7,1) * t664 + Ifges(7,4) * t663 + Ifges(7,5) * t689 + t669 * t695 - t670 * t729;
t685 = Ifges(6,5) * t715 + Ifges(6,6) * t714 + Ifges(6,3) * t731;
t687 = Ifges(6,1) * t715 + Ifges(6,4) * t714 + Ifges(6,5) * t731;
t625 = -mrSges(6,1) * t657 + mrSges(6,3) * t654 + Ifges(6,4) * t679 + Ifges(6,2) * t678 + Ifges(6,6) * t691 - pkin(5) * t789 + pkin(12) * t791 + t779 * t640 + t773 * t641 - t715 * t685 + t731 * t687;
t686 = Ifges(6,4) * t715 + Ifges(6,2) * t714 + Ifges(6,6) * t731;
t626 = mrSges(6,2) * t657 - mrSges(6,3) * t653 + Ifges(6,1) * t679 + Ifges(6,4) * t678 + Ifges(6,5) * t691 - pkin(12) * t639 - t640 * t773 + t641 * t779 + t685 * t714 - t686 * t731;
t701 = Ifges(5,5) * t733 + Ifges(5,6) * t732 + Ifges(5,3) * t761;
t702 = Ifges(5,4) * t733 + Ifges(5,2) * t732 + Ifges(5,6) * t761;
t610 = mrSges(5,2) * t680 - mrSges(5,3) * t665 + Ifges(5,1) * t693 + Ifges(5,4) * t692 + Ifges(5,5) * t750 - pkin(11) * t633 - t625 * t774 + t626 * t780 + t701 * t732 - t702 * t761;
t703 = Ifges(5,1) * t733 + Ifges(5,4) * t732 + Ifges(5,5) * t761;
t618 = Ifges(5,4) * t693 + Ifges(5,2) * t692 + Ifges(5,6) * t750 - t733 * t701 + t761 * t703 - mrSges(5,1) * t680 + mrSges(5,3) * t666 - Ifges(6,5) * t679 - Ifges(6,6) * t678 - Ifges(6,3) * t691 - t715 * t686 + t714 * t687 - mrSges(6,1) * t653 + mrSges(6,2) * t654 - Ifges(7,5) * t664 - Ifges(7,6) * t663 - Ifges(7,3) * t689 - t696 * t670 + t695 * t671 - mrSges(7,1) * t649 + mrSges(7,2) * t650 - pkin(5) * t639 - pkin(4) * t633;
t719 = Ifges(4,5) * t747 + Ifges(4,6) * t746 + Ifges(4,3) * t762;
t721 = Ifges(4,1) * t747 + Ifges(4,4) * t746 + Ifges(4,5) * t762;
t599 = -mrSges(4,1) * t710 + mrSges(4,3) * t684 + Ifges(4,4) * t726 + Ifges(4,2) * t725 + Ifges(4,6) * t751 - pkin(3) * t788 + pkin(10) * t793 + t775 * t610 + t781 * t618 - t747 * t719 + t762 * t721;
t720 = Ifges(4,4) * t747 + Ifges(4,2) * t746 + Ifges(4,6) * t762;
t601 = mrSges(4,2) * t710 - mrSges(4,3) * t683 + Ifges(4,1) * t726 + Ifges(4,4) * t725 + Ifges(4,5) * t751 - pkin(10) * t624 + t610 * t781 - t618 * t775 + t719 * t746 - t720 * t762;
t738 = Ifges(3,3) * t768 + (Ifges(3,5) * t777 + Ifges(3,6) * t783) * t800;
t739 = Ifges(3,6) * t768 + (Ifges(3,4) * t777 + Ifges(3,2) * t783) * t800;
t598 = mrSges(3,2) * t741 - mrSges(3,3) * t727 + Ifges(3,1) * t758 + Ifges(3,4) * t759 + Ifges(3,5) * t767 - pkin(9) * t617 - t599 * t776 + t601 * t782 + t738 * t796 - t739 * t768;
t740 = Ifges(3,5) * t768 + (Ifges(3,1) * t777 + Ifges(3,4) * t783) * t800;
t600 = -t738 * t797 - pkin(11) * t792 - pkin(2) * t617 - t774 * t626 - t780 * t625 + Ifges(3,6) * t767 + t768 * t740 - t747 * t720 - Ifges(5,3) * t750 - Ifges(4,3) * t751 + Ifges(3,4) * t758 + Ifges(3,2) * t759 - mrSges(3,1) * t741 + t746 * t721 + t732 * t703 - t733 * t702 - Ifges(4,6) * t725 - Ifges(4,5) * t726 + mrSges(3,3) * t728 - Ifges(5,6) * t692 - Ifges(5,5) * t693 - mrSges(4,1) * t683 + mrSges(4,2) * t684 - mrSges(5,1) * t665 + mrSges(5,2) * t666 - pkin(4) * t787 - pkin(3) * t624;
t790 = pkin(8) * t609 + t598 * t777 + t600 * t783;
t597 = Ifges(3,5) * t758 + Ifges(3,6) * t759 + Ifges(3,3) * t767 + mrSges(3,1) * t727 - mrSges(3,2) * t728 + t776 * t601 + t782 * t599 + pkin(2) * t786 + pkin(9) * t794 + (t739 * t777 - t740 * t783) * t800;
t596 = -mrSges(2,2) * g(3) - mrSges(2,3) * t763 + Ifges(2,5) * qJDD(1) - t785 * Ifges(2,6) + t783 * t598 - t777 * t600 + (-t604 * t771 - t605 * t772) * pkin(8);
t595 = mrSges(2,1) * g(3) + mrSges(2,3) * t764 + t785 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t604 - t771 * t597 + t772 * t790;
t1 = [-m(1) * g(1) + t795; -m(1) * g(2) + t802; (-m(1) - m(2)) * g(3) + t604; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t802 - t778 * t595 + t784 * t596; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t795 + t784 * t595 + t778 * t596; -mrSges(1,1) * g(2) + mrSges(2,1) * t763 + mrSges(1,2) * g(1) - mrSges(2,2) * t764 + Ifges(2,3) * qJDD(1) + pkin(1) * t605 + t772 * t597 + t771 * t790;];
tauB  = t1;
