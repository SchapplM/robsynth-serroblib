% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
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
% Datum: 2019-05-08 01:41
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRPR13_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR13_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR13_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR13_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR13_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR13_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR13_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 01:20:59
% EndTime: 2019-05-08 01:21:44
% DurationCPUTime: 19.97s
% Computational Cost: add. (320723->378), mult. (678573->474), div. (0->0), fcn. (535598->12), ass. (0->158)
t821 = Ifges(5,1) + Ifges(6,1);
t812 = Ifges(5,4) - Ifges(6,5);
t811 = Ifges(5,5) + Ifges(6,4);
t820 = Ifges(5,2) + Ifges(6,3);
t810 = Ifges(5,6) - Ifges(6,6);
t819 = -Ifges(5,3) - Ifges(6,2);
t769 = sin(pkin(6));
t774 = sin(qJ(2));
t778 = cos(qJ(2));
t796 = qJD(1) * qJD(2);
t756 = (-qJDD(1) * t778 + t774 * t796) * t769;
t770 = cos(pkin(6));
t766 = t770 * qJD(1) + qJD(2);
t773 = sin(qJ(3));
t777 = cos(qJ(3));
t798 = qJD(1) * t769;
t795 = t774 * t798;
t744 = t777 * t766 - t773 * t795;
t755 = (qJDD(1) * t774 + t778 * t796) * t769;
t765 = t770 * qJDD(1) + qJDD(2);
t723 = t744 * qJD(3) + t777 * t755 + t773 * t765;
t745 = t773 * t766 + t777 * t795;
t797 = qJD(1) * t778;
t794 = t769 * t797;
t761 = qJD(3) - t794;
t772 = sin(qJ(4));
t816 = cos(qJ(4));
t729 = t772 * t745 - t761 * t816;
t748 = qJDD(3) + t756;
t679 = -t729 * qJD(4) + t723 * t816 + t772 * t748;
t754 = (-pkin(2) * t778 - pkin(9) * t774) * t798;
t764 = t766 ^ 2;
t775 = sin(qJ(1));
t779 = cos(qJ(1));
t762 = t775 * g(1) - t779 * g(2);
t780 = qJD(1) ^ 2;
t815 = pkin(8) * t769;
t751 = qJDD(1) * pkin(1) + t780 * t815 + t762;
t763 = -t779 * g(1) - t775 * g(2);
t752 = -t780 * pkin(1) + qJDD(1) * t815 + t763;
t806 = t770 * t774;
t799 = t751 * t806 + t778 * t752;
t697 = -t764 * pkin(2) + t765 * pkin(9) + (-g(3) * t774 + t754 * t797) * t769 + t799;
t814 = t770 * g(3);
t698 = t756 * pkin(2) - t755 * pkin(9) - t814 + (-t751 + (pkin(2) * t774 - pkin(9) * t778) * t766 * qJD(1)) * t769;
t668 = -t773 * t697 + t777 * t698;
t727 = -t744 * pkin(3) - t745 * pkin(10);
t760 = t761 ^ 2;
t784 = t748 * pkin(3) + t760 * pkin(10) - t745 * t727 + t668;
t742 = qJD(4) - t744;
t809 = t729 * t742;
t818 = (-t679 + t809) * qJ(5) - t784;
t817 = 2 * qJD(5);
t813 = -mrSges(5,3) - mrSges(6,2);
t808 = t769 * t774;
t807 = t769 * t778;
t805 = t770 * t778;
t725 = -g(3) * t808 + t799;
t749 = t766 * mrSges(3,1) - mrSges(3,3) * t795;
t753 = (-mrSges(3,1) * t778 + mrSges(3,2) * t774) * t798;
t669 = t777 * t697 + t773 * t698;
t722 = -t745 * qJD(3) - t773 * t755 + t777 * t765;
t726 = -t744 * mrSges(4,1) + t745 * mrSges(4,2);
t732 = t761 * mrSges(4,1) - t745 * mrSges(4,3);
t665 = -t760 * pkin(3) + t748 * pkin(10) + t744 * t727 + t669;
t724 = -g(3) * t807 + t751 * t805 - t774 * t752;
t696 = -t765 * pkin(2) - t764 * pkin(9) + t754 * t795 - t724;
t667 = (-t744 * t761 - t723) * pkin(10) + (t745 * t761 - t722) * pkin(3) + t696;
t657 = t816 * t665 + t772 * t667;
t730 = t745 * t816 + t772 * t761;
t678 = t730 * qJD(4) + t772 * t723 - t748 * t816;
t709 = t742 * mrSges(5,1) - t730 * mrSges(5,3);
t720 = qJDD(4) - t722;
t703 = t729 * pkin(4) - t730 * qJ(5);
t741 = t742 ^ 2;
t653 = -t741 * pkin(4) + t720 * qJ(5) - t729 * t703 + t742 * t817 + t657;
t710 = -t742 * mrSges(6,1) + t730 * mrSges(6,2);
t656 = -t772 * t665 + t667 * t816;
t654 = -t720 * pkin(4) - t741 * qJ(5) + t730 * t703 + qJDD(5) - t656;
t648 = (-t679 - t809) * pkin(11) + (t729 * t730 - t720) * pkin(5) + t654;
t711 = -t742 * pkin(5) - t730 * pkin(11);
t728 = t729 ^ 2;
t649 = -t728 * pkin(5) + t678 * pkin(11) + t742 * t711 + t653;
t771 = sin(qJ(6));
t776 = cos(qJ(6));
t646 = t776 * t648 - t771 * t649;
t701 = t776 * t729 - t771 * t730;
t661 = t701 * qJD(6) + t771 * t678 + t776 * t679;
t702 = t771 * t729 + t776 * t730;
t675 = -t701 * mrSges(7,1) + t702 * mrSges(7,2);
t740 = qJD(6) - t742;
t682 = -t740 * mrSges(7,2) + t701 * mrSges(7,3);
t715 = qJDD(6) - t720;
t644 = m(7) * t646 + t715 * mrSges(7,1) - t661 * mrSges(7,3) - t702 * t675 + t740 * t682;
t647 = t771 * t648 + t776 * t649;
t660 = -t702 * qJD(6) + t776 * t678 - t771 * t679;
t683 = t740 * mrSges(7,1) - t702 * mrSges(7,3);
t645 = m(7) * t647 - t715 * mrSges(7,2) + t660 * mrSges(7,3) + t701 * t675 - t740 * t683;
t790 = -t771 * t644 + t776 * t645;
t787 = m(6) * t653 + t720 * mrSges(6,3) + t742 * t710 + t790;
t704 = t729 * mrSges(6,1) - t730 * mrSges(6,3);
t800 = -t729 * mrSges(5,1) - t730 * mrSges(5,2) - t704;
t635 = m(5) * t657 - t720 * mrSges(5,2) + t678 * t813 - t742 * t709 + t729 * t800 + t787;
t708 = -t742 * mrSges(5,2) - t729 * mrSges(5,3);
t637 = t776 * t644 + t771 * t645;
t707 = -t729 * mrSges(6,2) + t742 * mrSges(6,3);
t783 = -m(6) * t654 + t720 * mrSges(6,1) + t742 * t707 - t637;
t636 = m(5) * t656 + t720 * mrSges(5,1) + t679 * t813 + t742 * t708 + t730 * t800 + t783;
t791 = t816 * t635 - t772 * t636;
t632 = m(4) * t669 - t748 * mrSges(4,2) + t722 * mrSges(4,3) + t744 * t726 - t761 * t732 + t791;
t731 = -t761 * mrSges(4,2) + t744 * mrSges(4,3);
t655 = -0.2e1 * qJD(5) * t730 + (t730 * t742 + t678) * pkin(4) + t818;
t651 = -t728 * pkin(11) + (-pkin(4) - pkin(5)) * t678 + (-pkin(4) * t742 + t711 + t817) * t730 - t818;
t788 = -m(7) * t651 + t660 * mrSges(7,1) - t661 * mrSges(7,2) + t701 * t682 - t702 * t683;
t642 = m(6) * t655 + t678 * mrSges(6,1) - t679 * mrSges(6,3) + t729 * t707 - t730 * t710 + t788;
t781 = m(5) * t784 - t678 * mrSges(5,1) - t679 * mrSges(5,2) - t729 * t708 - t730 * t709 - t642;
t641 = m(4) * t668 + t748 * mrSges(4,1) - t723 * mrSges(4,3) - t745 * t726 + t761 * t731 + t781;
t792 = t777 * t632 - t773 * t641;
t623 = m(3) * t725 - t765 * mrSges(3,2) - t756 * mrSges(3,3) - t766 * t749 + t753 * t794 + t792;
t626 = t773 * t632 + t777 * t641;
t736 = -t769 * t751 - t814;
t750 = -t766 * mrSges(3,2) + mrSges(3,3) * t794;
t625 = m(3) * t736 + t756 * mrSges(3,1) + t755 * mrSges(3,2) + (t749 * t774 - t750 * t778) * t798 + t626;
t633 = t772 * t635 + t636 * t816;
t782 = -m(4) * t696 + t722 * mrSges(4,1) - t723 * mrSges(4,2) + t744 * t731 - t745 * t732 - t633;
t629 = m(3) * t724 + t765 * mrSges(3,1) - t755 * mrSges(3,3) + t766 * t750 - t753 * t795 + t782;
t612 = t623 * t806 - t769 * t625 + t629 * t805;
t610 = m(2) * t762 + qJDD(1) * mrSges(2,1) - t780 * mrSges(2,2) + t612;
t617 = t778 * t623 - t774 * t629;
t616 = m(2) * t763 - t780 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t617;
t804 = t779 * t610 + t775 * t616;
t803 = t729 * t820 - t730 * t812 - t742 * t810;
t802 = t729 * t810 - t730 * t811 + t742 * t819;
t801 = -t812 * t729 + t730 * t821 + t811 * t742;
t611 = t623 * t808 + t770 * t625 + t629 * t807;
t793 = -t775 * t610 + t779 * t616;
t670 = Ifges(7,5) * t702 + Ifges(7,6) * t701 + Ifges(7,3) * t740;
t672 = Ifges(7,1) * t702 + Ifges(7,4) * t701 + Ifges(7,5) * t740;
t638 = -mrSges(7,1) * t651 + mrSges(7,3) * t647 + Ifges(7,4) * t661 + Ifges(7,2) * t660 + Ifges(7,6) * t715 - t702 * t670 + t740 * t672;
t671 = Ifges(7,4) * t702 + Ifges(7,2) * t701 + Ifges(7,6) * t740;
t639 = mrSges(7,2) * t651 - mrSges(7,3) * t646 + Ifges(7,1) * t661 + Ifges(7,4) * t660 + Ifges(7,5) * t715 + t701 * t670 - t740 * t671;
t618 = mrSges(5,1) * t784 - mrSges(6,1) * t655 + mrSges(6,2) * t653 + mrSges(5,3) * t657 - pkin(4) * t642 - pkin(5) * t788 - pkin(11) * t790 - t776 * t638 - t771 * t639 - t678 * t820 + t812 * t679 + t810 * t720 + t802 * t730 + t801 * t742;
t619 = -mrSges(5,2) * t784 + mrSges(6,2) * t654 - mrSges(5,3) * t656 - mrSges(6,3) * t655 - pkin(11) * t637 - qJ(5) * t642 - t771 * t638 + t776 * t639 - t812 * t678 + t679 * t821 + t811 * t720 + t802 * t729 + t803 * t742;
t716 = Ifges(4,5) * t745 + Ifges(4,6) * t744 + Ifges(4,3) * t761;
t717 = Ifges(4,4) * t745 + Ifges(4,2) * t744 + Ifges(4,6) * t761;
t608 = mrSges(4,2) * t696 - mrSges(4,3) * t668 + Ifges(4,1) * t723 + Ifges(4,4) * t722 + Ifges(4,5) * t748 - pkin(10) * t633 - t772 * t618 + t619 * t816 + t744 * t716 - t761 * t717;
t718 = Ifges(4,1) * t745 + Ifges(4,4) * t744 + Ifges(4,5) * t761;
t613 = -qJ(5) * t787 + (qJ(5) * t704 - t801) * t729 + (pkin(4) * t704 + t803) * t730 + pkin(5) * t637 - pkin(4) * t783 - pkin(3) * t633 + t819 * t720 + mrSges(4,3) * t669 + Ifges(7,6) * t660 + Ifges(7,5) * t661 - mrSges(5,1) * t656 + mrSges(5,2) * t657 + Ifges(7,3) * t715 + Ifges(4,2) * t722 + Ifges(4,4) * t723 - t701 * t672 + t702 * t671 + mrSges(6,1) * t654 - mrSges(6,3) * t653 - mrSges(7,2) * t647 + mrSges(7,1) * t646 - mrSges(4,1) * t696 + (qJ(5) * mrSges(6,2) + t810) * t678 + (pkin(4) * mrSges(6,2) - t811) * t679 - t745 * t716 + Ifges(4,6) * t748 + t761 * t718;
t733 = Ifges(3,3) * t766 + (Ifges(3,5) * t774 + Ifges(3,6) * t778) * t798;
t734 = Ifges(3,6) * t766 + (Ifges(3,4) * t774 + Ifges(3,2) * t778) * t798;
t606 = mrSges(3,2) * t736 - mrSges(3,3) * t724 + Ifges(3,1) * t755 - Ifges(3,4) * t756 + Ifges(3,5) * t765 - pkin(9) * t626 + t777 * t608 - t773 * t613 + t733 * t794 - t766 * t734;
t735 = Ifges(3,5) * t766 + (Ifges(3,1) * t774 + Ifges(3,4) * t778) * t798;
t607 = Ifges(3,4) * t755 - Ifges(3,2) * t756 + Ifges(3,6) * t765 - t733 * t795 + t766 * t735 - mrSges(3,1) * t736 + mrSges(3,3) * t725 - Ifges(4,5) * t723 - Ifges(4,6) * t722 - Ifges(4,3) * t748 - t745 * t717 + t744 * t718 - mrSges(4,1) * t668 + mrSges(4,2) * t669 - t772 * t619 - t816 * t618 - pkin(3) * t781 - pkin(10) * t791 - pkin(2) * t626;
t785 = pkin(8) * t617 + t606 * t774 + t607 * t778;
t605 = Ifges(3,5) * t755 - Ifges(3,6) * t756 + Ifges(3,3) * t765 + mrSges(3,1) * t724 - mrSges(3,2) * t725 + t773 * t608 + t777 * t613 + pkin(2) * t782 + pkin(9) * t792 + (t734 * t774 - t735 * t778) * t798;
t604 = -mrSges(2,2) * g(3) - mrSges(2,3) * t762 + Ifges(2,5) * qJDD(1) - t780 * Ifges(2,6) + t778 * t606 - t774 * t607 + (-t611 * t769 - t612 * t770) * pkin(8);
t603 = mrSges(2,1) * g(3) + mrSges(2,3) * t763 + t780 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t611 - t769 * t605 + t770 * t785;
t1 = [-m(1) * g(1) + t793; -m(1) * g(2) + t804; (-m(1) - m(2)) * g(3) + t611; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t804 - t775 * t603 + t779 * t604; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t793 + t779 * t603 + t775 * t604; -mrSges(1,1) * g(2) + mrSges(2,1) * t762 + mrSges(1,2) * g(1) - mrSges(2,2) * t763 + Ifges(2,3) * qJDD(1) + pkin(1) * t612 + t770 * t605 + t769 * t785;];
tauB  = t1;
