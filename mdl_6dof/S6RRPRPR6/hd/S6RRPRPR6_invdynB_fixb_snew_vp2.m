% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-05-06 14:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRPR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR6_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR6_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR6_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR6_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR6_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 14:21:43
% EndTime: 2019-05-06 14:22:03
% DurationCPUTime: 17.78s
% Computational Cost: add. (255794->379), mult. (672317->473), div. (0->0), fcn. (523112->12), ass. (0->156)
t813 = Ifges(5,1) + Ifges(6,2);
t807 = Ifges(5,4) + Ifges(6,6);
t806 = Ifges(5,5) - Ifges(6,4);
t812 = Ifges(5,2) + Ifges(6,3);
t805 = Ifges(5,6) - Ifges(6,5);
t811 = -Ifges(5,3) - Ifges(6,1);
t763 = sin(pkin(6));
t768 = sin(qJ(2));
t771 = cos(qJ(2));
t791 = qJD(1) * qJD(2);
t749 = (qJDD(1) * t768 + t771 * t791) * t763;
t765 = cos(pkin(6));
t757 = qJDD(1) * t765 + qJDD(2);
t758 = qJD(1) * t765 + qJD(2);
t769 = sin(qJ(1));
t772 = cos(qJ(1));
t754 = t769 * g(1) - g(2) * t772;
t773 = qJD(1) ^ 2;
t808 = pkin(8) * t763;
t746 = qJDD(1) * pkin(1) + t773 * t808 + t754;
t755 = -g(1) * t772 - g(2) * t769;
t747 = -pkin(1) * t773 + qJDD(1) * t808 + t755;
t799 = t765 * t771;
t783 = t746 * t799 - t768 * t747;
t803 = t763 ^ 2 * t773;
t676 = t757 * pkin(2) - t749 * qJ(3) + (pkin(2) * t768 * t803 + (qJ(3) * qJD(1) * t758 - g(3)) * t763) * t771 + t783;
t800 = t765 * t768;
t802 = t763 * t768;
t714 = -g(3) * t802 + t746 * t800 + t771 * t747;
t793 = qJD(1) * t763;
t789 = t768 * t793;
t743 = pkin(2) * t758 - qJ(3) * t789;
t750 = (qJDD(1) * t771 - t768 * t791) * t763;
t790 = t771 ^ 2 * t803;
t679 = -pkin(2) * t790 + qJ(3) * t750 - t743 * t758 + t714;
t762 = sin(pkin(11));
t764 = cos(pkin(11));
t741 = (t762 * t771 + t764 * t768) * t793;
t659 = -0.2e1 * qJD(3) * t741 + t764 * t676 - t762 * t679;
t810 = -2 * qJD(5);
t809 = cos(qJ(4));
t767 = sin(qJ(4));
t725 = t741 * t767 - t758 * t809;
t788 = t771 * t793;
t740 = -t762 * t789 + t764 * t788;
t739 = qJD(4) - t740;
t804 = t725 * t739;
t801 = t763 * t771;
t660 = 0.2e1 * qJD(3) * t740 + t762 * t676 + t764 * t679;
t715 = -mrSges(4,1) * t740 + mrSges(4,2) * t741;
t720 = -t749 * t762 + t750 * t764;
t728 = mrSges(4,1) * t758 - mrSges(4,3) * t741;
t716 = -pkin(3) * t740 - pkin(9) * t741;
t756 = t758 ^ 2;
t658 = -pkin(3) * t756 + pkin(9) * t757 + t716 * t740 + t660;
t732 = -t765 * g(3) - t763 * t746;
t696 = -t750 * pkin(2) - qJ(3) * t790 + t743 * t789 + qJDD(3) + t732;
t721 = t749 * t764 + t750 * t762;
t662 = (-t740 * t758 - t721) * pkin(9) + (t741 * t758 - t720) * pkin(3) + t696;
t653 = -t767 * t658 + t662 * t809;
t693 = -t725 * qJD(4) + t721 * t809 + t767 * t757;
t726 = t741 * t809 + t767 * t758;
t698 = mrSges(5,1) * t725 + mrSges(5,2) * t726;
t703 = mrSges(6,1) * t725 - mrSges(6,3) * t739;
t705 = -mrSges(5,2) * t739 - mrSges(5,3) * t725;
t719 = qJDD(4) - t720;
t697 = pkin(4) * t725 - qJ(5) * t726;
t738 = t739 ^ 2;
t651 = -t719 * pkin(4) - t738 * qJ(5) + t726 * t697 + qJDD(5) - t653;
t646 = (t725 * t726 - t719) * pkin(10) + (t693 + t804) * pkin(5) + t651;
t692 = qJD(4) * t726 + t721 * t767 - t757 * t809;
t707 = pkin(5) * t726 - pkin(10) * t739;
t724 = t725 ^ 2;
t657 = -t757 * pkin(3) - t756 * pkin(9) + t741 * t716 - t659;
t775 = (-t693 + t804) * qJ(5) + t657 + (pkin(4) * t739 + t810) * t726;
t649 = -t724 * pkin(5) - t726 * t707 + (pkin(4) + pkin(10)) * t692 + t775;
t766 = sin(qJ(6));
t770 = cos(qJ(6));
t644 = t646 * t770 - t649 * t766;
t701 = t725 * t770 - t739 * t766;
t665 = qJD(6) * t701 + t692 * t766 + t719 * t770;
t702 = t725 * t766 + t739 * t770;
t670 = -mrSges(7,1) * t701 + mrSges(7,2) * t702;
t723 = qJD(6) + t726;
t677 = -mrSges(7,2) * t723 + mrSges(7,3) * t701;
t691 = qJDD(6) + t693;
t642 = m(7) * t644 + mrSges(7,1) * t691 - mrSges(7,3) * t665 - t670 * t702 + t677 * t723;
t645 = t646 * t766 + t649 * t770;
t664 = -qJD(6) * t702 + t692 * t770 - t719 * t766;
t678 = mrSges(7,1) * t723 - mrSges(7,3) * t702;
t643 = m(7) * t645 - mrSges(7,2) * t691 + mrSges(7,3) * t664 + t670 * t701 - t678 * t723;
t634 = t770 * t642 + t766 * t643;
t699 = -mrSges(6,2) * t725 - mrSges(6,3) * t726;
t779 = -m(6) * t651 - t693 * mrSges(6,1) - t726 * t699 - t634;
t632 = m(5) * t653 - t693 * mrSges(5,3) - t726 * t698 + (-t703 + t705) * t739 + (mrSges(5,1) - mrSges(6,2)) * t719 + t779;
t654 = t809 * t658 + t767 * t662;
t706 = mrSges(5,1) * t739 - mrSges(5,3) * t726;
t778 = -t738 * pkin(4) + t719 * qJ(5) - t725 * t697 + t654;
t650 = t739 * t810 - t778;
t704 = mrSges(6,1) * t726 + mrSges(6,2) * t739;
t648 = -t692 * pkin(5) - t724 * pkin(10) + ((2 * qJD(5)) + t707) * t739 + t778;
t780 = -m(7) * t648 + t664 * mrSges(7,1) - t665 * mrSges(7,2) + t701 * t677 - t702 * t678;
t776 = -m(6) * t650 + t719 * mrSges(6,3) + t739 * t704 - t780;
t639 = m(5) * t654 - t719 * mrSges(5,2) - t739 * t706 + (-t698 - t699) * t725 + (-mrSges(5,3) - mrSges(6,1)) * t692 + t776;
t785 = -t632 * t767 + t809 * t639;
t625 = m(4) * t660 - mrSges(4,2) * t757 + mrSges(4,3) * t720 + t715 * t740 - t728 * t758 + t785;
t727 = -mrSges(4,2) * t758 + mrSges(4,3) * t740;
t652 = t692 * pkin(4) + t775;
t797 = -t766 * t642 + t770 * t643;
t782 = -m(6) * t652 + t692 * mrSges(6,2) + t725 * t703 - t797;
t774 = -m(5) * t657 - t725 * t705 - t692 * mrSges(5,1) + (t704 - t706) * t726 + (-mrSges(5,2) + mrSges(6,3)) * t693 + t782;
t630 = m(4) * t659 + t757 * mrSges(4,1) - t721 * mrSges(4,3) - t741 * t715 + t758 * t727 + t774;
t622 = t762 * t625 + t764 * t630;
t713 = -g(3) * t801 + t783;
t745 = -mrSges(3,2) * t758 + mrSges(3,3) * t788;
t748 = (-mrSges(3,1) * t771 + mrSges(3,2) * t768) * t793;
t620 = m(3) * t713 + mrSges(3,1) * t757 - mrSges(3,3) * t749 + t745 * t758 - t748 * t789 + t622;
t744 = mrSges(3,1) * t758 - mrSges(3,3) * t789;
t786 = t764 * t625 - t630 * t762;
t621 = m(3) * t714 - mrSges(3,2) * t757 + mrSges(3,3) * t750 - t744 * t758 + t748 * t788 + t786;
t628 = t809 * t632 + t767 * t639;
t777 = m(4) * t696 - t720 * mrSges(4,1) + t721 * mrSges(4,2) - t740 * t727 + t741 * t728 + t628;
t627 = m(3) * t732 - t750 * mrSges(3,1) + t749 * mrSges(3,2) + (t744 * t768 - t745 * t771) * t793 + t777;
t607 = t620 * t799 + t621 * t800 - t627 * t763;
t605 = m(2) * t754 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t773 + t607;
t611 = -t620 * t768 + t771 * t621;
t610 = m(2) * t755 - mrSges(2,1) * t773 - qJDD(1) * mrSges(2,2) + t611;
t798 = t772 * t605 + t769 * t610;
t796 = t805 * t725 - t806 * t726 + t811 * t739;
t795 = t812 * t725 - t807 * t726 - t805 * t739;
t794 = -t807 * t725 + t813 * t726 + t806 * t739;
t606 = t620 * t801 + t621 * t802 + t765 * t627;
t787 = -t605 * t769 + t772 * t610;
t633 = -t693 * mrSges(6,3) - t726 * t704 - t782;
t666 = Ifges(7,5) * t702 + Ifges(7,6) * t701 + Ifges(7,3) * t723;
t668 = Ifges(7,1) * t702 + Ifges(7,4) * t701 + Ifges(7,5) * t723;
t635 = -mrSges(7,1) * t648 + mrSges(7,3) * t645 + Ifges(7,4) * t665 + Ifges(7,2) * t664 + Ifges(7,6) * t691 - t666 * t702 + t668 * t723;
t667 = Ifges(7,4) * t702 + Ifges(7,2) * t701 + Ifges(7,6) * t723;
t636 = mrSges(7,2) * t648 - mrSges(7,3) * t644 + Ifges(7,1) * t665 + Ifges(7,4) * t664 + Ifges(7,5) * t691 + t666 * t701 - t667 * t723;
t613 = -mrSges(5,1) * t657 - mrSges(6,1) * t650 + mrSges(6,2) * t652 + mrSges(5,3) * t654 - pkin(4) * t633 - pkin(5) * t780 - pkin(10) * t797 - t770 * t635 - t766 * t636 - t692 * t812 + t807 * t693 + t805 * t719 + t796 * t726 + t794 * t739;
t614 = mrSges(6,1) * t651 + mrSges(7,1) * t644 + mrSges(5,2) * t657 - mrSges(7,2) * t645 - mrSges(5,3) * t653 - mrSges(6,3) * t652 + Ifges(7,5) * t665 + Ifges(7,6) * t664 + Ifges(7,3) * t691 + pkin(5) * t634 - qJ(5) * t633 + t702 * t667 - t701 * t668 + t795 * t739 + t796 * t725 + t806 * t719 + t813 * t693 - t807 * t692;
t709 = Ifges(4,5) * t741 + Ifges(4,6) * t740 + Ifges(4,3) * t758;
t710 = Ifges(4,4) * t741 + Ifges(4,2) * t740 + Ifges(4,6) * t758;
t603 = mrSges(4,2) * t696 - mrSges(4,3) * t659 + Ifges(4,1) * t721 + Ifges(4,4) * t720 + Ifges(4,5) * t757 - pkin(9) * t628 - t767 * t613 + t614 * t809 + t740 * t709 - t758 * t710;
t711 = Ifges(4,1) * t741 + Ifges(4,4) * t740 + Ifges(4,5) * t758;
t612 = -pkin(3) * t628 + pkin(10) * t634 - qJ(5) * t776 + t766 * t635 - pkin(4) * (-t739 * t703 + t779) - t770 * t636 + t758 * t711 + Ifges(4,6) * t757 - t741 * t709 + Ifges(4,2) * t720 + Ifges(4,4) * t721 + mrSges(6,3) * t650 - mrSges(6,2) * t651 - mrSges(5,1) * t653 + mrSges(5,2) * t654 + mrSges(4,3) * t660 - mrSges(4,1) * t696 + t795 * t726 + (qJ(5) * t699 - t794) * t725 + (mrSges(6,2) * pkin(4) + t811) * t719 - t806 * t693 + (mrSges(6,1) * qJ(5) + t805) * t692;
t729 = Ifges(3,3) * t758 + (Ifges(3,5) * t768 + Ifges(3,6) * t771) * t793;
t731 = Ifges(3,5) * t758 + (Ifges(3,1) * t768 + Ifges(3,4) * t771) * t793;
t600 = -mrSges(3,1) * t732 + mrSges(3,3) * t714 + Ifges(3,4) * t749 + Ifges(3,2) * t750 + Ifges(3,6) * t757 - pkin(2) * t777 + qJ(3) * t786 + t762 * t603 + t764 * t612 - t729 * t789 + t758 * t731;
t730 = Ifges(3,6) * t758 + (Ifges(3,4) * t768 + Ifges(3,2) * t771) * t793;
t601 = mrSges(3,2) * t732 - mrSges(3,3) * t713 + Ifges(3,1) * t749 + Ifges(3,4) * t750 + Ifges(3,5) * t757 - qJ(3) * t622 + t603 * t764 - t612 * t762 + t729 * t788 - t730 * t758;
t781 = pkin(8) * t611 + t600 * t771 + t601 * t768;
t602 = Ifges(3,5) * t749 + Ifges(3,6) * t750 + mrSges(3,1) * t713 - mrSges(3,2) * t714 + Ifges(4,5) * t721 + Ifges(4,6) * t720 + t741 * t710 - t740 * t711 + mrSges(4,1) * t659 - mrSges(4,2) * t660 + t767 * t614 + t809 * t613 + pkin(3) * t774 + pkin(9) * t785 + pkin(2) * t622 + (Ifges(3,3) + Ifges(4,3)) * t757 + (t730 * t768 - t731 * t771) * t793;
t599 = -mrSges(2,2) * g(3) - mrSges(2,3) * t754 + Ifges(2,5) * qJDD(1) - t773 * Ifges(2,6) - t768 * t600 + t771 * t601 + (-t606 * t763 - t607 * t765) * pkin(8);
t598 = mrSges(2,1) * g(3) + mrSges(2,3) * t755 + t773 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t606 - t763 * t602 + t765 * t781;
t1 = [-m(1) * g(1) + t787; -m(1) * g(2) + t798; (-m(1) - m(2)) * g(3) + t606; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t798 - t769 * t598 + t772 * t599; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t787 + t772 * t598 + t769 * t599; -mrSges(1,1) * g(2) + mrSges(2,1) * t754 + mrSges(1,2) * g(1) - mrSges(2,2) * t755 + Ifges(2,3) * qJDD(1) + pkin(1) * t607 + t765 * t602 + t763 * t781;];
tauB  = t1;
