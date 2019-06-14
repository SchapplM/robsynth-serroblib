% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPPRP2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% Datum: 2019-05-06 09:13
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPPRP2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:09:25
% EndTime: 2019-05-06 09:09:37
% DurationCPUTime: 5.88s
% Computational Cost: add. (54957->352), mult. (126648->405), div. (0->0), fcn. (83214->8), ass. (0->138)
t812 = -2 * qJD(4);
t811 = Ifges(4,1) + Ifges(5,2);
t810 = -Ifges(5,1) - Ifges(4,3);
t809 = Ifges(6,1) + Ifges(7,1);
t800 = Ifges(6,4) + Ifges(7,4);
t799 = Ifges(4,5) - Ifges(5,4);
t798 = Ifges(6,5) + Ifges(7,5);
t808 = Ifges(4,2) + Ifges(5,3);
t807 = Ifges(6,2) + Ifges(7,2);
t797 = Ifges(4,6) - Ifges(5,5);
t796 = -Ifges(5,6) - Ifges(4,4);
t795 = Ifges(6,6) + Ifges(7,6);
t806 = Ifges(6,3) + Ifges(7,3);
t753 = sin(qJ(2));
t756 = cos(qJ(2));
t776 = qJD(1) * qJD(2);
t738 = qJDD(1) * t753 + t756 * t776;
t739 = qJDD(1) * t756 - t753 * t776;
t751 = sin(pkin(9));
t794 = cos(pkin(9));
t705 = t751 * t738 - t794 * t739;
t780 = qJD(1) * t756;
t781 = qJD(1) * t753;
t726 = t751 * t781 - t794 * t780;
t752 = sin(qJ(5));
t755 = cos(qJ(5));
t709 = qJD(2) * t755 + t726 * t752;
t668 = -qJD(5) * t709 - qJDD(2) * t752 + t705 * t755;
t805 = (mrSges(6,1) + mrSges(7,1)) * t668;
t727 = (t751 * t756 + t794 * t753) * qJD(1);
t696 = pkin(3) * t726 - qJ(4) * t727;
t758 = qJD(2) ^ 2;
t754 = sin(qJ(1));
t757 = cos(qJ(1));
t744 = -g(1) * t757 - g(2) * t754;
t759 = qJD(1) ^ 2;
t733 = -pkin(1) * t759 + qJDD(1) * pkin(7) + t744;
t792 = t753 * t733;
t802 = pkin(2) * t759;
t674 = qJDD(2) * pkin(2) - t738 * qJ(3) - t792 + (qJ(3) * t776 + t753 * t802 - g(3)) * t756;
t712 = -g(3) * t753 + t756 * t733;
t740 = qJD(2) * pkin(2) - qJ(3) * t781;
t750 = t756 ^ 2;
t675 = qJ(3) * t739 - qJD(2) * t740 - t750 * t802 + t712;
t786 = t751 * t674 + t794 * t675;
t804 = pkin(3) * t758 - qJDD(2) * qJ(4) + qJD(2) * t812 + t726 * t696 - t786;
t654 = -0.2e1 * qJD(3) * t727 + t794 * t674 - t751 * t675;
t708 = -qJD(2) * t752 + t726 * t755;
t724 = qJD(5) + t727;
t682 = -mrSges(7,2) * t724 + mrSges(7,3) * t708;
t683 = -mrSges(6,2) * t724 + mrSges(6,3) * t708;
t803 = -t805 - (t682 + t683) * t708;
t793 = t708 * t682;
t697 = mrSges(4,1) * t726 + mrSges(4,2) * t727;
t706 = t794 * t738 + t751 * t739;
t714 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t726;
t716 = mrSges(5,1) * t726 - qJD(2) * mrSges(5,3);
t651 = -qJDD(2) * pkin(3) - t758 * qJ(4) + t727 * t696 + qJDD(4) - t654;
t779 = qJD(2) * t726;
t645 = (t726 * t727 - qJDD(2)) * pkin(8) + (t706 + t779) * pkin(4) + t651;
t718 = pkin(4) * t727 - qJD(2) * pkin(8);
t725 = t726 ^ 2;
t743 = t754 * g(1) - t757 * g(2);
t767 = -qJDD(1) * pkin(1) - t743;
t681 = -t739 * pkin(2) + qJDD(3) + t740 * t781 + (-qJ(3) * t750 - pkin(7)) * t759 + t767;
t761 = (-t706 + t779) * qJ(4) + t681 + (qJD(2) * pkin(3) + t812) * t727;
t649 = -t725 * pkin(4) - t727 * t718 + (pkin(3) + pkin(8)) * t705 + t761;
t639 = t755 * t645 - t752 * t649;
t669 = qJD(5) * t708 + qJDD(2) * t755 + t705 * t752;
t676 = -mrSges(7,1) * t708 + mrSges(7,2) * t709;
t677 = -mrSges(6,1) * t708 + mrSges(6,2) * t709;
t704 = qJDD(5) + t706;
t636 = -0.2e1 * qJD(6) * t709 + (t708 * t724 - t669) * qJ(6) + (t708 * t709 + t704) * pkin(5) + t639;
t774 = m(7) * t636 + t704 * mrSges(7,1) + t724 * t682;
t631 = m(6) * t639 + t704 * mrSges(6,1) + t724 * t683 + (-t676 - t677) * t709 + (-mrSges(6,3) - mrSges(7,3)) * t669 + t774;
t640 = t752 * t645 + t755 * t649;
t685 = mrSges(7,1) * t724 - mrSges(7,3) * t709;
t686 = mrSges(6,1) * t724 - mrSges(6,3) * t709;
t684 = pkin(5) * t724 - qJ(6) * t709;
t707 = t708 ^ 2;
t638 = -pkin(5) * t707 + qJ(6) * t668 + 0.2e1 * qJD(6) * t708 - t684 * t724 + t640;
t773 = m(7) * t638 + t668 * mrSges(7,3) + t708 * t676;
t633 = m(6) * t640 + t668 * mrSges(6,3) + t708 * t677 + (-t685 - t686) * t724 + (-mrSges(6,2) - mrSges(7,2)) * t704 + t773;
t626 = t755 * t631 + t752 * t633;
t698 = -mrSges(5,2) * t726 - mrSges(5,3) * t727;
t764 = -m(5) * t651 - t706 * mrSges(5,1) - t727 * t698 - t626;
t623 = m(4) * t654 - t706 * mrSges(4,3) - t727 * t697 + (mrSges(4,1) - mrSges(5,2)) * qJDD(2) + (t714 - t716) * qJD(2) + t764;
t778 = qJD(3) * t726;
t721 = -0.2e1 * t778;
t655 = t721 + t786;
t715 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t727;
t650 = 0.2e1 * t778 + t804;
t717 = mrSges(5,1) * t727 + qJD(2) * mrSges(5,2);
t647 = -pkin(4) * t705 - pkin(8) * t725 + qJD(2) * t718 + t721 - t804;
t642 = -pkin(5) * t668 - qJ(6) * t707 + t684 * t709 + qJDD(6) + t647;
t772 = m(7) * t642 + t669 * mrSges(7,2) + t709 * t685;
t766 = -m(6) * t647 - t669 * mrSges(6,2) - t709 * t686 - t772;
t763 = -m(5) * t650 + qJDD(2) * mrSges(5,3) + qJD(2) * t717 - t766;
t629 = t763 - qJDD(2) * mrSges(4,2) + (-mrSges(4,3) - mrSges(5,1)) * t705 + (-t697 - t698) * t726 - qJD(2) * t715 + m(4) * t655 + t803;
t618 = t794 * t623 + t751 * t629;
t711 = -t756 * g(3) - t792;
t737 = (-mrSges(3,1) * t756 + mrSges(3,2) * t753) * qJD(1);
t742 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t780;
t616 = m(3) * t711 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t738 + qJD(2) * t742 - t737 * t781 + t618;
t741 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t781;
t769 = -t623 * t751 + t794 * t629;
t617 = m(3) * t712 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t739 - qJD(2) * t741 + t737 * t780 + t769;
t770 = -t616 * t753 + t756 * t617;
t611 = m(2) * t744 - mrSges(2,1) * t759 - qJDD(1) * mrSges(2,2) + t770;
t732 = -t759 * pkin(7) + t767;
t653 = t705 * pkin(3) + t761;
t790 = -t752 * t631 + t755 * t633;
t624 = m(5) * t653 - t705 * mrSges(5,2) - t706 * mrSges(5,3) - t726 * t716 - t727 * t717 + t790;
t762 = m(4) * t681 + t705 * mrSges(4,1) + mrSges(4,2) * t706 + t726 * t714 + t715 * t727 + t624;
t760 = -m(3) * t732 + t739 * mrSges(3,1) - mrSges(3,2) * t738 - t741 * t781 + t742 * t780 - t762;
t621 = m(2) * t743 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t759 + t760;
t791 = t754 * t611 + t757 * t621;
t612 = t756 * t616 + t753 * t617;
t789 = -t795 * t708 - t798 * t709 - t806 * t724;
t788 = t807 * t708 + t800 * t709 + t795 * t724;
t787 = -t800 * t708 - t809 * t709 - t798 * t724;
t784 = -t797 * qJD(2) + t808 * t726 + t796 * t727;
t783 = t810 * qJD(2) + t797 * t726 - t799 * t727;
t782 = t799 * qJD(2) + t796 * t726 + t811 * t727;
t771 = t757 * t611 - t621 * t754;
t730 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t753 + Ifges(3,4) * t756) * qJD(1);
t729 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t753 + Ifges(3,2) * t756) * qJD(1);
t728 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t753 + Ifges(3,6) * t756) * qJD(1);
t634 = -t669 * mrSges(7,3) - t709 * t676 + t774;
t625 = mrSges(6,2) * t647 + mrSges(7,2) * t642 - mrSges(6,3) * t639 - mrSges(7,3) * t636 - qJ(6) * t634 + t800 * t668 + t809 * t669 + t798 * t704 - t789 * t708 - t788 * t724;
t619 = -mrSges(6,1) * t647 + mrSges(6,3) * t640 - mrSges(7,1) * t642 + mrSges(7,3) * t638 - pkin(5) * (t772 - t793) + qJ(6) * t773 + (-qJ(6) * t685 - t787) * t724 + t789 * t709 + (-mrSges(7,2) * qJ(6) + t795) * t704 + t800 * t669 + (mrSges(7,1) * pkin(5) + t807) * t668;
t608 = mrSges(5,1) * t651 + mrSges(6,1) * t639 + mrSges(7,1) * t636 + mrSges(4,2) * t681 - mrSges(6,2) * t640 - mrSges(7,2) * t638 - mrSges(4,3) * t654 - mrSges(5,3) * t653 + pkin(4) * t626 + pkin(5) * t634 - qJ(4) * t624 + t783 * t726 + t788 * t709 + t787 * t708 + t811 * t706 + t796 * t705 + t806 * t704 + t798 * t669 + t795 * t668 + t799 * qJDD(2) + t784 * qJD(2);
t607 = -mrSges(4,1) * t681 + mrSges(4,3) * t655 - mrSges(5,1) * t650 + mrSges(5,2) * t653 - t752 * t625 - t755 * t619 - pkin(4) * (t766 - t803) - pkin(8) * t790 - pkin(3) * t624 + t783 * t727 - t796 * t706 - t808 * t705 + t797 * qJDD(2) + t782 * qJD(2);
t606 = (mrSges(5,2) * pkin(3) - Ifges(3,3) + t810) * qJDD(2) - qJ(4) * (-t708 * t683 + t763 - t793 - t805) - pkin(3) * (-qJD(2) * t716 + t764) + Ifges(2,6) * qJDD(1) + mrSges(2,1) * g(3) + (qJ(4) * t698 - t782) * t726 + t784 * t727 - t755 * t625 + t759 * Ifges(2,5) + t752 * t619 - Ifges(3,6) * t739 + mrSges(2,3) * t744 - Ifges(3,5) * t738 - mrSges(3,1) * t711 + mrSges(3,2) * t712 + mrSges(4,2) * t655 + mrSges(5,3) * t650 - mrSges(5,2) * t651 - mrSges(4,1) * t654 + pkin(8) * t626 - pkin(2) * t618 - pkin(1) * t612 + (-t729 * t753 + t730 * t756) * qJD(1) + (mrSges(5,1) * qJ(4) + t797) * t705 - t799 * t706;
t605 = mrSges(3,2) * t732 - mrSges(3,3) * t711 + Ifges(3,1) * t738 + Ifges(3,4) * t739 + Ifges(3,5) * qJDD(2) - qJ(3) * t618 - qJD(2) * t729 - t751 * t607 + t794 * t608 + t728 * t780;
t604 = -mrSges(3,1) * t732 + mrSges(3,3) * t712 + Ifges(3,4) * t738 + Ifges(3,2) * t739 + Ifges(3,6) * qJDD(2) - pkin(2) * t762 + qJ(3) * t769 + qJD(2) * t730 + t794 * t607 + t751 * t608 - t728 * t781;
t603 = -mrSges(2,2) * g(3) - mrSges(2,3) * t743 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t759 - pkin(7) * t612 - t604 * t753 + t605 * t756;
t1 = [-m(1) * g(1) + t771; -m(1) * g(2) + t791; (-m(1) - m(2)) * g(3) + t612; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t791 + t757 * t603 - t754 * t606; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t771 + t754 * t603 + t757 * t606; -mrSges(1,1) * g(2) + mrSges(2,1) * t743 + mrSges(1,2) * g(1) - mrSges(2,2) * t744 + Ifges(2,3) * qJDD(1) + pkin(1) * t760 + pkin(7) * t770 + t756 * t604 + t753 * t605;];
tauB  = t1;
