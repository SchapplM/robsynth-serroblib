% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PPRRRR2
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% tauJB [(6+6)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 20:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PPRRRR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRR2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_invdynJB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:55:46
% EndTime: 2019-05-04 20:56:08
% DurationCPUTime: 22.45s
% Computational Cost: add. (421478->295), mult. (768276->386), div. (0->0), fcn. (602083->16), ass. (0->138)
t762 = sin(pkin(12));
t766 = cos(pkin(12));
t753 = -g(1) * t766 - g(2) * t762;
t761 = sin(pkin(13));
t765 = cos(pkin(13));
t752 = g(1) * t762 - g(2) * t766;
t760 = -g(3) + qJDD(1);
t764 = sin(pkin(6));
t768 = cos(pkin(6));
t786 = t752 * t768 + t760 * t764;
t719 = -t761 * t753 + t765 * t786;
t720 = t765 * t753 + t761 * t786;
t735 = -t752 * t764 + t760 * t768 + qJDD(2);
t776 = cos(qJ(3));
t767 = cos(pkin(7));
t772 = sin(qJ(3));
t798 = t767 * t772;
t763 = sin(pkin(7));
t799 = t763 * t772;
t700 = t719 * t798 + t776 * t720 + t735 * t799;
t778 = qJD(3) ^ 2;
t698 = -pkin(3) * t778 + qJDD(3) * pkin(9) + t700;
t712 = -t719 * t763 + t735 * t767;
t771 = sin(qJ(4));
t775 = cos(qJ(4));
t691 = t775 * t698 + t771 * t712;
t749 = (-pkin(4) * t775 - pkin(10) * t771) * qJD(3);
t777 = qJD(4) ^ 2;
t795 = qJD(3) * t775;
t689 = -pkin(4) * t777 + qJDD(4) * pkin(10) + t749 * t795 + t691;
t699 = -t772 * t720 + (t719 * t767 + t735 * t763) * t776;
t697 = -qJDD(3) * pkin(3) - t778 * pkin(9) - t699;
t794 = qJD(3) * qJD(4);
t793 = t775 * t794;
t750 = qJDD(3) * t771 + t793;
t758 = t771 * t794;
t751 = qJDD(3) * t775 - t758;
t694 = (-t750 - t793) * pkin(10) + (-t751 + t758) * pkin(4) + t697;
t770 = sin(qJ(5));
t774 = cos(qJ(5));
t684 = -t770 * t689 + t774 * t694;
t796 = qJD(3) * t771;
t746 = qJD(4) * t774 - t770 * t796;
t727 = qJD(5) * t746 + qJDD(4) * t770 + t750 * t774;
t745 = qJDD(5) - t751;
t747 = qJD(4) * t770 + t774 * t796;
t757 = qJD(5) - t795;
t682 = (t746 * t757 - t727) * pkin(11) + (t746 * t747 + t745) * pkin(5) + t684;
t685 = t774 * t689 + t770 * t694;
t726 = -qJD(5) * t747 + qJDD(4) * t774 - t750 * t770;
t734 = pkin(5) * t757 - pkin(11) * t747;
t744 = t746 ^ 2;
t683 = -pkin(5) * t744 + pkin(11) * t726 - t734 * t757 + t685;
t769 = sin(qJ(6));
t773 = cos(qJ(6));
t680 = t682 * t773 - t683 * t769;
t728 = t746 * t773 - t747 * t769;
t706 = qJD(6) * t728 + t726 * t769 + t727 * t773;
t729 = t746 * t769 + t747 * t773;
t713 = -mrSges(7,1) * t728 + mrSges(7,2) * t729;
t756 = qJD(6) + t757;
t715 = -mrSges(7,2) * t756 + mrSges(7,3) * t728;
t741 = qJDD(6) + t745;
t676 = m(7) * t680 + mrSges(7,1) * t741 - mrSges(7,3) * t706 - t713 * t729 + t715 * t756;
t681 = t682 * t769 + t683 * t773;
t705 = -qJD(6) * t729 + t726 * t773 - t727 * t769;
t716 = mrSges(7,1) * t756 - mrSges(7,3) * t729;
t677 = m(7) * t681 - mrSges(7,2) * t741 + mrSges(7,3) * t705 + t713 * t728 - t716 * t756;
t668 = t773 * t676 + t769 * t677;
t730 = -mrSges(6,1) * t746 + mrSges(6,2) * t747;
t732 = -mrSges(6,2) * t757 + mrSges(6,3) * t746;
t666 = m(6) * t684 + mrSges(6,1) * t745 - mrSges(6,3) * t727 - t730 * t747 + t732 * t757 + t668;
t733 = mrSges(6,1) * t757 - mrSges(6,3) * t747;
t790 = -t676 * t769 + t773 * t677;
t667 = m(6) * t685 - mrSges(6,2) * t745 + mrSges(6,3) * t726 + t730 * t746 - t733 * t757 + t790;
t664 = -t666 * t770 + t774 * t667;
t748 = (-mrSges(5,1) * t775 + mrSges(5,2) * t771) * qJD(3);
t754 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t796;
t662 = m(5) * t691 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t751 - qJD(4) * t754 + t748 * t795 + t664;
t690 = -t771 * t698 + t712 * t775;
t688 = -qJDD(4) * pkin(4) - pkin(10) * t777 + t749 * t796 - t690;
t686 = -pkin(5) * t726 - pkin(11) * t744 + t734 * t747 + t688;
t783 = m(7) * t686 - t705 * mrSges(7,1) + mrSges(7,2) * t706 - t728 * t715 + t716 * t729;
t678 = -m(6) * t688 + t726 * mrSges(6,1) - mrSges(6,2) * t727 + t746 * t732 - t733 * t747 - t783;
t755 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t795;
t672 = m(5) * t690 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t750 + qJD(4) * t755 - t748 * t796 + t678;
t791 = t775 * t662 - t672 * t771;
t651 = m(4) * t700 - mrSges(4,1) * t778 - qJDD(3) * mrSges(4,2) + t791;
t654 = t771 * t662 + t775 * t672;
t653 = m(4) * t712 + t654;
t663 = t666 * t774 + t667 * t770;
t781 = -m(5) * t697 + t751 * mrSges(5,1) - mrSges(5,2) * t750 - t754 * t796 + t755 * t795 - t663;
t659 = m(4) * t699 + qJDD(3) * mrSges(4,1) - mrSges(4,2) * t778 + t781;
t800 = t659 * t776;
t640 = t651 * t798 - t653 * t763 + t767 * t800;
t636 = m(3) * t719 + t640;
t646 = t776 * t651 - t659 * t772;
t645 = m(3) * t720 + t646;
t804 = t636 * t765 + t645 * t761;
t707 = Ifges(7,5) * t729 + Ifges(7,6) * t728 + Ifges(7,3) * t756;
t709 = Ifges(7,1) * t729 + Ifges(7,4) * t728 + Ifges(7,5) * t756;
t669 = -mrSges(7,1) * t686 + mrSges(7,3) * t681 + Ifges(7,4) * t706 + Ifges(7,2) * t705 + Ifges(7,6) * t741 - t707 * t729 + t709 * t756;
t708 = Ifges(7,4) * t729 + Ifges(7,2) * t728 + Ifges(7,6) * t756;
t670 = mrSges(7,2) * t686 - mrSges(7,3) * t680 + Ifges(7,1) * t706 + Ifges(7,4) * t705 + Ifges(7,5) * t741 + t707 * t728 - t708 * t756;
t721 = Ifges(6,5) * t747 + Ifges(6,6) * t746 + Ifges(6,3) * t757;
t723 = Ifges(6,1) * t747 + Ifges(6,4) * t746 + Ifges(6,5) * t757;
t655 = -mrSges(6,1) * t688 + mrSges(6,3) * t685 + Ifges(6,4) * t727 + Ifges(6,2) * t726 + Ifges(6,6) * t745 - pkin(5) * t783 + pkin(11) * t790 + t773 * t669 + t769 * t670 - t747 * t721 + t757 * t723;
t722 = Ifges(6,4) * t747 + Ifges(6,2) * t746 + Ifges(6,6) * t757;
t656 = mrSges(6,2) * t688 - mrSges(6,3) * t684 + Ifges(6,1) * t727 + Ifges(6,4) * t726 + Ifges(6,5) * t745 - pkin(11) * t668 - t669 * t769 + t670 * t773 + t721 * t746 - t722 * t757;
t739 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t771 + Ifges(5,2) * t775) * qJD(3);
t740 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t771 + Ifges(5,4) * t775) * qJD(3);
t803 = mrSges(5,1) * t690 - mrSges(5,2) * t691 + Ifges(5,5) * t750 + Ifges(5,6) * t751 + Ifges(5,3) * qJDD(4) + pkin(4) * t678 + pkin(10) * t664 + t774 * t655 + t770 * t656 + (t739 * t771 - t740 * t775) * qJD(3);
t639 = t651 * t799 + t767 * t653 + t763 * t800;
t638 = m(3) * t735 + t639;
t626 = -t638 * t764 + t768 * t804;
t624 = m(2) * t752 + t626;
t632 = -t636 * t761 + t765 * t645;
t631 = m(2) * t753 + t632;
t797 = t766 * t624 + t762 * t631;
t625 = t768 * t638 + t764 * t804;
t792 = -t624 * t762 + t766 * t631;
t789 = m(2) * t760 + t625;
t738 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t771 + Ifges(5,6) * t775) * qJD(3);
t641 = mrSges(5,2) * t697 - mrSges(5,3) * t690 + Ifges(5,1) * t750 + Ifges(5,4) * t751 + Ifges(5,5) * qJDD(4) - pkin(10) * t663 - qJD(4) * t739 - t655 * t770 + t656 * t774 + t738 * t795;
t782 = -mrSges(7,1) * t680 + mrSges(7,2) * t681 - Ifges(7,5) * t706 - Ifges(7,6) * t705 - Ifges(7,3) * t741 - t729 * t708 + t728 * t709;
t779 = mrSges(6,1) * t684 - mrSges(6,2) * t685 + Ifges(6,5) * t727 + Ifges(6,6) * t726 + Ifges(6,3) * t745 + pkin(5) * t668 + t747 * t722 - t746 * t723 - t782;
t647 = -mrSges(5,1) * t697 + mrSges(5,3) * t691 + Ifges(5,4) * t750 + Ifges(5,2) * t751 + Ifges(5,6) * qJDD(4) - pkin(4) * t663 + qJD(4) * t740 - t738 * t796 - t779;
t628 = mrSges(4,2) * t712 - mrSges(4,3) * t699 + Ifges(4,5) * qJDD(3) - Ifges(4,6) * t778 - pkin(9) * t654 + t641 * t775 - t647 * t771;
t633 = -mrSges(4,1) * t712 + mrSges(4,3) * t700 + t778 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t654 - t803;
t785 = pkin(8) * t646 + t628 * t772 + t633 * t776;
t627 = mrSges(4,1) * t699 - mrSges(4,2) * t700 + Ifges(4,3) * qJDD(3) + pkin(3) * t781 + pkin(9) * t791 + t771 * t641 + t775 * t647;
t621 = -mrSges(3,1) * t735 + mrSges(3,3) * t720 - pkin(2) * t639 - t763 * t627 + t767 * t785;
t622 = mrSges(3,2) * t735 - mrSges(3,3) * t719 + t776 * t628 - t772 * t633 + (-t639 * t763 - t640 * t767) * pkin(8);
t784 = qJ(2) * t632 + t621 * t765 + t622 * t761;
t620 = mrSges(3,1) * t719 - mrSges(3,2) * t720 + pkin(2) * t640 + t767 * t627 + t763 * t785;
t619 = mrSges(2,2) * t760 - mrSges(2,3) * t752 - t761 * t621 + t765 * t622 + (-t625 * t764 - t626 * t768) * qJ(2);
t618 = -mrSges(2,1) * t760 + mrSges(2,3) * t753 - pkin(1) * t625 - t764 * t620 + t768 * t784;
t1 = [-m(1) * g(1) + t792; -m(1) * g(2) + t797; -m(1) * g(3) + t789; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t797 - t762 * t618 + t766 * t619; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t792 + t766 * t618 + t762 * t619; -mrSges(1,1) * g(2) + mrSges(2,1) * t752 + mrSges(1,2) * g(1) - mrSges(2,2) * t753 + pkin(1) * t626 + t768 * t620 + t764 * t784; t789; t638; t627; t803; t779; -t782;];
tauJB  = t1;
