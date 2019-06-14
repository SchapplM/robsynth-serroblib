% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
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
% Datum: 2019-05-04 23:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRPRPR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR6_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR6_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR6_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR6_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR6_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:08:20
% EndTime: 2019-05-04 23:08:26
% DurationCPUTime: 6.34s
% Computational Cost: add. (104855->291), mult. (207749->366), div. (0->0), fcn. (136787->12), ass. (0->128)
t766 = sin(pkin(10));
t769 = cos(pkin(10));
t751 = t766 * g(1) - t769 * g(2);
t752 = -t769 * g(1) - t766 * g(2);
t762 = -g(3) + qJDD(1);
t776 = cos(qJ(2));
t770 = cos(pkin(6));
t773 = sin(qJ(2));
t801 = t770 * t773;
t767 = sin(pkin(6));
t802 = t767 * t773;
t709 = t751 * t801 + t776 * t752 + t762 * t802;
t808 = -qJDD(2) * qJ(3) - (2 * qJD(3) * qJD(2)) - t709;
t708 = -t773 * t752 + (t751 * t770 + t762 * t767) * t776;
t807 = -pkin(2) - pkin(8);
t806 = mrSges(3,1) - mrSges(4,2);
t805 = (-Ifges(4,4) + Ifges(3,5));
t804 = Ifges(4,5) - Ifges(3,6);
t778 = qJD(2) ^ 2;
t782 = -t778 * qJ(3) + qJDD(3) - t708;
t703 = t807 * qJDD(2) + t782;
t726 = -t767 * t751 + t770 * t762;
t772 = sin(qJ(4));
t775 = cos(qJ(4));
t698 = t772 * t703 + t775 * t726;
t748 = (mrSges(5,1) * t772 + mrSges(5,2) * t775) * qJD(2);
t797 = qJD(2) * qJD(4);
t794 = t775 * t797;
t749 = t772 * qJDD(2) + t794;
t799 = qJD(2) * t775;
t754 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t799;
t747 = (pkin(4) * t772 - qJ(5) * t775) * qJD(2);
t777 = qJD(4) ^ 2;
t798 = t772 * qJD(2);
t685 = -t777 * pkin(4) + qJDD(4) * qJ(5) - t747 * t798 + t698;
t702 = t807 * t778 - t808;
t795 = t772 * t797;
t750 = t775 * qJDD(2) - t795;
t688 = (-t750 + t795) * qJ(5) + (t749 + t794) * pkin(4) + t702;
t765 = sin(pkin(11));
t768 = cos(pkin(11));
t740 = t765 * qJD(4) + t768 * t799;
t680 = -0.2e1 * qJD(5) * t740 - t765 * t685 + t768 * t688;
t724 = t765 * qJDD(4) + t768 * t750;
t739 = t768 * qJD(4) - t765 * t799;
t678 = (t739 * t798 - t724) * pkin(9) + (t739 * t740 + t749) * pkin(5) + t680;
t681 = 0.2e1 * qJD(5) * t739 + t768 * t685 + t765 * t688;
t723 = t768 * qJDD(4) - t765 * t750;
t725 = pkin(5) * t798 - t740 * pkin(9);
t738 = t739 ^ 2;
t679 = -t738 * pkin(5) + t723 * pkin(9) - t725 * t798 + t681;
t771 = sin(qJ(6));
t774 = cos(qJ(6));
t676 = t774 * t678 - t771 * t679;
t714 = t774 * t739 - t771 * t740;
t691 = t714 * qJD(6) + t771 * t723 + t774 * t724;
t715 = t771 * t739 + t774 * t740;
t699 = -t714 * mrSges(7,1) + t715 * mrSges(7,2);
t756 = qJD(6) + t798;
t706 = -t756 * mrSges(7,2) + t714 * mrSges(7,3);
t744 = qJDD(6) + t749;
t672 = m(7) * t676 + t744 * mrSges(7,1) - t691 * mrSges(7,3) - t715 * t699 + t756 * t706;
t677 = t771 * t678 + t774 * t679;
t690 = -t715 * qJD(6) + t774 * t723 - t771 * t724;
t707 = t756 * mrSges(7,1) - t715 * mrSges(7,3);
t673 = m(7) * t677 - t744 * mrSges(7,2) + t690 * mrSges(7,3) + t714 * t699 - t756 * t707;
t665 = t774 * t672 + t771 * t673;
t716 = -t739 * mrSges(6,1) + t740 * mrSges(6,2);
t721 = -mrSges(6,2) * t798 + t739 * mrSges(6,3);
t663 = m(6) * t680 + t749 * mrSges(6,1) - t724 * mrSges(6,3) - t740 * t716 + t721 * t798 + t665;
t722 = mrSges(6,1) * t798 - t740 * mrSges(6,3);
t790 = -t771 * t672 + t774 * t673;
t664 = m(6) * t681 - t749 * mrSges(6,2) + t723 * mrSges(6,3) + t739 * t716 - t722 * t798 + t790;
t791 = -t765 * t663 + t768 * t664;
t657 = m(5) * t698 - qJDD(4) * mrSges(5,2) - t749 * mrSges(5,3) - qJD(4) * t754 - t748 * t798 + t791;
t697 = t775 * t703 - t772 * t726;
t684 = -qJDD(4) * pkin(4) - t777 * qJ(5) + t747 * t799 + qJDD(5) - t697;
t682 = -t723 * pkin(5) - t738 * pkin(9) + t740 * t725 + t684;
t783 = m(7) * t682 - t690 * mrSges(7,1) + t691 * mrSges(7,2) - t714 * t706 + t715 * t707;
t675 = m(6) * t684 - t723 * mrSges(6,1) + t724 * mrSges(6,2) - t739 * t721 + t740 * t722 + t783;
t753 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t798;
t668 = m(5) * t697 + qJDD(4) * mrSges(5,1) - t750 * mrSges(5,3) + qJD(4) * t753 - t748 * t799 - t675;
t647 = t772 * t657 + t775 * t668;
t705 = -qJDD(2) * pkin(2) + t782;
t785 = -m(4) * t705 + (t778 * mrSges(4,3)) - t647;
t642 = m(3) * t708 - (t778 * mrSges(3,2)) + t806 * qJDD(2) + t785;
t803 = t642 * t776;
t792 = t775 * t657 - t772 * t668;
t646 = m(4) * t726 + t792;
t645 = m(3) * t726 + t646;
t704 = t778 * pkin(2) + t808;
t659 = t768 * t663 + t765 * t664;
t784 = -m(5) * t702 - t749 * mrSges(5,1) - t750 * mrSges(5,2) - t753 * t798 - t754 * t799 - t659;
t780 = -m(4) * t704 + (t778 * mrSges(4,2)) + qJDD(2) * mrSges(4,3) - t784;
t655 = m(3) * t709 - (t778 * mrSges(3,1)) - qJDD(2) * mrSges(3,2) + t780;
t633 = -t767 * t645 + t655 * t801 + t770 * t803;
t631 = m(2) * t751 + t633;
t638 = -t773 * t642 + t776 * t655;
t637 = m(2) * t752 + t638;
t800 = t769 * t631 + t766 * t637;
t632 = t770 * t645 + t655 * t802 + t767 * t803;
t793 = -t766 * t631 + t769 * t637;
t788 = m(2) * t762 + t632;
t692 = Ifges(7,5) * t715 + Ifges(7,6) * t714 + Ifges(7,3) * t756;
t694 = Ifges(7,1) * t715 + Ifges(7,4) * t714 + Ifges(7,5) * t756;
t666 = -mrSges(7,1) * t682 + mrSges(7,3) * t677 + Ifges(7,4) * t691 + Ifges(7,2) * t690 + Ifges(7,6) * t744 - t715 * t692 + t756 * t694;
t693 = Ifges(7,4) * t715 + Ifges(7,2) * t714 + Ifges(7,6) * t756;
t667 = mrSges(7,2) * t682 - mrSges(7,3) * t676 + Ifges(7,1) * t691 + Ifges(7,4) * t690 + Ifges(7,5) * t744 + t714 * t692 - t756 * t693;
t710 = Ifges(6,5) * t740 + Ifges(6,6) * t739 + Ifges(6,3) * t798;
t712 = Ifges(6,1) * t740 + Ifges(6,4) * t739 + Ifges(6,5) * t798;
t649 = -mrSges(6,1) * t684 + mrSges(6,3) * t681 + Ifges(6,4) * t724 + Ifges(6,2) * t723 + Ifges(6,6) * t749 - pkin(5) * t783 + pkin(9) * t790 + t774 * t666 + t771 * t667 - t740 * t710 + t712 * t798;
t711 = Ifges(6,4) * t740 + Ifges(6,2) * t739 + Ifges(6,6) * t798;
t651 = mrSges(6,2) * t684 - mrSges(6,3) * t680 + Ifges(6,1) * t724 + Ifges(6,4) * t723 + Ifges(6,5) * t749 - pkin(9) * t665 - t771 * t666 + t774 * t667 + t739 * t710 - t711 * t798;
t733 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t775 - Ifges(5,6) * t772) * qJD(2);
t734 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t775 - Ifges(5,2) * t772) * qJD(2);
t634 = mrSges(5,2) * t702 - mrSges(5,3) * t697 + Ifges(5,1) * t750 - Ifges(5,4) * t749 + Ifges(5,5) * qJDD(4) - qJ(5) * t659 - qJD(4) * t734 - t765 * t649 + t768 * t651 - t733 * t798;
t735 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t775 - Ifges(5,4) * t772) * qJD(2);
t779 = mrSges(7,1) * t676 - mrSges(7,2) * t677 + Ifges(7,5) * t691 + Ifges(7,6) * t690 + Ifges(7,3) * t744 + t715 * t693 - t714 * t694;
t639 = (-Ifges(5,2) - Ifges(6,3)) * t749 - t779 + Ifges(5,4) * t750 + t739 * t712 - t740 * t711 - Ifges(6,6) * t723 - Ifges(6,5) * t724 + qJD(4) * t735 + mrSges(5,3) * t698 - mrSges(5,1) * t702 - mrSges(6,1) * t680 + mrSges(6,2) * t681 - pkin(5) * t665 + Ifges(5,6) * qJDD(4) - pkin(4) * t659 - t733 * t799;
t628 = -mrSges(4,1) * t704 + mrSges(3,3) * t709 - pkin(2) * t646 - pkin(3) * t784 - pkin(8) * t792 - t804 * qJDD(2) - t772 * t634 - t775 * t639 - t806 * t726 + (t805 * t778);
t781 = mrSges(5,1) * t697 - mrSges(5,2) * t698 + Ifges(5,5) * t750 - Ifges(5,6) * t749 + Ifges(5,3) * qJDD(4) - pkin(4) * t675 + qJ(5) * t791 + t768 * t649 + t765 * t651 + t734 * t799 + t735 * t798;
t629 = t781 + mrSges(4,1) * t705 - mrSges(3,3) * t708 + pkin(3) * t647 - qJ(3) * t646 + t805 * qJDD(2) + (mrSges(3,2) - mrSges(4,3)) * t726 + t804 * t778;
t786 = pkin(7) * t638 + t628 * t776 + t629 * t773;
t643 = qJDD(2) * mrSges(4,2) - t785;
t627 = mrSges(3,1) * t708 - mrSges(3,2) * t709 + mrSges(4,2) * t705 - mrSges(4,3) * t704 + t775 * t634 - t772 * t639 - pkin(8) * t647 - pkin(2) * t643 + qJ(3) * t780 + (Ifges(3,3) + Ifges(4,1)) * qJDD(2);
t626 = mrSges(2,2) * t762 - mrSges(2,3) * t751 - t773 * t628 + t776 * t629 + (-t632 * t767 - t633 * t770) * pkin(7);
t625 = -mrSges(2,1) * t762 + mrSges(2,3) * t752 - pkin(1) * t632 - t767 * t627 + t786 * t770;
t1 = [-m(1) * g(1) + t793; -m(1) * g(2) + t800; -m(1) * g(3) + t788; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t800 - t766 * t625 + t769 * t626; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t793 + t769 * t625 + t766 * t626; -mrSges(1,1) * g(2) + mrSges(2,1) * t751 + mrSges(1,2) * g(1) - mrSges(2,2) * t752 + pkin(1) * t633 + t770 * t627 + t786 * t767; t788; t627; t643; t781; t675; t779;];
tauJB  = t1;
