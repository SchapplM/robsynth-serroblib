% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRRR7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR7_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR7_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR7_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR7_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR7_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:22:01
% EndTime: 2019-12-31 22:22:10
% DurationCPUTime: 8.76s
% Computational Cost: add. (132098->316), mult. (286748->399), div. (0->0), fcn. (206588->10), ass. (0->128)
t755 = sin(qJ(2));
t760 = cos(qJ(2));
t779 = qJD(1) * qJD(2);
t734 = t755 * qJDD(1) + t760 * t779;
t756 = sin(qJ(1));
t761 = cos(qJ(1));
t741 = -t761 * g(1) - t756 * g(2);
t762 = qJD(1) ^ 2;
t729 = -t762 * pkin(1) + qJDD(1) * pkin(6) + t741;
t783 = t755 * t729;
t784 = pkin(2) * t762;
t695 = qJDD(2) * pkin(2) - t734 * pkin(7) - t783 + (pkin(7) * t779 + t755 * t784 - g(3)) * t760;
t717 = -t755 * g(3) + t760 * t729;
t735 = t760 * qJDD(1) - t755 * t779;
t781 = qJD(1) * t755;
t739 = qJD(2) * pkin(2) - pkin(7) * t781;
t751 = t760 ^ 2;
t696 = t735 * pkin(7) - qJD(2) * t739 - t751 * t784 + t717;
t754 = sin(qJ(3));
t759 = cos(qJ(3));
t679 = t759 * t695 - t754 * t696;
t726 = (-t754 * t755 + t759 * t760) * qJD(1);
t703 = t726 * qJD(3) + t759 * t734 + t754 * t735;
t727 = (t754 * t760 + t755 * t759) * qJD(1);
t748 = qJDD(2) + qJDD(3);
t749 = qJD(2) + qJD(3);
t658 = (t726 * t749 - t703) * pkin(8) + (t726 * t727 + t748) * pkin(3) + t679;
t680 = t754 * t695 + t759 * t696;
t702 = -t727 * qJD(3) - t754 * t734 + t759 * t735;
t720 = t749 * pkin(3) - t727 * pkin(8);
t722 = t726 ^ 2;
t661 = -t722 * pkin(3) + t702 * pkin(8) - t749 * t720 + t680;
t753 = sin(qJ(4));
t758 = cos(qJ(4));
t656 = t753 * t658 + t758 * t661;
t714 = t753 * t726 + t758 * t727;
t675 = -t714 * qJD(4) + t758 * t702 - t753 * t703;
t713 = t758 * t726 - t753 * t727;
t689 = -t713 * mrSges(5,1) + t714 * mrSges(5,2);
t746 = qJD(4) + t749;
t706 = t746 * mrSges(5,1) - t714 * mrSges(5,3);
t745 = qJDD(4) + t748;
t690 = -t713 * pkin(4) - t714 * pkin(9);
t744 = t746 ^ 2;
t652 = -t744 * pkin(4) + t745 * pkin(9) + t713 * t690 + t656;
t740 = t756 * g(1) - t761 * g(2);
t772 = -qJDD(1) * pkin(1) - t740;
t704 = -t735 * pkin(2) + t739 * t781 + (-pkin(7) * t751 - pkin(6)) * t762 + t772;
t665 = -t702 * pkin(3) - t722 * pkin(8) + t727 * t720 + t704;
t676 = t713 * qJD(4) + t753 * t702 + t758 * t703;
t653 = (-t713 * t746 - t676) * pkin(9) + (t714 * t746 - t675) * pkin(4) + t665;
t752 = sin(qJ(5));
t757 = cos(qJ(5));
t649 = -t752 * t652 + t757 * t653;
t697 = -t752 * t714 + t757 * t746;
t663 = t697 * qJD(5) + t757 * t676 + t752 * t745;
t673 = qJDD(5) - t675;
t698 = t757 * t714 + t752 * t746;
t681 = -t697 * mrSges(6,1) + t698 * mrSges(6,2);
t710 = qJD(5) - t713;
t682 = -t710 * mrSges(6,2) + t697 * mrSges(6,3);
t645 = m(6) * t649 + t673 * mrSges(6,1) - t663 * mrSges(6,3) - t698 * t681 + t710 * t682;
t650 = t757 * t652 + t752 * t653;
t662 = -t698 * qJD(5) - t752 * t676 + t757 * t745;
t683 = t710 * mrSges(6,1) - t698 * mrSges(6,3);
t646 = m(6) * t650 - t673 * mrSges(6,2) + t662 * mrSges(6,3) + t697 * t681 - t710 * t683;
t774 = -t752 * t645 + t757 * t646;
t632 = m(5) * t656 - t745 * mrSges(5,2) + t675 * mrSges(5,3) + t713 * t689 - t746 * t706 + t774;
t655 = t758 * t658 - t753 * t661;
t705 = -t746 * mrSges(5,2) + t713 * mrSges(5,3);
t651 = -t745 * pkin(4) - t744 * pkin(9) + t714 * t690 - t655;
t769 = -m(6) * t651 + t662 * mrSges(6,1) - t663 * mrSges(6,2) + t697 * t682 - t698 * t683;
t641 = m(5) * t655 + t745 * mrSges(5,1) - t676 * mrSges(5,3) - t714 * t689 + t746 * t705 + t769;
t626 = t753 * t632 + t758 * t641;
t715 = -t726 * mrSges(4,1) + t727 * mrSges(4,2);
t718 = -t749 * mrSges(4,2) + t726 * mrSges(4,3);
t623 = m(4) * t679 + t748 * mrSges(4,1) - t703 * mrSges(4,3) - t727 * t715 + t749 * t718 + t626;
t719 = t749 * mrSges(4,1) - t727 * mrSges(4,3);
t775 = t758 * t632 - t753 * t641;
t624 = m(4) * t680 - t748 * mrSges(4,2) + t702 * mrSges(4,3) + t726 * t715 - t749 * t719 + t775;
t617 = t759 * t623 + t754 * t624;
t716 = -t760 * g(3) - t783;
t724 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t755 + Ifges(3,2) * t760) * qJD(1);
t725 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t755 + Ifges(3,4) * t760) * qJD(1);
t708 = Ifges(4,4) * t727 + Ifges(4,2) * t726 + Ifges(4,6) * t749;
t709 = Ifges(4,1) * t727 + Ifges(4,4) * t726 + Ifges(4,5) * t749;
t666 = Ifges(6,5) * t698 + Ifges(6,6) * t697 + Ifges(6,3) * t710;
t668 = Ifges(6,1) * t698 + Ifges(6,4) * t697 + Ifges(6,5) * t710;
t638 = -mrSges(6,1) * t651 + mrSges(6,3) * t650 + Ifges(6,4) * t663 + Ifges(6,2) * t662 + Ifges(6,6) * t673 - t698 * t666 + t710 * t668;
t667 = Ifges(6,4) * t698 + Ifges(6,2) * t697 + Ifges(6,6) * t710;
t639 = mrSges(6,2) * t651 - mrSges(6,3) * t649 + Ifges(6,1) * t663 + Ifges(6,4) * t662 + Ifges(6,5) * t673 + t697 * t666 - t710 * t667;
t685 = Ifges(5,4) * t714 + Ifges(5,2) * t713 + Ifges(5,6) * t746;
t686 = Ifges(5,1) * t714 + Ifges(5,4) * t713 + Ifges(5,5) * t746;
t768 = -mrSges(5,1) * t655 + mrSges(5,2) * t656 - Ifges(5,5) * t676 - Ifges(5,6) * t675 - Ifges(5,3) * t745 - pkin(4) * t769 - pkin(9) * t774 - t757 * t638 - t752 * t639 - t714 * t685 + t713 * t686;
t765 = -mrSges(4,1) * t679 + mrSges(4,2) * t680 - Ifges(4,5) * t703 - Ifges(4,6) * t702 - Ifges(4,3) * t748 - pkin(3) * t626 - t727 * t708 + t726 * t709 + t768;
t785 = mrSges(3,1) * t716 - mrSges(3,2) * t717 + Ifges(3,5) * t734 + Ifges(3,6) * t735 + Ifges(3,3) * qJDD(2) + pkin(2) * t617 + (t755 * t724 - t760 * t725) * qJD(1) - t765;
t733 = (-mrSges(3,1) * t760 + mrSges(3,2) * t755) * qJD(1);
t780 = qJD(1) * t760;
t738 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t780;
t615 = m(3) * t716 + qJDD(2) * mrSges(3,1) - t734 * mrSges(3,3) + qJD(2) * t738 - t733 * t781 + t617;
t737 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t781;
t776 = -t754 * t623 + t759 * t624;
t616 = m(3) * t717 - qJDD(2) * mrSges(3,2) + t735 * mrSges(3,3) - qJD(2) * t737 + t733 * t780 + t776;
t777 = -t755 * t615 + t760 * t616;
t608 = m(2) * t741 - t762 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t777;
t728 = -t762 * pkin(6) + t772;
t634 = t757 * t645 + t752 * t646;
t771 = m(5) * t665 - t675 * mrSges(5,1) + t676 * mrSges(5,2) - t713 * t705 + t714 * t706 + t634;
t767 = m(4) * t704 - t702 * mrSges(4,1) + t703 * mrSges(4,2) - t726 * t718 + t727 * t719 + t771;
t764 = -m(3) * t728 + t735 * mrSges(3,1) - t734 * mrSges(3,2) - t737 * t781 + t738 * t780 - t767;
t628 = m(2) * t740 + qJDD(1) * mrSges(2,1) - t762 * mrSges(2,2) + t764;
t782 = t756 * t608 + t761 * t628;
t610 = t760 * t615 + t755 * t616;
t778 = t761 * t608 - t756 * t628;
t684 = Ifges(5,5) * t714 + Ifges(5,6) * t713 + Ifges(5,3) * t746;
t618 = mrSges(5,2) * t665 - mrSges(5,3) * t655 + Ifges(5,1) * t676 + Ifges(5,4) * t675 + Ifges(5,5) * t745 - pkin(9) * t634 - t752 * t638 + t757 * t639 + t713 * t684 - t746 * t685;
t766 = mrSges(6,1) * t649 - mrSges(6,2) * t650 + Ifges(6,5) * t663 + Ifges(6,6) * t662 + Ifges(6,3) * t673 + t698 * t667 - t697 * t668;
t619 = -mrSges(5,1) * t665 + mrSges(5,3) * t656 + Ifges(5,4) * t676 + Ifges(5,2) * t675 + Ifges(5,6) * t745 - pkin(4) * t634 - t714 * t684 + t746 * t686 - t766;
t707 = Ifges(4,5) * t727 + Ifges(4,6) * t726 + Ifges(4,3) * t749;
t605 = -mrSges(4,1) * t704 + mrSges(4,3) * t680 + Ifges(4,4) * t703 + Ifges(4,2) * t702 + Ifges(4,6) * t748 - pkin(3) * t771 + pkin(8) * t775 + t753 * t618 + t758 * t619 - t727 * t707 + t749 * t709;
t611 = mrSges(4,2) * t704 - mrSges(4,3) * t679 + Ifges(4,1) * t703 + Ifges(4,4) * t702 + Ifges(4,5) * t748 - pkin(8) * t626 + t758 * t618 - t753 * t619 + t726 * t707 - t749 * t708;
t723 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t755 + Ifges(3,6) * t760) * qJD(1);
t601 = -mrSges(3,1) * t728 + mrSges(3,3) * t717 + Ifges(3,4) * t734 + Ifges(3,2) * t735 + Ifges(3,6) * qJDD(2) - pkin(2) * t767 + pkin(7) * t776 + qJD(2) * t725 + t759 * t605 + t754 * t611 - t723 * t781;
t603 = mrSges(3,2) * t728 - mrSges(3,3) * t716 + Ifges(3,1) * t734 + Ifges(3,4) * t735 + Ifges(3,5) * qJDD(2) - pkin(7) * t617 - qJD(2) * t724 - t754 * t605 + t759 * t611 + t723 * t780;
t770 = mrSges(2,1) * t740 - mrSges(2,2) * t741 + Ifges(2,3) * qJDD(1) + pkin(1) * t764 + pkin(6) * t777 + t760 * t601 + t755 * t603;
t604 = mrSges(2,1) * g(3) + mrSges(2,3) * t741 + t762 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t610 - t785;
t599 = -mrSges(2,2) * g(3) - mrSges(2,3) * t740 + Ifges(2,5) * qJDD(1) - t762 * Ifges(2,6) - pkin(6) * t610 - t755 * t601 + t760 * t603;
t1 = [-m(1) * g(1) + t778; -m(1) * g(2) + t782; (-m(1) - m(2)) * g(3) + t610; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t782 + t761 * t599 - t756 * t604; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t778 + t756 * t599 + t761 * t604; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t770; t770; t785; -t765; -t768; t766;];
tauJB = t1;
