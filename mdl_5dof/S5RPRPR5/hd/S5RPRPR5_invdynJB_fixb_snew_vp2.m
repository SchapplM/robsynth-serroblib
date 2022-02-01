% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m [6x1]
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
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRPR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR5_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR5_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR5_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:25:27
% EndTime: 2022-01-23 09:25:35
% DurationCPUTime: 7.98s
% Computational Cost: add. (80930->290), mult. (207955->381), div. (0->0), fcn. (141748->10), ass. (0->126)
t755 = sin(qJ(1));
t758 = cos(qJ(1));
t735 = -t758 * g(1) - t755 * g(2);
t759 = qJD(1) ^ 2;
t795 = -t759 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t735;
t750 = sin(pkin(8));
t752 = cos(pkin(8));
t704 = -t752 * g(3) - t795 * t750;
t783 = t752 * qJD(1);
t738 = qJD(3) - t783;
t754 = sin(qJ(3));
t784 = t750 * qJD(1);
t778 = t754 * t784;
t719 = -t738 * mrSges(4,2) - mrSges(4,3) * t778;
t757 = cos(qJ(3));
t777 = t757 * t784;
t721 = t738 * mrSges(4,1) - mrSges(4,3) * t777;
t794 = -t719 * t754 - t721 * t757;
t705 = -t750 * g(3) + t795 * t752;
t770 = -pkin(2) * t752 - pkin(6) * t750;
t730 = t770 * qJD(1);
t694 = t730 * t783 + t705;
t734 = t755 * g(1) - t758 * g(2);
t765 = -t759 * qJ(2) + qJDD(2) - t734;
t706 = (-pkin(1) + t770) * qJDD(1) + t765;
t703 = t757 * t706;
t781 = qJD(1) * qJD(3);
t724 = (qJDD(1) * t757 - t754 * t781) * t750;
t780 = t752 * qJDD(1);
t737 = qJDD(3) - t780;
t788 = t750 ^ 2 * t759;
t673 = t737 * pkin(3) - t724 * qJ(4) + t703 + (-pkin(3) * t757 * t788 - qJ(4) * t738 * t784 - t694) * t754;
t683 = t757 * t694 + t754 * t706;
t720 = t738 * pkin(3) - qJ(4) * t777;
t723 = (-qJDD(1) * t754 - t757 * t781) * t750;
t779 = t754 ^ 2 * t788;
t674 = -pkin(3) * t779 + t723 * qJ(4) - t738 * t720 + t683;
t749 = sin(pkin(9));
t751 = cos(pkin(9));
t715 = (-t749 * t754 + t751 * t757) * t784;
t659 = -0.2e1 * qJD(4) * t715 + t751 * t673 - t749 * t674;
t698 = t749 * t723 + t751 * t724;
t714 = (-t749 * t757 - t751 * t754) * t784;
t657 = (t714 * t738 - t698) * pkin(7) + (t714 * t715 + t737) * pkin(4) + t659;
t660 = 0.2e1 * qJD(4) * t714 + t749 * t673 + t751 * t674;
t697 = t751 * t723 - t749 * t724;
t701 = t738 * pkin(4) - t715 * pkin(7);
t713 = t714 ^ 2;
t658 = -t713 * pkin(4) + t697 * pkin(7) - t738 * t701 + t660;
t753 = sin(qJ(5));
t756 = cos(qJ(5));
t655 = t756 * t657 - t753 * t658;
t691 = t756 * t714 - t753 * t715;
t669 = t691 * qJD(5) + t753 * t697 + t756 * t698;
t692 = t753 * t714 + t756 * t715;
t680 = -t691 * mrSges(6,1) + t692 * mrSges(6,2);
t736 = qJD(5) + t738;
t684 = -t736 * mrSges(6,2) + t691 * mrSges(6,3);
t733 = qJDD(5) + t737;
t651 = m(6) * t655 + t733 * mrSges(6,1) - t669 * mrSges(6,3) - t692 * t680 + t736 * t684;
t656 = t753 * t657 + t756 * t658;
t668 = -t692 * qJD(5) + t756 * t697 - t753 * t698;
t685 = t736 * mrSges(6,1) - t692 * mrSges(6,3);
t652 = m(6) * t656 - t733 * mrSges(6,2) + t668 * mrSges(6,3) + t691 * t680 - t736 * t685;
t643 = t756 * t651 + t753 * t652;
t695 = -t714 * mrSges(5,1) + t715 * mrSges(5,2);
t699 = -t738 * mrSges(5,2) + t714 * mrSges(5,3);
t641 = m(5) * t659 + t737 * mrSges(5,1) - t698 * mrSges(5,3) - t715 * t695 + t738 * t699 + t643;
t700 = t738 * mrSges(5,1) - t715 * mrSges(5,3);
t771 = -t753 * t651 + t756 * t652;
t642 = m(5) * t660 - t737 * mrSges(5,2) + t697 * mrSges(5,3) + t714 * t695 - t738 * t700 + t771;
t637 = t751 * t641 + t749 * t642;
t682 = -t754 * t694 + t703;
t687 = Ifges(5,4) * t715 + Ifges(5,2) * t714 + Ifges(5,6) * t738;
t688 = Ifges(5,1) * t715 + Ifges(5,4) * t714 + Ifges(5,5) * t738;
t676 = Ifges(6,4) * t692 + Ifges(6,2) * t691 + Ifges(6,6) * t736;
t677 = Ifges(6,1) * t692 + Ifges(6,4) * t691 + Ifges(6,5) * t736;
t763 = -mrSges(6,1) * t655 + mrSges(6,2) * t656 - Ifges(6,5) * t669 - Ifges(6,6) * t668 - Ifges(6,3) * t733 - t692 * t676 + t691 * t677;
t793 = -mrSges(4,1) * t682 - mrSges(5,1) * t659 + mrSges(4,2) * t683 + mrSges(5,2) * t660 - Ifges(4,5) * t724 - Ifges(5,5) * t698 - Ifges(4,6) * t723 - Ifges(5,6) * t697 - pkin(3) * t637 - pkin(4) * t643 - t715 * t687 + t714 * t688 - (Ifges(4,3) + Ifges(5,3)) * t737 + t763;
t791 = mrSges(3,2) * t750;
t728 = (-mrSges(3,1) * t752 + t791) * qJD(1);
t722 = (mrSges(4,1) * t754 + mrSges(4,2) * t757) * t784;
t635 = m(4) * t682 + t737 * mrSges(4,1) - t724 * mrSges(4,3) + t738 * t719 - t722 * t777 + t637;
t772 = -t749 * t641 + t751 * t642;
t636 = m(4) * t683 - t737 * mrSges(4,2) + t723 * mrSges(4,3) - t738 * t721 - t722 * t778 + t772;
t773 = -t754 * t635 + t757 * t636;
t785 = qJDD(1) * mrSges(3,3);
t628 = m(3) * t705 + (qJD(1) * t728 + t785) * t752 + t773;
t693 = t730 * t784 - t704;
t681 = -t723 * pkin(3) - qJ(4) * t779 + t720 * t777 + qJDD(4) + t693;
t662 = -t697 * pkin(4) - t713 * pkin(7) + t715 * t701 + t681;
t767 = m(6) * t662 - t668 * mrSges(6,1) + t669 * mrSges(6,2) - t691 * t684 + t692 * t685;
t653 = m(5) * t681 - t697 * mrSges(5,1) + t698 * mrSges(5,2) - t714 * t699 + t715 * t700 + t767;
t761 = -m(4) * t693 + t723 * mrSges(4,1) - t724 * mrSges(4,2) - t653;
t650 = t761 + m(3) * t704 + (-t785 + (-t728 + t794) * qJD(1)) * t750;
t774 = t752 * t628 - t750 * t650;
t621 = m(2) * t735 - t759 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t774;
t631 = t757 * t635 + t754 * t636;
t726 = -qJDD(1) * pkin(1) + t765;
t762 = -m(3) * t726 + mrSges(3,1) * t780 - t631 + (t752 ^ 2 * t759 + t788) * mrSges(3,3);
t625 = m(2) * t734 - t759 * mrSges(2,2) + (mrSges(2,1) - t791) * qJDD(1) + t762;
t787 = t755 * t621 + t758 * t625;
t623 = t750 * t628 + t752 * t650;
t775 = t758 * t621 - t755 * t625;
t769 = Ifges(3,1) * t750 + Ifges(3,4) * t752;
t768 = Ifges(3,5) * t750 + Ifges(3,6) * t752;
t708 = Ifges(4,6) * t738 + (Ifges(4,4) * t757 - Ifges(4,2) * t754) * t784;
t709 = Ifges(4,5) * t738 + (Ifges(4,1) * t757 - Ifges(4,4) * t754) * t784;
t766 = t708 * t757 + t709 * t754;
t675 = Ifges(6,5) * t692 + Ifges(6,6) * t691 + Ifges(6,3) * t736;
t644 = -mrSges(6,1) * t662 + mrSges(6,3) * t656 + Ifges(6,4) * t669 + Ifges(6,2) * t668 + Ifges(6,6) * t733 - t692 * t675 + t736 * t677;
t645 = mrSges(6,2) * t662 - mrSges(6,3) * t655 + Ifges(6,1) * t669 + Ifges(6,4) * t668 + Ifges(6,5) * t733 + t691 * t675 - t736 * t676;
t686 = Ifges(5,5) * t715 + Ifges(5,6) * t714 + Ifges(5,3) * t738;
t632 = -mrSges(5,1) * t681 + mrSges(5,3) * t660 + Ifges(5,4) * t698 + Ifges(5,2) * t697 + Ifges(5,6) * t737 - pkin(4) * t767 + pkin(7) * t771 + t756 * t644 + t753 * t645 - t715 * t686 + t738 * t688;
t633 = mrSges(5,2) * t681 - mrSges(5,3) * t659 + Ifges(5,1) * t698 + Ifges(5,4) * t697 + Ifges(5,5) * t737 - pkin(7) * t643 - t753 * t644 + t756 * t645 + t714 * t686 - t738 * t687;
t707 = Ifges(4,3) * t738 + (Ifges(4,5) * t757 - Ifges(4,6) * t754) * t784;
t617 = -mrSges(4,1) * t693 + mrSges(4,3) * t683 + Ifges(4,4) * t724 + Ifges(4,2) * t723 + Ifges(4,6) * t737 - pkin(3) * t653 + qJ(4) * t772 + t751 * t632 + t749 * t633 - t707 * t777 + t738 * t709;
t618 = mrSges(4,2) * t693 - mrSges(4,3) * t682 + Ifges(4,1) * t724 + Ifges(4,4) * t723 + Ifges(4,5) * t737 - qJ(4) * t637 - t749 * t632 + t751 * t633 - t707 * t778 - t738 * t708;
t729 = t768 * qJD(1);
t614 = mrSges(3,2) * t726 - mrSges(3,3) * t704 - pkin(6) * t631 + t769 * qJDD(1) - t754 * t617 + t757 * t618 + t729 * t783;
t616 = (Ifges(3,4) * qJDD(1) + (-t729 - t766) * qJD(1)) * t750 + mrSges(3,3) * t705 + Ifges(3,2) * t780 - mrSges(3,1) * t726 - pkin(2) * t631 + t793;
t630 = qJDD(1) * t791 - t762;
t764 = mrSges(2,1) * t734 - mrSges(2,2) * t735 + Ifges(2,3) * qJDD(1) - pkin(1) * t630 + qJ(2) * t774 + t750 * t614 + t752 * t616;
t612 = t759 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t735 - mrSges(3,1) * t704 + mrSges(3,2) * t705 - t754 * t618 - t757 * t617 - pkin(2) * t761 - pkin(6) * t773 - pkin(1) * t623 + (Ifges(2,6) - t768) * qJDD(1) + (-pkin(2) * t794 * t750 + (-t750 * (Ifges(3,4) * t750 + Ifges(3,2) * t752) + t752 * t769) * qJD(1)) * qJD(1);
t611 = -mrSges(2,2) * g(3) - mrSges(2,3) * t734 + Ifges(2,5) * qJDD(1) - t759 * Ifges(2,6) - qJ(2) * t623 + t752 * t614 - t750 * t616;
t1 = [-m(1) * g(1) + t775; -m(1) * g(2) + t787; (-m(1) - m(2)) * g(3) + t623; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t787 + t758 * t611 - t755 * t612; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t775 + t755 * t611 + t758 * t612; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t764; t764; t630; t766 * t784 - t793; t653; -t763;];
tauJB = t1;
