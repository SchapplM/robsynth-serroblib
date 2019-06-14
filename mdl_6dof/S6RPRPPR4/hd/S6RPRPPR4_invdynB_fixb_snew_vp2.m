% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRPPR4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-05-05 16:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRPPR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:48:20
% EndTime: 2019-05-05 16:48:29
% DurationCPUTime: 7.93s
% Computational Cost: add. (102140->344), mult. (252882->417), div. (0->0), fcn. (185434->10), ass. (0->145)
t801 = -2 * qJD(4);
t800 = Ifges(5,1) + Ifges(6,1);
t791 = Ifges(5,4) - Ifges(6,5);
t799 = -Ifges(5,5) - Ifges(6,4);
t798 = Ifges(5,2) + Ifges(6,3);
t797 = -Ifges(6,2) - Ifges(5,3);
t789 = Ifges(5,6) - Ifges(6,6);
t750 = qJD(1) ^ 2;
t746 = sin(qJ(1));
t748 = cos(qJ(1));
t728 = -t748 * g(1) - t746 * g(2);
t724 = -t750 * pkin(1) + qJDD(1) * qJ(2) + t728;
t742 = sin(pkin(9));
t743 = cos(pkin(9));
t775 = qJD(1) * qJD(2);
t771 = -t743 * g(3) - 0.2e1 * t742 * t775;
t793 = pkin(2) * t743;
t683 = (-pkin(7) * qJDD(1) + t750 * t793 - t724) * t742 + t771;
t709 = -t742 * g(3) + (t724 + 0.2e1 * t775) * t743;
t773 = qJDD(1) * t743;
t740 = t743 ^ 2;
t785 = t740 * t750;
t693 = -pkin(2) * t785 + pkin(7) * t773 + t709;
t745 = sin(qJ(3));
t794 = cos(qJ(3));
t658 = t745 * t683 + t794 * t693;
t772 = t743 * t794;
t776 = t742 * qJD(1);
t722 = -qJD(1) * t772 + t745 * t776;
t757 = t794 * t742 + t743 * t745;
t723 = t757 * qJD(1);
t699 = t722 * pkin(3) - t723 * qJ(4);
t749 = qJD(3) ^ 2;
t648 = -t749 * pkin(3) + qJDD(3) * qJ(4) - t722 * t699 + t658;
t739 = t742 ^ 2;
t727 = t746 * g(1) - t748 * g(2);
t764 = qJDD(2) - t727;
t705 = (-pkin(1) - t793) * qJDD(1) + (-qJ(2) + (-t739 - t740) * pkin(7)) * t750 + t764;
t774 = qJDD(1) * t742;
t777 = t723 * qJD(3);
t706 = -qJDD(1) * t772 + t745 * t774 + t777;
t778 = t722 * qJD(3);
t707 = t757 * qJDD(1) - t778;
t650 = (-t707 + t778) * qJ(4) + (t706 + t777) * pkin(3) + t705;
t741 = sin(pkin(10));
t787 = cos(pkin(10));
t714 = t741 * qJD(3) + t787 * t723;
t639 = -t741 * t648 + t787 * t650 + t714 * t801;
t692 = t741 * qJDD(3) + t787 * t707;
t657 = t794 * t683 - t745 * t693;
t754 = qJDD(3) * pkin(3) + t749 * qJ(4) - t723 * t699 - qJDD(4) + t657;
t713 = -t787 * qJD(3) + t741 * t723;
t786 = t713 * t722;
t796 = (-t692 + t786) * qJ(5) - t754;
t795 = 2 * qJD(5);
t792 = -mrSges(5,3) - mrSges(6,2);
t788 = mrSges(3,2) * t742;
t700 = t722 * mrSges(4,1) + t723 * mrSges(4,2);
t716 = qJD(3) * mrSges(4,1) - t723 * mrSges(4,3);
t640 = t787 * t648 + t741 * t650 + t713 * t801;
t688 = t722 * mrSges(5,1) - t714 * mrSges(5,3);
t691 = -t787 * qJDD(3) + t741 * t707;
t675 = t713 * pkin(4) - t714 * qJ(5);
t721 = t722 ^ 2;
t636 = -t721 * pkin(4) + t706 * qJ(5) - t713 * t675 + t722 * t795 + t640;
t689 = -t722 * mrSges(6,1) + t714 * mrSges(6,2);
t637 = -t706 * pkin(4) - t721 * qJ(5) + t714 * t675 + qJDD(5) - t639;
t631 = (-t692 - t786) * pkin(8) + (t713 * t714 - t706) * pkin(5) + t637;
t690 = -t722 * pkin(5) - t714 * pkin(8);
t712 = t713 ^ 2;
t632 = -t712 * pkin(5) + t691 * pkin(8) + t722 * t690 + t636;
t744 = sin(qJ(6));
t747 = cos(qJ(6));
t629 = t747 * t631 - t744 * t632;
t673 = t747 * t713 - t744 * t714;
t647 = t673 * qJD(6) + t744 * t691 + t747 * t692;
t674 = t744 * t713 + t747 * t714;
t656 = -t673 * mrSges(7,1) + t674 * mrSges(7,2);
t719 = qJD(6) - t722;
t661 = -t719 * mrSges(7,2) + t673 * mrSges(7,3);
t704 = qJDD(6) - t706;
t626 = m(7) * t629 + t704 * mrSges(7,1) - t647 * mrSges(7,3) - t674 * t656 + t719 * t661;
t630 = t744 * t631 + t747 * t632;
t646 = -t674 * qJD(6) + t747 * t691 - t744 * t692;
t662 = t719 * mrSges(7,1) - t674 * mrSges(7,3);
t627 = m(7) * t630 - t704 * mrSges(7,2) + t646 * mrSges(7,3) + t673 * t656 - t719 * t662;
t766 = -t744 * t626 + t747 * t627;
t756 = m(6) * t636 + t706 * mrSges(6,3) + t722 * t689 + t766;
t676 = t713 * mrSges(6,1) - t714 * mrSges(6,3);
t780 = -t713 * mrSges(5,1) - t714 * mrSges(5,2) - t676;
t617 = m(5) * t640 - t706 * mrSges(5,2) - t722 * t688 + t792 * t691 + t780 * t713 + t756;
t687 = -t722 * mrSges(5,2) - t713 * mrSges(5,3);
t620 = t747 * t626 + t744 * t627;
t686 = -t713 * mrSges(6,2) + t722 * mrSges(6,3);
t755 = -m(6) * t637 + t706 * mrSges(6,1) + t722 * t686 - t620;
t619 = m(5) * t639 + t706 * mrSges(5,1) + t722 * t687 + t792 * t692 + t780 * t714 + t755;
t767 = t787 * t617 - t741 * t619;
t613 = m(4) * t658 - qJDD(3) * mrSges(4,2) - t706 * mrSges(4,3) - qJD(3) * t716 - t722 * t700 + t767;
t715 = -qJD(3) * mrSges(4,2) - t722 * mrSges(4,3);
t638 = -0.2e1 * qJD(5) * t714 + (t714 * t722 + t691) * pkin(4) + t796;
t634 = -t712 * pkin(8) + (-pkin(4) - pkin(5)) * t691 + (-pkin(4) * t722 + t690 + t795) * t714 - t796;
t759 = -m(7) * t634 + t646 * mrSges(7,1) - t647 * mrSges(7,2) + t673 * t661 - t674 * t662;
t628 = m(6) * t638 + t691 * mrSges(6,1) - t692 * mrSges(6,3) + t713 * t686 - t714 * t689 + t759;
t751 = m(5) * t754 - t691 * mrSges(5,1) - t692 * mrSges(5,2) - t713 * t687 - t714 * t688 - t628;
t624 = m(4) * t657 + qJDD(3) * mrSges(4,1) - t707 * mrSges(4,3) + qJD(3) * t715 - t723 * t700 + t751;
t608 = t745 * t613 + t794 * t624;
t708 = -t742 * t724 + t771;
t758 = mrSges(3,3) * qJDD(1) + t750 * (-mrSges(3,1) * t743 + t788);
t606 = m(3) * t708 - t758 * t742 + t608;
t768 = t794 * t613 - t745 * t624;
t607 = m(3) * t709 + t758 * t743 + t768;
t769 = -t742 * t606 + t743 * t607;
t599 = m(2) * t728 - t750 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t769;
t720 = -qJDD(1) * pkin(1) - t750 * qJ(2) + t764;
t614 = t741 * t617 + t787 * t619;
t753 = m(4) * t705 + t706 * mrSges(4,1) + t707 * mrSges(4,2) + t722 * t715 + t723 * t716 + t614;
t752 = -m(3) * t720 + mrSges(3,1) * t773 - t753 + (t739 * t750 + t785) * mrSges(3,3);
t610 = (mrSges(2,1) - t788) * qJDD(1) + t752 - t750 * mrSges(2,2) + m(2) * t727;
t784 = t746 * t599 + t748 * t610;
t600 = t743 * t606 + t742 * t607;
t783 = t713 * t798 - t714 * t791 - t722 * t789;
t782 = t713 * t789 + t714 * t799 + t722 * t797;
t781 = -t713 * t791 + t714 * t800 - t722 * t799;
t770 = t748 * t599 - t746 * t610;
t763 = Ifges(3,1) * t742 + Ifges(3,4) * t743;
t762 = Ifges(3,4) * t742 + Ifges(3,2) * t743;
t761 = Ifges(3,5) * t742 + Ifges(3,6) * t743;
t726 = t761 * qJD(1);
t696 = Ifges(4,1) * t723 - Ifges(4,4) * t722 + Ifges(4,5) * qJD(3);
t695 = Ifges(4,4) * t723 - Ifges(4,2) * t722 + Ifges(4,6) * qJD(3);
t694 = Ifges(4,5) * t723 - Ifges(4,6) * t722 + Ifges(4,3) * qJD(3);
t653 = Ifges(7,1) * t674 + Ifges(7,4) * t673 + Ifges(7,5) * t719;
t652 = Ifges(7,4) * t674 + Ifges(7,2) * t673 + Ifges(7,6) * t719;
t651 = Ifges(7,5) * t674 + Ifges(7,6) * t673 + Ifges(7,3) * t719;
t622 = mrSges(7,2) * t634 - mrSges(7,3) * t629 + Ifges(7,1) * t647 + Ifges(7,4) * t646 + Ifges(7,5) * t704 + t673 * t651 - t719 * t652;
t621 = -mrSges(7,1) * t634 + mrSges(7,3) * t630 + Ifges(7,4) * t647 + Ifges(7,2) * t646 + Ifges(7,6) * t704 - t674 * t651 + t719 * t653;
t602 = -mrSges(5,2) * t754 + mrSges(6,2) * t637 - mrSges(5,3) * t639 - mrSges(6,3) * t638 - pkin(8) * t620 - qJ(5) * t628 - t744 * t621 + t747 * t622 - t791 * t691 + t692 * t800 - t706 * t799 + t782 * t713 + t783 * t722;
t601 = mrSges(5,1) * t754 - mrSges(6,1) * t638 + mrSges(6,2) * t636 + mrSges(5,3) * t640 - pkin(4) * t628 - pkin(5) * t759 - pkin(8) * t766 - t747 * t621 - t744 * t622 - t691 * t798 + t791 * t692 + t789 * t706 + t782 * t714 + t781 * t722;
t596 = -pkin(4) * t755 - qJ(5) * t756 + (-Ifges(4,2) + t797) * t706 + (pkin(4) * mrSges(6,2) + t799) * t692 + Ifges(4,6) * qJDD(3) - t723 * t694 + Ifges(7,3) * t704 - mrSges(4,1) * t705 + Ifges(4,4) * t707 + qJD(3) * t696 - t673 * t653 + t674 * t652 + Ifges(7,5) * t647 + mrSges(4,3) * t658 + Ifges(7,6) * t646 + mrSges(6,1) * t637 - mrSges(5,1) * t639 + mrSges(5,2) * t640 - mrSges(6,3) * t636 - mrSges(7,2) * t630 + mrSges(7,1) * t629 + pkin(5) * t620 - pkin(3) * t614 + (qJ(5) * t676 - t781) * t713 + (pkin(4) * t676 + t783) * t714 + (qJ(5) * mrSges(6,2) + t789) * t691;
t595 = mrSges(4,2) * t705 - mrSges(4,3) * t657 + Ifges(4,1) * t707 - Ifges(4,4) * t706 + Ifges(4,5) * qJDD(3) - qJ(4) * t614 - qJD(3) * t695 - t741 * t601 + t787 * t602 - t722 * t694;
t594 = mrSges(2,1) * g(3) - pkin(1) * t600 + mrSges(2,3) * t728 - pkin(2) * t608 - mrSges(3,1) * t708 + mrSges(3,2) * t709 - qJ(4) * t767 - t741 * t602 - t787 * t601 - pkin(3) * t751 - Ifges(4,5) * t707 + Ifges(4,6) * t706 - Ifges(4,3) * qJDD(3) - t723 * t695 - t722 * t696 - mrSges(4,1) * t657 + mrSges(4,2) * t658 + (Ifges(2,6) - t761) * qJDD(1) + (-t742 * t762 + t743 * t763 + Ifges(2,5)) * t750;
t593 = t743 * qJD(1) * t726 + mrSges(3,2) * t720 - mrSges(3,3) * t708 - pkin(7) * t608 + t763 * qJDD(1) + t794 * t595 - t745 * t596;
t592 = -mrSges(3,1) * t720 + mrSges(3,3) * t709 - pkin(2) * t753 + pkin(7) * t768 + t762 * qJDD(1) + t745 * t595 + t794 * t596 - t726 * t776;
t591 = -mrSges(2,2) * g(3) - mrSges(2,3) * t727 + Ifges(2,5) * qJDD(1) - t750 * Ifges(2,6) - qJ(2) * t600 - t742 * t592 + t743 * t593;
t1 = [-m(1) * g(1) + t770; -m(1) * g(2) + t784; (-m(1) - m(2)) * g(3) + t600; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t784 + t748 * t591 - t746 * t594; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t770 + t746 * t591 + t748 * t594; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t727 - mrSges(2,2) * t728 + t742 * t593 + t743 * t592 + pkin(1) * (-mrSges(3,2) * t774 + t752) + qJ(2) * t769;];
tauB  = t1;
