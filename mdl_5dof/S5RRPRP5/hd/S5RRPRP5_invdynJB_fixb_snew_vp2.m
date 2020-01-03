% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPRP5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPRP5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP5_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP5_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP5_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP5_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:54:28
% EndTime: 2019-12-31 19:54:33
% DurationCPUTime: 5.86s
% Computational Cost: add. (58441->290), mult. (136067->358), div. (0->0), fcn. (91751->8), ass. (0->115)
t793 = Ifges(5,1) + Ifges(6,1);
t782 = Ifges(5,4) - Ifges(6,5);
t791 = Ifges(6,4) + Ifges(5,5);
t792 = Ifges(5,2) + Ifges(6,3);
t789 = Ifges(5,6) - Ifges(6,6);
t790 = -Ifges(6,2) - Ifges(5,3);
t753 = sin(pkin(8));
t754 = cos(pkin(8));
t756 = sin(qJ(2));
t758 = cos(qJ(2));
t724 = (-t753 * t756 + t754 * t758) * qJD(1);
t725 = (t753 * t758 + t754 * t756) * qJD(1);
t755 = sin(qJ(4));
t785 = cos(qJ(4));
t705 = -t785 * t724 + t755 * t725;
t706 = t755 * t724 + t785 * t725;
t750 = qJD(2) + qJD(4);
t788 = t792 * t705 - t782 * t706 - t789 * t750;
t787 = -t782 * t705 + t793 * t706 + t791 * t750;
t774 = qJD(1) * qJD(2);
t735 = t756 * qJDD(1) + t758 * t774;
t757 = sin(qJ(1));
t759 = cos(qJ(1));
t742 = -t759 * g(1) - t757 * g(2);
t760 = qJD(1) ^ 2;
t730 = -t760 * pkin(1) + qJDD(1) * pkin(6) + t742;
t780 = t756 * t730;
t784 = pkin(2) * t760;
t693 = qJDD(2) * pkin(2) - t735 * qJ(3) - t780 + (qJ(3) * t774 + t756 * t784 - g(3)) * t758;
t715 = -t756 * g(3) + t758 * t730;
t736 = t758 * qJDD(1) - t756 * t774;
t776 = qJD(1) * t756;
t738 = qJD(2) * pkin(2) - qJ(3) * t776;
t752 = t758 ^ 2;
t694 = t736 * qJ(3) - qJD(2) * t738 - t752 * t784 + t715;
t663 = -0.2e1 * qJD(3) * t725 + t754 * t693 - t753 * t694;
t713 = t754 * t735 + t753 * t736;
t658 = (qJD(2) * t724 - t713) * pkin(7) + (t724 * t725 + qJDD(2)) * pkin(3) + t663;
t664 = 0.2e1 * qJD(3) * t724 + t753 * t693 + t754 * t694;
t712 = -t753 * t735 + t754 * t736;
t718 = qJD(2) * pkin(3) - t725 * pkin(7);
t723 = t724 ^ 2;
t660 = -t723 * pkin(3) + t712 * pkin(7) - qJD(2) * t718 + t664;
t656 = t755 * t658 + t785 * t660;
t675 = t706 * qJD(4) - t785 * t712 + t755 * t713;
t699 = t750 * mrSges(5,1) - t706 * mrSges(5,3);
t749 = qJDD(2) + qJDD(4);
t687 = t705 * pkin(4) - t706 * qJ(5);
t748 = t750 ^ 2;
t650 = -t748 * pkin(4) + t749 * qJ(5) + 0.2e1 * qJD(5) * t750 - t705 * t687 + t656;
t700 = -t750 * mrSges(6,1) + t706 * mrSges(6,2);
t773 = m(6) * t650 + t749 * mrSges(6,3) + t750 * t700;
t688 = t705 * mrSges(6,1) - t706 * mrSges(6,3);
t777 = -t705 * mrSges(5,1) - t706 * mrSges(5,2) - t688;
t783 = -mrSges(5,3) - mrSges(6,2);
t640 = m(5) * t656 - t749 * mrSges(5,2) + t783 * t675 - t750 * t699 + t777 * t705 + t773;
t655 = t785 * t658 - t755 * t660;
t676 = -t705 * qJD(4) + t755 * t712 + t785 * t713;
t698 = -t750 * mrSges(5,2) - t705 * mrSges(5,3);
t651 = -t749 * pkin(4) - t748 * qJ(5) + t706 * t687 + qJDD(5) - t655;
t701 = -t705 * mrSges(6,2) + t750 * mrSges(6,3);
t768 = -m(6) * t651 + t749 * mrSges(6,1) + t750 * t701;
t642 = m(5) * t655 + t749 * mrSges(5,1) + t783 * t676 + t750 * t698 + t777 * t706 + t768;
t633 = t755 * t640 + t785 * t642;
t709 = -t724 * mrSges(4,1) + t725 * mrSges(4,2);
t716 = -qJD(2) * mrSges(4,2) + t724 * mrSges(4,3);
t629 = m(4) * t663 + qJDD(2) * mrSges(4,1) - t713 * mrSges(4,3) + qJD(2) * t716 - t725 * t709 + t633;
t717 = qJD(2) * mrSges(4,1) - t725 * mrSges(4,3);
t769 = t785 * t640 - t755 * t642;
t630 = m(4) * t664 - qJDD(2) * mrSges(4,2) + t712 * mrSges(4,3) - qJD(2) * t717 + t724 * t709 + t769;
t625 = t754 * t629 + t753 * t630;
t703 = Ifges(4,4) * t725 + Ifges(4,2) * t724 + Ifges(4,6) * qJD(2);
t704 = Ifges(4,1) * t725 + Ifges(4,4) * t724 + Ifges(4,5) * qJD(2);
t714 = -t758 * g(3) - t780;
t727 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t756 + Ifges(3,2) * t758) * qJD(1);
t728 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t756 + Ifges(3,4) * t758) * qJD(1);
t647 = t676 * mrSges(6,2) + t706 * t688 - t768;
t763 = -mrSges(5,1) * t655 + mrSges(6,1) * t651 + mrSges(5,2) * t656 - mrSges(6,3) * t650 + pkin(4) * t647 - qJ(5) * t773 + t790 * t749 + t788 * t706 + (qJ(5) * t688 - t787) * t705 - t791 * t676 + (qJ(5) * mrSges(6,2) + t789) * t675;
t786 = mrSges(3,1) * t714 + mrSges(4,1) * t663 - mrSges(3,2) * t715 - mrSges(4,2) * t664 + Ifges(3,5) * t735 + Ifges(4,5) * t713 + Ifges(3,6) * t736 + Ifges(4,6) * t712 + pkin(2) * t625 + pkin(3) * t633 + (t756 * t727 - t758 * t728) * qJD(1) + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + t725 * t703 - t724 * t704 - t763;
t734 = (-mrSges(3,1) * t758 + mrSges(3,2) * t756) * qJD(1);
t775 = qJD(1) * t758;
t740 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t775;
t623 = m(3) * t714 + qJDD(2) * mrSges(3,1) - t735 * mrSges(3,3) + qJD(2) * t740 - t734 * t776 + t625;
t739 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t776;
t770 = -t753 * t629 + t754 * t630;
t624 = m(3) * t715 - qJDD(2) * mrSges(3,2) + t736 * mrSges(3,3) - qJD(2) * t739 + t734 * t775 + t770;
t771 = -t756 * t623 + t758 * t624;
t615 = m(2) * t742 - t760 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t771;
t741 = t757 * g(1) - t759 * g(2);
t766 = -qJDD(1) * pkin(1) - t741;
t697 = -t736 * pkin(2) + qJDD(3) + t738 * t776 + (-qJ(3) * t752 - pkin(6)) * t760 + t766;
t662 = -t712 * pkin(3) - t723 * pkin(7) + t725 * t718 + t697;
t653 = -0.2e1 * qJD(5) * t706 + (t705 * t750 - t676) * qJ(5) + (t706 * t750 + t675) * pkin(4) + t662;
t643 = m(6) * t653 + t675 * mrSges(6,1) - t676 * mrSges(6,3) - t706 * t700 + t705 * t701;
t764 = m(5) * t662 + t675 * mrSges(5,1) + t676 * mrSges(5,2) + t705 * t698 + t706 * t699 + t643;
t637 = m(4) * t697 - t712 * mrSges(4,1) + t713 * mrSges(4,2) - t724 * t716 + t725 * t717 + t764;
t729 = -t760 * pkin(6) + t766;
t762 = -m(3) * t729 + t736 * mrSges(3,1) - t735 * mrSges(3,2) - t739 * t776 + t740 * t775 - t637;
t635 = m(2) * t741 + qJDD(1) * mrSges(2,1) - t760 * mrSges(2,2) + t762;
t779 = t757 * t615 + t759 * t635;
t617 = t758 * t623 + t756 * t624;
t778 = t789 * t705 - t791 * t706 + t790 * t750;
t772 = t759 * t615 - t757 * t635;
t631 = -mrSges(5,1) * t662 - mrSges(6,1) * t653 + mrSges(6,2) * t650 + mrSges(5,3) * t656 - pkin(4) * t643 - t792 * t675 + t782 * t676 + t778 * t706 + t789 * t749 + t787 * t750;
t632 = mrSges(5,2) * t662 + mrSges(6,2) * t651 - mrSges(5,3) * t655 - mrSges(6,3) * t653 - qJ(5) * t643 - t782 * t675 + t793 * t676 + t778 * t705 + t791 * t749 + t788 * t750;
t702 = Ifges(4,5) * t725 + Ifges(4,6) * t724 + Ifges(4,3) * qJD(2);
t618 = -mrSges(4,1) * t697 + mrSges(4,3) * t664 + Ifges(4,4) * t713 + Ifges(4,2) * t712 + Ifges(4,6) * qJDD(2) - pkin(3) * t764 + pkin(7) * t769 + qJD(2) * t704 + t785 * t631 + t755 * t632 - t725 * t702;
t619 = mrSges(4,2) * t697 - mrSges(4,3) * t663 + Ifges(4,1) * t713 + Ifges(4,4) * t712 + Ifges(4,5) * qJDD(2) - pkin(7) * t633 - qJD(2) * t703 - t755 * t631 + t785 * t632 + t724 * t702;
t726 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t756 + Ifges(3,6) * t758) * qJD(1);
t609 = -mrSges(3,1) * t729 + mrSges(3,3) * t715 + Ifges(3,4) * t735 + Ifges(3,2) * t736 + Ifges(3,6) * qJDD(2) - pkin(2) * t637 + qJ(3) * t770 + qJD(2) * t728 + t754 * t618 + t753 * t619 - t726 * t776;
t612 = mrSges(3,2) * t729 - mrSges(3,3) * t714 + Ifges(3,1) * t735 + Ifges(3,4) * t736 + Ifges(3,5) * qJDD(2) - qJ(3) * t625 - qJD(2) * t727 - t753 * t618 + t754 * t619 + t726 * t775;
t765 = mrSges(2,1) * t741 - mrSges(2,2) * t742 + Ifges(2,3) * qJDD(1) + pkin(1) * t762 + pkin(6) * t771 + t758 * t609 + t756 * t612;
t610 = mrSges(2,1) * g(3) + mrSges(2,3) * t742 + t760 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t617 - t786;
t607 = -mrSges(2,2) * g(3) - mrSges(2,3) * t741 + Ifges(2,5) * qJDD(1) - t760 * Ifges(2,6) - pkin(6) * t617 - t756 * t609 + t758 * t612;
t1 = [-m(1) * g(1) + t772; -m(1) * g(2) + t779; (-m(1) - m(2)) * g(3) + t617; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t779 + t759 * t607 - t757 * t610; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t772 + t757 * t607 + t759 * t610; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t765; t765; t786; t637; -t763; t647;];
tauJB = t1;
