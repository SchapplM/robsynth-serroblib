% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
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
% Datum: 2019-05-05 13:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPPPRR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:36:07
% EndTime: 2019-05-05 13:36:11
% DurationCPUTime: 4.55s
% Computational Cost: add. (55279->271), mult. (115006->327), div. (0->0), fcn. (73100->10), ass. (0->124)
t757 = sin(qJ(1));
t760 = cos(qJ(1));
t728 = t757 * g(1) - t760 * g(2);
t725 = qJDD(1) * pkin(1) + t728;
t729 = -t760 * g(1) - t757 * g(2);
t762 = qJD(1) ^ 2;
t726 = -t762 * pkin(1) + t729;
t752 = sin(pkin(9));
t754 = cos(pkin(9));
t711 = t754 * t725 - t752 * t726;
t770 = -t762 * qJ(3) + qJDD(3) - t711;
t795 = -pkin(2) - qJ(4);
t803 = -(2 * qJD(1) * qJD(4)) + qJDD(1) * t795 + t770;
t712 = t752 * t725 + t754 * t726;
t802 = qJDD(1) * qJ(3) + (2 * qJD(3) * qJD(1)) + t712;
t751 = sin(pkin(10));
t753 = cos(pkin(10));
t756 = sin(qJ(5));
t759 = cos(qJ(5));
t774 = t751 * t759 + t753 * t756;
t718 = t774 * qJD(1);
t695 = t762 * pkin(2) - t802;
t801 = -m(4) * t695 + t762 * mrSges(4,2) + qJDD(1) * mrSges(4,3);
t773 = -t751 * t756 + t753 * t759;
t719 = t773 * qJD(1);
t789 = t719 * qJD(5);
t708 = -qJDD(1) * t774 - t789;
t799 = pkin(4) * t762;
t798 = mrSges(3,1) - mrSges(4,2);
t797 = -Ifges(4,4) + Ifges(3,5);
t796 = Ifges(3,6) - Ifges(4,5);
t748 = -g(3) + qJDD(2);
t786 = qJDD(1) * t753;
t793 = t803 * t753;
t673 = -pkin(7) * t786 + (-t753 * t799 - t748) * t751 + t793;
t685 = t753 * t748 + t803 * t751;
t740 = t751 ^ 2;
t787 = qJDD(1) * t751;
t674 = -pkin(7) * t787 - t740 * t799 + t685;
t670 = t756 * t673 + t759 * t674;
t702 = t718 * mrSges(6,1) + t719 * mrSges(6,2);
t716 = qJD(5) * mrSges(6,1) - t719 * mrSges(6,3);
t707 = t718 * pkin(5) - t719 * pkin(8);
t761 = qJD(5) ^ 2;
t666 = -t761 * pkin(5) + qJDD(5) * pkin(8) - t718 * t707 + t670;
t772 = qJDD(4) + t802;
t792 = -t753 ^ 2 - t740;
t683 = pkin(4) * t787 + (pkin(7) * t792 + t795) * t762 + t772;
t790 = t718 * qJD(5);
t709 = qJDD(1) * t773 - t790;
t667 = (-t709 + t790) * pkin(8) + (-t708 + t789) * pkin(5) + t683;
t755 = sin(qJ(6));
t758 = cos(qJ(6));
t663 = -t755 * t666 + t758 * t667;
t713 = t758 * qJD(5) - t755 * t719;
t681 = t713 * qJD(6) + t755 * qJDD(5) + t758 * t709;
t714 = t755 * qJD(5) + t758 * t719;
t687 = -t713 * mrSges(7,1) + t714 * mrSges(7,2);
t717 = qJD(6) + t718;
t692 = -t717 * mrSges(7,2) + t713 * mrSges(7,3);
t706 = qJDD(6) - t708;
t660 = m(7) * t663 + t706 * mrSges(7,1) - t681 * mrSges(7,3) - t714 * t687 + t717 * t692;
t664 = t758 * t666 + t755 * t667;
t680 = -t714 * qJD(6) + t758 * qJDD(5) - t755 * t709;
t693 = t717 * mrSges(7,1) - t714 * mrSges(7,3);
t661 = m(7) * t664 - t706 * mrSges(7,2) + t680 * mrSges(7,3) + t713 * t687 - t717 * t693;
t778 = -t755 * t660 + t758 * t661;
t648 = m(6) * t670 - qJDD(5) * mrSges(6,2) + t708 * mrSges(6,3) - qJD(5) * t716 - t718 * t702 + t778;
t669 = t759 * t673 - t756 * t674;
t715 = -qJD(5) * mrSges(6,2) - t718 * mrSges(6,3);
t665 = -qJDD(5) * pkin(5) - t761 * pkin(8) + t719 * t707 - t669;
t768 = -m(7) * t665 + t680 * mrSges(7,1) - t681 * mrSges(7,2) + t713 * t692 - t714 * t693;
t656 = m(6) * t669 + qJDD(5) * mrSges(6,1) - t709 * mrSges(6,3) + qJD(5) * t715 - t719 * t702 + t768;
t640 = t756 * t648 + t759 * t656;
t684 = -t751 * t748 + t793;
t771 = -qJDD(1) * mrSges(5,3) - t762 * (mrSges(5,1) * t751 + mrSges(5,2) * t753);
t638 = m(5) * t684 + t753 * t771 + t640;
t779 = t759 * t648 - t756 * t656;
t639 = m(5) * t685 + t751 * t771 + t779;
t634 = t753 * t638 + t751 * t639;
t697 = -qJDD(1) * pkin(2) + t770;
t769 = -m(4) * t697 + t762 * mrSges(4,3) - t634;
t629 = m(3) * t711 - t762 * mrSges(3,2) + qJDD(1) * t798 + t769;
t691 = t762 * t795 + t772;
t650 = t758 * t660 + t755 * t661;
t767 = m(6) * t683 - t708 * mrSges(6,1) + t709 * mrSges(6,2) + t718 * t715 + t719 * t716 + t650;
t766 = m(5) * t691 + mrSges(5,1) * t787 + mrSges(5,2) * t786 + t767;
t783 = t792 * mrSges(5,3);
t643 = t766 + m(3) * t712 - qJDD(1) * mrSges(3,2) + (-mrSges(3,1) + t783) * t762 + t801;
t627 = t754 * t629 + t752 * t643;
t624 = m(2) * t728 + qJDD(1) * mrSges(2,1) - t762 * mrSges(2,2) + t627;
t781 = -t752 * t629 + t754 * t643;
t625 = m(2) * t729 - t762 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t781;
t794 = t760 * t624 + t757 * t625;
t775 = Ifges(5,5) * t753 - Ifges(5,6) * t751;
t791 = t762 * t775;
t782 = -t757 * t624 + t760 * t625;
t780 = -t751 * t638 + t753 * t639;
t633 = m(4) * t748 + t780;
t777 = Ifges(5,1) * t753 - Ifges(5,4) * t751;
t776 = Ifges(5,4) * t753 - Ifges(5,2) * t751;
t632 = m(3) * t748 + t633;
t675 = Ifges(7,5) * t714 + Ifges(7,6) * t713 + Ifges(7,3) * t717;
t677 = Ifges(7,1) * t714 + Ifges(7,4) * t713 + Ifges(7,5) * t717;
t653 = -mrSges(7,1) * t665 + mrSges(7,3) * t664 + Ifges(7,4) * t681 + Ifges(7,2) * t680 + Ifges(7,6) * t706 - t714 * t675 + t717 * t677;
t676 = Ifges(7,4) * t714 + Ifges(7,2) * t713 + Ifges(7,6) * t717;
t654 = mrSges(7,2) * t665 - mrSges(7,3) * t663 + Ifges(7,1) * t681 + Ifges(7,4) * t680 + Ifges(7,5) * t706 + t713 * t675 - t717 * t676;
t699 = Ifges(6,4) * t719 - Ifges(6,2) * t718 + Ifges(6,6) * qJD(5);
t700 = Ifges(6,1) * t719 - Ifges(6,4) * t718 + Ifges(6,5) * qJD(5);
t765 = mrSges(6,1) * t669 - mrSges(6,2) * t670 + Ifges(6,5) * t709 + Ifges(6,6) * t708 + Ifges(6,3) * qJDD(5) + pkin(5) * t768 + pkin(8) * t778 + t758 * t653 + t755 * t654 + t719 * t699 + t718 * t700;
t764 = mrSges(7,1) * t663 - mrSges(7,2) * t664 + Ifges(7,5) * t681 + Ifges(7,6) * t680 + Ifges(7,3) * t706 + t714 * t676 - t713 * t677;
t645 = t762 * t783 + t766;
t698 = Ifges(6,5) * t719 - Ifges(6,6) * t718 + Ifges(6,3) * qJD(5);
t635 = mrSges(6,2) * t683 - mrSges(6,3) * t669 + Ifges(6,1) * t709 + Ifges(6,4) * t708 + Ifges(6,5) * qJDD(5) - pkin(8) * t650 - qJD(5) * t699 - t755 * t653 + t758 * t654 - t718 * t698;
t636 = -mrSges(6,1) * t683 + mrSges(6,3) * t670 + Ifges(6,4) * t709 + Ifges(6,2) * t708 + Ifges(6,6) * qJDD(5) - pkin(5) * t650 + qJD(5) * t700 - t719 * t698 - t764;
t618 = -mrSges(5,1) * t691 + mrSges(5,3) * t685 - pkin(4) * t767 + pkin(7) * t779 + qJDD(1) * t776 + t756 * t635 + t759 * t636 - t753 * t791;
t620 = mrSges(5,2) * t691 - mrSges(5,3) * t684 - pkin(7) * t640 + qJDD(1) * t777 + t759 * t635 - t756 * t636 - t751 * t791;
t631 = qJDD(1) * mrSges(4,2) - t769;
t763 = -mrSges(2,2) * t729 - mrSges(3,2) * t712 - mrSges(4,3) * t695 + pkin(1) * t627 - pkin(2) * t631 - qJ(4) * t634 - t751 * t618 + t753 * t620 + qJ(3) * (t645 + t801) + mrSges(4,2) * t697 + mrSges(3,1) * t711 + mrSges(2,1) * t728 + (Ifges(4,1) + Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t617 = t765 + mrSges(4,1) * t697 + mrSges(5,1) * t684 - mrSges(5,2) * t685 + pkin(3) * t634 + pkin(4) * t640 - qJ(3) * t633 - mrSges(3,3) * t711 + (t775 + t797) * qJDD(1) + (mrSges(3,2) - mrSges(4,3)) * t748 + (t751 * t777 + t753 * t776 - t796) * t762;
t616 = -mrSges(4,1) * t695 + mrSges(3,3) * t712 - pkin(2) * t633 + pkin(3) * t645 - qJ(4) * t780 + qJDD(1) * t796 - t753 * t618 - t751 * t620 - t748 * t798 + t762 * t797;
t615 = -mrSges(2,2) * g(3) - mrSges(2,3) * t728 + Ifges(2,5) * qJDD(1) - t762 * Ifges(2,6) - qJ(2) * t627 - t752 * t616 + t754 * t617;
t614 = mrSges(2,1) * g(3) + mrSges(2,3) * t729 + t762 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t632 + qJ(2) * t781 + t754 * t616 + t752 * t617;
t1 = [-m(1) * g(1) + t782; -m(1) * g(2) + t794; (-m(1) - m(2)) * g(3) + t632; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t794 - t757 * t614 + t760 * t615; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t782 + t760 * t614 + t757 * t615; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t763; t763; t632; t631; t645; t765; t764;];
tauJB  = t1;
