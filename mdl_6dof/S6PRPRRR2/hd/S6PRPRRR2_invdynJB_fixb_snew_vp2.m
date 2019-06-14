% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRPRRR2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-05-05 00:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRPRRR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 00:22:24
% EndTime: 2019-05-05 00:22:37
% DurationCPUTime: 12.48s
% Computational Cost: add. (229389->299), mult. (426567->382), div. (0->0), fcn. (303513->14), ass. (0->134)
t769 = sin(pkin(11));
t772 = cos(pkin(11));
t756 = g(1) * t769 - g(2) * t772;
t757 = -g(1) * t772 - g(2) * t769;
t767 = -g(3) + qJDD(1);
t777 = sin(qJ(2));
t773 = cos(pkin(6));
t781 = cos(qJ(2));
t801 = t773 * t781;
t770 = sin(pkin(6));
t803 = t770 * t781;
t718 = t756 * t801 - t757 * t777 + t767 * t803;
t716 = qJDD(2) * pkin(2) + t718;
t802 = t773 * t777;
t804 = t770 * t777;
t719 = t756 * t802 + t781 * t757 + t767 * t804;
t783 = qJD(2) ^ 2;
t717 = -pkin(2) * t783 + t719;
t768 = sin(pkin(12));
t771 = cos(pkin(12));
t705 = t768 * t716 + t771 * t717;
t703 = -pkin(3) * t783 + qJDD(2) * pkin(8) + t705;
t736 = -t756 * t770 + t773 * t767;
t735 = qJDD(3) + t736;
t776 = sin(qJ(4));
t780 = cos(qJ(4));
t693 = t780 * t703 + t776 * t735;
t753 = (-pkin(4) * t780 - pkin(9) * t776) * qJD(2);
t782 = qJD(4) ^ 2;
t798 = qJD(2) * t780;
t688 = -pkin(4) * t782 + qJDD(4) * pkin(9) + t753 * t798 + t693;
t704 = t771 * t716 - t768 * t717;
t702 = -qJDD(2) * pkin(3) - t783 * pkin(8) - t704;
t797 = qJD(2) * qJD(4);
t796 = t780 * t797;
t754 = qJDD(2) * t776 + t796;
t765 = t776 * t797;
t755 = qJDD(2) * t780 - t765;
t691 = (-t754 - t796) * pkin(9) + (-t755 + t765) * pkin(4) + t702;
t775 = sin(qJ(5));
t779 = cos(qJ(5));
t683 = -t775 * t688 + t779 * t691;
t799 = qJD(2) * t776;
t750 = qJD(4) * t779 - t775 * t799;
t726 = qJD(5) * t750 + qJDD(4) * t775 + t754 * t779;
t748 = qJDD(5) - t755;
t751 = qJD(4) * t775 + t779 * t799;
t763 = qJD(5) - t798;
t681 = (t750 * t763 - t726) * pkin(10) + (t750 * t751 + t748) * pkin(5) + t683;
t684 = t779 * t688 + t775 * t691;
t725 = -qJD(5) * t751 + qJDD(4) * t779 - t754 * t775;
t734 = pkin(5) * t763 - pkin(10) * t751;
t747 = t750 ^ 2;
t682 = -pkin(5) * t747 + pkin(10) * t725 - t734 * t763 + t684;
t774 = sin(qJ(6));
t778 = cos(qJ(6));
t680 = t681 * t774 + t682 * t778;
t692 = -t776 * t703 + t735 * t780;
t687 = -qJDD(4) * pkin(4) - pkin(9) * t782 + t753 * t799 - t692;
t685 = -pkin(5) * t725 - pkin(10) * t747 + t734 * t751 + t687;
t728 = t750 * t774 + t751 * t778;
t698 = -qJD(6) * t728 + t725 * t778 - t726 * t774;
t727 = t750 * t778 - t751 * t774;
t699 = qJD(6) * t727 + t725 * t774 + t726 * t778;
t762 = qJD(6) + t763;
t706 = Ifges(7,5) * t728 + Ifges(7,6) * t727 + Ifges(7,3) * t762;
t708 = Ifges(7,1) * t728 + Ifges(7,4) * t727 + Ifges(7,5) * t762;
t744 = qJDD(6) + t748;
t668 = -mrSges(7,1) * t685 + mrSges(7,3) * t680 + Ifges(7,4) * t699 + Ifges(7,2) * t698 + Ifges(7,6) * t744 - t706 * t728 + t708 * t762;
t679 = t681 * t778 - t682 * t774;
t707 = Ifges(7,4) * t728 + Ifges(7,2) * t727 + Ifges(7,6) * t762;
t669 = mrSges(7,2) * t685 - mrSges(7,3) * t679 + Ifges(7,1) * t699 + Ifges(7,4) * t698 + Ifges(7,5) * t744 + t706 * t727 - t707 * t762;
t720 = Ifges(6,5) * t751 + Ifges(6,6) * t750 + Ifges(6,3) * t763;
t722 = Ifges(6,1) * t751 + Ifges(6,4) * t750 + Ifges(6,5) * t763;
t714 = -mrSges(7,2) * t762 + mrSges(7,3) * t727;
t715 = mrSges(7,1) * t762 - mrSges(7,3) * t728;
t788 = m(7) * t685 - t698 * mrSges(7,1) + mrSges(7,2) * t699 - t727 * t714 + t715 * t728;
t710 = -mrSges(7,1) * t727 + mrSges(7,2) * t728;
t675 = m(7) * t679 + mrSges(7,1) * t744 - mrSges(7,3) * t699 - t710 * t728 + t714 * t762;
t676 = m(7) * t680 - mrSges(7,2) * t744 + mrSges(7,3) * t698 + t710 * t727 - t715 * t762;
t792 = -t675 * t774 + t778 * t676;
t655 = -mrSges(6,1) * t687 + mrSges(6,3) * t684 + Ifges(6,4) * t726 + Ifges(6,2) * t725 + Ifges(6,6) * t748 - pkin(5) * t788 + pkin(10) * t792 + t778 * t668 + t774 * t669 - t751 * t720 + t763 * t722;
t667 = t778 * t675 + t774 * t676;
t721 = Ifges(6,4) * t751 + Ifges(6,2) * t750 + Ifges(6,6) * t763;
t656 = mrSges(6,2) * t687 - mrSges(6,3) * t683 + Ifges(6,1) * t726 + Ifges(6,4) * t725 + Ifges(6,5) * t748 - pkin(10) * t667 - t668 * t774 + t669 * t778 + t720 * t750 - t721 * t763;
t729 = -mrSges(6,1) * t750 + mrSges(6,2) * t751;
t732 = -mrSges(6,2) * t763 + mrSges(6,3) * t750;
t665 = m(6) * t683 + mrSges(6,1) * t748 - mrSges(6,3) * t726 - t729 * t751 + t732 * t763 + t667;
t733 = mrSges(6,1) * t763 - mrSges(6,3) * t751;
t666 = m(6) * t684 - mrSges(6,2) * t748 + mrSges(6,3) * t725 + t729 * t750 - t733 * t763 + t792;
t663 = -t665 * t775 + t779 * t666;
t677 = -m(6) * t687 + t725 * mrSges(6,1) - mrSges(6,2) * t726 + t750 * t732 - t733 * t751 - t788;
t742 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t776 + Ifges(5,2) * t780) * qJD(2);
t743 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t776 + Ifges(5,4) * t780) * qJD(2);
t805 = mrSges(5,1) * t692 - mrSges(5,2) * t693 + Ifges(5,5) * t754 + Ifges(5,6) * t755 + Ifges(5,3) * qJDD(4) + pkin(4) * t677 + pkin(9) * t663 + t779 * t655 + t775 * t656 + (t742 * t776 - t743 * t780) * qJD(2);
t752 = (-mrSges(5,1) * t780 + mrSges(5,2) * t776) * qJD(2);
t758 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t799;
t661 = m(5) * t693 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t755 - qJD(4) * t758 + t752 * t798 + t663;
t759 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t798;
t671 = m(5) * t692 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t754 + qJD(4) * t759 - t752 * t799 + t677;
t793 = t780 * t661 - t671 * t776;
t650 = m(4) * t705 - mrSges(4,1) * t783 - qJDD(2) * mrSges(4,2) + t793;
t662 = t665 * t779 + t666 * t775;
t786 = -m(5) * t702 + t755 * mrSges(5,1) - mrSges(5,2) * t754 - t758 * t799 + t759 * t798 - t662;
t658 = m(4) * t704 + qJDD(2) * mrSges(4,1) - mrSges(4,2) * t783 + t786;
t646 = t768 * t650 + t771 * t658;
t644 = m(3) * t718 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t783 + t646;
t794 = t771 * t650 - t658 * t768;
t645 = m(3) * t719 - mrSges(3,1) * t783 - qJDD(2) * mrSges(3,2) + t794;
t654 = t776 * t661 + t780 * t671;
t653 = m(4) * t735 + t654;
t652 = m(3) * t736 + t653;
t632 = t644 * t801 + t645 * t802 - t652 * t770;
t630 = m(2) * t756 + t632;
t637 = -t644 * t777 + t781 * t645;
t636 = m(2) * t757 + t637;
t800 = t772 * t630 + t769 * t636;
t631 = t644 * t803 + t645 * t804 + t773 * t652;
t795 = -t630 * t769 + t772 * t636;
t791 = m(2) * t767 + t631;
t741 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t776 + Ifges(5,6) * t780) * qJD(2);
t638 = mrSges(5,2) * t702 - mrSges(5,3) * t692 + Ifges(5,1) * t754 + Ifges(5,4) * t755 + Ifges(5,5) * qJDD(4) - pkin(9) * t662 - qJD(4) * t742 - t655 * t775 + t656 * t779 + t741 * t798;
t787 = -mrSges(7,1) * t679 + mrSges(7,2) * t680 - Ifges(7,5) * t699 - Ifges(7,6) * t698 - Ifges(7,3) * t744 - t728 * t707 + t727 * t708;
t784 = mrSges(6,1) * t683 - mrSges(6,2) * t684 + Ifges(6,5) * t726 + Ifges(6,6) * t725 + Ifges(6,3) * t748 + pkin(5) * t667 + t751 * t721 - t750 * t722 - t787;
t647 = -mrSges(5,1) * t702 + mrSges(5,3) * t693 + Ifges(5,4) * t754 + Ifges(5,2) * t755 + Ifges(5,6) * qJDD(4) - pkin(4) * t662 + qJD(4) * t743 - t741 * t799 - t784;
t628 = mrSges(4,2) * t735 - mrSges(4,3) * t704 + Ifges(4,5) * qJDD(2) - Ifges(4,6) * t783 - pkin(8) * t654 + t638 * t780 - t647 * t776;
t633 = -mrSges(4,1) * t735 + mrSges(4,3) * t705 + t783 * Ifges(4,5) + Ifges(4,6) * qJDD(2) - pkin(3) * t654 - t805;
t625 = -mrSges(3,1) * t736 + mrSges(3,3) * t719 + t783 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t653 + qJ(3) * t794 + t768 * t628 + t771 * t633;
t626 = mrSges(3,2) * t736 - mrSges(3,3) * t718 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t783 - qJ(3) * t646 + t628 * t771 - t633 * t768;
t789 = pkin(7) * t637 + t625 * t781 + t626 * t777;
t627 = mrSges(3,1) * t718 - mrSges(3,2) * t719 + mrSges(4,1) * t704 - mrSges(4,2) * t705 + t776 * t638 + t780 * t647 + pkin(3) * t786 + pkin(8) * t793 + pkin(2) * t646 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2);
t624 = mrSges(2,2) * t767 - mrSges(2,3) * t756 - t777 * t625 + t781 * t626 + (-t631 * t770 - t632 * t773) * pkin(7);
t623 = -mrSges(2,1) * t767 + mrSges(2,3) * t757 - pkin(1) * t631 - t770 * t627 + t773 * t789;
t1 = [-m(1) * g(1) + t795; -m(1) * g(2) + t800; -m(1) * g(3) + t791; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t800 - t769 * t623 + t772 * t624; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t795 + t772 * t623 + t769 * t624; -mrSges(1,1) * g(2) + mrSges(2,1) * t756 + mrSges(1,2) * g(1) - mrSges(2,2) * t757 + pkin(1) * t632 + t773 * t627 + t770 * t789; t791; t627; t653; t805; t784; -t787;];
tauJB  = t1;
