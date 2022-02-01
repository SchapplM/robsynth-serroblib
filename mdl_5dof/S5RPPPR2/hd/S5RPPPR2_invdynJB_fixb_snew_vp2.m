% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPPPR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR2_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR2_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR2_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 08:59:26
% EndTime: 2022-01-23 08:59:31
% DurationCPUTime: 5.29s
% Computational Cost: add. (44592->280), mult. (124905->382), div. (0->0), fcn. (84227->10), ass. (0->133)
t744 = sin(pkin(7));
t747 = cos(pkin(7));
t749 = sin(qJ(1));
t751 = cos(qJ(1));
t728 = -t751 * g(1) - t749 * g(2);
t752 = qJD(1) ^ 2;
t799 = -t752 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t728;
t698 = -t744 * g(3) + t799 * t747;
t765 = -pkin(2) * t747 - qJ(3) * t744;
t721 = t765 * qJD(1);
t785 = t747 * qJD(1);
t689 = t721 * t785 + t698;
t743 = sin(pkin(8));
t746 = cos(pkin(8));
t727 = t749 * g(1) - t751 * g(2);
t759 = -t752 * qJ(2) + qJDD(2) - t727;
t786 = t744 * qJD(1);
t798 = (-pkin(1) + t765) * qJDD(1) + t759 - 0.2e1 * qJD(3) * t786;
t669 = -t743 * t689 + t798 * t746;
t697 = -t747 * g(3) - t799 * t744;
t742 = sin(pkin(9));
t745 = cos(pkin(9));
t790 = t744 * t746;
t760 = t742 * t790 + t745 * t747;
t710 = t760 * qJD(1);
t708 = t760 * qJDD(1);
t797 = 2 * qJD(4);
t796 = mrSges(3,2) * t744;
t795 = Ifges(4,4) * t746;
t794 = Ifges(4,6) * t747;
t793 = t744 ^ 2 * t752;
t792 = t747 ^ 2 * t752;
t791 = t743 * t744;
t789 = t747 * t752;
t722 = (-mrSges(3,1) * t747 + t796) * qJD(1);
t670 = t746 * t689 + t798 * t743;
t769 = mrSges(4,1) * t743 + mrSges(4,2) * t746;
t715 = t769 * t786;
t762 = -mrSges(4,1) * t747 - mrSges(4,3) * t790;
t719 = t762 * qJD(1);
t761 = mrSges(4,2) * t747 - mrSges(4,3) * t791;
t714 = (pkin(3) * t743 - qJ(4) * t746) * t786;
t780 = t743 * t786;
t782 = qJDD(1) * t747;
t667 = -pkin(3) * t792 - qJ(4) * t782 - t714 * t780 + t670;
t688 = t721 * t786 + qJDD(3) - t697;
t675 = ((-qJDD(1) * t746 - t743 * t789) * qJ(4) + (qJDD(1) * t743 - t746 * t789) * pkin(3)) * t744 + t688;
t663 = t745 * t667 + t742 * t675 - t710 * t797;
t779 = t746 * t786;
t711 = -t742 * t785 + t745 * t779;
t690 = t710 * mrSges(5,1) + t711 * mrSges(5,2);
t696 = mrSges(5,1) * t780 - t711 * mrSges(5,3);
t691 = t710 * pkin(4) - t711 * pkin(6);
t783 = qJDD(1) * t744;
t776 = t743 * t783;
t781 = t743 ^ 2 * t793;
t661 = -pkin(4) * t781 + pkin(6) * t776 - t710 * t691 + t663;
t666 = pkin(3) * t782 - qJ(4) * t792 + t714 * t779 + qJDD(4) - t669;
t709 = (-t742 * t747 + t745 * t790) * qJDD(1);
t664 = (t710 * t780 - t709) * pkin(6) + (t711 * t780 + t708) * pkin(4) + t666;
t748 = sin(qJ(5));
t750 = cos(qJ(5));
t658 = -t748 * t661 + t750 * t664;
t692 = -t748 * t711 + t750 * t780;
t693 = t750 * t711 + t748 * t780;
t677 = -t692 * mrSges(6,1) + t693 * mrSges(6,2);
t679 = t692 * qJD(5) + t750 * t709 + t748 * t776;
t707 = qJD(5) + t710;
t680 = -t707 * mrSges(6,2) + t692 * mrSges(6,3);
t706 = qJDD(5) + t708;
t656 = m(6) * t658 + t706 * mrSges(6,1) - t679 * mrSges(6,3) - t693 * t677 + t707 * t680;
t659 = t750 * t661 + t748 * t664;
t678 = -t693 * qJD(5) - t748 * t709 + t750 * t776;
t681 = t707 * mrSges(6,1) - t693 * mrSges(6,3);
t657 = m(6) * t659 - t706 * mrSges(6,2) + t678 * mrSges(6,3) + t692 * t677 - t707 * t681;
t771 = -t748 * t656 + t750 * t657;
t649 = m(5) * t663 - t708 * mrSges(5,3) - t710 * t690 + (-mrSges(5,2) * qJDD(1) - qJD(1) * t696) * t791 + t771;
t764 = t742 * t667 - t745 * t675;
t662 = -0.2e1 * qJD(4) * t711 - t764;
t695 = -mrSges(5,2) * t780 - t710 * mrSges(5,3);
t660 = -pkin(4) * t776 - pkin(6) * t781 + (t797 + t691) * t711 + t764;
t756 = -m(6) * t660 + t678 * mrSges(6,1) - t679 * mrSges(6,2) + t692 * t680 - t693 * t681;
t654 = m(5) * t662 - t709 * mrSges(5,3) - t711 * t690 + (mrSges(5,1) * qJDD(1) + qJD(1) * t695) * t791 + t756;
t772 = t745 * t649 - t742 * t654;
t644 = m(4) * t670 + t761 * qJDD(1) + (-t715 * t791 + t719 * t747) * qJD(1) + t772;
t651 = t750 * t656 + t748 * t657;
t650 = m(5) * t666 + t708 * mrSges(5,1) + t709 * mrSges(5,2) + t710 * t695 + t711 * t696 + t651;
t718 = t761 * qJD(1);
t647 = m(4) * t669 + t762 * qJDD(1) + (-t715 * t790 - t718 * t747) * qJD(1) - t650;
t773 = t746 * t644 - t743 * t647;
t635 = m(3) * t698 + (qJDD(1) * mrSges(3,3) + qJD(1) * t722) * t747 + t773;
t646 = t742 * t649 + t745 * t654;
t758 = m(4) * t688 + t646;
t763 = t718 * t743 + t719 * t746;
t642 = m(3) * t697 + ((-mrSges(3,3) - t769) * qJDD(1) + (-t722 - t763) * qJD(1)) * t744 - t758;
t774 = t747 * t635 - t744 * t642;
t628 = m(2) * t728 - t752 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t774;
t638 = t743 * t644 + t746 * t647;
t717 = -qJDD(1) * pkin(1) + t759;
t754 = -m(3) * t717 + mrSges(3,1) * t782 - t638 + (t792 + t793) * mrSges(3,3);
t632 = m(2) * t727 - t752 * mrSges(2,2) + (mrSges(2,1) - t796) * qJDD(1) + t754;
t788 = t749 * t628 + t751 * t632;
t630 = t744 * t635 + t747 * t642;
t775 = t751 * t628 - t749 * t632;
t768 = Ifges(3,1) * t744 + Ifges(3,4) * t747;
t767 = Ifges(3,5) * t744 + Ifges(3,6) * t747;
t766 = Ifges(4,5) * t746 - Ifges(4,6) * t743;
t671 = Ifges(6,5) * t693 + Ifges(6,6) * t692 + Ifges(6,3) * t707;
t673 = Ifges(6,1) * t693 + Ifges(6,4) * t692 + Ifges(6,5) * t707;
t652 = -mrSges(6,1) * t660 + mrSges(6,3) * t659 + Ifges(6,4) * t679 + Ifges(6,2) * t678 + Ifges(6,6) * t706 - t693 * t671 + t707 * t673;
t672 = Ifges(6,4) * t693 + Ifges(6,2) * t692 + Ifges(6,6) * t707;
t653 = mrSges(6,2) * t660 - mrSges(6,3) * t658 + Ifges(6,1) * t679 + Ifges(6,4) * t678 + Ifges(6,5) * t706 + t692 * t671 - t707 * t672;
t682 = Ifges(5,5) * t711 - Ifges(5,6) * t710 + Ifges(5,3) * t780;
t683 = Ifges(5,4) * t711 - Ifges(5,2) * t710 + Ifges(5,6) * t780;
t639 = mrSges(5,2) * t666 - mrSges(5,3) * t662 + Ifges(5,1) * t709 - Ifges(5,4) * t708 - pkin(6) * t651 - t748 * t652 + t750 * t653 - t710 * t682 + (Ifges(5,5) * qJDD(1) - qJD(1) * t683) * t791;
t684 = Ifges(5,1) * t711 - Ifges(5,4) * t710 + Ifges(5,5) * t780;
t753 = mrSges(6,1) * t658 - mrSges(6,2) * t659 + Ifges(6,5) * t679 + Ifges(6,6) * t678 + Ifges(6,3) * t706 + t693 * t672 - t692 * t673;
t640 = -mrSges(5,1) * t666 + mrSges(5,3) * t663 + Ifges(5,4) * t709 - Ifges(5,2) * t708 - pkin(4) * t651 - t711 * t682 + (Ifges(5,6) * qJDD(1) + qJD(1) * t684) * t791 - t753;
t701 = (-Ifges(4,3) * t747 + t766 * t744) * qJD(1);
t702 = (-t794 + (-Ifges(4,2) * t743 + t795) * t744) * qJD(1);
t755 = -Ifges(4,5) * t747 + (Ifges(4,1) * t746 - Ifges(4,4) * t743) * t744;
t624 = mrSges(4,2) * t688 - mrSges(4,3) * t669 - qJ(4) * t646 + t745 * t639 - t742 * t640 + (-t701 * t791 + t702 * t747) * qJD(1) + t755 * qJDD(1);
t703 = t755 * qJD(1);
t625 = -mrSges(4,1) * t688 + mrSges(4,3) * t670 - Ifges(5,5) * t709 + Ifges(5,6) * t708 - t711 * t683 - t710 * t684 - mrSges(5,1) * t662 + mrSges(5,2) * t663 - t748 * t653 - t750 * t652 - pkin(4) * t756 - pkin(6) * t771 - pkin(3) * t646 + (-t701 * t790 - t747 * t703) * qJD(1) + (-t794 + (t795 + (-Ifges(4,2) - Ifges(5,3)) * t743) * t744) * qJDD(1);
t723 = t767 * qJD(1);
t621 = mrSges(3,2) * t717 - mrSges(3,3) * t697 - qJ(3) * t638 + t768 * qJDD(1) + t746 * t624 - t743 * t625 + t723 * t785;
t623 = -mrSges(3,1) * t717 + mrSges(3,3) * t698 - mrSges(4,1) * t669 + mrSges(4,2) * t670 - t742 * t639 - t745 * t640 + pkin(3) * t650 - qJ(4) * t772 - pkin(2) * t638 + (Ifges(3,2) + Ifges(4,3)) * t782 + ((Ifges(3,4) - t766) * qJDD(1) + (-t702 * t746 - t703 * t743 - t723) * qJD(1)) * t744;
t637 = mrSges(3,2) * t783 - t754;
t757 = mrSges(2,1) * t727 - mrSges(2,2) * t728 + Ifges(2,3) * qJDD(1) - pkin(1) * t637 + qJ(2) * t774 + t744 * t621 + t747 * t623;
t645 = (t763 * qJD(1) + t769 * qJDD(1)) * t744 + t758;
t619 = mrSges(2,1) * g(3) + mrSges(2,3) * t728 - mrSges(3,1) * t697 + mrSges(3,2) * t698 - t743 * t624 - t746 * t625 + pkin(2) * t645 - qJ(3) * t773 - pkin(1) * t630 + (Ifges(2,6) - t767) * qJDD(1) + (Ifges(2,5) - t744 * (Ifges(3,4) * t744 + Ifges(3,2) * t747) + t747 * t768) * t752;
t618 = -mrSges(2,2) * g(3) - mrSges(2,3) * t727 + Ifges(2,5) * qJDD(1) - t752 * Ifges(2,6) - qJ(2) * t630 + t747 * t621 - t744 * t623;
t1 = [-m(1) * g(1) + t775; -m(1) * g(2) + t788; (-m(1) - m(2)) * g(3) + t630; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t788 + t751 * t618 - t749 * t619; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t775 + t749 * t618 + t751 * t619; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t757; t757; t637; t645; t650; t753;];
tauJB = t1;
