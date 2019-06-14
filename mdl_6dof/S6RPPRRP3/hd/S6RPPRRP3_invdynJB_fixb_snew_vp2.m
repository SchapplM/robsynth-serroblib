% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
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
% Datum: 2019-05-05 14:52
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPPRRP3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:51:00
% EndTime: 2019-05-05 14:51:04
% DurationCPUTime: 2.99s
% Computational Cost: add. (31742->270), mult. (57582->314), div. (0->0), fcn. (31255->8), ass. (0->114)
t789 = Ifges(6,1) + Ifges(7,1);
t779 = Ifges(6,4) - Ifges(7,5);
t777 = -Ifges(6,5) - Ifges(7,4);
t788 = Ifges(6,2) + Ifges(7,3);
t776 = Ifges(6,6) - Ifges(7,6);
t787 = -Ifges(6,3) - Ifges(7,2);
t744 = sin(qJ(1));
t746 = cos(qJ(1));
t719 = t744 * g(1) - t746 * g(2);
t710 = qJDD(1) * pkin(1) + t719;
t720 = -t746 * g(1) - t744 * g(2);
t748 = qJD(1) ^ 2;
t712 = -t748 * pkin(1) + t720;
t740 = sin(pkin(9));
t741 = cos(pkin(9));
t681 = t740 * t710 + t741 * t712;
t786 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) - t681;
t742 = sin(qJ(5));
t745 = cos(qJ(4));
t769 = qJD(1) * t745;
t783 = cos(qJ(5));
t708 = -t783 * qJD(4) + t742 * t769;
t743 = sin(qJ(4));
t767 = qJD(1) * qJD(4);
t763 = t743 * t767;
t715 = t745 * qJDD(1) - t763;
t678 = -t708 * qJD(5) + t742 * qJDD(4) + t783 * t715;
t709 = t742 * qJD(4) + t783 * t769;
t685 = t708 * mrSges(7,1) - t709 * mrSges(7,3);
t784 = -pkin(2) - pkin(7);
t662 = t784 * t748 - t786;
t762 = t745 * t767;
t714 = -t743 * qJDD(1) - t762;
t653 = (-t715 + t763) * pkin(8) + (-t714 + t762) * pkin(4) + t662;
t680 = t741 * t710 - t740 * t712;
t756 = -t748 * qJ(3) + qJDD(3) - t680;
t663 = t784 * qJDD(1) + t756;
t737 = -g(3) + qJDD(2);
t659 = t743 * t663 + t745 * t737;
t713 = (pkin(4) * t743 - pkin(8) * t745) * qJD(1);
t747 = qJD(4) ^ 2;
t768 = t743 * qJD(1);
t656 = -t747 * pkin(4) + qJDD(4) * pkin(8) - t713 * t768 + t659;
t650 = t783 * t653 - t742 * t656;
t684 = t708 * pkin(5) - t709 * qJ(6);
t707 = qJDD(5) - t714;
t722 = qJD(5) + t768;
t721 = t722 ^ 2;
t648 = -t707 * pkin(5) - t721 * qJ(6) + t709 * t684 + qJDD(6) - t650;
t690 = -t708 * mrSges(7,2) + t722 * mrSges(7,3);
t757 = -m(7) * t648 + t707 * mrSges(7,1) + t722 * t690;
t644 = t678 * mrSges(7,2) + t709 * t685 - t757;
t651 = t742 * t653 + t783 * t656;
t647 = -t721 * pkin(5) + t707 * qJ(6) + 0.2e1 * qJD(6) * t722 - t708 * t684 + t651;
t677 = t709 * qJD(5) - t783 * qJDD(4) + t742 * t715;
t689 = -t722 * mrSges(7,1) + t709 * mrSges(7,2);
t765 = m(7) * t647 + t707 * mrSges(7,3) + t722 * t689;
t771 = t779 * t708 - t789 * t709 + t777 * t722;
t772 = t788 * t708 - t779 * t709 - t776 * t722;
t785 = -t776 * t677 - t777 * t678 - t787 * t707 - t771 * t708 - t772 * t709 + mrSges(6,1) * t650 - mrSges(7,1) * t648 - mrSges(6,2) * t651 + mrSges(7,3) * t647 - pkin(5) * t644 + qJ(6) * (-t677 * mrSges(7,2) - t708 * t685 + t765);
t782 = mrSges(3,1) - mrSges(4,2);
t781 = -mrSges(6,3) - mrSges(7,2);
t780 = -Ifges(4,4) + Ifges(3,5);
t778 = Ifges(4,5) - Ifges(3,6);
t711 = (mrSges(5,1) * t743 + mrSges(5,2) * t745) * qJD(1);
t718 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t769;
t688 = t722 * mrSges(6,1) - t709 * mrSges(6,3);
t770 = -t708 * mrSges(6,1) - t709 * mrSges(6,2) - t685;
t639 = m(6) * t651 - t707 * mrSges(6,2) + t781 * t677 - t722 * t688 + t770 * t708 + t765;
t687 = -t722 * mrSges(6,2) - t708 * mrSges(6,3);
t641 = m(6) * t650 + t707 * mrSges(6,1) + t781 * t678 + t722 * t687 + t770 * t709 + t757;
t758 = t783 * t639 - t742 * t641;
t630 = m(5) * t659 - qJDD(4) * mrSges(5,2) + t714 * mrSges(5,3) - qJD(4) * t718 - t711 * t768 + t758;
t658 = t745 * t663 - t743 * t737;
t717 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t768;
t655 = -qJDD(4) * pkin(4) - t747 * pkin(8) + t713 * t769 - t658;
t649 = -0.2e1 * qJD(6) * t709 + (t708 * t722 - t678) * qJ(6) + (t709 * t722 + t677) * pkin(5) + t655;
t645 = m(7) * t649 + t677 * mrSges(7,1) - t678 * mrSges(7,3) - t709 * t689 + t708 * t690;
t750 = -m(6) * t655 - t677 * mrSges(6,1) - t678 * mrSges(6,2) - t708 * t687 - t709 * t688 - t645;
t636 = m(5) * t658 + qJDD(4) * mrSges(5,1) - t715 * mrSges(5,3) + qJD(4) * t717 - t711 * t769 + t750;
t623 = t743 * t630 + t745 * t636;
t666 = -qJDD(1) * pkin(2) + t756;
t755 = -m(4) * t666 + t748 * mrSges(4,3) - t623;
t618 = m(3) * t680 - t748 * mrSges(3,2) + t782 * qJDD(1) + t755;
t664 = t748 * pkin(2) + t786;
t635 = t742 * t639 + t783 * t641;
t754 = -m(5) * t662 + t714 * mrSges(5,1) - t715 * mrSges(5,2) - t717 * t768 - t718 * t769 - t635;
t752 = -m(4) * t664 + t748 * mrSges(4,2) + qJDD(1) * mrSges(4,3) - t754;
t626 = m(3) * t681 - t748 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t752;
t615 = t741 * t618 + t740 * t626;
t612 = m(2) * t719 + qJDD(1) * mrSges(2,1) - t748 * mrSges(2,2) + t615;
t760 = -t740 * t618 + t741 * t626;
t613 = m(2) * t720 - t748 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t760;
t774 = t746 * t612 + t744 * t613;
t773 = t776 * t708 + t777 * t709 + t787 * t722;
t761 = -t744 * t612 + t746 * t613;
t759 = t745 * t630 - t743 * t636;
t622 = m(4) * t737 + t759;
t621 = m(3) * t737 + t622;
t631 = -mrSges(6,1) * t655 - mrSges(7,1) * t649 + mrSges(7,2) * t647 + mrSges(6,3) * t651 - pkin(5) * t645 - t788 * t677 + t779 * t678 + t776 * t707 + t773 * t709 - t771 * t722;
t633 = mrSges(6,2) * t655 + mrSges(7,2) * t648 - mrSges(6,3) * t650 - mrSges(7,3) * t649 - qJ(6) * t645 - t779 * t677 + t789 * t678 - t777 * t707 + t773 * t708 + t772 * t722;
t698 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t745 - Ifges(5,2) * t743) * qJD(1);
t699 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t745 - Ifges(5,4) * t743) * qJD(1);
t753 = mrSges(5,1) * t658 - mrSges(5,2) * t659 + Ifges(5,5) * t715 + Ifges(5,6) * t714 + Ifges(5,3) * qJDD(4) + pkin(4) * t750 + pkin(8) * t758 + t783 * t631 + t742 * t633 + t698 * t769 + t699 * t768;
t697 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t745 - Ifges(5,6) * t743) * qJD(1);
t608 = mrSges(5,2) * t662 - mrSges(5,3) * t658 + Ifges(5,1) * t715 + Ifges(5,4) * t714 + Ifges(5,5) * qJDD(4) - pkin(8) * t635 - qJD(4) * t698 - t742 * t631 + t783 * t633 - t697 * t768;
t616 = -mrSges(5,1) * t662 + mrSges(5,3) * t659 + Ifges(5,4) * t715 + Ifges(5,2) * t714 + Ifges(5,6) * qJDD(4) - pkin(4) * t635 + qJD(4) * t699 - t697 * t769 - t785;
t620 = qJDD(1) * mrSges(4,2) - t755;
t749 = mrSges(2,1) * t719 + mrSges(3,1) * t680 - mrSges(2,2) * t720 - mrSges(3,2) * t681 + mrSges(4,2) * t666 - mrSges(4,3) * t664 + pkin(1) * t615 - pkin(2) * t620 - pkin(7) * t623 + qJ(3) * t752 + t745 * t608 - t743 * t616 + (Ifges(2,3) + Ifges(3,3) + Ifges(4,1)) * qJDD(1);
t606 = mrSges(4,1) * t666 + t753 - mrSges(3,3) * t680 + t778 * t748 + (mrSges(3,2) - mrSges(4,3)) * t737 - qJ(3) * t622 + pkin(3) * t623 + t780 * qJDD(1);
t605 = -mrSges(4,1) * t664 + mrSges(3,3) * t681 - pkin(2) * t622 - pkin(3) * t754 - pkin(7) * t759 - t778 * qJDD(1) - t743 * t608 - t745 * t616 - t782 * t737 + t780 * t748;
t604 = -mrSges(2,2) * g(3) - mrSges(2,3) * t719 + Ifges(2,5) * qJDD(1) - t748 * Ifges(2,6) - qJ(2) * t615 - t740 * t605 + t741 * t606;
t603 = mrSges(2,1) * g(3) + mrSges(2,3) * t720 + t748 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t621 + qJ(2) * t760 + t741 * t605 + t740 * t606;
t1 = [-m(1) * g(1) + t761; -m(1) * g(2) + t774; (-m(1) - m(2)) * g(3) + t621; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t774 - t744 * t603 + t746 * t604; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t761 + t746 * t603 + t744 * t604; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t749; t749; t621; t620; t753; t785; t644;];
tauJB  = t1;
