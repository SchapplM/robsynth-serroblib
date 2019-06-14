% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2019-05-05 15:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPPRRR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR6_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR6_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR6_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR6_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR6_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:57:25
% EndTime: 2019-05-05 15:57:30
% DurationCPUTime: 3.47s
% Computational Cost: add. (42442->290), mult. (79958->340), div. (0->0), fcn. (45939->8), ass. (0->115)
t749 = sin(qJ(1));
t753 = cos(qJ(1));
t721 = t749 * g(1) - t753 * g(2);
t755 = qJD(1) ^ 2;
t697 = -qJDD(1) * pkin(1) - t755 * qJ(2) + qJDD(2) - t721;
t688 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) + t697;
t722 = -t753 * g(1) - t749 * g(2);
t783 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t722;
t782 = -m(3) - m(4);
t781 = mrSges(2,1) - mrSges(3,2);
t780 = mrSges(4,3) * t755;
t695 = pkin(1) * t755 - t783;
t689 = qJDD(3) + (-pkin(1) - qJ(3)) * t755 + t783;
t685 = -qJDD(1) * pkin(7) + t689;
t748 = sin(qJ(4));
t752 = cos(qJ(4));
t678 = -g(3) * t752 + t748 * t685;
t714 = (mrSges(5,1) * t748 + mrSges(5,2) * t752) * qJD(1);
t777 = qJD(1) * qJD(4);
t725 = t752 * t777;
t716 = -t748 * qJDD(1) - t725;
t778 = qJD(1) * t752;
t720 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t778;
t727 = t748 * qJD(1);
t684 = -t755 * pkin(7) - t688;
t772 = t748 * t777;
t717 = qJDD(1) * t752 - t772;
t663 = (-t717 + t772) * pkin(8) + (-t716 + t725) * pkin(4) + t684;
t715 = (pkin(4) * t748 - pkin(8) * t752) * qJD(1);
t754 = qJD(4) ^ 2;
t666 = -pkin(4) * t754 + qJDD(4) * pkin(8) - t715 * t727 + t678;
t747 = sin(qJ(5));
t751 = cos(qJ(5));
t648 = t751 * t663 - t747 * t666;
t712 = qJD(4) * t751 - t747 * t778;
t676 = qJD(5) * t712 + qJDD(4) * t747 + t717 * t751;
t711 = qJDD(5) - t716;
t713 = qJD(4) * t747 + t751 * t778;
t724 = t727 + qJD(5);
t645 = (t712 * t724 - t676) * pkin(9) + (t712 * t713 + t711) * pkin(5) + t648;
t649 = t747 * t663 + t751 * t666;
t675 = -qJD(5) * t713 + qJDD(4) * t751 - t717 * t747;
t692 = pkin(5) * t724 - pkin(9) * t713;
t710 = t712 ^ 2;
t646 = -pkin(5) * t710 + pkin(9) * t675 - t692 * t724 + t649;
t746 = sin(qJ(6));
t750 = cos(qJ(6));
t643 = t645 * t750 - t646 * t746;
t679 = t712 * t750 - t713 * t746;
t655 = qJD(6) * t679 + t675 * t746 + t676 * t750;
t680 = t712 * t746 + t713 * t750;
t662 = -mrSges(7,1) * t679 + mrSges(7,2) * t680;
t723 = qJD(6) + t724;
t667 = -mrSges(7,2) * t723 + mrSges(7,3) * t679;
t703 = qJDD(6) + t711;
t639 = m(7) * t643 + mrSges(7,1) * t703 - mrSges(7,3) * t655 - t662 * t680 + t667 * t723;
t644 = t645 * t746 + t646 * t750;
t654 = -qJD(6) * t680 + t675 * t750 - t676 * t746;
t668 = mrSges(7,1) * t723 - mrSges(7,3) * t680;
t640 = m(7) * t644 - mrSges(7,2) * t703 + mrSges(7,3) * t654 + t662 * t679 - t668 * t723;
t631 = t750 * t639 + t746 * t640;
t682 = -mrSges(6,1) * t712 + mrSges(6,2) * t713;
t690 = -mrSges(6,2) * t724 + mrSges(6,3) * t712;
t629 = m(6) * t648 + mrSges(6,1) * t711 - mrSges(6,3) * t676 - t682 * t713 + t690 * t724 + t631;
t691 = mrSges(6,1) * t724 - mrSges(6,3) * t713;
t768 = -t639 * t746 + t750 * t640;
t630 = m(6) * t649 - mrSges(6,2) * t711 + mrSges(6,3) * t675 + t682 * t712 - t691 * t724 + t768;
t769 = -t629 * t747 + t751 * t630;
t623 = m(5) * t678 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t716 - qJD(4) * t720 - t714 * t727 + t769;
t677 = g(3) * t748 + t685 * t752;
t719 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t727;
t665 = -qJDD(4) * pkin(4) - pkin(8) * t754 + t715 * t778 - t677;
t647 = -pkin(5) * t675 - pkin(9) * t710 + t692 * t713 + t665;
t762 = m(7) * t647 - t654 * mrSges(7,1) + mrSges(7,2) * t655 - t679 * t667 + t668 * t680;
t758 = -m(6) * t665 + t675 * mrSges(6,1) - mrSges(6,2) * t676 + t712 * t690 - t691 * t713 - t762;
t635 = m(5) * t677 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t717 + qJD(4) * t719 - t714 * t778 + t758;
t613 = t748 * t623 + t752 * t635;
t767 = m(4) * t689 + qJDD(1) * mrSges(4,2) + t613;
t763 = -m(3) * t695 + t755 * mrSges(3,2) + qJDD(1) * mrSges(3,3) + t767;
t609 = m(2) * t722 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - mrSges(4,3)) * t755 + t763;
t625 = t751 * t629 + t747 * t630;
t764 = -m(5) * t684 + t716 * mrSges(5,1) - t717 * mrSges(5,2) - t719 * t727 - t720 * t778 - t625;
t620 = m(4) * t688 - t755 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t764;
t759 = -m(3) * t697 + t755 * mrSges(3,3) - t620;
t617 = m(2) * t721 - t755 * mrSges(2,2) + t781 * qJDD(1) + t759;
t779 = t749 * t609 + t753 * t617;
t774 = Ifges(2,5) - Ifges(3,4) + Ifges(4,5);
t773 = Ifges(2,6) - Ifges(3,5) - Ifges(4,4);
t771 = t753 * t609 - t617 * t749;
t770 = t752 * t623 - t748 * t635;
t657 = Ifges(7,4) * t680 + Ifges(7,2) * t679 + Ifges(7,6) * t723;
t658 = Ifges(7,1) * t680 + Ifges(7,4) * t679 + Ifges(7,5) * t723;
t761 = -mrSges(7,1) * t643 + mrSges(7,2) * t644 - Ifges(7,5) * t655 - Ifges(7,6) * t654 - Ifges(7,3) * t703 - t680 * t657 + t679 * t658;
t656 = Ifges(7,5) * t680 + Ifges(7,6) * t679 + Ifges(7,3) * t723;
t632 = -mrSges(7,1) * t647 + mrSges(7,3) * t644 + Ifges(7,4) * t655 + Ifges(7,2) * t654 + Ifges(7,6) * t703 - t656 * t680 + t658 * t723;
t633 = mrSges(7,2) * t647 - mrSges(7,3) * t643 + Ifges(7,1) * t655 + Ifges(7,4) * t654 + Ifges(7,5) * t703 + t656 * t679 - t657 * t723;
t669 = Ifges(6,5) * t713 + Ifges(6,6) * t712 + Ifges(6,3) * t724;
t671 = Ifges(6,1) * t713 + Ifges(6,4) * t712 + Ifges(6,5) * t724;
t606 = -mrSges(6,1) * t665 + mrSges(6,3) * t649 + Ifges(6,4) * t676 + Ifges(6,2) * t675 + Ifges(6,6) * t711 - pkin(5) * t762 + pkin(9) * t768 + t750 * t632 + t746 * t633 - t713 * t669 + t724 * t671;
t670 = Ifges(6,4) * t713 + Ifges(6,2) * t712 + Ifges(6,6) * t724;
t615 = mrSges(6,2) * t665 - mrSges(6,3) * t648 + Ifges(6,1) * t676 + Ifges(6,4) * t675 + Ifges(6,5) * t711 - pkin(9) * t631 - t632 * t746 + t633 * t750 + t669 * t712 - t670 * t724;
t701 = (Ifges(5,6) * qJD(4)) + (Ifges(5,4) * t752 - Ifges(5,2) * t748) * qJD(1);
t702 = (Ifges(5,5) * qJD(4)) + (Ifges(5,1) * t752 - Ifges(5,4) * t748) * qJD(1);
t760 = mrSges(5,1) * t677 - mrSges(5,2) * t678 + Ifges(5,5) * t717 + Ifges(5,6) * t716 + Ifges(5,3) * qJDD(4) + pkin(4) * t758 + pkin(8) * t769 + t751 * t606 + t747 * t615 + t701 * t778 + t702 * t727;
t700 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t752 - Ifges(5,6) * t748) * qJD(1);
t603 = mrSges(5,2) * t684 - mrSges(5,3) * t677 + Ifges(5,1) * t717 + Ifges(5,4) * t716 + Ifges(5,5) * qJDD(4) - pkin(8) * t625 - qJD(4) * t701 - t606 * t747 + t615 * t751 - t700 * t727;
t756 = mrSges(6,1) * t648 - mrSges(6,2) * t649 + Ifges(6,5) * t676 + Ifges(6,6) * t675 + Ifges(6,3) * t711 + pkin(5) * t631 + t713 * t670 - t712 * t671 - t761;
t604 = -mrSges(5,1) * t684 + mrSges(5,3) * t678 + Ifges(5,4) * t717 + Ifges(5,2) * t716 + Ifges(5,6) * qJDD(4) - pkin(4) * t625 + qJD(4) * t702 - t700 * t778 - t756;
t619 = qJDD(1) * mrSges(3,2) - t759;
t757 = -mrSges(2,2) * t722 - mrSges(3,3) * t695 - mrSges(4,3) * t688 - pkin(7) * t613 - qJ(3) * t620 + t752 * t603 - t604 * t748 + qJ(2) * (t763 - t780) - pkin(1) * t619 + mrSges(4,2) * t689 + mrSges(3,2) * t697 + mrSges(2,1) * t721 + (Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1);
t612 = t782 * g(3) + t770;
t611 = t767 - t780;
t601 = t760 + (m(4) * qJ(3) + mrSges(4,3) + t781) * g(3) + t774 * t755 - qJ(3) * t770 + mrSges(2,3) * t722 - mrSges(3,1) * t695 + mrSges(4,1) * t689 + pkin(3) * t613 - pkin(1) * t612 + pkin(2) * t611 + t773 * qJDD(1);
t600 = -qJ(2) * t612 - mrSges(2,3) * t721 + pkin(2) * t620 + mrSges(3,1) * t697 + pkin(3) * t764 + pkin(7) * t770 + t748 * t603 + t752 * t604 + mrSges(4,1) * t688 - t773 * t755 + t774 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3) + mrSges(4,2)) * g(3);
t1 = [-m(1) * g(1) + t771; -m(1) * g(2) + t779; (-m(1) - m(2) + t782) * g(3) + t770; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t779 + t753 * t600 - t749 * t601; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t771 + t749 * t600 + t753 * t601; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t757; t757; t619; t611; t760; t756; -t761;];
tauJB  = t1;
