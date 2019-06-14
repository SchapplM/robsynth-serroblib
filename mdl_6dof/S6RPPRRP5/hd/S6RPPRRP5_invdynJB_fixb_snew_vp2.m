% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPPRRP5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2019-05-05 14:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPPRRP5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:57:35
% EndTime: 2019-05-05 14:57:38
% DurationCPUTime: 1.89s
% Computational Cost: add. (18066->267), mult. (33405->298), div. (0->0), fcn. (17027->6), ass. (0->106)
t783 = Ifges(6,1) + Ifges(7,1);
t775 = Ifges(6,4) + Ifges(7,4);
t774 = Ifges(6,5) + Ifges(7,5);
t782 = Ifges(6,2) + Ifges(7,2);
t773 = Ifges(6,6) + Ifges(7,6);
t781 = Ifges(6,3) + Ifges(7,3);
t735 = sin(qJ(1));
t738 = cos(qJ(1));
t712 = t735 * g(1) - t738 * g(2);
t740 = qJD(1) ^ 2;
t689 = -qJDD(1) * pkin(1) - t740 * qJ(2) + qJDD(2) - t712;
t678 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) + t689;
t713 = -t738 * g(1) - t735 * g(2);
t780 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t713;
t733 = sin(qJ(5));
t736 = cos(qJ(5));
t737 = cos(qJ(4));
t764 = qJD(1) * t737;
t703 = qJD(4) * t736 - t733 * t764;
t734 = sin(qJ(4));
t763 = qJD(1) * qJD(4);
t756 = t734 * t763;
t708 = qJDD(1) * t737 - t756;
t666 = t703 * qJD(5) + qJDD(4) * t733 + t708 * t736;
t704 = qJD(4) * t733 + t736 * t764;
t671 = -mrSges(7,1) * t703 + mrSges(7,2) * t704;
t674 = -t740 * pkin(7) - t678;
t755 = t737 * t763;
t707 = -qJDD(1) * t734 - t755;
t648 = (-t708 + t756) * pkin(8) + (-t707 + t755) * pkin(4) + t674;
t679 = qJDD(3) + (-pkin(1) - qJ(3)) * t740 + t780;
t675 = -qJDD(1) * pkin(7) + t679;
t668 = -g(3) * t737 + t734 * t675;
t706 = (pkin(4) * t734 - pkin(8) * t737) * qJD(1);
t739 = qJD(4) ^ 2;
t765 = qJD(1) * t734;
t651 = -pkin(4) * t739 + qJDD(4) * pkin(8) - t706 * t765 + t668;
t644 = t736 * t648 - t733 * t651;
t702 = qJDD(5) - t707;
t714 = qJD(5) + t765;
t640 = -0.2e1 * qJD(6) * t704 + (t703 * t714 - t666) * qJ(6) + (t703 * t704 + t702) * pkin(5) + t644;
t680 = -mrSges(7,2) * t714 + t703 * mrSges(7,3);
t758 = m(7) * t640 + t702 * mrSges(7,1) + t714 * t680;
t637 = -t666 * mrSges(7,3) - t704 * t671 + t758;
t645 = t733 * t648 + t736 * t651;
t665 = -t704 * qJD(5) + qJDD(4) * t736 - t708 * t733;
t682 = pkin(5) * t714 - t704 * qJ(6);
t701 = t703 ^ 2;
t642 = -t701 * pkin(5) + t665 * qJ(6) + 0.2e1 * qJD(6) * t703 - t682 * t714 + t645;
t767 = t775 * t703 + t783 * t704 + t774 * t714;
t768 = -t782 * t703 - t775 * t704 - t773 * t714;
t779 = mrSges(6,1) * t644 + mrSges(7,1) * t640 - mrSges(6,2) * t645 - mrSges(7,2) * t642 + pkin(5) * t637 + t773 * t665 + t774 * t666 + t781 * t702 - t767 * t703 - t768 * t704;
t778 = -m(3) - m(4);
t777 = mrSges(2,1) - mrSges(3,2);
t776 = -mrSges(6,2) - mrSges(7,2);
t771 = mrSges(4,3) * t740;
t687 = pkin(1) * t740 - t780;
t705 = (mrSges(5,1) * t734 + mrSges(5,2) * t737) * qJD(1);
t711 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t764;
t672 = -mrSges(6,1) * t703 + mrSges(6,2) * t704;
t681 = -mrSges(6,2) * t714 + t703 * mrSges(6,3);
t630 = m(6) * t644 + t702 * mrSges(6,1) + t714 * t681 + (-t671 - t672) * t704 + (-mrSges(6,3) - mrSges(7,3)) * t666 + t758;
t757 = m(7) * t642 + t665 * mrSges(7,3) + t703 * t671;
t683 = mrSges(7,1) * t714 - t704 * mrSges(7,3);
t766 = -mrSges(6,1) * t714 + t704 * mrSges(6,3) - t683;
t633 = m(6) * t645 + t665 * mrSges(6,3) + t703 * t672 + t776 * t702 + t766 * t714 + t757;
t752 = -t630 * t733 + t736 * t633;
t625 = m(5) * t668 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t707 - qJD(4) * t711 - t705 * t765 + t752;
t667 = g(3) * t734 + t675 * t737;
t710 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t765;
t650 = -qJDD(4) * pkin(4) - pkin(8) * t739 + t706 * t764 - t667;
t643 = -t665 * pkin(5) - t701 * qJ(6) + t704 * t682 + qJDD(6) + t650;
t750 = -m(7) * t643 + t665 * mrSges(7,1) + t703 * t680;
t742 = -m(6) * t650 + t665 * mrSges(6,1) + t776 * t666 + t703 * t681 + t766 * t704 + t750;
t635 = m(5) * t667 + qJDD(4) * mrSges(5,1) - t708 * mrSges(5,3) + qJD(4) * t710 - t705 * t764 + t742;
t614 = t734 * t625 + t737 * t635;
t751 = m(4) * t679 + qJDD(1) * mrSges(4,2) + t614;
t746 = -m(3) * t687 + t740 * mrSges(3,2) + qJDD(1) * mrSges(3,3) + t751;
t610 = m(2) * t713 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - mrSges(4,3)) * t740 + t746;
t628 = t736 * t630 + t733 * t633;
t747 = -m(5) * t674 + t707 * mrSges(5,1) - t708 * mrSges(5,2) - t710 * t765 - t711 * t764 - t628;
t621 = m(4) * t678 - t740 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t747;
t743 = -m(3) * t689 + t740 * mrSges(3,3) - t621;
t618 = m(2) * t712 - t740 * mrSges(2,2) + t777 * qJDD(1) + t743;
t770 = t735 * t610 + t738 * t618;
t769 = -t773 * t703 - t774 * t704 - t781 * t714;
t760 = Ifges(2,5) - Ifges(3,4) + Ifges(4,5);
t759 = Ifges(2,6) - Ifges(3,5) - Ifges(4,4);
t754 = t738 * t610 - t618 * t735;
t753 = t737 * t625 - t734 * t635;
t638 = t666 * mrSges(7,2) + t704 * t683 - t750;
t616 = -mrSges(6,1) * t650 + mrSges(6,3) * t645 - mrSges(7,1) * t643 + mrSges(7,3) * t642 - pkin(5) * t638 + qJ(6) * t757 + (-qJ(6) * t683 + t767) * t714 + t769 * t704 + (-mrSges(7,2) * qJ(6) + t773) * t702 + t775 * t666 + t782 * t665;
t626 = mrSges(6,2) * t650 + mrSges(7,2) * t643 - mrSges(6,3) * t644 - mrSges(7,3) * t640 - qJ(6) * t637 + t775 * t665 + t783 * t666 + t774 * t702 - t769 * t703 + t768 * t714;
t692 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t737 - Ifges(5,2) * t734) * qJD(1);
t693 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t737 - Ifges(5,4) * t734) * qJD(1);
t744 = mrSges(5,1) * t667 - mrSges(5,2) * t668 + Ifges(5,5) * t708 + Ifges(5,6) * t707 + Ifges(5,3) * qJDD(4) + pkin(4) * t742 + pkin(8) * t752 + t736 * t616 + t733 * t626 + t692 * t764 + t693 * t765;
t691 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t737 - Ifges(5,6) * t734) * qJD(1);
t606 = mrSges(5,2) * t674 - mrSges(5,3) * t667 + Ifges(5,1) * t708 + Ifges(5,4) * t707 + Ifges(5,5) * qJDD(4) - pkin(8) * t628 - qJD(4) * t692 - t616 * t733 + t626 * t736 - t691 * t765;
t607 = -mrSges(5,1) * t674 + mrSges(5,3) * t668 + Ifges(5,4) * t708 + Ifges(5,2) * t707 + Ifges(5,6) * qJDD(4) - pkin(4) * t628 + qJD(4) * t693 - t691 * t764 - t779;
t620 = qJDD(1) * mrSges(3,2) - t743;
t741 = -mrSges(2,2) * t713 - mrSges(3,3) * t687 - mrSges(4,3) * t678 - pkin(7) * t614 - qJ(3) * t621 + t737 * t606 - t607 * t734 + qJ(2) * (t746 - t771) - pkin(1) * t620 + mrSges(4,2) * t679 + mrSges(3,2) * t689 + mrSges(2,1) * t712 + (Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1);
t613 = t778 * g(3) + t753;
t612 = t751 - t771;
t604 = t744 + t760 * t740 - qJ(3) * t753 + t759 * qJDD(1) + (m(4) * qJ(3) + mrSges(4,3) + t777) * g(3) + mrSges(2,3) * t713 - mrSges(3,1) * t687 + mrSges(4,1) * t679 + pkin(3) * t614 + pkin(2) * t612 - pkin(1) * t613;
t603 = -qJ(2) * t613 - mrSges(2,3) * t712 + pkin(2) * t621 + mrSges(3,1) * t689 + t734 * t606 + t737 * t607 + pkin(3) * t747 + pkin(7) * t753 + mrSges(4,1) * t678 - t759 * t740 + t760 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3) + mrSges(4,2)) * g(3);
t1 = [-m(1) * g(1) + t754; -m(1) * g(2) + t770; (-m(1) - m(2) + t778) * g(3) + t753; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t770 + t738 * t603 - t735 * t604; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t754 + t735 * t603 + t738 * t604; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t741; t741; t620; t612; t744; t779; t638;];
tauJB  = t1;
