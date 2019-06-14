% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRPP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
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
% Datum: 2019-05-06 12:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRPP3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:29:30
% EndTime: 2019-05-06 12:29:40
% DurationCPUTime: 6.58s
% Computational Cost: add. (69220->339), mult. (148499->400), div. (0->0), fcn. (100477->8), ass. (0->128)
t746 = Ifges(5,1) + Ifges(6,2) + Ifges(7,3);
t733 = Ifges(5,4) + Ifges(6,6) - Ifges(7,6);
t732 = Ifges(5,5) - Ifges(6,4) + Ifges(7,5);
t745 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t731 = Ifges(5,6) - Ifges(6,5) - Ifges(7,4);
t744 = Ifges(5,3) + Ifges(6,1) + Ifges(7,1);
t743 = -2 * qJD(5);
t742 = 2 * qJD(6);
t741 = cos(qJ(4));
t740 = -mrSges(7,1) - mrSges(5,3);
t704 = sin(pkin(9));
t705 = cos(pkin(9));
t707 = sin(qJ(2));
t736 = qJD(1) * t707;
t686 = qJD(2) * t705 - t704 * t736;
t687 = qJD(2) * t704 + t705 * t736;
t706 = sin(qJ(4));
t658 = -t686 * t741 + t687 * t706;
t709 = cos(qJ(2));
t735 = qJD(1) * t709;
t699 = -qJD(4) + t735;
t739 = t658 * t699;
t708 = sin(qJ(1));
t710 = cos(qJ(1));
t697 = -g(1) * t710 - g(2) * t708;
t712 = qJD(1) ^ 2;
t680 = -pkin(1) * t712 + qJDD(1) * pkin(7) + t697;
t663 = -g(3) * t707 + t709 * t680;
t691 = (-mrSges(3,1) * t709 + mrSges(3,2) * t707) * qJD(1);
t734 = qJD(1) * qJD(2);
t701 = t707 * t734;
t693 = qJDD(1) * t709 - t701;
t694 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t736;
t696 = t708 * g(1) - t710 * g(2);
t679 = -qJDD(1) * pkin(1) - t712 * pkin(7) - t696;
t726 = t709 * t734;
t692 = qJDD(1) * t707 + t726;
t636 = (-t692 - t726) * qJ(3) + (-t693 + t701) * pkin(2) + t679;
t690 = (-pkin(2) * t709 - qJ(3) * t707) * qJD(1);
t711 = qJD(2) ^ 2;
t642 = -pkin(2) * t711 + qJDD(2) * qJ(3) + t690 * t735 + t663;
t603 = -0.2e1 * qJD(3) * t687 + t705 * t636 - t704 * t642;
t668 = qJDD(2) * t704 + t692 * t705;
t599 = (-t686 * t735 - t668) * pkin(8) + (t686 * t687 - t693) * pkin(3) + t603;
t604 = 0.2e1 * qJD(3) * t686 + t704 * t636 + t705 * t642;
t667 = qJDD(2) * t705 - t692 * t704;
t669 = -pkin(3) * t735 - pkin(8) * t687;
t685 = t686 ^ 2;
t602 = -pkin(3) * t685 + pkin(8) * t667 + t669 * t735 + t604;
t596 = t599 * t741 - t706 * t602;
t615 = -t658 * qJD(4) + t706 * t667 + t668 * t741;
t659 = t706 * t686 + t687 * t741;
t630 = -mrSges(7,2) * t659 + mrSges(7,3) * t658;
t632 = mrSges(5,1) * t658 + mrSges(5,2) * t659;
t643 = mrSges(5,2) * t699 - mrSges(5,3) * t658;
t647 = mrSges(6,1) * t658 + mrSges(6,3) * t699;
t689 = qJDD(4) - t693;
t631 = pkin(4) * t658 - qJ(5) * t659;
t698 = t699 ^ 2;
t593 = -t689 * pkin(4) - t698 * qJ(5) + t659 * t631 + qJDD(5) - t596;
t633 = -mrSges(6,2) * t658 - mrSges(6,3) * t659;
t587 = t699 * t742 + (t658 * t659 - t689) * qJ(6) + (t615 - t739) * pkin(5) + t593;
t648 = -mrSges(7,1) * t658 - mrSges(7,2) * t699;
t721 = m(7) * t587 - t689 * mrSges(7,3) + t699 * t648;
t718 = -m(6) * t593 - t615 * mrSges(6,1) - t659 * t633 - t721;
t579 = m(5) * t596 + (-t643 + t647) * t699 + (mrSges(5,1) - mrSges(6,2)) * t689 + (-t630 - t632) * t659 + t740 * t615 + t718;
t597 = t706 * t599 + t741 * t602;
t614 = qJD(4) * t659 - t667 * t741 + t668 * t706;
t644 = -mrSges(5,1) * t699 - mrSges(5,3) * t659;
t717 = -t698 * pkin(4) + t689 * qJ(5) - t658 * t631 + t597;
t592 = 0.2e1 * qJD(5) * t699 - t717;
t649 = mrSges(6,1) * t659 - mrSges(6,2) * t699;
t645 = pkin(5) * t659 + qJ(6) * t699;
t657 = t658 ^ 2;
t589 = -t614 * pkin(5) - t657 * qJ(6) + qJDD(6) + (t743 - t645) * t699 + t717;
t646 = mrSges(7,1) * t659 + mrSges(7,3) * t699;
t727 = -m(7) * t589 - t689 * mrSges(7,2) + t699 * t646;
t719 = -m(6) * t592 + t689 * mrSges(6,3) - t699 * t649 - t727;
t737 = -t630 - t633;
t582 = m(5) * t597 - mrSges(5,2) * t689 + t644 * t699 + (-t632 + t737) * t658 + (-mrSges(6,1) + t740) * t614 + t719;
t577 = t741 * t579 + t706 * t582;
t660 = -mrSges(4,1) * t686 + mrSges(4,2) * t687;
t665 = mrSges(4,2) * t735 + mrSges(4,3) * t686;
t575 = m(4) * t603 - mrSges(4,1) * t693 - mrSges(4,3) * t668 - t660 * t687 - t665 * t735 + t577;
t666 = -mrSges(4,1) * t735 - mrSges(4,3) * t687;
t722 = -t579 * t706 + t741 * t582;
t576 = m(4) * t604 + mrSges(4,2) * t693 + mrSges(4,3) * t667 + t660 * t686 + t666 * t735 + t722;
t723 = -t575 * t704 + t705 * t576;
t570 = m(3) * t663 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t693 - qJD(2) * t694 + t691 * t735 + t723;
t662 = -t709 * g(3) - t707 * t680;
t695 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t735;
t641 = -qJDD(2) * pkin(2) - qJ(3) * t711 + t690 * t736 + qJDD(3) - t662;
t605 = -pkin(3) * t667 - pkin(8) * t685 + t687 * t669 + t641;
t714 = (-t615 - t739) * qJ(5) + t605 + (-t699 * pkin(4) + t743) * t659;
t595 = pkin(4) * t614 + t714;
t591 = t714 + (pkin(4) + qJ(6)) * t614 - pkin(5) * t657 + t658 * t742 - t645 * t659;
t720 = m(7) * t591 - t615 * mrSges(7,2) + t614 * mrSges(7,3) - t659 * t646 + t658 * t648;
t585 = m(6) * t595 - t614 * mrSges(6,2) - t615 * mrSges(6,3) - t658 * t647 - t659 * t649 + t720;
t715 = m(5) * t605 + t614 * mrSges(5,1) + t615 * mrSges(5,2) + t658 * t643 + t659 * t644 + t585;
t713 = -m(4) * t641 + t667 * mrSges(4,1) - t668 * mrSges(4,2) + t686 * t665 - t687 * t666 - t715;
t584 = m(3) * t662 + qJDD(2) * mrSges(3,1) - t692 * mrSges(3,3) + qJD(2) * t695 - t691 * t736 + t713;
t724 = t709 * t570 - t584 * t707;
t564 = m(2) * t697 - mrSges(2,1) * t712 - qJDD(1) * mrSges(2,2) + t724;
t571 = t575 * t705 + t576 * t704;
t716 = -m(3) * t679 + t693 * mrSges(3,1) - mrSges(3,2) * t692 - t694 * t736 + t695 * t735 - t571;
t567 = m(2) * t696 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t712 + t716;
t738 = t708 * t564 + t710 * t567;
t565 = t707 * t570 + t709 * t584;
t730 = t658 * t731 - t659 * t732 + t699 * t744;
t729 = t658 * t745 + t659 * t733 - t699 * t731;
t728 = t733 * t658 - t659 * t746 + t732 * t699;
t725 = t710 * t564 - t567 * t708;
t678 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t707 + Ifges(3,4) * t709) * qJD(1);
t677 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t707 + Ifges(3,2) * t709) * qJD(1);
t676 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t707 + Ifges(3,6) * t709) * qJD(1);
t654 = Ifges(4,1) * t687 + Ifges(4,4) * t686 - Ifges(4,5) * t735;
t653 = Ifges(4,4) * t687 + Ifges(4,2) * t686 - Ifges(4,6) * t735;
t652 = Ifges(4,5) * t687 + Ifges(4,6) * t686 - Ifges(4,3) * t735;
t586 = mrSges(7,1) * t615 + t630 * t659 + t721;
t573 = mrSges(6,1) * t593 + mrSges(7,1) * t587 + mrSges(5,2) * t605 - mrSges(7,2) * t591 - mrSges(5,3) * t596 - mrSges(6,3) * t595 + pkin(5) * t586 - qJ(5) * t585 - t733 * t614 + t615 * t746 + t730 * t658 + t732 * t689 + t729 * t699;
t572 = -mrSges(5,1) * t605 + mrSges(5,3) * t597 - mrSges(6,1) * t592 + mrSges(6,2) * t595 + mrSges(7,1) * t589 - mrSges(7,3) * t591 - pkin(5) * (t658 * t630 + t727) - qJ(6) * t720 - pkin(4) * t585 + t728 * t699 + t731 * t689 + t730 * t659 + t733 * t615 + (-mrSges(7,1) * pkin(5) + t745) * t614;
t561 = mrSges(4,2) * t641 - mrSges(4,3) * t603 + Ifges(4,1) * t668 + Ifges(4,4) * t667 - Ifges(4,5) * t693 - pkin(8) * t577 - t706 * t572 + t573 * t741 + t686 * t652 + t653 * t735;
t560 = -mrSges(4,1) * t641 + mrSges(4,3) * t604 + Ifges(4,4) * t668 + Ifges(4,2) * t667 - Ifges(4,6) * t693 - pkin(3) * t715 + pkin(8) * t722 + t572 * t741 + t706 * t573 - t687 * t652 - t654 * t735;
t559 = -mrSges(5,1) * t596 + mrSges(5,2) * t597 - pkin(2) * t571 + mrSges(3,3) * t663 - pkin(3) * t577 - Ifges(4,6) * t667 - Ifges(4,5) * t668 + qJD(2) * t678 - mrSges(3,1) * t679 + mrSges(6,3) * t592 - mrSges(6,2) * t593 + qJ(6) * t586 - mrSges(4,1) * t603 + mrSges(4,2) * t604 + t686 * t654 - t687 * t653 + Ifges(3,4) * t692 + (-qJ(5) * t737 + t728) * t658 - t676 * t736 + (pkin(4) * t630 - t729) * t659 + (-qJ(5) * (-mrSges(6,1) - mrSges(7,1)) + t731) * t614 + (mrSges(7,1) * pkin(4) - t732) * t615 + (mrSges(6,2) * pkin(4) - t744) * t689 - pkin(4) * (t647 * t699 + t718) - qJ(5) * t719 + Ifges(3,6) * qJDD(2) + mrSges(7,3) * t587 - mrSges(7,2) * t589 + (Ifges(3,2) + Ifges(4,3)) * t693;
t558 = mrSges(3,2) * t679 - mrSges(3,3) * t662 + Ifges(3,1) * t692 + Ifges(3,4) * t693 + Ifges(3,5) * qJDD(2) - qJ(3) * t571 - qJD(2) * t677 - t560 * t704 + t561 * t705 + t676 * t735;
t557 = Ifges(2,6) * qJDD(1) + t712 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t697 - Ifges(3,5) * t692 - Ifges(3,6) * t693 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t662 + mrSges(3,2) * t663 - t704 * t561 - t705 * t560 - pkin(2) * t713 - qJ(3) * t723 - pkin(1) * t565 + (-t677 * t707 + t678 * t709) * qJD(1);
t556 = -mrSges(2,2) * g(3) - mrSges(2,3) * t696 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t712 - pkin(7) * t565 + t558 * t709 - t559 * t707;
t1 = [-m(1) * g(1) + t725; -m(1) * g(2) + t738; (-m(1) - m(2)) * g(3) + t565; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t738 + t710 * t556 - t708 * t557; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t725 + t708 * t556 + t710 * t557; -mrSges(1,1) * g(2) + mrSges(2,1) * t696 + mrSges(1,2) * g(1) - mrSges(2,2) * t697 + Ifges(2,3) * qJDD(1) + pkin(1) * t716 + pkin(7) * t724 + t707 * t558 + t709 * t559;];
tauB  = t1;
