% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPRP3
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
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPRP3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP3_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP3_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP3_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP3_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:51:07
% EndTime: 2019-12-31 19:51:10
% DurationCPUTime: 2.56s
% Computational Cost: add. (35593->226), mult. (49061->270), div. (0->0), fcn. (29563->8), ass. (0->105)
t714 = Ifges(5,1) + Ifges(6,1);
t706 = Ifges(5,4) - Ifges(6,5);
t705 = Ifges(5,5) + Ifges(6,4);
t713 = Ifges(5,2) + Ifges(6,3);
t704 = Ifges(5,6) - Ifges(6,6);
t712 = Ifges(5,3) + Ifges(6,2);
t658 = qJD(1) + qJD(2);
t654 = t658 ^ 2;
t662 = sin(pkin(8));
t663 = cos(pkin(8));
t664 = sin(qJ(4));
t709 = cos(qJ(4));
t711 = t662 * t664 - t663 * t709;
t633 = t711 * t658;
t676 = t709 * t662 + t663 * t664;
t634 = t676 * t658;
t618 = t633 * mrSges(6,1) - t634 * mrSges(6,3);
t655 = qJDD(1) + qJDD(2);
t691 = t633 * qJD(4);
t622 = t676 * t655 - t691;
t666 = sin(qJ(1));
t668 = cos(qJ(1));
t647 = t666 * g(1) - t668 * g(2);
t640 = qJDD(1) * pkin(1) + t647;
t648 = -t668 * g(1) - t666 * g(2);
t670 = qJD(1) ^ 2;
t641 = -t670 * pkin(1) + t648;
t665 = sin(qJ(2));
t667 = cos(qJ(2));
t626 = t665 * t640 + t667 * t641;
t623 = -t654 * pkin(2) + t655 * qJ(3) + t626;
t692 = qJD(3) * t658;
t687 = -t663 * g(3) - 0.2e1 * t662 * t692;
t708 = pkin(3) * t663;
t599 = (-pkin(7) * t655 + t654 * t708 - t623) * t662 + t687;
t603 = -t662 * g(3) + (t623 + 0.2e1 * t692) * t663;
t699 = t655 * t663;
t657 = t663 ^ 2;
t700 = t654 * t657;
t600 = -pkin(3) * t700 + pkin(7) * t699 + t603;
t595 = t709 * t599 - t664 * t600;
t617 = t633 * pkin(4) - t634 * qJ(5);
t669 = qJD(4) ^ 2;
t592 = -qJDD(4) * pkin(4) - t669 * qJ(5) + t634 * t617 + qJDD(5) - t595;
t632 = -t633 * mrSges(6,2) + qJD(4) * mrSges(6,3);
t682 = -m(6) * t592 + qJDD(4) * mrSges(6,1) + qJD(4) * t632;
t589 = t622 * mrSges(6,2) + t634 * t618 - t682;
t596 = t664 * t599 + t709 * t600;
t591 = -t669 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t633 * t617 + t596;
t690 = t634 * qJD(4);
t621 = t711 * t655 + t690;
t631 = -qJD(4) * mrSges(6,1) + t634 * mrSges(6,2);
t689 = m(6) * t591 + qJDD(4) * mrSges(6,3) + qJD(4) * t631;
t694 = t705 * qJD(4) - t706 * t633 + t714 * t634;
t696 = -t704 * qJD(4) + t713 * t633 - t706 * t634;
t710 = t712 * qJDD(4) - t704 * t621 + t705 * t622 + t694 * t633 - t696 * t634 + mrSges(5,1) * t595 - mrSges(6,1) * t592 - mrSges(5,2) * t596 + mrSges(6,3) * t591 - pkin(4) * t589 + qJ(5) * (-t621 * mrSges(6,2) - t633 * t618 + t689);
t707 = -mrSges(5,3) - mrSges(6,2);
t702 = mrSges(4,2) * t662;
t679 = Ifges(4,5) * t662 + Ifges(4,6) * t663;
t701 = t679 * t654;
t630 = qJD(4) * mrSges(5,1) - t634 * mrSges(5,3);
t693 = -t633 * mrSges(5,1) - t634 * mrSges(5,2) - t618;
t585 = m(5) * t596 - qJDD(4) * mrSges(5,2) - qJD(4) * t630 + t707 * t621 + t693 * t633 + t689;
t629 = -qJD(4) * mrSges(5,2) - t633 * mrSges(5,3);
t586 = m(5) * t595 + qJDD(4) * mrSges(5,1) + qJD(4) * t629 + t707 * t622 + t693 * t634 + t682;
t577 = t664 * t585 + t709 * t586;
t602 = -t662 * t623 + t687;
t677 = mrSges(4,3) * t655 + (-mrSges(4,1) * t663 + t702) * t654;
t575 = m(4) * t602 - t677 * t662 + t577;
t683 = t709 * t585 - t664 * t586;
t576 = m(4) * t603 + t677 * t663 + t683;
t684 = -t662 * t575 + t663 * t576;
t567 = m(3) * t626 - t654 * mrSges(3,1) - t655 * mrSges(3,2) + t684;
t625 = t667 * t640 - t665 * t641;
t678 = qJDD(3) - t625;
t620 = -t655 * pkin(2) - t654 * qJ(3) + t678;
t656 = t662 ^ 2;
t601 = (-pkin(2) - t708) * t655 + (-qJ(3) + (-t656 - t657) * pkin(7)) * t654 + t678;
t594 = -0.2e1 * qJD(5) * t634 + (-t622 + t691) * qJ(5) + (t621 + t690) * pkin(4) + t601;
t587 = m(6) * t594 + t621 * mrSges(6,1) - t622 * mrSges(6,3) - t634 * t631 + t633 * t632;
t673 = m(5) * t601 + t621 * mrSges(5,1) + t622 * mrSges(5,2) + t633 * t629 + t634 * t630 + t587;
t671 = -m(4) * t620 + mrSges(4,1) * t699 - t673 + (t654 * t656 + t700) * mrSges(4,3);
t579 = t671 - t654 * mrSges(3,2) + m(3) * t625 + (mrSges(3,1) - t702) * t655;
t564 = t665 * t567 + t667 * t579;
t561 = m(2) * t647 + qJDD(1) * mrSges(2,1) - t670 * mrSges(2,2) + t564;
t685 = t667 * t567 - t665 * t579;
t562 = m(2) * t648 - t670 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t685;
t697 = t668 * t561 + t666 * t562;
t569 = t663 * t575 + t662 * t576;
t695 = -t712 * qJD(4) + t704 * t633 - t705 * t634;
t686 = -t666 * t561 + t668 * t562;
t681 = Ifges(4,1) * t662 + Ifges(4,4) * t663;
t680 = Ifges(4,4) * t662 + Ifges(4,2) * t663;
t570 = -mrSges(5,1) * t601 - mrSges(6,1) * t594 + mrSges(6,2) * t591 + mrSges(5,3) * t596 - pkin(4) * t587 + t694 * qJD(4) + t704 * qJDD(4) - t713 * t621 + t706 * t622 + t695 * t634;
t571 = mrSges(5,2) * t601 + mrSges(6,2) * t592 - mrSges(5,3) * t595 - mrSges(6,3) * t594 - qJ(5) * t587 + t696 * qJD(4) + t705 * qJDD(4) - t706 * t621 + t714 * t622 + t695 * t633;
t555 = -mrSges(4,1) * t620 + mrSges(4,3) * t603 - pkin(3) * t673 + pkin(7) * t683 + t709 * t570 + t664 * t571 + t680 * t655 - t662 * t701;
t557 = mrSges(4,2) * t620 - mrSges(4,3) * t602 - pkin(7) * t577 - t664 * t570 + t709 * t571 + t681 * t655 + t663 * t701;
t581 = t655 * t702 - t671;
t675 = mrSges(3,1) * t625 - mrSges(3,2) * t626 + Ifges(3,3) * t655 - pkin(2) * t581 + qJ(3) * t684 + t663 * t555 + t662 * t557;
t674 = mrSges(2,1) * t647 - mrSges(2,2) * t648 + Ifges(2,3) * qJDD(1) + pkin(1) * t564 + t675;
t553 = mrSges(3,1) * g(3) + (Ifges(3,6) - t679) * t655 + mrSges(3,3) * t626 - mrSges(4,1) * t602 + mrSges(4,2) * t603 - pkin(3) * t577 - pkin(2) * t569 + (-t662 * t680 + t663 * t681 + Ifges(3,5)) * t654 - t710;
t552 = -mrSges(3,2) * g(3) - mrSges(3,3) * t625 + Ifges(3,5) * t655 - t654 * Ifges(3,6) - qJ(3) * t569 - t662 * t555 + t663 * t557;
t551 = -mrSges(2,2) * g(3) - mrSges(2,3) * t647 + Ifges(2,5) * qJDD(1) - t670 * Ifges(2,6) - pkin(6) * t564 + t667 * t552 - t665 * t553;
t550 = Ifges(2,6) * qJDD(1) + t670 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t648 + t665 * t552 + t667 * t553 - pkin(1) * (-m(3) * g(3) + t569) + pkin(6) * t685;
t1 = [-m(1) * g(1) + t686; -m(1) * g(2) + t697; (-m(1) - m(2) - m(3)) * g(3) + t569; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t697 - t666 * t550 + t668 * t551; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t686 + t668 * t550 + t666 * t551; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t674; t674; t675; t581; t710; t589;];
tauJB = t1;
