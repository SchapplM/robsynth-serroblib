% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRPR16
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
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
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRPR16_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR16_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR16_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR16_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR16_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR16_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR16_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:38:54
% EndTime: 2019-12-31 18:38:56
% DurationCPUTime: 1.48s
% Computational Cost: add. (10350->256), mult. (20172->298), div. (0->0), fcn. (9414->6), ass. (0->110)
t717 = Ifges(4,5) - Ifges(5,4);
t716 = Ifges(4,6) - Ifges(5,5);
t715 = (Ifges(4,3) + Ifges(5,1));
t667 = sin(qJ(1));
t670 = cos(qJ(1));
t647 = -t670 * g(1) - t667 * g(2);
t686 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t647;
t672 = qJD(1) ^ 2;
t711 = (-pkin(1) - pkin(6));
t696 = t711 * t672;
t607 = t696 + t686;
t666 = sin(qJ(3));
t669 = cos(qJ(3));
t698 = qJD(1) * qJD(3);
t695 = t669 * t698;
t637 = t666 * qJDD(1) + t695;
t652 = t666 * t698;
t638 = t669 * qJDD(1) - t652;
t700 = qJD(1) * t666;
t641 = -(qJD(3) * mrSges(4,2)) - mrSges(4,3) * t700;
t699 = t669 * qJD(1);
t642 = (qJD(3) * mrSges(4,1)) - mrSges(4,3) * t699;
t712 = -2 * qJD(4);
t678 = pkin(3) * t695 + t699 * t712 + t686 + (-t638 + t652) * qJ(4);
t585 = t637 * pkin(3) + t678 + t696;
t643 = mrSges(5,1) * t700 - (qJD(3) * mrSges(5,3));
t644 = mrSges(5,1) * t699 + (qJD(3) * mrSges(5,2));
t645 = pkin(4) * t699 - (qJD(3) * pkin(7));
t662 = t666 ^ 2;
t710 = pkin(3) + pkin(7);
t580 = -t645 * t699 + t710 * t637 + (-pkin(4) * t662 + t711) * t672 + t678;
t634 = (pkin(3) * t666 - qJ(4) * t669) * qJD(1);
t671 = qJD(3) ^ 2;
t646 = t667 * g(1) - t670 * g(2);
t685 = -t672 * qJ(2) + qJDD(2) - t646;
t608 = t711 * qJDD(1) + t685;
t702 = t669 * t608;
t684 = -t671 * qJ(4) + t634 * t699 + qJDD(4) - t702;
t709 = pkin(7) * t672;
t583 = t638 * pkin(4) - t710 * qJDD(3) + (pkin(4) * t698 + t669 * t709 - g(3)) * t666 + t684;
t665 = sin(qJ(5));
t668 = cos(qJ(5));
t578 = -t665 * t580 + t668 * t583;
t632 = -t665 * qJD(3) + t668 * t700;
t598 = t632 * qJD(5) + t668 * qJDD(3) + t665 * t637;
t633 = t668 * qJD(3) + t665 * t700;
t600 = -t632 * mrSges(6,1) + t633 * mrSges(6,2);
t650 = qJD(5) + t699;
t603 = -t650 * mrSges(6,2) + t632 * mrSges(6,3);
t631 = qJDD(5) + t638;
t575 = m(6) * t578 + t631 * mrSges(6,1) - t598 * mrSges(6,3) - t633 * t600 + t650 * t603;
t579 = t668 * t580 + t665 * t583;
t597 = -t633 * qJD(5) - t665 * qJDD(3) + t668 * t637;
t604 = t650 * mrSges(6,1) - t633 * mrSges(6,3);
t576 = m(6) * t579 - t631 * mrSges(6,2) + t597 * mrSges(6,3) + t632 * t600 - t650 * t604;
t692 = -t665 * t575 + t668 * t576;
t676 = m(5) * t585 - t638 * mrSges(5,3) - (t643 * t666 + t644 * t669) * qJD(1) + t692;
t706 = mrSges(4,1) - mrSges(5,2);
t714 = -m(4) * t607 - t638 * mrSges(4,2) - t706 * t637 - t641 * t700 - t642 * t699 - t676;
t708 = t666 * g(3);
t707 = mrSges(2,1) - mrSges(3,2);
t705 = Ifges(4,4) + Ifges(5,6);
t704 = Ifges(2,5) - Ifges(3,4);
t703 = (-Ifges(2,6) + Ifges(3,5));
t601 = t702 + t708;
t567 = t668 * t575 + t665 * t576;
t588 = -qJDD(3) * pkin(3) + t684 - t708;
t681 = -m(5) * t588 - t638 * mrSges(5,1) - t567;
t635 = (-mrSges(5,2) * t666 - mrSges(5,3) * t669) * qJD(1);
t690 = qJD(1) * (-t635 - (mrSges(4,1) * t666 + mrSges(4,2) * t669) * qJD(1));
t563 = m(4) * t601 - t638 * mrSges(4,3) + t706 * qJDD(3) + (t641 - t643) * qJD(3) + t669 * t690 + t681;
t602 = -t669 * g(3) + t666 * t608;
t679 = -t671 * pkin(3) + qJDD(3) * qJ(4) - t634 * t700 + t602;
t586 = (qJD(3) * t712) - t679;
t582 = -t662 * t709 - t637 * pkin(4) + ((2 * qJD(4)) + t645) * qJD(3) + t679;
t683 = -m(6) * t582 + t597 * mrSges(6,1) - t598 * mrSges(6,2) + t632 * t603 - t633 * t604;
t677 = -m(5) * t586 + qJDD(3) * mrSges(5,3) + qJD(3) * t644 - t683;
t572 = m(4) * t602 - qJDD(3) * mrSges(4,2) - qJD(3) * t642 + (-mrSges(4,3) - mrSges(5,1)) * t637 + t666 * t690 + t677;
t558 = t669 * t563 + t666 * t572;
t613 = -qJDD(1) * pkin(1) + t685;
t682 = -m(3) * t613 + (t672 * mrSges(3,3)) - t558;
t554 = m(2) * t646 - (t672 * mrSges(2,2)) + t707 * qJDD(1) + t682;
t611 = t672 * pkin(1) - t686;
t674 = -m(3) * t611 + (t672 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t714;
t561 = m(2) * t647 - (t672 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t674;
t701 = t670 * t554 + t667 * t561;
t694 = -t667 * t554 + t670 * t561;
t693 = -t666 * t563 + t669 * t572;
t691 = qJD(1) * (-(t715 * qJD(3)) + (t716 * t666 - t717 * t669) * qJD(1));
t591 = Ifges(6,4) * t633 + Ifges(6,2) * t632 + Ifges(6,6) * t650;
t592 = Ifges(6,1) * t633 + Ifges(6,4) * t632 + Ifges(6,5) * t650;
t680 = mrSges(6,1) * t578 - mrSges(6,2) * t579 + Ifges(6,5) * t598 + Ifges(6,6) * t597 + Ifges(6,3) * t631 + t633 * t591 - t632 * t592;
t564 = -t637 * mrSges(5,2) + t676;
t590 = Ifges(6,5) * t633 + Ifges(6,6) * t632 + Ifges(6,3) * t650;
t569 = -mrSges(6,1) * t582 + mrSges(6,3) * t579 + Ifges(6,4) * t598 + Ifges(6,2) * t597 + Ifges(6,6) * t631 - t633 * t590 + t650 * t592;
t570 = mrSges(6,2) * t582 - mrSges(6,3) * t578 + Ifges(6,1) * t598 + Ifges(6,4) * t597 + Ifges(6,5) * t631 + t632 * t590 - t650 * t591;
t617 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t669 - Ifges(4,4) * t666) * qJD(1);
t619 = (Ifges(5,4) * qJD(3)) + (-Ifges(5,2) * t669 + Ifges(5,6) * t666) * qJD(1);
t550 = -mrSges(4,1) * t607 + mrSges(4,3) * t602 - mrSges(5,1) * t586 + mrSges(5,2) * t585 - t665 * t570 - t668 * t569 - pkin(4) * t683 - pkin(7) * t692 - pkin(3) * t564 + t705 * t638 + (-Ifges(4,2) - Ifges(5,3)) * t637 + t716 * qJDD(3) + (t617 - t619) * qJD(3) + t669 * t691;
t616 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t669 - Ifges(4,2) * t666) * qJD(1);
t618 = (Ifges(5,5) * qJD(3)) + (-Ifges(5,6) * t669 + Ifges(5,3) * t666) * qJD(1);
t552 = t680 + (Ifges(4,1) + Ifges(5,2)) * t638 - t705 * t637 + (-t616 + t618) * qJD(3) + t717 * qJDD(3) + mrSges(4,2) * t607 + mrSges(5,1) * t588 - mrSges(4,3) * t601 - mrSges(5,3) * t585 + pkin(4) * t567 - qJ(4) * t564 + t666 * t691;
t556 = qJDD(1) * mrSges(3,2) - t682;
t675 = mrSges(2,1) * t646 - mrSges(2,2) * t647 + mrSges(3,2) * t613 - mrSges(3,3) * t611 - pkin(1) * t556 - pkin(6) * t558 + qJ(2) * t674 - t666 * t550 + t669 * t552 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t566 = qJDD(3) * mrSges(5,2) + qJD(3) * t643 + t635 * t699 - t681;
t673 = -mrSges(4,2) * t602 - mrSges(5,3) * t586 - pkin(7) * t567 - t665 * t569 - pkin(3) * t566 + t668 * t570 + qJ(4) * (-t635 * t700 + t677) + mrSges(5,2) * t588 + mrSges(4,1) * t601 + t617 * t700 + t616 * t699 + (-t669 * t618 - t666 * t619) * qJD(1) + t717 * t638 + (-qJ(4) * mrSges(5,1) - t716) * t637 + t715 * qJDD(3);
t557 = -m(3) * g(3) + t693;
t549 = t673 + (t703 * t672) + t704 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t646 + mrSges(3,1) * t613 - qJ(2) * t557 + pkin(2) * t558;
t548 = -mrSges(3,1) * t611 + mrSges(2,3) * t647 - pkin(1) * t557 - pkin(2) * t714 - pkin(6) * t693 + t707 * g(3) - t703 * qJDD(1) - t669 * t550 - t666 * t552 + t704 * t672;
t1 = [-m(1) * g(1) + t694; -m(1) * g(2) + t701; (-m(1) - m(2) - m(3)) * g(3) + t693; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t701 - t667 * t548 + t670 * t549; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t694 + t670 * t548 + t667 * t549; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t675; t675; t556; t673; t566; t680;];
tauJB = t1;
