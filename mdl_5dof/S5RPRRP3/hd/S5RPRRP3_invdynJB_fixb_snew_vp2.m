% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRRP3
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRRP3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP3_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP3_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP3_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:29:49
% EndTime: 2022-01-23 09:29:52
% DurationCPUTime: 2.99s
% Computational Cost: add. (28140->248), mult. (55929->301), div. (0->0), fcn. (33035->8), ass. (0->101)
t702 = Ifges(5,4) + Ifges(6,4);
t711 = Ifges(5,2) + Ifges(6,2);
t707 = Ifges(5,6) + Ifges(6,6);
t674 = sin(qJ(4));
t675 = sin(qJ(3));
t677 = cos(qJ(4));
t678 = cos(qJ(3));
t643 = (-t674 * t675 + t677 * t678) * qJD(1);
t696 = qJD(1) * qJD(3);
t693 = t678 * t696;
t651 = t675 * qJDD(1) + t693;
t652 = t678 * qJDD(1) - t675 * t696;
t612 = t643 * qJD(4) + t677 * t651 + t674 * t652;
t644 = (t674 * t678 + t675 * t677) * qJD(1);
t625 = -t643 * mrSges(6,1) + t644 * mrSges(6,2);
t676 = sin(qJ(1));
t679 = cos(qJ(1));
t657 = t676 * g(1) - t679 * g(2);
t648 = qJDD(1) * pkin(1) + t657;
t658 = -t679 * g(1) - t676 * g(2);
t680 = qJD(1) ^ 2;
t650 = -t680 * pkin(1) + t658;
t672 = sin(pkin(8));
t673 = cos(pkin(8));
t629 = t672 * t648 + t673 * t650;
t624 = -t680 * pkin(2) + qJDD(1) * pkin(6) + t629;
t671 = -g(3) + qJDD(2);
t609 = -t675 * t624 + t678 * t671;
t595 = (-t651 + t693) * pkin(7) + (t675 * t678 * t680 + qJDD(3)) * pkin(3) + t609;
t610 = t678 * t624 + t675 * t671;
t698 = qJD(1) * t675;
t656 = qJD(3) * pkin(3) - pkin(7) * t698;
t670 = t678 ^ 2;
t596 = -t670 * t680 * pkin(3) + t652 * pkin(7) - qJD(3) * t656 + t610;
t590 = t677 * t595 - t674 * t596;
t666 = qJDD(3) + qJDD(4);
t667 = qJD(3) + qJD(4);
t584 = -0.2e1 * qJD(5) * t644 + (t643 * t667 - t612) * qJ(5) + (t643 * t644 + t666) * pkin(4) + t590;
t631 = -t667 * mrSges(6,2) + t643 * mrSges(6,3);
t695 = m(6) * t584 + t666 * mrSges(6,1) + t667 * t631;
t580 = -t612 * mrSges(6,3) - t644 * t625 + t695;
t591 = t674 * t595 + t677 * t596;
t611 = -t644 * qJD(4) - t674 * t651 + t677 * t652;
t633 = t667 * pkin(4) - t644 * qJ(5);
t636 = t643 ^ 2;
t586 = -t636 * pkin(4) + t611 * qJ(5) + 0.2e1 * qJD(5) * t643 - t667 * t633 + t591;
t708 = Ifges(5,5) + Ifges(6,5);
t709 = Ifges(5,1) + Ifges(6,1);
t699 = -t702 * t643 - t709 * t644 - t708 * t667;
t705 = t711 * t643 + t644 * t702 + t707 * t667;
t706 = Ifges(5,3) + Ifges(6,3);
t710 = mrSges(5,1) * t590 + mrSges(6,1) * t584 - mrSges(5,2) * t591 - mrSges(6,2) * t586 + pkin(4) * t580 + t707 * t611 + t708 * t612 + t699 * t643 + t705 * t644 + t706 * t666;
t626 = -t643 * mrSges(5,1) + t644 * mrSges(5,2);
t632 = -t667 * mrSges(5,2) + t643 * mrSges(5,3);
t574 = m(5) * t590 + t666 * mrSges(5,1) + t667 * t632 + (-t625 - t626) * t644 + (-mrSges(5,3) - mrSges(6,3)) * t612 + t695;
t634 = t667 * mrSges(6,1) - t644 * mrSges(6,3);
t635 = t667 * mrSges(5,1) - t644 * mrSges(5,3);
t694 = m(6) * t586 + t611 * mrSges(6,3) + t643 * t625;
t577 = m(5) * t591 + t611 * mrSges(5,3) + t643 * t626 + (-t634 - t635) * t667 + (-mrSges(5,2) - mrSges(6,2)) * t666 + t694;
t570 = t677 * t574 + t674 * t577;
t641 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t675 + Ifges(4,2) * t678) * qJD(1);
t642 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t675 + Ifges(4,4) * t678) * qJD(1);
t704 = mrSges(4,1) * t609 - mrSges(4,2) * t610 + Ifges(4,5) * t651 + Ifges(4,6) * t652 + Ifges(4,3) * qJDD(3) + pkin(3) * t570 + (t675 * t641 - t678 * t642) * qJD(1) + t710;
t649 = (-mrSges(4,1) * t678 + mrSges(4,2) * t675) * qJD(1);
t697 = qJD(1) * t678;
t655 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t697;
t568 = m(4) * t609 + qJDD(3) * mrSges(4,1) - t651 * mrSges(4,3) + qJD(3) * t655 - t649 * t698 + t570;
t654 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t698;
t688 = -t674 * t574 + t677 * t577;
t569 = m(4) * t610 - qJDD(3) * mrSges(4,2) + t652 * mrSges(4,3) - qJD(3) * t654 + t649 * t697 + t688;
t689 = -t675 * t568 + t678 * t569;
t559 = m(3) * t629 - t680 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t689;
t628 = t673 * t648 - t672 * t650;
t686 = -qJDD(1) * pkin(2) - t628;
t623 = -t680 * pkin(6) + t686;
t597 = -t652 * pkin(3) + t656 * t698 + (-pkin(7) * t670 - pkin(6)) * t680 + t686;
t588 = -t611 * pkin(4) - t636 * qJ(5) + t644 * t633 + qJDD(5) + t597;
t581 = m(6) * t588 - t611 * mrSges(6,1) + t612 * mrSges(6,2) - t643 * t631 + t644 * t634;
t684 = m(5) * t597 - t611 * mrSges(5,1) + t612 * mrSges(5,2) - t643 * t632 + t644 * t635 + t581;
t682 = -m(4) * t623 + t652 * mrSges(4,1) - t651 * mrSges(4,2) - t654 * t698 + t655 * t697 - t684;
t572 = m(3) * t628 + qJDD(1) * mrSges(3,1) - t680 * mrSges(3,2) + t682;
t556 = t672 * t559 + t673 * t572;
t553 = m(2) * t657 + qJDD(1) * mrSges(2,1) - t680 * mrSges(2,2) + t556;
t690 = t673 * t559 - t672 * t572;
t554 = m(2) * t658 - t680 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t690;
t701 = t679 * t553 + t676 * t554;
t562 = t678 * t568 + t675 * t569;
t700 = -t707 * t643 - t708 * t644 - t706 * t667;
t560 = m(3) * t671 + t562;
t691 = -t676 * t553 + t679 * t554;
t563 = -mrSges(5,1) * t597 + mrSges(5,3) * t591 - mrSges(6,1) * t588 + mrSges(6,3) * t586 - pkin(4) * t581 + qJ(5) * t694 + (-qJ(5) * t634 - t699) * t667 + (-qJ(5) * mrSges(6,2) + t707) * t666 + t700 * t644 + t702 * t612 + t711 * t611;
t564 = mrSges(5,2) * t597 + mrSges(6,2) * t588 - mrSges(5,3) * t590 - mrSges(6,3) * t584 - qJ(5) * t580 + t702 * t611 + t709 * t612 - t700 * t643 + t708 * t666 - t705 * t667;
t640 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t675 + Ifges(4,6) * t678) * qJD(1);
t547 = -mrSges(4,1) * t623 + mrSges(4,3) * t610 + Ifges(4,4) * t651 + Ifges(4,2) * t652 + Ifges(4,6) * qJDD(3) - pkin(3) * t684 + pkin(7) * t688 + qJD(3) * t642 + t677 * t563 + t674 * t564 - t640 * t698;
t549 = mrSges(4,2) * t623 - mrSges(4,3) * t609 + Ifges(4,1) * t651 + Ifges(4,4) * t652 + Ifges(4,5) * qJDD(3) - pkin(7) * t570 - qJD(3) * t641 - t674 * t563 + t677 * t564 + t640 * t697;
t685 = mrSges(2,1) * t657 + mrSges(3,1) * t628 - mrSges(2,2) * t658 - mrSges(3,2) * t629 + pkin(1) * t556 + pkin(2) * t682 + pkin(6) * t689 + t678 * t547 + t675 * t549 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t545 = -mrSges(3,1) * t671 + mrSges(3,3) * t629 + t680 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t562 - t704;
t544 = mrSges(3,2) * t671 - mrSges(3,3) * t628 + Ifges(3,5) * qJDD(1) - t680 * Ifges(3,6) - pkin(6) * t562 - t675 * t547 + t678 * t549;
t543 = -mrSges(2,2) * g(3) - mrSges(2,3) * t657 + Ifges(2,5) * qJDD(1) - t680 * Ifges(2,6) - qJ(2) * t556 + t673 * t544 - t672 * t545;
t542 = mrSges(2,1) * g(3) + mrSges(2,3) * t658 + t680 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t560 + qJ(2) * t690 + t672 * t544 + t673 * t545;
t1 = [-m(1) * g(1) + t691; -m(1) * g(2) + t701; (-m(1) - m(2)) * g(3) + t560; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t701 - t676 * t542 + t679 * t543; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t691 + t679 * t542 + t676 * t543; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t685; t685; t560; t704; t710; t581;];
tauJB = t1;
