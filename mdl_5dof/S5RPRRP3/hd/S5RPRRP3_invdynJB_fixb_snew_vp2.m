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
% Datum: 2019-12-05 18:04
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:03:30
% EndTime: 2019-12-05 18:03:34
% DurationCPUTime: 2.88s
% Computational Cost: add. (28140->248), mult. (55929->301), div. (0->0), fcn. (33035->8), ass. (0->101)
t708 = Ifges(5,4) + Ifges(6,4);
t717 = Ifges(5,2) + Ifges(6,2);
t713 = Ifges(5,6) + Ifges(6,6);
t680 = sin(qJ(4));
t681 = sin(qJ(3));
t683 = cos(qJ(4));
t684 = cos(qJ(3));
t647 = (-t680 * t681 + t683 * t684) * qJD(1);
t703 = qJD(1) * qJD(3);
t700 = t684 * t703;
t655 = t681 * qJDD(1) + t700;
t656 = t684 * qJDD(1) - t681 * t703;
t616 = t647 * qJD(4) + t683 * t655 + t680 * t656;
t648 = (t680 * t684 + t681 * t683) * qJD(1);
t629 = -t647 * mrSges(6,1) + t648 * mrSges(6,2);
t682 = sin(qJ(1));
t685 = cos(qJ(1));
t662 = t685 * g(2) + t682 * g(3);
t652 = qJDD(1) * pkin(1) + t662;
t661 = t682 * g(2) - t685 * g(3);
t686 = qJD(1) ^ 2;
t654 = -t686 * pkin(1) + t661;
t678 = sin(pkin(8));
t679 = cos(pkin(8));
t633 = t678 * t652 + t679 * t654;
t628 = -t686 * pkin(2) + qJDD(1) * pkin(6) + t633;
t677 = -g(1) + qJDD(2);
t613 = -t681 * t628 + t684 * t677;
t599 = (-t655 + t700) * pkin(7) + (t681 * t684 * t686 + qJDD(3)) * pkin(3) + t613;
t614 = t684 * t628 + t681 * t677;
t705 = qJD(1) * t681;
t660 = qJD(3) * pkin(3) - pkin(7) * t705;
t676 = t684 ^ 2;
t600 = -t676 * t686 * pkin(3) + t656 * pkin(7) - qJD(3) * t660 + t614;
t594 = t683 * t599 - t680 * t600;
t672 = qJDD(3) + qJDD(4);
t673 = qJD(3) + qJD(4);
t588 = -0.2e1 * qJD(5) * t648 + (t647 * t673 - t616) * qJ(5) + (t647 * t648 + t672) * pkin(4) + t594;
t635 = -t673 * mrSges(6,2) + t647 * mrSges(6,3);
t702 = m(6) * t588 + t672 * mrSges(6,1) + t673 * t635;
t584 = -t616 * mrSges(6,3) - t648 * t629 + t702;
t595 = t680 * t599 + t683 * t600;
t615 = -t648 * qJD(4) - t680 * t655 + t683 * t656;
t637 = t673 * pkin(4) - t648 * qJ(5);
t640 = t647 ^ 2;
t590 = -t640 * pkin(4) + t615 * qJ(5) + 0.2e1 * qJD(5) * t647 - t673 * t637 + t595;
t714 = Ifges(5,5) + Ifges(6,5);
t715 = Ifges(5,1) + Ifges(6,1);
t706 = -t708 * t647 - t715 * t648 - t714 * t673;
t711 = t717 * t647 + t708 * t648 + t713 * t673;
t712 = Ifges(5,3) + Ifges(6,3);
t716 = mrSges(5,1) * t594 + mrSges(6,1) * t588 - mrSges(5,2) * t595 - mrSges(6,2) * t590 + pkin(4) * t584 + t713 * t615 + t714 * t616 + t706 * t647 + t711 * t648 + t712 * t672;
t630 = -t647 * mrSges(5,1) + t648 * mrSges(5,2);
t636 = -t673 * mrSges(5,2) + t647 * mrSges(5,3);
t578 = m(5) * t594 + t672 * mrSges(5,1) + t673 * t636 + (-t629 - t630) * t648 + (-mrSges(5,3) - mrSges(6,3)) * t616 + t702;
t638 = t673 * mrSges(6,1) - t648 * mrSges(6,3);
t639 = t673 * mrSges(5,1) - t648 * mrSges(5,3);
t701 = m(6) * t590 + t615 * mrSges(6,3) + t647 * t629;
t581 = m(5) * t595 + t615 * mrSges(5,3) + t647 * t630 + (-t638 - t639) * t673 + (-mrSges(5,2) - mrSges(6,2)) * t672 + t701;
t574 = t683 * t578 + t680 * t581;
t645 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t681 + Ifges(4,2) * t684) * qJD(1);
t646 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t681 + Ifges(4,4) * t684) * qJD(1);
t710 = mrSges(4,1) * t613 - mrSges(4,2) * t614 + Ifges(4,5) * t655 + Ifges(4,6) * t656 + Ifges(4,3) * qJDD(3) + pkin(3) * t574 + (t681 * t645 - t684 * t646) * qJD(1) + t716;
t653 = (-mrSges(4,1) * t684 + mrSges(4,2) * t681) * qJD(1);
t704 = qJD(1) * t684;
t659 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t704;
t572 = m(4) * t613 + qJDD(3) * mrSges(4,1) - t655 * mrSges(4,3) + qJD(3) * t659 - t653 * t705 + t574;
t658 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t705;
t695 = -t680 * t578 + t683 * t581;
t573 = m(4) * t614 - qJDD(3) * mrSges(4,2) + t656 * mrSges(4,3) - qJD(3) * t658 + t653 * t704 + t695;
t696 = -t681 * t572 + t684 * t573;
t563 = m(3) * t633 - t686 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t696;
t632 = t679 * t652 - t678 * t654;
t692 = -qJDD(1) * pkin(2) - t632;
t627 = -t686 * pkin(6) + t692;
t601 = -t656 * pkin(3) + t660 * t705 + (-pkin(7) * t676 - pkin(6)) * t686 + t692;
t592 = -t615 * pkin(4) - t640 * qJ(5) + t648 * t637 + qJDD(5) + t601;
t585 = m(6) * t592 - t615 * mrSges(6,1) + t616 * mrSges(6,2) - t647 * t635 + t648 * t638;
t690 = m(5) * t601 - t615 * mrSges(5,1) + t616 * mrSges(5,2) - t647 * t636 + t648 * t639 + t585;
t688 = -m(4) * t627 + t656 * mrSges(4,1) - t655 * mrSges(4,2) - t658 * t705 + t659 * t704 - t690;
t576 = m(3) * t632 + qJDD(1) * mrSges(3,1) - t686 * mrSges(3,2) + t688;
t560 = t678 * t563 + t679 * t576;
t566 = t684 * t572 + t681 * t573;
t707 = -t713 * t647 - t714 * t648 - t712 * t673;
t564 = m(3) * t677 + t566;
t697 = t679 * t563 - t678 * t576;
t557 = m(2) * t661 - t686 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t697;
t558 = m(2) * t662 + qJDD(1) * mrSges(2,1) - t686 * mrSges(2,2) + t560;
t698 = t685 * t557 - t682 * t558;
t694 = -t682 * t557 - t685 * t558;
t567 = -mrSges(5,1) * t601 + mrSges(5,3) * t595 - mrSges(6,1) * t592 + mrSges(6,3) * t590 - pkin(4) * t585 + qJ(5) * t701 + (-qJ(5) * t638 - t706) * t673 + (-qJ(5) * mrSges(6,2) + t713) * t672 + t707 * t648 + t708 * t616 + t717 * t615;
t568 = mrSges(5,2) * t601 + mrSges(6,2) * t592 - mrSges(5,3) * t594 - mrSges(6,3) * t588 - qJ(5) * t584 + t708 * t615 + t715 * t616 - t707 * t647 + t714 * t672 - t711 * t673;
t644 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t681 + Ifges(4,6) * t684) * qJD(1);
t553 = -mrSges(4,1) * t627 + mrSges(4,3) * t614 + Ifges(4,4) * t655 + Ifges(4,2) * t656 + Ifges(4,6) * qJDD(3) - pkin(3) * t690 + pkin(7) * t695 + qJD(3) * t646 + t683 * t567 + t680 * t568 - t644 * t705;
t555 = mrSges(4,2) * t627 - mrSges(4,3) * t613 + Ifges(4,1) * t655 + Ifges(4,4) * t656 + Ifges(4,5) * qJDD(3) - pkin(7) * t574 - qJD(3) * t645 - t680 * t567 + t683 * t568 + t644 * t704;
t691 = mrSges(2,1) * t662 + mrSges(3,1) * t632 - mrSges(2,2) * t661 - mrSges(3,2) * t633 + pkin(1) * t560 + pkin(2) * t688 + pkin(6) * t696 + t684 * t553 + t681 * t555 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t551 = -mrSges(3,1) * t677 + mrSges(3,3) * t633 + t686 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t566 - t710;
t550 = mrSges(3,2) * t677 - mrSges(3,3) * t632 + Ifges(3,5) * qJDD(1) - t686 * Ifges(3,6) - pkin(6) * t566 - t681 * t553 + t684 * t555;
t549 = -mrSges(2,2) * g(1) - mrSges(2,3) * t662 + Ifges(2,5) * qJDD(1) - t686 * Ifges(2,6) - qJ(2) * t560 + t679 * t550 - t678 * t551;
t548 = mrSges(2,1) * g(1) + mrSges(2,3) * t661 + t686 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t564 + qJ(2) * t697 + t678 * t550 + t679 * t551;
t1 = [(-m(1) - m(2)) * g(1) + t564; -m(1) * g(2) + t694; -m(1) * g(3) + t698; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t691; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t698 - t685 * t548 - t682 * t549; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t694 - t682 * t548 + t685 * t549; t691; t564; t710; t716; t585;];
tauJB = t1;
