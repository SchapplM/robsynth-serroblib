% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPPPR1
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
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPPPR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR1_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR1_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR1_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:20:27
% EndTime: 2020-01-03 11:20:31
% DurationCPUTime: 3.73s
% Computational Cost: add. (35992->234), mult. (83592->324), div. (0->0), fcn. (50610->10), ass. (0->118)
t680 = sin(qJ(1));
t682 = cos(qJ(1));
t657 = -t682 * g(2) - t680 * g(3);
t651 = qJDD(1) * pkin(1) + t657;
t656 = -t680 * g(2) + t682 * g(3);
t683 = qJD(1) ^ 2;
t652 = -t683 * pkin(1) + t656;
t675 = sin(pkin(7));
t678 = cos(pkin(7));
t634 = t675 * t651 + t678 * t652;
t724 = -pkin(2) * t683 + qJDD(1) * qJ(3) + (2 * qJD(1) * qJD(3)) + t634;
t633 = t678 * t651 - t675 * t652;
t691 = -t683 * qJ(3) + qJDD(3) - t633;
t674 = sin(pkin(8));
t677 = cos(pkin(8));
t698 = -pkin(3) * t677 - qJ(4) * t674;
t716 = qJD(1) * t674;
t723 = (-pkin(2) + t698) * qJDD(1) + t691 - 0.2e1 * qJD(4) * t716;
t672 = -g(1) + qJDD(2);
t618 = t677 * t672 - t724 * t674;
t722 = mrSges(4,2) * t674;
t669 = t674 ^ 2;
t721 = t669 * t683;
t673 = sin(pkin(9));
t720 = t673 * t674;
t676 = cos(pkin(9));
t719 = t674 * t676;
t619 = t674 * t672 + t724 * t677;
t649 = (-mrSges(4,1) * t677 + t722) * qJD(1);
t648 = t698 * qJD(1);
t715 = t677 * qJD(1);
t612 = t648 * t715 + t619;
t696 = -pkin(4) * t677 - pkin(6) * t719;
t717 = t723 * t676;
t605 = t696 * qJDD(1) + (-t612 + (-pkin(4) * t669 * t676 + pkin(6) * t674 * t677) * t683) * t673 + t717;
t608 = t676 * t612 + t723 * t673;
t647 = t696 * qJD(1);
t711 = t673 ^ 2 * t721;
t713 = qJDD(1) * t674;
t606 = -pkin(6) * t673 * t713 - pkin(4) * t711 + t647 * t715 + t608;
t679 = sin(qJ(5));
t681 = cos(qJ(5));
t603 = t605 * t681 - t606 * t679;
t693 = (-t673 * t681 - t676 * t679) * t674;
t638 = qJD(1) * t693;
t692 = (-t673 * t679 + t676 * t681) * t674;
t639 = qJD(1) * t692;
t622 = -mrSges(6,1) * t638 + mrSges(6,2) * t639;
t625 = t638 * qJD(5) + qJDD(1) * t692;
t659 = qJD(5) - t715;
t631 = -mrSges(6,2) * t659 + mrSges(6,3) * t638;
t712 = t677 * qJDD(1);
t658 = qJDD(5) - t712;
t600 = m(6) * t603 + mrSges(6,1) * t658 - t625 * mrSges(6,3) - t622 * t639 + t631 * t659;
t604 = t605 * t679 + t606 * t681;
t624 = -t639 * qJD(5) + qJDD(1) * t693;
t632 = mrSges(6,1) * t659 - mrSges(6,3) * t639;
t601 = m(6) * t604 - mrSges(6,2) * t658 + t624 * mrSges(6,3) + t622 * t638 - t632 * t659;
t592 = t681 * t600 + t679 * t601;
t607 = -t673 * t612 + t717;
t702 = mrSges(5,1) * t673 + mrSges(5,2) * t676;
t640 = t702 * t716;
t694 = mrSges(5,2) * t677 - mrSges(5,3) * t720;
t642 = t694 * qJD(1);
t695 = -mrSges(5,1) * t677 - mrSges(5,3) * t719;
t590 = m(5) * t607 + t695 * qJDD(1) + (-t640 * t719 - t642 * t677) * qJD(1) + t592;
t643 = t695 * qJD(1);
t705 = -t679 * t600 + t681 * t601;
t591 = m(5) * t608 + t694 * qJDD(1) + (-t640 * t720 + t643 * t677) * qJD(1) + t705;
t706 = -t673 * t590 + t676 * t591;
t585 = m(4) * t619 + (qJDD(1) * mrSges(4,3) + qJD(1) * t649) * t677 + t706;
t611 = t648 * t716 + qJDD(4) - t618;
t609 = -pkin(6) * t711 + (pkin(4) * qJDD(1) * t673 + qJD(1) * t647 * t676) * t674 + t611;
t690 = m(6) * t609 - t624 * mrSges(6,1) + t625 * mrSges(6,2) - t638 * t631 + t639 * t632;
t686 = m(5) * t611 + t690;
t697 = t642 * t673 + t643 * t676;
t596 = m(4) * t618 + ((-mrSges(4,3) - t702) * qJDD(1) + (-t649 - t697) * qJD(1)) * t674 - t686;
t707 = t677 * t585 - t596 * t674;
t575 = m(3) * t634 - mrSges(3,1) * t683 - qJDD(1) * mrSges(3,2) + t707;
t588 = t676 * t590 + t673 * t591;
t628 = -qJDD(1) * pkin(2) + t691;
t687 = -m(4) * t628 + mrSges(4,1) * t712 - t588 + (t677 ^ 2 * t683 + t721) * mrSges(4,3);
t582 = m(3) * t633 - t683 * mrSges(3,2) + (mrSges(3,1) - t722) * qJDD(1) + t687;
t708 = t678 * t575 - t582 * t675;
t567 = m(2) * t656 - mrSges(2,1) * t683 - qJDD(1) * mrSges(2,2) + t708;
t570 = t675 * t575 + t678 * t582;
t568 = m(2) * t657 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t683 + t570;
t718 = t680 * t567 + t682 * t568;
t578 = t674 * t585 + t677 * t596;
t576 = m(3) * t672 + t578;
t709 = -t567 * t682 + t680 * t568;
t701 = Ifges(4,1) * t674 + Ifges(4,4) * t677;
t700 = Ifges(4,5) * t674 + Ifges(4,6) * t677;
t699 = Ifges(5,5) * t676 - Ifges(5,6) * t673;
t689 = -Ifges(5,5) * t677 + (Ifges(5,1) * t676 - Ifges(5,4) * t673) * t674;
t688 = -Ifges(5,6) * t677 + (Ifges(5,4) * t676 - Ifges(5,2) * t673) * t674;
t613 = Ifges(6,5) * t639 + Ifges(6,6) * t638 + Ifges(6,3) * t659;
t615 = Ifges(6,1) * t639 + Ifges(6,4) * t638 + Ifges(6,5) * t659;
t593 = -mrSges(6,1) * t609 + mrSges(6,3) * t604 + Ifges(6,4) * t625 + Ifges(6,2) * t624 + Ifges(6,6) * t658 - t613 * t639 + t615 * t659;
t614 = Ifges(6,4) * t639 + Ifges(6,2) * t638 + Ifges(6,6) * t659;
t594 = mrSges(6,2) * t609 - mrSges(6,3) * t603 + Ifges(6,1) * t625 + Ifges(6,4) * t624 + Ifges(6,5) * t658 + t613 * t638 - t614 * t659;
t635 = (-Ifges(5,3) * t677 + t699 * t674) * qJD(1);
t637 = t689 * qJD(1);
t579 = -mrSges(5,1) * t611 + mrSges(5,3) * t608 + t679 * t594 + t681 * t593 - pkin(4) * t690 + pkin(6) * t705 + (-t635 * t719 - t637 * t677) * qJD(1) + t688 * qJDD(1);
t636 = t688 * qJD(1);
t580 = mrSges(5,2) * t611 - mrSges(5,3) * t607 - pkin(6) * t592 - t679 * t593 + t681 * t594 + (-t635 * t720 + t636 * t677) * qJD(1) + t689 * qJDD(1);
t650 = t700 * qJD(1);
t563 = mrSges(4,2) * t628 - mrSges(4,3) * t618 - qJ(4) * t588 + t701 * qJDD(1) - t673 * t579 + t676 * t580 + t650 * t715;
t684 = mrSges(6,1) * t603 - mrSges(6,2) * t604 + Ifges(6,5) * t625 + Ifges(6,6) * t624 + Ifges(6,3) * t658 + t639 * t614 - t638 * t615;
t572 = -mrSges(4,1) * t628 - mrSges(5,1) * t607 + mrSges(5,2) * t608 + mrSges(4,3) * t619 - pkin(3) * t588 - pkin(4) * t592 + (Ifges(4,2) + Ifges(5,3)) * t712 + ((Ifges(4,4) - t699) * qJDD(1) + (-t636 * t676 - t637 * t673 - t650) * qJD(1)) * t674 - t684;
t587 = mrSges(4,2) * t713 - t687;
t685 = mrSges(2,1) * t657 + mrSges(3,1) * t633 - mrSges(2,2) * t656 - mrSges(3,2) * t634 + pkin(1) * t570 - pkin(2) * t587 + qJ(3) * t707 + t674 * t563 + t677 * t572 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t602 = (t697 * qJD(1) + t702 * qJDD(1)) * t674 + t686;
t561 = -mrSges(3,1) * t672 + mrSges(3,3) * t634 - mrSges(4,1) * t618 + mrSges(4,2) * t619 - t673 * t580 - t676 * t579 + pkin(3) * t602 - qJ(4) * t706 - pkin(2) * t578 + (Ifges(3,6) - t700) * qJDD(1) + (Ifges(3,5) - t674 * (Ifges(4,4) * t674 + Ifges(4,2) * t677) + t677 * t701) * t683;
t560 = mrSges(3,2) * t672 - mrSges(3,3) * t633 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t683 - qJ(3) * t578 + t563 * t677 - t572 * t674;
t559 = -mrSges(2,2) * g(1) - mrSges(2,3) * t657 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t683 - qJ(2) * t570 + t560 * t678 - t561 * t675;
t558 = mrSges(2,1) * g(1) + mrSges(2,3) * t656 + t683 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t576 + qJ(2) * t708 + t675 * t560 + t678 * t561;
t1 = [(-m(1) - m(2)) * g(1) + t576; -m(1) * g(2) + t718; -m(1) * g(3) + t709; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t685; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t709 + t682 * t558 + t680 * t559; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t718 + t680 * t558 - t682 * t559; t685; t576; t587; t602; t684;];
tauJB = t1;
