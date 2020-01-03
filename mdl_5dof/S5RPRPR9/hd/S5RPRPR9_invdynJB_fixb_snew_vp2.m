% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRPR9
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRPR9_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR9_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR9_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR9_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR9_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR9_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR9_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:23:55
% EndTime: 2019-12-31 18:23:58
% DurationCPUTime: 2.08s
% Computational Cost: add. (18312->255), mult. (35572->303), div. (0->0), fcn. (17890->8), ass. (0->112)
t716 = Ifges(4,1) + Ifges(5,2);
t706 = Ifges(4,4) + Ifges(5,6);
t705 = Ifges(4,5) - Ifges(5,4);
t715 = Ifges(4,2) + Ifges(5,3);
t704 = Ifges(4,6) - Ifges(5,5);
t714 = Ifges(4,3) + Ifges(5,1);
t670 = sin(qJ(3));
t673 = cos(qJ(3));
t641 = (mrSges(5,2) * t673 - mrSges(5,3) * t670) * qJD(1);
t696 = qJD(1) * t673;
t650 = -mrSges(5,1) * t696 - qJD(3) * mrSges(5,3);
t694 = qJD(1) * qJD(3);
t693 = t670 * t694;
t645 = t673 * qJDD(1) - t693;
t695 = t670 * qJD(1);
t652 = pkin(4) * t695 - qJD(3) * pkin(7);
t665 = t673 ^ 2;
t676 = qJD(1) ^ 2;
t692 = t673 * t694;
t644 = t670 * qJDD(1) + t692;
t671 = sin(qJ(1));
t674 = cos(qJ(1));
t653 = t671 * g(1) - t674 * g(2);
t639 = qJDD(1) * pkin(1) + t653;
t654 = -t674 * g(1) - t671 * g(2);
t643 = -t676 * pkin(1) + t654;
t667 = sin(pkin(8));
t668 = cos(pkin(8));
t614 = t668 * t639 - t667 * t643;
t687 = -qJDD(1) * pkin(2) - t614;
t710 = -2 * qJD(4);
t679 = pkin(3) * t693 + t695 * t710 + (-t644 - t692) * qJ(4) + t687;
t709 = -pkin(3) - pkin(7);
t590 = -t652 * t695 + (-pkin(4) * t665 - pkin(6)) * t676 + t709 * t645 + t679;
t666 = -g(3) + qJDD(2);
t615 = t667 * t639 + t668 * t643;
t603 = -t676 * pkin(2) + qJDD(1) * pkin(6) + t615;
t600 = t670 * t603;
t640 = (-pkin(3) * t673 - qJ(4) * t670) * qJD(1);
t675 = qJD(3) ^ 2;
t688 = -t675 * qJ(4) + t640 * t695 + qJDD(4) + t600;
t708 = pkin(7) * t676;
t593 = t644 * pkin(4) + t709 * qJDD(3) + (-pkin(4) * t694 - t670 * t708 - t666) * t673 + t688;
t669 = sin(qJ(5));
t672 = cos(qJ(5));
t588 = -t669 * t590 + t672 * t593;
t637 = -t669 * qJD(3) - t672 * t696;
t612 = t637 * qJD(5) + t672 * qJDD(3) - t669 * t645;
t638 = t672 * qJD(3) - t669 * t696;
t616 = -t637 * mrSges(6,1) + t638 * mrSges(6,2);
t656 = qJD(5) + t695;
t617 = -t656 * mrSges(6,2) + t637 * mrSges(6,3);
t636 = qJDD(5) + t644;
t585 = m(6) * t588 + t636 * mrSges(6,1) - t612 * mrSges(6,3) - t638 * t616 + t656 * t617;
t589 = t672 * t590 + t669 * t593;
t611 = -t638 * qJD(5) - t669 * qJDD(3) - t672 * t645;
t618 = t656 * mrSges(6,1) - t638 * mrSges(6,3);
t586 = m(6) * t589 - t636 * mrSges(6,2) + t611 * mrSges(6,3) + t637 * t616 - t656 * t618;
t576 = t672 * t585 + t669 * t586;
t702 = t673 * t666;
t596 = -qJDD(3) * pkin(3) + t688 - t702;
t684 = -m(5) * t596 - t644 * mrSges(5,1) - t576;
t575 = qJDD(3) * mrSges(5,2) + qJD(3) * t650 + t641 * t695 - t684;
t599 = t673 * t603 + t670 * t666;
t682 = -t675 * pkin(3) + qJDD(3) * qJ(4) + t640 * t696 + t599;
t592 = -t665 * t708 + t645 * pkin(4) + ((2 * qJD(4)) + t652) * qJD(3) + t682;
t604 = Ifges(6,5) * t638 + Ifges(6,6) * t637 + Ifges(6,3) * t656;
t606 = Ifges(6,1) * t638 + Ifges(6,4) * t637 + Ifges(6,5) * t656;
t577 = -mrSges(6,1) * t592 + mrSges(6,3) * t589 + Ifges(6,4) * t612 + Ifges(6,2) * t611 + Ifges(6,6) * t636 - t638 * t604 + t656 * t606;
t605 = Ifges(6,4) * t638 + Ifges(6,2) * t637 + Ifges(6,6) * t656;
t578 = mrSges(6,2) * t592 - mrSges(6,3) * t588 + Ifges(6,1) * t612 + Ifges(6,4) * t611 + Ifges(6,5) * t636 + t637 * t604 - t656 * t605;
t595 = qJD(3) * t710 - t682;
t598 = -t600 + t702;
t651 = mrSges(5,1) * t695 + qJD(3) * mrSges(5,2);
t685 = -m(6) * t592 + t611 * mrSges(6,1) - t612 * mrSges(6,2) + t637 * t617 - t638 * t618;
t681 = -m(5) * t595 + qJDD(3) * mrSges(5,3) + qJD(3) * t651 + t641 * t696 - t685;
t697 = t705 * qJD(3) + (t716 * t670 + t706 * t673) * qJD(1);
t698 = t704 * qJD(3) + (t706 * t670 + t715 * t673) * qJD(1);
t713 = (t698 * t670 - t697 * t673) * qJD(1) + t714 * qJDD(3) + t705 * t644 + t704 * t645 + mrSges(4,1) * t598 - mrSges(4,2) * t599 + mrSges(5,2) * t596 - mrSges(5,3) * t595 - pkin(3) * t575 - pkin(7) * t576 + qJ(4) * (t645 * mrSges(5,1) + t681) - t669 * t577 + t672 * t578;
t707 = t676 * pkin(6);
t642 = (-mrSges(4,1) * t673 + mrSges(4,2) * t670) * qJD(1);
t649 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t696;
t573 = m(4) * t598 - t644 * mrSges(4,3) + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t649 - t650) * qJD(3) + (-t641 - t642) * t695 + t684;
t648 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t695;
t581 = t642 * t696 + m(4) * t599 - qJDD(3) * mrSges(4,2) - qJD(3) * t648 + (mrSges(4,3) + mrSges(5,1)) * t645 + t681;
t689 = -t670 * t573 + t673 * t581;
t565 = m(3) * t615 - t676 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t689;
t602 = t687 - t707;
t594 = -t645 * pkin(3) + t679 - t707;
t700 = -t669 * t585 + t672 * t586;
t686 = -m(5) * t594 - t645 * mrSges(5,2) + t651 * t695 - t700;
t678 = -m(4) * t602 + t649 * t696 + t645 * mrSges(4,1) + (-mrSges(4,2) + mrSges(5,3)) * t644 + (-t648 * t670 - t650 * t673) * qJD(1) + t686;
t570 = m(3) * t614 + qJDD(1) * mrSges(3,1) - t676 * mrSges(3,2) + t678;
t562 = t667 * t565 + t668 * t570;
t559 = m(2) * t653 + qJDD(1) * mrSges(2,1) - t676 * mrSges(2,2) + t562;
t690 = t668 * t565 - t667 * t570;
t560 = m(2) * t654 - t676 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t690;
t701 = t674 * t559 + t671 * t560;
t568 = t673 * t573 + t670 * t581;
t699 = t714 * qJD(3) + (t705 * t670 + t704 * t673) * qJD(1);
t566 = m(3) * t666 + t568;
t691 = -t671 * t559 + t674 * t560;
t683 = mrSges(6,1) * t588 - mrSges(6,2) * t589 + Ifges(6,5) * t612 + Ifges(6,6) * t611 + Ifges(6,3) * t636 + t638 * t605 - t637 * t606;
t574 = -t644 * mrSges(5,3) + t650 * t696 - t686;
t553 = -mrSges(4,1) * t602 - mrSges(5,1) * t595 + mrSges(5,2) * t594 + mrSges(4,3) * t599 - pkin(3) * t574 - pkin(4) * t685 - pkin(7) * t700 + t697 * qJD(3) + t704 * qJDD(3) - t672 * t577 - t669 * t578 + t706 * t644 + t715 * t645 - t699 * t695;
t555 = mrSges(5,1) * t596 + mrSges(4,2) * t602 - mrSges(4,3) * t598 - mrSges(5,3) * t594 + pkin(4) * t576 - qJ(4) * t574 - t698 * qJD(3) + t705 * qJDD(3) + t716 * t644 + t706 * t645 + t699 * t696 + t683;
t680 = mrSges(2,1) * t653 + mrSges(3,1) * t614 - mrSges(2,2) * t654 - mrSges(3,2) * t615 + pkin(1) * t562 + pkin(2) * t678 + pkin(6) * t689 + t673 * t553 + t670 * t555 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t551 = -mrSges(3,1) * t666 + mrSges(3,3) * t615 + t676 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t568 - t713;
t550 = mrSges(3,2) * t666 - mrSges(3,3) * t614 + Ifges(3,5) * qJDD(1) - t676 * Ifges(3,6) - pkin(6) * t568 - t670 * t553 + t673 * t555;
t549 = -mrSges(2,2) * g(3) - mrSges(2,3) * t653 + Ifges(2,5) * qJDD(1) - t676 * Ifges(2,6) - qJ(2) * t562 + t668 * t550 - t667 * t551;
t548 = mrSges(2,1) * g(3) + mrSges(2,3) * t654 + t676 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t566 + qJ(2) * t690 + t667 * t550 + t668 * t551;
t1 = [-m(1) * g(1) + t691; -m(1) * g(2) + t701; (-m(1) - m(2)) * g(3) + t566; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t701 - t671 * t548 + t674 * t549; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t691 + t674 * t548 + t671 * t549; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t680; t680; t566; t713; t575; t683;];
tauJB = t1;
