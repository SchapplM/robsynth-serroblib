% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRRR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2020-01-03 12:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRRR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR6_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR6_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR6_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:15:11
% EndTime: 2020-01-03 12:15:17
% DurationCPUTime: 6.17s
% Computational Cost: add. (117086->275), mult. (150157->346), div. (0->0), fcn. (97178->10), ass. (0->115)
t674 = qJDD(1) + qJDD(2);
t681 = sin(qJ(3));
t686 = cos(qJ(3));
t676 = qJD(1) + qJD(2);
t705 = qJD(3) * t676;
t653 = t681 * t674 + t686 * t705;
t683 = sin(qJ(1));
t688 = cos(qJ(1));
t665 = -t688 * g(2) - t683 * g(3);
t658 = qJDD(1) * pkin(1) + t665;
t664 = -t683 * g(2) + t688 * g(3);
t689 = qJD(1) ^ 2;
t659 = -t689 * pkin(1) + t664;
t682 = sin(qJ(2));
t687 = cos(qJ(2));
t639 = t682 * t658 + t687 * t659;
t672 = t676 ^ 2;
t636 = -t672 * pkin(2) + t674 * pkin(7) + t639;
t707 = t681 * t636;
t710 = pkin(3) * t672;
t614 = qJDD(3) * pkin(3) - t653 * pkin(8) - t707 + (pkin(8) * t705 + t681 * t710 - g(1)) * t686;
t626 = -t681 * g(1) + t686 * t636;
t654 = t686 * t674 - t681 * t705;
t709 = t676 * t681;
t662 = qJD(3) * pkin(3) - pkin(8) * t709;
t678 = t686 ^ 2;
t615 = t654 * pkin(8) - qJD(3) * t662 - t678 * t710 + t626;
t680 = sin(qJ(4));
t685 = cos(qJ(4));
t596 = t685 * t614 - t680 * t615;
t647 = (-t680 * t681 + t685 * t686) * t676;
t622 = t647 * qJD(4) + t685 * t653 + t680 * t654;
t648 = (t680 * t686 + t681 * t685) * t676;
t673 = qJDD(3) + qJDD(4);
t675 = qJD(3) + qJD(4);
t591 = (t647 * t675 - t622) * pkin(9) + (t647 * t648 + t673) * pkin(4) + t596;
t597 = t680 * t614 + t685 * t615;
t621 = -t648 * qJD(4) - t680 * t653 + t685 * t654;
t642 = t675 * pkin(4) - t648 * pkin(9);
t643 = t647 ^ 2;
t592 = -t643 * pkin(4) + t621 * pkin(9) - t675 * t642 + t597;
t679 = sin(qJ(5));
t684 = cos(qJ(5));
t589 = t684 * t591 - t679 * t592;
t631 = t684 * t647 - t679 * t648;
t603 = t631 * qJD(5) + t679 * t621 + t684 * t622;
t632 = t679 * t647 + t684 * t648;
t610 = -t631 * mrSges(6,1) + t632 * mrSges(6,2);
t670 = qJD(5) + t675;
t623 = -t670 * mrSges(6,2) + t631 * mrSges(6,3);
t669 = qJDD(5) + t673;
t586 = m(6) * t589 + t669 * mrSges(6,1) - t603 * mrSges(6,3) - t632 * t610 + t670 * t623;
t590 = t679 * t591 + t684 * t592;
t602 = -t632 * qJD(5) + t684 * t621 - t679 * t622;
t624 = t670 * mrSges(6,1) - t632 * mrSges(6,3);
t587 = m(6) * t590 - t669 * mrSges(6,2) + t602 * mrSges(6,3) + t631 * t610 - t670 * t624;
t577 = t684 * t586 + t679 * t587;
t634 = -t647 * mrSges(5,1) + t648 * mrSges(5,2);
t640 = -t675 * mrSges(5,2) + t647 * mrSges(5,3);
t574 = m(5) * t596 + t673 * mrSges(5,1) - t622 * mrSges(5,3) - t648 * t634 + t675 * t640 + t577;
t641 = t675 * mrSges(5,1) - t648 * mrSges(5,3);
t700 = -t679 * t586 + t684 * t587;
t575 = m(5) * t597 - t673 * mrSges(5,2) + t621 * mrSges(5,3) + t647 * t634 - t675 * t641 + t700;
t570 = t685 * t574 + t680 * t575;
t625 = -t686 * g(1) - t707;
t645 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t681 + Ifges(4,2) * t686) * t676;
t646 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t681 + Ifges(4,4) * t686) * t676;
t628 = Ifges(5,4) * t648 + Ifges(5,2) * t647 + Ifges(5,6) * t675;
t629 = Ifges(5,1) * t648 + Ifges(5,4) * t647 + Ifges(5,5) * t675;
t606 = Ifges(6,4) * t632 + Ifges(6,2) * t631 + Ifges(6,6) * t670;
t607 = Ifges(6,1) * t632 + Ifges(6,4) * t631 + Ifges(6,5) * t670;
t695 = -mrSges(6,1) * t589 + mrSges(6,2) * t590 - Ifges(6,5) * t603 - Ifges(6,6) * t602 - Ifges(6,3) * t669 - t632 * t606 + t631 * t607;
t692 = -mrSges(5,1) * t596 + mrSges(5,2) * t597 - Ifges(5,5) * t622 - Ifges(5,6) * t621 - Ifges(5,3) * t673 - pkin(4) * t577 - t648 * t628 + t647 * t629 + t695;
t711 = mrSges(4,1) * t625 - mrSges(4,2) * t626 + Ifges(4,5) * t653 + Ifges(4,6) * t654 + Ifges(4,3) * qJDD(3) + pkin(3) * t570 + (t681 * t645 - t686 * t646) * t676 - t692;
t708 = t676 * t686;
t652 = (-mrSges(4,1) * t686 + mrSges(4,2) * t681) * t676;
t661 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t708;
t568 = m(4) * t625 + qJDD(3) * mrSges(4,1) - t653 * mrSges(4,3) + qJD(3) * t661 - t652 * t709 + t570;
t660 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t709;
t701 = -t680 * t574 + t685 * t575;
t569 = m(4) * t626 - qJDD(3) * mrSges(4,2) + t654 * mrSges(4,3) - qJD(3) * t660 + t652 * t708 + t701;
t702 = -t681 * t568 + t686 * t569;
t560 = m(3) * t639 - t672 * mrSges(3,1) - t674 * mrSges(3,2) + t702;
t638 = t687 * t658 - t682 * t659;
t697 = -t674 * pkin(2) - t638;
t635 = -t672 * pkin(7) + t697;
t616 = -t654 * pkin(3) + t662 * t709 + (-pkin(8) * t678 - pkin(7)) * t672 + t697;
t594 = -t621 * pkin(4) - t643 * pkin(9) + t648 * t642 + t616;
t699 = m(6) * t594 - t602 * mrSges(6,1) + t603 * mrSges(6,2) - t631 * t623 + t632 * t624;
t693 = m(5) * t616 - t621 * mrSges(5,1) + t622 * mrSges(5,2) - t647 * t640 + t648 * t641 + t699;
t691 = -m(4) * t635 + t654 * mrSges(4,1) - t653 * mrSges(4,2) - t660 * t709 + t661 * t708 - t693;
t581 = m(3) * t638 + t674 * mrSges(3,1) - t672 * mrSges(3,2) + t691;
t703 = t687 * t560 - t682 * t581;
t554 = m(2) * t664 - t689 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t703;
t557 = t682 * t560 + t687 * t581;
t555 = m(2) * t665 + qJDD(1) * mrSges(2,1) - t689 * mrSges(2,2) + t557;
t706 = t683 * t554 + t688 * t555;
t562 = t686 * t568 + t681 * t569;
t704 = -t688 * t554 + t683 * t555;
t605 = Ifges(6,5) * t632 + Ifges(6,6) * t631 + Ifges(6,3) * t670;
t578 = -mrSges(6,1) * t594 + mrSges(6,3) * t590 + Ifges(6,4) * t603 + Ifges(6,2) * t602 + Ifges(6,6) * t669 - t632 * t605 + t670 * t607;
t579 = mrSges(6,2) * t594 - mrSges(6,3) * t589 + Ifges(6,1) * t603 + Ifges(6,4) * t602 + Ifges(6,5) * t669 + t631 * t605 - t670 * t606;
t627 = Ifges(5,5) * t648 + Ifges(5,6) * t647 + Ifges(5,3) * t675;
t563 = -mrSges(5,1) * t616 + mrSges(5,3) * t597 + Ifges(5,4) * t622 + Ifges(5,2) * t621 + Ifges(5,6) * t673 - pkin(4) * t699 + pkin(9) * t700 + t684 * t578 + t679 * t579 - t648 * t627 + t675 * t629;
t564 = mrSges(5,2) * t616 - mrSges(5,3) * t596 + Ifges(5,1) * t622 + Ifges(5,4) * t621 + Ifges(5,5) * t673 - pkin(9) * t577 - t679 * t578 + t684 * t579 + t647 * t627 - t675 * t628;
t644 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t681 + Ifges(4,6) * t686) * t676;
t548 = -mrSges(4,1) * t635 + mrSges(4,3) * t626 + Ifges(4,4) * t653 + Ifges(4,2) * t654 + Ifges(4,6) * qJDD(3) - pkin(3) * t693 + pkin(8) * t701 + qJD(3) * t646 + t685 * t563 + t680 * t564 - t644 * t709;
t550 = mrSges(4,2) * t635 - mrSges(4,3) * t625 + Ifges(4,1) * t653 + Ifges(4,4) * t654 + Ifges(4,5) * qJDD(3) - pkin(8) * t570 - qJD(3) * t645 - t680 * t563 + t685 * t564 + t644 * t708;
t696 = mrSges(3,1) * t638 - mrSges(3,2) * t639 + Ifges(3,3) * t674 + pkin(2) * t691 + pkin(7) * t702 + t686 * t548 + t681 * t550;
t694 = mrSges(2,1) * t665 - mrSges(2,2) * t664 + Ifges(2,3) * qJDD(1) + pkin(1) * t557 + t696;
t546 = mrSges(3,1) * g(1) + mrSges(3,3) * t639 + t672 * Ifges(3,5) + Ifges(3,6) * t674 - pkin(2) * t562 - t711;
t545 = -mrSges(3,2) * g(1) - mrSges(3,3) * t638 + Ifges(3,5) * t674 - t672 * Ifges(3,6) - pkin(7) * t562 - t681 * t548 + t686 * t550;
t544 = -mrSges(2,2) * g(1) - mrSges(2,3) * t665 + Ifges(2,5) * qJDD(1) - t689 * Ifges(2,6) - pkin(6) * t557 + t687 * t545 - t682 * t546;
t543 = Ifges(2,6) * qJDD(1) + t689 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t664 + t682 * t545 + t687 * t546 - pkin(1) * (-m(3) * g(1) + t562) + pkin(6) * t703;
t1 = [(-m(1) - m(2) - m(3)) * g(1) + t562; -m(1) * g(2) + t706; -m(1) * g(3) + t704; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t694; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t704 + t688 * t543 + t683 * t544; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t706 + t683 * t543 - t688 * t544; t694; t696; t711; -t692; -t695;];
tauJB = t1;
