% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2022-01-20 11:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRRP2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP2_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP2_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP2_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:49:01
% EndTime: 2022-01-20 11:49:05
% DurationCPUTime: 3.49s
% Computational Cost: add. (46882->251), mult. (59433->304), div. (0->0), fcn. (35123->8), ass. (0->104)
t707 = Ifges(5,4) + Ifges(6,4);
t717 = Ifges(5,2) + Ifges(6,2);
t713 = Ifges(5,6) + Ifges(6,6);
t673 = qJD(1) + qJD(2);
t676 = sin(qJ(4));
t677 = sin(qJ(3));
t680 = cos(qJ(4));
t681 = cos(qJ(3));
t645 = (-t676 * t677 + t680 * t681) * t673;
t671 = qJDD(1) + qJDD(2);
t700 = qJD(3) * t673;
t651 = t677 * t671 + t681 * t700;
t652 = t681 * t671 - t677 * t700;
t615 = t645 * qJD(4) + t680 * t651 + t676 * t652;
t646 = (t676 * t681 + t677 * t680) * t673;
t628 = -t645 * mrSges(6,1) + t646 * mrSges(6,2);
t679 = sin(qJ(1));
t683 = cos(qJ(1));
t662 = t679 * g(1) - t683 * g(2);
t656 = qJDD(1) * pkin(1) + t662;
t663 = -t683 * g(1) - t679 * g(2);
t684 = qJD(1) ^ 2;
t657 = -t684 * pkin(1) + t663;
t678 = sin(qJ(2));
t682 = cos(qJ(2));
t634 = t678 * t656 + t682 * t657;
t669 = t673 ^ 2;
t631 = -t669 * pkin(2) + t671 * pkin(7) + t634;
t704 = t677 * t631;
t708 = pkin(3) * t669;
t600 = qJDD(3) * pkin(3) - t651 * pkin(8) - t704 + (pkin(8) * t700 + t677 * t708 - g(3)) * t681;
t617 = -t677 * g(3) + t681 * t631;
t706 = t673 * t677;
t660 = qJD(3) * pkin(3) - pkin(8) * t706;
t675 = t681 ^ 2;
t601 = t652 * pkin(8) - qJD(3) * t660 - t675 * t708 + t617;
t595 = t680 * t600 - t676 * t601;
t670 = qJDD(3) + qJDD(4);
t672 = qJD(3) + qJD(4);
t589 = -0.2e1 * qJD(5) * t646 + (t645 * t672 - t615) * qJ(5) + (t645 * t646 + t670) * pkin(4) + t595;
t636 = -t672 * mrSges(6,2) + t645 * mrSges(6,3);
t699 = m(6) * t589 + t670 * mrSges(6,1) + t672 * t636;
t585 = -t615 * mrSges(6,3) - t646 * t628 + t699;
t596 = t676 * t600 + t680 * t601;
t614 = -t646 * qJD(4) - t676 * t651 + t680 * t652;
t638 = t672 * pkin(4) - t646 * qJ(5);
t641 = t645 ^ 2;
t591 = -t641 * pkin(4) + t614 * qJ(5) + 0.2e1 * qJD(5) * t645 - t672 * t638 + t596;
t714 = Ifges(5,5) + Ifges(6,5);
t715 = Ifges(5,1) + Ifges(6,1);
t701 = -t707 * t645 - t715 * t646 - t714 * t672;
t711 = t717 * t645 + t707 * t646 + t713 * t672;
t712 = Ifges(5,3) + Ifges(6,3);
t716 = mrSges(5,1) * t595 + mrSges(6,1) * t589 - mrSges(5,2) * t596 - mrSges(6,2) * t591 + pkin(4) * t585 + t713 * t614 + t714 * t615 + t701 * t645 + t711 * t646 + t712 * t670;
t629 = -t645 * mrSges(5,1) + t646 * mrSges(5,2);
t637 = -t672 * mrSges(5,2) + t645 * mrSges(5,3);
t580 = m(5) * t595 + t670 * mrSges(5,1) + t672 * t637 + (-t628 - t629) * t646 + (-mrSges(5,3) - mrSges(6,3)) * t615 + t699;
t639 = t672 * mrSges(6,1) - t646 * mrSges(6,3);
t640 = t672 * mrSges(5,1) - t646 * mrSges(5,3);
t698 = m(6) * t591 + t614 * mrSges(6,3) + t645 * t628;
t583 = m(5) * t596 + t614 * mrSges(5,3) + t645 * t629 + (-t639 - t640) * t672 + (-mrSges(5,2) - mrSges(6,2)) * t670 + t698;
t575 = t680 * t580 + t676 * t583;
t616 = -t681 * g(3) - t704;
t643 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t677 + Ifges(4,2) * t681) * t673;
t644 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t677 + Ifges(4,4) * t681) * t673;
t710 = mrSges(4,1) * t616 - mrSges(4,2) * t617 + Ifges(4,5) * t651 + Ifges(4,6) * t652 + Ifges(4,3) * qJDD(3) + pkin(3) * t575 + (t677 * t643 - t681 * t644) * t673 + t716;
t705 = t673 * t681;
t650 = (-mrSges(4,1) * t681 + mrSges(4,2) * t677) * t673;
t659 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t705;
t573 = m(4) * t616 + qJDD(3) * mrSges(4,1) - t651 * mrSges(4,3) + qJD(3) * t659 - t650 * t706 + t575;
t658 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t706;
t693 = -t676 * t580 + t680 * t583;
t574 = m(4) * t617 - qJDD(3) * mrSges(4,2) + t652 * mrSges(4,3) - qJD(3) * t658 + t650 * t705 + t693;
t694 = -t677 * t573 + t681 * t574;
t565 = m(3) * t634 - t669 * mrSges(3,1) - t671 * mrSges(3,2) + t694;
t633 = t682 * t656 - t678 * t657;
t691 = -t671 * pkin(2) - t633;
t630 = -t669 * pkin(7) + t691;
t602 = -t652 * pkin(3) + t660 * t706 + (-pkin(8) * t675 - pkin(7)) * t669 + t691;
t593 = -t614 * pkin(4) - t641 * qJ(5) + t646 * t638 + qJDD(5) + t602;
t586 = m(6) * t593 - t614 * mrSges(6,1) + t615 * mrSges(6,2) - t645 * t636 + t646 * t639;
t688 = m(5) * t602 - t614 * mrSges(5,1) + t615 * mrSges(5,2) - t645 * t637 + t646 * t640 + t586;
t686 = -m(4) * t630 + t652 * mrSges(4,1) - t651 * mrSges(4,2) - t658 * t706 + t659 * t705 - t688;
t577 = m(3) * t633 + t671 * mrSges(3,1) - t669 * mrSges(3,2) + t686;
t562 = t678 * t565 + t682 * t577;
t559 = m(2) * t662 + qJDD(1) * mrSges(2,1) - t684 * mrSges(2,2) + t562;
t695 = t682 * t565 - t678 * t577;
t560 = m(2) * t663 - t684 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t695;
t703 = t683 * t559 + t679 * t560;
t567 = t681 * t573 + t677 * t574;
t702 = -t713 * t645 - t714 * t646 - t712 * t672;
t696 = -t679 * t559 + t683 * t560;
t568 = -mrSges(5,1) * t602 + mrSges(5,3) * t596 - mrSges(6,1) * t593 + mrSges(6,3) * t591 - pkin(4) * t586 + qJ(5) * t698 + (-qJ(5) * t639 - t701) * t672 + (-qJ(5) * mrSges(6,2) + t713) * t670 + t702 * t646 + t707 * t615 + t717 * t614;
t569 = mrSges(5,2) * t602 + mrSges(6,2) * t593 - mrSges(5,3) * t595 - mrSges(6,3) * t589 - qJ(5) * t585 + t707 * t614 + t715 * t615 - t702 * t645 + t714 * t670 - t711 * t672;
t642 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t677 + Ifges(4,6) * t681) * t673;
t553 = -mrSges(4,1) * t630 + mrSges(4,3) * t617 + Ifges(4,4) * t651 + Ifges(4,2) * t652 + Ifges(4,6) * qJDD(3) - pkin(3) * t688 + pkin(8) * t693 + qJD(3) * t644 + t680 * t568 + t676 * t569 - t642 * t706;
t555 = mrSges(4,2) * t630 - mrSges(4,3) * t616 + Ifges(4,1) * t651 + Ifges(4,4) * t652 + Ifges(4,5) * qJDD(3) - pkin(8) * t575 - qJD(3) * t643 - t676 * t568 + t680 * t569 + t642 * t705;
t690 = mrSges(3,1) * t633 - mrSges(3,2) * t634 + Ifges(3,3) * t671 + pkin(2) * t686 + pkin(7) * t694 + t681 * t553 + t677 * t555;
t689 = mrSges(2,1) * t662 - mrSges(2,2) * t663 + Ifges(2,3) * qJDD(1) + pkin(1) * t562 + t690;
t551 = mrSges(3,1) * g(3) + mrSges(3,3) * t634 + t669 * Ifges(3,5) + Ifges(3,6) * t671 - pkin(2) * t567 - t710;
t550 = -mrSges(3,2) * g(3) - mrSges(3,3) * t633 + Ifges(3,5) * t671 - t669 * Ifges(3,6) - pkin(7) * t567 - t677 * t553 + t681 * t555;
t549 = -mrSges(2,2) * g(3) - mrSges(2,3) * t662 + Ifges(2,5) * qJDD(1) - t684 * Ifges(2,6) - pkin(6) * t562 + t682 * t550 - t678 * t551;
t548 = Ifges(2,6) * qJDD(1) + t684 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t663 + t678 * t550 + t682 * t551 - pkin(1) * (-m(3) * g(3) + t567) + pkin(6) * t695;
t1 = [-m(1) * g(1) + t696; -m(1) * g(2) + t703; (-m(1) - m(2) - m(3)) * g(3) + t567; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t703 - t679 * t548 + t683 * t549; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t696 + t683 * t548 + t679 * t549; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t689; t689; t690; t710; t716; t586;];
tauJB = t1;
