% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPPRP2
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPPRP2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP2_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP2_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP2_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP2_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP2_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:15
% EndTime: 2019-12-31 17:49:17
% DurationCPUTime: 2.22s
% Computational Cost: add. (21445->225), mult. (46239->272), div. (0->0), fcn. (27835->8), ass. (0->105)
t708 = Ifges(5,1) + Ifges(6,1);
t701 = Ifges(5,4) - Ifges(6,5);
t700 = Ifges(5,5) + Ifges(6,4);
t707 = Ifges(5,2) + Ifges(6,3);
t699 = Ifges(5,6) - Ifges(6,6);
t706 = Ifges(5,3) + Ifges(6,2);
t666 = qJD(1) ^ 2;
t662 = sin(qJ(4));
t660 = cos(pkin(8));
t704 = cos(qJ(4));
t682 = t660 * t704;
t658 = sin(pkin(8));
t687 = t658 * qJD(1);
t628 = -qJD(1) * t682 + t662 * t687;
t671 = t704 * t658 + t660 * t662;
t629 = t671 * qJD(1);
t612 = t628 * mrSges(6,1) - t629 * mrSges(6,3);
t689 = t628 * qJD(4);
t618 = t671 * qJDD(1) - t689;
t663 = sin(qJ(1));
t664 = cos(qJ(1));
t638 = t663 * g(1) - t664 * g(2);
t635 = qJDD(1) * pkin(1) + t638;
t639 = -t664 * g(1) - t663 * g(2);
t636 = -t666 * pkin(1) + t639;
t659 = sin(pkin(7));
t661 = cos(pkin(7));
t621 = t659 * t635 + t661 * t636;
t607 = -t666 * pkin(2) + qJDD(1) * qJ(3) + t621;
t657 = -g(3) + qJDD(2);
t686 = qJD(1) * qJD(3);
t690 = t660 * t657 - 0.2e1 * t658 * t686;
t703 = pkin(3) * t660;
t594 = (-pkin(6) * qJDD(1) + t666 * t703 - t607) * t658 + t690;
t598 = t658 * t657 + (t607 + 0.2e1 * t686) * t660;
t684 = qJDD(1) * t660;
t652 = t660 ^ 2;
t696 = t652 * t666;
t595 = -pkin(3) * t696 + pkin(6) * t684 + t598;
t590 = t704 * t594 - t662 * t595;
t611 = t628 * pkin(4) - t629 * qJ(5);
t665 = qJD(4) ^ 2;
t587 = -qJDD(4) * pkin(4) - t665 * qJ(5) + t629 * t611 + qJDD(5) - t590;
t627 = -t628 * mrSges(6,2) + qJD(4) * mrSges(6,3);
t677 = -m(6) * t587 + qJDD(4) * mrSges(6,1) + qJD(4) * t627;
t584 = t618 * mrSges(6,2) + t629 * t612 - t677;
t591 = t662 * t594 + t704 * t595;
t586 = -t665 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t628 * t611 + t591;
t685 = qJDD(1) * t658;
t688 = t629 * qJD(4);
t617 = -qJDD(1) * t682 + t662 * t685 + t688;
t626 = -qJD(4) * mrSges(6,1) + t629 * mrSges(6,2);
t683 = m(6) * t586 + qJDD(4) * mrSges(6,3) + qJD(4) * t626;
t692 = -t700 * qJD(4) + t701 * t628 - t708 * t629;
t694 = -t699 * qJD(4) + t707 * t628 - t701 * t629;
t705 = t706 * qJDD(4) - t699 * t617 + t700 * t618 - t692 * t628 - t694 * t629 + mrSges(5,1) * t590 - mrSges(6,1) * t587 - mrSges(5,2) * t591 + mrSges(6,3) * t586 - pkin(4) * t584 + qJ(5) * (-t617 * mrSges(6,2) - t628 * t612 + t683);
t702 = -mrSges(5,3) - mrSges(6,2);
t697 = mrSges(4,2) * t658;
t625 = qJD(4) * mrSges(5,1) - t629 * mrSges(5,3);
t691 = -t628 * mrSges(5,1) - t629 * mrSges(5,2) - t612;
t580 = m(5) * t591 - qJDD(4) * mrSges(5,2) - qJD(4) * t625 + t702 * t617 + t691 * t628 + t683;
t624 = -qJD(4) * mrSges(5,2) - t628 * mrSges(5,3);
t581 = m(5) * t590 + qJDD(4) * mrSges(5,1) + qJD(4) * t624 + t702 * t618 + t691 * t629 + t677;
t572 = t662 * t580 + t704 * t581;
t597 = -t658 * t607 + t690;
t672 = mrSges(4,3) * qJDD(1) + t666 * (-mrSges(4,1) * t660 + t697);
t570 = m(4) * t597 - t672 * t658 + t572;
t678 = t704 * t580 - t662 * t581;
t571 = m(4) * t598 + t672 * t660 + t678;
t679 = -t658 * t570 + t660 * t571;
t561 = m(3) * t621 - t666 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t679;
t620 = t661 * t635 - t659 * t636;
t673 = qJDD(3) - t620;
t600 = -qJDD(1) * pkin(2) - t666 * qJ(3) + t673;
t651 = t658 ^ 2;
t596 = (-pkin(2) - t703) * qJDD(1) + (-qJ(3) + (-t651 - t652) * pkin(6)) * t666 + t673;
t589 = -0.2e1 * qJD(5) * t629 + (-t618 + t689) * qJ(5) + (t617 + t688) * pkin(4) + t596;
t582 = m(6) * t589 + t617 * mrSges(6,1) - t618 * mrSges(6,3) - t629 * t626 + t628 * t627;
t669 = m(5) * t596 + t617 * mrSges(5,1) + t618 * mrSges(5,2) + t628 * t624 + t629 * t625 + t582;
t667 = -m(4) * t600 + mrSges(4,1) * t684 - t669 + (t651 * t666 + t696) * mrSges(4,3);
t574 = (mrSges(3,1) - t697) * qJDD(1) - t666 * mrSges(3,2) + m(3) * t620 + t667;
t558 = t659 * t561 + t661 * t574;
t555 = m(2) * t638 + qJDD(1) * mrSges(2,1) - t666 * mrSges(2,2) + t558;
t680 = t661 * t561 - t659 * t574;
t556 = m(2) * t639 - t666 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t680;
t695 = t664 * t555 + t663 * t556;
t564 = t660 * t570 + t658 * t571;
t693 = -t706 * qJD(4) + t699 * t628 - t700 * t629;
t562 = m(3) * t657 + t564;
t681 = -t663 * t555 + t664 * t556;
t676 = Ifges(4,1) * t658 + Ifges(4,4) * t660;
t675 = Ifges(4,4) * t658 + Ifges(4,2) * t660;
t674 = Ifges(4,5) * t658 + Ifges(4,6) * t660;
t565 = -mrSges(5,1) * t596 - mrSges(6,1) * t589 + mrSges(6,2) * t586 + mrSges(5,3) * t591 - pkin(4) * t582 - t692 * qJD(4) + t699 * qJDD(4) - t707 * t617 + t701 * t618 + t693 * t629;
t566 = mrSges(5,2) * t596 + mrSges(6,2) * t587 - mrSges(5,3) * t590 - mrSges(6,3) * t589 - qJ(5) * t582 + t694 * qJD(4) + t700 * qJDD(4) - t701 * t617 + t708 * t618 + t693 * t628;
t634 = t674 * qJD(1);
t549 = -mrSges(4,1) * t600 + mrSges(4,3) * t598 - pkin(3) * t669 + pkin(6) * t678 + t675 * qJDD(1) + t704 * t565 + t662 * t566 - t634 * t687;
t551 = t660 * qJD(1) * t634 + mrSges(4,2) * t600 - mrSges(4,3) * t597 - pkin(6) * t572 + t676 * qJDD(1) - t662 * t565 + t704 * t566;
t576 = mrSges(4,2) * t685 - t667;
t670 = mrSges(2,1) * t638 + mrSges(3,1) * t620 - mrSges(2,2) * t639 - mrSges(3,2) * t621 + pkin(1) * t558 - pkin(2) * t576 + qJ(3) * t679 + t660 * t549 + t658 * t551 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t547 = (Ifges(3,6) - t674) * qJDD(1) - mrSges(3,1) * t657 + mrSges(3,3) * t621 + mrSges(4,2) * t598 - mrSges(4,1) * t597 - pkin(3) * t572 - pkin(2) * t564 + (-t658 * t675 + t660 * t676 + Ifges(3,5)) * t666 - t705;
t546 = mrSges(3,2) * t657 - mrSges(3,3) * t620 + Ifges(3,5) * qJDD(1) - t666 * Ifges(3,6) - qJ(3) * t564 - t658 * t549 + t660 * t551;
t545 = -mrSges(2,2) * g(3) - mrSges(2,3) * t638 + Ifges(2,5) * qJDD(1) - t666 * Ifges(2,6) - qJ(2) * t558 + t661 * t546 - t659 * t547;
t544 = mrSges(2,1) * g(3) + mrSges(2,3) * t639 + t666 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t562 + qJ(2) * t680 + t659 * t546 + t661 * t547;
t1 = [-m(1) * g(1) + t681; -m(1) * g(2) + t695; (-m(1) - m(2)) * g(3) + t562; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t695 - t663 * t544 + t664 * t545; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t681 + t664 * t544 + t663 * t545; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t670; t670; t562; t576; t705; t584;];
tauJB = t1;
