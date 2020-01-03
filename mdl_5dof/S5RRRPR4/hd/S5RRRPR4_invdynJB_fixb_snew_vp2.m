% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRPR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRPR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR4_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR4_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR4_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR4_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR4_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:14
% EndTime: 2019-12-31 21:11:17
% DurationCPUTime: 2.39s
% Computational Cost: add. (30538->254), mult. (38780->310), div. (0->0), fcn. (19896->8), ass. (0->107)
t724 = Ifges(4,1) + Ifges(5,1);
t715 = Ifges(4,4) - Ifges(5,5);
t714 = Ifges(4,5) + Ifges(5,4);
t723 = Ifges(4,2) + Ifges(5,3);
t713 = Ifges(4,6) - Ifges(5,6);
t722 = Ifges(4,3) + Ifges(5,2);
t674 = qJD(1) + qJD(2);
t681 = sin(qJ(3));
t685 = cos(qJ(3));
t644 = (-mrSges(5,1) * t685 - mrSges(5,3) * t681) * t674;
t672 = qJDD(1) + qJDD(2);
t705 = qJD(3) * t685;
t704 = t674 * t705;
t646 = t681 * t672 + t704;
t683 = sin(qJ(1));
t687 = cos(qJ(1));
t662 = t683 * g(1) - t687 * g(2);
t654 = qJDD(1) * pkin(1) + t662;
t663 = -t687 * g(1) - t683 * g(2);
t689 = qJD(1) ^ 2;
t655 = -t689 * pkin(1) + t663;
t682 = sin(qJ(2));
t686 = cos(qJ(2));
t624 = t682 * t654 + t686 * t655;
t670 = t674 ^ 2;
t621 = -t670 * pkin(2) + t672 * pkin(7) + t624;
t611 = -t681 * g(3) + t685 * t621;
t643 = (-pkin(3) * t685 - qJ(4) * t681) * t674;
t688 = qJD(3) ^ 2;
t710 = t674 * t685;
t718 = 2 * qJD(4);
t601 = -t688 * pkin(3) + qJDD(3) * qJ(4) + qJD(3) * t718 + t643 * t710 + t611;
t711 = t674 * t681;
t647 = -qJD(3) * t711 + t685 * t672;
t660 = -qJD(3) * pkin(4) - pkin(8) * t711;
t679 = t685 ^ 2;
t597 = -t679 * t670 * pkin(4) - t647 * pkin(8) + qJD(3) * t660 + t601;
t610 = -t685 * g(3) - t681 * t621;
t602 = -qJDD(3) * pkin(3) - t688 * qJ(4) + t643 * t711 + qJDD(4) - t610;
t598 = (-t646 + t704) * pkin(8) + (-t670 * t681 * t685 - qJDD(3)) * pkin(4) + t602;
t680 = sin(qJ(5));
t684 = cos(qJ(5));
t593 = -t680 * t597 + t684 * t598;
t636 = (-t680 * t681 - t684 * t685) * t674;
t609 = t636 * qJD(5) + t684 * t646 - t680 * t647;
t637 = (-t680 * t685 + t681 * t684) * t674;
t619 = -t636 * mrSges(6,1) + t637 * mrSges(6,2);
t673 = -qJD(3) + qJD(5);
t625 = -t673 * mrSges(6,2) + t636 * mrSges(6,3);
t671 = -qJDD(3) + qJDD(5);
t590 = m(6) * t593 + t671 * mrSges(6,1) - t609 * mrSges(6,3) - t637 * t619 + t673 * t625;
t594 = t684 * t597 + t680 * t598;
t608 = -t637 * qJD(5) - t680 * t646 - t684 * t647;
t626 = t673 * mrSges(6,1) - t637 * mrSges(6,3);
t591 = m(6) * t594 - t671 * mrSges(6,2) + t608 * mrSges(6,3) + t636 * t619 - t673 * t626;
t582 = t684 * t590 + t680 * t591;
t659 = mrSges(5,2) * t710 + qJD(3) * mrSges(5,3);
t694 = -m(5) * t602 + qJDD(3) * mrSges(5,1) + qJD(3) * t659 - t582;
t581 = t646 * mrSges(5,2) + t644 * t711 - t694;
t613 = Ifges(6,4) * t637 + Ifges(6,2) * t636 + Ifges(6,6) * t673;
t614 = Ifges(6,1) * t637 + Ifges(6,4) * t636 + Ifges(6,5) * t673;
t693 = -mrSges(6,1) * t593 + mrSges(6,2) * t594 - Ifges(6,5) * t609 - Ifges(6,6) * t608 - Ifges(6,3) * t671 - t637 * t613 + t636 * t614;
t657 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t711;
t700 = -t680 * t590 + t684 * t591;
t696 = m(5) * t601 + qJDD(3) * mrSges(5,3) + qJD(3) * t657 + t644 * t710 + t700;
t706 = (t724 * t681 + t715 * t685) * t674 + t714 * qJD(3);
t707 = (-t715 * t681 - t723 * t685) * t674 - t713 * qJD(3);
t721 = -(t707 * t681 + t706 * t685) * t674 + t722 * qJDD(3) + t714 * t646 + t713 * t647 + mrSges(4,1) * t610 - mrSges(5,1) * t602 - mrSges(4,2) * t611 + mrSges(5,3) * t601 - pkin(3) * t581 - pkin(4) * t582 + qJ(4) * (t647 * mrSges(5,2) + t696) + t693;
t717 = t670 * pkin(7);
t716 = mrSges(4,3) + mrSges(5,2);
t645 = (-mrSges(4,1) * t685 + mrSges(4,2) * t681) * t674;
t656 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t711;
t578 = m(4) * t611 - qJDD(3) * mrSges(4,2) - qJD(3) * t656 + t645 * t710 + t716 * t647 + t696;
t658 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t710;
t579 = m(4) * t610 + qJDD(3) * mrSges(4,1) + qJD(3) * t658 + (-t644 - t645) * t711 - t716 * t646 + t694;
t701 = t685 * t578 - t681 * t579;
t572 = m(3) * t624 - t670 * mrSges(3,1) - t672 * mrSges(3,2) + t701;
t623 = t686 * t654 - t682 * t655;
t699 = t672 * pkin(2) + t623;
t697 = -t646 * qJ(4) - t699;
t599 = -t647 * pkin(3) - t717 + (-0.2e1 * qJD(4) * t681 + (pkin(3) * t681 - qJ(4) * t685) * qJD(3)) * t674 + t697;
t596 = (-pkin(8) * t679 + pkin(7)) * t670 + (pkin(3) + pkin(4)) * t647 + (qJ(4) * t705 + (-pkin(3) * qJD(3) + t660 + t718) * t681) * t674 - t697;
t698 = -m(6) * t596 + t608 * mrSges(6,1) - t609 * mrSges(6,2) + t636 * t625 - t637 * t626;
t588 = m(5) * t599 - t647 * mrSges(5,1) - t646 * mrSges(5,3) - t657 * t711 - t659 * t710 + t698;
t620 = -t699 - t717;
t691 = -m(4) * t620 + t647 * mrSges(4,1) - t646 * mrSges(4,2) - t656 * t711 + t658 * t710 - t588;
t584 = m(3) * t623 + t672 * mrSges(3,1) - t670 * mrSges(3,2) + t691;
t569 = t682 * t572 + t686 * t584;
t566 = m(2) * t662 + qJDD(1) * mrSges(2,1) - t689 * mrSges(2,2) + t569;
t702 = t686 * t572 - t682 * t584;
t567 = m(2) * t663 - t689 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t702;
t709 = t687 * t566 + t683 * t567;
t574 = t681 * t578 + t685 * t579;
t708 = (t714 * t681 + t713 * t685) * t674 + t722 * qJD(3);
t703 = -t683 * t566 + t687 * t567;
t612 = Ifges(6,5) * t637 + Ifges(6,6) * t636 + Ifges(6,3) * t673;
t586 = -mrSges(6,1) * t596 + mrSges(6,3) * t594 + Ifges(6,4) * t609 + Ifges(6,2) * t608 + Ifges(6,6) * t671 - t637 * t612 + t673 * t614;
t587 = mrSges(6,2) * t596 - mrSges(6,3) * t593 + Ifges(6,1) * t609 + Ifges(6,4) * t608 + Ifges(6,5) * t671 + t636 * t612 - t673 * t613;
t560 = -mrSges(4,1) * t620 - mrSges(5,1) * t599 + mrSges(5,2) * t601 + mrSges(4,3) * t611 - pkin(3) * t588 - pkin(4) * t698 - pkin(8) * t700 + t706 * qJD(3) + t713 * qJDD(3) - t684 * t586 - t680 * t587 + t715 * t646 + t723 * t647 - t708 * t711;
t562 = mrSges(4,2) * t620 + mrSges(5,2) * t602 - mrSges(4,3) * t610 - mrSges(5,3) * t599 - pkin(8) * t582 - qJ(4) * t588 + t707 * qJD(3) + t714 * qJDD(3) - t680 * t586 + t684 * t587 + t724 * t646 + t715 * t647 + t708 * t710;
t695 = mrSges(3,1) * t623 - mrSges(3,2) * t624 + Ifges(3,3) * t672 + pkin(2) * t691 + pkin(7) * t701 + t685 * t560 + t681 * t562;
t692 = mrSges(2,1) * t662 - mrSges(2,2) * t663 + Ifges(2,3) * qJDD(1) + pkin(1) * t569 + t695;
t558 = mrSges(3,1) * g(3) + mrSges(3,3) * t624 + t670 * Ifges(3,5) + Ifges(3,6) * t672 - pkin(2) * t574 - t721;
t557 = -mrSges(3,2) * g(3) - mrSges(3,3) * t623 + Ifges(3,5) * t672 - t670 * Ifges(3,6) - pkin(7) * t574 - t681 * t560 + t685 * t562;
t556 = -mrSges(2,2) * g(3) - mrSges(2,3) * t662 + Ifges(2,5) * qJDD(1) - t689 * Ifges(2,6) - pkin(6) * t569 + t686 * t557 - t682 * t558;
t555 = Ifges(2,6) * qJDD(1) + t689 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t663 + t682 * t557 + t686 * t558 - pkin(1) * (-m(3) * g(3) + t574) + pkin(6) * t702;
t1 = [-m(1) * g(1) + t703; -m(1) * g(2) + t709; (-m(1) - m(2) - m(3)) * g(3) + t574; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t709 - t683 * t555 + t687 * t556; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t703 + t687 * t555 + t683 * t556; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t692; t692; t695; t721; t581; -t693;];
tauJB = t1;
