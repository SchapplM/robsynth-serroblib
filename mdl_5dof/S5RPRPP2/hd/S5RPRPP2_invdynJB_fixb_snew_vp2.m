% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRPP2
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRPP2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP2_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP2_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP2_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP2_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP2_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:10:49
% EndTime: 2019-12-31 18:10:51
% DurationCPUTime: 1.36s
% Computational Cost: add. (9193->223), mult. (17682->268), div. (0->0), fcn. (7625->6), ass. (0->93)
t718 = Ifges(4,1) + Ifges(5,1) + Ifges(6,1);
t701 = Ifges(4,4) - Ifges(5,5) - Ifges(6,4);
t700 = Ifges(4,5) + Ifges(5,4) - Ifges(6,5);
t717 = Ifges(4,2) + Ifges(5,3) + Ifges(6,2);
t699 = Ifges(4,6) - Ifges(5,6) + Ifges(6,6);
t716 = Ifges(4,3) + Ifges(5,2) + Ifges(6,3);
t676 = sin(qJ(3));
t678 = cos(qJ(3));
t640 = (-mrSges(5,1) * t678 - mrSges(5,3) * t676) * qJD(1);
t702 = qJD(1) * qJD(3);
t644 = t676 * qJDD(1) + t678 * t702;
t673 = -g(3) + qJDD(2);
t677 = sin(qJ(1));
t679 = cos(qJ(1));
t658 = t677 * g(1) - t679 * g(2);
t638 = qJDD(1) * pkin(1) + t658;
t659 = -t679 * g(1) - t677 * g(2);
t681 = qJD(1) ^ 2;
t643 = -t681 * pkin(1) + t659;
t674 = sin(pkin(7));
t675 = cos(pkin(7));
t609 = t674 * t638 + t675 * t643;
t606 = -t681 * pkin(2) + qJDD(1) * pkin(6) + t609;
t603 = t676 * t606;
t639 = (-pkin(3) * t678 - qJ(4) * t676) * qJD(1);
t680 = qJD(3) ^ 2;
t704 = qJD(1) * t676;
t688 = -t680 * qJ(4) + t639 * t704 + qJDD(4) + t603;
t697 = -0.2e1 * qJD(1) * qJD(5);
t710 = pkin(4) * t681;
t711 = pkin(3) + pkin(4);
t596 = t676 * t697 - t644 * qJ(5) - t711 * qJDD(3) + (qJ(5) * t702 - t676 * t710 - t673) * t678 + t688;
t641 = (mrSges(6,1) * t678 + mrSges(6,2) * t676) * qJD(1);
t703 = qJD(1) * t678;
t655 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t703;
t591 = m(6) * t596 - qJDD(3) * mrSges(6,1) - t644 * mrSges(6,3) - qJD(3) * t655 - t641 * t704;
t706 = t678 * t673;
t600 = -qJDD(3) * pkin(3) + t688 - t706;
t657 = mrSges(5,2) * t703 + qJD(3) * mrSges(5,3);
t685 = -m(5) * t600 + qJDD(3) * mrSges(5,1) + qJD(3) * t657 - t591;
t589 = t644 * mrSges(5,2) + t640 * t704 - t685;
t602 = t678 * t606 + t676 * t673;
t712 = 2 * qJD(4);
t599 = -t680 * pkin(3) + qJDD(3) * qJ(4) + qJD(3) * t712 + t639 * t703 + t602;
t645 = t678 * qJDD(1) - t676 * t702;
t651 = -qJD(3) * pkin(4) - qJ(5) * t704;
t672 = t678 ^ 2;
t595 = -t645 * qJ(5) + qJD(3) * t651 - t672 * t710 + t678 * t697 + t599;
t601 = -t603 + t706;
t654 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t704;
t652 = -qJD(3) * mrSges(6,1) - mrSges(6,3) * t704;
t690 = m(6) * t595 + qJDD(3) * mrSges(6,2) - t645 * mrSges(6,3) + qJD(3) * t652;
t686 = m(5) * t599 + qJDD(3) * mrSges(5,3) + qJD(3) * t654 + t640 * t703 + t690;
t694 = t700 * qJD(3) + (t718 * t676 + t701 * t678) * qJD(1);
t695 = -t699 * qJD(3) + (-t701 * t676 - t717 * t678) * qJD(1);
t715 = -(t695 * t676 + t694 * t678) * qJD(1) + t716 * qJDD(3) + t700 * t644 + t699 * t645 + mrSges(4,1) * t601 - mrSges(5,1) * t600 - mrSges(6,1) * t596 - mrSges(4,2) * t602 + mrSges(6,2) * t595 + mrSges(5,3) * t599 - pkin(3) * t589 - pkin(4) * t591 + qJ(4) * (t645 * mrSges(5,2) - t641 * t703 + t686);
t709 = t681 * pkin(6);
t708 = mrSges(4,3) + mrSges(5,2);
t707 = qJ(4) * t678;
t642 = (-mrSges(4,1) * t678 + mrSges(4,2) * t676) * qJD(1);
t653 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t704;
t585 = m(4) * t602 - qJDD(3) * mrSges(4,2) - qJD(3) * t653 + t708 * t645 + (-t641 + t642) * t703 + t686;
t656 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t703;
t586 = m(4) * t601 + qJDD(3) * mrSges(4,1) + qJD(3) * t656 - t708 * t644 + (-t640 - t642) * t704 + t685;
t691 = t678 * t585 - t676 * t586;
t575 = m(3) * t609 - t681 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t691;
t608 = t675 * t638 - t674 * t643;
t689 = qJDD(1) * pkin(2) + t608;
t687 = -t644 * qJ(4) - t689;
t593 = qJDD(5) + (-qJ(5) * t672 + pkin(6)) * t681 + t711 * t645 + (qJD(3) * t707 + (-pkin(3) * qJD(3) + t651 + t712) * t676) * qJD(1) - t687;
t590 = m(6) * t593 + t645 * mrSges(6,1) + t644 * mrSges(6,2) + t652 * t704 + t655 * t703;
t597 = -t645 * pkin(3) - t709 + (-0.2e1 * qJD(4) * t676 + (pkin(3) * t676 - t707) * qJD(3)) * qJD(1) + t687;
t587 = m(5) * t597 - t645 * mrSges(5,1) - t644 * mrSges(5,3) - t654 * t704 - t657 * t703 - t590;
t605 = -t689 - t709;
t682 = -m(4) * t605 + t645 * mrSges(4,1) - t644 * mrSges(4,2) - t653 * t704 + t656 * t703 - t587;
t580 = m(3) * t608 + qJDD(1) * mrSges(3,1) - t681 * mrSges(3,2) + t682;
t570 = t674 * t575 + t675 * t580;
t567 = m(2) * t658 + qJDD(1) * mrSges(2,1) - t681 * mrSges(2,2) + t570;
t692 = t675 * t575 - t674 * t580;
t568 = m(2) * t659 - t681 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t692;
t705 = t679 * t567 + t677 * t568;
t578 = t676 * t585 + t678 * t586;
t696 = -t716 * qJD(3) + (-t700 * t676 - t699 * t678) * qJD(1);
t576 = m(3) * t673 + t578;
t693 = -t677 * t567 + t679 * t568;
t563 = -mrSges(4,1) * t605 + mrSges(4,3) * t602 - mrSges(5,1) * t597 + mrSges(5,2) * t599 + mrSges(6,1) * t593 - mrSges(6,3) * t595 + pkin(4) * t590 - qJ(5) * t690 - pkin(3) * t587 + t717 * t645 + t701 * t644 + t699 * qJDD(3) + t694 * qJD(3) + (qJ(5) * t641 * t678 + t696 * t676) * qJD(1);
t572 = mrSges(4,2) * t605 + mrSges(5,2) * t600 + mrSges(6,2) * t593 - mrSges(4,3) * t601 - mrSges(5,3) * t597 - mrSges(6,3) * t596 - qJ(4) * t587 - qJ(5) * t591 + t695 * qJD(3) + t700 * qJDD(3) + t718 * t644 + t701 * t645 - t696 * t703;
t684 = mrSges(2,1) * t658 + mrSges(3,1) * t608 - mrSges(2,2) * t659 - mrSges(3,2) * t609 + pkin(1) * t570 + pkin(2) * t682 + pkin(6) * t691 + t678 * t563 + t676 * t572 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t561 = -mrSges(3,1) * t673 + mrSges(3,3) * t609 + t681 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t578 - t715;
t560 = mrSges(3,2) * t673 - mrSges(3,3) * t608 + Ifges(3,5) * qJDD(1) - t681 * Ifges(3,6) - pkin(6) * t578 - t676 * t563 + t678 * t572;
t559 = -mrSges(2,2) * g(3) - mrSges(2,3) * t658 + Ifges(2,5) * qJDD(1) - t681 * Ifges(2,6) - qJ(2) * t570 + t675 * t560 - t674 * t561;
t558 = mrSges(2,1) * g(3) + mrSges(2,3) * t659 + t681 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t576 + qJ(2) * t692 + t674 * t560 + t675 * t561;
t1 = [-m(1) * g(1) + t693; -m(1) * g(2) + t705; (-m(1) - m(2)) * g(3) + t576; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t705 - t677 * t558 + t679 * t559; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t693 + t679 * t558 + t677 * t559; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t684; t684; t576; t715; t589; t590;];
tauJB = t1;
