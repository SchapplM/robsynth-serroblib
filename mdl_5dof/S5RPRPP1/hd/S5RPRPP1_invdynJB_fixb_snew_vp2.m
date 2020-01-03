% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRPP1
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
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
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRPP1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP1_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP1_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP1_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP1_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:08:51
% EndTime: 2019-12-31 18:08:53
% DurationCPUTime: 2.49s
% Computational Cost: add. (24385->246), mult. (51173->301), div. (0->0), fcn. (29289->8), ass. (0->101)
t686 = Ifges(5,1) + Ifges(6,1);
t680 = Ifges(5,4) - Ifges(6,5);
t679 = Ifges(5,5) + Ifges(6,4);
t685 = -Ifges(5,2) - Ifges(6,3);
t684 = -Ifges(6,2) - Ifges(5,3);
t678 = Ifges(5,6) - Ifges(6,6);
t650 = sin(qJ(1));
t652 = cos(qJ(1));
t633 = t650 * g(1) - t652 * g(2);
t624 = qJDD(1) * pkin(1) + t633;
t634 = -t652 * g(1) - t650 * g(2);
t654 = qJD(1) ^ 2;
t626 = -t654 * pkin(1) + t634;
t647 = sin(pkin(7));
t648 = cos(pkin(7));
t602 = t647 * t624 + t648 * t626;
t590 = -t654 * pkin(2) + qJDD(1) * pkin(6) + t602;
t645 = -g(3) + qJDD(2);
t649 = sin(qJ(3));
t651 = cos(qJ(3));
t580 = -t649 * t590 + t651 * t645;
t669 = qJD(1) * qJD(3);
t666 = t651 * t669;
t627 = t649 * qJDD(1) + t666;
t577 = (-t627 + t666) * qJ(4) + (t649 * t651 * t654 + qJDD(3)) * pkin(3) + t580;
t581 = t651 * t590 + t649 * t645;
t628 = t651 * qJDD(1) - t649 * t669;
t671 = qJD(1) * t649;
t630 = qJD(3) * pkin(3) - qJ(4) * t671;
t644 = t651 ^ 2;
t578 = -t644 * t654 * pkin(3) + t628 * qJ(4) - qJD(3) * t630 + t581;
t646 = sin(pkin(8));
t670 = qJD(1) * t651;
t677 = cos(pkin(8));
t612 = t646 * t671 - t677 * t670;
t682 = -2 * qJD(4);
t574 = t646 * t577 + t677 * t578 + t612 * t682;
t603 = t646 * t627 - t677 * t628;
t613 = (t646 * t651 + t677 * t649) * qJD(1);
t608 = qJD(3) * mrSges(5,1) - t613 * mrSges(5,3);
t594 = t612 * pkin(4) - t613 * qJ(5);
t653 = qJD(3) ^ 2;
t569 = -t653 * pkin(4) + qJDD(3) * qJ(5) + 0.2e1 * qJD(5) * qJD(3) - t612 * t594 + t574;
t609 = -qJD(3) * mrSges(6,1) + t613 * mrSges(6,2);
t667 = m(6) * t569 + qJDD(3) * mrSges(6,3) + qJD(3) * t609;
t595 = t612 * mrSges(6,1) - t613 * mrSges(6,3);
t672 = -t612 * mrSges(5,1) - t613 * mrSges(5,2) - t595;
t681 = -mrSges(5,3) - mrSges(6,2);
t561 = m(5) * t574 - qJDD(3) * mrSges(5,2) - qJD(3) * t608 + t681 * t603 + t672 * t612 + t667;
t658 = t677 * t577 - t646 * t578;
t573 = t613 * t682 + t658;
t604 = t677 * t627 + t646 * t628;
t607 = -qJD(3) * mrSges(5,2) - t612 * mrSges(5,3);
t570 = -qJDD(3) * pkin(4) - t653 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t594) * t613 - t658;
t610 = -t612 * mrSges(6,2) + qJD(3) * mrSges(6,3);
t661 = -m(6) * t570 + qJDD(3) * mrSges(6,1) + qJD(3) * t610;
t562 = m(5) * t573 + qJDD(3) * mrSges(5,1) + qJD(3) * t607 + t681 * t604 + t672 * t613 + t661;
t555 = t646 * t561 + t677 * t562;
t566 = t604 * mrSges(6,2) + t613 * t595 - t661;
t618 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t649 + Ifges(4,2) * t651) * qJD(1);
t619 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t649 + Ifges(4,4) * t651) * qJD(1);
t673 = t679 * qJD(3) - t680 * t612 + t686 * t613;
t674 = t678 * qJD(3) + t685 * t612 + t680 * t613;
t683 = (t649 * t618 - t651 * t619) * qJD(1) + (Ifges(4,3) - t684) * qJDD(3) - t678 * t603 + t679 * t604 + t673 * t612 + t674 * t613 + mrSges(4,1) * t580 + mrSges(5,1) * t573 - mrSges(6,1) * t570 - mrSges(4,2) * t581 - mrSges(5,2) * t574 + mrSges(6,3) * t569 + Ifges(4,5) * t627 + Ifges(4,6) * t628 + pkin(3) * t555 - pkin(4) * t566 + qJ(5) * (-t603 * mrSges(6,2) - t612 * t595 + t667);
t625 = (-mrSges(4,1) * t651 + mrSges(4,2) * t649) * qJD(1);
t632 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t670;
t551 = m(4) * t580 + qJDD(3) * mrSges(4,1) - t627 * mrSges(4,3) + qJD(3) * t632 - t625 * t671 + t555;
t631 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t671;
t662 = t677 * t561 - t646 * t562;
t552 = m(4) * t581 - qJDD(3) * mrSges(4,2) + t628 * mrSges(4,3) - qJD(3) * t631 + t625 * t670 + t662;
t663 = -t649 * t551 + t651 * t552;
t544 = m(3) * t602 - t654 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t663;
t601 = t648 * t624 - t647 * t626;
t659 = -qJDD(1) * pkin(2) - t601;
t579 = -t628 * pkin(3) + qJDD(4) + t630 * t671 + (-qJ(4) * t644 - pkin(6)) * t654 + t659;
t572 = -0.2e1 * qJD(5) * t613 + (qJD(3) * t612 - t604) * qJ(5) + (qJD(3) * t613 + t603) * pkin(4) + t579;
t567 = m(6) * t572 + t603 * mrSges(6,1) - t604 * mrSges(6,3) - t613 * t609 + t612 * t610;
t564 = m(5) * t579 + t603 * mrSges(5,1) + t604 * mrSges(5,2) + t612 * t607 + t613 * t608 + t567;
t589 = -t654 * pkin(6) + t659;
t656 = -m(4) * t589 + t628 * mrSges(4,1) - t627 * mrSges(4,2) - t631 * t671 + t632 * t670 - t564;
t557 = m(3) * t601 + qJDD(1) * mrSges(3,1) - t654 * mrSges(3,2) + t656;
t541 = t647 * t544 + t648 * t557;
t538 = m(2) * t633 + qJDD(1) * mrSges(2,1) - t654 * mrSges(2,2) + t541;
t664 = t648 * t544 - t647 * t557;
t539 = m(2) * t634 - t654 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t664;
t676 = t652 * t538 + t650 * t539;
t547 = t651 * t551 + t649 * t552;
t675 = t684 * qJD(3) + t678 * t612 - t679 * t613;
t545 = m(3) * t645 + t547;
t665 = -t650 * t538 + t652 * t539;
t553 = -mrSges(5,1) * t579 - mrSges(6,1) * t572 + mrSges(6,2) * t569 + mrSges(5,3) * t574 - pkin(4) * t567 + t673 * qJD(3) + t678 * qJDD(3) + t685 * t603 + t680 * t604 + t675 * t613;
t554 = mrSges(5,2) * t579 + mrSges(6,2) * t570 - mrSges(5,3) * t573 - mrSges(6,3) * t572 - qJ(5) * t567 - t674 * qJD(3) + t679 * qJDD(3) - t680 * t603 + t686 * t604 + t675 * t612;
t617 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t649 + Ifges(4,6) * t651) * qJD(1);
t532 = -mrSges(4,1) * t589 + mrSges(4,3) * t581 + Ifges(4,4) * t627 + Ifges(4,2) * t628 + Ifges(4,6) * qJDD(3) - pkin(3) * t564 + qJ(4) * t662 + qJD(3) * t619 + t677 * t553 + t646 * t554 - t617 * t671;
t534 = mrSges(4,2) * t589 - mrSges(4,3) * t580 + Ifges(4,1) * t627 + Ifges(4,4) * t628 + Ifges(4,5) * qJDD(3) - qJ(4) * t555 - qJD(3) * t618 - t646 * t553 + t677 * t554 + t617 * t670;
t657 = mrSges(2,1) * t633 + mrSges(3,1) * t601 - mrSges(2,2) * t634 - mrSges(3,2) * t602 + pkin(1) * t541 + pkin(2) * t656 + pkin(6) * t663 + t651 * t532 + t649 * t534 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t530 = -mrSges(3,1) * t645 + mrSges(3,3) * t602 + t654 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t547 - t683;
t529 = mrSges(3,2) * t645 - mrSges(3,3) * t601 + Ifges(3,5) * qJDD(1) - t654 * Ifges(3,6) - pkin(6) * t547 - t649 * t532 + t651 * t534;
t528 = -mrSges(2,2) * g(3) - mrSges(2,3) * t633 + Ifges(2,5) * qJDD(1) - t654 * Ifges(2,6) - qJ(2) * t541 + t648 * t529 - t647 * t530;
t527 = mrSges(2,1) * g(3) + mrSges(2,3) * t634 + t654 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t545 + qJ(2) * t664 + t647 * t529 + t648 * t530;
t1 = [-m(1) * g(1) + t665; -m(1) * g(2) + t676; (-m(1) - m(2)) * g(3) + t545; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t676 - t650 * t527 + t652 * t528; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t665 + t652 * t527 + t650 * t528; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t657; t657; t545; t683; t564; t566;];
tauJB = t1;
