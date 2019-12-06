% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRPR1
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRPR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR1_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR1_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR1_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR1_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:24
% EndTime: 2019-12-05 17:47:27
% DurationCPUTime: 3.43s
% Computational Cost: add. (36245->264), mult. (78766->327), div. (0->0), fcn. (49520->8), ass. (0->106)
t660 = sin(qJ(1));
t663 = cos(qJ(1));
t640 = -t663 * g(1) - t660 * g(2);
t673 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t640;
t688 = -pkin(1) - pkin(6);
t687 = mrSges(2,1) - mrSges(3,2);
t686 = Ifges(2,5) - Ifges(3,4);
t685 = (-Ifges(2,6) + Ifges(3,5));
t639 = t660 * g(1) - t663 * g(2);
t664 = qJD(1) ^ 2;
t672 = -t664 * qJ(2) + qJDD(2) - t639;
t612 = t688 * qJDD(1) + t672;
t659 = sin(qJ(3));
t662 = cos(qJ(3));
t603 = t659 * g(3) + t662 * t612;
t681 = qJD(1) * qJD(3);
t679 = t659 * t681;
t634 = qJDD(1) * t662 - t679;
t583 = (-t634 - t679) * qJ(4) + (-t659 * t662 * t664 + qJDD(3)) * pkin(3) + t603;
t604 = -g(3) * t662 + t659 * t612;
t633 = -qJDD(1) * t659 - t662 * t681;
t682 = qJD(1) * t662;
t637 = qJD(3) * pkin(3) - qJ(4) * t682;
t653 = t659 ^ 2;
t584 = -pkin(3) * t653 * t664 + qJ(4) * t633 - qJD(3) * t637 + t604;
t656 = sin(pkin(8));
t657 = cos(pkin(8));
t622 = (-t656 * t659 + t657 * t662) * qJD(1);
t566 = -0.2e1 * qJD(4) * t622 + t657 * t583 - t656 * t584;
t601 = t633 * t656 + t634 * t657;
t621 = (-t656 * t662 - t657 * t659) * qJD(1);
t562 = (qJD(3) * t621 - t601) * pkin(7) + (t621 * t622 + qJDD(3)) * pkin(4) + t566;
t567 = 0.2e1 * qJD(4) * t621 + t656 * t583 + t657 * t584;
t600 = t633 * t657 - t634 * t656;
t611 = qJD(3) * pkin(4) - pkin(7) * t622;
t620 = t621 ^ 2;
t563 = -pkin(4) * t620 + pkin(7) * t600 - qJD(3) * t611 + t567;
t658 = sin(qJ(5));
t661 = cos(qJ(5));
t560 = t562 * t661 - t563 * t658;
t593 = t621 * t661 - t622 * t658;
t574 = qJD(5) * t593 + t600 * t658 + t601 * t661;
t594 = t621 * t658 + t622 * t661;
t579 = -mrSges(6,1) * t593 + mrSges(6,2) * t594;
t646 = qJD(3) + qJD(5);
t588 = -mrSges(6,2) * t646 + mrSges(6,3) * t593;
t645 = qJDD(3) + qJDD(5);
t556 = m(6) * t560 + mrSges(6,1) * t645 - mrSges(6,3) * t574 - t579 * t594 + t588 * t646;
t561 = t562 * t658 + t563 * t661;
t573 = -qJD(5) * t594 + t600 * t661 - t601 * t658;
t589 = mrSges(6,1) * t646 - mrSges(6,3) * t594;
t557 = m(6) * t561 - mrSges(6,2) * t645 + mrSges(6,3) * t573 + t579 * t593 - t589 * t646;
t546 = t661 * t556 + t658 * t557;
t596 = -mrSges(5,1) * t621 + mrSges(5,2) * t622;
t609 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t621;
t543 = m(5) * t566 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t601 + qJD(3) * t609 - t596 * t622 + t546;
t610 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t622;
t675 = -t556 * t658 + t661 * t557;
t544 = m(5) * t567 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t600 - qJD(3) * t610 + t596 * t621 + t675;
t539 = t657 * t543 + t656 * t544;
t632 = (mrSges(4,1) * t659 + mrSges(4,2) * t662) * qJD(1);
t683 = qJD(1) * t659;
t636 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t683;
t536 = m(4) * t603 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t634 + qJD(3) * t636 - t632 * t682 + t539;
t638 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t682;
t676 = -t543 * t656 + t657 * t544;
t537 = m(4) * t604 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t633 - qJD(3) * t638 - t632 * t683 + t676;
t532 = t662 * t536 + t659 * t537;
t619 = -qJDD(1) * pkin(1) + t672;
t670 = -m(3) * t619 + (t664 * mrSges(3,3)) - t532;
t528 = m(2) * t639 - (t664 * mrSges(2,2)) + t687 * qJDD(1) + t670;
t615 = pkin(1) * t664 - t673;
t586 = -t633 * pkin(3) + qJDD(4) + t637 * t682 + (-qJ(4) * t653 + t688) * t664 + t673;
t568 = -t600 * pkin(4) - t620 * pkin(7) + t622 * t611 + t586;
t671 = m(6) * t568 - mrSges(6,1) * t573 + t574 * mrSges(6,2) - t588 * t593 + t594 * t589;
t558 = m(5) * t586 - mrSges(5,1) * t600 + t601 * mrSges(5,2) - t609 * t621 + t622 * t610 + t671;
t608 = t688 * t664 + t673;
t667 = -m(4) * t608 + mrSges(4,1) * t633 - t634 * mrSges(4,2) - t636 * t683 - t638 * t682 - t558;
t666 = -m(3) * t615 + (t664 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t667;
t551 = m(2) * t640 - (mrSges(2,1) * t664) - qJDD(1) * mrSges(2,2) + t666;
t684 = t663 * t528 + t660 * t551;
t678 = -t528 * t660 + t663 * t551;
t677 = -t659 * t536 + t662 * t537;
t576 = Ifges(6,4) * t594 + Ifges(6,2) * t593 + Ifges(6,6) * t646;
t577 = Ifges(6,1) * t594 + Ifges(6,4) * t593 + Ifges(6,5) * t646;
t669 = mrSges(6,1) * t560 - mrSges(6,2) * t561 + Ifges(6,5) * t574 + Ifges(6,6) * t573 + Ifges(6,3) * t645 + t594 * t576 - t593 * t577;
t575 = Ifges(6,5) * t594 + Ifges(6,6) * t593 + Ifges(6,3) * t646;
t547 = -mrSges(6,1) * t568 + mrSges(6,3) * t561 + Ifges(6,4) * t574 + Ifges(6,2) * t573 + Ifges(6,6) * t645 - t575 * t594 + t577 * t646;
t548 = mrSges(6,2) * t568 - mrSges(6,3) * t560 + Ifges(6,1) * t574 + Ifges(6,4) * t573 + Ifges(6,5) * t645 + t575 * t593 - t576 * t646;
t590 = Ifges(5,5) * t622 + Ifges(5,6) * t621 + (Ifges(5,3) * qJD(3));
t592 = Ifges(5,1) * t622 + Ifges(5,4) * t621 + Ifges(5,5) * qJD(3);
t533 = -mrSges(5,1) * t586 + mrSges(5,3) * t567 + Ifges(5,4) * t601 + Ifges(5,2) * t600 + Ifges(5,6) * qJDD(3) - pkin(4) * t671 + pkin(7) * t675 + qJD(3) * t592 + t661 * t547 + t658 * t548 - t622 * t590;
t591 = Ifges(5,4) * t622 + Ifges(5,2) * t621 + Ifges(5,6) * qJD(3);
t534 = mrSges(5,2) * t586 - mrSges(5,3) * t566 + Ifges(5,1) * t601 + Ifges(5,4) * t600 + Ifges(5,5) * qJDD(3) - pkin(7) * t546 - qJD(3) * t591 - t547 * t658 + t548 * t661 + t590 * t621;
t623 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t662 - Ifges(4,6) * t659) * qJD(1);
t625 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t662 - Ifges(4,4) * t659) * qJD(1);
t524 = -mrSges(4,1) * t608 + mrSges(4,3) * t604 + Ifges(4,4) * t634 + Ifges(4,2) * t633 + Ifges(4,6) * qJDD(3) - pkin(3) * t558 + qJ(4) * t676 + qJD(3) * t625 + t657 * t533 + t656 * t534 - t623 * t682;
t624 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t662 - Ifges(4,2) * t659) * qJD(1);
t526 = mrSges(4,2) * t608 - mrSges(4,3) * t603 + Ifges(4,1) * t634 + Ifges(4,4) * t633 + Ifges(4,5) * qJDD(3) - qJ(4) * t539 - qJD(3) * t624 - t533 * t656 + t534 * t657 - t623 * t683;
t530 = qJDD(1) * mrSges(3,2) - t670;
t668 = mrSges(2,1) * t639 - mrSges(2,2) * t640 + mrSges(3,2) * t619 - mrSges(3,3) * t615 - pkin(1) * t530 - pkin(6) * t532 + qJ(2) * t666 - t524 * t659 + t662 * t526 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t665 = mrSges(4,1) * t603 + mrSges(5,1) * t566 - mrSges(4,2) * t604 - mrSges(5,2) * t567 + Ifges(4,5) * t634 + Ifges(5,5) * t601 + Ifges(4,6) * t633 + Ifges(5,6) * t600 + pkin(3) * t539 + pkin(4) * t546 + t622 * t591 - t621 * t592 + t624 * t682 + t625 * t683 + t669 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3);
t531 = -m(3) * g(3) + t677;
t523 = (t685 * t664) + t665 + pkin(2) * t532 - mrSges(2,3) * t639 + mrSges(3,1) * t619 - qJ(2) * t531 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + t686 * qJDD(1);
t522 = -mrSges(3,1) * t615 + mrSges(2,3) * t640 - pkin(1) * t531 - pkin(2) * t667 - pkin(6) * t677 + t687 * g(3) - t685 * qJDD(1) - t662 * t524 - t659 * t526 + t686 * t664;
t1 = [-m(1) * g(1) + t678; -m(1) * g(2) + t684; (-m(1) - m(2) - m(3)) * g(3) + t677; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t684 - t660 * t522 + t663 * t523; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t678 + t663 * t522 + t660 * t523; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t668; t668; t530; t665; t558; t669;];
tauJB = t1;
