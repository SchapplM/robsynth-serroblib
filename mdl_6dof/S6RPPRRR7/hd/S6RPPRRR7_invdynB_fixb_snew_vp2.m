% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 16:07
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPPRRR7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR7_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR7_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR7_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR7_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR7_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:05:45
% EndTime: 2019-05-05 16:05:52
% DurationCPUTime: 6.82s
% Computational Cost: add. (102006->314), mult. (234537->383), div. (0->0), fcn. (173510->10), ass. (0->131)
t650 = sin(qJ(1));
t654 = cos(qJ(1));
t626 = t650 * g(1) - t654 * g(2);
t655 = qJD(1) ^ 2;
t663 = -t655 * qJ(2) + qJDD(2) - t626;
t685 = -pkin(1) - qJ(3);
t691 = -(2 * qJD(1) * qJD(3)) + t685 * qJDD(1) + t663;
t645 = sin(pkin(10));
t638 = t645 ^ 2;
t646 = cos(pkin(10));
t682 = t646 ^ 2 + t638;
t676 = t682 * mrSges(4,3);
t627 = -t654 * g(1) - t650 * g(2);
t690 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t627;
t689 = pkin(3) * t655;
t688 = mrSges(2,1) - mrSges(3,2);
t687 = -Ifges(3,4) + Ifges(2,5);
t686 = -Ifges(2,6) + Ifges(3,5);
t684 = mrSges(4,2) * t646;
t607 = t645 * g(3) + t691 * t646;
t589 = (-pkin(7) * qJDD(1) - t645 * t689) * t646 + t607;
t608 = -t646 * g(3) + t691 * t645;
t678 = qJDD(1) * t645;
t592 = -pkin(7) * t678 - t638 * t689 + t608;
t649 = sin(qJ(4));
t653 = cos(qJ(4));
t575 = t653 * t589 - t649 * t592;
t666 = -t645 * t649 + t646 * t653;
t667 = -t645 * t653 - t646 * t649;
t622 = t667 * qJD(1);
t680 = t622 * qJD(4);
t610 = t666 * qJDD(1) + t680;
t623 = t666 * qJD(1);
t559 = (-t610 + t680) * pkin(8) + (t622 * t623 + qJDD(4)) * pkin(4) + t575;
t576 = t649 * t589 + t653 * t592;
t609 = -t623 * qJD(4) + t667 * qJDD(1);
t617 = qJD(4) * pkin(4) - t623 * pkin(8);
t621 = t622 ^ 2;
t561 = -t621 * pkin(4) + t609 * pkin(8) - qJD(4) * t617 + t576;
t648 = sin(qJ(5));
t652 = cos(qJ(5));
t557 = t648 * t559 + t652 * t561;
t602 = t648 * t622 + t652 * t623;
t573 = -t602 * qJD(5) + t652 * t609 - t648 * t610;
t601 = t652 * t622 - t648 * t623;
t584 = -t601 * mrSges(6,1) + t602 * mrSges(6,2);
t640 = qJD(4) + qJD(5);
t594 = t640 * mrSges(6,1) - t602 * mrSges(6,3);
t637 = qJDD(4) + qJDD(5);
t585 = -t601 * pkin(5) - t602 * pkin(9);
t636 = t640 ^ 2;
t554 = -t636 * pkin(5) + t637 * pkin(9) + t601 * t585 + t557;
t662 = qJDD(3) + t690;
t596 = pkin(3) * t678 + (-t682 * pkin(7) + t685) * t655 + t662;
t570 = -t609 * pkin(4) - t621 * pkin(8) + t623 * t617 + t596;
t574 = t601 * qJD(5) + t648 * t609 + t652 * t610;
t555 = (-t601 * t640 - t574) * pkin(9) + (t602 * t640 - t573) * pkin(5) + t570;
t647 = sin(qJ(6));
t651 = cos(qJ(6));
t551 = -t647 * t554 + t651 * t555;
t590 = -t647 * t602 + t651 * t640;
t564 = t590 * qJD(6) + t651 * t574 + t647 * t637;
t572 = qJDD(6) - t573;
t591 = t651 * t602 + t647 * t640;
t577 = -t590 * mrSges(7,1) + t591 * mrSges(7,2);
t597 = qJD(6) - t601;
t578 = -t597 * mrSges(7,2) + t590 * mrSges(7,3);
t549 = m(7) * t551 + t572 * mrSges(7,1) - t564 * mrSges(7,3) - t591 * t577 + t597 * t578;
t552 = t651 * t554 + t647 * t555;
t563 = -t591 * qJD(6) - t647 * t574 + t651 * t637;
t579 = t597 * mrSges(7,1) - t591 * mrSges(7,3);
t550 = m(7) * t552 - t572 * mrSges(7,2) + t563 * mrSges(7,3) + t590 * t577 - t597 * t579;
t671 = -t647 * t549 + t651 * t550;
t540 = m(6) * t557 - t637 * mrSges(6,2) + t573 * mrSges(6,3) + t601 * t584 - t640 * t594 + t671;
t556 = t652 * t559 - t648 * t561;
t593 = -t640 * mrSges(6,2) + t601 * mrSges(6,3);
t553 = -t637 * pkin(5) - t636 * pkin(9) + t602 * t585 - t556;
t660 = -m(7) * t553 + t563 * mrSges(7,1) - t564 * mrSges(7,2) + t590 * t578 - t591 * t579;
t545 = m(6) * t556 + t637 * mrSges(6,1) - t574 * mrSges(6,3) - t602 * t584 + t640 * t593 + t660;
t534 = t648 * t540 + t652 * t545;
t605 = -t622 * mrSges(5,1) + t623 * mrSges(5,2);
t615 = -qJD(4) * mrSges(5,2) + t622 * mrSges(5,3);
t532 = m(5) * t575 + qJDD(4) * mrSges(5,1) - t610 * mrSges(5,3) + qJD(4) * t615 - t623 * t605 + t534;
t616 = qJD(4) * mrSges(5,1) - t623 * mrSges(5,3);
t672 = t652 * t540 - t648 * t545;
t533 = m(5) * t576 - qJDD(4) * mrSges(5,2) + t609 * mrSges(5,3) - qJD(4) * t616 + t622 * t605 + t672;
t526 = t653 * t532 + t649 * t533;
t665 = -mrSges(4,3) * qJDD(1) - t655 * (mrSges(4,1) * t645 + t684);
t524 = m(4) * t607 + t665 * t646 + t526;
t673 = -t649 * t532 + t653 * t533;
t525 = m(4) * t608 + t665 * t645 + t673;
t521 = t646 * t524 + t645 * t525;
t620 = -qJDD(1) * pkin(1) + t663;
t661 = -m(3) * t620 + t655 * mrSges(3,3) - t521;
t519 = m(2) * t626 - t655 * mrSges(2,2) + t688 * qJDD(1) + t661;
t619 = t655 * pkin(1) - t690;
t614 = t685 * t655 + t662;
t541 = t651 * t549 + t647 * t550;
t659 = m(6) * t570 - t573 * mrSges(6,1) + t574 * mrSges(6,2) - t601 * t593 + t602 * t594 + t541;
t658 = m(5) * t596 - t609 * mrSges(5,1) + t610 * mrSges(5,2) - t622 * t615 + t623 * t616 + t659;
t657 = -m(4) * t614 - mrSges(4,1) * t678 - qJDD(1) * t684 - t658;
t656 = -m(3) * t619 + t655 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t657;
t537 = t656 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - t676) * t655 + m(2) * t627;
t683 = t654 * t519 + t650 * t537;
t668 = Ifges(4,5) * t646 - Ifges(4,6) * t645;
t681 = t655 * t668;
t675 = -t650 * t519 + t654 * t537;
t674 = -t645 * t524 + t646 * t525;
t670 = Ifges(4,1) * t646 - Ifges(4,4) * t645;
t669 = Ifges(4,4) * t646 - Ifges(4,2) * t645;
t600 = Ifges(5,1) * t623 + Ifges(5,4) * t622 + Ifges(5,5) * qJD(4);
t599 = Ifges(5,4) * t623 + Ifges(5,2) * t622 + Ifges(5,6) * qJD(4);
t598 = Ifges(5,5) * t623 + Ifges(5,6) * t622 + Ifges(5,3) * qJD(4);
t582 = Ifges(6,1) * t602 + Ifges(6,4) * t601 + Ifges(6,5) * t640;
t581 = Ifges(6,4) * t602 + Ifges(6,2) * t601 + Ifges(6,6) * t640;
t580 = Ifges(6,5) * t602 + Ifges(6,6) * t601 + Ifges(6,3) * t640;
t567 = Ifges(7,1) * t591 + Ifges(7,4) * t590 + Ifges(7,5) * t597;
t566 = Ifges(7,4) * t591 + Ifges(7,2) * t590 + Ifges(7,6) * t597;
t565 = Ifges(7,5) * t591 + Ifges(7,6) * t590 + Ifges(7,3) * t597;
t543 = mrSges(7,2) * t553 - mrSges(7,3) * t551 + Ifges(7,1) * t564 + Ifges(7,4) * t563 + Ifges(7,5) * t572 + t590 * t565 - t597 * t566;
t542 = -mrSges(7,1) * t553 + mrSges(7,3) * t552 + Ifges(7,4) * t564 + Ifges(7,2) * t563 + Ifges(7,6) * t572 - t591 * t565 + t597 * t567;
t528 = -mrSges(6,1) * t570 - mrSges(7,1) * t551 + mrSges(7,2) * t552 + mrSges(6,3) * t557 + Ifges(6,4) * t574 - Ifges(7,5) * t564 + Ifges(6,2) * t573 + Ifges(6,6) * t637 - Ifges(7,6) * t563 - Ifges(7,3) * t572 - pkin(5) * t541 - t591 * t566 + t590 * t567 - t602 * t580 + t640 * t582;
t527 = mrSges(6,2) * t570 - mrSges(6,3) * t556 + Ifges(6,1) * t574 + Ifges(6,4) * t573 + Ifges(6,5) * t637 - pkin(9) * t541 - t647 * t542 + t651 * t543 + t601 * t580 - t640 * t581;
t522 = mrSges(5,2) * t596 - mrSges(5,3) * t575 + Ifges(5,1) * t610 + Ifges(5,4) * t609 + Ifges(5,5) * qJDD(4) - pkin(8) * t534 - qJD(4) * t599 + t652 * t527 - t648 * t528 + t622 * t598;
t520 = -m(3) * g(3) + t674;
t517 = -mrSges(5,1) * t596 + mrSges(5,3) * t576 + Ifges(5,4) * t610 + Ifges(5,2) * t609 + Ifges(5,6) * qJDD(4) - pkin(4) * t659 + pkin(8) * t672 + qJD(4) * t600 + t648 * t527 + t652 * t528 - t623 * t598;
t516 = mrSges(4,2) * t614 - mrSges(4,3) * t607 - pkin(7) * t526 + t670 * qJDD(1) - t649 * t517 + t653 * t522 - t645 * t681;
t515 = -mrSges(4,1) * t614 + mrSges(4,3) * t608 - pkin(3) * t658 + pkin(7) * t673 + t669 * qJDD(1) + t653 * t517 + t649 * t522 - t646 * t681;
t514 = (t668 + t687) * qJDD(1) + pkin(3) * t526 + pkin(5) * t660 + pkin(9) * t671 + Ifges(5,3) * qJDD(4) + pkin(2) * t521 - qJ(2) * t520 + t647 * t543 + t651 * t542 + Ifges(6,3) * t637 + t623 * t599 - mrSges(2,3) * t626 + mrSges(3,1) * t620 - t622 * t600 + t602 * t581 + mrSges(4,1) * t607 - mrSges(4,2) * t608 + Ifges(5,6) * t609 + Ifges(5,5) * t610 - t601 * t582 + Ifges(6,5) * t574 + mrSges(5,1) * t575 - mrSges(5,2) * t576 + Ifges(6,6) * t573 + mrSges(6,1) * t556 - mrSges(6,2) * t557 + (t645 * t670 + t646 * t669 + t686) * t655 + pkin(4) * t534 + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t513 = mrSges(2,3) * t627 - mrSges(3,1) * t619 - t645 * t516 - t646 * t515 - pkin(2) * t657 - qJ(3) * t674 - pkin(1) * t520 - t686 * qJDD(1) + t688 * g(3) + (-pkin(2) * t676 + t687) * t655;
t1 = [-m(1) * g(1) + t675; -m(1) * g(2) + t683; (-m(1) - m(2) - m(3)) * g(3) + t674; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t683 - t650 * t513 + t654 * t514; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t675 + t654 * t513 + t650 * t514; pkin(1) * t661 + qJ(2) * (-t655 * t676 + t656) + t646 * t516 - t645 * t515 - qJ(3) * t521 + mrSges(2,1) * t626 - mrSges(2,2) * t627 + mrSges(3,2) * t620 - mrSges(3,3) * t619 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;
