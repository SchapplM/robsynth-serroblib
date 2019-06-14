% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-05-06 04:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRRR8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR8_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR8_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR8_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR8_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR8_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 04:14:57
% EndTime: 2019-05-06 04:15:04
% DurationCPUTime: 7.37s
% Computational Cost: add. (122339->339), mult. (239441->416), div. (0->0), fcn. (164122->10), ass. (0->132)
t659 = sin(qJ(1));
t664 = cos(qJ(1));
t641 = -t664 * g(1) - t659 * g(2);
t673 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t641;
t689 = -pkin(1) - pkin(7);
t688 = mrSges(2,1) - mrSges(3,2);
t687 = Ifges(2,5) - Ifges(3,4);
t686 = (-Ifges(2,6) + Ifges(3,5));
t640 = t659 * g(1) - t664 * g(2);
t665 = qJD(1) ^ 2;
t672 = -t665 * qJ(2) + qJDD(2) - t640;
t618 = qJDD(1) * t689 + t672;
t658 = sin(qJ(3));
t663 = cos(qJ(3));
t607 = t658 * g(3) + t663 * t618;
t682 = qJD(1) * qJD(3);
t680 = t658 * t682;
t636 = t663 * qJDD(1) - t680;
t585 = (-t636 - t680) * pkin(8) + (-t658 * t663 * t665 + qJDD(3)) * pkin(3) + t607;
t608 = -t663 * g(3) + t658 * t618;
t635 = -t658 * qJDD(1) - t663 * t682;
t683 = qJD(1) * t663;
t639 = qJD(3) * pkin(3) - pkin(8) * t683;
t652 = t658 ^ 2;
t588 = -t652 * t665 * pkin(3) + t635 * pkin(8) - qJD(3) * t639 + t608;
t657 = sin(qJ(4));
t662 = cos(qJ(4));
t571 = t657 * t585 + t662 * t588;
t628 = (-t657 * t658 + t662 * t663) * qJD(1);
t596 = -t628 * qJD(4) + t662 * t635 - t657 * t636;
t684 = qJD(1) * t658;
t627 = -t657 * t683 - t662 * t684;
t605 = -t627 * mrSges(5,1) + t628 * mrSges(5,2);
t649 = qJD(3) + qJD(4);
t616 = t649 * mrSges(5,1) - t628 * mrSges(5,3);
t648 = qJDD(3) + qJDD(4);
t591 = -t635 * pkin(3) + t639 * t683 + (-pkin(8) * t652 + t689) * t665 + t673;
t597 = t627 * qJD(4) + t657 * t635 + t662 * t636;
t562 = (-t627 * t649 - t597) * pkin(9) + (t628 * t649 - t596) * pkin(4) + t591;
t606 = -t627 * pkin(4) - t628 * pkin(9);
t647 = t649 ^ 2;
t565 = -t647 * pkin(4) + t648 * pkin(9) + t627 * t606 + t571;
t656 = sin(qJ(5));
t661 = cos(qJ(5));
t554 = t661 * t562 - t656 * t565;
t610 = -t656 * t628 + t661 * t649;
t575 = t610 * qJD(5) + t661 * t597 + t656 * t648;
t595 = qJDD(5) - t596;
t611 = t661 * t628 + t656 * t649;
t623 = qJD(5) - t627;
t552 = (t610 * t623 - t575) * pkin(10) + (t610 * t611 + t595) * pkin(5) + t554;
t555 = t656 * t562 + t661 * t565;
t574 = -t611 * qJD(5) - t656 * t597 + t661 * t648;
t600 = t623 * pkin(5) - t611 * pkin(10);
t609 = t610 ^ 2;
t553 = -t609 * pkin(5) + t574 * pkin(10) - t623 * t600 + t555;
t655 = sin(qJ(6));
t660 = cos(qJ(6));
t550 = t660 * t552 - t655 * t553;
t586 = t660 * t610 - t655 * t611;
t559 = t586 * qJD(6) + t655 * t574 + t660 * t575;
t587 = t655 * t610 + t660 * t611;
t572 = -t586 * mrSges(7,1) + t587 * mrSges(7,2);
t620 = qJD(6) + t623;
t576 = -t620 * mrSges(7,2) + t586 * mrSges(7,3);
t592 = qJDD(6) + t595;
t548 = m(7) * t550 + t592 * mrSges(7,1) - t559 * mrSges(7,3) - t587 * t572 + t620 * t576;
t551 = t655 * t552 + t660 * t553;
t558 = -t587 * qJD(6) + t660 * t574 - t655 * t575;
t577 = t620 * mrSges(7,1) - t587 * mrSges(7,3);
t549 = m(7) * t551 - t592 * mrSges(7,2) + t558 * mrSges(7,3) + t586 * t572 - t620 * t577;
t540 = t660 * t548 + t655 * t549;
t589 = -t610 * mrSges(6,1) + t611 * mrSges(6,2);
t598 = -t623 * mrSges(6,2) + t610 * mrSges(6,3);
t538 = m(6) * t554 + t595 * mrSges(6,1) - t575 * mrSges(6,3) - t611 * t589 + t623 * t598 + t540;
t599 = t623 * mrSges(6,1) - t611 * mrSges(6,3);
t675 = -t655 * t548 + t660 * t549;
t539 = m(6) * t555 - t595 * mrSges(6,2) + t574 * mrSges(6,3) + t610 * t589 - t623 * t599 + t675;
t676 = -t656 * t538 + t661 * t539;
t533 = m(5) * t571 - t648 * mrSges(5,2) + t596 * mrSges(5,3) + t627 * t605 - t649 * t616 + t676;
t570 = t662 * t585 - t657 * t588;
t615 = -t649 * mrSges(5,2) + t627 * mrSges(5,3);
t564 = -t648 * pkin(4) - t647 * pkin(9) + t628 * t606 - t570;
t556 = -t574 * pkin(5) - t609 * pkin(10) + t611 * t600 + t564;
t670 = m(7) * t556 - t558 * mrSges(7,1) + t559 * mrSges(7,2) - t586 * t576 + t587 * t577;
t667 = -m(6) * t564 + t574 * mrSges(6,1) - t575 * mrSges(6,2) + t610 * t598 - t611 * t599 - t670;
t544 = m(5) * t570 + t648 * mrSges(5,1) - t597 * mrSges(5,3) - t628 * t605 + t649 * t615 + t667;
t526 = t657 * t533 + t662 * t544;
t634 = (mrSges(4,1) * t658 + mrSges(4,2) * t663) * qJD(1);
t637 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t684;
t523 = m(4) * t607 + qJDD(3) * mrSges(4,1) - t636 * mrSges(4,3) + qJD(3) * t637 - t634 * t683 + t526;
t638 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t683;
t677 = t662 * t533 - t657 * t544;
t524 = m(4) * t608 - qJDD(3) * mrSges(4,2) + t635 * mrSges(4,3) - qJD(3) * t638 - t634 * t684 + t677;
t520 = t663 * t523 + t658 * t524;
t622 = -qJDD(1) * pkin(1) + t672;
t671 = -m(3) * t622 + (t665 * mrSges(3,3)) - t520;
t518 = m(2) * t640 - (t665 * mrSges(2,2)) + qJDD(1) * t688 + t671;
t619 = t665 * pkin(1) - t673;
t617 = t665 * t689 + t673;
t534 = t661 * t538 + t656 * t539;
t669 = m(5) * t591 - t596 * mrSges(5,1) + t597 * mrSges(5,2) - t627 * t615 + t628 * t616 + t534;
t668 = -m(4) * t617 + t635 * mrSges(4,1) - t636 * mrSges(4,2) - t637 * t684 - t638 * t683 - t669;
t666 = -m(3) * t619 + (t665 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t668;
t530 = m(2) * t641 - (t665 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t666;
t685 = t664 * t518 + t659 * t530;
t679 = -t659 * t518 + t664 * t530;
t678 = -t658 * t523 + t663 * t524;
t626 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t663 - Ifges(4,4) * t658) * qJD(1);
t625 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t663 - Ifges(4,2) * t658) * qJD(1);
t624 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t663 - Ifges(4,6) * t658) * qJD(1);
t603 = Ifges(5,1) * t628 + Ifges(5,4) * t627 + Ifges(5,5) * t649;
t602 = Ifges(5,4) * t628 + Ifges(5,2) * t627 + Ifges(5,6) * t649;
t601 = Ifges(5,5) * t628 + Ifges(5,6) * t627 + Ifges(5,3) * t649;
t580 = Ifges(6,1) * t611 + Ifges(6,4) * t610 + Ifges(6,5) * t623;
t579 = Ifges(6,4) * t611 + Ifges(6,2) * t610 + Ifges(6,6) * t623;
t578 = Ifges(6,5) * t611 + Ifges(6,6) * t610 + Ifges(6,3) * t623;
t568 = Ifges(7,1) * t587 + Ifges(7,4) * t586 + Ifges(7,5) * t620;
t567 = Ifges(7,4) * t587 + Ifges(7,2) * t586 + Ifges(7,6) * t620;
t566 = Ifges(7,5) * t587 + Ifges(7,6) * t586 + Ifges(7,3) * t620;
t542 = mrSges(7,2) * t556 - mrSges(7,3) * t550 + Ifges(7,1) * t559 + Ifges(7,4) * t558 + Ifges(7,5) * t592 + t586 * t566 - t620 * t567;
t541 = -mrSges(7,1) * t556 + mrSges(7,3) * t551 + Ifges(7,4) * t559 + Ifges(7,2) * t558 + Ifges(7,6) * t592 - t587 * t566 + t620 * t568;
t527 = mrSges(6,2) * t564 - mrSges(6,3) * t554 + Ifges(6,1) * t575 + Ifges(6,4) * t574 + Ifges(6,5) * t595 - pkin(10) * t540 - t655 * t541 + t660 * t542 + t610 * t578 - t623 * t579;
t525 = -mrSges(6,1) * t564 + mrSges(6,3) * t555 + Ifges(6,4) * t575 + Ifges(6,2) * t574 + Ifges(6,6) * t595 - pkin(5) * t670 + pkin(10) * t675 + t660 * t541 + t655 * t542 - t611 * t578 + t623 * t580;
t521 = Ifges(5,4) * t597 + Ifges(5,2) * t596 + Ifges(5,6) * t648 - t628 * t601 + t649 * t603 - mrSges(5,1) * t591 + mrSges(5,3) * t571 - Ifges(6,5) * t575 - Ifges(6,6) * t574 - Ifges(6,3) * t595 - t611 * t579 + t610 * t580 - mrSges(6,1) * t554 + mrSges(6,2) * t555 - Ifges(7,5) * t559 - Ifges(7,6) * t558 - Ifges(7,3) * t592 - t587 * t567 + t586 * t568 - mrSges(7,1) * t550 + mrSges(7,2) * t551 - pkin(5) * t540 - pkin(4) * t534;
t519 = -m(3) * g(3) + t678;
t516 = mrSges(5,2) * t591 - mrSges(5,3) * t570 + Ifges(5,1) * t597 + Ifges(5,4) * t596 + Ifges(5,5) * t648 - pkin(9) * t534 - t656 * t525 + t661 * t527 + t627 * t601 - t649 * t602;
t515 = mrSges(4,2) * t617 - mrSges(4,3) * t607 + Ifges(4,1) * t636 + Ifges(4,4) * t635 + Ifges(4,5) * qJDD(3) - pkin(8) * t526 - qJD(3) * t625 + t662 * t516 - t657 * t521 - t624 * t684;
t514 = -mrSges(4,1) * t617 + mrSges(4,3) * t608 + Ifges(4,4) * t636 + Ifges(4,2) * t635 + Ifges(4,6) * qJDD(3) - pkin(3) * t669 + pkin(8) * t677 + qJD(3) * t626 + t657 * t516 + t662 * t521 - t624 * t683;
t513 = (t686 * t665) + t687 * qJDD(1) + pkin(4) * t667 - qJ(2) * t519 + pkin(9) * t676 + pkin(2) * t520 + Ifges(4,3) * qJDD(3) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + (t663 * t625 + t658 * t626) * qJD(1) + t661 * t525 + t656 * t527 + Ifges(4,5) * t636 - mrSges(2,3) * t640 + Ifges(5,3) * t648 - t627 * t603 + t628 * t602 + Ifges(4,6) * t635 + mrSges(3,1) * t622 + mrSges(4,1) * t607 - mrSges(4,2) * t608 + Ifges(5,6) * t596 + Ifges(5,5) * t597 + mrSges(5,1) * t570 - mrSges(5,2) * t571 + pkin(3) * t526;
t512 = -mrSges(3,1) * t619 + mrSges(2,3) * t641 - pkin(1) * t519 - pkin(2) * t668 - pkin(7) * t678 + g(3) * t688 - qJDD(1) * t686 - t663 * t514 - t658 * t515 + t665 * t687;
t1 = [-m(1) * g(1) + t679; -m(1) * g(2) + t685; (-m(1) - m(2) - m(3)) * g(3) + t678; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t685 - t659 * t512 + t664 * t513; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t679 + t664 * t512 + t659 * t513; pkin(1) * t671 + qJ(2) * t666 + t663 * t515 - t658 * t514 - pkin(7) * t520 + mrSges(2,1) * t640 - mrSges(2,2) * t641 + mrSges(3,2) * t622 - mrSges(3,3) * t619 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;
