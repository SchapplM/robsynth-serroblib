% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRRR10
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRRR10_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR10_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR10_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR10_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR10_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR10_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR10_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:33:03
% EndTime: 2019-12-31 22:33:16
% DurationCPUTime: 13.09s
% Computational Cost: add. (212740->325), mult. (456802->425), div. (0->0), fcn. (360449->12), ass. (0->136)
t609 = cos(pkin(5));
t640 = t609 * g(3);
t608 = sin(pkin(5));
t613 = sin(qJ(2));
t639 = t608 * t613;
t618 = cos(qJ(2));
t638 = t608 * t618;
t637 = t609 * t613;
t636 = t609 * t618;
t614 = sin(qJ(1));
t619 = cos(qJ(1));
t600 = t614 * g(1) - t619 * g(2);
t620 = qJD(1) ^ 2;
t591 = t620 * t608 * pkin(7) + qJDD(1) * pkin(1) + t600;
t601 = -t619 * g(1) - t614 * g(2);
t631 = qJDD(1) * t608;
t592 = -t620 * pkin(1) + pkin(7) * t631 + t601;
t634 = t591 * t637 + t618 * t592;
t568 = -g(3) * t639 + t634;
t605 = t609 * qJD(1) + qJD(2);
t633 = qJD(1) * t608;
t630 = t613 * t633;
t589 = t605 * mrSges(3,1) - mrSges(3,3) * t630;
t593 = (-mrSges(3,1) * t618 + mrSges(3,2) * t613) * t633;
t596 = -qJD(2) * t630 + t618 * t631;
t604 = t609 * qJDD(1) + qJDD(2);
t594 = (-pkin(2) * t618 - pkin(8) * t613) * t633;
t603 = t605 ^ 2;
t632 = qJD(1) * t618;
t553 = -t603 * pkin(2) + t604 * pkin(8) + (-g(3) * t613 + t594 * t632) * t608 + t634;
t595 = (qJD(2) * t632 + qJDD(1) * t613) * t608;
t554 = -t596 * pkin(2) - t595 * pkin(8) - t640 + (-t591 + (pkin(2) * t613 - pkin(8) * t618) * t605 * qJD(1)) * t608;
t612 = sin(qJ(3));
t617 = cos(qJ(3));
t530 = -t612 * t553 + t617 * t554;
t583 = t617 * t605 - t612 * t630;
t566 = t583 * qJD(3) + t617 * t595 + t612 * t604;
t584 = t612 * t605 + t617 * t630;
t588 = qJDD(3) - t596;
t629 = t608 * t632;
t599 = qJD(3) - t629;
t523 = (t583 * t599 - t566) * pkin(9) + (t583 * t584 + t588) * pkin(3) + t530;
t531 = t617 * t553 + t612 * t554;
t565 = -t584 * qJD(3) - t612 * t595 + t617 * t604;
t575 = t599 * pkin(3) - t584 * pkin(9);
t582 = t583 ^ 2;
t525 = -t582 * pkin(3) + t565 * pkin(9) - t599 * t575 + t531;
t611 = sin(qJ(4));
t616 = cos(qJ(4));
t521 = t611 * t523 + t616 * t525;
t571 = t611 * t583 + t616 * t584;
t538 = -t571 * qJD(4) + t616 * t565 - t611 * t566;
t570 = t616 * t583 - t611 * t584;
t547 = -t570 * mrSges(5,1) + t571 * mrSges(5,2);
t598 = qJD(4) + t599;
t558 = t598 * mrSges(5,1) - t571 * mrSges(5,3);
t587 = qJDD(4) + t588;
t548 = -t570 * pkin(4) - t571 * pkin(10);
t597 = t598 ^ 2;
t518 = -t597 * pkin(4) + t587 * pkin(10) + t570 * t548 + t521;
t567 = -g(3) * t638 + t591 * t636 - t613 * t592;
t552 = -t604 * pkin(2) - t603 * pkin(8) + t594 * t630 - t567;
t529 = -t565 * pkin(3) - t582 * pkin(9) + t584 * t575 + t552;
t539 = t570 * qJD(4) + t611 * t565 + t616 * t566;
t519 = (-t570 * t598 - t539) * pkin(10) + (t571 * t598 - t538) * pkin(4) + t529;
t610 = sin(qJ(5));
t615 = cos(qJ(5));
t515 = -t610 * t518 + t615 * t519;
t555 = -t610 * t571 + t615 * t598;
t528 = t555 * qJD(5) + t615 * t539 + t610 * t587;
t537 = qJDD(5) - t538;
t556 = t615 * t571 + t610 * t598;
t540 = -t555 * mrSges(6,1) + t556 * mrSges(6,2);
t569 = qJD(5) - t570;
t541 = -t569 * mrSges(6,2) + t555 * mrSges(6,3);
t513 = m(6) * t515 + t537 * mrSges(6,1) - t528 * mrSges(6,3) - t556 * t540 + t569 * t541;
t516 = t615 * t518 + t610 * t519;
t527 = -t556 * qJD(5) - t610 * t539 + t615 * t587;
t542 = t569 * mrSges(6,1) - t556 * mrSges(6,3);
t514 = m(6) * t516 - t537 * mrSges(6,2) + t527 * mrSges(6,3) + t555 * t540 - t569 * t542;
t625 = -t610 * t513 + t615 * t514;
t504 = m(5) * t521 - t587 * mrSges(5,2) + t538 * mrSges(5,3) + t570 * t547 - t598 * t558 + t625;
t520 = t616 * t523 - t611 * t525;
t557 = -t598 * mrSges(5,2) + t570 * mrSges(5,3);
t517 = -t587 * pkin(4) - t597 * pkin(10) + t571 * t548 - t520;
t623 = -m(6) * t517 + t527 * mrSges(6,1) - t528 * mrSges(6,2) + t555 * t541 - t556 * t542;
t509 = m(5) * t520 + t587 * mrSges(5,1) - t539 * mrSges(5,3) - t571 * t547 + t598 * t557 + t623;
t498 = t611 * t504 + t616 * t509;
t572 = -t583 * mrSges(4,1) + t584 * mrSges(4,2);
t573 = -t599 * mrSges(4,2) + t583 * mrSges(4,3);
t496 = m(4) * t530 + t588 * mrSges(4,1) - t566 * mrSges(4,3) - t584 * t572 + t599 * t573 + t498;
t574 = t599 * mrSges(4,1) - t584 * mrSges(4,3);
t626 = t616 * t504 - t611 * t509;
t497 = m(4) * t531 - t588 * mrSges(4,2) + t565 * mrSges(4,3) + t583 * t572 - t599 * t574 + t626;
t627 = -t612 * t496 + t617 * t497;
t487 = m(3) * t568 - t604 * mrSges(3,2) + t596 * mrSges(3,3) - t605 * t589 + t593 * t629 + t627;
t490 = t617 * t496 + t612 * t497;
t579 = -t608 * t591 - t640;
t590 = -t605 * mrSges(3,2) + mrSges(3,3) * t629;
t489 = m(3) * t579 - t596 * mrSges(3,1) + t595 * mrSges(3,2) + (t589 * t613 - t590 * t618) * t633 + t490;
t505 = t615 * t513 + t610 * t514;
t622 = m(5) * t529 - t538 * mrSges(5,1) + t539 * mrSges(5,2) - t570 * t557 + t571 * t558 + t505;
t621 = -m(4) * t552 + t565 * mrSges(4,1) - t566 * mrSges(4,2) + t583 * t573 - t584 * t574 - t622;
t501 = m(3) * t567 + t604 * mrSges(3,1) - t595 * mrSges(3,3) + t605 * t590 - t593 * t630 + t621;
t477 = t487 * t637 - t608 * t489 + t501 * t636;
t475 = m(2) * t600 + qJDD(1) * mrSges(2,1) - t620 * mrSges(2,2) + t477;
t483 = t618 * t487 - t613 * t501;
t482 = m(2) * t601 - t620 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t483;
t635 = t619 * t475 + t614 * t482;
t476 = t487 * t639 + t609 * t489 + t501 * t638;
t628 = -t614 * t475 + t619 * t482;
t532 = Ifges(6,5) * t556 + Ifges(6,6) * t555 + Ifges(6,3) * t569;
t534 = Ifges(6,1) * t556 + Ifges(6,4) * t555 + Ifges(6,5) * t569;
t506 = -mrSges(6,1) * t517 + mrSges(6,3) * t516 + Ifges(6,4) * t528 + Ifges(6,2) * t527 + Ifges(6,6) * t537 - t556 * t532 + t569 * t534;
t533 = Ifges(6,4) * t556 + Ifges(6,2) * t555 + Ifges(6,6) * t569;
t507 = mrSges(6,2) * t517 - mrSges(6,3) * t515 + Ifges(6,1) * t528 + Ifges(6,4) * t527 + Ifges(6,5) * t537 + t555 * t532 - t569 * t533;
t543 = Ifges(5,5) * t571 + Ifges(5,6) * t570 + Ifges(5,3) * t598;
t544 = Ifges(5,4) * t571 + Ifges(5,2) * t570 + Ifges(5,6) * t598;
t491 = mrSges(5,2) * t529 - mrSges(5,3) * t520 + Ifges(5,1) * t539 + Ifges(5,4) * t538 + Ifges(5,5) * t587 - pkin(10) * t505 - t610 * t506 + t615 * t507 + t570 * t543 - t598 * t544;
t545 = Ifges(5,1) * t571 + Ifges(5,4) * t570 + Ifges(5,5) * t598;
t492 = -mrSges(5,1) * t529 - mrSges(6,1) * t515 + mrSges(6,2) * t516 + mrSges(5,3) * t521 + Ifges(5,4) * t539 - Ifges(6,5) * t528 + Ifges(5,2) * t538 + Ifges(5,6) * t587 - Ifges(6,6) * t527 - Ifges(6,3) * t537 - pkin(4) * t505 - t556 * t533 + t555 * t534 - t571 * t543 + t598 * t545;
t559 = Ifges(4,5) * t584 + Ifges(4,6) * t583 + Ifges(4,3) * t599;
t561 = Ifges(4,1) * t584 + Ifges(4,4) * t583 + Ifges(4,5) * t599;
t478 = -mrSges(4,1) * t552 + mrSges(4,3) * t531 + Ifges(4,4) * t566 + Ifges(4,2) * t565 + Ifges(4,6) * t588 - pkin(3) * t622 + pkin(9) * t626 + t611 * t491 + t616 * t492 - t584 * t559 + t599 * t561;
t560 = Ifges(4,4) * t584 + Ifges(4,2) * t583 + Ifges(4,6) * t599;
t479 = mrSges(4,2) * t552 - mrSges(4,3) * t530 + Ifges(4,1) * t566 + Ifges(4,4) * t565 + Ifges(4,5) * t588 - pkin(9) * t498 + t616 * t491 - t611 * t492 + t583 * t559 - t599 * t560;
t576 = Ifges(3,3) * t605 + (Ifges(3,5) * t613 + Ifges(3,6) * t618) * t633;
t577 = Ifges(3,6) * t605 + (Ifges(3,4) * t613 + Ifges(3,2) * t618) * t633;
t472 = mrSges(3,2) * t579 - mrSges(3,3) * t567 + Ifges(3,1) * t595 + Ifges(3,4) * t596 + Ifges(3,5) * t604 - pkin(8) * t490 - t612 * t478 + t617 * t479 + t576 * t629 - t605 * t577;
t578 = Ifges(3,5) * t605 + (Ifges(3,1) * t613 + Ifges(3,4) * t618) * t633;
t473 = -pkin(4) * t623 - t576 * t630 + t583 * t561 - t584 * t560 - Ifges(5,3) * t587 - Ifges(4,3) * t588 + Ifges(3,4) * t595 + Ifges(3,2) * t596 - pkin(3) * t498 + Ifges(3,6) * t604 - pkin(2) * t490 - pkin(10) * t625 - t571 * t544 - mrSges(3,1) * t579 - Ifges(4,5) * t566 + mrSges(3,3) * t568 + t570 * t545 - Ifges(5,6) * t538 - Ifges(5,5) * t539 - Ifges(4,6) * t565 - mrSges(5,1) * t520 + mrSges(5,2) * t521 - mrSges(4,1) * t530 + mrSges(4,2) * t531 + t605 * t578 - t610 * t507 - t615 * t506;
t624 = pkin(7) * t483 + t472 * t613 + t473 * t618;
t471 = Ifges(3,5) * t595 + Ifges(3,6) * t596 + Ifges(3,3) * t604 + mrSges(3,1) * t567 - mrSges(3,2) * t568 + t612 * t479 + t617 * t478 + pkin(2) * t621 + pkin(8) * t627 + (t577 * t613 - t578 * t618) * t633;
t470 = -mrSges(2,2) * g(3) - mrSges(2,3) * t600 + Ifges(2,5) * qJDD(1) - t620 * Ifges(2,6) + t618 * t472 - t613 * t473 + (-t476 * t608 - t477 * t609) * pkin(7);
t469 = mrSges(2,1) * g(3) + mrSges(2,3) * t601 + t620 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t476 - t608 * t471 + t609 * t624;
t1 = [-m(1) * g(1) + t628; -m(1) * g(2) + t635; (-m(1) - m(2)) * g(3) + t476; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t635 - t614 * t469 + t619 * t470; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t628 + t619 * t469 + t614 * t470; -mrSges(1,1) * g(2) + mrSges(2,1) * t600 + mrSges(1,2) * g(1) - mrSges(2,2) * t601 + Ifges(2,3) * qJDD(1) + pkin(1) * t477 + t609 * t471 + t608 * t624;];
tauB = t1;
