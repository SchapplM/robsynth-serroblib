% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRPRPR6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
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
% Datum: 2019-05-04 23:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRPRPR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR6_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR6_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR6_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR6_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR6_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:08:11
% EndTime: 2019-05-04 23:08:20
% DurationCPUTime: 7.18s
% Computational Cost: add. (94372->291), mult. (186971->367), div. (0->0), fcn. (123079->12), ass. (0->124)
t584 = sin(pkin(10));
t587 = cos(pkin(10));
t571 = t584 * g(1) - t587 * g(2);
t572 = -t587 * g(1) - t584 * g(2);
t580 = -g(3) + qJDD(1);
t594 = cos(qJ(2));
t588 = cos(pkin(6));
t591 = sin(qJ(2));
t617 = t588 * t591;
t585 = sin(pkin(6));
t618 = t585 * t591;
t533 = t571 * t617 + t594 * t572 + t580 * t618;
t624 = -qJDD(2) * qJ(3) - (2 * qJD(3) * qJD(2)) - t533;
t532 = -t591 * t572 + (t571 * t588 + t580 * t585) * t594;
t623 = -pkin(2) - pkin(8);
t622 = mrSges(3,1) - mrSges(4,2);
t621 = (-Ifges(4,4) + Ifges(3,5));
t620 = Ifges(4,5) - Ifges(3,6);
t596 = qJD(2) ^ 2;
t599 = -t596 * qJ(3) + qJDD(3) - t532;
t527 = t623 * qJDD(2) + t599;
t550 = -t585 * t571 + t588 * t580;
t590 = sin(qJ(4));
t593 = cos(qJ(4));
t522 = t590 * t527 + t593 * t550;
t568 = (mrSges(5,1) * t590 + mrSges(5,2) * t593) * qJD(2);
t613 = qJD(2) * qJD(4);
t610 = t593 * t613;
t569 = t590 * qJDD(2) + t610;
t615 = qJD(2) * t593;
t574 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t615;
t567 = (pkin(4) * t590 - qJ(5) * t593) * qJD(2);
t595 = qJD(4) ^ 2;
t614 = t590 * qJD(2);
t510 = -t595 * pkin(4) + qJDD(4) * qJ(5) - t567 * t614 + t522;
t526 = t623 * t596 - t624;
t611 = t590 * t613;
t570 = t593 * qJDD(2) - t611;
t513 = (-t570 + t611) * qJ(5) + (t569 + t610) * pkin(4) + t526;
t583 = sin(pkin(11));
t586 = cos(pkin(11));
t562 = t583 * qJD(4) + t586 * t615;
t505 = -0.2e1 * qJD(5) * t562 - t583 * t510 + t586 * t513;
t548 = t583 * qJDD(4) + t586 * t570;
t561 = t586 * qJD(4) - t583 * t615;
t503 = (t561 * t614 - t548) * pkin(9) + (t561 * t562 + t569) * pkin(5) + t505;
t506 = 0.2e1 * qJD(5) * t561 + t586 * t510 + t583 * t513;
t547 = t586 * qJDD(4) - t583 * t570;
t549 = pkin(5) * t614 - t562 * pkin(9);
t560 = t561 ^ 2;
t504 = -t560 * pkin(5) + t547 * pkin(9) - t549 * t614 + t506;
t589 = sin(qJ(6));
t592 = cos(qJ(6));
t501 = t592 * t503 - t589 * t504;
t538 = t592 * t561 - t589 * t562;
t516 = t538 * qJD(6) + t589 * t547 + t592 * t548;
t539 = t589 * t561 + t592 * t562;
t523 = -t538 * mrSges(7,1) + t539 * mrSges(7,2);
t576 = qJD(6) + t614;
t530 = -t576 * mrSges(7,2) + t538 * mrSges(7,3);
t564 = qJDD(6) + t569;
t499 = m(7) * t501 + t564 * mrSges(7,1) - t516 * mrSges(7,3) - t539 * t523 + t576 * t530;
t502 = t589 * t503 + t592 * t504;
t515 = -t539 * qJD(6) + t592 * t547 - t589 * t548;
t531 = t576 * mrSges(7,1) - t539 * mrSges(7,3);
t500 = m(7) * t502 - t564 * mrSges(7,2) + t515 * mrSges(7,3) + t538 * t523 - t576 * t531;
t492 = t592 * t499 + t589 * t500;
t540 = -t561 * mrSges(6,1) + t562 * mrSges(6,2);
t545 = -mrSges(6,2) * t614 + t561 * mrSges(6,3);
t490 = m(6) * t505 + t569 * mrSges(6,1) - t548 * mrSges(6,3) - t562 * t540 + t545 * t614 + t492;
t546 = mrSges(6,1) * t614 - t562 * mrSges(6,3);
t606 = -t589 * t499 + t592 * t500;
t491 = m(6) * t506 - t569 * mrSges(6,2) + t547 * mrSges(6,3) + t561 * t540 - t546 * t614 + t606;
t607 = -t583 * t490 + t586 * t491;
t485 = m(5) * t522 - qJDD(4) * mrSges(5,2) - t569 * mrSges(5,3) - qJD(4) * t574 - t568 * t614 + t607;
t521 = t593 * t527 - t590 * t550;
t573 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t614;
t509 = -qJDD(4) * pkin(4) - t595 * qJ(5) + t567 * t615 + qJDD(5) - t521;
t507 = -t547 * pkin(5) - t560 * pkin(9) + t562 * t549 + t509;
t600 = m(7) * t507 - t515 * mrSges(7,1) + t516 * mrSges(7,2) - t538 * t530 + t539 * t531;
t597 = -m(6) * t509 + t547 * mrSges(6,1) - t548 * mrSges(6,2) + t561 * t545 - t562 * t546 - t600;
t495 = m(5) * t521 + qJDD(4) * mrSges(5,1) - t570 * mrSges(5,3) + qJD(4) * t573 - t568 * t615 + t597;
t477 = t590 * t485 + t593 * t495;
t529 = -qJDD(2) * pkin(2) + t599;
t602 = -m(4) * t529 + (t596 * mrSges(4,3)) - t477;
t473 = m(3) * t532 - (t596 * mrSges(3,2)) + qJDD(2) * t622 + t602;
t619 = t473 * t594;
t608 = t593 * t485 - t590 * t495;
t476 = m(4) * t550 + t608;
t475 = m(3) * t550 + t476;
t528 = t596 * pkin(2) + t624;
t486 = t586 * t490 + t583 * t491;
t601 = -m(5) * t526 - t569 * mrSges(5,1) - t570 * mrSges(5,2) - t573 * t614 - t574 * t615 - t486;
t598 = -m(4) * t528 + (t596 * mrSges(4,2)) + qJDD(2) * mrSges(4,3) - t601;
t483 = m(3) * t533 - (t596 * mrSges(3,1)) - qJDD(2) * mrSges(3,2) + t598;
t464 = -t585 * t475 + t483 * t617 + t588 * t619;
t462 = m(2) * t571 + t464;
t469 = -t591 * t473 + t594 * t483;
t468 = m(2) * t572 + t469;
t616 = t587 * t462 + t584 * t468;
t463 = t588 * t475 + t483 * t618 + t585 * t619;
t609 = -t584 * t462 + t587 * t468;
t517 = Ifges(7,5) * t539 + Ifges(7,6) * t538 + Ifges(7,3) * t576;
t519 = Ifges(7,1) * t539 + Ifges(7,4) * t538 + Ifges(7,5) * t576;
t493 = -mrSges(7,1) * t507 + mrSges(7,3) * t502 + Ifges(7,4) * t516 + Ifges(7,2) * t515 + Ifges(7,6) * t564 - t539 * t517 + t576 * t519;
t518 = Ifges(7,4) * t539 + Ifges(7,2) * t538 + Ifges(7,6) * t576;
t494 = mrSges(7,2) * t507 - mrSges(7,3) * t501 + Ifges(7,1) * t516 + Ifges(7,4) * t515 + Ifges(7,5) * t564 + t538 * t517 - t576 * t518;
t534 = Ifges(6,5) * t562 + Ifges(6,6) * t561 + Ifges(6,3) * t614;
t536 = Ifges(6,1) * t562 + Ifges(6,4) * t561 + Ifges(6,5) * t614;
t478 = -mrSges(6,1) * t509 + mrSges(6,3) * t506 + Ifges(6,4) * t548 + Ifges(6,2) * t547 + Ifges(6,6) * t569 - pkin(5) * t600 + pkin(9) * t606 + t592 * t493 + t589 * t494 - t562 * t534 + t536 * t614;
t535 = Ifges(6,4) * t562 + Ifges(6,2) * t561 + Ifges(6,6) * t614;
t479 = mrSges(6,2) * t509 - mrSges(6,3) * t505 + Ifges(6,1) * t548 + Ifges(6,4) * t547 + Ifges(6,5) * t569 - pkin(9) * t492 - t589 * t493 + t592 * t494 + t561 * t534 - t535 * t614;
t555 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t593 - Ifges(5,6) * t590) * qJD(2);
t556 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t593 - Ifges(5,2) * t590) * qJD(2);
t465 = mrSges(5,2) * t526 - mrSges(5,3) * t521 + Ifges(5,1) * t570 - Ifges(5,4) * t569 + Ifges(5,5) * qJDD(4) - qJ(5) * t486 - qJD(4) * t556 - t583 * t478 + t586 * t479 - t555 * t614;
t557 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t593 - Ifges(5,4) * t590) * qJD(2);
t470 = Ifges(5,4) * t570 + Ifges(5,6) * qJDD(4) - t555 * t615 + qJD(4) * t557 - mrSges(5,1) * t526 + mrSges(5,3) * t522 - Ifges(6,5) * t548 - Ifges(6,6) * t547 - t562 * t535 + t561 * t536 - mrSges(6,1) * t505 + mrSges(6,2) * t506 - Ifges(7,5) * t516 - Ifges(7,6) * t515 - Ifges(7,3) * t564 - t539 * t518 + t538 * t519 - mrSges(7,1) * t501 + mrSges(7,2) * t502 - pkin(5) * t492 - pkin(4) * t486 + (-Ifges(5,2) - Ifges(6,3)) * t569;
t459 = -mrSges(4,1) * t528 + mrSges(3,3) * t533 - pkin(2) * t476 - pkin(3) * t601 - pkin(8) * t608 - qJDD(2) * t620 - t590 * t465 - t593 * t470 - t550 * t622 + (t596 * t621);
t460 = -qJ(3) * t476 - mrSges(3,3) * t532 + pkin(3) * t477 + mrSges(4,1) * t529 + qJ(5) * t607 + t583 * t479 + t586 * t478 + pkin(4) * t597 + mrSges(5,1) * t521 - mrSges(5,2) * t522 + Ifges(5,5) * t570 - Ifges(5,6) * t569 + Ifges(5,3) * qJDD(4) + t620 * t596 + (mrSges(3,2) - mrSges(4,3)) * t550 + t621 * qJDD(2) + (t593 * t556 + t590 * t557) * qJD(2);
t603 = pkin(7) * t469 + t459 * t594 + t460 * t591;
t458 = mrSges(3,1) * t532 - mrSges(3,2) * t533 + mrSges(4,2) * t529 - mrSges(4,3) * t528 + t593 * t465 - t590 * t470 - pkin(8) * t477 + pkin(2) * t602 + qJ(3) * t598 + (-pkin(2) * mrSges(4,2) + Ifges(4,1) + Ifges(3,3)) * qJDD(2);
t457 = mrSges(2,2) * t580 - mrSges(2,3) * t571 - t591 * t459 + t594 * t460 + (-t463 * t585 - t464 * t588) * pkin(7);
t456 = -mrSges(2,1) * t580 + mrSges(2,3) * t572 - pkin(1) * t463 - t585 * t458 + t603 * t588;
t1 = [-m(1) * g(1) + t609; -m(1) * g(2) + t616; -m(1) * g(3) + m(2) * t580 + t463; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t616 - t584 * t456 + t587 * t457; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t609 + t587 * t456 + t584 * t457; -mrSges(1,1) * g(2) + mrSges(2,1) * t571 + mrSges(1,2) * g(1) - mrSges(2,2) * t572 + pkin(1) * t464 + t588 * t458 + t585 * t603;];
tauB  = t1;
