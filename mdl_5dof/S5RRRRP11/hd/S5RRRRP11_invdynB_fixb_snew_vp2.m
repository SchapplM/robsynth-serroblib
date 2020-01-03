% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRRP11
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRRP11_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP11_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP11_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP11_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP11_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP11_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP11_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:15:01
% EndTime: 2019-12-31 22:15:11
% DurationCPUTime: 6.40s
% Computational Cost: add. (92283->301), mult. (197250->384), div. (0->0), fcn. (149782->10), ass. (0->127)
t627 = Ifges(5,1) + Ifges(6,1);
t619 = Ifges(5,4) - Ifges(6,5);
t626 = -Ifges(5,5) - Ifges(6,4);
t625 = Ifges(5,2) + Ifges(6,3);
t617 = Ifges(5,6) - Ifges(6,6);
t624 = -Ifges(5,3) - Ifges(6,2);
t583 = sin(pkin(5));
t587 = sin(qJ(2));
t590 = cos(qJ(2));
t604 = qJD(1) * qJD(2);
t572 = (-qJDD(1) * t590 + t587 * t604) * t583;
t623 = cos(qJ(4));
t622 = pkin(7) * t583;
t584 = cos(pkin(5));
t621 = t584 * g(3);
t620 = -mrSges(5,3) - mrSges(6,2);
t616 = t583 * t587;
t615 = t583 * t590;
t614 = t584 * t587;
t613 = t584 * t590;
t588 = sin(qJ(1));
t591 = cos(qJ(1));
t576 = t588 * g(1) - t591 * g(2);
t592 = qJD(1) ^ 2;
t567 = qJDD(1) * pkin(1) + t592 * t622 + t576;
t577 = -t591 * g(1) - t588 * g(2);
t568 = -t592 * pkin(1) + qJDD(1) * t622 + t577;
t607 = t567 * t614 + t590 * t568;
t544 = -g(3) * t616 + t607;
t580 = t584 * qJD(1) + qJD(2);
t606 = qJD(1) * t583;
t602 = t587 * t606;
t565 = t580 * mrSges(3,1) - mrSges(3,3) * t602;
t569 = (-mrSges(3,1) * t590 + mrSges(3,2) * t587) * t606;
t579 = t584 * qJDD(1) + qJDD(2);
t570 = (-pkin(2) * t590 - pkin(8) * t587) * t606;
t578 = t580 ^ 2;
t605 = qJD(1) * t590;
t521 = -t578 * pkin(2) + t579 * pkin(8) + (-g(3) * t587 + t570 * t605) * t583 + t607;
t571 = (qJDD(1) * t587 + t590 * t604) * t583;
t522 = t572 * pkin(2) - t571 * pkin(8) - t621 + (-t567 + (pkin(2) * t587 - pkin(8) * t590) * t580 * qJD(1)) * t583;
t586 = sin(qJ(3));
t589 = cos(qJ(3));
t504 = t589 * t521 + t586 * t522;
t561 = t586 * t580 + t589 * t602;
t541 = -t561 * qJD(3) - t586 * t571 + t589 * t579;
t560 = t589 * t580 - t586 * t602;
t545 = -t560 * mrSges(4,1) + t561 * mrSges(4,2);
t601 = t583 * t605;
t575 = qJD(3) - t601;
t550 = t575 * mrSges(4,1) - t561 * mrSges(4,3);
t564 = qJDD(3) + t572;
t546 = -t560 * pkin(3) - t561 * pkin(9);
t574 = t575 ^ 2;
t500 = -t574 * pkin(3) + t564 * pkin(9) + t560 * t546 + t504;
t543 = -g(3) * t615 + t567 * t613 - t587 * t568;
t520 = -t579 * pkin(2) - t578 * pkin(8) + t570 * t602 - t543;
t542 = t560 * qJD(3) + t589 * t571 + t586 * t579;
t502 = (-t560 * t575 - t542) * pkin(9) + (t561 * t575 - t541) * pkin(3) + t520;
t585 = sin(qJ(4));
t497 = t623 * t500 + t585 * t502;
t548 = t623 * t561 + t585 * t575;
t507 = t548 * qJD(4) + t585 * t542 - t623 * t564;
t558 = qJD(4) - t560;
t531 = t558 * mrSges(5,1) - t548 * mrSges(5,3);
t539 = qJDD(4) - t541;
t547 = t585 * t561 - t623 * t575;
t525 = t547 * pkin(4) - t548 * qJ(5);
t557 = t558 ^ 2;
t493 = -t557 * pkin(4) + t539 * qJ(5) + 0.2e1 * qJD(5) * t558 - t547 * t525 + t497;
t532 = -t558 * mrSges(6,1) + t548 * mrSges(6,2);
t603 = m(6) * t493 + t539 * mrSges(6,3) + t558 * t532;
t526 = t547 * mrSges(6,1) - t548 * mrSges(6,3);
t608 = -t547 * mrSges(5,1) - t548 * mrSges(5,2) - t526;
t489 = m(5) * t497 - t539 * mrSges(5,2) + t620 * t507 - t558 * t531 + t608 * t547 + t603;
t496 = -t585 * t500 + t623 * t502;
t508 = -t547 * qJD(4) + t623 * t542 + t585 * t564;
t530 = -t558 * mrSges(5,2) - t547 * mrSges(5,3);
t494 = -t539 * pkin(4) - t557 * qJ(5) + t548 * t525 + qJDD(5) - t496;
t529 = -t547 * mrSges(6,2) + t558 * mrSges(6,3);
t597 = -m(6) * t494 + t539 * mrSges(6,1) + t558 * t529;
t490 = m(5) * t496 + t539 * mrSges(5,1) + t620 * t508 + t558 * t530 + t608 * t548 + t597;
t598 = t623 * t489 - t585 * t490;
t483 = m(4) * t504 - t564 * mrSges(4,2) + t541 * mrSges(4,3) + t560 * t545 - t575 * t550 + t598;
t503 = -t586 * t521 + t589 * t522;
t549 = -t575 * mrSges(4,2) + t560 * mrSges(4,3);
t499 = -t564 * pkin(3) - t574 * pkin(9) + t561 * t546 - t503;
t495 = -0.2e1 * qJD(5) * t548 + (t547 * t558 - t508) * qJ(5) + (t548 * t558 + t507) * pkin(4) + t499;
t491 = m(6) * t495 + t507 * mrSges(6,1) - t508 * mrSges(6,3) + t547 * t529 - t548 * t532;
t593 = -m(5) * t499 - t507 * mrSges(5,1) - t508 * mrSges(5,2) - t547 * t530 - t548 * t531 - t491;
t487 = m(4) * t503 + t564 * mrSges(4,1) - t542 * mrSges(4,3) - t561 * t545 + t575 * t549 + t593;
t599 = t589 * t483 - t586 * t487;
t473 = m(3) * t544 - t579 * mrSges(3,2) - t572 * mrSges(3,3) - t580 * t565 + t569 * t601 + t599;
t476 = t586 * t483 + t589 * t487;
t554 = -t583 * t567 - t621;
t566 = -t580 * mrSges(3,2) + mrSges(3,3) * t601;
t475 = m(3) * t554 + t572 * mrSges(3,1) + t571 * mrSges(3,2) + (t565 * t587 - t566 * t590) * t606 + t476;
t485 = t585 * t489 + t623 * t490;
t594 = -m(4) * t520 + t541 * mrSges(4,1) - t542 * mrSges(4,2) + t560 * t549 - t561 * t550 - t485;
t479 = m(3) * t543 + t579 * mrSges(3,1) - t571 * mrSges(3,3) + t580 * t566 - t569 * t602 + t594;
t463 = t473 * t614 - t583 * t475 + t479 * t613;
t461 = m(2) * t576 + qJDD(1) * mrSges(2,1) - t592 * mrSges(2,2) + t463;
t468 = t590 * t473 - t587 * t479;
t467 = m(2) * t577 - t592 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t468;
t612 = t591 * t461 + t588 * t467;
t611 = t625 * t547 - t619 * t548 - t617 * t558;
t610 = t617 * t547 + t626 * t548 + t624 * t558;
t609 = -t619 * t547 + t627 * t548 - t626 * t558;
t462 = t473 * t616 + t584 * t475 + t479 * t615;
t600 = -t588 * t461 + t591 * t467;
t480 = -mrSges(5,1) * t499 - mrSges(6,1) * t495 + mrSges(6,2) * t493 + mrSges(5,3) * t497 - pkin(4) * t491 - t625 * t507 + t619 * t508 + t617 * t539 + t610 * t548 + t609 * t558;
t484 = mrSges(5,2) * t499 + mrSges(6,2) * t494 - mrSges(5,3) * t496 - mrSges(6,3) * t495 - qJ(5) * t491 - t619 * t507 + t627 * t508 - t539 * t626 + t610 * t547 + t611 * t558;
t535 = Ifges(4,5) * t561 + Ifges(4,6) * t560 + Ifges(4,3) * t575;
t536 = Ifges(4,4) * t561 + Ifges(4,2) * t560 + Ifges(4,6) * t575;
t464 = mrSges(4,2) * t520 - mrSges(4,3) * t503 + Ifges(4,1) * t542 + Ifges(4,4) * t541 + Ifges(4,5) * t564 - pkin(9) * t485 - t585 * t480 + t623 * t484 + t560 * t535 - t575 * t536;
t537 = Ifges(4,1) * t561 + Ifges(4,4) * t560 + Ifges(4,5) * t575;
t469 = Ifges(4,4) * t542 + Ifges(4,2) * t541 + Ifges(4,6) * t564 - t561 * t535 + t575 * t537 - mrSges(4,1) * t520 + mrSges(4,3) * t504 - mrSges(5,1) * t496 + mrSges(5,2) * t497 + mrSges(6,1) * t494 - mrSges(6,3) * t493 - pkin(4) * t597 - qJ(5) * t603 - pkin(3) * t485 + (pkin(4) * t526 + t611) * t548 + (qJ(5) * t526 - t609) * t547 + t624 * t539 + (pkin(4) * mrSges(6,2) + t626) * t508 + (qJ(5) * mrSges(6,2) + t617) * t507;
t551 = Ifges(3,3) * t580 + (Ifges(3,5) * t587 + Ifges(3,6) * t590) * t606;
t552 = Ifges(3,6) * t580 + (Ifges(3,4) * t587 + Ifges(3,2) * t590) * t606;
t458 = mrSges(3,2) * t554 - mrSges(3,3) * t543 + Ifges(3,1) * t571 - Ifges(3,4) * t572 + Ifges(3,5) * t579 - pkin(8) * t476 + t589 * t464 - t586 * t469 + t551 * t601 - t580 * t552;
t553 = Ifges(3,5) * t580 + (Ifges(3,1) * t587 + Ifges(3,4) * t590) * t606;
t459 = Ifges(3,4) * t571 - Ifges(3,2) * t572 + Ifges(3,6) * t579 - t551 * t602 + t580 * t553 - mrSges(3,1) * t554 + mrSges(3,3) * t544 - Ifges(4,5) * t542 - Ifges(4,6) * t541 - Ifges(4,3) * t564 - t561 * t536 + t560 * t537 - mrSges(4,1) * t503 + mrSges(4,2) * t504 - t585 * t484 - t623 * t480 - pkin(3) * t593 - pkin(9) * t598 - pkin(2) * t476;
t595 = pkin(7) * t468 + t458 * t587 + t459 * t590;
t457 = Ifges(3,5) * t571 - Ifges(3,6) * t572 + Ifges(3,3) * t579 + mrSges(3,1) * t543 - mrSges(3,2) * t544 + t586 * t464 + t589 * t469 + pkin(2) * t594 + pkin(8) * t599 + (t552 * t587 - t553 * t590) * t606;
t456 = -mrSges(2,2) * g(3) - mrSges(2,3) * t576 + Ifges(2,5) * qJDD(1) - t592 * Ifges(2,6) + t590 * t458 - t587 * t459 + (-t462 * t583 - t463 * t584) * pkin(7);
t455 = mrSges(2,1) * g(3) + mrSges(2,3) * t577 + t592 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t462 - t583 * t457 + t595 * t584;
t1 = [-m(1) * g(1) + t600; -m(1) * g(2) + t612; (-m(1) - m(2)) * g(3) + t462; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t612 - t588 * t455 + t591 * t456; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t600 + t591 * t455 + t588 * t456; -mrSges(1,1) * g(2) + mrSges(2,1) * t576 + mrSges(1,2) * g(1) - mrSges(2,2) * t577 + Ifges(2,3) * qJDD(1) + pkin(1) * t463 + t584 * t457 + t595 * t583;];
tauB = t1;
