% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-05-04 23:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRPRRP1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:29:29
% EndTime: 2019-05-04 23:29:38
% DurationCPUTime: 7.29s
% Computational Cost: add. (93838->276), mult. (171776->340), div. (0->0), fcn. (116727->12), ass. (0->121)
t629 = Ifges(6,1) + Ifges(7,1);
t624 = Ifges(6,4) + Ifges(7,4);
t623 = Ifges(6,5) + Ifges(7,5);
t628 = Ifges(6,2) + Ifges(7,2);
t627 = Ifges(6,6) + Ifges(7,6);
t626 = Ifges(6,3) + Ifges(7,3);
t625 = -mrSges(6,2) - mrSges(7,2);
t585 = sin(pkin(6));
t591 = sin(qJ(2));
t621 = t585 * t591;
t594 = cos(qJ(2));
t620 = t585 * t594;
t588 = cos(pkin(6));
t619 = t588 * t591;
t618 = t588 * t594;
t584 = sin(pkin(10));
t587 = cos(pkin(10));
t574 = g(1) * t584 - g(2) * t587;
t575 = -g(1) * t587 - g(2) * t584;
t582 = -g(3) + qJDD(1);
t529 = t574 * t618 - t575 * t591 + t582 * t620;
t527 = qJDD(2) * pkin(2) + t529;
t530 = t574 * t619 + t594 * t575 + t582 * t621;
t596 = qJD(2) ^ 2;
t528 = -pkin(2) * t596 + t530;
t583 = sin(pkin(11));
t586 = cos(pkin(11));
t522 = t583 * t527 + t586 * t528;
t520 = -pkin(3) * t596 + qJDD(2) * pkin(8) + t522;
t555 = -t574 * t585 + t588 * t582;
t554 = qJDD(3) + t555;
t590 = sin(qJ(4));
t593 = cos(qJ(4));
t516 = t593 * t520 + t590 * t554;
t570 = (-mrSges(5,1) * t593 + mrSges(5,2) * t590) * qJD(2);
t610 = qJD(2) * qJD(4);
t606 = t590 * t610;
t573 = qJDD(2) * t593 - t606;
t612 = qJD(2) * t590;
t576 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t612;
t571 = (-pkin(4) * t593 - pkin(9) * t590) * qJD(2);
t595 = qJD(4) ^ 2;
t611 = qJD(2) * t593;
t511 = -pkin(4) * t595 + qJDD(4) * pkin(9) + t571 * t611 + t516;
t521 = t586 * t527 - t583 * t528;
t519 = -qJDD(2) * pkin(3) - t596 * pkin(8) - t521;
t605 = t593 * t610;
t572 = qJDD(2) * t590 + t605;
t514 = (-t572 - t605) * pkin(9) + (-t573 + t606) * pkin(4) + t519;
t589 = sin(qJ(5));
t592 = cos(qJ(5));
t506 = -t589 * t511 + t592 * t514;
t568 = qJD(4) * t592 - t589 * t612;
t543 = qJD(5) * t568 + qJDD(4) * t589 + t572 * t592;
t569 = qJD(4) * t589 + t592 * t612;
t545 = -mrSges(7,1) * t568 + mrSges(7,2) * t569;
t546 = -mrSges(6,1) * t568 + mrSges(6,2) * t569;
t580 = qJD(5) - t611;
t550 = -mrSges(6,2) * t580 + mrSges(6,3) * t568;
t566 = qJDD(5) - t573;
t503 = -0.2e1 * qJD(6) * t569 + (t568 * t580 - t543) * qJ(6) + (t568 * t569 + t566) * pkin(5) + t506;
t549 = -mrSges(7,2) * t580 + mrSges(7,3) * t568;
t609 = m(7) * t503 + t566 * mrSges(7,1) + t580 * t549;
t496 = m(6) * t506 + t566 * mrSges(6,1) + t580 * t550 + (-t545 - t546) * t569 + (-mrSges(6,3) - mrSges(7,3)) * t543 + t609;
t507 = t592 * t511 + t589 * t514;
t542 = -qJD(5) * t569 + qJDD(4) * t592 - t572 * t589;
t551 = pkin(5) * t580 - qJ(6) * t569;
t565 = t568 ^ 2;
t505 = -pkin(5) * t565 + qJ(6) * t542 + 0.2e1 * qJD(6) * t568 - t551 * t580 + t507;
t608 = m(7) * t505 + t542 * mrSges(7,3) + t568 * t545;
t552 = mrSges(7,1) * t580 - mrSges(7,3) * t569;
t613 = -mrSges(6,1) * t580 + mrSges(6,3) * t569 - t552;
t498 = m(6) * t507 + t542 * mrSges(6,3) + t568 * t546 + t625 * t566 + t613 * t580 + t608;
t601 = -t496 * t589 + t592 * t498;
t493 = m(5) * t516 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t573 - qJD(4) * t576 + t570 * t611 + t601;
t515 = -t590 * t520 + t554 * t593;
t577 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t611;
t510 = -qJDD(4) * pkin(4) - pkin(9) * t595 + t571 * t612 - t515;
t508 = -pkin(5) * t542 - qJ(6) * t565 + t551 * t569 + qJDD(6) + t510;
t600 = m(7) * t508 - t542 * mrSges(7,1) - t568 * t549;
t597 = -m(6) * t510 + t542 * mrSges(6,1) + t625 * t543 + t568 * t550 + t613 * t569 - t600;
t500 = m(5) * t515 + qJDD(4) * mrSges(5,1) - t572 * mrSges(5,3) + qJD(4) * t577 - t570 * t612 + t597;
t602 = t593 * t493 - t500 * t590;
t484 = m(4) * t522 - mrSges(4,1) * t596 - qJDD(2) * mrSges(4,2) + t602;
t495 = t496 * t592 + t498 * t589;
t598 = -m(5) * t519 + t573 * mrSges(5,1) - mrSges(5,2) * t572 - t576 * t612 + t577 * t611 - t495;
t490 = m(4) * t521 + qJDD(2) * mrSges(4,1) - mrSges(4,2) * t596 + t598;
t480 = t583 * t484 + t586 * t490;
t478 = m(3) * t529 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t596 + t480;
t603 = t586 * t484 - t490 * t583;
t479 = m(3) * t530 - mrSges(3,1) * t596 - qJDD(2) * mrSges(3,2) + t603;
t487 = t590 * t493 + t593 * t500;
t607 = m(4) * t554 + t487;
t486 = m(3) * t555 + t607;
t466 = t478 * t618 + t479 * t619 - t486 * t585;
t464 = m(2) * t574 + t466;
t470 = -t478 * t591 + t594 * t479;
t469 = m(2) * t575 + t470;
t617 = t587 * t464 + t584 * t469;
t616 = t568 * t627 + t569 * t623 + t580 * t626;
t615 = -t568 * t628 - t569 * t624 - t580 * t627;
t614 = t624 * t568 + t569 * t629 + t623 * t580;
t465 = t478 * t620 + t479 * t621 + t588 * t486;
t604 = -t464 * t584 + t587 * t469;
t488 = -mrSges(6,1) * t510 + mrSges(6,3) * t507 - mrSges(7,1) * t508 + mrSges(7,3) * t505 - pkin(5) * t600 + qJ(6) * t608 + (-qJ(6) * t552 + t614) * t580 + (-pkin(5) * t552 - t616) * t569 + (-mrSges(7,2) * qJ(6) + t627) * t566 + (-mrSges(7,2) * pkin(5) + t624) * t543 + t628 * t542;
t501 = -t543 * mrSges(7,3) - t569 * t545 + t609;
t494 = mrSges(6,2) * t510 + mrSges(7,2) * t508 - mrSges(6,3) * t506 - mrSges(7,3) * t503 - qJ(6) * t501 + t624 * t542 + t543 * t629 + t623 * t566 + t616 * t568 + t615 * t580;
t559 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t590 + Ifges(5,6) * t593) * qJD(2);
t560 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t590 + Ifges(5,2) * t593) * qJD(2);
t472 = mrSges(5,2) * t519 - mrSges(5,3) * t515 + Ifges(5,1) * t572 + Ifges(5,4) * t573 + Ifges(5,5) * qJDD(4) - pkin(9) * t495 - qJD(4) * t560 - t488 * t589 + t494 * t592 + t559 * t611;
t561 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t590 + Ifges(5,4) * t593) * qJD(2);
t481 = -t559 * t612 - mrSges(5,1) * t519 - mrSges(6,1) * t506 - mrSges(7,1) * t503 + mrSges(6,2) * t507 + mrSges(7,2) * t505 + mrSges(5,3) * t516 + Ifges(5,4) * t572 + Ifges(5,2) * t573 + Ifges(5,6) * qJDD(4) - pkin(4) * t495 - pkin(5) * t501 + qJD(4) * t561 + t615 * t569 + t614 * t568 - t626 * t566 - t623 * t543 - t627 * t542;
t462 = mrSges(4,2) * t554 - mrSges(4,3) * t521 + Ifges(4,5) * qJDD(2) - Ifges(4,6) * t596 - pkin(8) * t487 + t472 * t593 - t481 * t590;
t471 = Ifges(4,6) * qJDD(2) + t596 * Ifges(4,5) - mrSges(4,1) * t554 + mrSges(4,3) * t522 - Ifges(5,5) * t572 - Ifges(5,6) * t573 - Ifges(5,3) * qJDD(4) - mrSges(5,1) * t515 + mrSges(5,2) * t516 - t589 * t494 - t592 * t488 - pkin(4) * t597 - pkin(9) * t601 - pkin(3) * t487 + (-t560 * t590 + t561 * t593) * qJD(2);
t459 = -mrSges(3,1) * t555 + mrSges(3,3) * t530 + t596 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t607 + qJ(3) * t603 + t583 * t462 + t586 * t471;
t460 = mrSges(3,2) * t555 - mrSges(3,3) * t529 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t596 - qJ(3) * t480 + t462 * t586 - t471 * t583;
t599 = pkin(7) * t470 + t459 * t594 + t460 * t591;
t461 = mrSges(3,1) * t529 - mrSges(3,2) * t530 + mrSges(4,1) * t521 - mrSges(4,2) * t522 + t590 * t472 + t593 * t481 + pkin(3) * t598 + pkin(8) * t602 + pkin(2) * t480 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2);
t458 = mrSges(2,2) * t582 - mrSges(2,3) * t574 - t591 * t459 + t594 * t460 + (-t465 * t585 - t466 * t588) * pkin(7);
t457 = -mrSges(2,1) * t582 + mrSges(2,3) * t575 - pkin(1) * t465 - t585 * t461 + t599 * t588;
t1 = [-m(1) * g(1) + t604; -m(1) * g(2) + t617; -m(1) * g(3) + m(2) * t582 + t465; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t617 - t584 * t457 + t587 * t458; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t604 + t587 * t457 + t584 * t458; -mrSges(1,1) * g(2) + mrSges(2,1) * t574 + mrSges(1,2) * g(1) - mrSges(2,2) * t575 + pkin(1) * t466 + t588 * t461 + t599 * t585;];
tauB  = t1;
