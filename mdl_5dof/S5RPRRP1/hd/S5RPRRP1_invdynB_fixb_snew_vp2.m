% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRRP1
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRRP1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP1_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP1_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP1_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP1_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:28
% EndTime: 2019-12-05 17:59:32
% DurationCPUTime: 1.66s
% Computational Cost: add. (12573->245), mult. (24983->285), div. (0->0), fcn. (14356->6), ass. (0->94)
t545 = Ifges(5,1) + Ifges(6,1);
t537 = Ifges(5,4) + Ifges(6,4);
t535 = Ifges(5,5) + Ifges(6,5);
t544 = Ifges(5,2) + Ifges(6,2);
t533 = Ifges(5,6) + Ifges(6,6);
t543 = Ifges(5,3) + Ifges(6,3);
t508 = qJD(1) ^ 2;
t504 = sin(qJ(1));
t507 = cos(qJ(1));
t490 = -g(1) * t507 - g(2) * t504;
t513 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t490;
t540 = -pkin(1) - pkin(6);
t469 = t540 * t508 + t513;
t503 = sin(qJ(3));
t506 = cos(qJ(3));
t524 = qJD(1) * qJD(3);
t484 = -qJDD(1) * t503 - t506 * t524;
t519 = t503 * t524;
t485 = qJDD(1) * t506 - t519;
t526 = qJD(1) * t503;
t486 = -(qJD(3) * mrSges(4,2)) - mrSges(4,3) * t526;
t525 = qJD(1) * t506;
t487 = (qJD(3) * mrSges(4,1)) - mrSges(4,3) * t525;
t488 = (qJD(3) * pkin(3)) - pkin(7) * t525;
t499 = t503 ^ 2;
t441 = -pkin(3) * t484 + t488 * t525 + (-pkin(7) * t499 + t540) * t508 + t513;
t502 = sin(qJ(4));
t505 = cos(qJ(4));
t478 = (-t502 * t503 + t505 * t506) * qJD(1);
t446 = -qJD(4) * t478 + t484 * t505 - t485 * t502;
t477 = (-t502 * t506 - t503 * t505) * qJD(1);
t447 = qJD(4) * t477 + t484 * t502 + t485 * t505;
t496 = qJD(3) + qJD(4);
t464 = -mrSges(6,2) * t496 + mrSges(6,3) * t477;
t465 = -mrSges(5,2) * t496 + mrSges(5,3) * t477;
t468 = mrSges(5,1) * t496 - mrSges(5,3) * t478;
t466 = pkin(4) * t496 - qJ(5) * t478;
t473 = t477 ^ 2;
t432 = -pkin(4) * t446 - qJ(5) * t473 + t466 * t478 + qJDD(5) + t441;
t467 = mrSges(6,1) * t496 - mrSges(6,3) * t478;
t520 = m(6) * t432 + t447 * mrSges(6,2) + t478 * t467;
t510 = m(5) * t441 + t447 * mrSges(5,2) - (mrSges(5,1) + mrSges(6,1)) * t446 + t478 * t468 - (t464 + t465) * t477 + t520;
t542 = -m(4) * t469 + mrSges(4,1) * t484 - t485 * mrSges(4,2) - t486 * t526 - t487 * t525 - t510;
t539 = mrSges(2,1) - mrSges(3,2);
t536 = Ifges(2,5) - Ifges(3,4);
t534 = -Ifges(2,6) + Ifges(3,5);
t489 = g(1) * t504 - t507 * g(2);
t512 = -qJ(2) * t508 + qJDD(2) - t489;
t470 = t540 * qJDD(1) + t512;
t458 = t503 * g(3) + t506 * t470;
t438 = (-t485 - t519) * pkin(7) + (-t503 * t506 * t508 + qJDD(3)) * pkin(3) + t458;
t459 = -g(3) * t506 + t503 * t470;
t439 = -pkin(3) * t499 * t508 + pkin(7) * t484 - qJD(3) * t488 + t459;
t433 = t505 * t438 - t439 * t502;
t456 = -mrSges(6,1) * t477 + mrSges(6,2) * t478;
t457 = -mrSges(5,1) * t477 + mrSges(5,2) * t478;
t495 = qJDD(3) + qJDD(4);
t428 = -0.2e1 * qJD(5) * t478 + (t477 * t496 - t447) * qJ(5) + (t477 * t478 + t495) * pkin(4) + t433;
t522 = m(6) * t428 + t495 * mrSges(6,1) + t496 * t464;
t422 = m(5) * t433 + mrSges(5,1) * t495 + t465 * t496 + (-t456 - t457) * t478 + (-mrSges(5,3) - mrSges(6,3)) * t447 + t522;
t434 = t502 * t438 + t505 * t439;
t430 = -pkin(4) * t473 + qJ(5) * t446 + 0.2e1 * qJD(5) * t477 - t466 * t496 + t434;
t521 = m(6) * t430 + t446 * mrSges(6,3) + t477 * t456;
t425 = m(5) * t434 + mrSges(5,3) * t446 + t457 * t477 + (-t467 - t468) * t496 + (-mrSges(5,2) - mrSges(6,2)) * t495 + t521;
t417 = t505 * t422 + t502 * t425;
t483 = (mrSges(4,1) * t503 + mrSges(4,2) * t506) * qJD(1);
t415 = m(4) * t458 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t485 + qJD(3) * t486 - t483 * t525 + t417;
t516 = -t422 * t502 + t505 * t425;
t416 = m(4) * t459 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t484 - qJD(3) * t487 - t483 * t526 + t516;
t411 = t415 * t506 + t416 * t503;
t472 = -qJDD(1) * pkin(1) + t512;
t511 = -m(3) * t472 + t508 * mrSges(3,3) - t411;
t409 = m(2) * t489 - mrSges(2,2) * t508 + t539 * qJDD(1) + t511;
t471 = pkin(1) * t508 - t513;
t509 = -m(3) * t471 + t508 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t542;
t420 = m(2) * t490 - mrSges(2,1) * t508 - qJDD(1) * mrSges(2,2) + t509;
t531 = t507 * t409 + t504 * t420;
t530 = -t533 * t477 - t535 * t478 - t543 * t496;
t529 = t544 * t477 + t537 * t478 + t533 * t496;
t528 = -t537 * t477 - t545 * t478 - t535 * t496;
t518 = -t409 * t504 + t507 * t420;
t517 = -t503 * t415 + t506 * t416;
t476 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t506 - Ifges(4,4) * t503) * qJD(1);
t475 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t506 - Ifges(4,2) * t503) * qJD(1);
t474 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t506 - Ifges(4,6) * t503) * qJD(1);
t426 = -mrSges(6,3) * t447 - t456 * t478 + t522;
t413 = mrSges(5,2) * t441 + mrSges(6,2) * t432 - mrSges(5,3) * t433 - mrSges(6,3) * t428 - qJ(5) * t426 + t537 * t446 + t447 * t545 - t530 * t477 + t535 * t495 - t529 * t496;
t412 = -mrSges(5,1) * t441 + mrSges(5,3) * t434 - mrSges(6,1) * t432 + mrSges(6,3) * t430 - pkin(4) * (-t464 * t477 + t520) + qJ(5) * t521 + (-qJ(5) * t467 - t528) * t496 + (-mrSges(6,2) * qJ(5) + t533) * t495 + t530 * t478 + t537 * t447 + (mrSges(6,1) * pkin(4) + t544) * t446;
t410 = -m(3) * g(3) + t517;
t407 = mrSges(4,2) * t469 - mrSges(4,3) * t458 + Ifges(4,1) * t485 + Ifges(4,4) * t484 + Ifges(4,5) * qJDD(3) - pkin(7) * t417 - qJD(3) * t475 - t412 * t502 + t413 * t505 - t474 * t526;
t406 = -mrSges(4,1) * t469 + mrSges(4,3) * t459 + Ifges(4,4) * t485 + Ifges(4,2) * t484 + Ifges(4,6) * qJDD(3) - pkin(3) * t510 + pkin(7) * t516 + qJD(3) * t476 + t505 * t412 + t502 * t413 - t474 * t525;
t405 = mrSges(3,1) * t472 + mrSges(4,1) * t458 + mrSges(5,1) * t433 + mrSges(6,1) * t428 - mrSges(4,2) * t459 - mrSges(5,2) * t434 - mrSges(6,2) * t430 - mrSges(2,3) * t489 + Ifges(4,5) * t485 + Ifges(4,6) * t484 + Ifges(4,3) * qJDD(3) + pkin(2) * t411 + pkin(3) * t417 + pkin(4) * t426 - qJ(2) * t410 + t534 * t508 + t543 * t495 + t529 * t478 + t528 * t477 + t535 * t447 + t533 * t446 + t536 * qJDD(1) + (t475 * t506 + t476 * t503) * qJD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t404 = -mrSges(3,1) * t471 + mrSges(2,3) * t490 - pkin(1) * t410 - pkin(2) * t542 - pkin(6) * t517 + t539 * g(3) - t534 * qJDD(1) - t506 * t406 - t503 * t407 + t536 * t508;
t1 = [-m(1) * g(1) + t518; -m(1) * g(2) + t531; (-m(1) - m(2) - m(3)) * g(3) + t517; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t531 - t504 * t404 + t507 * t405; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t518 + t507 * t404 + t504 * t405; pkin(1) * t511 + qJ(2) * t509 + t506 * t407 - t503 * t406 - pkin(6) * t411 + mrSges(2,1) * t489 - mrSges(2,2) * t490 + mrSges(3,2) * t472 - mrSges(3,3) * t471 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-mrSges(3,2) * pkin(1) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB = t1;
