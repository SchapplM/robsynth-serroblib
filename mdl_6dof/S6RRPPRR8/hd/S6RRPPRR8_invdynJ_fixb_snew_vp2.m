% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPPRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 11:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPPRR8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR8_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR8_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR8_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR8_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR8_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 11:18:19
% EndTime: 2019-05-06 11:18:27
% DurationCPUTime: 4.72s
% Computational Cost: add. (42242->329), mult. (93148->401), div. (0->0), fcn. (63201->10), ass. (0->131)
t597 = -2 * qJD(3);
t596 = -2 * qJD(4);
t595 = Ifges(4,1) + Ifges(5,1);
t589 = Ifges(4,4) - Ifges(5,5);
t594 = Ifges(4,5) + Ifges(5,4);
t593 = Ifges(4,2) + Ifges(5,3);
t592 = Ifges(5,2) + Ifges(4,3);
t591 = Ifges(4,6) - Ifges(5,6);
t561 = qJD(1) ^ 2;
t555 = sin(qJ(1));
t559 = cos(qJ(1));
t570 = -g(1) * t559 - g(2) * t555;
t523 = -pkin(1) * t561 + qJDD(1) * pkin(7) + t570;
t554 = sin(qJ(2));
t558 = cos(qJ(2));
t498 = -t558 * g(3) - t554 * t523;
t533 = (-pkin(2) * t558 - qJ(3) * t554) * qJD(1);
t560 = qJD(2) ^ 2;
t580 = qJD(1) * t554;
t476 = -qJDD(2) * pkin(2) - t560 * qJ(3) + t533 * t580 + qJDD(3) - t498;
t577 = qJD(1) * qJD(2);
t575 = t558 * t577;
t535 = qJDD(1) * t554 + t575;
t551 = sin(pkin(10));
t586 = cos(pkin(10));
t508 = -qJDD(2) * t586 + t535 * t551;
t509 = t551 * qJDD(2) + t535 * t586;
t529 = t551 * qJD(2) + t580 * t586;
t545 = t558 * qJD(1);
t528 = -qJD(2) * t586 + t551 * t580;
t576 = t528 * t545;
t452 = t476 + (-t509 - t576) * qJ(4) + (-t529 * t545 + t508) * pkin(3) + t529 * t596;
t574 = g(1) * t555 - t559 * g(2);
t522 = -qJDD(1) * pkin(1) - pkin(7) * t561 - t574;
t542 = t554 * t577;
t536 = qJDD(1) * t558 - t542;
t473 = (-t535 - t575) * qJ(3) + (-t536 + t542) * pkin(2) + t522;
t499 = -g(3) * t554 + t558 * t523;
t477 = -pkin(2) * t560 + qJDD(2) * qJ(3) + t533 * t545 + t499;
t453 = t473 * t586 - t551 * t477 + t529 * t597;
t590 = -mrSges(4,3) - mrSges(5,2);
t585 = t558 ^ 2 * t561;
t495 = pkin(3) * t528 - qJ(4) * t529;
t445 = t536 * pkin(3) - qJ(4) * t585 + t529 * t495 + qJDD(4) - t453;
t431 = (-t509 + t576) * pkin(8) + (t528 * t529 + t536) * pkin(4) + t445;
t454 = t551 * t473 + t586 * t477 + t528 * t597;
t444 = -pkin(3) * t585 - t536 * qJ(4) - t528 * t495 + t545 * t596 + t454;
t510 = pkin(4) * t545 - pkin(8) * t529;
t526 = t528 ^ 2;
t433 = -pkin(4) * t526 + pkin(8) * t508 - t510 * t545 + t444;
t553 = sin(qJ(5));
t557 = cos(qJ(5));
t425 = t557 * t431 - t433 * t553;
t493 = t528 * t557 - t529 * t553;
t463 = qJD(5) * t493 + t508 * t553 + t509 * t557;
t494 = t528 * t553 + t529 * t557;
t532 = qJDD(5) + t536;
t540 = t545 + qJD(5);
t422 = (t493 * t540 - t463) * pkin(9) + (t493 * t494 + t532) * pkin(5) + t425;
t426 = t553 * t431 + t557 * t433;
t462 = -qJD(5) * t494 + t508 * t557 - t509 * t553;
t482 = pkin(5) * t540 - pkin(9) * t494;
t492 = t493 ^ 2;
t423 = -pkin(5) * t492 + pkin(9) * t462 - t482 * t540 + t426;
t552 = sin(qJ(6));
t556 = cos(qJ(6));
t420 = t422 * t556 - t423 * t552;
t469 = t493 * t556 - t494 * t552;
t439 = qJD(6) * t469 + t462 * t552 + t463 * t556;
t470 = t493 * t552 + t494 * t556;
t451 = -mrSges(7,1) * t469 + mrSges(7,2) * t470;
t539 = qJD(6) + t540;
t456 = -mrSges(7,2) * t539 + mrSges(7,3) * t469;
t525 = qJDD(6) + t532;
t417 = m(7) * t420 + mrSges(7,1) * t525 - mrSges(7,3) * t439 - t451 * t470 + t456 * t539;
t421 = t422 * t552 + t423 * t556;
t438 = -qJD(6) * t470 + t462 * t556 - t463 * t552;
t457 = mrSges(7,1) * t539 - mrSges(7,3) * t470;
t418 = m(7) * t421 - mrSges(7,2) * t525 + mrSges(7,3) * t438 + t451 * t469 - t457 * t539;
t409 = t556 * t417 + t552 * t418;
t584 = t593 * t528 - t589 * t529 + t591 * t545;
t583 = t591 * t528 - t594 * t529 + t592 * t545;
t582 = t589 * t528 - t595 * t529 + t594 * t545;
t496 = mrSges(5,1) * t528 - mrSges(5,3) * t529;
t581 = -mrSges(4,1) * t528 - mrSges(4,2) * t529 - t496;
t506 = -mrSges(4,1) * t545 - mrSges(4,3) * t529;
t507 = mrSges(5,1) * t545 + mrSges(5,2) * t529;
t471 = -mrSges(6,1) * t493 + mrSges(6,2) * t494;
t478 = -mrSges(6,2) * t540 + mrSges(6,3) * t493;
t406 = m(6) * t425 + mrSges(6,1) * t532 - mrSges(6,3) * t463 - t471 * t494 + t478 * t540 + t409;
t479 = mrSges(6,1) * t540 - mrSges(6,3) * t494;
t571 = -t417 * t552 + t556 * t418;
t407 = m(6) * t426 - mrSges(6,2) * t532 + mrSges(6,3) * t462 + t471 * t493 - t479 * t540 + t571;
t572 = -t406 * t553 + t557 * t407;
t567 = m(5) * t444 - t536 * mrSges(5,3) + t572;
t401 = m(4) * t454 + mrSges(4,2) * t536 + t581 * t528 + t590 * t508 + (t506 - t507) * t545 + t567;
t504 = -mrSges(5,2) * t528 - mrSges(5,3) * t545;
t505 = mrSges(4,2) * t545 - mrSges(4,3) * t528;
t404 = t406 * t557 + t407 * t553;
t566 = -m(5) * t445 - t536 * mrSges(5,1) - t404;
t402 = m(4) * t453 - mrSges(4,1) * t536 + t581 * t529 + t590 * t509 + (-t504 - t505) * t545 + t566;
t573 = t586 * t401 - t402 * t551;
t442 = -pkin(4) * t508 - pkin(8) * t526 + t529 * t510 - t452;
t428 = -pkin(5) * t462 - pkin(9) * t492 + t482 * t494 + t442;
t568 = m(7) * t428 - t438 * mrSges(7,1) + t439 * mrSges(7,2) - t469 * t456 + t470 * t457;
t397 = t551 * t401 + t402 * t586;
t565 = -m(6) * t442 + t462 * mrSges(6,1) - t463 * mrSges(6,2) + t493 * t478 - t494 * t479 - t568;
t447 = Ifges(7,4) * t470 + Ifges(7,2) * t469 + Ifges(7,6) * t539;
t448 = Ifges(7,1) * t470 + Ifges(7,4) * t469 + Ifges(7,5) * t539;
t564 = mrSges(7,1) * t420 - mrSges(7,2) * t421 + Ifges(7,5) * t439 + Ifges(7,6) * t438 + Ifges(7,3) * t525 + t470 * t447 - t469 * t448;
t413 = m(5) * t452 + t508 * mrSges(5,1) - t509 * mrSges(5,3) + t528 * t504 - t529 * t507 + t565;
t465 = Ifges(6,4) * t494 + Ifges(6,2) * t493 + Ifges(6,6) * t540;
t466 = Ifges(6,1) * t494 + Ifges(6,4) * t493 + Ifges(6,5) * t540;
t562 = mrSges(6,1) * t425 - mrSges(6,2) * t426 + Ifges(6,5) * t463 + Ifges(6,6) * t462 + Ifges(6,3) * t532 + pkin(5) * t409 + t494 * t465 - t493 * t466 + t564;
t412 = m(4) * t476 + t508 * mrSges(4,1) + t509 * mrSges(4,2) + t528 * t505 + t529 * t506 + t413;
t538 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t545;
t537 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t580;
t534 = (-mrSges(3,1) * t558 + mrSges(3,2) * t554) * qJD(1);
t521 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t554 + Ifges(3,4) * t558) * qJD(1);
t520 = Ifges(3,6) * qJD(2) + (t554 * Ifges(3,4) + Ifges(3,2) * t558) * qJD(1);
t519 = Ifges(3,3) * qJD(2) + (t554 * Ifges(3,5) + t558 * Ifges(3,6)) * qJD(1);
t464 = Ifges(6,5) * t494 + Ifges(6,6) * t493 + Ifges(6,3) * t540;
t446 = Ifges(7,5) * t470 + Ifges(7,6) * t469 + Ifges(7,3) * t539;
t411 = mrSges(7,2) * t428 - mrSges(7,3) * t420 + Ifges(7,1) * t439 + Ifges(7,4) * t438 + Ifges(7,5) * t525 + t446 * t469 - t447 * t539;
t410 = -mrSges(7,1) * t428 + mrSges(7,3) * t421 + Ifges(7,4) * t439 + Ifges(7,2) * t438 + Ifges(7,6) * t525 - t446 * t470 + t448 * t539;
t403 = mrSges(5,2) * t509 + t496 * t529 + t504 * t545 - t566;
t399 = mrSges(6,2) * t442 - mrSges(6,3) * t425 + Ifges(6,1) * t463 + Ifges(6,4) * t462 + Ifges(6,5) * t532 - pkin(9) * t409 - t410 * t552 + t411 * t556 + t464 * t493 - t465 * t540;
t398 = -mrSges(6,1) * t442 + mrSges(6,3) * t426 + Ifges(6,4) * t463 + Ifges(6,2) * t462 + Ifges(6,6) * t532 - pkin(5) * t568 + pkin(9) * t571 + t556 * t410 + t552 * t411 - t494 * t464 + t540 * t466;
t396 = mrSges(4,2) * t476 + mrSges(5,2) * t445 - mrSges(4,3) * t453 - mrSges(5,3) * t452 - pkin(8) * t404 - qJ(4) * t413 - t553 * t398 + t557 * t399 - t589 * t508 + t595 * t509 + t583 * t528 - t536 * t594 - t584 * t545;
t395 = -mrSges(4,1) * t476 - mrSges(5,1) * t452 + mrSges(5,2) * t444 + mrSges(4,3) * t454 - pkin(3) * t413 - pkin(4) * t565 - pkin(8) * t572 - t557 * t398 - t553 * t399 - t593 * t508 + t589 * t509 + t583 * t529 - t536 * t591 + t582 * t545;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t574 - mrSges(2,2) * t570 + t554 * (mrSges(3,2) * t522 - mrSges(3,3) * t498 + Ifges(3,1) * t535 + Ifges(3,4) * t536 + Ifges(3,5) * qJDD(2) - qJ(3) * t397 - qJD(2) * t520 - t551 * t395 + t586 * t396 + t519 * t545) + t558 * (Ifges(3,4) * t535 + qJD(2) * t521 - mrSges(3,1) * t522 + mrSges(3,3) * t499 + Ifges(3,6) * qJDD(2) - mrSges(5,3) * t444 + mrSges(5,1) * t445 - mrSges(4,1) * t453 + mrSges(4,2) * t454 - qJ(4) * (-t507 * t545 + t567) + pkin(4) * t404 + pkin(3) * t403 + t562 - pkin(2) * t397 + (Ifges(3,2) + t592) * t536 + t584 * t529 + (qJ(4) * t496 + t582) * t528 - t594 * t509 + (qJ(4) * mrSges(5,2) + t591) * t508 - t519 * t580) + pkin(1) * (-m(3) * t522 + t536 * mrSges(3,1) - t535 * mrSges(3,2) + (-t537 * t554 + t538 * t558) * qJD(1) - t397) + pkin(7) * (t558 * (m(3) * t499 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t536 - qJD(2) * t537 + t534 * t545 + t573) - t554 * (m(3) * t498 + qJDD(2) * mrSges(3,1) - t535 * mrSges(3,3) + qJD(2) * t538 - t534 * t580 - t412)); Ifges(3,5) * t535 + Ifges(3,6) * t536 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t498 - mrSges(3,2) * t499 + t551 * t396 + t586 * t395 - pkin(2) * t412 + qJ(3) * t573 + (t554 * t520 - t558 * t521) * qJD(1); t412; t403; t562; t564;];
tauJ  = t1;
