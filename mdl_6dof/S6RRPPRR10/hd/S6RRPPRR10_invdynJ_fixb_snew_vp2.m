% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-05-06 11:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPPRR10_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR10_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR10_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR10_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR10_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR10_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR10_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 11:47:16
% EndTime: 2019-05-06 11:47:24
% DurationCPUTime: 5.13s
% Computational Cost: add. (49497->334), mult. (109011->408), div. (0->0), fcn. (69490->10), ass. (0->133)
t573 = -2 * qJD(3);
t572 = Ifges(3,1) + Ifges(4,2);
t567 = Ifges(3,4) + Ifges(4,6);
t566 = Ifges(3,5) - Ifges(4,4);
t571 = Ifges(3,2) + Ifges(4,3);
t565 = Ifges(3,6) - Ifges(4,5);
t570 = (Ifges(3,3) + Ifges(4,1));
t540 = qJD(1) ^ 2;
t534 = sin(qJ(1));
t538 = cos(qJ(1));
t551 = -g(1) * t538 - g(2) * t534;
t497 = -pkin(1) * t540 + qJDD(1) * pkin(7) + t551;
t533 = sin(qJ(2));
t537 = cos(qJ(2));
t475 = -g(3) * t533 + t537 * t497;
t505 = (-pkin(2) * t537 - qJ(3) * t533) * qJD(1);
t539 = qJD(2) ^ 2;
t560 = qJD(1) * t537;
t460 = pkin(2) * t539 - qJDD(2) * qJ(3) + (qJD(2) * t573) - t505 * t560 - t475;
t569 = t540 * pkin(7);
t568 = mrSges(3,1) - mrSges(4,2);
t559 = qJD(1) * qJD(2);
t556 = t533 * t559;
t509 = qJDD(1) * t537 - t556;
t522 = t533 * qJD(1);
t514 = pkin(3) * t522 - (qJD(2) * qJ(4));
t528 = t537 ^ 2;
t557 = t537 * t559;
t508 = qJDD(1) * t533 + t557;
t555 = t534 * g(1) - t538 * g(2);
t549 = -qJDD(1) * pkin(1) - t555;
t544 = pkin(2) * t556 + t522 * t573 + (-t508 - t557) * qJ(3) + t549;
t439 = -t514 * t522 + (-pkin(3) * t528 - pkin(7)) * t540 + (-pkin(2) - qJ(4)) * t509 + t544;
t474 = -t537 * g(3) - t533 * t497;
t461 = -qJDD(2) * pkin(2) - t539 * qJ(3) + t505 * t522 + qJDD(3) - t474;
t454 = (-t533 * t537 * t540 - qJDD(2)) * qJ(4) + (t508 - t557) * pkin(3) + t461;
t529 = sin(pkin(10));
t530 = cos(pkin(10));
t502 = qJD(2) * t530 - t529 * t560;
t428 = -0.2e1 * qJD(4) * t502 - t529 * t439 + t530 * t454;
t480 = qJDD(2) * t530 - t509 * t529;
t501 = -qJD(2) * t529 - t530 * t560;
t419 = (t501 * t522 - t480) * pkin(8) + (t501 * t502 + t508) * pkin(4) + t428;
t429 = 0.2e1 * qJD(4) * t501 + t530 * t439 + t529 * t454;
t479 = -qJDD(2) * t529 - t509 * t530;
t481 = pkin(4) * t522 - pkin(8) * t502;
t500 = t501 ^ 2;
t421 = -pkin(4) * t500 + pkin(8) * t479 - t481 * t522 + t429;
t532 = sin(qJ(5));
t536 = cos(qJ(5));
t413 = t536 * t419 - t532 * t421;
t471 = t501 * t536 - t502 * t532;
t446 = qJD(5) * t471 + t479 * t532 + t480 * t536;
t472 = t501 * t532 + t502 * t536;
t504 = qJDD(5) + t508;
t519 = t522 + qJD(5);
t410 = (t471 * t519 - t446) * pkin(9) + (t471 * t472 + t504) * pkin(5) + t413;
t414 = t532 * t419 + t536 * t421;
t445 = -qJD(5) * t472 + t479 * t536 - t480 * t532;
t464 = pkin(5) * t519 - pkin(9) * t472;
t470 = t471 ^ 2;
t411 = -pkin(5) * t470 + pkin(9) * t445 - t464 * t519 + t414;
t531 = sin(qJ(6));
t535 = cos(qJ(6));
t408 = t410 * t535 - t411 * t531;
t456 = t471 * t535 - t472 * t531;
t426 = qJD(6) * t456 + t445 * t531 + t446 * t535;
t457 = t471 * t531 + t472 * t535;
t436 = -mrSges(7,1) * t456 + mrSges(7,2) * t457;
t517 = qJD(6) + t519;
t440 = -mrSges(7,2) * t517 + mrSges(7,3) * t456;
t499 = qJDD(6) + t504;
t404 = m(7) * t408 + mrSges(7,1) * t499 - mrSges(7,3) * t426 - t436 * t457 + t440 * t517;
t409 = t410 * t531 + t411 * t535;
t425 = -qJD(6) * t457 + t445 * t535 - t446 * t531;
t441 = mrSges(7,1) * t517 - mrSges(7,3) * t457;
t405 = m(7) * t409 - mrSges(7,2) * t499 + mrSges(7,3) * t425 + t436 * t456 - t441 * t517;
t398 = t535 * t404 + t531 * t405;
t458 = -mrSges(6,1) * t471 + mrSges(6,2) * t472;
t462 = -mrSges(6,2) * t519 + mrSges(6,3) * t471;
t395 = m(6) * t413 + mrSges(6,1) * t504 - mrSges(6,3) * t446 - t458 * t472 + t462 * t519 + t398;
t463 = mrSges(6,1) * t519 - mrSges(6,3) * t472;
t552 = -t404 * t531 + t535 * t405;
t396 = m(6) * t414 - mrSges(6,2) * t504 + mrSges(6,3) * t445 + t458 * t471 - t463 * t519 + t552;
t391 = t536 * t395 + t532 * t396;
t564 = (t570 * qJD(2)) + (t533 * t566 + t537 * t565) * qJD(1);
t563 = t565 * qJD(2) + (t533 * t567 + t537 * t571) * qJD(1);
t562 = t566 * qJD(2) + (t572 * t533 + t537 * t567) * qJD(1);
t515 = -mrSges(4,1) * t560 - (qJD(2) * mrSges(4,3));
t561 = -(qJD(2) * mrSges(3,2)) + mrSges(3,3) * t560 - t515;
t473 = -mrSges(5,1) * t501 + mrSges(5,2) * t502;
t477 = -mrSges(5,2) * t522 + mrSges(5,3) * t501;
t389 = m(5) * t428 + mrSges(5,1) * t508 - mrSges(5,3) * t480 - t473 * t502 + t477 * t522 + t391;
t478 = mrSges(5,1) * t522 - mrSges(5,3) * t502;
t553 = -t395 * t532 + t536 * t396;
t390 = m(5) * t429 - mrSges(5,2) * t508 + mrSges(5,3) * t479 + t473 * t501 - t478 * t522 + t553;
t554 = -t529 * t389 + t530 * t390;
t385 = t530 * t389 + t529 * t390;
t459 = -t509 * pkin(2) + t544 - t569;
t550 = m(4) * t459 + t554;
t450 = -qJ(4) * t528 * t540 + pkin(3) * t509 + qJD(2) * t514 + qJDD(4) - t460;
t435 = -pkin(4) * t479 - pkin(8) * t500 + t502 * t481 + t450;
t416 = -pkin(5) * t445 - pkin(9) * t470 + t464 * t472 + t435;
t547 = m(7) * t416 - t425 * mrSges(7,1) + t426 * mrSges(7,2) - t456 * t440 + t457 * t441;
t546 = m(4) * t461 + t508 * mrSges(4,1) + t385;
t431 = Ifges(7,4) * t457 + Ifges(7,2) * t456 + Ifges(7,6) * t517;
t432 = Ifges(7,1) * t457 + Ifges(7,4) * t456 + Ifges(7,5) * t517;
t545 = mrSges(7,1) * t408 - mrSges(7,2) * t409 + Ifges(7,5) * t426 + Ifges(7,6) * t425 + Ifges(7,3) * t499 + t457 * t431 - t456 * t432;
t543 = m(6) * t435 - t445 * mrSges(6,1) + t446 * mrSges(6,2) - t471 * t462 + t472 * t463 + t547;
t452 = Ifges(6,4) * t472 + Ifges(6,2) * t471 + Ifges(6,6) * t519;
t453 = Ifges(6,1) * t472 + Ifges(6,4) * t471 + Ifges(6,5) * t519;
t542 = mrSges(6,1) * t413 - mrSges(6,2) * t414 + Ifges(6,5) * t446 + Ifges(6,6) * t445 + Ifges(6,3) * t504 + pkin(5) * t398 + t472 * t452 - t471 * t453 + t545;
t406 = m(5) * t450 - t479 * mrSges(5,1) + t480 * mrSges(5,2) - t501 * t477 + t502 * t478 + t543;
t506 = (mrSges(4,2) * t537 - mrSges(4,3) * t533) * qJD(1);
t516 = mrSges(4,1) * t522 + qJD(2) * mrSges(4,2);
t541 = -m(4) * t460 + qJDD(2) * mrSges(4,3) + qJD(2) * t516 + t506 * t560 + t406;
t512 = (qJD(2) * mrSges(3,1)) - mrSges(3,3) * t522;
t507 = (-mrSges(3,1) * t537 + mrSges(3,2) * t533) * qJD(1);
t496 = t549 - t569;
t467 = Ifges(5,1) * t502 + Ifges(5,4) * t501 + Ifges(5,5) * t522;
t466 = Ifges(5,4) * t502 + Ifges(5,2) * t501 + Ifges(5,6) * t522;
t465 = Ifges(5,5) * t502 + Ifges(5,6) * t501 + Ifges(5,3) * t522;
t451 = Ifges(6,5) * t472 + Ifges(6,6) * t471 + Ifges(6,3) * t519;
t430 = Ifges(7,5) * t457 + Ifges(7,6) * t456 + Ifges(7,3) * t517;
t400 = mrSges(7,2) * t416 - mrSges(7,3) * t408 + Ifges(7,1) * t426 + Ifges(7,4) * t425 + Ifges(7,5) * t499 + t430 * t456 - t431 * t517;
t399 = -mrSges(7,1) * t416 + mrSges(7,3) * t409 + Ifges(7,4) * t426 + Ifges(7,2) * t425 + Ifges(7,6) * t499 - t430 * t457 + t432 * t517;
t387 = mrSges(6,2) * t435 - mrSges(6,3) * t413 + Ifges(6,1) * t446 + Ifges(6,4) * t445 + Ifges(6,5) * t504 - pkin(9) * t398 - t399 * t531 + t400 * t535 + t451 * t471 - t452 * t519;
t386 = -mrSges(6,1) * t435 + mrSges(6,3) * t414 + Ifges(6,4) * t446 + Ifges(6,2) * t445 + Ifges(6,6) * t504 - pkin(5) * t547 + pkin(9) * t552 + t535 * t399 + t531 * t400 - t472 * t451 + t519 * t453;
t384 = qJDD(2) * mrSges(4,2) + qJD(2) * t515 + t506 * t522 + t546;
t383 = t509 * mrSges(4,2) - t508 * mrSges(4,3) + (t515 * t537 - t516 * t533) * qJD(1) + t550;
t382 = mrSges(5,2) * t450 - mrSges(5,3) * t428 + Ifges(5,1) * t480 + Ifges(5,4) * t479 + Ifges(5,5) * t508 - pkin(8) * t391 - t386 * t532 + t387 * t536 + t465 * t501 - t466 * t522;
t381 = -mrSges(5,1) * t450 + mrSges(5,3) * t429 + Ifges(5,4) * t480 + Ifges(5,2) * t479 + Ifges(5,6) * t508 - pkin(4) * t543 + pkin(8) * t553 + t536 * t386 + t532 * t387 - t502 * t465 + t467 * t522;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t555 - mrSges(2,2) * t551 + t533 * (t542 + t567 * t509 - t563 * qJD(2) + t564 * t560 + (Ifges(5,3) + t572) * t508 + t502 * t466 + mrSges(3,2) * t496 - t501 * t467 + Ifges(5,6) * t479 + Ifges(5,5) * t480 - mrSges(3,3) * t474 + mrSges(4,1) * t461 + t566 * qJDD(2) - mrSges(4,3) * t459 - mrSges(5,2) * t429 + mrSges(5,1) * t428 + pkin(4) * t391 + pkin(3) * t385 - qJ(3) * t383) + t537 * (-mrSges(3,1) * t496 - mrSges(4,1) * t460 + mrSges(4,2) * t459 + mrSges(3,3) * t475 - pkin(2) * t383 + pkin(3) * t406 - qJ(4) * t554 + t562 * qJD(2) + t565 * qJDD(2) - t530 * t381 - t529 * t382 + t567 * t508 + t571 * t509 - t564 * t522) + pkin(1) * (-m(3) * t496 + t568 * t509 + (-mrSges(3,2) + mrSges(4,3)) * t508 + (t561 * t537 + (-t512 + t516) * t533) * qJD(1) - t550) + pkin(7) * (t537 * (t541 - qJDD(2) * mrSges(3,2) + (mrSges(3,3) + mrSges(4,1)) * t509 - qJD(2) * t512 + m(3) * t475 + t507 * t560) + (-m(3) * t474 + t508 * mrSges(3,3) - t568 * qJDD(2) - t561 * qJD(2) + (t506 + t507) * t522 + t546) * t533); mrSges(3,1) * t474 - mrSges(3,2) * t475 + mrSges(4,2) * t461 - mrSges(4,3) * t460 + t530 * t382 - t529 * t381 - qJ(4) * t385 - pkin(2) * t384 + qJ(3) * t541 + (mrSges(4,1) * qJ(3) + t565) * t509 + t566 * t508 + t570 * qJDD(2) + (t563 * t533 - t562 * t537) * qJD(1); t384; t406; t542; t545;];
tauJ  = t1;
