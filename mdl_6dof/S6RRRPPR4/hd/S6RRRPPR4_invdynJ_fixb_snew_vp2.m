% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-05-07 04:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPPR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:41:39
% EndTime: 2019-05-07 04:41:53
% DurationCPUTime: 5.34s
% Computational Cost: add. (51071->328), mult. (106230->398), div. (0->0), fcn. (74660->10), ass. (0->132)
t579 = Ifges(5,1) + Ifges(6,1);
t568 = Ifges(5,4) - Ifges(6,5);
t567 = Ifges(5,5) + Ifges(6,4);
t578 = -Ifges(5,2) - Ifges(6,3);
t566 = Ifges(6,6) - Ifges(5,6);
t577 = -2 * qJD(4);
t576 = -Ifges(6,2) - Ifges(5,3);
t537 = qJD(1) ^ 2;
t531 = sin(qJ(1));
t535 = cos(qJ(1));
t552 = t531 * g(1) - t535 * g(2);
t504 = -qJDD(1) * pkin(1) - t537 * pkin(7) - t552;
t530 = sin(qJ(2));
t534 = cos(qJ(2));
t556 = qJD(1) * qJD(2);
t553 = t534 * t556;
t516 = qJDD(1) * t530 + t553;
t554 = t530 * t556;
t517 = t534 * qJDD(1) - t554;
t465 = (-t516 - t553) * pkin(8) + (-t517 + t554) * pkin(2) + t504;
t547 = -g(1) * t535 - g(2) * t531;
t505 = -pkin(1) * t537 + qJDD(1) * pkin(7) + t547;
t495 = -g(3) * t530 + t534 * t505;
t515 = (-pkin(2) * t534 - pkin(8) * t530) * qJD(1);
t536 = qJD(2) ^ 2;
t558 = qJD(1) * t534;
t470 = -pkin(2) * t536 + qJDD(2) * pkin(8) + t515 * t558 + t495;
t529 = sin(qJ(3));
t533 = cos(qJ(3));
t438 = t533 * t465 - t529 * t470;
t559 = qJD(1) * t530;
t512 = qJD(2) * t533 - t529 * t559;
t486 = qJD(3) * t512 + qJDD(2) * t529 + t516 * t533;
t511 = qJDD(3) - t517;
t513 = qJD(2) * t529 + t533 * t559;
t523 = qJD(3) - t558;
t427 = (t512 * t523 - t486) * qJ(4) + (t512 * t513 + t511) * pkin(3) + t438;
t439 = t529 * t465 + t533 * t470;
t485 = -qJD(3) * t513 + qJDD(2) * t533 - t516 * t529;
t492 = pkin(3) * t523 - qJ(4) * t513;
t510 = t512 ^ 2;
t430 = -pkin(3) * t510 + qJ(4) * t485 - t492 * t523 + t439;
t527 = sin(pkin(10));
t565 = cos(pkin(10));
t488 = -t512 * t565 + t513 * t527;
t416 = t527 * t427 + t565 * t430 + t488 * t577;
t453 = -t485 * t565 + t486 * t527;
t489 = t527 * t512 + t513 * t565;
t472 = mrSges(5,1) * t523 - mrSges(5,3) * t489;
t460 = pkin(4) * t488 - qJ(5) * t489;
t522 = t523 ^ 2;
t570 = 2 * qJD(5);
t413 = -pkin(4) * t522 + t511 * qJ(5) - t488 * t460 + t523 * t570 + t416;
t473 = -mrSges(6,1) * t523 + mrSges(6,2) * t489;
t415 = t427 * t565 - t527 * t430 + t489 * t577;
t414 = -t511 * pkin(4) - t522 * qJ(5) + t489 * t460 + qJDD(5) - t415;
t454 = t527 * t485 + t486 * t565;
t564 = t488 * t523;
t408 = (-t454 - t564) * pkin(9) + (t488 * t489 - t511) * pkin(5) + t414;
t475 = -pkin(5) * t523 - pkin(9) * t489;
t487 = t488 ^ 2;
t409 = -pkin(5) * t487 + pkin(9) * t453 + t475 * t523 + t413;
t528 = sin(qJ(6));
t532 = cos(qJ(6));
t406 = t408 * t532 - t409 * t528;
t458 = t488 * t532 - t489 * t528;
t424 = qJD(6) * t458 + t453 * t528 + t454 * t532;
t459 = t488 * t528 + t489 * t532;
t436 = -mrSges(7,1) * t458 + mrSges(7,2) * t459;
t521 = qJD(6) - t523;
t442 = -mrSges(7,2) * t521 + mrSges(7,3) * t458;
t509 = qJDD(6) - t511;
t403 = m(7) * t406 + mrSges(7,1) * t509 - mrSges(7,3) * t424 - t436 * t459 + t442 * t521;
t407 = t408 * t528 + t409 * t532;
t423 = -qJD(6) * t459 + t453 * t532 - t454 * t528;
t443 = mrSges(7,1) * t521 - mrSges(7,3) * t459;
t404 = m(7) * t407 - mrSges(7,2) * t509 + mrSges(7,3) * t423 + t436 * t458 - t443 * t521;
t549 = -t528 * t403 + t532 * t404;
t544 = m(6) * t413 + t511 * mrSges(6,3) + t523 * t473 + t549;
t461 = mrSges(6,1) * t488 - mrSges(6,3) * t489;
t560 = -mrSges(5,1) * t488 - mrSges(5,2) * t489 - t461;
t569 = -mrSges(5,3) - mrSges(6,2);
t392 = m(5) * t416 - t511 * mrSges(5,2) + t453 * t569 - t523 * t472 + t488 * t560 + t544;
t471 = -mrSges(5,2) * t523 - mrSges(5,3) * t488;
t397 = t532 * t403 + t528 * t404;
t474 = -mrSges(6,2) * t488 + mrSges(6,3) * t523;
t542 = -m(6) * t414 + t511 * mrSges(6,1) + t523 * t474 - t397;
t394 = m(5) * t415 + t511 * mrSges(5,1) + t454 * t569 + t523 * t471 + t489 * t560 + t542;
t389 = t527 * t392 + t565 * t394;
t396 = t454 * mrSges(6,2) + t489 * t461 - t542;
t478 = Ifges(4,4) * t513 + Ifges(4,2) * t512 + Ifges(4,6) * t523;
t479 = Ifges(4,1) * t513 + Ifges(4,4) * t512 + Ifges(4,5) * t523;
t432 = Ifges(7,4) * t459 + Ifges(7,2) * t458 + Ifges(7,6) * t521;
t433 = Ifges(7,1) * t459 + Ifges(7,4) * t458 + Ifges(7,5) * t521;
t541 = -mrSges(7,1) * t406 + mrSges(7,2) * t407 - Ifges(7,5) * t424 - Ifges(7,6) * t423 - Ifges(7,3) * t509 - t459 * t432 + t458 * t433;
t561 = -t568 * t488 + t579 * t489 + t567 * t523;
t562 = t578 * t488 + t568 * t489 - t566 * t523;
t575 = t453 * t566 + t454 * t567 + mrSges(4,1) * t438 + mrSges(5,1) * t415 - mrSges(6,1) * t414 - mrSges(4,2) * t439 - mrSges(5,2) * t416 + mrSges(6,3) * t413 + Ifges(4,5) * t486 + Ifges(4,6) * t485 + pkin(3) * t389 - pkin(4) * t396 - pkin(5) * t397 + qJ(5) * (-t453 * mrSges(6,2) - t488 * t461 + t544) + t513 * t478 - t512 * t479 + t541 + t489 * t562 + t488 * t561;
t572 = (Ifges(4,3) - t576) * t511;
t494 = -t534 * g(3) - t530 * t505;
t543 = qJDD(2) * pkin(2) + t536 * pkin(8) - t515 * t559 + t494;
t540 = t485 * pkin(3) + t510 * qJ(4) - t513 * t492 - qJDD(4) + t543;
t571 = (-t454 + t564) * qJ(5) - t540;
t563 = -t566 * t488 - t567 * t489 + t576 * t523;
t490 = -mrSges(4,1) * t512 + mrSges(4,2) * t513;
t491 = -mrSges(4,2) * t523 + mrSges(4,3) * t512;
t387 = m(4) * t438 + mrSges(4,1) * t511 - mrSges(4,3) * t486 - t490 * t513 + t491 * t523 + t389;
t493 = mrSges(4,1) * t523 - mrSges(4,3) * t513;
t550 = t565 * t392 - t394 * t527;
t388 = m(4) * t439 - mrSges(4,2) * t511 + mrSges(4,3) * t485 + t490 * t512 - t493 * t523 + t550;
t551 = -t387 * t529 + t533 * t388;
t411 = -t487 * pkin(9) + (-pkin(4) - pkin(5)) * t453 + (-pkin(4) * t523 + t475 + t570) * t489 - t571;
t545 = -m(7) * t411 + t423 * mrSges(7,1) - t424 * mrSges(7,2) + t458 * t442 - t459 * t443;
t383 = t533 * t387 + t529 * t388;
t418 = -0.2e1 * qJD(5) * t489 + (t489 * t523 + t453) * pkin(4) + t571;
t401 = m(6) * t418 + t453 * mrSges(6,1) - t454 * mrSges(6,3) - t489 * t473 + t488 * t474 + t545;
t400 = -m(5) * t540 + t453 * mrSges(5,1) + t454 * mrSges(5,2) + t488 * t471 + t489 * t472 + t401;
t539 = m(4) * t543 + t485 * mrSges(4,1) - t486 * mrSges(4,2) + t512 * t491 - t513 * t493 - t400;
t519 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t558;
t518 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t559;
t514 = (-mrSges(3,1) * t534 + mrSges(3,2) * t530) * qJD(1);
t503 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t530 + Ifges(3,4) * t534) * qJD(1);
t502 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t530 + Ifges(3,2) * t534) * qJD(1);
t477 = Ifges(4,5) * t513 + Ifges(4,6) * t512 + Ifges(4,3) * t523;
t431 = Ifges(7,5) * t459 + Ifges(7,6) * t458 + Ifges(7,3) * t521;
t399 = mrSges(7,2) * t411 - mrSges(7,3) * t406 + Ifges(7,1) * t424 + Ifges(7,4) * t423 + Ifges(7,5) * t509 + t431 * t458 - t432 * t521;
t398 = -mrSges(7,1) * t411 + mrSges(7,3) * t407 + Ifges(7,4) * t424 + Ifges(7,2) * t423 + Ifges(7,6) * t509 - t431 * t459 + t433 * t521;
t385 = -mrSges(5,2) * t540 + mrSges(6,2) * t414 - mrSges(5,3) * t415 - mrSges(6,3) * t418 - pkin(9) * t397 - qJ(5) * t401 - t528 * t398 + t532 * t399 - t568 * t453 + t579 * t454 + t563 * t488 + t567 * t511 - t562 * t523;
t384 = mrSges(5,1) * t540 - mrSges(6,1) * t418 + mrSges(6,2) * t413 + mrSges(5,3) * t416 - pkin(4) * t401 - pkin(5) * t545 - pkin(9) * t549 - t532 * t398 - t528 * t399 + t578 * t453 + t568 * t454 + t563 * t489 - t566 * t511 + t561 * t523;
t382 = -mrSges(4,2) * t543 - mrSges(4,3) * t438 + Ifges(4,1) * t486 + Ifges(4,4) * t485 + Ifges(4,5) * t511 - qJ(4) * t389 - t527 * t384 + t385 * t565 + t512 * t477 - t523 * t478;
t381 = mrSges(4,1) * t543 + mrSges(4,3) * t439 + Ifges(4,4) * t486 + Ifges(4,2) * t485 + Ifges(4,6) * t511 - pkin(3) * t400 + qJ(4) * t550 + t384 * t565 + t527 * t385 - t513 * t477 + t523 * t479;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t552 - mrSges(2,2) * t547 + t530 * (mrSges(3,2) * t504 - mrSges(3,3) * t494 + Ifges(3,1) * t516 + Ifges(3,4) * t517 + Ifges(3,5) * qJDD(2) - pkin(8) * t383 - qJD(2) * t502 - t529 * t381 + t533 * t382) + t534 * (-mrSges(3,1) * t504 + mrSges(3,3) * t495 + Ifges(3,4) * t516 + Ifges(3,2) * t517 + Ifges(3,6) * qJDD(2) - pkin(2) * t383 + qJD(2) * t503 - t575) + pkin(1) * (-m(3) * t504 + t517 * mrSges(3,1) - t516 * mrSges(3,2) + (-t518 * t530 + t519 * t534) * qJD(1) - t383) + pkin(7) * (t534 * (m(3) * t495 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t517 - qJD(2) * t518 + t514 * t558 + t551) - t530 * (m(3) * t494 + qJDD(2) * mrSges(3,1) - t516 * mrSges(3,3) + qJD(2) * t519 - t514 * t559 + t539)) - t534 * t572; Ifges(3,5) * t516 + Ifges(3,6) * t517 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t494 - mrSges(3,2) * t495 + t529 * t382 + t533 * t381 + pkin(2) * t539 + pkin(8) * t551 + (t530 * t502 - t534 * t503) * qJD(1); t572 + t575; t400; t396; -t541;];
tauJ  = t1;
