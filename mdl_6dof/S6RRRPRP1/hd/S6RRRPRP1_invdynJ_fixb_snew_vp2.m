% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-05-07 07:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPRP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:25:15
% EndTime: 2019-05-07 07:25:24
% DurationCPUTime: 5.93s
% Computational Cost: add. (60842->328), mult. (136733->404), div. (0->0), fcn. (100354->10), ass. (0->129)
t567 = Ifges(6,1) + Ifges(7,1);
t561 = Ifges(6,4) + Ifges(7,4);
t560 = Ifges(6,5) + Ifges(7,5);
t566 = Ifges(6,2) + Ifges(7,2);
t559 = Ifges(6,6) + Ifges(7,6);
t565 = Ifges(6,3) + Ifges(7,3);
t527 = sin(qJ(2));
t531 = cos(qJ(2));
t549 = qJD(1) * qJD(2);
t511 = qJDD(1) * t527 + t531 * t549;
t533 = qJD(1) ^ 2;
t528 = sin(qJ(1));
t532 = cos(qJ(1));
t540 = -g(1) * t532 - g(2) * t528;
t508 = -pkin(1) * t533 + qJDD(1) * pkin(7) + t540;
t557 = t508 * t527;
t563 = pkin(2) * t533;
t473 = qJDD(2) * pkin(2) - pkin(8) * t511 - t557 + (pkin(8) * t549 + t527 * t563 - g(3)) * t531;
t496 = -g(3) * t527 + t531 * t508;
t512 = qJDD(1) * t531 - t527 * t549;
t552 = qJD(1) * t527;
t515 = qJD(2) * pkin(2) - pkin(8) * t552;
t522 = t531 ^ 2;
t474 = pkin(8) * t512 - qJD(2) * t515 - t522 * t563 + t496;
t526 = sin(qJ(3));
t530 = cos(qJ(3));
t447 = t530 * t473 - t474 * t526;
t505 = (-t526 * t527 + t530 * t531) * qJD(1);
t481 = qJD(3) * t505 + t511 * t530 + t512 * t526;
t506 = (t526 * t531 + t527 * t530) * qJD(1);
t520 = qJDD(2) + qJDD(3);
t521 = qJD(2) + qJD(3);
t423 = (t505 * t521 - t481) * qJ(4) + (t505 * t506 + t520) * pkin(3) + t447;
t448 = t526 * t473 + t530 * t474;
t480 = -qJD(3) * t506 - t511 * t526 + t512 * t530;
t498 = pkin(3) * t521 - qJ(4) * t506;
t501 = t505 ^ 2;
t426 = -pkin(3) * t501 + qJ(4) * t480 - t498 * t521 + t448;
t523 = sin(pkin(10));
t524 = cos(pkin(10));
t493 = t505 * t523 + t506 * t524;
t417 = -0.2e1 * qJD(4) * t493 + t423 * t524 - t523 * t426;
t455 = t480 * t523 + t481 * t524;
t525 = sin(qJ(5));
t529 = cos(qJ(5));
t478 = -t493 * t525 + t521 * t529;
t431 = qJD(5) * t478 + t455 * t529 + t520 * t525;
t479 = t493 * t529 + t521 * t525;
t456 = -mrSges(7,1) * t478 + mrSges(7,2) * t479;
t492 = t505 * t524 - t506 * t523;
t418 = 0.2e1 * qJD(4) * t492 + t523 * t423 + t524 * t426;
t468 = -pkin(4) * t492 - pkin(9) * t493;
t519 = t521 ^ 2;
t415 = -pkin(4) * t519 + pkin(9) * t520 + t468 * t492 + t418;
t546 = g(1) * t528 - t532 * g(2);
t539 = -qJDD(1) * pkin(1) - t546;
t482 = -pkin(2) * t512 + t515 * t552 + (-pkin(8) * t522 - pkin(7)) * t533 + t539;
t433 = -pkin(3) * t480 - qJ(4) * t501 + t506 * t498 + qJDD(4) + t482;
t454 = t480 * t524 - t481 * t523;
t421 = (-t492 * t521 - t455) * pkin(9) + (t493 * t521 - t454) * pkin(4) + t433;
t410 = -t415 * t525 + t529 * t421;
t453 = qJDD(5) - t454;
t487 = qJD(5) - t492;
t407 = -0.2e1 * qJD(6) * t479 + (t478 * t487 - t431) * qJ(6) + (t478 * t479 + t453) * pkin(5) + t410;
t458 = -mrSges(7,2) * t487 + mrSges(7,3) * t478;
t548 = m(7) * t407 + t453 * mrSges(7,1) + t487 * t458;
t404 = -mrSges(7,3) * t431 - t456 * t479 + t548;
t411 = t529 * t415 + t525 * t421;
t430 = -qJD(5) * t479 - t455 * t525 + t520 * t529;
t460 = pkin(5) * t487 - qJ(6) * t479;
t475 = t478 ^ 2;
t409 = -pkin(5) * t475 + qJ(6) * t430 + 0.2e1 * qJD(6) * t478 - t460 * t487 + t411;
t554 = t561 * t478 + t479 * t567 + t560 * t487;
t555 = -t478 * t566 - t479 * t561 - t487 * t559;
t564 = mrSges(6,1) * t410 + mrSges(7,1) * t407 - mrSges(6,2) * t411 - mrSges(7,2) * t409 + pkin(5) * t404 + t430 * t559 + t431 * t560 + t453 * t565 - t478 * t554 - t479 * t555;
t562 = -mrSges(6,2) - mrSges(7,2);
t467 = -mrSges(5,1) * t492 + mrSges(5,2) * t493;
t484 = mrSges(5,1) * t521 - mrSges(5,3) * t493;
t457 = -mrSges(6,1) * t478 + mrSges(6,2) * t479;
t459 = -mrSges(6,2) * t487 + mrSges(6,3) * t478;
t397 = m(6) * t410 + mrSges(6,1) * t453 + t459 * t487 + (-t456 - t457) * t479 + (-mrSges(6,3) - mrSges(7,3)) * t431 + t548;
t547 = m(7) * t409 + t430 * mrSges(7,3) + t478 * t456;
t461 = mrSges(7,1) * t487 - mrSges(7,3) * t479;
t553 = -mrSges(6,1) * t487 + mrSges(6,3) * t479 - t461;
t400 = m(6) * t411 + mrSges(6,3) * t430 + t453 * t562 + t457 * t478 + t487 * t553 + t547;
t543 = -t397 * t525 + t529 * t400;
t390 = m(5) * t418 - mrSges(5,2) * t520 + mrSges(5,3) * t454 + t467 * t492 - t484 * t521 + t543;
t483 = -mrSges(5,2) * t521 + mrSges(5,3) * t492;
t414 = -pkin(4) * t520 - pkin(9) * t519 + t493 * t468 - t417;
t412 = -pkin(5) * t430 - qJ(6) * t475 + t460 * t479 + qJDD(6) + t414;
t541 = -m(7) * t412 + t430 * mrSges(7,1) + t478 * t458;
t536 = -m(6) * t414 + t430 * mrSges(6,1) + t431 * t562 + t478 * t459 + t479 * t553 + t541;
t402 = m(5) * t417 + mrSges(5,1) * t520 - mrSges(5,3) * t455 - t467 * t493 + t483 * t521 + t536;
t385 = t523 * t390 + t524 * t402;
t494 = -mrSges(4,1) * t505 + mrSges(4,2) * t506;
t497 = -mrSges(4,2) * t521 + mrSges(4,3) * t505;
t382 = m(4) * t447 + mrSges(4,1) * t520 - mrSges(4,3) * t481 - t494 * t506 + t497 * t521 + t385;
t499 = mrSges(4,1) * t521 - mrSges(4,3) * t506;
t544 = t524 * t390 - t402 * t523;
t383 = m(4) * t448 - mrSges(4,2) * t520 + mrSges(4,3) * t480 + t494 * t505 - t499 * t521 + t544;
t377 = t530 * t382 + t526 * t383;
t395 = t529 * t397 + t525 * t400;
t556 = -t478 * t559 - t479 * t560 - t487 * t565;
t551 = qJD(1) * t531;
t545 = -t382 * t526 + t530 * t383;
t538 = -m(5) * t433 + mrSges(5,1) * t454 - t455 * mrSges(5,2) + t483 * t492 - t493 * t484 - t395;
t535 = m(4) * t482 - mrSges(4,1) * t480 + mrSges(4,2) * t481 - t497 * t505 + t499 * t506 - t538;
t405 = mrSges(7,2) * t431 + t461 * t479 - t541;
t387 = -mrSges(6,1) * t414 + mrSges(6,3) * t411 - mrSges(7,1) * t412 + mrSges(7,3) * t409 - pkin(5) * t405 + qJ(6) * t547 + (-qJ(6) * t461 + t554) * t487 + t556 * t479 + (-mrSges(7,2) * qJ(6) + t559) * t453 + t561 * t431 + t566 * t430;
t392 = mrSges(6,2) * t414 + mrSges(7,2) * t412 - mrSges(6,3) * t410 - mrSges(7,3) * t407 - qJ(6) * t404 + t561 * t430 + t431 * t567 + t560 * t453 - t556 * t478 + t555 * t487;
t464 = Ifges(5,4) * t493 + Ifges(5,2) * t492 + Ifges(5,6) * t521;
t465 = Ifges(5,1) * t493 + Ifges(5,4) * t492 + Ifges(5,5) * t521;
t489 = Ifges(4,4) * t506 + Ifges(4,2) * t505 + Ifges(4,6) * t521;
t490 = Ifges(4,1) * t506 + Ifges(4,4) * t505 + Ifges(4,5) * t521;
t534 = mrSges(4,1) * t447 + mrSges(5,1) * t417 - mrSges(4,2) * t448 - mrSges(5,2) * t418 + pkin(3) * t385 + pkin(4) * t536 + pkin(9) * t543 + t529 * t387 + t525 * t392 + t493 * t464 - t492 * t465 - t505 * t490 + Ifges(5,6) * t454 + Ifges(5,5) * t455 + t506 * t489 + Ifges(4,6) * t480 + Ifges(4,5) * t481 + (Ifges(5,3) + Ifges(4,3)) * t520;
t514 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t551;
t513 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t552;
t510 = (-mrSges(3,1) * t531 + mrSges(3,2) * t527) * qJD(1);
t507 = -pkin(7) * t533 + t539;
t504 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t527 + Ifges(3,4) * t531) * qJD(1);
t503 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t527 + Ifges(3,2) * t531) * qJD(1);
t495 = -g(3) * t531 - t557;
t488 = Ifges(4,5) * t506 + Ifges(4,6) * t505 + Ifges(4,3) * t521;
t463 = Ifges(5,5) * t493 + Ifges(5,6) * t492 + Ifges(5,3) * t521;
t378 = -mrSges(5,1) * t433 + mrSges(5,3) * t418 + Ifges(5,4) * t455 + Ifges(5,2) * t454 + Ifges(5,6) * t520 - pkin(4) * t395 - t493 * t463 + t521 * t465 - t564;
t376 = mrSges(5,2) * t433 - mrSges(5,3) * t417 + Ifges(5,1) * t455 + Ifges(5,4) * t454 + Ifges(5,5) * t520 - pkin(9) * t395 - t387 * t525 + t392 * t529 + t463 * t492 - t464 * t521;
t375 = mrSges(4,2) * t482 - mrSges(4,3) * t447 + Ifges(4,1) * t481 + Ifges(4,4) * t480 + Ifges(4,5) * t520 - qJ(4) * t385 + t376 * t524 - t378 * t523 + t488 * t505 - t489 * t521;
t374 = -mrSges(4,1) * t482 + mrSges(4,3) * t448 + Ifges(4,4) * t481 + Ifges(4,2) * t480 + Ifges(4,6) * t520 + pkin(3) * t538 + qJ(4) * t544 + t523 * t376 + t524 * t378 - t506 * t488 + t521 * t490;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t546 - mrSges(2,2) * t540 + t527 * (mrSges(3,2) * t507 - mrSges(3,3) * t495 + Ifges(3,1) * t511 + Ifges(3,4) * t512 + Ifges(3,5) * qJDD(2) - pkin(8) * t377 - qJD(2) * t503 - t526 * t374 + t530 * t375) + t531 * (-mrSges(3,1) * t507 + mrSges(3,3) * t496 + Ifges(3,4) * t511 + Ifges(3,2) * t512 + Ifges(3,6) * qJDD(2) - pkin(2) * t535 + pkin(8) * t545 + qJD(2) * t504 + t530 * t374 + t526 * t375) + pkin(1) * (-t535 + (-t513 * t527 + t514 * t531) * qJD(1) - m(3) * t507 + mrSges(3,1) * t512 - mrSges(3,2) * t511) + pkin(7) * (t531 * (m(3) * t496 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t512 - qJD(2) * t513 + t510 * t551 + t545) - t527 * (m(3) * t495 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t511 + qJD(2) * t514 - t510 * t552 + t377)); t534 + (t527 * t503 - t531 * t504) * qJD(1) + Ifges(3,3) * qJDD(2) + Ifges(3,5) * t511 + Ifges(3,6) * t512 + mrSges(3,1) * t495 - mrSges(3,2) * t496 + pkin(2) * t377; t534; -t538; t564; t405;];
tauJ  = t1;
