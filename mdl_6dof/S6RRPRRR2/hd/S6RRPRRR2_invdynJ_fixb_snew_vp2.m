% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-06 19:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRRR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 19:53:38
% EndTime: 2019-05-06 19:53:51
% DurationCPUTime: 8.86s
% Computational Cost: add. (127752->350), mult. (295211->446), div. (0->0), fcn. (223604->12), ass. (0->138)
t542 = qJD(1) ^ 2;
t560 = pkin(2) * t542;
t536 = sin(qJ(1));
t541 = cos(qJ(1));
t550 = -g(1) * t541 - g(2) * t536;
t516 = -pkin(1) * t542 + qJDD(1) * pkin(7) + t550;
t535 = sin(qJ(2));
t559 = t535 * t516;
t540 = cos(qJ(2));
t556 = qJD(1) * qJD(2);
t519 = qJDD(1) * t535 + t540 * t556;
t479 = qJDD(2) * pkin(2) - t519 * qJ(3) - t559 + (qJ(3) * t556 + t535 * t560 - g(3)) * t540;
t501 = -g(3) * t535 + t540 * t516;
t520 = qJDD(1) * t540 - t535 * t556;
t558 = qJD(1) * t535;
t521 = qJD(2) * pkin(2) - qJ(3) * t558;
t529 = t540 ^ 2;
t480 = qJ(3) * t520 - qJD(2) * t521 - t529 * t560 + t501;
t530 = sin(pkin(11));
t531 = cos(pkin(11));
t511 = (t530 * t540 + t531 * t535) * qJD(1);
t452 = -0.2e1 * qJD(3) * t511 + t531 * t479 - t530 * t480;
t499 = t519 * t531 + t520 * t530;
t510 = (-t530 * t535 + t531 * t540) * qJD(1);
t436 = (qJD(2) * t510 - t499) * pkin(8) + (t510 * t511 + qJDD(2)) * pkin(3) + t452;
t453 = 0.2e1 * qJD(3) * t510 + t530 * t479 + t531 * t480;
t498 = -t519 * t530 + t520 * t531;
t504 = qJD(2) * pkin(3) - pkin(8) * t511;
t509 = t510 ^ 2;
t440 = -pkin(3) * t509 + pkin(8) * t498 - qJD(2) * t504 + t453;
t534 = sin(qJ(4));
t539 = cos(qJ(4));
t429 = t534 * t436 + t539 * t440;
t494 = t510 * t534 + t511 * t539;
t462 = -t494 * qJD(4) + t498 * t539 - t534 * t499;
t493 = t510 * t539 - t534 * t511;
t474 = -mrSges(5,1) * t493 + mrSges(5,2) * t494;
t528 = qJD(2) + qJD(4);
t486 = mrSges(5,1) * t528 - mrSges(5,3) * t494;
t527 = qJDD(2) + qJDD(4);
t475 = -pkin(4) * t493 - pkin(9) * t494;
t526 = t528 ^ 2;
t418 = -pkin(4) * t526 + pkin(9) * t527 + t475 * t493 + t429;
t555 = t536 * g(1) - t541 * g(2);
t549 = -qJDD(1) * pkin(1) - t555;
t482 = -t520 * pkin(2) + qJDD(3) + t521 * t558 + (-qJ(3) * t529 - pkin(7)) * t542 + t549;
t451 = -t498 * pkin(3) - t509 * pkin(8) + t511 * t504 + t482;
t463 = qJD(4) * t493 + t498 * t534 + t499 * t539;
t426 = (-t493 * t528 - t463) * pkin(9) + (t494 * t528 - t462) * pkin(4) + t451;
t533 = sin(qJ(5));
t538 = cos(qJ(5));
t413 = -t533 * t418 + t538 * t426;
t483 = -t494 * t533 + t528 * t538;
t443 = qJD(5) * t483 + t463 * t538 + t527 * t533;
t461 = qJDD(5) - t462;
t484 = t494 * t538 + t528 * t533;
t489 = qJD(5) - t493;
t411 = (t483 * t489 - t443) * pkin(10) + (t483 * t484 + t461) * pkin(5) + t413;
t414 = t538 * t418 + t533 * t426;
t442 = -qJD(5) * t484 - t463 * t533 + t527 * t538;
t469 = pkin(5) * t489 - pkin(10) * t484;
t481 = t483 ^ 2;
t412 = -pkin(5) * t481 + pkin(10) * t442 - t469 * t489 + t414;
t532 = sin(qJ(6));
t537 = cos(qJ(6));
t409 = t411 * t537 - t412 * t532;
t464 = t483 * t537 - t484 * t532;
t423 = qJD(6) * t464 + t442 * t532 + t443 * t537;
t465 = t483 * t532 + t484 * t537;
t437 = -mrSges(7,1) * t464 + mrSges(7,2) * t465;
t487 = qJD(6) + t489;
t444 = -mrSges(7,2) * t487 + mrSges(7,3) * t464;
t457 = qJDD(6) + t461;
t405 = m(7) * t409 + mrSges(7,1) * t457 - mrSges(7,3) * t423 - t437 * t465 + t444 * t487;
t410 = t411 * t532 + t412 * t537;
t422 = -qJD(6) * t465 + t442 * t537 - t443 * t532;
t445 = mrSges(7,1) * t487 - mrSges(7,3) * t465;
t406 = m(7) * t410 - mrSges(7,2) * t457 + mrSges(7,3) * t422 + t437 * t464 - t445 * t487;
t397 = t537 * t405 + t532 * t406;
t466 = -mrSges(6,1) * t483 + mrSges(6,2) * t484;
t467 = -mrSges(6,2) * t489 + mrSges(6,3) * t483;
t395 = m(6) * t413 + mrSges(6,1) * t461 - mrSges(6,3) * t443 - t466 * t484 + t467 * t489 + t397;
t468 = mrSges(6,1) * t489 - mrSges(6,3) * t484;
t551 = -t405 * t532 + t537 * t406;
t396 = m(6) * t414 - mrSges(6,2) * t461 + mrSges(6,3) * t442 + t466 * t483 - t468 * t489 + t551;
t552 = -t395 * t533 + t538 * t396;
t388 = m(5) * t429 - mrSges(5,2) * t527 + mrSges(5,3) * t462 + t474 * t493 - t486 * t528 + t552;
t428 = t436 * t539 - t534 * t440;
t485 = -mrSges(5,2) * t528 + mrSges(5,3) * t493;
t417 = -pkin(4) * t527 - pkin(9) * t526 + t494 * t475 - t428;
t415 = -pkin(5) * t442 - pkin(10) * t481 + t469 * t484 + t417;
t548 = m(7) * t415 - t422 * mrSges(7,1) + mrSges(7,2) * t423 - t464 * t444 + t445 * t465;
t544 = -m(6) * t417 + t442 * mrSges(6,1) - mrSges(6,2) * t443 + t483 * t467 - t468 * t484 - t548;
t401 = m(5) * t428 + mrSges(5,1) * t527 - mrSges(5,3) * t463 - t474 * t494 + t485 * t528 + t544;
t381 = t534 * t388 + t539 * t401;
t496 = -mrSges(4,1) * t510 + mrSges(4,2) * t511;
t502 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t510;
t379 = m(4) * t452 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t499 + qJD(2) * t502 - t496 * t511 + t381;
t503 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t511;
t553 = t539 * t388 - t401 * t534;
t380 = m(4) * t453 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t498 - qJD(2) * t503 + t496 * t510 + t553;
t374 = t531 * t379 + t530 * t380;
t391 = t538 * t395 + t533 * t396;
t557 = qJD(1) * t540;
t554 = -t379 * t530 + t531 * t380;
t547 = -m(5) * t451 + t462 * mrSges(5,1) - t463 * mrSges(5,2) + t493 * t485 - t494 * t486 - t391;
t432 = Ifges(7,4) * t465 + Ifges(7,2) * t464 + Ifges(7,6) * t487;
t433 = Ifges(7,1) * t465 + Ifges(7,4) * t464 + Ifges(7,5) * t487;
t546 = -mrSges(7,1) * t409 + mrSges(7,2) * t410 - Ifges(7,5) * t423 - Ifges(7,6) * t422 - Ifges(7,3) * t457 - t465 * t432 + t464 * t433;
t431 = Ifges(7,5) * t465 + Ifges(7,6) * t464 + Ifges(7,3) * t487;
t398 = -mrSges(7,1) * t415 + mrSges(7,3) * t410 + Ifges(7,4) * t423 + Ifges(7,2) * t422 + Ifges(7,6) * t457 - t431 * t465 + t433 * t487;
t399 = mrSges(7,2) * t415 - mrSges(7,3) * t409 + Ifges(7,1) * t423 + Ifges(7,4) * t422 + Ifges(7,5) * t457 + t431 * t464 - t432 * t487;
t446 = Ifges(6,5) * t484 + Ifges(6,6) * t483 + Ifges(6,3) * t489;
t448 = Ifges(6,1) * t484 + Ifges(6,4) * t483 + Ifges(6,5) * t489;
t383 = -mrSges(6,1) * t417 + mrSges(6,3) * t414 + Ifges(6,4) * t443 + Ifges(6,2) * t442 + Ifges(6,6) * t461 - pkin(5) * t548 + pkin(10) * t551 + t537 * t398 + t532 * t399 - t484 * t446 + t489 * t448;
t447 = Ifges(6,4) * t484 + Ifges(6,2) * t483 + Ifges(6,6) * t489;
t385 = mrSges(6,2) * t417 - mrSges(6,3) * t413 + Ifges(6,1) * t443 + Ifges(6,4) * t442 + Ifges(6,5) * t461 - pkin(10) * t397 - t398 * t532 + t399 * t537 + t446 * t483 - t447 * t489;
t471 = Ifges(5,4) * t494 + Ifges(5,2) * t493 + Ifges(5,6) * t528;
t472 = Ifges(5,1) * t494 + Ifges(5,4) * t493 + Ifges(5,5) * t528;
t545 = mrSges(5,1) * t428 - mrSges(5,2) * t429 + Ifges(5,5) * t463 + Ifges(5,6) * t462 + Ifges(5,3) * t527 + pkin(4) * t544 + pkin(9) * t552 + t538 * t383 + t533 * t385 + t494 * t471 - t493 * t472;
t389 = m(4) * t482 - t498 * mrSges(4,1) + t499 * mrSges(4,2) - t510 * t502 + t511 * t503 - t547;
t543 = mrSges(6,1) * t413 - mrSges(6,2) * t414 + Ifges(6,5) * t443 + Ifges(6,6) * t442 + Ifges(6,3) * t461 + pkin(5) * t397 + t484 * t447 - t483 * t448 - t546;
t523 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t557;
t522 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t558;
t518 = (-mrSges(3,1) * t540 + mrSges(3,2) * t535) * qJD(1);
t515 = -t542 * pkin(7) + t549;
t514 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t535 + Ifges(3,4) * t540) * qJD(1);
t513 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t535 + Ifges(3,2) * t540) * qJD(1);
t500 = -t540 * g(3) - t559;
t492 = Ifges(4,1) * t511 + Ifges(4,4) * t510 + Ifges(4,5) * qJD(2);
t491 = Ifges(4,4) * t511 + Ifges(4,2) * t510 + Ifges(4,6) * qJD(2);
t490 = Ifges(4,5) * t511 + Ifges(4,6) * t510 + Ifges(4,3) * qJD(2);
t470 = Ifges(5,5) * t494 + Ifges(5,6) * t493 + Ifges(5,3) * t528;
t375 = -mrSges(5,1) * t451 + mrSges(5,3) * t429 + Ifges(5,4) * t463 + Ifges(5,2) * t462 + Ifges(5,6) * t527 - pkin(4) * t391 - t494 * t470 + t528 * t472 - t543;
t373 = mrSges(5,2) * t451 - mrSges(5,3) * t428 + Ifges(5,1) * t463 + Ifges(5,4) * t462 + Ifges(5,5) * t527 - pkin(9) * t391 - t383 * t533 + t385 * t538 + t470 * t493 - t471 * t528;
t372 = mrSges(4,2) * t482 - mrSges(4,3) * t452 + Ifges(4,1) * t499 + Ifges(4,4) * t498 + Ifges(4,5) * qJDD(2) - pkin(8) * t381 - qJD(2) * t491 + t373 * t539 - t375 * t534 + t490 * t510;
t371 = -mrSges(4,1) * t482 + mrSges(4,3) * t453 + Ifges(4,4) * t499 + Ifges(4,2) * t498 + Ifges(4,6) * qJDD(2) + pkin(3) * t547 + pkin(8) * t553 + qJD(2) * t492 + t534 * t373 + t539 * t375 - t511 * t490;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t555 - mrSges(2,2) * t550 + t535 * (mrSges(3,2) * t515 - mrSges(3,3) * t500 + Ifges(3,1) * t519 + Ifges(3,4) * t520 + Ifges(3,5) * qJDD(2) - qJ(3) * t374 - qJD(2) * t513 - t530 * t371 + t531 * t372) + t540 * (-mrSges(3,1) * t515 + mrSges(3,3) * t501 + Ifges(3,4) * t519 + Ifges(3,2) * t520 + Ifges(3,6) * qJDD(2) - pkin(2) * t389 + qJ(3) * t554 + qJD(2) * t514 + t531 * t371 + t530 * t372) + pkin(1) * (-t389 - m(3) * t515 - t519 * mrSges(3,2) + t520 * mrSges(3,1) + (-t522 * t535 + t523 * t540) * qJD(1)) + pkin(7) * (t540 * (m(3) * t501 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t520 - qJD(2) * t522 + t518 * t557 + t554) - t535 * (m(3) * t500 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t519 + qJD(2) * t523 - t518 * t558 + t374)); pkin(3) * t381 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + (t535 * t513 - t540 * t514) * qJD(1) + pkin(2) * t374 + t545 + Ifges(3,5) * t519 + Ifges(3,6) * t520 - t510 * t492 + t511 * t491 + Ifges(4,6) * t498 + Ifges(4,5) * t499 + mrSges(3,1) * t500 - mrSges(3,2) * t501 + mrSges(4,1) * t452 - mrSges(4,2) * t453; t389; t545; t543; -t546;];
tauJ  = t1;
