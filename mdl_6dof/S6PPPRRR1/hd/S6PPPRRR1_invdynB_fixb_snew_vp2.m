% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PPPRRR1
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
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
% Datum: 2019-05-04 19:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PPPRRR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPPRRR1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPPRRR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_invdynB_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPPRRR1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPPRRR1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPPRRR1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:32:54
% EndTime: 2019-05-04 19:33:25
% DurationCPUTime: 30.30s
% Computational Cost: add. (580061->248), mult. (982366->336), div. (0->0), fcn. (872348->18), ass. (0->126)
t499 = sin(pkin(12));
t505 = cos(pkin(12));
t492 = -g(1) * t505 - g(2) * t499;
t498 = sin(pkin(13));
t504 = cos(pkin(13));
t491 = g(1) * t499 - g(2) * t505;
t496 = -g(3) + qJDD(1);
t502 = sin(pkin(6));
t508 = cos(pkin(6));
t522 = t491 * t508 + t496 * t502;
t467 = t504 * t492 + t498 * t522;
t497 = sin(pkin(14));
t503 = cos(pkin(14));
t466 = -t498 * t492 + t504 * t522;
t478 = -t491 * t502 + t496 * t508 + qJDD(2);
t501 = sin(pkin(7));
t507 = cos(pkin(7));
t523 = t466 * t507 + t478 * t501;
t462 = -t497 * t467 + t503 * t523;
t463 = t503 * t467 + t497 * t523;
t465 = -t466 * t501 + t478 * t507 + qJDD(3);
t514 = cos(qJ(4));
t506 = cos(pkin(8));
t511 = sin(qJ(4));
t535 = t506 * t511;
t500 = sin(pkin(8));
t536 = t500 * t511;
t456 = t462 * t535 + t463 * t514 + t465 * t536;
t516 = qJD(4) ^ 2;
t454 = -pkin(4) * t516 + qJDD(4) * pkin(10) + t456;
t458 = -t462 * t500 + t465 * t506;
t510 = sin(qJ(5));
t513 = cos(qJ(5));
t450 = t454 * t513 + t458 * t510;
t487 = (-mrSges(6,1) * t513 + mrSges(6,2) * t510) * qJD(4);
t530 = qJD(4) * qJD(5);
t529 = t510 * t530;
t490 = qJDD(4) * t513 - t529;
t532 = qJD(4) * t510;
t493 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t532;
t488 = (-pkin(5) * t513 - pkin(11) * t510) * qJD(4);
t515 = qJD(5) ^ 2;
t531 = t513 * qJD(4);
t448 = -pkin(5) * t515 + qJDD(5) * pkin(11) + t488 * t531 + t450;
t455 = -t511 * t463 + (t462 * t506 + t465 * t500) * t514;
t453 = -qJDD(4) * pkin(4) - t516 * pkin(10) - t455;
t528 = t513 * t530;
t489 = qJDD(4) * t510 + t528;
t451 = (-t489 - t528) * pkin(11) + (-t490 + t529) * pkin(5) + t453;
t509 = sin(qJ(6));
t512 = cos(qJ(6));
t445 = -t448 * t509 + t451 * t512;
t485 = qJD(5) * t512 - t509 * t532;
t474 = qJD(6) * t485 + qJDD(5) * t509 + t489 * t512;
t486 = qJD(5) * t509 + t512 * t532;
t475 = -mrSges(7,1) * t485 + mrSges(7,2) * t486;
t495 = qJD(6) - t531;
t476 = -mrSges(7,2) * t495 + mrSges(7,3) * t485;
t484 = qJDD(6) - t490;
t443 = m(7) * t445 + mrSges(7,1) * t484 - mrSges(7,3) * t474 - t475 * t486 + t476 * t495;
t446 = t448 * t512 + t451 * t509;
t473 = -qJD(6) * t486 + qJDD(5) * t512 - t489 * t509;
t477 = mrSges(7,1) * t495 - mrSges(7,3) * t486;
t444 = m(7) * t446 - mrSges(7,2) * t484 + mrSges(7,3) * t473 + t475 * t485 - t477 * t495;
t525 = -t443 * t509 + t444 * t512;
t436 = m(6) * t450 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t490 - qJD(5) * t493 + t487 * t531 + t525;
t534 = t513 * t458;
t449 = -t454 * t510 + t534;
t494 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t531;
t447 = -qJDD(5) * pkin(5) - t515 * pkin(11) - t534 + (qJD(4) * t488 + t454) * t510;
t518 = -m(7) * t447 + mrSges(7,1) * t473 - mrSges(7,2) * t474 + t476 * t485 - t477 * t486;
t441 = m(6) * t449 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t489 + qJD(5) * t494 - t487 * t532 + t518;
t526 = t436 * t513 - t441 * t510;
t427 = m(5) * t456 - mrSges(5,1) * t516 - qJDD(4) * mrSges(5,2) + t526;
t430 = t436 * t510 + t441 * t513;
t429 = m(5) * t458 + t430;
t437 = t443 * t512 + t444 * t509;
t517 = -m(6) * t453 + mrSges(6,1) * t490 - mrSges(6,2) * t489 - t493 * t532 + t494 * t531 - t437;
t433 = m(5) * t455 + qJDD(4) * mrSges(5,1) - mrSges(5,2) * t516 + t517;
t537 = t433 * t514;
t416 = t427 * t535 - t429 * t500 + t506 * t537;
t412 = m(4) * t462 + t416;
t421 = t427 * t514 - t433 * t511;
t420 = m(4) * t463 + t421;
t543 = t412 * t503 + t420 * t497;
t415 = t427 * t536 + t429 * t506 + t500 * t537;
t414 = m(4) * t465 + t415;
t401 = -t501 * t414 + t507 * t543;
t397 = m(3) * t466 + t401;
t406 = -t412 * t497 + t420 * t503;
t405 = m(3) * t467 + t406;
t542 = t397 * t504 + t405 * t498;
t400 = t507 * t414 + t501 * t543;
t399 = m(3) * t478 + t400;
t387 = -t502 * t399 + t508 * t542;
t385 = m(2) * t491 + t387;
t394 = -t397 * t498 + t405 * t504;
t393 = m(2) * t492 + t394;
t533 = t385 * t505 + t393 * t499;
t386 = t508 * t399 + t502 * t542;
t527 = -t385 * t499 + t393 * t505;
t468 = Ifges(7,5) * t486 + Ifges(7,6) * t485 + Ifges(7,3) * t495;
t470 = Ifges(7,1) * t486 + Ifges(7,4) * t485 + Ifges(7,5) * t495;
t438 = -mrSges(7,1) * t447 + mrSges(7,3) * t446 + Ifges(7,4) * t474 + Ifges(7,2) * t473 + Ifges(7,6) * t484 - t468 * t486 + t470 * t495;
t469 = Ifges(7,4) * t486 + Ifges(7,2) * t485 + Ifges(7,6) * t495;
t439 = mrSges(7,2) * t447 - mrSges(7,3) * t445 + Ifges(7,1) * t474 + Ifges(7,4) * t473 + Ifges(7,5) * t484 + t468 * t485 - t469 * t495;
t479 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t510 + Ifges(6,6) * t513) * qJD(4);
t480 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t510 + Ifges(6,2) * t513) * qJD(4);
t422 = mrSges(6,2) * t453 - mrSges(6,3) * t449 + Ifges(6,1) * t489 + Ifges(6,4) * t490 + Ifges(6,5) * qJDD(5) - pkin(11) * t437 - qJD(5) * t480 - t438 * t509 + t439 * t512 + t479 * t531;
t481 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t510 + Ifges(6,4) * t513) * qJD(4);
t423 = -mrSges(6,1) * t453 - mrSges(7,1) * t445 + mrSges(7,2) * t446 + mrSges(6,3) * t450 + Ifges(6,4) * t489 - Ifges(7,5) * t474 + Ifges(6,2) * t490 + Ifges(6,6) * qJDD(5) - Ifges(7,6) * t473 - Ifges(7,3) * t484 - pkin(5) * t437 + qJD(5) * t481 - t469 * t486 + t470 * t485 - t479 * t532;
t408 = mrSges(5,2) * t458 - mrSges(5,3) * t455 + Ifges(5,5) * qJDD(4) - Ifges(5,6) * t516 - pkin(10) * t430 + t422 * t513 - t423 * t510;
t409 = Ifges(5,6) * qJDD(4) + t516 * Ifges(5,5) - mrSges(5,1) * t458 + mrSges(5,3) * t456 - Ifges(6,5) * t489 - Ifges(6,6) * t490 - Ifges(6,3) * qJDD(5) - mrSges(6,1) * t449 + mrSges(6,2) * t450 - t509 * t439 - t512 * t438 - pkin(5) * t518 - pkin(11) * t525 - pkin(4) * t430 + (-t480 * t510 + t481 * t513) * qJD(4);
t521 = pkin(9) * t421 + t408 * t511 + t409 * t514;
t407 = mrSges(5,1) * t455 - mrSges(5,2) * t456 + Ifges(5,3) * qJDD(4) + pkin(4) * t517 + pkin(10) * t526 + t510 * t422 + t513 * t423;
t388 = mrSges(4,1) * t462 - mrSges(4,2) * t463 + pkin(3) * t416 + t506 * t407 + t500 * t521;
t389 = -mrSges(4,1) * t465 + mrSges(4,3) * t463 - pkin(3) * t415 - t500 * t407 + t506 * t521;
t390 = mrSges(4,2) * t465 - mrSges(4,3) * t462 + t514 * t408 - t511 * t409 + (-t415 * t500 - t416 * t506) * pkin(9);
t519 = qJ(3) * t406 + t389 * t503 + t390 * t497;
t382 = -mrSges(3,1) * t478 + mrSges(3,3) * t467 - pkin(2) * t400 - t501 * t388 + t507 * t519;
t383 = mrSges(3,2) * t478 - mrSges(3,3) * t466 - t497 * t389 + t503 * t390 + (-t400 * t501 - t401 * t507) * qJ(3);
t520 = qJ(2) * t394 + t382 * t504 + t383 * t498;
t381 = mrSges(3,1) * t466 - mrSges(3,2) * t467 + pkin(2) * t401 + t507 * t388 + t501 * t519;
t380 = mrSges(2,2) * t496 - mrSges(2,3) * t491 - t498 * t382 + t504 * t383 + (-t386 * t502 - t387 * t508) * qJ(2);
t379 = -mrSges(2,1) * t496 + mrSges(2,3) * t492 - pkin(1) * t386 - t502 * t381 + t508 * t520;
t1 = [-m(1) * g(1) + t527; -m(1) * g(2) + t533; -m(1) * g(3) + m(2) * t496 + t386; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t533 - t379 * t499 + t380 * t505; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t527 + t505 * t379 + t499 * t380; -mrSges(1,1) * g(2) + mrSges(2,1) * t491 + mrSges(1,2) * g(1) - mrSges(2,2) * t492 + pkin(1) * t387 + t508 * t381 + t502 * t520;];
tauB  = t1;
