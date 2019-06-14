% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRPPR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
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
% Datum: 2019-05-05 02:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRPPR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 02:31:03
% EndTime: 2019-05-05 02:31:06
% DurationCPUTime: 3.44s
% Computational Cost: add. (33526->289), mult. (74665->375), div. (0->0), fcn. (53894->14), ass. (0->119)
t516 = -2 * qJD(4);
t484 = sin(pkin(10));
t487 = cos(pkin(10));
t473 = g(1) * t484 - g(2) * t487;
t481 = -g(3) + qJDD(1);
t485 = sin(pkin(6));
t488 = cos(pkin(6));
t515 = t473 * t488 + t481 * t485;
t474 = -g(1) * t487 - g(2) * t484;
t491 = sin(qJ(2));
t494 = cos(qJ(2));
t437 = t494 * t474 + t515 * t491;
t496 = qJD(2) ^ 2;
t428 = -pkin(2) * t496 + qJDD(2) * pkin(8) + t437;
t454 = -t473 * t485 + t481 * t488;
t490 = sin(qJ(3));
t493 = cos(qJ(3));
t414 = -t428 * t490 + t493 * t454;
t508 = qJD(2) * qJD(3);
t507 = t493 * t508;
t471 = qJDD(2) * t490 + t507;
t409 = (-t471 + t507) * qJ(4) + (t490 * t493 * t496 + qJDD(3)) * pkin(3) + t414;
t415 = t493 * t428 + t490 * t454;
t472 = qJDD(2) * t493 - t490 * t508;
t511 = qJD(2) * t490;
t475 = qJD(3) * pkin(3) - qJ(4) * t511;
t480 = t493 ^ 2;
t410 = -pkin(3) * t480 * t496 + qJ(4) * t472 - qJD(3) * t475 + t415;
t483 = sin(pkin(11));
t514 = cos(pkin(11));
t459 = (t483 * t493 + t514 * t490) * qJD(2);
t393 = t514 * t409 - t483 * t410 + t459 * t516;
t436 = -t491 * t474 + t515 * t494;
t510 = qJD(2) * t493;
t458 = t483 * t511 - t514 * t510;
t394 = t483 * t409 + t514 * t410 + t458 * t516;
t440 = mrSges(5,1) * t458 + mrSges(5,2) * t459;
t443 = t471 * t483 - t514 * t472;
t453 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t459;
t439 = pkin(4) * t458 - qJ(5) * t459;
t495 = qJD(3) ^ 2;
t392 = -pkin(4) * t495 + qJDD(3) * qJ(5) - t439 * t458 + t394;
t499 = -qJDD(2) * pkin(2) - t436;
t411 = -pkin(3) * t472 + qJDD(4) + t475 * t511 + (-qJ(4) * t480 - pkin(8)) * t496 + t499;
t444 = t514 * t471 + t483 * t472;
t397 = (qJD(3) * t458 - t444) * qJ(5) + (qJD(3) * t459 + t443) * pkin(4) + t411;
t482 = sin(pkin(12));
t486 = cos(pkin(12));
t449 = qJD(3) * t482 + t459 * t486;
t387 = -0.2e1 * qJD(5) * t449 - t392 * t482 + t486 * t397;
t432 = qJDD(3) * t482 + t444 * t486;
t448 = qJD(3) * t486 - t459 * t482;
t385 = (t448 * t458 - t432) * pkin(9) + (t448 * t449 + t443) * pkin(5) + t387;
t388 = 0.2e1 * qJD(5) * t448 + t486 * t392 + t482 * t397;
t429 = pkin(5) * t458 - pkin(9) * t449;
t431 = qJDD(3) * t486 - t444 * t482;
t447 = t448 ^ 2;
t386 = -pkin(5) * t447 + pkin(9) * t431 - t429 * t458 + t388;
t489 = sin(qJ(6));
t492 = cos(qJ(6));
t383 = t385 * t492 - t386 * t489;
t420 = t448 * t492 - t449 * t489;
t400 = qJD(6) * t420 + t431 * t489 + t432 * t492;
t421 = t448 * t489 + t449 * t492;
t405 = -mrSges(7,1) * t420 + mrSges(7,2) * t421;
t457 = qJD(6) + t458;
t412 = -mrSges(7,2) * t457 + mrSges(7,3) * t420;
t442 = qJDD(6) + t443;
t380 = m(7) * t383 + mrSges(7,1) * t442 - mrSges(7,3) * t400 - t405 * t421 + t412 * t457;
t384 = t385 * t489 + t386 * t492;
t399 = -qJD(6) * t421 + t431 * t492 - t432 * t489;
t413 = mrSges(7,1) * t457 - mrSges(7,3) * t421;
t381 = m(7) * t384 - mrSges(7,2) * t442 + mrSges(7,3) * t399 + t405 * t420 - t413 * t457;
t372 = t492 * t380 + t489 * t381;
t422 = -mrSges(6,1) * t448 + mrSges(6,2) * t449;
t425 = -mrSges(6,2) * t458 + mrSges(6,3) * t448;
t370 = m(6) * t387 + mrSges(6,1) * t443 - mrSges(6,3) * t432 - t422 * t449 + t425 * t458 + t372;
t426 = mrSges(6,1) * t458 - mrSges(6,3) * t449;
t503 = -t380 * t489 + t492 * t381;
t371 = m(6) * t388 - mrSges(6,2) * t443 + mrSges(6,3) * t431 + t422 * t448 - t426 * t458 + t503;
t504 = -t370 * t482 + t486 * t371;
t364 = m(5) * t394 - qJDD(3) * mrSges(5,2) - mrSges(5,3) * t443 - qJD(3) * t453 - t440 * t458 + t504;
t391 = -qJDD(3) * pkin(4) - t495 * qJ(5) + t459 * t439 + qJDD(5) - t393;
t389 = -t431 * pkin(5) - t447 * pkin(9) + t449 * t429 + t391;
t500 = m(7) * t389 - t399 * mrSges(7,1) + mrSges(7,2) * t400 - t420 * t412 + t413 * t421;
t382 = m(6) * t391 - t431 * mrSges(6,1) + mrSges(6,2) * t432 - t448 * t425 + t426 * t449 + t500;
t452 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t458;
t376 = m(5) * t393 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t444 + qJD(3) * t452 - t440 * t459 - t382;
t359 = t483 * t364 + t514 * t376;
t366 = t486 * t370 + t482 * t371;
t470 = (-mrSges(4,1) * t493 + mrSges(4,2) * t490) * qJD(2);
t477 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t510;
t357 = m(4) * t414 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t471 + qJD(3) * t477 - t470 * t511 + t359;
t476 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t511;
t505 = t514 * t364 - t376 * t483;
t358 = m(4) * t415 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t472 - qJD(3) * t476 + t470 * t510 + t505;
t506 = -t357 * t490 + t493 * t358;
t365 = m(5) * t411 + t443 * mrSges(5,1) + mrSges(5,2) * t444 + t458 * t452 + t453 * t459 + t366;
t402 = Ifges(7,4) * t421 + Ifges(7,2) * t420 + Ifges(7,6) * t457;
t403 = Ifges(7,1) * t421 + Ifges(7,4) * t420 + Ifges(7,5) * t457;
t498 = mrSges(7,1) * t383 - mrSges(7,2) * t384 + Ifges(7,5) * t400 + Ifges(7,6) * t399 + Ifges(7,3) * t442 + t421 * t402 - t420 * t403;
t427 = -pkin(8) * t496 + t499;
t497 = -m(4) * t427 + t472 * mrSges(4,1) - mrSges(4,2) * t471 - t476 * t511 + t477 * t510 - t365;
t463 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t490 + Ifges(4,4) * t493) * qJD(2);
t462 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t490 + Ifges(4,2) * t493) * qJD(2);
t435 = Ifges(5,1) * t459 - Ifges(5,4) * t458 + Ifges(5,5) * qJD(3);
t434 = Ifges(5,4) * t459 - Ifges(5,2) * t458 + Ifges(5,6) * qJD(3);
t433 = Ifges(5,5) * t459 - Ifges(5,6) * t458 + Ifges(5,3) * qJD(3);
t418 = Ifges(6,1) * t449 + Ifges(6,4) * t448 + Ifges(6,5) * t458;
t417 = Ifges(6,4) * t449 + Ifges(6,2) * t448 + Ifges(6,6) * t458;
t416 = Ifges(6,5) * t449 + Ifges(6,6) * t448 + Ifges(6,3) * t458;
t401 = Ifges(7,5) * t421 + Ifges(7,6) * t420 + Ifges(7,3) * t457;
t374 = mrSges(7,2) * t389 - mrSges(7,3) * t383 + Ifges(7,1) * t400 + Ifges(7,4) * t399 + Ifges(7,5) * t442 + t401 * t420 - t402 * t457;
t373 = -mrSges(7,1) * t389 + mrSges(7,3) * t384 + Ifges(7,4) * t400 + Ifges(7,2) * t399 + Ifges(7,6) * t442 - t401 * t421 + t403 * t457;
t361 = mrSges(6,2) * t391 - mrSges(6,3) * t387 + Ifges(6,1) * t432 + Ifges(6,4) * t431 + Ifges(6,5) * t443 - pkin(9) * t372 - t373 * t489 + t374 * t492 + t416 * t448 - t417 * t458;
t360 = -mrSges(6,1) * t391 + mrSges(6,3) * t388 + Ifges(6,4) * t432 + Ifges(6,2) * t431 + Ifges(6,6) * t443 - pkin(5) * t500 + pkin(9) * t503 + t492 * t373 + t489 * t374 - t449 * t416 + t458 * t418;
t355 = Ifges(5,6) * qJDD(3) - t459 * t433 + t448 * t418 - t449 * t417 + qJD(3) * t435 + Ifges(5,4) * t444 - Ifges(6,6) * t431 - Ifges(6,5) * t432 - mrSges(5,1) * t411 + mrSges(5,3) * t394 - mrSges(6,1) * t387 + mrSges(6,2) * t388 - pkin(5) * t372 - pkin(4) * t366 - t498 + (-Ifges(5,2) - Ifges(6,3)) * t443;
t354 = mrSges(5,2) * t411 - mrSges(5,3) * t393 + Ifges(5,1) * t444 - Ifges(5,4) * t443 + Ifges(5,5) * qJDD(3) - qJ(5) * t366 - qJD(3) * t434 - t360 * t482 + t361 * t486 - t433 * t458;
t1 = [m(2) * t481 + t488 * (m(3) * t454 + t357 * t493 + t358 * t490) + (t491 * (m(3) * t437 - mrSges(3,1) * t496 - qJDD(2) * mrSges(3,2) + t506) + t494 * (m(3) * t436 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t496 + t497)) * t485; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t436 - mrSges(3,2) * t437 + t490 * (mrSges(4,2) * t427 - mrSges(4,3) * t414 + Ifges(4,1) * t471 + Ifges(4,4) * t472 + Ifges(4,5) * qJDD(3) - qJ(4) * t359 - qJD(3) * t462 + t514 * t354 - t483 * t355) + t493 * (-mrSges(4,1) * t427 + mrSges(4,3) * t415 + Ifges(4,4) * t471 + Ifges(4,2) * t472 + Ifges(4,6) * qJDD(3) - pkin(3) * t365 + qJ(4) * t505 + qJD(3) * t463 + t483 * t354 + t514 * t355) + pkin(2) * t497 + pkin(8) * t506; Ifges(4,5) * t471 + Ifges(4,6) * t472 + mrSges(4,1) * t414 - mrSges(4,2) * t415 + Ifges(5,5) * t444 - Ifges(5,6) * t443 + t459 * t434 + t458 * t435 + mrSges(5,1) * t393 - mrSges(5,2) * t394 + t482 * t361 + t486 * t360 - pkin(4) * t382 + qJ(5) * t504 + pkin(3) * t359 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (t462 * t490 - t463 * t493) * qJD(2); t365; t382; t498;];
tauJ  = t1;
