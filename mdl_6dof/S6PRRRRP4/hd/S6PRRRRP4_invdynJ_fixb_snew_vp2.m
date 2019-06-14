% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-05-05 09:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRRRP4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:51:28
% EndTime: 2019-05-05 09:51:32
% DurationCPUTime: 2.62s
% Computational Cost: add. (21377->267), mult. (41189->328), div. (0->0), fcn. (29204->12), ass. (0->114)
t510 = Ifges(6,1) + Ifges(7,1);
t500 = Ifges(6,4) - Ifges(7,5);
t508 = Ifges(7,4) + Ifges(6,5);
t509 = Ifges(6,2) + Ifges(7,3);
t507 = Ifges(6,6) - Ifges(7,6);
t506 = -Ifges(6,3) - Ifges(7,2);
t473 = sin(qJ(4));
t476 = cos(qJ(4));
t474 = sin(qJ(3));
t495 = qJD(2) * t474;
t452 = qJD(3) * t476 - t473 * t495;
t453 = qJD(3) * t473 + t476 * t495;
t472 = sin(qJ(5));
t502 = cos(qJ(5));
t427 = -t452 * t502 + t472 * t453;
t428 = t472 * t452 + t453 * t502;
t477 = cos(qJ(3));
t494 = qJD(2) * t477;
t465 = qJD(4) - t494;
t464 = qJD(5) + t465;
t505 = t427 * t509 - t428 * t500 - t464 * t507;
t504 = -t500 * t427 + t428 * t510 + t508 * t464;
t468 = sin(pkin(11));
t470 = cos(pkin(11));
t458 = g(1) * t468 - g(2) * t470;
t467 = -g(3) + qJDD(1);
t469 = sin(pkin(6));
t471 = cos(pkin(6));
t503 = t458 * t471 + t467 * t469;
t459 = -g(1) * t470 - g(2) * t468;
t475 = sin(qJ(2));
t478 = cos(qJ(2));
t416 = -t475 * t459 + t478 * t503;
t501 = -mrSges(6,3) - mrSges(7,2);
t417 = t478 * t459 + t503 * t475;
t480 = qJD(2) ^ 2;
t411 = -pkin(2) * t480 + qJDD(2) * pkin(8) + t417;
t434 = -t458 * t469 + t467 * t471;
t402 = t477 * t411 + t474 * t434;
t455 = (-pkin(3) * t477 - pkin(9) * t474) * qJD(2);
t479 = qJD(3) ^ 2;
t388 = -pkin(3) * t479 + qJDD(3) * pkin(9) + t455 * t494 + t402;
t410 = -qJDD(2) * pkin(2) - t480 * pkin(8) - t416;
t493 = qJD(2) * qJD(3);
t491 = t477 * t493;
t456 = qJDD(2) * t474 + t491;
t466 = t474 * t493;
t457 = qJDD(2) * t477 - t466;
t391 = (-t456 - t491) * pkin(9) + (-t457 + t466) * pkin(3) + t410;
t372 = -t473 * t388 + t476 * t391;
t426 = qJD(4) * t452 + qJDD(3) * t473 + t456 * t476;
t449 = qJDD(4) - t457;
t369 = (t452 * t465 - t426) * pkin(10) + (t452 * t453 + t449) * pkin(4) + t372;
t373 = t476 * t388 + t473 * t391;
t425 = -qJD(4) * t453 + qJDD(3) * t476 - t456 * t473;
t433 = pkin(4) * t465 - pkin(10) * t453;
t448 = t452 ^ 2;
t371 = -pkin(4) * t448 + pkin(10) * t425 - t433 * t465 + t373;
t365 = t472 * t369 + t371 * t502;
t384 = t428 * qJD(5) - t425 * t502 + t472 * t426;
t414 = mrSges(6,1) * t464 - mrSges(6,3) * t428;
t445 = qJDD(5) + t449;
t403 = pkin(5) * t427 - qJ(6) * t428;
t463 = t464 ^ 2;
t361 = -pkin(5) * t463 + qJ(6) * t445 + 0.2e1 * qJD(6) * t464 - t403 * t427 + t365;
t415 = -mrSges(7,1) * t464 + mrSges(7,2) * t428;
t492 = m(7) * t361 + t445 * mrSges(7,3) + t464 * t415;
t404 = mrSges(7,1) * t427 - mrSges(7,3) * t428;
t496 = -mrSges(6,1) * t427 - mrSges(6,2) * t428 - t404;
t351 = m(6) * t365 - t445 * mrSges(6,2) + t384 * t501 - t464 * t414 + t427 * t496 + t492;
t364 = t369 * t502 - t472 * t371;
t385 = -t427 * qJD(5) + t472 * t425 + t426 * t502;
t413 = -mrSges(6,2) * t464 - mrSges(6,3) * t427;
t362 = -t445 * pkin(5) - t463 * qJ(6) + t428 * t403 + qJDD(6) - t364;
t412 = -mrSges(7,2) * t427 + mrSges(7,3) * t464;
t487 = -m(7) * t362 + t445 * mrSges(7,1) + t464 * t412;
t353 = m(6) * t364 + t445 * mrSges(6,1) + t385 * t501 + t464 * t413 + t428 * t496 + t487;
t347 = t472 * t351 + t353 * t502;
t497 = t507 * t427 - t508 * t428 + t506 * t464;
t454 = (-mrSges(4,1) * t477 + mrSges(4,2) * t474) * qJD(2);
t460 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t495;
t429 = -mrSges(5,1) * t452 + mrSges(5,2) * t453;
t431 = -mrSges(5,2) * t465 + mrSges(5,3) * t452;
t343 = m(5) * t372 + mrSges(5,1) * t449 - mrSges(5,3) * t426 - t429 * t453 + t431 * t465 + t347;
t432 = mrSges(5,1) * t465 - mrSges(5,3) * t453;
t488 = t351 * t502 - t353 * t472;
t344 = m(5) * t373 - mrSges(5,2) * t449 + mrSges(5,3) * t425 + t429 * t452 - t432 * t465 + t488;
t489 = -t343 * t473 + t476 * t344;
t340 = m(4) * t402 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t457 - qJD(3) * t460 + t454 * t494 + t489;
t401 = -t474 * t411 + t477 * t434;
t461 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t494;
t387 = -qJDD(3) * pkin(3) - t479 * pkin(9) + t455 * t495 - t401;
t374 = -t425 * pkin(4) - t448 * pkin(10) + t453 * t433 + t387;
t367 = -0.2e1 * qJD(6) * t428 + (t427 * t464 - t385) * qJ(6) + (t428 * t464 + t384) * pkin(5) + t374;
t358 = m(7) * t367 + t384 * mrSges(7,1) - t385 * mrSges(7,3) + t427 * t412 - t428 * t415;
t485 = m(6) * t374 + t384 * mrSges(6,1) + t385 * mrSges(6,2) + t427 * t413 + t428 * t414 + t358;
t482 = -m(5) * t387 + t425 * mrSges(5,1) - t426 * mrSges(5,2) + t452 * t431 - t453 * t432 - t485;
t348 = m(4) * t401 + qJDD(3) * mrSges(4,1) - t456 * mrSges(4,3) + qJD(3) * t461 - t454 * t495 + t482;
t490 = t477 * t340 - t348 * t474;
t341 = t343 * t476 + t344 * t473;
t484 = -m(4) * t410 + t457 * mrSges(4,1) - mrSges(4,2) * t456 - t460 * t495 + t461 * t494 - t341;
t357 = t385 * mrSges(7,2) + t428 * t404 - t487;
t483 = -mrSges(6,1) * t364 + mrSges(7,1) * t362 + mrSges(6,2) * t365 - mrSges(7,3) * t361 + pkin(5) * t357 - qJ(6) * t492 + t506 * t445 + t505 * t428 + (qJ(6) * t404 - t504) * t427 - t508 * t385 + (qJ(6) * mrSges(7,2) + t507) * t384;
t420 = Ifges(5,4) * t453 + Ifges(5,2) * t452 + Ifges(5,6) * t465;
t421 = Ifges(5,1) * t453 + Ifges(5,4) * t452 + Ifges(5,5) * t465;
t481 = mrSges(5,1) * t372 - mrSges(5,2) * t373 + Ifges(5,5) * t426 + Ifges(5,6) * t425 + Ifges(5,3) * t449 + pkin(4) * t347 + t453 * t420 - t452 * t421 - t483;
t444 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t474 + Ifges(4,4) * t477) * qJD(2);
t443 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t474 + Ifges(4,2) * t477) * qJD(2);
t419 = Ifges(5,5) * t453 + Ifges(5,6) * t452 + Ifges(5,3) * t465;
t346 = mrSges(6,2) * t374 + mrSges(7,2) * t362 - mrSges(6,3) * t364 - mrSges(7,3) * t367 - qJ(6) * t358 - t500 * t384 + t385 * t510 + t497 * t427 + t508 * t445 + t505 * t464;
t345 = -mrSges(6,1) * t374 - mrSges(7,1) * t367 + mrSges(7,2) * t361 + mrSges(6,3) * t365 - pkin(5) * t358 - t384 * t509 + t500 * t385 + t497 * t428 + t507 * t445 + t504 * t464;
t338 = mrSges(5,2) * t387 - mrSges(5,3) * t372 + Ifges(5,1) * t426 + Ifges(5,4) * t425 + Ifges(5,5) * t449 - pkin(10) * t347 - t472 * t345 + t346 * t502 + t452 * t419 - t465 * t420;
t337 = -mrSges(5,1) * t387 + mrSges(5,3) * t373 + Ifges(5,4) * t426 + Ifges(5,2) * t425 + Ifges(5,6) * t449 - pkin(4) * t485 + pkin(10) * t488 + t345 * t502 + t472 * t346 - t453 * t419 + t465 * t421;
t1 = [m(2) * t467 + t471 * (m(3) * t434 + t340 * t474 + t348 * t477) + (t475 * (m(3) * t417 - mrSges(3,1) * t480 - qJDD(2) * mrSges(3,2) + t490) + t478 * (m(3) * t416 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t480 + t484)) * t469; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t416 - mrSges(3,2) * t417 + t474 * (mrSges(4,2) * t410 - mrSges(4,3) * t401 + Ifges(4,1) * t456 + Ifges(4,4) * t457 + Ifges(4,5) * qJDD(3) - pkin(9) * t341 - qJD(3) * t443 - t337 * t473 + t338 * t476) + t477 * (-mrSges(4,1) * t410 + mrSges(4,3) * t402 + Ifges(4,4) * t456 + Ifges(4,2) * t457 + Ifges(4,6) * qJDD(3) - pkin(3) * t341 + qJD(3) * t444 - t481) + pkin(2) * t484 + pkin(8) * t490; Ifges(4,5) * t456 + Ifges(4,6) * t457 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t401 - mrSges(4,2) * t402 + t473 * t338 + t476 * t337 + pkin(3) * t482 + pkin(9) * t489 + (t443 * t474 - t444 * t477) * qJD(2); t481; -t483; t357;];
tauJ  = t1;
