% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-05-05 03:41
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRPRP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:36:55
% EndTime: 2019-05-05 03:36:58
% DurationCPUTime: 2.08s
% Computational Cost: add. (16105->268), mult. (34550->333), div. (0->0), fcn. (24336->12), ass. (0->113)
t522 = -2 * qJD(4);
t521 = Ifges(6,1) + Ifges(7,1);
t515 = Ifges(6,4) + Ifges(7,4);
t514 = Ifges(6,5) + Ifges(7,5);
t520 = Ifges(6,2) + Ifges(7,2);
t513 = Ifges(6,6) + Ifges(7,6);
t519 = Ifges(6,3) + Ifges(7,3);
t475 = sin(pkin(10));
t478 = cos(pkin(10));
t466 = g(1) * t475 - g(2) * t478;
t473 = -g(3) + qJDD(1);
t476 = sin(pkin(6));
t479 = cos(pkin(6));
t518 = t466 * t479 + t473 * t476;
t467 = -g(1) * t478 - g(2) * t475;
t482 = sin(qJ(2));
t485 = cos(qJ(2));
t432 = t485 * t467 + t482 * t518;
t487 = qJD(2) ^ 2;
t427 = -pkin(2) * t487 + qJDD(2) * pkin(8) + t432;
t448 = -t466 * t476 + t473 * t479;
t481 = sin(qJ(3));
t484 = cos(qJ(3));
t401 = -t481 * t427 + t484 * t448;
t502 = qJD(2) * qJD(3);
t499 = t484 * t502;
t464 = qJDD(2) * t481 + t499;
t397 = (-t464 + t499) * qJ(4) + (t481 * t484 * t487 + qJDD(3)) * pkin(3) + t401;
t402 = t484 * t427 + t481 * t448;
t465 = qJDD(2) * t484 - t481 * t502;
t505 = qJD(2) * t481;
t468 = qJD(3) * pkin(3) - qJ(4) * t505;
t472 = t484 ^ 2;
t398 = -pkin(3) * t472 * t487 + qJ(4) * t465 - qJD(3) * t468 + t402;
t474 = sin(pkin(11));
t477 = cos(pkin(11));
t453 = (t474 * t484 + t477 * t481) * qJD(2);
t389 = t397 * t477 - t474 * t398 + t453 * t522;
t452 = (t474 * t481 - t477 * t484) * qJD(2);
t431 = -t482 * t467 + t485 * t518;
t440 = t464 * t477 + t465 * t474;
t480 = sin(qJ(5));
t483 = cos(qJ(5));
t442 = qJD(3) * t483 - t453 * t480;
t415 = qJD(5) * t442 + qJDD(3) * t480 + t440 * t483;
t443 = qJD(3) * t480 + t453 * t483;
t417 = -mrSges(7,1) * t442 + mrSges(7,2) * t443;
t390 = t474 * t397 + t477 * t398 + t452 * t522;
t435 = pkin(4) * t452 - pkin(9) * t453;
t486 = qJD(3) ^ 2;
t388 = -pkin(4) * t486 + qJDD(3) * pkin(9) - t435 * t452 + t390;
t491 = -qJDD(2) * pkin(2) - t431;
t399 = -t465 * pkin(3) + qJDD(4) + t468 * t505 + (-qJ(4) * t472 - pkin(8)) * t487 + t491;
t439 = -t464 * t474 + t465 * t477;
t393 = (qJD(3) * t452 - t440) * pkin(9) + (qJD(3) * t453 - t439) * pkin(4) + t399;
t383 = -t480 * t388 + t483 * t393;
t438 = qJDD(5) - t439;
t451 = qJD(5) + t452;
t380 = -0.2e1 * qJD(6) * t443 + (t442 * t451 - t415) * qJ(6) + (t442 * t443 + t438) * pkin(5) + t383;
t421 = -mrSges(7,2) * t451 + mrSges(7,3) * t442;
t501 = m(7) * t380 + t438 * mrSges(7,1) + t451 * t421;
t377 = -t415 * mrSges(7,3) - t443 * t417 + t501;
t384 = t483 * t388 + t480 * t393;
t414 = -qJD(5) * t443 + qJDD(3) * t483 - t440 * t480;
t423 = pkin(5) * t451 - qJ(6) * t443;
t441 = t442 ^ 2;
t382 = -pkin(5) * t441 + qJ(6) * t414 + 0.2e1 * qJD(6) * t442 - t423 * t451 + t384;
t507 = t442 * t515 + t443 * t521 + t451 * t514;
t508 = -t442 * t520 - t443 * t515 - t451 * t513;
t517 = mrSges(6,1) * t383 + mrSges(7,1) * t380 - mrSges(6,2) * t384 - mrSges(7,2) * t382 + pkin(5) * t377 + t513 * t414 + t514 * t415 + t438 * t519 - t507 * t442 - t508 * t443;
t516 = -mrSges(6,2) - mrSges(7,2);
t434 = mrSges(5,1) * t452 + mrSges(5,2) * t453;
t447 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t453;
t418 = -mrSges(6,1) * t442 + mrSges(6,2) * t443;
t422 = -mrSges(6,2) * t451 + mrSges(6,3) * t442;
t371 = m(6) * t383 + t438 * mrSges(6,1) + t451 * t422 + (-t417 - t418) * t443 + (-mrSges(6,3) - mrSges(7,3)) * t415 + t501;
t500 = m(7) * t382 + t414 * mrSges(7,3) + t442 * t417;
t424 = mrSges(7,1) * t451 - mrSges(7,3) * t443;
t506 = -mrSges(6,1) * t451 + mrSges(6,3) * t443 - t424;
t376 = m(6) * t384 + t414 * mrSges(6,3) + t442 * t418 + t516 * t438 + t506 * t451 + t500;
t496 = -t371 * t480 + t483 * t376;
t366 = m(5) * t390 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t439 - qJD(3) * t447 - t434 * t452 + t496;
t446 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t452;
t387 = -qJDD(3) * pkin(4) - pkin(9) * t486 + t453 * t435 - t389;
t385 = -pkin(5) * t414 - qJ(6) * t441 + t423 * t443 + qJDD(6) + t387;
t494 = -m(7) * t385 + t414 * mrSges(7,1) + t442 * t421;
t489 = -m(6) * t387 + t414 * mrSges(6,1) + t516 * t415 + t442 * t422 + t506 * t443 + t494;
t373 = m(5) * t389 + qJDD(3) * mrSges(5,1) - t440 * mrSges(5,3) + qJD(3) * t446 - t453 * t434 + t489;
t362 = t474 * t366 + t477 * t373;
t369 = t483 * t371 + t480 * t376;
t509 = -t442 * t513 - t443 * t514 - t451 * t519;
t504 = qJD(2) * t484;
t463 = (-mrSges(4,1) * t484 + mrSges(4,2) * t481) * qJD(2);
t470 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t504;
t360 = m(4) * t401 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t464 + qJD(3) * t470 - t463 * t505 + t362;
t469 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t505;
t497 = t477 * t366 - t373 * t474;
t361 = m(4) * t402 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t465 - qJD(3) * t469 + t463 * t504 + t497;
t498 = -t360 * t481 + t484 * t361;
t368 = m(5) * t399 - t439 * mrSges(5,1) + mrSges(5,2) * t440 + t452 * t446 + t447 * t453 + t369;
t426 = -t487 * pkin(8) + t491;
t488 = -m(4) * t426 + t465 * mrSges(4,1) - mrSges(4,2) * t464 - t469 * t505 + t470 * t504 - t368;
t457 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t481 + Ifges(4,4) * t484) * qJD(2);
t456 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t481 + Ifges(4,2) * t484) * qJD(2);
t430 = Ifges(5,1) * t453 - Ifges(5,4) * t452 + Ifges(5,5) * qJD(3);
t429 = Ifges(5,4) * t453 - Ifges(5,2) * t452 + Ifges(5,6) * qJD(3);
t428 = Ifges(5,5) * t453 - Ifges(5,6) * t452 + Ifges(5,3) * qJD(3);
t378 = t415 * mrSges(7,2) + t443 * t424 - t494;
t367 = mrSges(6,2) * t387 + mrSges(7,2) * t385 - mrSges(6,3) * t383 - mrSges(7,3) * t380 - qJ(6) * t377 + t515 * t414 + t415 * t521 + t514 * t438 - t509 * t442 + t508 * t451;
t363 = -mrSges(6,1) * t387 + mrSges(6,3) * t384 - mrSges(7,1) * t385 + mrSges(7,3) * t382 - pkin(5) * t378 + qJ(6) * t500 + (-qJ(6) * t424 + t507) * t451 + t509 * t443 + (-mrSges(7,2) * qJ(6) + t513) * t438 + t515 * t415 + t520 * t414;
t358 = -mrSges(5,1) * t399 + mrSges(5,3) * t390 + Ifges(5,4) * t440 + Ifges(5,2) * t439 + Ifges(5,6) * qJDD(3) - pkin(4) * t369 + qJD(3) * t430 - t453 * t428 - t517;
t357 = mrSges(5,2) * t399 - mrSges(5,3) * t389 + Ifges(5,1) * t440 + Ifges(5,4) * t439 + Ifges(5,5) * qJDD(3) - pkin(9) * t369 - qJD(3) * t429 - t363 * t480 + t367 * t483 - t428 * t452;
t1 = [m(2) * t473 + t479 * (m(3) * t448 + t360 * t484 + t361 * t481) + (t482 * (m(3) * t432 - mrSges(3,1) * t487 - qJDD(2) * mrSges(3,2) + t498) + t485 * (m(3) * t431 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t487 + t488)) * t476; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t431 - mrSges(3,2) * t432 + t481 * (mrSges(4,2) * t426 - mrSges(4,3) * t401 + Ifges(4,1) * t464 + Ifges(4,4) * t465 + Ifges(4,5) * qJDD(3) - qJ(4) * t362 - qJD(3) * t456 + t357 * t477 - t358 * t474) + t484 * (-mrSges(4,1) * t426 + mrSges(4,3) * t402 + Ifges(4,4) * t464 + Ifges(4,2) * t465 + Ifges(4,6) * qJDD(3) - pkin(3) * t368 + qJ(4) * t497 + qJD(3) * t457 + t474 * t357 + t477 * t358) + pkin(2) * t488 + pkin(8) * t498; Ifges(4,5) * t464 + Ifges(4,6) * t465 + mrSges(4,1) * t401 - mrSges(4,2) * t402 + Ifges(5,5) * t440 + Ifges(5,6) * t439 + t453 * t429 + t452 * t430 + mrSges(5,1) * t389 - mrSges(5,2) * t390 + t480 * t367 + t483 * t363 + pkin(4) * t489 + pkin(9) * t496 + pkin(3) * t362 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (t456 * t481 - t457 * t484) * qJD(2); t368; t517; t378;];
tauJ  = t1;
