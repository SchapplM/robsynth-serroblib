% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-05-05 14:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPPRRP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:43:52
% EndTime: 2019-05-05 14:43:55
% DurationCPUTime: 1.80s
% Computational Cost: add. (13538->244), mult. (30030->295), div. (0->0), fcn. (20624->10), ass. (0->111)
t501 = Ifges(6,1) + Ifges(7,1);
t495 = Ifges(6,4) + Ifges(7,4);
t494 = Ifges(6,5) + Ifges(7,5);
t500 = Ifges(6,2) + Ifges(7,2);
t493 = Ifges(6,6) + Ifges(7,6);
t499 = Ifges(6,3) + Ifges(7,3);
t462 = qJD(1) ^ 2;
t451 = sin(pkin(10));
t453 = cos(pkin(10));
t456 = sin(qJ(4));
t459 = cos(qJ(4));
t468 = t451 * t456 - t453 * t459;
t431 = t468 * qJD(1);
t469 = t451 * t459 + t453 * t456;
t432 = t469 * qJD(1);
t481 = t432 * qJD(4);
t421 = -qJDD(1) * t468 - t481;
t482 = t431 * qJD(4);
t422 = qJDD(1) * t469 - t482;
t455 = sin(qJ(5));
t458 = cos(qJ(5));
t426 = t458 * qJD(4) - t455 * t432;
t395 = qJD(5) * t426 + qJDD(4) * t455 + t422 * t458;
t427 = t455 * qJD(4) + t458 * t432;
t401 = -mrSges(7,1) * t426 + mrSges(7,2) * t427;
t457 = sin(qJ(1));
t460 = cos(qJ(1));
t476 = t457 * g(1) - t460 * g(2);
t438 = qJDD(1) * pkin(1) + t476;
t471 = -t460 * g(1) - t457 * g(2);
t439 = -t462 * pkin(1) + t471;
t452 = sin(pkin(9));
t454 = cos(pkin(9));
t424 = t452 * t438 + t454 * t439;
t414 = -pkin(2) * t462 + qJDD(1) * qJ(3) + t424;
t450 = -g(3) + qJDD(2);
t480 = qJD(1) * qJD(3);
t484 = t453 * t450 - 0.2e1 * t451 * t480;
t497 = pkin(3) * t453;
t393 = (-pkin(7) * qJDD(1) + t462 * t497 - t414) * t451 + t484;
t400 = t451 * t450 + (t414 + 0.2e1 * t480) * t453;
t479 = qJDD(1) * t453;
t449 = t453 ^ 2;
t490 = t449 * t462;
t396 = -pkin(3) * t490 + pkin(7) * t479 + t400;
t377 = t456 * t393 + t459 * t396;
t420 = t431 * pkin(4) - t432 * pkin(8);
t461 = qJD(4) ^ 2;
t372 = -pkin(4) * t461 + qJDD(4) * pkin(8) - t420 * t431 + t377;
t448 = t451 ^ 2;
t423 = t454 * t438 - t452 * t439;
t470 = qJDD(3) - t423;
t398 = (-pkin(2) - t497) * qJDD(1) + (-qJ(3) + (-t448 - t449) * pkin(7)) * t462 + t470;
t375 = (-t422 + t482) * pkin(8) + (-t421 + t481) * pkin(4) + t398;
t367 = -t455 * t372 + t458 * t375;
t419 = qJDD(5) - t421;
t430 = qJD(5) + t431;
t364 = -0.2e1 * qJD(6) * t427 + (t426 * t430 - t395) * qJ(6) + (t426 * t427 + t419) * pkin(5) + t367;
t404 = -mrSges(7,2) * t430 + mrSges(7,3) * t426;
t478 = m(7) * t364 + t419 * mrSges(7,1) + t430 * t404;
t361 = -t395 * mrSges(7,3) - t427 * t401 + t478;
t368 = t458 * t372 + t455 * t375;
t394 = -qJD(5) * t427 + qJDD(4) * t458 - t422 * t455;
t406 = pkin(5) * t430 - qJ(6) * t427;
t425 = t426 ^ 2;
t366 = -pkin(5) * t425 + qJ(6) * t394 + 0.2e1 * qJD(6) * t426 - t406 * t430 + t368;
t486 = t495 * t426 + t501 * t427 + t494 * t430;
t487 = -t500 * t426 - t495 * t427 - t493 * t430;
t498 = mrSges(6,1) * t367 + mrSges(7,1) * t364 - mrSges(6,2) * t368 - mrSges(7,2) * t366 + pkin(5) * t361 + t394 * t493 + t395 * t494 + t499 * t419 - t426 * t486 - t427 * t487;
t496 = -mrSges(6,2) - mrSges(7,2);
t491 = mrSges(4,2) * t451;
t417 = t431 * mrSges(5,1) + t432 * mrSges(5,2);
t429 = qJD(4) * mrSges(5,1) - t432 * mrSges(5,3);
t402 = -mrSges(6,1) * t426 + mrSges(6,2) * t427;
t405 = -mrSges(6,2) * t430 + mrSges(6,3) * t426;
t355 = m(6) * t367 + t419 * mrSges(6,1) + t430 * t405 + (-t401 - t402) * t427 + (-mrSges(6,3) - mrSges(7,3)) * t395 + t478;
t477 = m(7) * t366 + t394 * mrSges(7,3) + t426 * t401;
t407 = mrSges(7,1) * t430 - mrSges(7,3) * t427;
t485 = -mrSges(6,1) * t430 + mrSges(6,3) * t427 - t407;
t358 = m(6) * t368 + t394 * mrSges(6,3) + t426 * t402 + t419 * t496 + t430 * t485 + t477;
t473 = -t355 * t455 + t458 * t358;
t351 = m(5) * t377 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t421 - qJD(4) * t429 - t417 * t431 + t473;
t376 = t393 * t459 - t456 * t396;
t428 = -qJD(4) * mrSges(5,2) - t431 * mrSges(5,3);
t371 = -qJDD(4) * pkin(4) - pkin(8) * t461 + t432 * t420 - t376;
t369 = -pkin(5) * t394 - qJ(6) * t425 + t406 * t427 + qJDD(6) + t371;
t472 = -m(7) * t369 + t394 * mrSges(7,1) + t426 * t404;
t464 = -m(6) * t371 + t394 * mrSges(6,1) + t395 * t496 + t426 * t405 + t427 * t485 + t472;
t360 = m(5) * t376 + qJDD(4) * mrSges(5,1) - t422 * mrSges(5,3) + qJD(4) * t428 - t432 * t417 + t464;
t489 = t456 * t351 + t459 * t360;
t353 = t458 * t355 + t455 * t358;
t488 = -t493 * t426 - t494 * t427 - t499 * t430;
t399 = -t414 * t451 + t484;
t467 = mrSges(4,3) * qJDD(1) + t462 * (-t453 * mrSges(4,1) + t491);
t345 = m(4) * t399 - t451 * t467 + t489;
t474 = t459 * t351 - t456 * t360;
t346 = m(4) * t400 + t453 * t467 + t474;
t475 = -t345 * t451 + t453 * t346;
t466 = m(5) * t398 - t421 * mrSges(5,1) + t422 * mrSges(5,2) + t431 * t428 + t432 * t429 + t353;
t410 = -qJDD(1) * pkin(2) - t462 * qJ(3) + t470;
t463 = -m(4) * t410 + mrSges(4,1) * t479 - t466 + (t448 * t462 + t490) * mrSges(4,3);
t413 = Ifges(5,1) * t432 - Ifges(5,4) * t431 + Ifges(5,5) * qJD(4);
t412 = Ifges(5,4) * t432 - Ifges(5,2) * t431 + Ifges(5,6) * qJD(4);
t411 = Ifges(5,5) * t432 - Ifges(5,6) * t431 + Ifges(5,3) * qJD(4);
t362 = t395 * mrSges(7,2) + t427 * t407 - t472;
t352 = mrSges(6,2) * t371 + mrSges(7,2) * t369 - mrSges(6,3) * t367 - mrSges(7,3) * t364 - qJ(6) * t361 + t495 * t394 + t501 * t395 + t494 * t419 - t488 * t426 + t487 * t430;
t348 = qJDD(1) * t491 - t463;
t347 = -mrSges(6,1) * t371 + mrSges(6,3) * t368 - mrSges(7,1) * t369 + mrSges(7,3) * t366 - pkin(5) * t362 + qJ(6) * t477 + (-qJ(6) * t407 + t486) * t430 + t488 * t427 + (-mrSges(7,2) * qJ(6) + t493) * t419 + t495 * t395 + t500 * t394;
t343 = -mrSges(5,1) * t398 + mrSges(5,3) * t377 + Ifges(5,4) * t422 + Ifges(5,2) * t421 + Ifges(5,6) * qJDD(4) - pkin(4) * t353 + qJD(4) * t413 - t432 * t411 - t498;
t342 = mrSges(5,2) * t398 - mrSges(5,3) * t376 + Ifges(5,1) * t422 + Ifges(5,4) * t421 + Ifges(5,5) * qJDD(4) - pkin(8) * t353 - qJD(4) * t412 - t347 * t455 + t352 * t458 - t411 * t431;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t476 - mrSges(2,2) * t471 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t423 - mrSges(3,2) * t424 + t451 * (mrSges(4,2) * t410 - mrSges(4,3) * t399 + t459 * t342 - t456 * t343 - pkin(7) * t489 + (Ifges(4,1) * t451 + Ifges(4,4) * t453) * qJDD(1)) + t453 * (-mrSges(4,1) * t410 + mrSges(4,3) * t400 + t456 * t342 + t459 * t343 - pkin(3) * t466 + pkin(7) * t474 + (Ifges(4,4) * t451 + Ifges(4,2) * t453) * qJDD(1)) - pkin(2) * t348 + qJ(3) * t475 + pkin(1) * (t452 * (m(3) * t424 - mrSges(3,1) * t462 - qJDD(1) * mrSges(3,2) + t475) + t454 * (t463 - t462 * mrSges(3,2) + m(3) * t423 + (mrSges(3,1) - t491) * qJDD(1))); m(3) * t450 + t345 * t453 + t346 * t451; t348; mrSges(5,1) * t376 - mrSges(5,2) * t377 + Ifges(5,5) * t422 + Ifges(5,6) * t421 + Ifges(5,3) * qJDD(4) + pkin(4) * t464 + pkin(8) * t473 + t458 * t347 + t455 * t352 + t432 * t412 + t431 * t413; t498; t362;];
tauJ  = t1;
