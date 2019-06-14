% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPPRRP7
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-05-05 15:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPPRRP7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP7_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP7_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP7_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP7_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP7_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:03:32
% EndTime: 2019-05-05 15:03:34
% DurationCPUTime: 1.72s
% Computational Cost: add. (11948->236), mult. (26849->283), div. (0->0), fcn. (18261->8), ass. (0->106)
t513 = Ifges(6,1) + Ifges(7,1);
t504 = Ifges(6,4) + Ifges(7,4);
t503 = Ifges(6,5) + Ifges(7,5);
t512 = Ifges(6,2) + Ifges(7,2);
t502 = Ifges(6,6) + Ifges(7,6);
t511 = Ifges(6,3) + Ifges(7,3);
t467 = qJD(1) ^ 2;
t462 = sin(qJ(1));
t465 = cos(qJ(1));
t484 = g(1) * t462 - t465 * g(2);
t473 = -qJ(2) * t467 + qJDD(2) - t484;
t500 = -pkin(1) - qJ(3);
t510 = -(2 * qJD(1) * qJD(3)) + t500 * qJDD(1) + t473;
t459 = cos(pkin(9));
t509 = t459 ^ 2;
t479 = -g(1) * t465 - g(2) * t462;
t508 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t479;
t458 = sin(pkin(9));
t461 = sin(qJ(4));
t464 = cos(qJ(4));
t477 = t458 * t464 + t459 * t461;
t443 = t477 * qJD(1);
t476 = -t458 * t461 + t459 * t464;
t444 = t476 * qJD(1);
t490 = qJD(4) * t444;
t428 = -t477 * qJDD(1) - t490;
t491 = qJD(4) * t443;
t429 = t476 * qJDD(1) - t491;
t460 = sin(qJ(5));
t463 = cos(qJ(5));
t432 = qJD(4) * t463 - t444 * t460;
t400 = qJD(5) * t432 + qJDD(4) * t460 + t429 * t463;
t433 = qJD(4) * t460 + t444 * t463;
t403 = -mrSges(7,1) * t432 + mrSges(7,2) * t433;
t425 = t458 * g(3) + t510 * t459;
t506 = pkin(3) * t467;
t408 = (-pkin(7) * qJDD(1) - t458 * t506) * t459 + t425;
t426 = -g(3) * t459 + t510 * t458;
t455 = t458 ^ 2;
t488 = t458 * qJDD(1);
t409 = -pkin(7) * t488 - t455 * t506 + t426;
t386 = t461 * t408 + t464 * t409;
t427 = pkin(4) * t443 - pkin(8) * t444;
t466 = qJD(4) ^ 2;
t381 = -pkin(4) * t466 + qJDD(4) * pkin(8) - t427 * t443 + t386;
t472 = qJDD(3) + t508;
t493 = -t455 - t509;
t416 = pkin(3) * t488 + (t493 * pkin(7) + t500) * t467 + t472;
t384 = (-t429 + t491) * pkin(8) + (-t428 + t490) * pkin(4) + t416;
t377 = -t381 * t460 + t463 * t384;
t424 = qJDD(5) - t428;
t441 = qJD(5) + t443;
t373 = -0.2e1 * qJD(6) * t433 + (t432 * t441 - t400) * qJ(6) + (t432 * t433 + t424) * pkin(5) + t377;
t410 = -mrSges(7,2) * t441 + mrSges(7,3) * t432;
t486 = m(7) * t373 + t424 * mrSges(7,1) + t441 * t410;
t370 = -mrSges(7,3) * t400 - t403 * t433 + t486;
t378 = t463 * t381 + t460 * t384;
t399 = -qJD(5) * t433 + qJDD(4) * t463 - t429 * t460;
t412 = pkin(5) * t441 - qJ(6) * t433;
t431 = t432 ^ 2;
t375 = -pkin(5) * t431 + qJ(6) * t399 + 0.2e1 * qJD(6) * t432 - t412 * t441 + t378;
t495 = t504 * t432 + t513 * t433 + t503 * t441;
t496 = -t512 * t432 - t504 * t433 - t502 * t441;
t507 = mrSges(6,1) * t377 + mrSges(7,1) * t373 - mrSges(6,2) * t378 - mrSges(7,2) * t375 + pkin(5) * t370 + t502 * t399 + t503 * t400 + t511 * t424 - t495 * t432 - t496 * t433;
t505 = -mrSges(6,2) - mrSges(7,2);
t499 = t459 * mrSges(4,2);
t422 = mrSges(5,1) * t443 + mrSges(5,2) * t444;
t438 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t444;
t404 = -mrSges(6,1) * t432 + mrSges(6,2) * t433;
t411 = -mrSges(6,2) * t441 + mrSges(6,3) * t432;
t364 = m(6) * t377 + mrSges(6,1) * t424 + t411 * t441 + (-t403 - t404) * t433 + (-mrSges(6,3) - mrSges(7,3)) * t400 + t486;
t485 = m(7) * t375 + t399 * mrSges(7,3) + t432 * t403;
t413 = mrSges(7,1) * t441 - mrSges(7,3) * t433;
t494 = -mrSges(6,1) * t441 + mrSges(6,3) * t433 - t413;
t367 = m(6) * t378 + mrSges(6,3) * t399 + t404 * t432 + t505 * t424 + t494 * t441 + t485;
t481 = -t364 * t460 + t463 * t367;
t360 = m(5) * t386 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t428 - qJD(4) * t438 - t422 * t443 + t481;
t385 = t408 * t464 - t461 * t409;
t437 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t443;
t380 = -qJDD(4) * pkin(4) - pkin(8) * t466 + t444 * t427 - t385;
t376 = -pkin(5) * t399 - qJ(6) * t431 + t412 * t433 + qJDD(6) + t380;
t480 = -m(7) * t376 + t399 * mrSges(7,1) + t432 * t410;
t468 = -m(6) * t380 + t399 * mrSges(6,1) + t505 * t400 + t432 * t411 + t494 * t433 + t480;
t369 = m(5) * t385 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t429 + qJD(4) * t437 - t422 * t444 + t468;
t498 = t461 * t360 + t464 * t369;
t362 = t463 * t364 + t460 * t367;
t497 = -t502 * t432 - t503 * t433 - t511 * t441;
t483 = t493 * mrSges(4,3);
t482 = t464 * t360 - t369 * t461;
t475 = -qJDD(1) * mrSges(4,3) - t467 * (t458 * mrSges(4,1) + t499);
t478 = (m(4) * t425 + t475 * t459 + t498) * t459 + (m(4) * t426 + t475 * t458 + t482) * t458;
t471 = m(5) * t416 - mrSges(5,1) * t428 + t429 * mrSges(5,2) + t437 * t443 + t444 * t438 + t362;
t436 = t500 * t467 + t472;
t469 = m(4) * t436 + mrSges(4,1) * t488 + qJDD(1) * t499 + t471;
t442 = -qJDD(1) * pkin(1) + t473;
t440 = pkin(1) * t467 - t508;
t419 = Ifges(5,1) * t444 - Ifges(5,4) * t443 + Ifges(5,5) * qJD(4);
t418 = Ifges(5,4) * t444 - Ifges(5,2) * t443 + Ifges(5,6) * qJD(4);
t417 = Ifges(5,5) * t444 - Ifges(5,6) * t443 + Ifges(5,3) * qJD(4);
t371 = mrSges(7,2) * t400 + t413 * t433 - t480;
t361 = mrSges(6,2) * t380 + mrSges(7,2) * t376 - mrSges(6,3) * t377 - mrSges(7,3) * t373 - qJ(6) * t370 + t504 * t399 + t513 * t400 + t503 * t424 - t497 * t432 + t496 * t441;
t357 = -mrSges(6,1) * t380 + mrSges(6,3) * t378 - mrSges(7,1) * t376 + mrSges(7,3) * t375 - pkin(5) * t371 + qJ(6) * t485 + (-qJ(6) * t413 + t495) * t441 + t497 * t433 + (-qJ(6) * mrSges(7,2) + t502) * t424 + t504 * t400 + t512 * t399;
t354 = -mrSges(5,1) * t416 + mrSges(5,3) * t386 + Ifges(5,4) * t429 + Ifges(5,2) * t428 + Ifges(5,6) * qJDD(4) - pkin(4) * t362 + qJD(4) * t419 - t444 * t417 - t507;
t353 = m(3) * t442 + qJDD(1) * mrSges(3,2) - mrSges(3,3) * t467 + t478;
t352 = mrSges(5,2) * t416 - mrSges(5,3) * t385 + Ifges(5,1) * t429 + Ifges(5,4) * t428 + Ifges(5,5) * qJDD(4) - pkin(8) * t362 - qJD(4) * t418 - t357 * t460 + t361 * t463 - t417 * t443;
t1 = [mrSges(2,1) * t484 - mrSges(2,2) * t479 + mrSges(3,2) * t442 - mrSges(3,3) * t440 + t459 * (mrSges(4,2) * t436 - mrSges(4,3) * t425 - pkin(7) * t498 + t464 * t352 - t461 * t354) - t458 * (-mrSges(4,1) * t436 + mrSges(4,3) * t426 - pkin(3) * t471 + pkin(7) * t482 + t461 * t352 + t464 * t354) - qJ(3) * t478 - pkin(1) * t353 + qJ(2) * (-m(3) * t440 + (mrSges(3,2) + t483) * t467 + t469) + (Ifges(4,1) * t509 + qJ(2) * mrSges(3,3) + Ifges(3,1) + Ifges(2,3) + (-0.2e1 * Ifges(4,4) * t459 + Ifges(4,2) * t458) * t458) * qJDD(1); t353; t467 * t483 + t469; mrSges(5,1) * t385 - mrSges(5,2) * t386 + Ifges(5,5) * t429 + Ifges(5,6) * t428 + Ifges(5,3) * qJDD(4) + pkin(4) * t468 + pkin(8) * t481 + t463 * t357 + t460 * t361 + t444 * t418 + t443 * t419; t507; t371;];
tauJ  = t1;
