% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRPRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-05-05 17:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRPRP7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP7_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP7_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP7_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP7_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP7_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:57:02
% EndTime: 2019-05-05 17:57:06
% DurationCPUTime: 1.92s
% Computational Cost: add. (13739->267), mult. (29844->324), div. (0->0), fcn. (19147->8), ass. (0->106)
t507 = -2 * qJD(4);
t506 = Ifges(6,1) + Ifges(7,1);
t500 = Ifges(6,4) + Ifges(7,4);
t499 = Ifges(6,5) + Ifges(7,5);
t505 = Ifges(6,2) + Ifges(7,2);
t498 = Ifges(6,6) + Ifges(7,6);
t504 = Ifges(6,3) + Ifges(7,3);
t471 = qJD(1) ^ 2;
t466 = sin(qJ(1));
t469 = cos(qJ(1));
t484 = t466 * g(1) - t469 * g(2);
t474 = -qJ(2) * t471 + qJDD(2) - t484;
t502 = -pkin(1) - pkin(7);
t438 = t502 * qJDD(1) + t474;
t465 = sin(qJ(3));
t468 = cos(qJ(3));
t428 = t465 * g(3) + t468 * t438;
t489 = qJD(1) * qJD(3);
t485 = t465 * t489;
t452 = qJDD(1) * t468 - t485;
t405 = (-t452 - t485) * qJ(4) + (-t465 * t468 * t471 + qJDD(3)) * pkin(3) + t428;
t429 = -g(3) * t468 + t465 * t438;
t451 = -qJDD(1) * t465 - t468 * t489;
t491 = qJD(1) * t468;
t454 = qJD(3) * pkin(3) - qJ(4) * t491;
t461 = t465 ^ 2;
t406 = -pkin(3) * t461 * t471 + qJ(4) * t451 - qJD(3) * t454 + t429;
t462 = sin(pkin(9));
t463 = cos(pkin(9));
t492 = qJD(1) * t465;
t445 = -t462 * t492 + t463 * t491;
t385 = t405 * t463 - t462 * t406 + t445 * t507;
t479 = -g(1) * t469 - g(2) * t466;
t475 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t479;
t444 = (t468 * t462 + t465 * t463) * qJD(1);
t427 = t451 * t462 + t452 * t463;
t464 = sin(qJ(5));
t467 = cos(qJ(5));
t431 = qJD(3) * t467 - t445 * t464;
t400 = qJD(5) * t431 + qJDD(3) * t464 + t427 * t467;
t432 = qJD(3) * t464 + t445 * t467;
t409 = -mrSges(7,1) * t431 + mrSges(7,2) * t432;
t386 = t462 * t405 + t463 * t406 + t444 * t507;
t422 = pkin(4) * t444 - pkin(8) * t445;
t470 = qJD(3) ^ 2;
t381 = -pkin(4) * t470 + qJDD(3) * pkin(8) - t422 * t444 + t386;
t408 = -pkin(3) * t451 + qJDD(4) + t454 * t491 + (-qJ(4) * t461 + t502) * t471 + t475;
t426 = t451 * t463 - t452 * t462;
t384 = (qJD(3) * t444 - t427) * pkin(8) + (qJD(3) * t445 - t426) * pkin(4) + t408;
t377 = -t381 * t464 + t467 * t384;
t425 = qJDD(5) - t426;
t442 = qJD(5) + t444;
t373 = -0.2e1 * qJD(6) * t432 + (t431 * t442 - t400) * qJ(6) + (t431 * t432 + t425) * pkin(5) + t377;
t412 = -mrSges(7,2) * t442 + mrSges(7,3) * t431;
t487 = m(7) * t373 + t425 * mrSges(7,1) + t442 * t412;
t370 = -mrSges(7,3) * t400 - t409 * t432 + t487;
t378 = t467 * t381 + t464 * t384;
t399 = -qJD(5) * t432 + qJDD(3) * t467 - t427 * t464;
t414 = pkin(5) * t442 - qJ(6) * t432;
t430 = t431 ^ 2;
t375 = -pkin(5) * t430 + qJ(6) * t399 + 0.2e1 * qJD(6) * t431 - t414 * t442 + t378;
t494 = t500 * t431 + t506 * t432 + t499 * t442;
t495 = -t505 * t431 - t500 * t432 - t498 * t442;
t503 = mrSges(6,1) * t377 + mrSges(7,1) * t373 - mrSges(6,2) * t378 - mrSges(7,2) * t375 + pkin(5) * t370 + t498 * t399 + t499 * t400 + t504 * t425 - t494 * t431 - t495 * t432;
t501 = -mrSges(6,2) - mrSges(7,2);
t421 = mrSges(5,1) * t444 + mrSges(5,2) * t445;
t437 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t445;
t410 = -mrSges(6,1) * t431 + mrSges(6,2) * t432;
t413 = -mrSges(6,2) * t442 + mrSges(6,3) * t431;
t364 = m(6) * t377 + mrSges(6,1) * t425 + t413 * t442 + (-t409 - t410) * t432 + (-mrSges(6,3) - mrSges(7,3)) * t400 + t487;
t486 = m(7) * t375 + t399 * mrSges(7,3) + t431 * t409;
t415 = mrSges(7,1) * t442 - mrSges(7,3) * t432;
t493 = -mrSges(6,1) * t442 + mrSges(6,3) * t432 - t415;
t367 = m(6) * t378 + mrSges(6,3) * t399 + t410 * t431 + t501 * t425 + t493 * t442 + t486;
t482 = -t364 * t464 + t467 * t367;
t359 = m(5) * t386 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t426 - qJD(3) * t437 - t421 * t444 + t482;
t436 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t444;
t380 = -qJDD(3) * pkin(4) - pkin(8) * t470 + t445 * t422 - t385;
t376 = -pkin(5) * t399 - qJ(6) * t430 + t414 * t432 + qJDD(6) + t380;
t480 = -m(7) * t376 + t399 * mrSges(7,1) + t431 * t412;
t472 = -m(6) * t380 + t399 * mrSges(6,1) + t501 * t400 + t431 * t413 + t493 * t432 + t480;
t369 = m(5) * t385 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t427 + qJD(3) * t436 - t421 * t445 + t472;
t355 = t462 * t359 + t463 * t369;
t362 = t467 * t364 + t464 * t367;
t496 = -t498 * t431 - t499 * t432 - t504 * t442;
t483 = t463 * t359 - t369 * t462;
t450 = (t465 * mrSges(4,1) + t468 * mrSges(4,2)) * qJD(1);
t453 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t492;
t455 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t491;
t478 = (m(4) * t428 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t452 + qJD(3) * t453 - t450 * t491 + t355) * t468 + t465 * (m(4) * t429 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t451 - qJD(3) * t455 - t450 * t492 + t483);
t361 = m(5) * t408 - mrSges(5,1) * t426 + t427 * mrSges(5,2) + t436 * t444 + t445 * t437 + t362;
t448 = Ifges(4,5) * qJD(3) + (t468 * Ifges(4,1) - t465 * Ifges(4,4)) * qJD(1);
t447 = Ifges(4,6) * qJD(3) + (t468 * Ifges(4,4) - t465 * Ifges(4,2)) * qJD(1);
t443 = -qJDD(1) * pkin(1) + t474;
t439 = pkin(1) * t471 - t475;
t435 = t502 * t471 + t475;
t419 = Ifges(5,1) * t445 - Ifges(5,4) * t444 + Ifges(5,5) * qJD(3);
t418 = Ifges(5,4) * t445 - Ifges(5,2) * t444 + Ifges(5,6) * qJD(3);
t417 = Ifges(5,5) * t445 - Ifges(5,6) * t444 + Ifges(5,3) * qJD(3);
t371 = mrSges(7,2) * t400 + t415 * t432 - t480;
t360 = mrSges(6,2) * t380 + mrSges(7,2) * t376 - mrSges(6,3) * t377 - mrSges(7,3) * t373 - qJ(6) * t370 + t500 * t399 + t506 * t400 + t499 * t425 - t496 * t431 + t495 * t442;
t356 = -mrSges(6,1) * t380 + mrSges(6,3) * t378 - mrSges(7,1) * t376 + mrSges(7,3) * t375 - pkin(5) * t371 + qJ(6) * t486 + (-qJ(6) * t415 + t494) * t442 + t496 * t432 + (-qJ(6) * mrSges(7,2) + t498) * t425 + t500 * t400 + t505 * t399;
t352 = -mrSges(5,1) * t408 + mrSges(5,3) * t386 + Ifges(5,4) * t427 + Ifges(5,2) * t426 + Ifges(5,6) * qJDD(3) - pkin(4) * t362 + qJD(3) * t419 - t445 * t417 - t503;
t351 = m(3) * t443 + qJDD(1) * mrSges(3,2) - mrSges(3,3) * t471 + t478;
t350 = mrSges(5,2) * t408 - mrSges(5,3) * t385 + Ifges(5,1) * t427 + Ifges(5,4) * t426 + Ifges(5,5) * qJDD(3) - pkin(8) * t362 - qJD(3) * t418 - t356 * t464 + t360 * t467 - t417 * t444;
t1 = [mrSges(2,1) * t484 - mrSges(2,2) * t479 + mrSges(3,2) * t443 - mrSges(3,3) * t439 + t468 * (mrSges(4,2) * t435 - mrSges(4,3) * t428 + Ifges(4,1) * t452 + Ifges(4,4) * t451 + Ifges(4,5) * qJDD(3) - qJ(4) * t355 - qJD(3) * t447 + t350 * t463 - t352 * t462) - t465 * (-mrSges(4,1) * t435 + mrSges(4,3) * t429 + Ifges(4,4) * t452 + Ifges(4,2) * t451 + Ifges(4,6) * qJDD(3) - pkin(3) * t361 + qJ(4) * t483 + qJD(3) * t448 + t462 * t350 + t463 * t352) - pkin(7) * t478 - pkin(1) * t351 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-m(3) * t439 + m(4) * t435 - mrSges(4,1) * t451 + mrSges(3,2) * t471 + mrSges(4,2) * t452 + t361 + qJDD(1) * mrSges(3,3) + (t453 * t465 + t455 * t468) * qJD(1)) * qJ(2); t351; Ifges(4,5) * t452 + Ifges(4,6) * t451 + mrSges(4,1) * t428 - mrSges(4,2) * t429 + Ifges(5,5) * t427 + Ifges(5,6) * t426 + t445 * t418 + t444 * t419 + mrSges(5,1) * t385 - mrSges(5,2) * t386 + t464 * t360 + t467 * t356 + pkin(4) * t472 + pkin(8) * t482 + pkin(3) * t355 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (t468 * t447 + t465 * t448) * qJD(1); t361; t503; t371;];
tauJ  = t1;
