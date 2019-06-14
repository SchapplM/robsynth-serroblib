% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-05-05 21:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRPP2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:20:53
% EndTime: 2019-05-05 21:20:57
% DurationCPUTime: 1.52s
% Computational Cost: add. (7248->246), mult. (13682->287), div. (0->0), fcn. (8001->8), ass. (0->98)
t501 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t486 = Ifges(5,4) - Ifges(6,5) - Ifges(7,4);
t485 = Ifges(5,5) + Ifges(6,4) - Ifges(7,5);
t500 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t484 = Ifges(5,6) - Ifges(6,6) + Ifges(7,6);
t499 = Ifges(7,3) + Ifges(5,3) + Ifges(6,2);
t462 = sin(qJ(1));
t464 = cos(qJ(1));
t477 = t462 * g(1) - g(2) * t464;
t439 = qJDD(1) * pkin(1) + t477;
t466 = qJD(1) ^ 2;
t472 = -g(1) * t464 - g(2) * t462;
t441 = -pkin(1) * t466 + t472;
t458 = sin(pkin(9));
t459 = cos(pkin(9));
t408 = t439 * t459 - t458 * t441;
t385 = -qJDD(1) * pkin(2) - pkin(7) * t466 - t408;
t461 = sin(qJ(3));
t463 = cos(qJ(3));
t487 = qJD(1) * qJD(3);
t478 = t463 * t487;
t443 = qJDD(1) * t461 + t478;
t453 = t461 * t487;
t444 = qJDD(1) * t463 - t453;
t372 = (-t443 - t478) * pkin(8) + (-t444 + t453) * pkin(3) + t385;
t409 = t458 * t439 + t459 * t441;
t386 = -pkin(2) * t466 + qJDD(1) * pkin(7) + t409;
t457 = -g(3) + qJDD(2);
t378 = t463 * t386 + t461 * t457;
t442 = (-t463 * pkin(3) - t461 * pkin(8)) * qJD(1);
t465 = qJD(3) ^ 2;
t488 = qJD(1) * t463;
t376 = -pkin(3) * t465 + qJDD(3) * pkin(8) + t442 * t488 + t378;
t460 = sin(qJ(4));
t494 = cos(qJ(4));
t369 = t494 * t372 - t460 * t376;
t489 = qJD(1) * t461;
t437 = -t494 * qJD(3) + t460 * t489;
t438 = t460 * qJD(3) + t494 * t489;
t412 = pkin(4) * t437 - qJ(5) * t438;
t436 = -qJDD(4) + t444;
t449 = -qJD(4) + t488;
t448 = t449 ^ 2;
t367 = t436 * pkin(4) - t448 * qJ(5) + t438 * t412 + qJDD(5) - t369;
t422 = -mrSges(6,2) * t437 - mrSges(6,3) * t449;
t498 = -m(6) * t367 - t436 * mrSges(6,1) - t449 * t422;
t407 = -t437 * qJD(4) + t460 * qJDD(3) + t494 * t443;
t377 = -t461 * t386 + t463 * t457;
t470 = qJDD(3) * pkin(3) + pkin(8) * t465 - t442 * t489 + t377;
t491 = t437 * t449;
t497 = -(t407 + t491) * qJ(5) - t470;
t416 = -mrSges(7,2) * t449 + mrSges(7,3) * t437;
t360 = -0.2e1 * qJD(6) * t438 + (-t407 + t491) * qJ(6) + (t437 * t438 + t436) * pkin(5) + t367;
t414 = -mrSges(7,1) * t437 + mrSges(7,2) * t438;
t473 = -m(7) * t360 + t407 * mrSges(7,3) + t438 * t414;
t358 = mrSges(7,1) * t436 + t416 * t449 - t473;
t413 = mrSges(6,1) * t437 - mrSges(6,3) * t438;
t355 = mrSges(6,2) * t407 + t413 * t438 + t358 - t498;
t370 = t460 * t372 + t494 * t376;
t495 = -2 * qJD(5);
t366 = -pkin(4) * t448 - t436 * qJ(5) - t437 * t412 + t449 * t495 + t370;
t406 = qJD(4) * t438 - t494 * qJDD(3) + t443 * t460;
t418 = pkin(5) * t449 - qJ(6) * t438;
t435 = t437 ^ 2;
t362 = -pkin(5) * t435 + qJ(6) * t406 + 0.2e1 * qJD(6) * t437 - t418 * t449 + t366;
t419 = mrSges(7,1) * t449 - mrSges(7,3) * t438;
t421 = mrSges(6,1) * t449 + mrSges(6,2) * t438;
t482 = m(7) * t362 + t406 * mrSges(7,3) + t437 * t414;
t471 = m(6) * t366 - t436 * mrSges(6,3) - t449 * t421 + t482;
t479 = -t486 * t437 + t501 * t438 - t485 * t449;
t480 = t500 * t437 + t486 * t438 - t484 * t449;
t496 = -t406 * t484 + t407 * t485 - t499 * t436 + t437 * t479 + t438 * t480 + mrSges(5,1) * t369 - mrSges(6,1) * t367 - mrSges(7,1) * t360 - mrSges(5,2) * t370 + mrSges(7,2) * t362 + mrSges(6,3) * t366 - pkin(4) * t355 - pkin(5) * t358 + qJ(5) * (-mrSges(6,2) * t406 - mrSges(7,2) * t436 - t413 * t437 - t419 * t449 + t471);
t492 = -mrSges(5,3) - mrSges(6,2);
t490 = -mrSges(5,1) * t437 - mrSges(5,2) * t438 - t413;
t481 = t484 * t437 - t485 * t438 + t499 * t449;
t440 = (-t463 * mrSges(4,1) + t461 * mrSges(4,2)) * qJD(1);
t445 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t489;
t417 = mrSges(5,2) * t449 - mrSges(5,3) * t437;
t353 = m(5) * t369 + (-t416 - t417) * t449 + t490 * t438 + (-mrSges(5,1) - mrSges(7,1)) * t436 + t492 * t407 + t473 + t498;
t420 = -mrSges(5,1) * t449 - mrSges(5,3) * t438;
t354 = m(5) * t370 + (-t419 + t420) * t449 + t490 * t437 + (mrSges(5,2) - mrSges(7,2)) * t436 + t492 * t406 + t471;
t475 = -t353 * t460 + t494 * t354;
t349 = m(4) * t378 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t444 - qJD(3) * t445 + t440 * t488 + t475;
t446 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t488;
t364 = -qJ(6) * t435 + qJDD(6) + (-pkin(4) - pkin(5)) * t406 + (pkin(4) * t449 + (2 * qJD(5)) + t418) * t438 - t497;
t359 = m(7) * t364 - t406 * mrSges(7,1) + t407 * mrSges(7,2) - t437 * t416 + t438 * t419;
t368 = t438 * t495 + (-t438 * t449 + t406) * pkin(4) + t497;
t357 = m(6) * t368 + t406 * mrSges(6,1) - t407 * mrSges(6,3) - t438 * t421 + t437 * t422 - t359;
t467 = m(5) * t470 - t406 * mrSges(5,1) - t407 * mrSges(5,2) - t437 * t417 - t438 * t420 - t357;
t351 = m(4) * t377 + qJDD(3) * mrSges(4,1) - t443 * mrSges(4,3) + qJD(3) * t446 - t440 * t489 + t467;
t476 = t463 * t349 - t351 * t461;
t350 = t494 * t353 + t460 * t354;
t469 = -m(4) * t385 + t444 * mrSges(4,1) - t443 * mrSges(4,2) - t445 * t489 + t446 * t488 - t350;
t429 = Ifges(4,5) * qJD(3) + (t461 * Ifges(4,1) + t463 * Ifges(4,4)) * qJD(1);
t428 = Ifges(4,6) * qJD(3) + (t461 * Ifges(4,4) + t463 * Ifges(4,2)) * qJD(1);
t347 = -mrSges(5,2) * t470 + mrSges(6,2) * t367 + mrSges(7,2) * t364 - mrSges(5,3) * t369 - mrSges(6,3) * t368 - mrSges(7,3) * t360 - qJ(5) * t357 - qJ(6) * t358 - t486 * t406 + t501 * t407 - t485 * t436 + t481 * t437 + t480 * t449;
t346 = mrSges(5,1) * t470 + mrSges(5,3) * t370 - mrSges(6,1) * t368 + mrSges(6,2) * t366 + mrSges(7,1) * t364 - mrSges(7,3) * t362 + pkin(5) * t359 - qJ(6) * t482 - pkin(4) * t357 + (qJ(6) * t419 - t479) * t449 + t481 * t438 + (qJ(6) * mrSges(7,2) - t484) * t436 + t486 * t407 + t500 * t406;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t477 - mrSges(2,2) * t472 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t408 - mrSges(3,2) * t409 + t461 * (mrSges(4,2) * t385 - mrSges(4,3) * t377 + Ifges(4,1) * t443 + Ifges(4,4) * t444 + Ifges(4,5) * qJDD(3) - pkin(8) * t350 - qJD(3) * t428 - t460 * t346 + t347 * t494) + pkin(2) * t469 + pkin(7) * t476 + pkin(1) * (t458 * (m(3) * t409 - mrSges(3,1) * t466 - qJDD(1) * mrSges(3,2) + t476) + t459 * (m(3) * t408 + qJDD(1) * mrSges(3,1) - t466 * mrSges(3,2) + t469)) + (-mrSges(4,1) * t385 + mrSges(4,3) * t378 + Ifges(4,4) * t443 + Ifges(4,2) * t444 + Ifges(4,6) * qJDD(3) - pkin(3) * t350 + qJD(3) * t429 - t496) * t463; m(3) * t457 + t349 * t461 + t351 * t463; Ifges(4,5) * t443 + Ifges(4,6) * t444 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t377 - mrSges(4,2) * t378 + t460 * t347 + t494 * t346 + pkin(3) * t467 + pkin(8) * t475 + (t461 * t428 - t463 * t429) * qJD(1); t496; t355; t359;];
tauJ  = t1;
