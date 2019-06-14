% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRPRP10
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
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
% Datum: 2019-05-05 18:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRPRP10_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP10_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP10_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP10_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP10_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP10_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP10_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:09:45
% EndTime: 2019-05-05 18:09:49
% DurationCPUTime: 1.73s
% Computational Cost: add. (5314->247), mult. (10386->283), div. (0->0), fcn. (5167->6), ass. (0->106)
t533 = Ifges(6,1) + Ifges(7,1);
t513 = Ifges(6,4) - Ifges(7,5);
t528 = Ifges(7,4) + Ifges(6,5);
t532 = Ifges(6,2) + Ifges(7,3);
t527 = Ifges(6,6) - Ifges(7,6);
t473 = sin(qJ(3));
t475 = cos(qJ(3));
t510 = -Ifges(5,6) - Ifges(4,4);
t531 = t475 * (Ifges(5,2) + Ifges(4,1)) + t473 * t510;
t530 = t475 * t510 - t473 * (-Ifges(4,2) - Ifges(5,3));
t529 = -2 * qJD(4);
t512 = Ifges(4,5) - Ifges(5,4);
t511 = Ifges(4,6) - Ifges(5,5);
t526 = Ifges(6,3) + Ifges(7,2);
t472 = sin(qJ(5));
t503 = qJD(1) * t473;
t518 = cos(qJ(5));
t444 = qJD(3) * t472 - t518 * t503;
t445 = t518 * qJD(3) + t472 * t503;
t502 = qJD(1) * t475;
t462 = qJD(5) + t502;
t525 = t444 * t532 - t445 * t513 - t462 * t527;
t524 = -t513 * t444 + t445 * t533 + t528 * t462;
t521 = (-qJD(1) * t530 + t511 * qJD(3)) * t475;
t478 = qJD(1) ^ 2;
t474 = sin(qJ(1));
t476 = cos(qJ(1));
t495 = g(1) * t474 - t476 * g(2);
t486 = -qJ(2) * t478 + qJDD(2) - t495;
t520 = -pkin(1) - pkin(7);
t426 = t520 * qJDD(1) + t486;
t419 = -g(3) * t475 + t473 * t426;
t446 = (pkin(3) * t473 - qJ(4) * t475) * qJD(1);
t477 = qJD(3) ^ 2;
t393 = pkin(3) * t477 - qJDD(3) * qJ(4) + qJD(3) * t529 + t446 * t503 - t419;
t491 = -g(1) * t476 - g(2) * t474;
t487 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t491;
t519 = pkin(3) + pkin(8);
t517 = pkin(8) * t478;
t516 = g(3) * t473;
t515 = mrSges(4,1) - mrSges(5,2);
t514 = -mrSges(6,3) - mrSges(7,2);
t509 = t426 * t475;
t501 = qJD(1) * qJD(3);
t496 = t475 * t501;
t449 = qJDD(1) * t473 + t496;
t457 = pkin(4) * t502 - qJD(3) * pkin(8);
t471 = t473 ^ 2;
t464 = t473 * t501;
t450 = qJDD(1) * t475 - t464;
t482 = pkin(3) * t496 + t502 * t529 + t487 + (-t450 + t464) * qJ(4);
t386 = -t457 * t502 + t519 * t449 + (-pkin(4) * t471 + t520) * t478 + t482;
t485 = -qJ(4) * t477 + t446 * t502 + qJDD(4) - t509;
t390 = pkin(4) * t450 - t519 * qJDD(3) + (pkin(4) * t501 + t475 * t517 - g(3)) * t473 + t485;
t384 = t518 * t386 + t472 * t390;
t508 = t527 * t444 - t528 * t445 - t526 * t462;
t416 = mrSges(7,1) * t444 - mrSges(7,3) * t445;
t507 = -mrSges(6,1) * t444 - mrSges(6,2) * t445 - t416;
t505 = qJD(1) * t531 + t512 * qJD(3);
t455 = mrSges(5,1) * t503 - qJD(3) * mrSges(5,3);
t504 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t503 - t455;
t498 = t520 * t478;
t415 = pkin(5) * t444 - qJ(6) * t445;
t443 = qJDD(5) + t450;
t458 = t462 ^ 2;
t378 = -pkin(5) * t458 + qJ(6) * t443 + 0.2e1 * qJD(6) * t462 - t415 * t444 + t384;
t422 = -mrSges(7,1) * t462 + mrSges(7,2) * t445;
t497 = m(7) * t378 + t443 * mrSges(7,3) + t462 * t422;
t411 = qJD(5) * t445 + qJDD(3) * t472 - t518 * t449;
t421 = mrSges(6,1) * t462 - mrSges(6,3) * t445;
t370 = m(6) * t384 - mrSges(6,2) * t443 + t514 * t411 - t421 * t462 + t507 * t444 + t497;
t383 = -t472 * t386 + t518 * t390;
t412 = -t444 * qJD(5) + t518 * qJDD(3) + t472 * t449;
t420 = -mrSges(6,2) * t462 - mrSges(6,3) * t444;
t379 = -t443 * pkin(5) - t458 * qJ(6) + t445 * t415 + qJDD(6) - t383;
t423 = -mrSges(7,2) * t444 + mrSges(7,3) * t462;
t492 = -m(7) * t379 + t443 * mrSges(7,1) + t462 * t423;
t371 = m(6) * t383 + mrSges(6,1) * t443 + t514 * t412 + t420 * t462 + t507 * t445 + t492;
t494 = t518 * t370 - t371 * t472;
t447 = (-mrSges(5,2) * t473 - mrSges(5,3) * t475) * qJD(1);
t493 = qJD(1) * (-t447 - (mrSges(4,1) * t473 + mrSges(4,2) * t475) * qJD(1));
t418 = t509 + t516;
t454 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t502;
t456 = mrSges(5,1) * t502 + qJD(3) * mrSges(5,2);
t389 = -pkin(4) * t449 + qJD(3) * t457 - t471 * t517 - t393;
t381 = -0.2e1 * qJD(6) * t445 + (t444 * t462 - t412) * qJ(6) + (t445 * t462 + t411) * pkin(5) + t389;
t375 = m(7) * t381 + t411 * mrSges(7,1) - mrSges(7,3) * t412 - t422 * t445 + t444 * t423;
t481 = m(6) * t389 + mrSges(6,1) * t411 + t412 * mrSges(6,2) + t420 * t444 + t445 * t421 + t375;
t480 = -m(5) * t393 + qJDD(3) * mrSges(5,3) + qJD(3) * t456 + t481;
t367 = t472 * t370 + t518 * t371;
t394 = -qJDD(3) * pkin(3) + t485 - t516;
t484 = m(5) * t394 + t450 * mrSges(5,1) + t367;
t490 = (m(4) * t418 - t450 * mrSges(4,3) + t504 * qJD(3) + t515 * qJDD(3) + t475 * t493 - t484) * t475 + (-qJDD(3) * mrSges(4,2) + (-mrSges(4,3) - mrSges(5,1)) * t449 + m(4) * t419 - qJD(3) * t454 + t480 + t473 * t493) * t473;
t392 = pkin(3) * t449 + t482 + t498;
t488 = m(5) * t392 - t450 * mrSges(5,3) + t494;
t374 = mrSges(7,2) * t412 + t416 * t445 - t492;
t479 = mrSges(6,1) * t383 - mrSges(7,1) * t379 - mrSges(6,2) * t384 + mrSges(7,3) * t378 - pkin(5) * t374 + qJ(6) * t497 - t525 * t445 + (-qJ(6) * t416 + t524) * t444 + t526 * t443 + t528 * t412 + (-qJ(6) * mrSges(7,2) - t527) * t411;
t428 = -qJDD(1) * pkin(1) + t486;
t427 = pkin(1) * t478 - t487;
t425 = t498 + t487;
t366 = mrSges(6,2) * t389 + mrSges(7,2) * t379 - mrSges(6,3) * t383 - mrSges(7,3) * t381 - qJ(6) * t375 - t513 * t411 + t412 * t533 + t528 * t443 + t508 * t444 + t525 * t462;
t365 = -mrSges(6,1) * t389 - mrSges(7,1) * t381 + mrSges(7,2) * t378 + mrSges(6,3) * t384 - pkin(5) * t375 - t411 * t532 + t513 * t412 + t527 * t443 + t508 * t445 + t524 * t462;
t364 = qJDD(3) * mrSges(5,2) + qJD(3) * t455 + t447 * t502 + t484;
t363 = -mrSges(5,2) * t449 + (-t455 * t473 - t456 * t475) * qJD(1) + t488;
t361 = m(3) * t428 + qJDD(1) * mrSges(3,2) - mrSges(3,3) * t478 + t490;
t1 = [mrSges(2,1) * t495 - mrSges(2,2) * t491 + mrSges(3,2) * t428 - mrSges(3,3) * t427 + t475 * (mrSges(5,1) * t394 + mrSges(4,2) * t425 - mrSges(4,3) * t418 - mrSges(5,3) * t392 + pkin(4) * t367 - qJ(4) * t363 + t479) - t473 * (-mrSges(4,1) * t425 - mrSges(5,1) * t393 + mrSges(5,2) * t392 + mrSges(4,3) * t419 - pkin(3) * t363 + pkin(4) * t481 - pkin(8) * t494 - t518 * t365 - t472 * t366) - pkin(7) * t490 - pkin(1) * t361 + t531 * t450 + (-t473 * t511 + t475 * t512) * qJDD(3) + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-t473 * t505 - t521) * qJD(3) + t530 * t449 + (-m(3) * t427 + m(4) * t425 + mrSges(3,2) * t478 + t488 + mrSges(4,2) * t450 + qJDD(1) * mrSges(3,3) + t515 * t449 + ((t454 - t456) * t475 + t504 * t473) * qJD(1)) * qJ(2); t361; mrSges(4,1) * t418 - mrSges(4,2) * t419 + mrSges(5,2) * t394 - mrSges(5,3) * t393 + t518 * t366 - t472 * t365 - pkin(8) * t367 - pkin(3) * t364 + qJ(4) * t480 + t512 * t450 + (-qJ(4) * mrSges(5,1) - t511) * t449 + (Ifges(4,3) + Ifges(5,1)) * qJDD(3) + (t521 + (-qJ(4) * t447 + t505) * t473) * qJD(1); t364; t479; t374;];
tauJ  = t1;
