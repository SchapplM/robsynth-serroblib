% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-05-05 06:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRRPP2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPP2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:42:24
% EndTime: 2019-05-05 06:42:26
% DurationCPUTime: 1.39s
% Computational Cost: add. (7513->244), mult. (14214->288), div. (0->0), fcn. (9223->10), ass. (0->100)
t473 = sin(pkin(10));
t475 = cos(pkin(10));
t460 = g(1) * t473 - g(2) * t475;
t472 = -g(3) + qJDD(1);
t474 = sin(pkin(6));
t476 = cos(pkin(6));
t520 = t460 * t476 + t472 * t474;
t519 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t502 = Ifges(5,4) - Ifges(6,5) - Ifges(7,4);
t501 = Ifges(5,5) + Ifges(6,4) - Ifges(7,5);
t518 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t500 = -Ifges(6,6) + Ifges(7,6) + Ifges(5,6);
t517 = Ifges(7,3) + Ifges(5,3) + Ifges(6,2);
t461 = -g(1) * t475 - g(2) * t473;
t479 = sin(qJ(2));
t481 = cos(qJ(2));
t402 = t481 * t461 + t520 * t479;
t483 = qJD(2) ^ 2;
t396 = -pkin(2) * t483 + qJDD(2) * pkin(8) + t402;
t439 = -t460 * t474 + t472 * t476;
t478 = sin(qJ(3));
t480 = cos(qJ(3));
t392 = t480 * t396 + t478 * t439;
t457 = (-t480 * pkin(3) - t478 * pkin(9)) * qJD(2);
t482 = qJD(3) ^ 2;
t504 = qJD(2) * t480;
t388 = -pkin(3) * t482 + qJDD(3) * pkin(9) + t457 * t504 + t392;
t401 = -t479 * t461 + t520 * t481;
t395 = -qJDD(2) * pkin(2) - pkin(8) * t483 - t401;
t503 = qJD(2) * qJD(3);
t494 = t480 * t503;
t458 = qJDD(2) * t478 + t494;
t469 = t478 * t503;
t459 = qJDD(2) * t480 - t469;
t390 = (-t458 - t494) * pkin(9) + (-t459 + t469) * pkin(3) + t395;
t477 = sin(qJ(4));
t512 = cos(qJ(4));
t383 = -t477 * t388 + t512 * t390;
t505 = qJD(2) * t478;
t454 = -t512 * qJD(3) + t477 * t505;
t455 = t477 * qJD(3) + t512 * t505;
t426 = pkin(4) * t454 - qJ(5) * t455;
t451 = -qJDD(4) + t459;
t467 = -qJD(4) + t504;
t466 = t467 ^ 2;
t381 = t451 * pkin(4) - t466 * qJ(5) + t455 * t426 + qJDD(5) - t383;
t438 = -mrSges(6,2) * t454 - mrSges(6,3) * t467;
t516 = -m(6) * t381 - t451 * mrSges(6,1) - t467 * t438;
t423 = -t454 * qJD(4) + t477 * qJDD(3) + t512 * t458;
t391 = -t478 * t396 + t480 * t439;
t487 = qJDD(3) * pkin(3) + pkin(9) * t482 - t457 * t505 + t391;
t509 = t454 * t467;
t515 = -(t423 + t509) * qJ(5) - t487;
t432 = -mrSges(7,2) * t467 + mrSges(7,3) * t454;
t374 = -0.2e1 * qJD(6) * t455 + (-t423 + t509) * qJ(6) + (t454 * t455 + t451) * pkin(5) + t381;
t428 = -mrSges(7,1) * t454 + mrSges(7,2) * t455;
t490 = -m(7) * t374 + t423 * mrSges(7,3) + t455 * t428;
t372 = mrSges(7,1) * t451 + t432 * t467 - t490;
t427 = mrSges(6,1) * t454 - mrSges(6,3) * t455;
t369 = mrSges(6,2) * t423 + t427 * t455 + t372 - t516;
t384 = t512 * t388 + t477 * t390;
t513 = -2 * qJD(5);
t380 = -pkin(4) * t466 - t451 * qJ(5) - t454 * t426 + t467 * t513 + t384;
t422 = qJD(4) * t455 - t512 * qJDD(3) + t458 * t477;
t434 = pkin(5) * t467 - qJ(6) * t455;
t450 = t454 ^ 2;
t376 = -pkin(5) * t450 + qJ(6) * t422 + 0.2e1 * qJD(6) * t454 - t434 * t467 + t380;
t435 = mrSges(7,1) * t467 - mrSges(7,3) * t455;
t437 = mrSges(6,1) * t467 + mrSges(6,2) * t455;
t498 = m(7) * t376 + t422 * mrSges(7,3) + t454 * t428;
t488 = m(6) * t380 - t451 * mrSges(6,3) - t467 * t437 + t498;
t495 = -t502 * t454 + t519 * t455 - t501 * t467;
t496 = t518 * t454 + t502 * t455 - t500 * t467;
t514 = -t500 * t422 + t501 * t423 - t517 * t451 + t495 * t454 + t496 * t455 + mrSges(5,1) * t383 - mrSges(6,1) * t381 - mrSges(7,1) * t374 - mrSges(5,2) * t384 + mrSges(7,2) * t376 + mrSges(6,3) * t380 - pkin(4) * t369 - pkin(5) * t372 + qJ(5) * (-mrSges(6,2) * t422 - mrSges(7,2) * t451 - t427 * t454 - t435 * t467 + t488);
t510 = -mrSges(5,3) - mrSges(6,2);
t506 = -mrSges(5,1) * t454 - mrSges(5,2) * t455 - t427;
t497 = t500 * t454 - t501 * t455 + t517 * t467;
t456 = (-mrSges(4,1) * t480 + mrSges(4,2) * t478) * qJD(2);
t462 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t505;
t433 = mrSges(5,2) * t467 - mrSges(5,3) * t454;
t367 = m(5) * t383 + (-t432 - t433) * t467 + t506 * t455 + (-mrSges(5,1) - mrSges(7,1)) * t451 + t510 * t423 + t490 + t516;
t436 = -mrSges(5,1) * t467 - mrSges(5,3) * t455;
t368 = m(5) * t384 + (-t435 + t436) * t467 + t506 * t454 + (mrSges(5,2) - mrSges(7,2)) * t451 + t510 * t422 + t488;
t492 = -t367 * t477 + t512 * t368;
t363 = m(4) * t392 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t459 - qJD(3) * t462 + t456 * t504 + t492;
t463 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t504;
t378 = -qJ(6) * t450 + qJDD(6) + (-pkin(4) - pkin(5)) * t422 + (pkin(4) * t467 + (2 * qJD(5)) + t434) * t455 - t515;
t373 = m(7) * t378 - t422 * mrSges(7,1) + t423 * mrSges(7,2) - t454 * t432 + t455 * t435;
t382 = t455 * t513 + (-t455 * t467 + t422) * pkin(4) + t515;
t371 = m(6) * t382 + t422 * mrSges(6,1) - t423 * mrSges(6,3) - t455 * t437 + t454 * t438 - t373;
t484 = m(5) * t487 - t422 * mrSges(5,1) - t423 * mrSges(5,2) - t454 * t433 - t455 * t436 - t371;
t365 = m(4) * t391 + qJDD(3) * mrSges(4,1) - t458 * mrSges(4,3) + qJD(3) * t463 - t456 * t505 + t484;
t493 = t480 * t363 - t365 * t478;
t364 = t512 * t367 + t477 * t368;
t486 = -m(4) * t395 + t459 * mrSges(4,1) - t458 * mrSges(4,2) - t462 * t505 + t463 * t504 - t364;
t444 = Ifges(4,5) * qJD(3) + (t478 * Ifges(4,1) + t480 * Ifges(4,4)) * qJD(2);
t443 = Ifges(4,6) * qJD(3) + (t478 * Ifges(4,4) + t480 * Ifges(4,2)) * qJD(2);
t361 = -mrSges(5,2) * t487 + mrSges(6,2) * t381 + mrSges(7,2) * t378 - mrSges(5,3) * t383 - mrSges(6,3) * t382 - mrSges(7,3) * t374 - qJ(5) * t371 - qJ(6) * t372 - t502 * t422 + t519 * t423 - t501 * t451 + t497 * t454 + t496 * t467;
t360 = mrSges(5,1) * t487 + mrSges(5,3) * t384 - mrSges(6,1) * t382 + mrSges(6,2) * t380 + mrSges(7,1) * t378 - mrSges(7,3) * t376 + pkin(5) * t373 - qJ(6) * t498 - pkin(4) * t371 + (qJ(6) * t435 - t495) * t467 + t497 * t455 + (qJ(6) * mrSges(7,2) - t500) * t451 + t502 * t423 + t518 * t422;
t1 = [m(2) * t472 + t476 * (m(3) * t439 + t363 * t478 + t365 * t480) + (t479 * (m(3) * t402 - mrSges(3,1) * t483 - qJDD(2) * mrSges(3,2) + t493) + t481 * (m(3) * t401 + qJDD(2) * mrSges(3,1) - t483 * mrSges(3,2) + t486)) * t474; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t401 - mrSges(3,2) * t402 + t478 * (mrSges(4,2) * t395 - mrSges(4,3) * t391 + Ifges(4,1) * t458 + Ifges(4,4) * t459 + Ifges(4,5) * qJDD(3) - pkin(9) * t364 - qJD(3) * t443 - t477 * t360 + t512 * t361) + pkin(2) * t486 + pkin(8) * t493 + (-mrSges(4,1) * t395 + mrSges(4,3) * t392 + Ifges(4,4) * t458 + Ifges(4,2) * t459 + Ifges(4,6) * qJDD(3) - pkin(3) * t364 + qJD(3) * t444 - t514) * t480; Ifges(4,5) * t458 + Ifges(4,6) * t459 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t391 - mrSges(4,2) * t392 + t477 * t361 + t512 * t360 + pkin(3) * t484 + pkin(9) * t492 + (t443 * t478 - t444 * t480) * qJD(2); t514; t369; t373;];
tauJ  = t1;
