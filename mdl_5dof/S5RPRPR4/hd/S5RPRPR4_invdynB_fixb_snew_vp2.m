% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRPR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR4_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR4_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR4_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:38:45
% EndTime: 2020-01-03 11:38:52
% DurationCPUTime: 4.92s
% Computational Cost: add. (53872->269), mult. (115951->344), div. (0->0), fcn. (73622->10), ass. (0->106)
t510 = sin(qJ(1));
t513 = cos(qJ(1));
t494 = -t510 * g(2) + t513 * g(3);
t514 = qJD(1) ^ 2;
t495 = -t513 * g(2) - t510 * g(3);
t486 = qJDD(1) * pkin(1) + t495;
t488 = -t514 * pkin(1) + t494;
t505 = sin(pkin(8));
t507 = cos(pkin(8));
t466 = t505 * t486 + t507 * t488;
t461 = -t514 * pkin(2) + qJDD(1) * pkin(6) + t466;
t503 = -g(1) + qJDD(2);
t509 = sin(qJ(3));
t512 = cos(qJ(3));
t450 = -t509 * t461 + t512 * t503;
t526 = qJD(1) * qJD(3);
t524 = t512 * t526;
t489 = t509 * qJDD(1) + t524;
t447 = (-t489 + t524) * qJ(4) + (t509 * t512 * t514 + qJDD(3)) * pkin(3) + t450;
t451 = t512 * t461 + t509 * t503;
t490 = t512 * qJDD(1) - t509 * t526;
t528 = qJD(1) * t509;
t491 = qJD(3) * pkin(3) - qJ(4) * t528;
t502 = t512 ^ 2;
t448 = -t502 * t514 * pkin(3) + t490 * qJ(4) - qJD(3) * t491 + t451;
t504 = sin(pkin(9));
t506 = cos(pkin(9));
t476 = (t504 * t512 + t506 * t509) * qJD(1);
t430 = -0.2e1 * qJD(4) * t476 + t506 * t447 - t504 * t448;
t468 = t506 * t489 + t504 * t490;
t475 = (-t504 * t509 + t506 * t512) * qJD(1);
t428 = (qJD(3) * t475 - t468) * pkin(7) + (t475 * t476 + qJDD(3)) * pkin(4) + t430;
t431 = 0.2e1 * qJD(4) * t475 + t504 * t447 + t506 * t448;
t467 = -t504 * t489 + t506 * t490;
t471 = qJD(3) * pkin(4) - t476 * pkin(7);
t474 = t475 ^ 2;
t429 = -t474 * pkin(4) + t467 * pkin(7) - qJD(3) * t471 + t431;
t508 = sin(qJ(5));
t511 = cos(qJ(5));
t426 = t511 * t428 - t508 * t429;
t458 = t511 * t475 - t508 * t476;
t437 = t458 * qJD(5) + t508 * t467 + t511 * t468;
t459 = t508 * t475 + t511 * t476;
t446 = -t458 * mrSges(6,1) + t459 * mrSges(6,2);
t501 = qJD(3) + qJD(5);
t452 = -t501 * mrSges(6,2) + t458 * mrSges(6,3);
t500 = qJDD(3) + qJDD(5);
t424 = m(6) * t426 + t500 * mrSges(6,1) - t437 * mrSges(6,3) - t459 * t446 + t501 * t452;
t427 = t508 * t428 + t511 * t429;
t436 = -t459 * qJD(5) + t511 * t467 - t508 * t468;
t453 = t501 * mrSges(6,1) - t459 * mrSges(6,3);
t425 = m(6) * t427 - t500 * mrSges(6,2) + t436 * mrSges(6,3) + t458 * t446 - t501 * t453;
t416 = t511 * t424 + t508 * t425;
t463 = -t475 * mrSges(5,1) + t476 * mrSges(5,2);
t469 = -qJD(3) * mrSges(5,2) + t475 * mrSges(5,3);
t414 = m(5) * t430 + qJDD(3) * mrSges(5,1) - t468 * mrSges(5,3) + qJD(3) * t469 - t476 * t463 + t416;
t470 = qJD(3) * mrSges(5,1) - t476 * mrSges(5,3);
t519 = -t508 * t424 + t511 * t425;
t415 = m(5) * t431 - qJDD(3) * mrSges(5,2) + t467 * mrSges(5,3) - qJD(3) * t470 + t475 * t463 + t519;
t410 = t506 * t414 + t504 * t415;
t487 = (-mrSges(4,1) * t512 + mrSges(4,2) * t509) * qJD(1);
t527 = qJD(1) * t512;
t493 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t527;
t408 = m(4) * t450 + qJDD(3) * mrSges(4,1) - t489 * mrSges(4,3) + qJD(3) * t493 - t487 * t528 + t410;
t492 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t528;
t520 = -t504 * t414 + t506 * t415;
t409 = m(4) * t451 - qJDD(3) * mrSges(4,2) + t490 * mrSges(4,3) - qJD(3) * t492 + t487 * t527 + t520;
t521 = -t509 * t408 + t512 * t409;
t401 = m(3) * t466 - t514 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t521;
t465 = t507 * t486 - t505 * t488;
t517 = -qJDD(1) * pkin(2) - t465;
t460 = -t514 * pkin(6) + t517;
t449 = -t490 * pkin(3) + qJDD(4) + t491 * t528 + (-qJ(4) * t502 - pkin(6)) * t514 + t517;
t433 = -t467 * pkin(4) - t474 * pkin(7) + t476 * t471 + t449;
t518 = m(6) * t433 - t436 * mrSges(6,1) + t437 * mrSges(6,2) - t458 * t452 + t459 * t453;
t516 = m(5) * t449 - t467 * mrSges(5,1) + t468 * mrSges(5,2) - t475 * t469 + t476 * t470 + t518;
t515 = -m(4) * t460 + t490 * mrSges(4,1) - t489 * mrSges(4,2) - t492 * t528 + t493 * t527 - t516;
t420 = m(3) * t465 + qJDD(1) * mrSges(3,1) - t514 * mrSges(3,2) + t515;
t522 = t507 * t401 - t505 * t420;
t396 = m(2) * t494 - t514 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t522;
t398 = t505 * t401 + t507 * t420;
t397 = m(2) * t495 + qJDD(1) * mrSges(2,1) - t514 * mrSges(2,2) + t398;
t529 = t510 * t396 + t513 * t397;
t402 = t512 * t408 + t509 * t409;
t525 = m(3) * t503 + t402;
t523 = -t513 * t396 + t510 * t397;
t482 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t509 + Ifges(4,4) * t512) * qJD(1);
t481 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t509 + Ifges(4,2) * t512) * qJD(1);
t480 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t509 + Ifges(4,6) * t512) * qJD(1);
t457 = Ifges(5,1) * t476 + Ifges(5,4) * t475 + Ifges(5,5) * qJD(3);
t456 = Ifges(5,4) * t476 + Ifges(5,2) * t475 + Ifges(5,6) * qJD(3);
t455 = Ifges(5,5) * t476 + Ifges(5,6) * t475 + Ifges(5,3) * qJD(3);
t440 = Ifges(6,1) * t459 + Ifges(6,4) * t458 + Ifges(6,5) * t501;
t439 = Ifges(6,4) * t459 + Ifges(6,2) * t458 + Ifges(6,6) * t501;
t438 = Ifges(6,5) * t459 + Ifges(6,6) * t458 + Ifges(6,3) * t501;
t418 = mrSges(6,2) * t433 - mrSges(6,3) * t426 + Ifges(6,1) * t437 + Ifges(6,4) * t436 + Ifges(6,5) * t500 + t458 * t438 - t501 * t439;
t417 = -mrSges(6,1) * t433 + mrSges(6,3) * t427 + Ifges(6,4) * t437 + Ifges(6,2) * t436 + Ifges(6,6) * t500 - t459 * t438 + t501 * t440;
t404 = mrSges(5,2) * t449 - mrSges(5,3) * t430 + Ifges(5,1) * t468 + Ifges(5,4) * t467 + Ifges(5,5) * qJDD(3) - pkin(7) * t416 - qJD(3) * t456 - t508 * t417 + t511 * t418 + t475 * t455;
t403 = -mrSges(5,1) * t449 + mrSges(5,3) * t431 + Ifges(5,4) * t468 + Ifges(5,2) * t467 + Ifges(5,6) * qJDD(3) - pkin(4) * t518 + pkin(7) * t519 + qJD(3) * t457 + t511 * t417 + t508 * t418 - t476 * t455;
t392 = mrSges(4,2) * t460 - mrSges(4,3) * t450 + Ifges(4,1) * t489 + Ifges(4,4) * t490 + Ifges(4,5) * qJDD(3) - qJ(4) * t410 - qJD(3) * t481 - t504 * t403 + t506 * t404 + t480 * t527;
t391 = -mrSges(4,1) * t460 + mrSges(4,3) * t451 + Ifges(4,4) * t489 + Ifges(4,2) * t490 + Ifges(4,6) * qJDD(3) - pkin(3) * t516 + qJ(4) * t520 + qJD(3) * t482 + t506 * t403 + t504 * t404 - t480 * t528;
t390 = -pkin(3) * t410 - mrSges(3,1) * t503 + t514 * Ifges(3,5) - Ifges(6,3) * t500 - Ifges(4,6) * t490 - Ifges(4,5) * t489 + t475 * t457 - t476 * t456 - Ifges(6,6) * t436 - Ifges(6,5) * t437 - mrSges(5,1) * t430 + mrSges(5,2) * t431 + Ifges(3,6) * qJDD(1) - pkin(4) * t416 - pkin(2) * t402 - t459 * t439 + mrSges(3,3) * t466 - Ifges(5,6) * t467 - Ifges(5,5) * t468 - mrSges(4,1) * t450 + mrSges(4,2) * t451 + t458 * t440 - mrSges(6,1) * t426 + mrSges(6,2) * t427 + (-t509 * t481 + t512 * t482) * qJD(1) + (-Ifges(4,3) - Ifges(5,3)) * qJDD(3);
t389 = mrSges(3,2) * t503 - mrSges(3,3) * t465 + Ifges(3,5) * qJDD(1) - t514 * Ifges(3,6) - pkin(6) * t402 - t509 * t391 + t512 * t392;
t388 = -mrSges(2,2) * g(1) - mrSges(2,3) * t495 + Ifges(2,5) * qJDD(1) - t514 * Ifges(2,6) - qJ(2) * t398 + t507 * t389 - t505 * t390;
t387 = mrSges(2,1) * g(1) + mrSges(2,3) * t494 + t514 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t525 + qJ(2) * t522 + t505 * t389 + t507 * t390;
t1 = [(-m(1) - m(2)) * g(1) + t525; -m(1) * g(2) + t529; -m(1) * g(3) + t523; pkin(1) * t398 + pkin(6) * t521 + t509 * t392 + t512 * t391 + pkin(2) * t515 + mrSges(3,1) * t465 - mrSges(3,2) * t466 + mrSges(2,1) * t495 - mrSges(2,2) * t494 - mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1); mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t523 + t513 * t387 + t510 * t388; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t529 + t510 * t387 - t513 * t388;];
tauB = t1;
