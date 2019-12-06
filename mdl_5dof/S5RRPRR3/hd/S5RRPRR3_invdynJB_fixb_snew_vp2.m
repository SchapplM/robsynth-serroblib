% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPRR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPRR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR3_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR3_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR3_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:30:29
% EndTime: 2019-12-05 18:30:31
% DurationCPUTime: 2.50s
% Computational Cost: add. (44455->186), mult. (53565->233), div. (0->0), fcn. (28670->10), ass. (0->85)
t503 = sin(qJ(1));
t507 = cos(qJ(1));
t481 = t507 * g(2) + t503 * g(3);
t477 = qJDD(1) * pkin(1) + t481;
t480 = t503 * g(2) - t507 * g(3);
t508 = qJD(1) ^ 2;
t478 = -t508 * pkin(1) + t480;
t502 = sin(qJ(2));
t506 = cos(qJ(2));
t462 = t506 * t477 - t502 * t478;
t494 = qJDD(1) + qJDD(2);
t459 = t494 * pkin(2) + t462;
t463 = t502 * t477 + t506 * t478;
t495 = qJD(1) + qJD(2);
t493 = t495 ^ 2;
t460 = -t493 * pkin(2) + t463;
t498 = sin(pkin(9));
t499 = cos(pkin(9));
t454 = t499 * t459 - t498 * t460;
t451 = t494 * pkin(3) + t454;
t455 = t498 * t459 + t499 * t460;
t452 = -t493 * pkin(3) + t455;
t501 = sin(qJ(4));
t505 = cos(qJ(4));
t448 = t501 * t451 + t505 * t452;
t487 = qJD(4) + t495;
t485 = t487 ^ 2;
t486 = qJDD(4) + t494;
t445 = -t485 * pkin(4) + t486 * pkin(8) + t448;
t497 = -g(1) + qJDD(3);
t500 = sin(qJ(5));
t504 = cos(qJ(5));
t442 = -t500 * t445 + t504 * t497;
t443 = t504 * t445 + t500 * t497;
t465 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t500 + Ifges(6,2) * t504) * t487;
t466 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t500 + Ifges(6,4) * t504) * t487;
t522 = qJD(5) * t487;
t470 = t500 * t486 + t504 * t522;
t471 = t504 * t486 - t500 * t522;
t525 = mrSges(6,1) * t442 - mrSges(6,2) * t443 + Ifges(6,5) * t470 + Ifges(6,6) * t471 + Ifges(6,3) * qJDD(5) + (t465 * t500 - t466 * t504) * t487;
t524 = t487 * t500;
t523 = t487 * t504;
t469 = (-mrSges(6,1) * t504 + mrSges(6,2) * t500) * t487;
t476 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t523;
t440 = m(6) * t442 + qJDD(5) * mrSges(6,1) - t470 * mrSges(6,3) + qJD(5) * t476 - t469 * t524;
t475 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t524;
t441 = m(6) * t443 - qJDD(5) * mrSges(6,2) + t471 * mrSges(6,3) - qJD(5) * t475 + t469 * t523;
t516 = -t500 * t440 + t504 * t441;
t426 = m(5) * t448 - t485 * mrSges(5,1) - t486 * mrSges(5,2) + t516;
t447 = t505 * t451 - t501 * t452;
t444 = -t486 * pkin(4) - t485 * pkin(8) - t447;
t512 = -m(6) * t444 + t471 * mrSges(6,1) - t470 * mrSges(6,2) - t475 * t524 + t476 * t523;
t435 = m(5) * t447 + t486 * mrSges(5,1) - t485 * mrSges(5,2) + t512;
t423 = t501 * t426 + t505 * t435;
t419 = m(4) * t454 + t494 * mrSges(4,1) - t493 * mrSges(4,2) + t423;
t517 = t505 * t426 - t501 * t435;
t420 = m(4) * t455 - t493 * mrSges(4,1) - t494 * mrSges(4,2) + t517;
t414 = t499 * t419 + t498 * t420;
t411 = m(3) * t462 + t494 * mrSges(3,1) - t493 * mrSges(3,2) + t414;
t518 = -t498 * t419 + t499 * t420;
t412 = m(3) * t463 - t493 * mrSges(3,1) - t494 * mrSges(3,2) + t518;
t405 = t506 * t411 + t502 * t412;
t429 = t504 * t440 + t500 * t441;
t521 = m(5) * t497 + t429;
t519 = -t502 * t411 + t506 * t412;
t402 = m(2) * t480 - t508 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t519;
t403 = m(2) * t481 + qJDD(1) * mrSges(2,1) - t508 * mrSges(2,2) + t405;
t520 = t507 * t402 - t503 * t403;
t427 = m(4) * t497 + t521;
t515 = -t503 * t402 - t507 * t403;
t464 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t500 + Ifges(6,6) * t504) * t487;
t432 = -mrSges(6,1) * t444 + mrSges(6,3) * t443 + Ifges(6,4) * t470 + Ifges(6,2) * t471 + Ifges(6,6) * qJDD(5) + qJD(5) * t466 - t464 * t524;
t433 = mrSges(6,2) * t444 - mrSges(6,3) * t442 + Ifges(6,1) * t470 + Ifges(6,4) * t471 + Ifges(6,5) * qJDD(5) - qJD(5) * t465 + t464 * t523;
t513 = mrSges(5,1) * t447 - mrSges(5,2) * t448 + Ifges(5,3) * t486 + pkin(4) * t512 + pkin(8) * t516 + t504 * t432 + t500 * t433;
t510 = mrSges(3,1) * t462 + mrSges(4,1) * t454 - mrSges(3,2) * t463 - mrSges(4,2) * t455 + pkin(2) * t414 + pkin(3) * t423 + t513 + (Ifges(3,3) + Ifges(4,3)) * t494;
t509 = mrSges(2,1) * t481 - mrSges(2,2) * t480 + Ifges(2,3) * qJDD(1) + pkin(1) * t405 + t510;
t421 = -mrSges(5,1) * t497 + mrSges(5,3) * t448 + t485 * Ifges(5,5) + Ifges(5,6) * t486 - pkin(4) * t429 - t525;
t415 = mrSges(5,2) * t497 - mrSges(5,3) * t447 + Ifges(5,5) * t486 - t485 * Ifges(5,6) - pkin(8) * t429 - t500 * t432 + t504 * t433;
t407 = mrSges(4,2) * t497 - mrSges(4,3) * t454 + Ifges(4,5) * t494 - t493 * Ifges(4,6) - pkin(7) * t423 + t505 * t415 - t501 * t421;
t406 = -mrSges(4,1) * t497 + mrSges(4,3) * t455 + t493 * Ifges(4,5) + Ifges(4,6) * t494 - pkin(3) * t521 + pkin(7) * t517 + t501 * t415 + t505 * t421;
t400 = -mrSges(3,2) * g(1) - mrSges(3,3) * t462 + Ifges(3,5) * t494 - t493 * Ifges(3,6) - qJ(3) * t414 - t498 * t406 + t499 * t407;
t399 = mrSges(3,1) * g(1) + mrSges(3,3) * t463 + t493 * Ifges(3,5) + Ifges(3,6) * t494 - pkin(2) * t427 + qJ(3) * t518 + t499 * t406 + t498 * t407;
t398 = -mrSges(2,2) * g(1) - mrSges(2,3) * t481 + Ifges(2,5) * qJDD(1) - t508 * Ifges(2,6) - pkin(6) * t405 - t502 * t399 + t506 * t400;
t397 = Ifges(2,6) * qJDD(1) + t508 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t480 + t502 * t400 + t506 * t399 - pkin(1) * (-m(3) * g(1) + t427) + pkin(6) * t519;
t1 = [(-m(1) - m(2) - m(3)) * g(1) + t427; -m(1) * g(2) + t515; -m(1) * g(3) + t520; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t509; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t520 - t507 * t397 - t503 * t398; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t515 - t503 * t397 + t507 * t398; t509; t510; t427; t513; t525;];
tauJB = t1;
