% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRPR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR2_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR2_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR2_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:40:55
% EndTime: 2019-12-05 18:40:58
% DurationCPUTime: 2.52s
% Computational Cost: add. (46436->187), mult. (53945->233), div. (0->0), fcn. (28876->10), ass. (0->85)
t498 = sin(qJ(1));
t502 = cos(qJ(1));
t477 = t502 * g(2) + t498 * g(3);
t473 = qJDD(1) * pkin(1) + t477;
t476 = t498 * g(2) - t502 * g(3);
t503 = qJD(1) ^ 2;
t474 = -t503 * pkin(1) + t476;
t497 = sin(qJ(2));
t501 = cos(qJ(2));
t458 = t501 * t473 - t497 * t474;
t489 = qJDD(1) + qJDD(2);
t455 = t489 * pkin(2) + t458;
t459 = t497 * t473 + t501 * t474;
t490 = qJD(1) + qJD(2);
t488 = t490 ^ 2;
t456 = -t488 * pkin(2) + t459;
t496 = sin(qJ(3));
t500 = cos(qJ(3));
t450 = t500 * t455 - t496 * t456;
t482 = qJDD(3) + t489;
t447 = t482 * pkin(3) + t450;
t451 = t496 * t455 + t500 * t456;
t483 = qJD(3) + t490;
t481 = t483 ^ 2;
t448 = -t481 * pkin(3) + t451;
t493 = sin(pkin(9));
t494 = cos(pkin(9));
t444 = t493 * t447 + t494 * t448;
t441 = -t481 * pkin(4) + t482 * pkin(8) + t444;
t492 = -g(1) + qJDD(4);
t495 = sin(qJ(5));
t499 = cos(qJ(5));
t438 = -t495 * t441 + t499 * t492;
t439 = t499 * t441 + t495 * t492;
t461 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t495 + Ifges(6,2) * t499) * t483;
t462 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t495 + Ifges(6,4) * t499) * t483;
t516 = qJD(5) * t483;
t466 = t495 * t482 + t499 * t516;
t467 = t499 * t482 - t495 * t516;
t520 = mrSges(6,1) * t438 - mrSges(6,2) * t439 + Ifges(6,5) * t466 + Ifges(6,6) * t467 + Ifges(6,3) * qJDD(5) + (t461 * t495 - t462 * t499) * t483;
t519 = -m(3) - m(4);
t518 = t483 * t495;
t517 = t483 * t499;
t465 = (-mrSges(6,1) * t499 + mrSges(6,2) * t495) * t483;
t472 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t517;
t436 = m(6) * t438 + qJDD(5) * mrSges(6,1) - t466 * mrSges(6,3) + qJD(5) * t472 - t465 * t518;
t471 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t518;
t437 = m(6) * t439 - qJDD(5) * mrSges(6,2) + t467 * mrSges(6,3) - qJD(5) * t471 + t465 * t517;
t511 = -t495 * t436 + t499 * t437;
t422 = m(5) * t444 - t481 * mrSges(5,1) - t482 * mrSges(5,2) + t511;
t443 = t494 * t447 - t493 * t448;
t440 = -t482 * pkin(4) - t481 * pkin(8) - t443;
t508 = -m(6) * t440 + t467 * mrSges(6,1) - t466 * mrSges(6,2) - t471 * t518 + t472 * t517;
t431 = m(5) * t443 + t482 * mrSges(5,1) - t481 * mrSges(5,2) + t508;
t419 = t493 * t422 + t494 * t431;
t415 = m(4) * t450 + t482 * mrSges(4,1) - t481 * mrSges(4,2) + t419;
t512 = t494 * t422 - t493 * t431;
t416 = m(4) * t451 - t481 * mrSges(4,1) - t482 * mrSges(4,2) + t512;
t410 = t500 * t415 + t496 * t416;
t407 = m(3) * t458 + t489 * mrSges(3,1) - t488 * mrSges(3,2) + t410;
t513 = -t496 * t415 + t500 * t416;
t408 = m(3) * t459 - t488 * mrSges(3,1) - t489 * mrSges(3,2) + t513;
t401 = t501 * t407 + t497 * t408;
t425 = t499 * t436 + t495 * t437;
t423 = m(5) * t492 + t425;
t514 = -t497 * t407 + t501 * t408;
t398 = m(2) * t476 - t503 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t514;
t399 = m(2) * t477 + qJDD(1) * mrSges(2,1) - t503 * mrSges(2,2) + t401;
t515 = t502 * t398 - t498 * t399;
t510 = -t498 * t398 - t502 * t399;
t460 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t495 + Ifges(6,6) * t499) * t483;
t428 = -mrSges(6,1) * t440 + mrSges(6,3) * t439 + Ifges(6,4) * t466 + Ifges(6,2) * t467 + Ifges(6,6) * qJDD(5) + qJD(5) * t462 - t460 * t518;
t429 = mrSges(6,2) * t440 - mrSges(6,3) * t438 + Ifges(6,1) * t466 + Ifges(6,4) * t467 + Ifges(6,5) * qJDD(5) - qJD(5) * t461 + t460 * t517;
t506 = mrSges(4,1) * t450 + mrSges(5,1) * t443 - mrSges(4,2) * t451 - mrSges(5,2) * t444 + pkin(3) * t419 + pkin(4) * t508 + pkin(8) * t511 + t499 * t428 + t495 * t429 + (Ifges(5,3) + Ifges(4,3)) * t482;
t505 = mrSges(3,1) * t458 - mrSges(3,2) * t459 + Ifges(3,3) * t489 + pkin(2) * t410 + t506;
t504 = mrSges(2,1) * t477 - mrSges(2,2) * t476 + Ifges(2,3) * qJDD(1) + pkin(1) * t401 + t505;
t417 = -mrSges(5,1) * t492 + mrSges(5,3) * t444 + t481 * Ifges(5,5) + Ifges(5,6) * t482 - pkin(4) * t425 - t520;
t411 = mrSges(5,2) * t492 - mrSges(5,3) * t443 + Ifges(5,5) * t482 - t481 * Ifges(5,6) - pkin(8) * t425 - t495 * t428 + t499 * t429;
t403 = -mrSges(4,2) * g(1) - mrSges(4,3) * t450 + Ifges(4,5) * t482 - t481 * Ifges(4,6) - qJ(4) * t419 + t494 * t411 - t493 * t417;
t402 = mrSges(4,1) * g(1) + mrSges(4,3) * t451 + t481 * Ifges(4,5) + Ifges(4,6) * t482 - pkin(3) * t423 + qJ(4) * t512 + t493 * t411 + t494 * t417;
t396 = -mrSges(3,2) * g(1) - mrSges(3,3) * t458 + Ifges(3,5) * t489 - t488 * Ifges(3,6) - pkin(7) * t410 - t496 * t402 + t500 * t403;
t395 = Ifges(3,6) * t489 + t488 * Ifges(3,5) + mrSges(3,1) * g(1) + mrSges(3,3) * t459 + t496 * t403 + t500 * t402 - pkin(2) * (-m(4) * g(1) + t423) + pkin(7) * t513;
t394 = -mrSges(2,2) * g(1) - mrSges(2,3) * t477 + Ifges(2,5) * qJDD(1) - t503 * Ifges(2,6) - pkin(6) * t401 - t497 * t395 + t501 * t396;
t393 = Ifges(2,6) * qJDD(1) + t503 * Ifges(2,5) + mrSges(2,3) * t476 + t497 * t396 + t501 * t395 - pkin(1) * t423 + pkin(6) * t514 + (-pkin(1) * t519 + mrSges(2,1)) * g(1);
t1 = [(-m(1) - m(2) + t519) * g(1) + t423; -m(1) * g(2) + t510; -m(1) * g(3) + t515; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t504; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t515 - t502 * t393 - t498 * t394; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t510 - t498 * t393 + t502 * t394; t504; t505; t506; t423; t520;];
tauJB = t1;
