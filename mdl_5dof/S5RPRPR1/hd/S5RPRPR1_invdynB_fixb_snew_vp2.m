% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRPR1
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRPR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR1_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR1_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR1_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR1_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:21
% EndTime: 2019-12-05 17:47:24
% DurationCPUTime: 2.86s
% Computational Cost: add. (28361->264), mult. (61605->328), div. (0->0), fcn. (38710->8), ass. (0->102)
t495 = sin(qJ(1));
t498 = cos(qJ(1));
t480 = -t498 * g(1) - t495 * g(2);
t506 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t480;
t521 = -pkin(1) - pkin(6);
t520 = mrSges(2,1) - mrSges(3,2);
t519 = Ifges(2,5) - Ifges(3,4);
t518 = (-Ifges(2,6) + Ifges(3,5));
t479 = t495 * g(1) - t498 * g(2);
t499 = qJD(1) ^ 2;
t505 = -t499 * qJ(2) + qJDD(2) - t479;
t458 = t521 * qJDD(1) + t505;
t494 = sin(qJ(3));
t497 = cos(qJ(3));
t449 = t494 * g(3) + t497 * t458;
t514 = qJD(1) * qJD(3);
t512 = t494 * t514;
t475 = t497 * qJDD(1) - t512;
t433 = (-t475 - t512) * qJ(4) + (-t494 * t497 * t499 + qJDD(3)) * pkin(3) + t449;
t450 = -t497 * g(3) + t494 * t458;
t474 = -t494 * qJDD(1) - t497 * t514;
t515 = qJD(1) * t497;
t477 = qJD(3) * pkin(3) - qJ(4) * t515;
t488 = t494 ^ 2;
t434 = -t488 * t499 * pkin(3) + t474 * qJ(4) - qJD(3) * t477 + t450;
t491 = sin(pkin(8));
t492 = cos(pkin(8));
t465 = (-t491 * t494 + t492 * t497) * qJD(1);
t419 = -0.2e1 * qJD(4) * t465 + t492 * t433 - t491 * t434;
t448 = t491 * t474 + t492 * t475;
t464 = (-t491 * t497 - t492 * t494) * qJD(1);
t416 = (qJD(3) * t464 - t448) * pkin(7) + (t464 * t465 + qJDD(3)) * pkin(4) + t419;
t420 = 0.2e1 * qJD(4) * t464 + t491 * t433 + t492 * t434;
t447 = t492 * t474 - t491 * t475;
t457 = qJD(3) * pkin(4) - t465 * pkin(7);
t463 = t464 ^ 2;
t417 = -t463 * pkin(4) + t447 * pkin(7) - qJD(3) * t457 + t420;
t493 = sin(qJ(5));
t496 = cos(qJ(5));
t414 = t496 * t416 - t493 * t417;
t442 = t496 * t464 - t493 * t465;
t424 = t442 * qJD(5) + t493 * t447 + t496 * t448;
t443 = t493 * t464 + t496 * t465;
t429 = -t442 * mrSges(6,1) + t443 * mrSges(6,2);
t485 = qJD(3) + qJD(5);
t437 = -t485 * mrSges(6,2) + t442 * mrSges(6,3);
t484 = qJDD(3) + qJDD(5);
t412 = m(6) * t414 + t484 * mrSges(6,1) - t424 * mrSges(6,3) - t443 * t429 + t485 * t437;
t415 = t493 * t416 + t496 * t417;
t423 = -t443 * qJD(5) + t496 * t447 - t493 * t448;
t438 = t485 * mrSges(6,1) - t443 * mrSges(6,3);
t413 = m(6) * t415 - t484 * mrSges(6,2) + t423 * mrSges(6,3) + t442 * t429 - t485 * t438;
t403 = t496 * t412 + t493 * t413;
t445 = -t464 * mrSges(5,1) + t465 * mrSges(5,2);
t455 = -qJD(3) * mrSges(5,2) + t464 * mrSges(5,3);
t401 = m(5) * t419 + qJDD(3) * mrSges(5,1) - t448 * mrSges(5,3) + qJD(3) * t455 - t465 * t445 + t403;
t456 = qJD(3) * mrSges(5,1) - t465 * mrSges(5,3);
t508 = -t493 * t412 + t496 * t413;
t402 = m(5) * t420 - qJDD(3) * mrSges(5,2) + t447 * mrSges(5,3) - qJD(3) * t456 + t464 * t445 + t508;
t397 = t492 * t401 + t491 * t402;
t473 = (mrSges(4,1) * t494 + mrSges(4,2) * t497) * qJD(1);
t516 = qJD(1) * t494;
t476 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t516;
t395 = m(4) * t449 + qJDD(3) * mrSges(4,1) - t475 * mrSges(4,3) + qJD(3) * t476 - t473 * t515 + t397;
t478 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t515;
t509 = -t491 * t401 + t492 * t402;
t396 = m(4) * t450 - qJDD(3) * mrSges(4,2) + t474 * mrSges(4,3) - qJD(3) * t478 - t473 * t516 + t509;
t391 = t497 * t395 + t494 * t396;
t462 = -qJDD(1) * pkin(1) + t505;
t503 = -m(3) * t462 + (t499 * mrSges(3,3)) - t391;
t389 = m(2) * t479 - (t499 * mrSges(2,2)) + t520 * qJDD(1) + t503;
t459 = t499 * pkin(1) - t506;
t454 = t521 * t499 + t506;
t436 = -t474 * pkin(3) + qJDD(4) + t477 * t515 + (-qJ(4) * t488 + t521) * t499 + t506;
t421 = -t447 * pkin(4) - t463 * pkin(7) + t465 * t457 + t436;
t504 = m(6) * t421 - t423 * mrSges(6,1) + t424 * mrSges(6,2) - t442 * t437 + t443 * t438;
t502 = m(5) * t436 - t447 * mrSges(5,1) + t448 * mrSges(5,2) - t464 * t455 + t465 * t456 + t504;
t501 = -m(4) * t454 + t474 * mrSges(4,1) - t475 * mrSges(4,2) - t476 * t516 - t478 * t515 - t502;
t500 = -m(3) * t459 + (t499 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t501;
t408 = m(2) * t480 - (t499 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t500;
t517 = t498 * t389 + t495 * t408;
t511 = -t495 * t389 + t498 * t408;
t510 = -t494 * t395 + t497 * t396;
t468 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t497 - Ifges(4,4) * t494) * qJD(1);
t467 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t497 - Ifges(4,2) * t494) * qJD(1);
t466 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t497 - Ifges(4,6) * t494) * qJD(1);
t441 = Ifges(5,1) * t465 + Ifges(5,4) * t464 + (Ifges(5,5) * qJD(3));
t440 = Ifges(5,4) * t465 + Ifges(5,2) * t464 + (Ifges(5,6) * qJD(3));
t439 = Ifges(5,5) * t465 + Ifges(5,6) * t464 + (Ifges(5,3) * qJD(3));
t427 = Ifges(6,1) * t443 + Ifges(6,4) * t442 + Ifges(6,5) * t485;
t426 = Ifges(6,4) * t443 + Ifges(6,2) * t442 + Ifges(6,6) * t485;
t425 = Ifges(6,5) * t443 + Ifges(6,6) * t442 + Ifges(6,3) * t485;
t405 = mrSges(6,2) * t421 - mrSges(6,3) * t414 + Ifges(6,1) * t424 + Ifges(6,4) * t423 + Ifges(6,5) * t484 + t442 * t425 - t485 * t426;
t404 = -mrSges(6,1) * t421 + mrSges(6,3) * t415 + Ifges(6,4) * t424 + Ifges(6,2) * t423 + Ifges(6,6) * t484 - t443 * t425 + t485 * t427;
t393 = mrSges(5,2) * t436 - mrSges(5,3) * t419 + Ifges(5,1) * t448 + Ifges(5,4) * t447 + Ifges(5,5) * qJDD(3) - pkin(7) * t403 - qJD(3) * t440 - t493 * t404 + t496 * t405 + t464 * t439;
t392 = -mrSges(5,1) * t436 + mrSges(5,3) * t420 + Ifges(5,4) * t448 + Ifges(5,2) * t447 + Ifges(5,6) * qJDD(3) - pkin(4) * t504 + pkin(7) * t508 + qJD(3) * t441 + t496 * t404 + t493 * t405 - t465 * t439;
t390 = -m(3) * g(3) + t510;
t387 = mrSges(4,2) * t454 - mrSges(4,3) * t449 + Ifges(4,1) * t475 + Ifges(4,4) * t474 + Ifges(4,5) * qJDD(3) - qJ(4) * t397 - qJD(3) * t467 - t491 * t392 + t492 * t393 - t466 * t516;
t386 = -mrSges(4,1) * t454 + mrSges(4,3) * t450 + Ifges(4,4) * t475 + Ifges(4,2) * t474 + Ifges(4,6) * qJDD(3) - pkin(3) * t502 + qJ(4) * t509 + qJD(3) * t468 + t492 * t392 + t491 * t393 - t466 * t515;
t385 = -t442 * t427 + t443 * t426 + Ifges(5,6) * t447 + Ifges(5,5) * t448 + pkin(2) * t391 - qJ(2) * t390 + (t518 * t499) + t519 * qJDD(1) + Ifges(4,5) * t475 - mrSges(2,3) * t479 + Ifges(6,3) * t484 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + mrSges(4,1) * t449 - mrSges(4,2) * t450 + mrSges(3,1) * t462 + pkin(4) * t403 - t464 * t441 + t465 * t440 + Ifges(4,6) * t474 + Ifges(6,6) * t423 + Ifges(6,5) * t424 + pkin(3) * t397 + mrSges(5,1) * t419 - mrSges(5,2) * t420 + mrSges(6,1) * t414 - mrSges(6,2) * t415 + (t497 * t467 + t494 * t468) * qJD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t384 = -mrSges(3,1) * t459 + mrSges(2,3) * t480 - pkin(1) * t390 - pkin(2) * t501 - pkin(6) * t510 + t520 * g(3) - t518 * qJDD(1) - t497 * t386 - t494 * t387 + t519 * t499;
t1 = [-m(1) * g(1) + t511; -m(1) * g(2) + t517; (-m(1) - m(2) - m(3)) * g(3) + t510; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t517 - t495 * t384 + t498 * t385; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t511 + t498 * t384 + t495 * t385; pkin(1) * t503 + qJ(2) * t500 - t494 * t386 - pkin(6) * t391 + mrSges(2,1) * t479 - mrSges(2,2) * t480 + t497 * t387 + mrSges(3,2) * t462 - mrSges(3,3) * t459 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB = t1;
