% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRPR15
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
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRPR15_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR15_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR15_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR15_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR15_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR15_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR15_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:36:50
% EndTime: 2019-12-31 18:36:54
% DurationCPUTime: 2.32s
% Computational Cost: add. (23410->263), mult. (48358->323), div. (0->0), fcn. (28723->8), ass. (0->103)
t476 = sin(qJ(1));
t479 = cos(qJ(1));
t464 = -t479 * g(1) - t476 * g(2);
t504 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t464;
t503 = -pkin(1) - pkin(6);
t502 = mrSges(2,1) - mrSges(3,2);
t501 = -Ifges(3,4) + Ifges(2,5);
t500 = (Ifges(3,5) - Ifges(2,6));
t463 = t476 * g(1) - t479 * g(2);
t481 = qJD(1) ^ 2;
t487 = -t481 * qJ(2) + qJDD(2) - t463;
t441 = t503 * qJDD(1) + t487;
t475 = sin(qJ(3));
t478 = cos(qJ(3));
t431 = -t478 * g(3) + t475 * t441;
t458 = (mrSges(4,1) * t475 + mrSges(4,2) * t478) * qJD(1);
t496 = qJD(1) * qJD(3);
t493 = t478 * t496;
t459 = t475 * qJDD(1) + t493;
t498 = qJD(1) * t478;
t462 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t498;
t435 = t503 * t481 - t504;
t494 = t475 * t496;
t460 = t478 * qJDD(1) - t494;
t417 = (-t460 + t494) * qJ(4) + (t459 + t493) * pkin(3) + t435;
t457 = (pkin(3) * t475 - qJ(4) * t478) * qJD(1);
t480 = qJD(3) ^ 2;
t497 = t475 * qJD(1);
t420 = -t480 * pkin(3) + qJDD(3) * qJ(4) - t457 * t497 + t431;
t472 = sin(pkin(8));
t473 = cos(pkin(8));
t454 = t472 * qJD(3) + t473 * t498;
t404 = -0.2e1 * qJD(4) * t454 + t473 * t417 - t472 * t420;
t439 = t472 * qJDD(3) + t473 * t460;
t453 = t473 * qJD(3) - t472 * t498;
t402 = (t453 * t497 - t439) * pkin(7) + (t453 * t454 + t459) * pkin(4) + t404;
t405 = 0.2e1 * qJD(4) * t453 + t472 * t417 + t473 * t420;
t438 = t473 * qJDD(3) - t472 * t460;
t440 = pkin(4) * t497 - t454 * pkin(7);
t452 = t453 ^ 2;
t403 = -t452 * pkin(4) + t438 * pkin(7) - t440 * t497 + t405;
t474 = sin(qJ(5));
t477 = cos(qJ(5));
t400 = t477 * t402 - t474 * t403;
t427 = t477 * t453 - t474 * t454;
t409 = t427 * qJD(5) + t474 * t438 + t477 * t439;
t428 = t474 * t453 + t477 * t454;
t414 = -t427 * mrSges(6,1) + t428 * mrSges(6,2);
t465 = qJD(5) + t497;
t421 = -t465 * mrSges(6,2) + t427 * mrSges(6,3);
t456 = qJDD(5) + t459;
t398 = m(6) * t400 + t456 * mrSges(6,1) - t409 * mrSges(6,3) - t428 * t414 + t465 * t421;
t401 = t474 * t402 + t477 * t403;
t408 = -t428 * qJD(5) + t477 * t438 - t474 * t439;
t422 = t465 * mrSges(6,1) - t428 * mrSges(6,3);
t399 = m(6) * t401 - t456 * mrSges(6,2) + t408 * mrSges(6,3) + t427 * t414 - t465 * t422;
t391 = t477 * t398 + t474 * t399;
t429 = -t453 * mrSges(5,1) + t454 * mrSges(5,2);
t436 = -mrSges(5,2) * t497 + t453 * mrSges(5,3);
t389 = m(5) * t404 + t459 * mrSges(5,1) - t439 * mrSges(5,3) - t454 * t429 + t436 * t497 + t391;
t437 = mrSges(5,1) * t497 - t454 * mrSges(5,3);
t489 = -t474 * t398 + t477 * t399;
t390 = m(5) * t405 - t459 * mrSges(5,2) + t438 * mrSges(5,3) + t453 * t429 - t437 * t497 + t489;
t490 = -t472 * t389 + t473 * t390;
t384 = m(4) * t431 - qJDD(3) * mrSges(4,2) - t459 * mrSges(4,3) - qJD(3) * t462 - t458 * t497 + t490;
t430 = t475 * g(3) + t478 * t441;
t461 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t497;
t419 = -qJDD(3) * pkin(3) - t480 * qJ(4) + t457 * t498 + qJDD(4) - t430;
t406 = -t438 * pkin(4) - t452 * pkin(7) + t454 * t440 + t419;
t484 = m(6) * t406 - t408 * mrSges(6,1) + t409 * mrSges(6,2) - t427 * t421 + t428 * t422;
t482 = -m(5) * t419 + t438 * mrSges(5,1) - t439 * mrSges(5,2) + t453 * t436 - t454 * t437 - t484;
t394 = m(4) * t430 + qJDD(3) * mrSges(4,1) - t460 * mrSges(4,3) + qJD(3) * t461 - t458 * t498 + t482;
t378 = t475 * t384 + t478 * t394;
t443 = -qJDD(1) * pkin(1) + t487;
t486 = -m(3) * t443 + (t481 * mrSges(3,3)) - t378;
t375 = m(2) * t463 - (t481 * mrSges(2,2)) + t502 * qJDD(1) + t486;
t442 = t481 * pkin(1) + t504;
t385 = t473 * t389 + t472 * t390;
t485 = -m(4) * t435 - t459 * mrSges(4,1) - t460 * mrSges(4,2) - t461 * t497 - t462 * t498 - t385;
t483 = -m(3) * t442 + (t481 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t485;
t382 = m(2) * t464 - (t481 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t483;
t499 = t479 * t375 + t476 * t382;
t492 = -t476 * t375 + t479 * t382;
t491 = t478 * t384 - t475 * t394;
t449 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t478 - Ifges(4,4) * t475) * qJD(1);
t448 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t478 - Ifges(4,2) * t475) * qJD(1);
t447 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t478 - Ifges(4,6) * t475) * qJD(1);
t425 = Ifges(5,1) * t454 + Ifges(5,4) * t453 + Ifges(5,5) * t497;
t424 = Ifges(5,4) * t454 + Ifges(5,2) * t453 + Ifges(5,6) * t497;
t423 = Ifges(5,5) * t454 + Ifges(5,6) * t453 + Ifges(5,3) * t497;
t412 = Ifges(6,1) * t428 + Ifges(6,4) * t427 + Ifges(6,5) * t465;
t411 = Ifges(6,4) * t428 + Ifges(6,2) * t427 + Ifges(6,6) * t465;
t410 = Ifges(6,5) * t428 + Ifges(6,6) * t427 + Ifges(6,3) * t465;
t393 = mrSges(6,2) * t406 - mrSges(6,3) * t400 + Ifges(6,1) * t409 + Ifges(6,4) * t408 + Ifges(6,5) * t456 + t427 * t410 - t465 * t411;
t392 = -mrSges(6,1) * t406 + mrSges(6,3) * t401 + Ifges(6,4) * t409 + Ifges(6,2) * t408 + Ifges(6,6) * t456 - t428 * t410 + t465 * t412;
t379 = mrSges(5,2) * t419 - mrSges(5,3) * t404 + Ifges(5,1) * t439 + Ifges(5,4) * t438 + Ifges(5,5) * t459 - pkin(7) * t391 - t474 * t392 + t477 * t393 + t453 * t423 - t424 * t497;
t377 = -m(3) * g(3) + t491;
t376 = -mrSges(5,1) * t419 + mrSges(5,3) * t405 + Ifges(5,4) * t439 + Ifges(5,2) * t438 + Ifges(5,6) * t459 - pkin(4) * t484 + pkin(7) * t489 + t477 * t392 + t474 * t393 - t454 * t423 + t425 * t497;
t373 = Ifges(4,4) * t460 + Ifges(4,6) * qJDD(3) - t447 * t498 + qJD(3) * t449 - mrSges(4,1) * t435 + mrSges(4,3) * t431 - Ifges(5,5) * t439 - Ifges(5,6) * t438 - t454 * t424 + t453 * t425 - mrSges(5,1) * t404 + mrSges(5,2) * t405 - Ifges(6,5) * t409 - Ifges(6,6) * t408 - Ifges(6,3) * t456 - t428 * t411 + t427 * t412 - mrSges(6,1) * t400 + mrSges(6,2) * t401 - pkin(4) * t391 - pkin(3) * t385 + (-Ifges(4,2) - Ifges(5,3)) * t459;
t372 = mrSges(4,2) * t435 - mrSges(4,3) * t430 + Ifges(4,1) * t460 - Ifges(4,4) * t459 + Ifges(4,5) * qJDD(3) - qJ(4) * t385 - qJD(3) * t448 - t472 * t376 + t473 * t379 - t447 * t497;
t371 = -qJ(2) * t377 - mrSges(2,3) * t463 + pkin(2) * t378 + mrSges(3,1) * t443 + t472 * t379 + t473 * t376 + pkin(3) * t482 + qJ(4) * t490 + Ifges(4,5) * t460 - Ifges(4,6) * t459 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t430 - mrSges(4,2) * t431 + (t500 * t481) + t501 * qJDD(1) + (t478 * t448 + t475 * t449) * qJD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t370 = -mrSges(3,1) * t442 + mrSges(2,3) * t464 - pkin(1) * t377 - pkin(2) * t485 - pkin(6) * t491 + t502 * g(3) - t500 * qJDD(1) - t475 * t372 - t478 * t373 + t501 * t481;
t1 = [-m(1) * g(1) + t492; -m(1) * g(2) + t499; (-m(1) - m(2) - m(3)) * g(3) + t491; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t499 - t476 * t370 + t479 * t371; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t492 + t479 * t370 + t476 * t371; pkin(1) * t486 + qJ(2) * t483 + mrSges(2,1) * t463 - mrSges(2,2) * t464 + t478 * t372 - t475 * t373 - pkin(6) * t378 + mrSges(3,2) * t443 - mrSges(3,3) * t442 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB = t1;
