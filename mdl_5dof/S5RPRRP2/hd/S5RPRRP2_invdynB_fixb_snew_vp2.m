% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRRP2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:45:16
% EndTime: 2020-01-03 11:45:18
% DurationCPUTime: 1.75s
% Computational Cost: add. (19107->209), mult. (26070->250), div. (0->0), fcn. (12664->8), ass. (0->91)
t475 = Ifges(5,1) + Ifges(6,1);
t468 = Ifges(5,4) + Ifges(6,4);
t467 = Ifges(5,5) + Ifges(6,5);
t474 = Ifges(5,2) + Ifges(6,2);
t473 = Ifges(5,6) + Ifges(6,6);
t472 = Ifges(5,3) + Ifges(6,3);
t432 = qJD(1) + qJD(3);
t430 = t432 ^ 2;
t440 = sin(qJ(1));
t443 = cos(qJ(1));
t424 = -t443 * g(2) - t440 * g(3);
t416 = qJDD(1) * pkin(1) + t424;
t423 = -t440 * g(2) + t443 * g(3);
t444 = qJD(1) ^ 2;
t417 = -t444 * pkin(1) + t423;
t436 = sin(pkin(8));
t437 = cos(pkin(8));
t394 = t437 * t416 - t436 * t417;
t392 = qJDD(1) * pkin(2) + t394;
t395 = t436 * t416 + t437 * t417;
t393 = -t444 * pkin(2) + t395;
t439 = sin(qJ(3));
t442 = cos(qJ(3));
t387 = t442 * t392 - t439 * t393;
t431 = qJDD(1) + qJDD(3);
t446 = -t431 * pkin(3) - t387;
t385 = -t430 * pkin(7) + t446;
t438 = sin(qJ(4));
t441 = cos(qJ(4));
t458 = qJD(4) * t432;
t453 = t441 * t458;
t410 = t438 * t431 + t453;
t411 = t441 * t431 - t438 * t458;
t464 = t432 * t441;
t422 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t464;
t465 = t432 * t438;
t418 = qJD(4) * pkin(4) - qJ(5) * t465;
t434 = t441 ^ 2;
t381 = t418 * t465 - t411 * pkin(4) + qJDD(5) + (-qJ(5) * t434 - pkin(7)) * t430 + t446;
t421 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t464;
t447 = m(6) * t381 - t411 * mrSges(6,1) - t421 * t464;
t469 = -mrSges(5,2) - mrSges(6,2);
t471 = -m(5) * t385 + t411 * mrSges(5,1) + t410 * t469 + t422 * t464 - t447;
t470 = pkin(4) * t430;
t388 = t439 * t392 + t442 * t393;
t386 = -t430 * pkin(3) + t431 * pkin(7) + t388;
t435 = -g(1) + qJDD(2);
t383 = t441 * t386 + t438 * t435;
t409 = (-mrSges(5,1) * t441 + mrSges(5,2) * t438) * t432;
t457 = qJD(5) * t432;
t380 = t411 * qJ(5) - qJD(4) * t418 - t434 * t470 + 0.2e1 * t441 * t457 + t383;
t408 = (-mrSges(6,1) * t441 + mrSges(6,2) * t438) * t432;
t455 = m(6) * t380 + t411 * mrSges(6,3) + t408 * t464;
t419 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t465;
t459 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t465 - t419;
t375 = m(5) * t383 + t411 * mrSges(5,3) + qJD(4) * t459 + qJDD(4) * t469 + t409 * t464 + t455;
t373 = t441 * t375;
t426 = t441 * t435;
t382 = -t438 * t386 + t426;
t379 = qJDD(4) * pkin(4) + t426 + (-t410 + t453) * qJ(5) + (t441 * t470 - t386 - 0.2e1 * t457) * t438;
t456 = m(6) * t379 + qJDD(4) * mrSges(6,1) + qJD(4) * t421;
t374 = m(5) * t382 + qJDD(4) * mrSges(5,1) + qJD(4) * t422 + (-t408 - t409) * t465 + (-mrSges(5,3) - mrSges(6,3)) * t410 + t456;
t367 = m(4) * t388 - t430 * mrSges(4,1) - t431 * mrSges(4,2) - t438 * t374 + t373;
t452 = t432 * t459;
t370 = m(4) * t387 + t431 * mrSges(4,1) - t430 * mrSges(4,2) + t438 * t452 + t471;
t362 = t439 * t367 + t442 * t370;
t360 = m(3) * t394 + qJDD(1) * mrSges(3,1) - t444 * mrSges(3,2) + t362;
t449 = t442 * t367 - t439 * t370;
t361 = m(3) * t395 - t444 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t449;
t450 = -t436 * t360 + t437 * t361;
t353 = m(2) * t423 - t444 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t450;
t355 = t437 * t360 + t436 * t361;
t354 = m(2) * t424 + qJDD(1) * mrSges(2,1) - t444 * mrSges(2,2) + t355;
t463 = t440 * t353 + t443 * t354;
t368 = t441 * t374 + t438 * t375;
t462 = (t467 * t438 + t473 * t441) * t432 + t472 * qJD(4);
t461 = (-t468 * t438 - t474 * t441) * t432 - t473 * qJD(4);
t460 = (t475 * t438 + t468 * t441) * t432 + t467 * qJD(4);
t454 = m(4) * t435 + t368;
t451 = -t443 * t353 + t440 * t354;
t448 = m(3) * t435 + t454;
t376 = -t410 * mrSges(6,3) - t408 * t465 + t456;
t364 = mrSges(5,2) * t385 + mrSges(6,2) * t381 - mrSges(5,3) * t382 - mrSges(6,3) * t379 - qJ(5) * t376 + t461 * qJD(4) + t467 * qJDD(4) + t475 * t410 + t468 * t411 + t462 * t464;
t363 = -mrSges(5,1) * t385 + mrSges(5,3) * t383 - mrSges(6,1) * t381 + mrSges(6,3) * t380 - pkin(4) * t447 + qJ(5) * t455 + (-pkin(4) * t419 - t462) * t465 + t474 * t411 + (-pkin(4) * mrSges(6,2) + t468) * t410 + (-qJ(5) * mrSges(6,2) + t473) * qJDD(4) + (-qJ(5) * t419 + t460) * qJD(4);
t356 = -mrSges(4,1) * t435 - mrSges(5,1) * t382 - mrSges(6,1) * t379 + mrSges(5,2) * t383 + mrSges(6,2) * t380 + mrSges(4,3) * t388 + t430 * Ifges(4,5) + Ifges(4,6) * t431 - pkin(3) * t368 - pkin(4) * t376 - t473 * t411 - t467 * t410 - t472 * qJDD(4) + (t438 * t461 + t460 * t441) * t432;
t349 = mrSges(4,2) * t435 - mrSges(4,3) * t387 + Ifges(4,5) * t431 - t430 * Ifges(4,6) - pkin(7) * t368 - t438 * t363 + t441 * t364;
t348 = mrSges(3,2) * t435 - mrSges(3,3) * t394 + Ifges(3,5) * qJDD(1) - t444 * Ifges(3,6) - pkin(6) * t362 + t442 * t349 - t439 * t356;
t347 = -mrSges(3,1) * t435 + mrSges(3,3) * t395 + t444 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t454 + pkin(6) * t449 + t439 * t349 + t442 * t356;
t346 = -mrSges(2,2) * g(1) - mrSges(2,3) * t424 + Ifges(2,5) * qJDD(1) - t444 * Ifges(2,6) - qJ(2) * t355 - t436 * t347 + t437 * t348;
t345 = mrSges(2,1) * g(1) + mrSges(2,3) * t423 + t444 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t448 + qJ(2) * t450 + t437 * t347 + t436 * t348;
t1 = [(-m(1) - m(2)) * g(1) + t448; -m(1) * g(2) + t463; -m(1) * g(3) + t451; pkin(1) * t355 + pkin(2) * t362 + mrSges(3,1) * t394 - mrSges(3,2) * t395 + t441 * t363 + pkin(3) * t471 + pkin(7) * t373 + mrSges(4,1) * t387 - mrSges(4,2) * t388 + mrSges(2,1) * t424 - mrSges(2,2) * t423 + Ifges(4,3) * t431 - mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + (pkin(3) * t452 - pkin(7) * t374 + t364) * t438 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1); mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t451 + t443 * t345 + t440 * t346; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t463 + t440 * t345 - t443 * t346;];
tauB = t1;
