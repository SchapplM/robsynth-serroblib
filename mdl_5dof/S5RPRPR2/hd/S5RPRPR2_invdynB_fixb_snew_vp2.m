% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRPR2
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
% Datum: 2020-01-03 11:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRPR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:33:53
% EndTime: 2020-01-03 11:33:56
% DurationCPUTime: 3.04s
% Computational Cost: add. (37791->206), mult. (54393->259), div. (0->0), fcn. (31463->10), ass. (0->95)
t432 = qJD(1) + qJD(3);
t428 = t432 ^ 2;
t436 = cos(pkin(9));
t468 = pkin(4) * t436;
t434 = sin(pkin(9));
t467 = mrSges(5,2) * t434;
t431 = t436 ^ 2;
t466 = t428 * t431;
t429 = qJDD(1) + qJDD(3);
t465 = t429 * t436;
t451 = Ifges(5,5) * t434 + Ifges(5,6) * t436;
t464 = t428 * t451;
t440 = sin(qJ(1));
t443 = cos(qJ(1));
t418 = -t440 * g(2) + t443 * g(3);
t444 = qJD(1) ^ 2;
t419 = -t443 * g(2) - t440 * g(3);
t414 = qJDD(1) * pkin(1) + t419;
t415 = -t444 * pkin(1) + t418;
t435 = sin(pkin(8));
t437 = cos(pkin(8));
t403 = t437 * t414 - t435 * t415;
t401 = qJDD(1) * pkin(2) + t403;
t404 = t435 * t414 + t437 * t415;
t402 = -t444 * pkin(2) + t404;
t439 = sin(qJ(3));
t442 = cos(qJ(3));
t389 = t439 * t401 + t442 * t402;
t387 = -t428 * pkin(3) + t429 * qJ(4) + t389;
t433 = -g(1) + qJDD(2);
t461 = qJD(4) * t432;
t462 = t436 * t433 - 0.2e1 * t434 * t461;
t380 = (-pkin(7) * t429 + t428 * t468 - t387) * t434 + t462;
t384 = t434 * t433 + (t387 + 0.2e1 * t461) * t436;
t381 = -pkin(4) * t466 + pkin(7) * t465 + t384;
t438 = sin(qJ(5));
t441 = cos(qJ(5));
t378 = t441 * t380 - t438 * t381;
t447 = -t434 * t438 + t436 * t441;
t407 = t447 * t432;
t448 = t434 * t441 + t436 * t438;
t408 = t448 * t432;
t395 = -t407 * mrSges(6,1) + t408 * mrSges(6,2);
t397 = t407 * qJD(5) + t448 * t429;
t405 = -qJD(5) * mrSges(6,2) + t407 * mrSges(6,3);
t376 = m(6) * t378 + qJDD(5) * mrSges(6,1) - t397 * mrSges(6,3) + qJD(5) * t405 - t408 * t395;
t379 = t438 * t380 + t441 * t381;
t396 = -t408 * qJD(5) + t447 * t429;
t406 = qJD(5) * mrSges(6,1) - t408 * mrSges(6,3);
t377 = m(6) * t379 - qJDD(5) * mrSges(6,2) + t396 * mrSges(6,3) - qJD(5) * t406 + t407 * t395;
t368 = t441 * t376 + t438 * t377;
t383 = -t434 * t387 + t462;
t449 = mrSges(5,3) * t429 + (-mrSges(5,1) * t436 + t467) * t428;
t366 = m(5) * t383 - t449 * t434 + t368;
t455 = -t438 * t376 + t441 * t377;
t367 = m(5) * t384 + t449 * t436 + t455;
t456 = -t434 * t366 + t436 * t367;
t361 = m(4) * t389 - t428 * mrSges(4,1) - t429 * mrSges(4,2) + t456;
t388 = t442 * t401 - t439 * t402;
t450 = qJDD(4) - t388;
t386 = -t429 * pkin(3) - t428 * qJ(4) + t450;
t430 = t434 ^ 2;
t382 = (-pkin(3) - t468) * t429 + (-qJ(4) + (-t430 - t431) * pkin(7)) * t428 + t450;
t446 = m(6) * t382 - t396 * mrSges(6,1) + t397 * mrSges(6,2) - t407 * t405 + t408 * t406;
t445 = -m(5) * t386 + mrSges(5,1) * t465 - t446 + (t428 * t430 + t466) * mrSges(5,3);
t372 = m(4) * t388 - t428 * mrSges(4,2) + (mrSges(4,1) - t467) * t429 + t445;
t357 = t439 * t361 + t442 * t372;
t354 = m(3) * t403 + qJDD(1) * mrSges(3,1) - t444 * mrSges(3,2) + t357;
t457 = t442 * t361 - t439 * t372;
t355 = m(3) * t404 - t444 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t457;
t458 = -t435 * t354 + t437 * t355;
t347 = m(2) * t418 - t444 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t458;
t349 = t437 * t354 + t435 * t355;
t348 = m(2) * t419 + qJDD(1) * mrSges(2,1) - t444 * mrSges(2,2) + t349;
t463 = t440 * t347 + t443 * t348;
t362 = t436 * t366 + t434 * t367;
t460 = m(4) * t433 + t362;
t459 = -t443 * t347 + t440 * t348;
t454 = m(3) * t433 + t460;
t453 = Ifges(5,1) * t434 + Ifges(5,4) * t436;
t452 = Ifges(5,4) * t434 + Ifges(5,2) * t436;
t392 = Ifges(6,1) * t408 + Ifges(6,4) * t407 + Ifges(6,5) * qJD(5);
t391 = Ifges(6,4) * t408 + Ifges(6,2) * t407 + Ifges(6,6) * qJD(5);
t390 = Ifges(6,5) * t408 + Ifges(6,6) * t407 + Ifges(6,3) * qJD(5);
t370 = mrSges(6,2) * t382 - mrSges(6,3) * t378 + Ifges(6,1) * t397 + Ifges(6,4) * t396 + Ifges(6,5) * qJDD(5) - qJD(5) * t391 + t407 * t390;
t369 = -mrSges(6,1) * t382 + mrSges(6,3) * t379 + Ifges(6,4) * t397 + Ifges(6,2) * t396 + Ifges(6,6) * qJDD(5) + qJD(5) * t392 - t408 * t390;
t358 = mrSges(5,2) * t386 - mrSges(5,3) * t383 - pkin(7) * t368 - t438 * t369 + t441 * t370 + t453 * t429 + t436 * t464;
t356 = -mrSges(5,1) * t386 + mrSges(5,3) * t384 - pkin(4) * t446 + pkin(7) * t455 + t441 * t369 + t438 * t370 + t452 * t429 - t434 * t464;
t350 = -mrSges(4,1) * t433 - mrSges(5,1) * t383 - mrSges(6,1) * t378 + mrSges(5,2) * t384 + mrSges(6,2) * t379 + mrSges(4,3) * t389 - Ifges(6,5) * t397 - Ifges(6,6) * t396 - Ifges(6,3) * qJDD(5) - pkin(3) * t362 - pkin(4) * t368 - t408 * t391 + t407 * t392 + (Ifges(4,6) - t451) * t429 + (-t434 * t452 + t436 * t453 + Ifges(4,5)) * t428;
t343 = mrSges(4,2) * t433 - mrSges(4,3) * t388 + Ifges(4,5) * t429 - t428 * Ifges(4,6) - qJ(4) * t362 - t434 * t356 + t436 * t358;
t342 = mrSges(3,2) * t433 - mrSges(3,3) * t403 + Ifges(3,5) * qJDD(1) - t444 * Ifges(3,6) - pkin(6) * t357 + t442 * t343 - t439 * t350;
t341 = -mrSges(3,1) * t433 + mrSges(3,3) * t404 + t444 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t460 + pkin(6) * t457 + t439 * t343 + t442 * t350;
t340 = -mrSges(2,2) * g(1) - mrSges(2,3) * t419 + Ifges(2,5) * qJDD(1) - t444 * Ifges(2,6) - qJ(2) * t349 - t435 * t341 + t437 * t342;
t339 = mrSges(2,1) * g(1) + mrSges(2,3) * t418 + t444 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t454 + qJ(2) * t458 + t437 * t341 + t435 * t342;
t1 = [(-m(1) - m(2)) * g(1) + t454; -m(1) * g(2) + t463; -m(1) * g(3) + t459; pkin(1) * t349 + pkin(2) * t357 + mrSges(3,1) * t403 - mrSges(3,2) * t404 + qJ(4) * t456 + t434 * t358 + t436 * t356 + pkin(3) * (-t429 * t467 + t445) + mrSges(4,1) * t388 - mrSges(4,2) * t389 + mrSges(2,1) * t419 - mrSges(2,2) * t418 + Ifges(4,3) * t429 - mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1); mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t459 + t443 * t339 + t440 * t340; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t463 + t440 * t339 - t443 * t340;];
tauB = t1;
