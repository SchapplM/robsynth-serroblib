% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2020-01-03 11:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPPR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR1_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR1_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR1_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:55:53
% EndTime: 2020-01-03 11:55:55
% DurationCPUTime: 2.45s
% Computational Cost: add. (39959->207), mult. (54393->259), div. (0->0), fcn. (31463->10), ass. (0->94)
t432 = qJD(1) + qJD(2);
t428 = t432 ^ 2;
t436 = cos(pkin(9));
t467 = pkin(4) * t436;
t434 = sin(pkin(9));
t466 = mrSges(5,2) * t434;
t451 = Ifges(5,5) * t434 + Ifges(5,6) * t436;
t465 = t451 * t428;
t431 = t436 ^ 2;
t464 = t428 * t431;
t429 = qJDD(1) + qJDD(2);
t463 = t429 * t436;
t440 = sin(qJ(1));
t443 = cos(qJ(1));
t419 = -g(2) * t440 + t443 * g(3);
t444 = qJD(1) ^ 2;
t420 = -g(2) * t443 - g(3) * t440;
t415 = qJDD(1) * pkin(1) + t420;
t416 = -pkin(1) * t444 + t419;
t439 = sin(qJ(2));
t442 = cos(qJ(2));
t404 = t442 * t415 - t416 * t439;
t402 = pkin(2) * t429 + t404;
t405 = t439 * t415 + t442 * t416;
t403 = -pkin(2) * t428 + t405;
t435 = sin(pkin(8));
t437 = cos(pkin(8));
t390 = t435 * t402 + t437 * t403;
t388 = -pkin(3) * t428 + qJ(4) * t429 + t390;
t433 = -g(1) + qJDD(3);
t460 = qJD(4) * t432;
t461 = t436 * t433 - 0.2e1 * t434 * t460;
t381 = (-pkin(7) * t429 + t428 * t467 - t388) * t434 + t461;
t385 = t434 * t433 + (t388 + 0.2e1 * t460) * t436;
t382 = -pkin(4) * t464 + pkin(7) * t463 + t385;
t438 = sin(qJ(5));
t441 = cos(qJ(5));
t379 = t381 * t441 - t382 * t438;
t447 = -t434 * t438 + t436 * t441;
t408 = t447 * t432;
t448 = t434 * t441 + t436 * t438;
t409 = t448 * t432;
t396 = -mrSges(6,1) * t408 + mrSges(6,2) * t409;
t398 = t408 * qJD(5) + t448 * t429;
t406 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t408;
t377 = m(6) * t379 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t398 + qJD(5) * t406 - t396 * t409;
t380 = t381 * t438 + t382 * t441;
t397 = -t409 * qJD(5) + t447 * t429;
t407 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t409;
t378 = m(6) * t380 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t397 - qJD(5) * t407 + t396 * t408;
t369 = t441 * t377 + t438 * t378;
t384 = -t388 * t434 + t461;
t449 = mrSges(5,3) * t429 + (-mrSges(5,1) * t436 + t466) * t428;
t367 = m(5) * t384 - t449 * t434 + t369;
t454 = -t438 * t377 + t441 * t378;
t368 = m(5) * t385 + t449 * t436 + t454;
t455 = -t367 * t434 + t436 * t368;
t362 = m(4) * t390 - mrSges(4,1) * t428 - mrSges(4,2) * t429 + t455;
t389 = t437 * t402 - t435 * t403;
t450 = qJDD(4) - t389;
t387 = -t429 * pkin(3) - t428 * qJ(4) + t450;
t430 = t434 ^ 2;
t383 = (-pkin(3) - t467) * t429 + (-qJ(4) + (-t430 - t431) * pkin(7)) * t428 + t450;
t446 = m(6) * t383 - t397 * mrSges(6,1) + t398 * mrSges(6,2) - t408 * t406 + t409 * t407;
t445 = -m(5) * t387 + mrSges(5,1) * t463 - t446 + (t428 * t430 + t464) * mrSges(5,3);
t373 = m(4) * t389 - t428 * mrSges(4,2) + (mrSges(4,1) - t466) * t429 + t445;
t358 = t435 * t362 + t437 * t373;
t355 = m(3) * t404 + mrSges(3,1) * t429 - mrSges(3,2) * t428 + t358;
t456 = t437 * t362 - t373 * t435;
t356 = m(3) * t405 - mrSges(3,1) * t428 - mrSges(3,2) * t429 + t456;
t457 = -t439 * t355 + t442 * t356;
t348 = m(2) * t419 - mrSges(2,1) * t444 - qJDD(1) * mrSges(2,2) + t457;
t350 = t442 * t355 + t439 * t356;
t349 = m(2) * t420 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t444 + t350;
t462 = t440 * t348 + t443 * t349;
t363 = t436 * t367 + t434 * t368;
t459 = m(4) * t433 + t363;
t458 = -t348 * t443 + t440 * t349;
t453 = Ifges(5,1) * t434 + Ifges(5,4) * t436;
t452 = Ifges(5,4) * t434 + Ifges(5,2) * t436;
t393 = Ifges(6,1) * t409 + Ifges(6,4) * t408 + Ifges(6,5) * qJD(5);
t392 = Ifges(6,4) * t409 + Ifges(6,2) * t408 + Ifges(6,6) * qJD(5);
t391 = Ifges(6,5) * t409 + Ifges(6,6) * t408 + Ifges(6,3) * qJD(5);
t371 = mrSges(6,2) * t383 - mrSges(6,3) * t379 + Ifges(6,1) * t398 + Ifges(6,4) * t397 + Ifges(6,5) * qJDD(5) - qJD(5) * t392 + t391 * t408;
t370 = -mrSges(6,1) * t383 + mrSges(6,3) * t380 + Ifges(6,4) * t398 + Ifges(6,2) * t397 + Ifges(6,6) * qJDD(5) + qJD(5) * t393 - t391 * t409;
t359 = mrSges(5,2) * t387 - mrSges(5,3) * t384 - pkin(7) * t369 - t438 * t370 + t441 * t371 + t453 * t429 + t436 * t465;
t357 = -mrSges(5,1) * t387 + mrSges(5,3) * t385 - pkin(4) * t446 + pkin(7) * t454 + t441 * t370 + t438 * t371 + t452 * t429 - t434 * t465;
t351 = -mrSges(4,1) * t433 - mrSges(5,1) * t384 - mrSges(6,1) * t379 + mrSges(5,2) * t385 + mrSges(6,2) * t380 + mrSges(4,3) * t390 - Ifges(6,5) * t398 - Ifges(6,6) * t397 - Ifges(6,3) * qJDD(5) - pkin(3) * t363 - pkin(4) * t369 - t409 * t392 + t408 * t393 + (Ifges(4,6) - t451) * t429 + (-t434 * t452 + t436 * t453 + Ifges(4,5)) * t428;
t344 = mrSges(4,2) * t433 - mrSges(4,3) * t389 + Ifges(4,5) * t429 - Ifges(4,6) * t428 - qJ(4) * t363 - t357 * t434 + t359 * t436;
t343 = -mrSges(3,2) * g(1) - mrSges(3,3) * t404 + Ifges(3,5) * t429 - Ifges(3,6) * t428 - qJ(3) * t358 + t344 * t437 - t351 * t435;
t342 = mrSges(3,1) * g(1) + mrSges(3,3) * t405 + t428 * Ifges(3,5) + Ifges(3,6) * t429 - pkin(2) * t459 + qJ(3) * t456 + t435 * t344 + t437 * t351;
t341 = -mrSges(2,2) * g(1) - mrSges(2,3) * t420 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t444 - pkin(6) * t350 - t342 * t439 + t343 * t442;
t340 = Ifges(2,6) * qJDD(1) + t444 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t419 + t439 * t343 + t442 * t342 - pkin(1) * (-m(3) * g(1) + t459) + pkin(6) * t457;
t1 = [(-m(1) - m(2) - m(3)) * g(1) + t459; -m(1) * g(2) + t462; -m(1) * g(3) + t458; pkin(1) * t350 + pkin(2) * t358 + mrSges(3,1) * t404 - mrSges(3,2) * t405 + t436 * t357 + pkin(3) * t445 + qJ(4) * t455 + t434 * t359 + mrSges(4,1) * t389 - mrSges(4,2) * t390 + mrSges(2,1) * t420 - mrSges(2,2) * t419 + Ifges(2,3) * qJDD(1) - mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + (-pkin(3) * t466 + Ifges(3,3) + Ifges(4,3)) * t429; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t458 + t443 * t340 + t440 * t341; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t462 + t440 * t340 - t443 * t341;];
tauB = t1;
