% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPRR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPRR7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR7_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR7_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR7_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR7_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR7_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:36
% EndTime: 2019-12-31 20:15:37
% DurationCPUTime: 1.53s
% Computational Cost: add. (21366->224), mult. (25930->274), div. (0->0), fcn. (14043->8), ass. (0->92)
t435 = sin(qJ(1));
t439 = cos(qJ(1));
t417 = t435 * g(1) - g(2) * t439;
t412 = qJDD(1) * pkin(1) + t417;
t418 = -g(1) * t439 - g(2) * t435;
t440 = qJD(1) ^ 2;
t413 = -pkin(1) * t440 + t418;
t434 = sin(qJ(2));
t438 = cos(qJ(2));
t395 = t434 * t412 + t438 * t413;
t428 = qJDD(1) + qJDD(2);
t430 = (qJD(1) + qJD(2));
t446 = t428 * qJ(3) + (2 * qJD(3) * t430) + t395;
t462 = -m(3) - m(4);
t461 = -pkin(2) - pkin(7);
t460 = mrSges(3,1) - mrSges(4,2);
t459 = Ifges(3,5) - Ifges(4,4);
t458 = (-Ifges(3,6) + Ifges(4,5));
t433 = sin(qJ(4));
t457 = t430 * t433;
t437 = cos(qJ(4));
t456 = t430 * t437;
t394 = t412 * t438 - t434 * t413;
t426 = t430 ^ 2;
t445 = -qJ(3) * t426 + qJDD(3) - t394;
t386 = t461 * t428 + t445;
t380 = t433 * g(3) + t437 * t386;
t453 = qJD(4) * t430;
t451 = t433 * t453;
t408 = t428 * t437 - t451;
t374 = (-t408 - t451) * pkin(8) + (-t426 * t433 * t437 + qJDD(4)) * pkin(4) + t380;
t381 = -g(3) * t437 + t433 * t386;
t407 = -t428 * t433 - t437 * t453;
t416 = qJD(4) * pkin(4) - pkin(8) * t456;
t431 = t433 ^ 2;
t376 = -pkin(4) * t426 * t431 + pkin(8) * t407 - qJD(4) * t416 + t381;
t432 = sin(qJ(5));
t436 = cos(qJ(5));
t371 = t374 * t436 - t376 * t432;
t401 = (-t432 * t437 - t433 * t436) * t430;
t379 = qJD(5) * t401 + t407 * t432 + t408 * t436;
t402 = (-t432 * t433 + t436 * t437) * t430;
t393 = -mrSges(6,1) * t401 + mrSges(6,2) * t402;
t429 = qJD(4) + qJD(5);
t396 = -mrSges(6,2) * t429 + mrSges(6,3) * t401;
t427 = qJDD(4) + qJDD(5);
t369 = m(6) * t371 + mrSges(6,1) * t427 - mrSges(6,3) * t379 - t393 * t402 + t396 * t429;
t372 = t374 * t432 + t376 * t436;
t378 = -qJD(5) * t402 + t407 * t436 - t408 * t432;
t397 = mrSges(6,1) * t429 - mrSges(6,3) * t402;
t370 = m(6) * t372 - mrSges(6,2) * t427 + mrSges(6,3) * t378 + t393 * t401 - t397 * t429;
t360 = t436 * t369 + t432 * t370;
t406 = (mrSges(5,1) * t433 + mrSges(5,2) * t437) * t430;
t414 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t457;
t358 = m(5) * t380 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t408 + qJD(4) * t414 - t406 * t456 + t360;
t415 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t456;
t447 = -t369 * t432 + t436 * t370;
t359 = m(5) * t381 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t407 - qJD(4) * t415 - t406 * t457 + t447;
t356 = t358 * t437 + t359 * t433;
t392 = -pkin(2) * t428 + t445;
t443 = -m(4) * t392 + (t426 * mrSges(4,3)) - t356;
t354 = m(3) * t394 - (mrSges(3,2) * t426) + t460 * t428 + t443;
t390 = pkin(2) * t426 - t446;
t383 = t461 * t426 + t446;
t375 = t416 * t456 - pkin(4) * t407 + (-pkin(8) * t431 + t461) * t426 + t446;
t444 = m(6) * t375 - mrSges(6,1) * t378 + t379 * mrSges(6,2) - t396 * t401 + t402 * t397;
t442 = -m(5) * t383 + mrSges(5,1) * t407 - t408 * mrSges(5,2) - t414 * t457 - t415 * t456 - t444;
t441 = -m(4) * t390 + (t426 * mrSges(4,2)) + t428 * mrSges(4,3) - t442;
t365 = m(3) * t395 - (mrSges(3,1) * t426) - mrSges(3,2) * t428 + t441;
t350 = t438 * t354 + t434 * t365;
t348 = m(2) * t417 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t440 + t350;
t449 = -t434 * t354 + t438 * t365;
t349 = m(2) * t418 - mrSges(2,1) * t440 - qJDD(1) * mrSges(2,2) + t449;
t455 = t439 * t348 + t435 * t349;
t450 = -t348 * t435 + t439 * t349;
t448 = t358 * t433 - t437 * t359;
t400 = (Ifges(5,5) * qJD(4)) + (Ifges(5,1) * t437 - Ifges(5,4) * t433) * t430;
t399 = (Ifges(5,6) * qJD(4)) + (Ifges(5,4) * t437 - Ifges(5,2) * t433) * t430;
t398 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t437 - Ifges(5,6) * t433) * t430;
t389 = Ifges(6,1) * t402 + Ifges(6,4) * t401 + Ifges(6,5) * t429;
t388 = Ifges(6,4) * t402 + Ifges(6,2) * t401 + Ifges(6,6) * t429;
t387 = Ifges(6,5) * t402 + Ifges(6,6) * t401 + Ifges(6,3) * t429;
t362 = mrSges(6,2) * t375 - mrSges(6,3) * t371 + Ifges(6,1) * t379 + Ifges(6,4) * t378 + Ifges(6,5) * t427 + t387 * t401 - t388 * t429;
t361 = -mrSges(6,1) * t375 + mrSges(6,3) * t372 + Ifges(6,4) * t379 + Ifges(6,2) * t378 + Ifges(6,6) * t427 - t387 * t402 + t389 * t429;
t355 = -m(4) * g(3) - t448;
t352 = mrSges(5,2) * t383 - mrSges(5,3) * t380 + Ifges(5,1) * t408 + Ifges(5,4) * t407 + Ifges(5,5) * qJDD(4) - pkin(8) * t360 - qJD(4) * t399 - t361 * t432 + t362 * t436 - t398 * t457;
t351 = -mrSges(5,1) * t383 + mrSges(5,3) * t381 + Ifges(5,4) * t408 + Ifges(5,2) * t407 + Ifges(5,6) * qJDD(4) - pkin(4) * t444 + pkin(8) * t447 + qJD(4) * t400 + t436 * t361 + t432 * t362 - t398 * t456;
t344 = Ifges(5,3) * qJDD(4) + Ifges(6,3) * t427 + Ifges(5,6) * t407 + Ifges(5,5) * t408 - t401 * t389 + t402 * t388 + mrSges(5,1) * t380 - mrSges(5,2) * t381 + mrSges(4,1) * t392 - mrSges(3,3) * t394 + mrSges(6,1) * t371 - mrSges(6,2) * t372 + Ifges(6,6) * t378 + Ifges(6,5) * t379 + pkin(4) * t360 - qJ(3) * t355 + pkin(3) * t356 + (t399 * t437 + t400 * t433) * t430 + t459 * t428 + (t458 * t426) + (-mrSges(3,2) + mrSges(4,3)) * g(3);
t343 = -mrSges(4,1) * t390 + mrSges(3,3) * t395 - pkin(2) * t355 - pkin(3) * t442 + pkin(7) * t448 + t460 * g(3) - t437 * t351 - t433 * t352 + t459 * t426 - t458 * t428;
t342 = -mrSges(2,2) * g(3) - mrSges(2,3) * t417 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t440 - pkin(6) * t350 - t343 * t434 + t344 * t438;
t341 = Ifges(2,6) * qJDD(1) + t440 * Ifges(2,5) + mrSges(2,3) * t418 + t434 * t344 + t438 * t343 + pkin(1) * t448 + pkin(6) * t449 + (-pkin(1) * t462 + mrSges(2,1)) * g(3);
t1 = [-m(1) * g(1) + t450; -m(1) * g(2) + t455; (-m(1) - m(2) + t462) * g(3) - t448; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t455 - t435 * t341 + t439 * t342; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t450 + t439 * t341 + t435 * t342; pkin(1) * t350 + mrSges(2,1) * t417 - mrSges(2,2) * t418 + pkin(2) * t443 + qJ(3) * t441 + t437 * t352 - t433 * t351 - pkin(7) * t356 + mrSges(3,1) * t394 - mrSges(3,2) * t395 + mrSges(4,2) * t392 - mrSges(4,3) * t390 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + (-mrSges(4,2) * pkin(2) + Ifges(4,1) + Ifges(3,3)) * t428;];
tauB = t1;
