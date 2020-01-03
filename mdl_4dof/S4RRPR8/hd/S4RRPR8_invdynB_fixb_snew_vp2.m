% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RRPR8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR8_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR8_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR8_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR8_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR8_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:07:54
% EndTime: 2019-12-31 17:07:56
% DurationCPUTime: 1.08s
% Computational Cost: add. (6062->220), mult. (12539->273), div. (0->0), fcn. (6315->6), ass. (0->88)
t466 = Ifges(3,1) + Ifges(4,1);
t460 = Ifges(3,4) - Ifges(4,5);
t459 = Ifges(3,5) + Ifges(4,4);
t465 = Ifges(3,2) + Ifges(4,3);
t458 = Ifges(3,6) - Ifges(4,6);
t464 = Ifges(3,3) + Ifges(4,2);
t463 = 2 * qJD(3);
t439 = qJD(1) ^ 2;
t462 = pkin(5) * t439;
t461 = mrSges(3,3) + mrSges(4,2);
t436 = cos(qJ(2));
t457 = qJ(3) * t436;
t434 = sin(qJ(1));
t437 = cos(qJ(1));
t420 = -g(1) * t437 - g(2) * t434;
t401 = -pkin(1) * t439 + qJDD(1) * pkin(5) + t420;
t433 = sin(qJ(2));
t384 = -g(3) * t433 + t436 * t401;
t409 = (-mrSges(3,1) * t436 + mrSges(3,2) * t433) * qJD(1);
t450 = qJD(1) * qJD(2);
t411 = qJDD(1) * t436 - t433 * t450;
t452 = qJD(1) * t433;
t414 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t452;
t407 = (-pkin(2) * t436 - qJ(3) * t433) * qJD(1);
t438 = qJD(2) ^ 2;
t451 = qJD(1) * t436;
t371 = -pkin(2) * t438 + qJDD(2) * qJ(3) + qJD(2) * t463 + t407 * t451 + t384;
t408 = (-mrSges(4,1) * t436 - mrSges(4,3) * t433) * qJD(1);
t415 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t452;
t418 = -qJD(2) * pkin(3) - pkin(6) * t452;
t431 = t436 ^ 2;
t367 = -pkin(3) * t431 * t439 - pkin(6) * t411 + qJD(2) * t418 + t371;
t383 = -t436 * g(3) - t433 * t401;
t374 = -qJDD(2) * pkin(2) - qJ(3) * t438 + t407 * t452 + qJDD(3) - t383;
t449 = t436 * t450;
t410 = qJDD(1) * t433 + t449;
t368 = (-t410 + t449) * pkin(6) + (-t433 * t436 * t439 - qJDD(2)) * pkin(3) + t374;
t432 = sin(qJ(4));
t435 = cos(qJ(4));
t363 = -t367 * t432 + t368 * t435;
t398 = (-t432 * t433 - t435 * t436) * qJD(1);
t376 = qJD(4) * t398 + t410 * t435 - t411 * t432;
t399 = (-t432 * t436 + t433 * t435) * qJD(1);
t382 = -mrSges(5,1) * t398 + mrSges(5,2) * t399;
t427 = -qJD(2) + qJD(4);
t385 = -mrSges(5,2) * t427 + mrSges(5,3) * t398;
t426 = -qJDD(2) + qJDD(4);
t361 = m(5) * t363 + mrSges(5,1) * t426 - mrSges(5,3) * t376 - t382 * t399 + t385 * t427;
t364 = t367 * t435 + t368 * t432;
t375 = -qJD(4) * t399 - t410 * t432 - t411 * t435;
t386 = mrSges(5,1) * t427 - mrSges(5,3) * t399;
t362 = m(5) * t364 - mrSges(5,2) * t426 + mrSges(5,3) * t375 + t382 * t398 - t386 * t427;
t446 = -t361 * t432 + t435 * t362;
t442 = m(4) * t371 + qJDD(2) * mrSges(4,3) + qJD(2) * t415 + t408 * t451 + t446;
t352 = m(3) * t384 - qJDD(2) * mrSges(3,2) - qJD(2) * t414 + t409 * t451 + t461 * t411 + t442;
t416 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t451;
t354 = t361 * t435 + t432 * t362;
t417 = mrSges(4,2) * t451 + qJD(2) * mrSges(4,3);
t441 = -m(4) * t374 + qJDD(2) * mrSges(4,1) + qJD(2) * t417 - t354;
t353 = m(3) * t383 + qJDD(2) * mrSges(3,1) + qJD(2) * t416 - t461 * t410 + (-t408 - t409) * t452 + t441;
t447 = t436 * t352 - t353 * t433;
t347 = m(2) * t420 - mrSges(2,1) * t439 - qJDD(1) * mrSges(2,2) + t447;
t419 = t434 * g(1) - t437 * g(2);
t445 = qJDD(1) * pkin(1) + t419;
t443 = -qJ(3) * t410 - t445;
t369 = -pkin(2) * t411 - t462 + (-0.2e1 * qJD(3) * t433 + (pkin(2) * t433 - t457) * qJD(2)) * qJD(1) + t443;
t366 = (-pkin(6) * t431 + pkin(5)) * t439 + (pkin(2) + pkin(3)) * t411 + (qJD(2) * t457 + (-pkin(2) * qJD(2) + t418 + t463) * t433) * qJD(1) - t443;
t444 = -m(5) * t366 + t375 * mrSges(5,1) - t376 * mrSges(5,2) + t398 * t385 - t399 * t386;
t359 = m(4) * t369 - mrSges(4,1) * t411 - t410 * mrSges(4,3) - t415 * t452 - t417 * t451 + t444;
t400 = -t445 - t462;
t440 = -m(3) * t400 + t411 * mrSges(3,1) - mrSges(3,2) * t410 - t414 * t452 + t416 * t451 - t359;
t356 = m(2) * t419 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t439 + t440;
t456 = t434 * t347 + t437 * t356;
t348 = t433 * t352 + t436 * t353;
t455 = t464 * qJD(2) + (t459 * t433 + t458 * t436) * qJD(1);
t454 = -t458 * qJD(2) + (-t460 * t433 - t465 * t436) * qJD(1);
t453 = t459 * qJD(2) + (t466 * t433 + t460 * t436) * qJD(1);
t448 = t437 * t347 - t356 * t434;
t379 = Ifges(5,1) * t399 + Ifges(5,4) * t398 + Ifges(5,5) * t427;
t378 = Ifges(5,4) * t399 + Ifges(5,2) * t398 + Ifges(5,6) * t427;
t377 = Ifges(5,5) * t399 + Ifges(5,6) * t398 + Ifges(5,3) * t427;
t358 = mrSges(5,2) * t366 - mrSges(5,3) * t363 + Ifges(5,1) * t376 + Ifges(5,4) * t375 + Ifges(5,5) * t426 + t377 * t398 - t378 * t427;
t357 = -mrSges(5,1) * t366 + mrSges(5,3) * t364 + Ifges(5,4) * t376 + Ifges(5,2) * t375 + Ifges(5,6) * t426 - t377 * t399 + t379 * t427;
t344 = mrSges(3,2) * t400 + mrSges(4,2) * t374 - mrSges(3,3) * t383 - mrSges(4,3) * t369 - pkin(6) * t354 - qJ(3) * t359 + t454 * qJD(2) + t459 * qJDD(2) - t432 * t357 + t435 * t358 + t466 * t410 + t460 * t411 + t455 * t451;
t343 = -mrSges(3,1) * t400 - mrSges(4,1) * t369 + mrSges(4,2) * t371 + mrSges(3,3) * t384 - pkin(2) * t359 - pkin(3) * t444 - pkin(6) * t446 + t453 * qJD(2) + t458 * qJDD(2) - t435 * t357 - t432 * t358 + t460 * t410 + t465 * t411 - t455 * t452;
t342 = -mrSges(4,3) * t371 + mrSges(4,1) * t374 + Ifges(5,6) * t375 + Ifges(5,5) * t376 - mrSges(3,1) * t383 + mrSges(3,2) * t384 - t398 * t379 + t399 * t378 + mrSges(2,3) * t420 + Ifges(5,3) * t426 - qJ(3) * t442 + t439 * Ifges(2,5) - pkin(1) * t348 + pkin(3) * t354 + mrSges(5,1) * t363 - mrSges(5,2) * t364 + Ifges(2,6) * qJDD(1) + mrSges(2,1) * g(3) - pkin(2) * t441 + (-qJ(3) * mrSges(4,2) - t458) * t411 + (pkin(2) * mrSges(4,2) - t459) * t410 - t464 * qJDD(2) + (t453 * t436 + (pkin(2) * t408 + t454) * t433) * qJD(1);
t341 = -mrSges(2,2) * g(3) - mrSges(2,3) * t419 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t439 - pkin(5) * t348 - t343 * t433 + t344 * t436;
t1 = [-m(1) * g(1) + t448; -m(1) * g(2) + t456; (-m(1) - m(2)) * g(3) + t348; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t456 + t437 * t341 - t434 * t342; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t448 + t434 * t341 + t437 * t342; -mrSges(1,1) * g(2) + mrSges(2,1) * t419 + mrSges(1,2) * g(1) - mrSges(2,2) * t420 + Ifges(2,3) * qJDD(1) + pkin(1) * t440 + pkin(5) * t447 + t436 * t343 + t433 * t344;];
tauB = t1;
