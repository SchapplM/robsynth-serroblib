% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPPRR11
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPPRR11_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR11_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR11_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR11_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR11_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR11_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR11_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:37
% EndTime: 2019-12-31 18:05:38
% DurationCPUTime: 0.85s
% Computational Cost: add. (6306->216), mult. (11439->253), div. (0->0), fcn. (5333->6), ass. (0->81)
t418 = sin(qJ(1));
t421 = cos(qJ(1));
t399 = t418 * g(1) - t421 * g(2);
t423 = qJD(1) ^ 2;
t382 = -qJDD(1) * pkin(1) - t423 * qJ(2) + qJDD(2) - t399;
t377 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) + t382;
t400 = -t421 * g(1) - t418 * g(2);
t447 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t400;
t446 = -m(3) - m(4);
t417 = sin(qJ(4));
t445 = t417 * g(3);
t444 = mrSges(2,1) - mrSges(3,2);
t381 = t423 * pkin(1) - t447;
t378 = qJDD(3) + (-pkin(1) - qJ(3)) * t423 + t447;
t375 = -qJDD(1) * pkin(6) + t378;
t420 = cos(qJ(4));
t370 = -t420 * g(3) + t417 * t375;
t393 = (mrSges(5,1) * t417 + mrSges(5,2) * t420) * qJD(1);
t440 = qJD(1) * qJD(4);
t434 = t420 * t440;
t395 = -t417 * qJDD(1) - t434;
t442 = qJD(1) * t420;
t398 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t442;
t374 = -t423 * pkin(6) - t377;
t435 = t417 * t440;
t396 = t420 * qJDD(1) - t435;
t359 = (-t396 + t435) * pkin(7) + (-t395 + t434) * pkin(4) + t374;
t394 = (pkin(4) * t417 - pkin(7) * t420) * qJD(1);
t422 = qJD(4) ^ 2;
t441 = t417 * qJD(1);
t361 = -t422 * pkin(4) + qJDD(4) * pkin(7) - t394 * t441 + t370;
t416 = sin(qJ(5));
t419 = cos(qJ(5));
t357 = t419 * t359 - t416 * t361;
t391 = t419 * qJD(4) - t416 * t442;
t368 = t391 * qJD(5) + t416 * qJDD(4) + t419 * t396;
t392 = t416 * qJD(4) + t419 * t442;
t372 = -t391 * mrSges(6,1) + t392 * mrSges(6,2);
t401 = qJD(5) + t441;
t379 = -t401 * mrSges(6,2) + t391 * mrSges(6,3);
t390 = qJDD(5) - t395;
t355 = m(6) * t357 + t390 * mrSges(6,1) - t368 * mrSges(6,3) - t392 * t372 + t401 * t379;
t358 = t416 * t359 + t419 * t361;
t367 = -t392 * qJD(5) + t419 * qJDD(4) - t416 * t396;
t380 = t401 * mrSges(6,1) - t392 * mrSges(6,3);
t356 = m(6) * t358 - t390 * mrSges(6,2) + t367 * mrSges(6,3) + t391 * t372 - t401 * t380;
t431 = -t416 * t355 + t419 * t356;
t346 = m(5) * t370 - qJDD(4) * mrSges(5,2) + t395 * mrSges(5,3) - qJD(4) * t398 - t393 * t441 + t431;
t369 = t420 * t375 + t445;
t397 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t441;
t360 = -qJDD(4) * pkin(4) - t422 * pkin(7) - t445 + (qJD(1) * t394 - t375) * t420;
t425 = -m(6) * t360 + t367 * mrSges(6,1) - t368 * mrSges(6,2) + t391 * t379 - t392 * t380;
t351 = m(5) * t369 + qJDD(4) * mrSges(5,1) - t396 * mrSges(5,3) + qJD(4) * t397 - t393 * t442 + t425;
t340 = t417 * t346 + t420 * t351;
t430 = -m(4) * t378 - qJDD(1) * mrSges(4,2) - t340;
t426 = -m(3) * t381 + t423 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t430;
t338 = m(2) * t400 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - mrSges(4,3)) * t423 + t426;
t347 = t419 * t355 + t416 * t356;
t427 = -m(5) * t374 + t395 * mrSges(5,1) - t396 * mrSges(5,2) - t397 * t441 - t398 * t442 - t347;
t343 = m(4) * t377 - t423 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t427;
t424 = -m(3) * t382 + t423 * mrSges(3,3) - t343;
t342 = m(2) * t399 - t423 * mrSges(2,2) + qJDD(1) * t444 + t424;
t443 = t418 * t338 + t421 * t342;
t437 = Ifges(2,5) - Ifges(3,4) + Ifges(4,5);
t436 = Ifges(2,6) - Ifges(3,5) - Ifges(4,4);
t433 = t421 * t338 - t418 * t342;
t432 = t420 * t346 - t417 * t351;
t385 = (Ifges(5,5) * qJD(4)) + (Ifges(5,1) * t420 - Ifges(5,4) * t417) * qJD(1);
t384 = (Ifges(5,6) * qJD(4)) + (Ifges(5,4) * t420 - Ifges(5,2) * t417) * qJD(1);
t383 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t420 - Ifges(5,6) * t417) * qJD(1);
t364 = Ifges(6,1) * t392 + Ifges(6,4) * t391 + Ifges(6,5) * t401;
t363 = Ifges(6,4) * t392 + Ifges(6,2) * t391 + Ifges(6,6) * t401;
t362 = Ifges(6,5) * t392 + Ifges(6,6) * t391 + Ifges(6,3) * t401;
t349 = mrSges(6,2) * t360 - mrSges(6,3) * t357 + Ifges(6,1) * t368 + Ifges(6,4) * t367 + Ifges(6,5) * t390 + t391 * t362 - t401 * t363;
t348 = -mrSges(6,1) * t360 + mrSges(6,3) * t358 + Ifges(6,4) * t368 + Ifges(6,2) * t367 + Ifges(6,6) * t390 - t392 * t362 + t401 * t364;
t339 = g(3) * t446 + t432;
t335 = -mrSges(5,1) * t374 - mrSges(6,1) * t357 + mrSges(6,2) * t358 + mrSges(5,3) * t370 + Ifges(5,4) * t396 - Ifges(6,5) * t368 + Ifges(5,2) * t395 + Ifges(5,6) * qJDD(4) - Ifges(6,6) * t367 - Ifges(6,3) * t390 - pkin(4) * t347 + qJD(4) * t385 - t392 * t363 + t391 * t364 - t383 * t442;
t334 = mrSges(5,2) * t374 - mrSges(5,3) * t369 + Ifges(5,1) * t396 + Ifges(5,4) * t395 + Ifges(5,5) * qJDD(4) - pkin(7) * t347 - qJD(4) * t384 - t416 * t348 + t419 * t349 - t383 * t441;
t333 = Ifges(5,3) * qJDD(4) - pkin(1) * t339 + pkin(3) * t340 + mrSges(5,1) * t369 - mrSges(5,2) * t370 + mrSges(4,1) * t378 - mrSges(3,1) * t381 + pkin(4) * t425 + Ifges(5,6) * t395 + Ifges(5,5) * t396 + mrSges(2,3) * t400 + t416 * t349 + pkin(7) * t431 - qJ(3) * t432 + t419 * t348 - pkin(2) * t430 + (t420 * t384 + t417 * t385) * qJD(1) + (-pkin(2) * mrSges(4,3) + t437) * t423 + t436 * qJDD(1) + (qJ(3) * m(4) + mrSges(4,3) + t444) * g(3);
t332 = -qJ(2) * t339 - mrSges(2,3) * t399 + pkin(2) * t343 + mrSges(3,1) * t382 + t417 * t334 + t420 * t335 + pkin(3) * t427 + pkin(6) * t432 + mrSges(4,1) * t377 - t436 * t423 + t437 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3) + mrSges(4,2)) * g(3);
t1 = [-m(1) * g(1) + t433; -m(1) * g(2) + t443; (-m(1) - m(2) + t446) * g(3) + t432; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t443 + t421 * t332 - t418 * t333; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t433 + t418 * t332 + t421 * t333; qJ(2) * (-t423 * mrSges(4,3) + t426) + pkin(1) * t424 + mrSges(2,1) * t399 - mrSges(2,2) * t400 - qJ(3) * t343 + mrSges(3,2) * t382 - mrSges(3,3) * t381 - pkin(6) * t340 + t420 * t334 - t417 * t335 + mrSges(4,2) * t378 - mrSges(4,3) * t377 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1);];
tauB = t1;
