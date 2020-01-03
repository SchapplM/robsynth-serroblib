% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPRP2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPRP2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:33
% EndTime: 2019-12-31 19:49:35
% DurationCPUTime: 1.58s
% Computational Cost: add. (19910->205), mult. (25449->249), div. (0->0), fcn. (12352->8), ass. (0->84)
t453 = Ifges(5,1) + Ifges(6,1);
t449 = Ifges(5,4) - Ifges(6,5);
t448 = Ifges(5,5) + Ifges(6,4);
t452 = Ifges(5,2) + Ifges(6,3);
t447 = Ifges(5,6) - Ifges(6,6);
t451 = Ifges(5,3) + Ifges(6,2);
t450 = mrSges(5,3) + mrSges(6,2);
t417 = qJD(1) + qJD(2);
t423 = sin(qJ(4));
t446 = t417 * t423;
t426 = cos(qJ(4));
t445 = t417 * t426;
t420 = -g(3) + qJDD(3);
t444 = t426 * t420;
t425 = sin(qJ(1));
t428 = cos(qJ(1));
t409 = t425 * g(1) - t428 * g(2);
t403 = qJDD(1) * pkin(1) + t409;
t410 = -t428 * g(1) - t425 * g(2);
t430 = qJD(1) ^ 2;
t404 = -t430 * pkin(1) + t410;
t424 = sin(qJ(2));
t427 = cos(qJ(2));
t379 = t427 * t403 - t424 * t404;
t416 = qJDD(1) + qJDD(2);
t377 = t416 * pkin(2) + t379;
t380 = t424 * t403 + t427 * t404;
t415 = t417 ^ 2;
t378 = -t415 * pkin(2) + t380;
t421 = sin(pkin(8));
t422 = cos(pkin(8));
t373 = t421 * t377 + t422 * t378;
t371 = -t415 * pkin(3) + t416 * pkin(7) + t373;
t368 = t426 * t371 + t423 * t420;
t395 = (-mrSges(5,1) * t426 + mrSges(5,2) * t423) * t417;
t439 = qJD(4) * t417;
t397 = t426 * t416 - t423 * t439;
t405 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t446;
t393 = (-pkin(4) * t426 - qJ(5) * t423) * t417;
t429 = qJD(4) ^ 2;
t365 = -t429 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t393 * t445 + t368;
t394 = (-mrSges(6,1) * t426 - mrSges(6,3) * t423) * t417;
t406 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t446;
t433 = m(6) * t365 + qJDD(4) * mrSges(6,3) + qJD(4) * t406 + t394 * t445;
t360 = m(5) * t368 - qJDD(4) * mrSges(5,2) - qJD(4) * t405 + t395 * t445 + t450 * t397 + t433;
t367 = -t423 * t371 + t444;
t396 = t423 * t416 + t426 * t439;
t407 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t445;
t366 = -qJDD(4) * pkin(4) - t429 * qJ(5) - t444 + qJDD(5) + (t393 * t417 + t371) * t423;
t408 = mrSges(6,2) * t445 + qJD(4) * mrSges(6,3);
t432 = -m(6) * t366 + qJDD(4) * mrSges(6,1) + qJD(4) * t408;
t361 = m(5) * t367 + qJDD(4) * mrSges(5,1) + qJD(4) * t407 + (-t394 - t395) * t446 - t450 * t396 + t432;
t434 = t426 * t360 - t423 * t361;
t353 = m(4) * t373 - t415 * mrSges(4,1) - t416 * mrSges(4,2) + t434;
t372 = t422 * t377 - t421 * t378;
t370 = -t416 * pkin(3) - t415 * pkin(7) - t372;
t363 = -t397 * pkin(4) - t396 * qJ(5) + (-0.2e1 * qJD(5) * t423 + (pkin(4) * t423 - qJ(5) * t426) * qJD(4)) * t417 + t370;
t362 = m(6) * t363 - t397 * mrSges(6,1) - t396 * mrSges(6,3) - t406 * t446 - t408 * t445;
t431 = -m(5) * t370 + t397 * mrSges(5,1) - t396 * mrSges(5,2) - t405 * t446 + t407 * t445 - t362;
t356 = m(4) * t372 + t416 * mrSges(4,1) - t415 * mrSges(4,2) + t431;
t348 = t421 * t353 + t422 * t356;
t346 = m(3) * t379 + t416 * mrSges(3,1) - t415 * mrSges(3,2) + t348;
t435 = t422 * t353 - t421 * t356;
t347 = m(3) * t380 - t415 * mrSges(3,1) - t416 * mrSges(3,2) + t435;
t340 = t427 * t346 + t424 * t347;
t338 = m(2) * t409 + qJDD(1) * mrSges(2,1) - t430 * mrSges(2,2) + t340;
t436 = -t424 * t346 + t427 * t347;
t339 = m(2) * t410 - t430 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t436;
t443 = t428 * t338 + t425 * t339;
t354 = t423 * t360 + t426 * t361;
t442 = (-t449 * t423 - t452 * t426) * t417 - t447 * qJD(4);
t441 = (t448 * t423 + t447 * t426) * t417 + t451 * qJD(4);
t440 = (t453 * t423 + t449 * t426) * t417 + t448 * qJD(4);
t438 = m(4) * t420 + t354;
t437 = -t425 * t338 + t428 * t339;
t350 = mrSges(5,2) * t370 + mrSges(6,2) * t366 - mrSges(5,3) * t367 - mrSges(6,3) * t363 - qJ(5) * t362 + t442 * qJD(4) + t448 * qJDD(4) + t453 * t396 + t449 * t397 + t441 * t445;
t349 = -mrSges(5,1) * t370 - mrSges(6,1) * t363 + mrSges(6,2) * t365 + mrSges(5,3) * t368 - pkin(4) * t362 + t440 * qJD(4) + t447 * qJDD(4) + t449 * t396 + t452 * t397 - t441 * t446;
t342 = Ifges(4,6) * t416 + t415 * Ifges(4,5) - mrSges(4,1) * t420 + mrSges(4,3) * t373 - mrSges(5,1) * t367 + mrSges(5,2) * t368 + mrSges(6,1) * t366 - mrSges(6,3) * t365 - pkin(4) * t432 - qJ(5) * t433 - pkin(3) * t354 + (-qJ(5) * mrSges(6,2) - t447) * t397 + (pkin(4) * mrSges(6,2) - t448) * t396 - t451 * qJDD(4) + (t440 * t426 + (pkin(4) * t394 + t442) * t423) * t417;
t341 = mrSges(4,2) * t420 - mrSges(4,3) * t372 + Ifges(4,5) * t416 - t415 * Ifges(4,6) - pkin(7) * t354 - t423 * t349 + t426 * t350;
t334 = -mrSges(3,2) * g(3) - mrSges(3,3) * t379 + Ifges(3,5) * t416 - t415 * Ifges(3,6) - qJ(3) * t348 + t422 * t341 - t421 * t342;
t333 = mrSges(3,1) * g(3) + mrSges(3,3) * t380 + t415 * Ifges(3,5) + Ifges(3,6) * t416 - pkin(2) * t438 + qJ(3) * t435 + t421 * t341 + t422 * t342;
t332 = -mrSges(2,2) * g(3) - mrSges(2,3) * t409 + Ifges(2,5) * qJDD(1) - t430 * Ifges(2,6) - pkin(6) * t340 - t424 * t333 + t427 * t334;
t331 = Ifges(2,6) * qJDD(1) + t430 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t410 + t424 * t334 + t427 * t333 - pkin(1) * (-m(3) * g(3) + t438) + pkin(6) * t436;
t1 = [-m(1) * g(1) + t437; -m(1) * g(2) + t443; (-m(1) - m(2) - m(3)) * g(3) + t438; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t443 - t425 * t331 + t428 * t332; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t437 + t428 * t331 + t425 * t332; pkin(1) * t340 + mrSges(2,1) * t409 - mrSges(2,2) * t410 + pkin(2) * t348 + mrSges(3,1) * t379 - mrSges(3,2) * t380 + pkin(3) * t431 + pkin(7) * t434 + t423 * t350 + t426 * t349 + mrSges(4,1) * t372 - mrSges(4,2) * t373 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + (Ifges(4,3) + Ifges(3,3)) * t416;];
tauB = t1;
