% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRRP3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRRP3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP3_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP3_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP3_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP3_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:15
% EndTime: 2019-12-31 21:49:17
% DurationCPUTime: 1.60s
% Computational Cost: add. (24356->207), mult. (25449->251), div. (0->0), fcn. (12352->8), ass. (0->86)
t454 = Ifges(5,1) + Ifges(6,1);
t448 = Ifges(5,4) - Ifges(6,5);
t447 = Ifges(5,5) + Ifges(6,4);
t453 = Ifges(5,2) + Ifges(6,3);
t446 = Ifges(5,6) - Ifges(6,6);
t452 = Ifges(5,3) + Ifges(6,2);
t451 = -m(3) - m(4);
t426 = cos(qJ(4));
t450 = t426 * g(3);
t449 = mrSges(5,3) + mrSges(6,2);
t419 = qJD(1) + qJD(2);
t415 = qJD(3) + t419;
t422 = sin(qJ(4));
t445 = t415 * t422;
t444 = t415 * t426;
t425 = sin(qJ(1));
t429 = cos(qJ(1));
t411 = t425 * g(1) - t429 * g(2);
t408 = qJDD(1) * pkin(1) + t411;
t412 = -t429 * g(1) - t425 * g(2);
t431 = qJD(1) ^ 2;
t409 = -t431 * pkin(1) + t412;
t424 = sin(qJ(2));
t428 = cos(qJ(2));
t380 = t428 * t408 - t424 * t409;
t418 = qJDD(1) + qJDD(2);
t378 = t418 * pkin(2) + t380;
t381 = t424 * t408 + t428 * t409;
t417 = t419 ^ 2;
t379 = -t417 * pkin(2) + t381;
t423 = sin(qJ(3));
t427 = cos(qJ(3));
t374 = t423 * t378 + t427 * t379;
t413 = t415 ^ 2;
t414 = qJDD(3) + t418;
t372 = -t413 * pkin(3) + t414 * pkin(8) + t374;
t369 = -t422 * g(3) + t426 * t372;
t396 = (-mrSges(5,1) * t426 + mrSges(5,2) * t422) * t415;
t439 = qJD(4) * t415;
t398 = t426 * t414 - t422 * t439;
t404 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t445;
t394 = (-pkin(4) * t426 - qJ(5) * t422) * t415;
t430 = qJD(4) ^ 2;
t366 = -t430 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t394 * t444 + t369;
t395 = (-mrSges(6,1) * t426 - mrSges(6,3) * t422) * t415;
t405 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t445;
t434 = m(6) * t366 + qJDD(4) * mrSges(6,3) + qJD(4) * t405 + t395 * t444;
t361 = m(5) * t369 - qJDD(4) * mrSges(5,2) - qJD(4) * t404 + t396 * t444 + t449 * t398 + t434;
t368 = -t422 * t372 - t450;
t397 = t422 * t414 + t426 * t439;
t406 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t444;
t367 = -qJDD(4) * pkin(4) + t450 - t430 * qJ(5) + qJDD(5) + (t394 * t415 + t372) * t422;
t407 = mrSges(6,2) * t444 + qJD(4) * mrSges(6,3);
t433 = -m(6) * t367 + qJDD(4) * mrSges(6,1) + qJD(4) * t407;
t362 = m(5) * t368 + qJDD(4) * mrSges(5,1) + qJD(4) * t406 + (-t395 - t396) * t445 - t449 * t397 + t433;
t435 = t426 * t361 - t422 * t362;
t354 = m(4) * t374 - t413 * mrSges(4,1) - t414 * mrSges(4,2) + t435;
t373 = t427 * t378 - t423 * t379;
t371 = -t414 * pkin(3) - t413 * pkin(8) - t373;
t364 = -t398 * pkin(4) - t397 * qJ(5) + (-0.2e1 * qJD(5) * t422 + (pkin(4) * t422 - qJ(5) * t426) * qJD(4)) * t415 + t371;
t363 = m(6) * t364 - t398 * mrSges(6,1) - t397 * mrSges(6,3) - t405 * t445 - t407 * t444;
t432 = -m(5) * t371 + t398 * mrSges(5,1) - t397 * mrSges(5,2) - t404 * t445 + t406 * t444 - t363;
t357 = m(4) * t373 + t414 * mrSges(4,1) - t413 * mrSges(4,2) + t432;
t349 = t423 * t354 + t427 * t357;
t347 = m(3) * t380 + t418 * mrSges(3,1) - t417 * mrSges(3,2) + t349;
t436 = t427 * t354 - t423 * t357;
t348 = m(3) * t381 - t417 * mrSges(3,1) - t418 * mrSges(3,2) + t436;
t341 = t428 * t347 + t424 * t348;
t339 = m(2) * t411 + qJDD(1) * mrSges(2,1) - t431 * mrSges(2,2) + t341;
t437 = -t424 * t347 + t428 * t348;
t340 = m(2) * t412 - t431 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t437;
t443 = t429 * t339 + t425 * t340;
t355 = t422 * t361 + t426 * t362;
t442 = (-t448 * t422 - t453 * t426) * t415 - t446 * qJD(4);
t441 = (t447 * t422 + t446 * t426) * t415 + t452 * qJD(4);
t440 = (t454 * t422 + t448 * t426) * t415 + t447 * qJD(4);
t438 = -t425 * t339 + t429 * t340;
t351 = mrSges(5,2) * t371 + mrSges(6,2) * t367 - mrSges(5,3) * t368 - mrSges(6,3) * t364 - qJ(5) * t363 + t442 * qJD(4) + t447 * qJDD(4) + t454 * t397 + t448 * t398 + t441 * t444;
t350 = -mrSges(5,1) * t371 - mrSges(6,1) * t364 + mrSges(6,2) * t366 + mrSges(5,3) * t369 - pkin(4) * t363 + t440 * qJD(4) + t446 * qJDD(4) + t448 * t397 + t453 * t398 - t441 * t445;
t343 = Ifges(4,6) * t414 + t413 * Ifges(4,5) + mrSges(4,1) * g(3) + mrSges(4,3) * t374 - mrSges(5,1) * t368 + mrSges(5,2) * t369 + mrSges(6,1) * t367 - mrSges(6,3) * t366 - pkin(4) * t433 - qJ(5) * t434 - pkin(3) * t355 + (-qJ(5) * mrSges(6,2) - t446) * t398 + (pkin(4) * mrSges(6,2) - t447) * t397 - t452 * qJDD(4) + (t440 * t426 + (pkin(4) * t395 + t442) * t422) * t415;
t342 = -mrSges(4,2) * g(3) - mrSges(4,3) * t373 + Ifges(4,5) * t414 - t413 * Ifges(4,6) - pkin(8) * t355 - t422 * t350 + t426 * t351;
t335 = -mrSges(3,2) * g(3) - mrSges(3,3) * t380 + Ifges(3,5) * t418 - t417 * Ifges(3,6) - pkin(7) * t349 + t427 * t342 - t423 * t343;
t334 = Ifges(3,6) * t418 + t417 * Ifges(3,5) + mrSges(3,1) * g(3) + mrSges(3,3) * t381 + t423 * t342 + t427 * t343 - pkin(2) * (-m(4) * g(3) + t355) + pkin(7) * t436;
t333 = -mrSges(2,2) * g(3) - mrSges(2,3) * t411 + Ifges(2,5) * qJDD(1) - t431 * Ifges(2,6) - pkin(6) * t341 - t424 * t334 + t428 * t335;
t332 = Ifges(2,6) * qJDD(1) + t431 * Ifges(2,5) + mrSges(2,3) * t412 + t424 * t335 + t428 * t334 - pkin(1) * t355 + pkin(6) * t437 + (-pkin(1) * t451 + mrSges(2,1)) * g(3);
t1 = [-m(1) * g(1) + t438; -m(1) * g(2) + t443; (-m(1) - m(2) + t451) * g(3) + t355; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t443 - t425 * t332 + t429 * t333; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t438 + t429 * t332 + t425 * t333; pkin(1) * t341 + mrSges(2,1) * t411 - mrSges(2,2) * t412 + pkin(2) * t349 + mrSges(3,1) * t380 - mrSges(3,2) * t381 + t426 * t350 + pkin(3) * t432 + pkin(8) * t435 + mrSges(4,1) * t373 - mrSges(4,2) * t374 + t422 * t351 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(4,3) * t414 + Ifges(3,3) * t418 + Ifges(2,3) * qJDD(1);];
tauB = t1;
