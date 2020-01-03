% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPPR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPPR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR5_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR5_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR5_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR5_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:29:09
% EndTime: 2019-12-31 19:29:11
% DurationCPUTime: 1.58s
% Computational Cost: add. (8418->256), mult. (19627->312), div. (0->0), fcn. (12602->8), ass. (0->104)
t468 = -2 * qJD(3);
t467 = Ifges(4,1) + Ifges(5,1);
t458 = Ifges(4,4) - Ifges(5,5);
t457 = Ifges(4,5) + Ifges(5,4);
t466 = Ifges(4,2) + Ifges(5,3);
t465 = -Ifges(5,2) - Ifges(4,3);
t456 = Ifges(4,6) - Ifges(5,6);
t422 = sin(qJ(2));
t425 = cos(qJ(2));
t443 = qJD(1) * qJD(2);
t406 = qJDD(1) * t425 - t422 * t443;
t447 = qJD(1) * t422;
t407 = qJD(2) * pkin(2) - qJ(3) * t447;
t419 = t425 ^ 2;
t428 = qJD(1) ^ 2;
t423 = sin(qJ(1));
t426 = cos(qJ(1));
t448 = t423 * g(1) - t426 * g(2);
t437 = qJDD(1) * pkin(1) + t448;
t364 = -t406 * pkin(2) + t407 * t447 - (qJ(3) * t419 + pkin(6)) * t428 + qJDD(3) - t437;
t405 = qJDD(1) * t422 + t425 * t443;
t420 = sin(pkin(8));
t455 = cos(pkin(8));
t381 = t455 * t405 + t420 * t406;
t446 = qJD(1) * t425;
t395 = t420 * t447 - t455 * t446;
t445 = qJD(2) * t395;
t464 = t364 + (-t381 + t445) * qJ(4);
t438 = -g(1) * t426 - g(2) * t423;
t402 = -pkin(1) * t428 + qJDD(1) * pkin(6) + t438;
t454 = t422 * t402;
t461 = pkin(2) * t428;
t360 = qJDD(2) * pkin(2) - t405 * qJ(3) - t454 + (qJ(3) * t443 + t422 * t461 - g(3)) * t425;
t385 = -g(3) * t422 + t425 * t402;
t361 = qJ(3) * t406 - qJD(2) * t407 - t419 * t461 + t385;
t396 = (t420 * t425 + t455 * t422) * qJD(1);
t347 = t455 * t360 - t420 * t361 + t396 * t468;
t380 = t405 * t420 - t455 * t406;
t386 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t395;
t387 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t396;
t388 = -qJD(2) * mrSges(5,1) + mrSges(5,2) * t396;
t346 = -0.2e1 * qJD(4) * t396 + (qJD(2) * t396 + t380) * pkin(3) + t464;
t389 = -mrSges(5,2) * t395 + qJD(2) * mrSges(5,3);
t390 = -qJD(2) * pkin(4) - pkin(7) * t396;
t394 = t395 ^ 2;
t462 = 0.2e1 * qJD(4);
t341 = -t394 * pkin(7) + (-pkin(3) - pkin(4)) * t380 + (-pkin(3) * qJD(2) + t390 + t462) * t396 - t464;
t421 = sin(qJ(5));
t424 = cos(qJ(5));
t374 = t395 * t421 + t396 * t424;
t350 = -qJD(5) * t374 + t380 * t424 - t381 * t421;
t373 = t395 * t424 - t396 * t421;
t351 = qJD(5) * t373 + t380 * t421 + t381 * t424;
t415 = -qJD(2) + qJD(5);
t365 = -mrSges(6,2) * t415 + mrSges(6,3) * t373;
t366 = mrSges(6,1) * t415 - mrSges(6,3) * t374;
t433 = m(6) * t341 - t350 * mrSges(6,1) + t351 * mrSges(6,2) - t373 * t365 + t374 * t366;
t431 = -m(5) * t346 - t380 * mrSges(5,1) - t395 * t389 + t433;
t332 = m(4) * t364 + t380 * mrSges(4,1) + (mrSges(4,2) - mrSges(5,3)) * t381 + t395 * t386 + (t387 - t388) * t396 - t431;
t459 = -mrSges(4,3) - mrSges(5,2);
t348 = t420 * t360 + t455 * t361 + t395 * t468;
t376 = pkin(3) * t395 - qJ(4) * t396;
t427 = qJD(2) ^ 2;
t343 = -pkin(3) * t427 + qJDD(2) * qJ(4) + qJD(2) * t462 - t395 * t376 + t348;
t344 = -qJDD(2) * pkin(3) - t427 * qJ(4) + t396 * t376 + qJDD(4) - t347;
t339 = (-t381 - t445) * pkin(7) + (t395 * t396 - qJDD(2)) * pkin(4) + t344;
t340 = -pkin(4) * t394 + pkin(7) * t380 + qJD(2) * t390 + t343;
t337 = t339 * t424 - t340 * t421;
t356 = -mrSges(6,1) * t373 + mrSges(6,2) * t374;
t414 = -qJDD(2) + qJDD(5);
t335 = m(6) * t337 + mrSges(6,1) * t414 - mrSges(6,3) * t351 - t356 * t374 + t365 * t415;
t338 = t339 * t421 + t340 * t424;
t336 = m(6) * t338 - mrSges(6,2) * t414 + mrSges(6,3) * t350 + t356 * t373 - t366 * t415;
t440 = -t421 * t335 + t424 * t336;
t435 = m(5) * t343 + qJDD(2) * mrSges(5,3) + qJD(2) * t388 + t440;
t377 = mrSges(5,1) * t395 - mrSges(5,3) * t396;
t450 = -mrSges(4,1) * t395 - mrSges(4,2) * t396 - t377;
t326 = m(4) * t348 - qJDD(2) * mrSges(4,2) - qJD(2) * t387 + t459 * t380 + t450 * t395 + t435;
t329 = t424 * t335 + t421 * t336;
t432 = -m(5) * t344 + qJDD(2) * mrSges(5,1) + qJD(2) * t389 - t329;
t327 = m(4) * t347 + qJDD(2) * mrSges(4,1) + qJD(2) * t386 + t459 * t381 + t450 * t396 + t432;
t322 = t420 * t326 + t455 * t327;
t453 = -t456 * qJD(2) + t466 * t395 - t458 * t396;
t452 = t465 * qJD(2) + t456 * t395 - t457 * t396;
t451 = t457 * qJD(2) - t458 * t395 + t467 * t396;
t441 = t455 * t326 - t327 * t420;
t353 = Ifges(6,4) * t374 + Ifges(6,2) * t373 + Ifges(6,6) * t415;
t354 = Ifges(6,1) * t374 + Ifges(6,4) * t373 + Ifges(6,5) * t415;
t430 = mrSges(6,1) * t337 - mrSges(6,2) * t338 + Ifges(6,5) * t351 + Ifges(6,6) * t350 + Ifges(6,3) * t414 + t374 * t353 - t373 * t354;
t409 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t446;
t408 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t447;
t404 = (-mrSges(3,1) * t425 + mrSges(3,2) * t422) * qJD(1);
t401 = -t428 * pkin(6) - t437;
t399 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t422 + Ifges(3,4) * t425) * qJD(1);
t398 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t422 + Ifges(3,2) * t425) * qJD(1);
t384 = -t425 * g(3) - t454;
t352 = Ifges(6,5) * t374 + Ifges(6,6) * t373 + Ifges(6,3) * t415;
t333 = -t381 * mrSges(5,3) - t396 * t388 - t431;
t331 = mrSges(6,2) * t341 - mrSges(6,3) * t337 + Ifges(6,1) * t351 + Ifges(6,4) * t350 + Ifges(6,5) * t414 + t352 * t373 - t353 * t415;
t330 = -mrSges(6,1) * t341 + mrSges(6,3) * t338 + Ifges(6,4) * t351 + Ifges(6,2) * t350 + Ifges(6,6) * t414 - t352 * t374 + t354 * t415;
t328 = t381 * mrSges(5,2) + t396 * t377 - t432;
t321 = mrSges(4,2) * t364 + mrSges(5,2) * t344 - mrSges(4,3) * t347 - mrSges(5,3) * t346 - pkin(7) * t329 - qJ(4) * t333 + t453 * qJD(2) + t457 * qJDD(2) - t421 * t330 + t424 * t331 - t458 * t380 + t467 * t381 + t452 * t395;
t320 = -mrSges(4,1) * t364 - mrSges(5,1) * t346 + mrSges(5,2) * t343 + mrSges(4,3) * t348 - pkin(3) * t333 + pkin(4) * t433 - pkin(7) * t440 + t451 * qJD(2) + t456 * qJDD(2) - t424 * t330 - t421 * t331 - t466 * t380 + t458 * t381 + t452 * t396;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t448 - mrSges(2,2) * t438 + t422 * (mrSges(3,2) * t401 - mrSges(3,3) * t384 + Ifges(3,1) * t405 + Ifges(3,4) * t406 + Ifges(3,5) * qJDD(2) - qJ(3) * t322 - qJD(2) * t398 - t420 * t320 + t455 * t321) + t425 * (-mrSges(3,1) * t401 + mrSges(3,3) * t385 + Ifges(3,4) * t405 + Ifges(3,2) * t406 + Ifges(3,6) * qJDD(2) - pkin(2) * t332 + qJ(3) * t441 + qJD(2) * t399 + t455 * t320 + t420 * t321) + pkin(1) * (-m(3) * t401 - t405 * mrSges(3,2) + t406 * mrSges(3,1) + (-t408 * t422 + t409 * t425) * qJD(1) - t332) + pkin(6) * (t425 * (m(3) * t385 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t406 - qJD(2) * t408 + t404 * t446 + t441) - t422 * (m(3) * t384 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t405 + qJD(2) * t409 - t404 * t447 + t322)); -t430 - t453 * t396 + (t422 * t398 - t425 * t399) * qJD(1) + (Ifges(3,3) - t465) * qJDD(2) + (-qJ(4) * t377 + t451) * t395 + (-mrSges(5,2) * qJ(4) - t456) * t380 + t457 * t381 + Ifges(3,5) * t405 + Ifges(3,6) * t406 - mrSges(3,2) * t385 + mrSges(3,1) * t384 - mrSges(4,2) * t348 + mrSges(5,3) * t343 - mrSges(5,1) * t344 + mrSges(4,1) * t347 + qJ(4) * t435 - pkin(4) * t329 - pkin(3) * t328 + pkin(2) * t322; t332; t328; t430;];
tauJ = t1;
