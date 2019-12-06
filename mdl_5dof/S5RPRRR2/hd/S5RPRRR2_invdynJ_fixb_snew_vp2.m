% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRRR2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRRR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR2_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR2_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR2_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR2_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR2_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:11:51
% EndTime: 2019-12-05 18:11:54
% DurationCPUTime: 2.90s
% Computational Cost: add. (27779->245), mult. (69076->315), div. (0->0), fcn. (52617->10), ass. (0->107)
t425 = qJD(1) ^ 2;
t448 = pkin(2) * t425;
t447 = pkin(6) * qJDD(1);
t420 = sin(qJ(1));
t424 = cos(qJ(1));
t436 = -t424 * g(1) - t420 * g(2);
t401 = -t425 * pkin(1) + qJDD(1) * qJ(2) + t436;
t415 = sin(pkin(9));
t416 = cos(pkin(9));
t442 = qJD(1) * qJD(2);
t440 = -t416 * g(3) - 0.2e1 * t415 * t442;
t378 = (t416 * t448 - t401 - t447) * t415 + t440;
t392 = -t415 * g(3) + (t401 + 0.2e1 * t442) * t416;
t413 = t416 ^ 2;
t379 = -t413 * t448 + t416 * t447 + t392;
t419 = sin(qJ(3));
t423 = cos(qJ(3));
t361 = t423 * t378 - t419 * t379;
t433 = t415 * t423 + t416 * t419;
t432 = -t415 * t419 + t416 * t423;
t399 = t432 * qJD(1);
t443 = t399 * qJD(3);
t390 = t433 * qJDD(1) + t443;
t400 = t433 * qJD(1);
t345 = (-t390 + t443) * pkin(7) + (t399 * t400 + qJDD(3)) * pkin(3) + t361;
t362 = t419 * t378 + t423 * t379;
t389 = -t400 * qJD(3) + t432 * qJDD(1);
t395 = qJD(3) * pkin(3) - t400 * pkin(7);
t398 = t399 ^ 2;
t352 = -t398 * pkin(3) + t389 * pkin(7) - qJD(3) * t395 + t362;
t418 = sin(qJ(4));
t422 = cos(qJ(4));
t335 = t422 * t345 - t418 * t352;
t384 = t422 * t399 - t418 * t400;
t360 = t384 * qJD(4) + t418 * t389 + t422 * t390;
t385 = t418 * t399 + t422 * t400;
t411 = qJDD(3) + qJDD(4);
t414 = qJD(3) + qJD(4);
t330 = (t384 * t414 - t360) * pkin(8) + (t384 * t385 + t411) * pkin(4) + t335;
t336 = t418 * t345 + t422 * t352;
t359 = -t385 * qJD(4) + t422 * t389 - t418 * t390;
t376 = t414 * pkin(4) - t385 * pkin(8);
t380 = t384 ^ 2;
t331 = -t380 * pkin(4) + t359 * pkin(8) - t414 * t376 + t336;
t417 = sin(qJ(5));
t421 = cos(qJ(5));
t328 = t421 * t330 - t417 * t331;
t369 = t421 * t384 - t417 * t385;
t341 = t369 * qJD(5) + t417 * t359 + t421 * t360;
t370 = t417 * t384 + t421 * t385;
t350 = -t369 * mrSges(6,1) + t370 * mrSges(6,2);
t409 = qJD(5) + t414;
t363 = -t409 * mrSges(6,2) + t369 * mrSges(6,3);
t408 = qJDD(5) + t411;
t325 = m(6) * t328 + t408 * mrSges(6,1) - t341 * mrSges(6,3) - t370 * t350 + t409 * t363;
t329 = t417 * t330 + t421 * t331;
t340 = -t370 * qJD(5) + t421 * t359 - t417 * t360;
t364 = t409 * mrSges(6,1) - t370 * mrSges(6,3);
t326 = m(6) * t329 - t408 * mrSges(6,2) + t340 * mrSges(6,3) + t369 * t350 - t409 * t364;
t318 = t421 * t325 + t417 * t326;
t371 = -t384 * mrSges(5,1) + t385 * mrSges(5,2);
t374 = -t414 * mrSges(5,2) + t384 * mrSges(5,3);
t315 = m(5) * t335 + t411 * mrSges(5,1) - t360 * mrSges(5,3) - t385 * t371 + t414 * t374 + t318;
t375 = t414 * mrSges(5,1) - t385 * mrSges(5,3);
t437 = -t417 * t325 + t421 * t326;
t316 = m(5) * t336 - t411 * mrSges(5,2) + t359 * mrSges(5,3) + t384 * t371 - t414 * t375 + t437;
t311 = t422 * t315 + t418 * t316;
t387 = -t399 * mrSges(4,1) + t400 * mrSges(4,2);
t393 = -qJD(3) * mrSges(4,2) + t399 * mrSges(4,3);
t309 = m(4) * t361 + qJDD(3) * mrSges(4,1) - t390 * mrSges(4,3) + qJD(3) * t393 - t400 * t387 + t311;
t394 = qJD(3) * mrSges(4,1) - t400 * mrSges(4,3);
t438 = -t418 * t315 + t422 * t316;
t310 = m(4) * t362 - qJDD(3) * mrSges(4,2) + t389 * mrSges(4,3) - qJD(3) * t394 + t399 * t387 + t438;
t446 = t423 * t309 + t419 * t310;
t445 = -t415 ^ 2 - t413;
t441 = t420 * g(1) - t424 * g(2);
t439 = -t419 * t309 + t423 * t310;
t435 = qJDD(2) - t441;
t434 = -t416 * mrSges(3,1) + t415 * mrSges(3,2);
t431 = mrSges(3,3) * qJDD(1) + t425 * t434;
t388 = (-pkin(2) * t416 - pkin(1)) * qJDD(1) + (t445 * pkin(6) - qJ(2)) * t425 + t435;
t355 = -t389 * pkin(3) - t398 * pkin(7) + t400 * t395 + t388;
t333 = -t359 * pkin(4) - t380 * pkin(8) + t385 * t376 + t355;
t430 = m(6) * t333 - t340 * mrSges(6,1) + t341 * mrSges(6,2) - t369 * t363 + t370 * t364;
t347 = Ifges(6,4) * t370 + Ifges(6,2) * t369 + Ifges(6,6) * t409;
t348 = Ifges(6,1) * t370 + Ifges(6,4) * t369 + Ifges(6,5) * t409;
t429 = mrSges(6,1) * t328 - mrSges(6,2) * t329 + Ifges(6,5) * t341 + Ifges(6,6) * t340 + Ifges(6,3) * t408 + t370 * t347 - t369 * t348;
t428 = m(5) * t355 - t359 * mrSges(5,1) + t360 * mrSges(5,2) - t384 * t374 + t385 * t375 + t430;
t366 = Ifges(5,4) * t385 + Ifges(5,2) * t384 + Ifges(5,6) * t414;
t367 = Ifges(5,1) * t385 + Ifges(5,4) * t384 + Ifges(5,5) * t414;
t427 = mrSges(5,1) * t335 - mrSges(5,2) * t336 + Ifges(5,5) * t360 + Ifges(5,6) * t359 + Ifges(5,3) * t411 + pkin(4) * t318 + t385 * t366 - t384 * t367 + t429;
t426 = m(4) * t388 - t389 * mrSges(4,1) + t390 * mrSges(4,2) - t399 * t393 + t400 * t394 + t428;
t397 = -qJDD(1) * pkin(1) - t425 * qJ(2) + t435;
t391 = -t415 * t401 + t440;
t383 = Ifges(4,1) * t400 + Ifges(4,4) * t399 + Ifges(4,5) * qJD(3);
t382 = Ifges(4,4) * t400 + Ifges(4,2) * t399 + Ifges(4,6) * qJD(3);
t381 = Ifges(4,5) * t400 + Ifges(4,6) * t399 + Ifges(4,3) * qJD(3);
t365 = Ifges(5,5) * t385 + Ifges(5,6) * t384 + Ifges(5,3) * t414;
t346 = Ifges(6,5) * t370 + Ifges(6,6) * t369 + Ifges(6,3) * t409;
t321 = t445 * t425 * mrSges(3,3) + m(3) * t397 + t434 * qJDD(1) + t426;
t320 = mrSges(6,2) * t333 - mrSges(6,3) * t328 + Ifges(6,1) * t341 + Ifges(6,4) * t340 + Ifges(6,5) * t408 + t369 * t346 - t409 * t347;
t319 = -mrSges(6,1) * t333 + mrSges(6,3) * t329 + Ifges(6,4) * t341 + Ifges(6,2) * t340 + Ifges(6,6) * t408 - t370 * t346 + t409 * t348;
t305 = mrSges(5,2) * t355 - mrSges(5,3) * t335 + Ifges(5,1) * t360 + Ifges(5,4) * t359 + Ifges(5,5) * t411 - pkin(8) * t318 - t417 * t319 + t421 * t320 + t384 * t365 - t414 * t366;
t304 = -mrSges(5,1) * t355 + mrSges(5,3) * t336 + Ifges(5,4) * t360 + Ifges(5,2) * t359 + Ifges(5,6) * t411 - pkin(4) * t430 + pkin(8) * t437 + t421 * t319 + t417 * t320 - t385 * t365 + t414 * t367;
t303 = mrSges(4,2) * t388 - mrSges(4,3) * t361 + Ifges(4,1) * t390 + Ifges(4,4) * t389 + Ifges(4,5) * qJDD(3) - pkin(7) * t311 - qJD(3) * t382 - t418 * t304 + t422 * t305 + t399 * t381;
t302 = -mrSges(4,1) * t388 + mrSges(4,3) * t362 + Ifges(4,4) * t390 + Ifges(4,2) * t389 + Ifges(4,6) * qJDD(3) - pkin(3) * t428 + pkin(7) * t438 + qJD(3) * t383 + t422 * t304 + t418 * t305 - t400 * t381;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t441 - mrSges(2,2) * t436 + t415 * (mrSges(3,2) * t397 - mrSges(3,3) * t391 + t423 * t303 - t419 * t302 - pkin(6) * t446 + (Ifges(3,1) * t415 + Ifges(3,4) * t416) * qJDD(1)) + t416 * (-mrSges(3,1) * t397 + mrSges(3,3) * t392 + t419 * t303 + t423 * t302 - pkin(2) * t426 + pkin(6) * t439 + (Ifges(3,4) * t415 + Ifges(3,2) * t416) * qJDD(1)) - pkin(1) * t321 + qJ(2) * ((m(3) * t392 + t431 * t416 + t439) * t416 + (-m(3) * t391 + t431 * t415 - t446) * t415); t321; mrSges(4,1) * t361 - mrSges(4,2) * t362 + Ifges(4,5) * t390 + Ifges(4,6) * t389 + Ifges(4,3) * qJDD(3) + pkin(3) * t311 + t400 * t382 - t399 * t383 + t427; t427; t429;];
tauJ = t1;
