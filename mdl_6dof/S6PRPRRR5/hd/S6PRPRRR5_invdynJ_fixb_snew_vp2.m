% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 01:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRPRRR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 01:15:59
% EndTime: 2019-05-05 01:16:01
% DurationCPUTime: 1.52s
% Computational Cost: add. (11299->235), mult. (21385->299), div. (0->0), fcn. (14559->12), ass. (0->102)
t429 = sin(pkin(11));
t431 = cos(pkin(11));
t412 = g(1) * t429 - g(2) * t431;
t426 = -g(3) + qJDD(1);
t430 = sin(pkin(6));
t432 = cos(pkin(6));
t464 = t412 * t432 + t426 * t430;
t413 = -g(1) * t431 - g(2) * t429;
t436 = sin(qJ(2));
t440 = cos(qJ(2));
t386 = t440 * t413 + t464 * t436;
t452 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t386;
t434 = sin(qJ(5));
t435 = sin(qJ(4));
t438 = cos(qJ(5));
t439 = cos(qJ(4));
t401 = (t434 * t439 + t435 * t438) * qJD(2);
t385 = -t436 * t413 + t464 * t440;
t463 = -pkin(2) - pkin(8);
t441 = qJD(2) ^ 2;
t462 = pkin(4) * t441;
t445 = -t441 * qJ(3) + qJDD(3) - t385;
t371 = t463 * qJDD(2) + t445;
t368 = t439 * t371;
t395 = -t412 * t430 + t426 * t432;
t457 = qJD(2) * qJD(4);
t411 = qJDD(2) * t439 - t435 * t457;
t352 = (qJDD(4) * pkin(4)) - t411 * pkin(9) + t368 + (-pkin(9) * t457 - t439 * t462 - t395) * t435;
t360 = t435 * t371 + t439 * t395;
t410 = -qJDD(2) * t435 - t439 * t457;
t458 = qJD(2) * t439;
t417 = (qJD(4) * pkin(4)) - pkin(9) * t458;
t425 = t435 ^ 2;
t353 = pkin(9) * t410 - qJD(4) * t417 - t425 * t462 + t360;
t348 = t434 * t352 + t438 * t353;
t402 = (-t434 * t435 + t438 * t439) * qJD(2);
t377 = -qJD(5) * t402 + t410 * t438 - t411 * t434;
t388 = mrSges(6,1) * t401 + mrSges(6,2) * t402;
t422 = qJD(4) + qJD(5);
t394 = mrSges(6,1) * t422 - mrSges(6,3) * t402;
t421 = qJDD(4) + qJDD(5);
t389 = pkin(5) * t401 - pkin(10) * t402;
t420 = t422 ^ 2;
t345 = -pkin(5) * t420 + pkin(10) * t421 - t389 * t401 + t348;
t358 = -t410 * pkin(4) + t417 * t458 + (-pkin(9) * t425 + t463) * t441 + t452;
t378 = -qJD(5) * t401 + t410 * t434 + t411 * t438;
t349 = (t401 * t422 - t378) * pkin(10) + (t402 * t422 - t377) * pkin(5) + t358;
t433 = sin(qJ(6));
t437 = cos(qJ(6));
t342 = -t345 * t433 + t349 * t437;
t390 = -t402 * t433 + t422 * t437;
t356 = qJD(6) * t390 + t378 * t437 + t421 * t433;
t391 = t402 * t437 + t422 * t433;
t365 = -mrSges(7,1) * t390 + mrSges(7,2) * t391;
t375 = qJDD(6) - t377;
t396 = qJD(6) + t401;
t380 = -mrSges(7,2) * t396 + mrSges(7,3) * t390;
t339 = m(7) * t342 + mrSges(7,1) * t375 - mrSges(7,3) * t356 - t365 * t391 + t380 * t396;
t343 = t345 * t437 + t349 * t433;
t355 = -qJD(6) * t391 - t378 * t433 + t421 * t437;
t381 = mrSges(7,1) * t396 - mrSges(7,3) * t391;
t340 = m(7) * t343 - mrSges(7,2) * t375 + mrSges(7,3) * t355 + t365 * t390 - t381 * t396;
t454 = -t339 * t433 + t437 * t340;
t327 = m(6) * t348 - mrSges(6,2) * t421 + mrSges(6,3) * t377 - t388 * t401 - t394 * t422 + t454;
t347 = t352 * t438 - t353 * t434;
t393 = -mrSges(6,2) * t422 - mrSges(6,3) * t401;
t344 = -pkin(5) * t421 - pkin(10) * t420 + t389 * t402 - t347;
t447 = -m(7) * t344 + t355 * mrSges(7,1) - mrSges(7,2) * t356 + t390 * t380 - t381 * t391;
t335 = m(6) * t347 + mrSges(6,1) * t421 - mrSges(6,3) * t378 - t388 * t402 + t393 * t422 + t447;
t324 = t434 * t327 + t438 * t335;
t329 = t437 * t339 + t433 * t340;
t459 = qJD(2) * t435;
t455 = t438 * t327 - t335 * t434;
t359 = -t435 * t395 + t368;
t409 = (mrSges(5,1) * t435 + mrSges(5,2) * t439) * qJD(2);
t414 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t459;
t322 = m(5) * t359 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t411 + qJD(4) * t414 - t409 * t458 + t324;
t415 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t458;
t323 = m(5) * t360 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t410 - qJD(4) * t415 - t409 * t459 + t455;
t451 = t439 * t322 + t435 * t323;
t379 = -qJDD(2) * pkin(2) + t445;
t448 = -m(4) * t379 + (t441 * mrSges(4,3)) - t451;
t446 = m(6) * t358 - mrSges(6,1) * t377 + t378 * mrSges(6,2) + t393 * t401 + t402 * t394 + t329;
t361 = Ifges(7,5) * t391 + Ifges(7,6) * t390 + Ifges(7,3) * t396;
t363 = Ifges(7,1) * t391 + Ifges(7,4) * t390 + Ifges(7,5) * t396;
t332 = -mrSges(7,1) * t344 + mrSges(7,3) * t343 + Ifges(7,4) * t356 + Ifges(7,2) * t355 + Ifges(7,6) * t375 - t361 * t391 + t363 * t396;
t362 = Ifges(7,4) * t391 + Ifges(7,2) * t390 + Ifges(7,6) * t396;
t333 = mrSges(7,2) * t344 - mrSges(7,3) * t342 + Ifges(7,1) * t356 + Ifges(7,4) * t355 + Ifges(7,5) * t375 + t361 * t390 - t362 * t396;
t383 = Ifges(6,4) * t402 - Ifges(6,2) * t401 + Ifges(6,6) * t422;
t384 = Ifges(6,1) * t402 - Ifges(6,4) * t401 + Ifges(6,5) * t422;
t444 = mrSges(6,1) * t347 - mrSges(6,2) * t348 + Ifges(6,5) * t378 + Ifges(6,6) * t377 + (Ifges(6,3) * t421) + pkin(5) * t447 + pkin(10) * t454 + t437 * t332 + t433 * t333 + t402 * t383 + t401 * t384;
t443 = mrSges(7,1) * t342 - mrSges(7,2) * t343 + Ifges(7,5) * t356 + Ifges(7,6) * t355 + Ifges(7,3) * t375 + t362 * t391 - t363 * t390;
t370 = t463 * t441 + t452;
t376 = t441 * pkin(2) - t452;
t442 = -m(4) * t376 + m(5) * t370 - mrSges(5,1) * t410 + (t441 * mrSges(4,2)) + t411 * mrSges(5,2) + qJDD(2) * mrSges(4,3) + t414 * t459 + t415 * t458 + t446;
t400 = (Ifges(5,5) * qJD(4)) + (Ifges(5,1) * t439 - Ifges(5,4) * t435) * qJD(2);
t399 = (Ifges(5,6) * qJD(4)) + (Ifges(5,4) * t439 - Ifges(5,2) * t435) * qJD(2);
t382 = Ifges(6,5) * t402 - Ifges(6,6) * t401 + Ifges(6,3) * t422;
t321 = -mrSges(6,1) * t358 + mrSges(6,3) * t348 + Ifges(6,4) * t378 + Ifges(6,2) * t377 + Ifges(6,6) * t421 - pkin(5) * t329 - t382 * t402 + t384 * t422 - t443;
t320 = mrSges(6,2) * t358 - mrSges(6,3) * t347 + Ifges(6,1) * t378 + Ifges(6,4) * t377 + Ifges(6,5) * t421 - pkin(10) * t329 - t332 * t433 + t333 * t437 - t382 * t401 - t383 * t422;
t319 = qJDD(2) * mrSges(4,2) - t448;
t1 = [m(2) * t426 + t432 * (-t435 * t322 + t439 * t323 + (m(3) + m(4)) * t395) + (t436 * (m(3) * t386 - (mrSges(3,1) * t441) - qJDD(2) * mrSges(3,2) + t442) + t440 * (m(3) * t385 - t441 * mrSges(3,2) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + t448)) * t430; mrSges(3,1) * t385 - mrSges(3,2) * t386 + mrSges(4,2) * t379 - mrSges(4,3) * t376 + t439 * (mrSges(5,2) * t370 - mrSges(5,3) * t359 + Ifges(5,1) * t411 + Ifges(5,4) * t410 + (Ifges(5,5) * qJDD(4)) - pkin(9) * t324 - qJD(4) * t399 + t320 * t438 - t321 * t434) - t435 * (-mrSges(5,1) * t370 + mrSges(5,3) * t360 + Ifges(5,4) * t411 + Ifges(5,2) * t410 + (Ifges(5,6) * qJDD(4)) - pkin(4) * t446 + pkin(9) * t455 + qJD(4) * t400 + t434 * t320 + t438 * t321) - pkin(8) * t451 - pkin(2) * t319 + qJ(3) * t442 + (Ifges(3,3) + Ifges(4,1)) * qJDD(2); t319; t444 + (t399 * t439 + t400 * t435) * qJD(2) + (Ifges(5,3) * qJDD(4)) + Ifges(5,6) * t410 + Ifges(5,5) * t411 + mrSges(5,1) * t359 - mrSges(5,2) * t360 + pkin(4) * t324; t444; t443;];
tauJ  = t1;
