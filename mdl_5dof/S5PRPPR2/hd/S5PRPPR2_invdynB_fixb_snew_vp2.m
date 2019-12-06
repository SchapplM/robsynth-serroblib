% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRPPR2
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRPPR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:24:11
% EndTime: 2019-12-05 15:24:14
% DurationCPUTime: 1.96s
% Computational Cost: add. (20355->195), mult. (38418->248), div. (0->0), fcn. (24486->10), ass. (0->90)
t418 = qJD(2) ^ 2;
t411 = cos(pkin(9));
t441 = pkin(4) * t411;
t408 = sin(pkin(9));
t440 = mrSges(5,2) * t408;
t406 = t411 ^ 2;
t439 = t406 * t418;
t410 = sin(pkin(7));
t413 = cos(pkin(7));
t396 = -g(1) * t413 - g(2) * t410;
t407 = -g(3) + qJDD(1);
t415 = sin(qJ(2));
t417 = cos(qJ(2));
t385 = -t396 * t415 + t417 * t407;
t381 = qJDD(2) * pkin(2) + t385;
t386 = t417 * t396 + t415 * t407;
t382 = -pkin(2) * t418 + t386;
t409 = sin(pkin(8));
t412 = cos(pkin(8));
t369 = t409 * t381 + t412 * t382;
t367 = -pkin(3) * t418 + qJDD(2) * qJ(4) + t369;
t395 = g(1) * t410 - t413 * g(2);
t394 = qJDD(3) - t395;
t435 = qJD(2) * qJD(4);
t437 = t411 * t394 - 0.2e1 * t408 * t435;
t360 = (-pkin(6) * qJDD(2) + t418 * t441 - t367) * t408 + t437;
t363 = t408 * t394 + (t367 + 0.2e1 * t435) * t411;
t434 = qJDD(2) * t411;
t361 = -pkin(4) * t439 + pkin(6) * t434 + t363;
t414 = sin(qJ(5));
t416 = cos(qJ(5));
t358 = t360 * t416 - t361 * t414;
t422 = -t408 * t414 + t411 * t416;
t387 = t422 * qJD(2);
t423 = t408 * t416 + t411 * t414;
t388 = t423 * qJD(2);
t374 = -mrSges(6,1) * t387 + mrSges(6,2) * t388;
t377 = qJD(5) * t387 + t423 * qJDD(2);
t383 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t387;
t356 = m(6) * t358 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t377 + qJD(5) * t383 - t374 * t388;
t359 = t360 * t414 + t361 * t416;
t376 = -qJD(5) * t388 + t422 * qJDD(2);
t384 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t388;
t357 = m(6) * t359 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t376 - qJD(5) * t384 + t374 * t387;
t348 = t416 * t356 + t414 * t357;
t362 = -t367 * t408 + t437;
t421 = mrSges(5,3) * qJDD(2) + t418 * (-mrSges(5,1) * t411 + t440);
t346 = m(5) * t362 - t421 * t408 + t348;
t429 = -t356 * t414 + t416 * t357;
t347 = m(5) * t363 + t421 * t411 + t429;
t430 = -t346 * t408 + t411 * t347;
t339 = m(4) * t369 - mrSges(4,1) * t418 - qJDD(2) * mrSges(4,2) + t430;
t368 = t381 * t412 - t409 * t382;
t424 = qJDD(4) - t368;
t366 = -qJDD(2) * pkin(3) - qJ(4) * t418 + t424;
t405 = t408 ^ 2;
t364 = (-pkin(3) - t441) * qJDD(2) + (-qJ(4) + (-t405 - t406) * pkin(6)) * t418 + t424;
t420 = m(6) * t364 - t376 * mrSges(6,1) + mrSges(6,2) * t377 - t387 * t383 + t384 * t388;
t419 = -m(5) * t366 + mrSges(5,1) * t434 - t420 + (t405 * t418 + t439) * mrSges(5,3);
t352 = m(4) * t368 - mrSges(4,2) * t418 + (mrSges(4,1) - t440) * qJDD(2) + t419;
t335 = t409 * t339 + t412 * t352;
t332 = m(3) * t385 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t418 + t335;
t431 = t412 * t339 - t409 * t352;
t333 = m(3) * t386 - mrSges(3,1) * t418 - qJDD(2) * mrSges(3,2) + t431;
t432 = -t332 * t415 + t417 * t333;
t326 = m(2) * t396 + t432;
t342 = t411 * t346 + t408 * t347;
t428 = m(4) * t394 + t342;
t341 = (m(2) + m(3)) * t395 - t428;
t438 = t410 * t326 + t413 * t341;
t327 = t417 * t332 + t415 * t333;
t425 = Ifges(5,5) * t408 + Ifges(5,6) * t411;
t436 = t418 * t425;
t433 = t413 * t326 - t341 * t410;
t427 = Ifges(5,1) * t408 + Ifges(5,4) * t411;
t426 = Ifges(5,4) * t408 + Ifges(5,2) * t411;
t372 = Ifges(6,1) * t388 + Ifges(6,4) * t387 + Ifges(6,5) * qJD(5);
t371 = Ifges(6,4) * t388 + Ifges(6,2) * t387 + Ifges(6,6) * qJD(5);
t370 = Ifges(6,5) * t388 + Ifges(6,6) * t387 + Ifges(6,3) * qJD(5);
t350 = mrSges(6,2) * t364 - mrSges(6,3) * t358 + Ifges(6,1) * t377 + Ifges(6,4) * t376 + Ifges(6,5) * qJDD(5) - qJD(5) * t371 + t370 * t387;
t349 = -mrSges(6,1) * t364 + mrSges(6,3) * t359 + Ifges(6,4) * t377 + Ifges(6,2) * t376 + Ifges(6,6) * qJDD(5) + qJD(5) * t372 - t370 * t388;
t336 = mrSges(5,2) * t366 - mrSges(5,3) * t362 - pkin(6) * t348 + t427 * qJDD(2) - t349 * t414 + t350 * t416 + t411 * t436;
t334 = -mrSges(5,1) * t366 + mrSges(5,3) * t363 - pkin(4) * t420 + pkin(6) * t429 + t426 * qJDD(2) + t416 * t349 + t414 * t350 - t408 * t436;
t328 = -mrSges(4,1) * t394 - mrSges(5,1) * t362 - mrSges(6,1) * t358 + mrSges(5,2) * t363 + mrSges(6,2) * t359 + mrSges(4,3) * t369 - Ifges(6,5) * t377 - Ifges(6,6) * t376 - Ifges(6,3) * qJDD(5) - pkin(3) * t342 - pkin(4) * t348 - t388 * t371 + t387 * t372 + (Ifges(4,6) - t425) * qJDD(2) + (-t408 * t426 + t411 * t427 + Ifges(4,5)) * t418;
t323 = mrSges(4,2) * t394 - mrSges(4,3) * t368 + Ifges(4,5) * qJDD(2) - Ifges(4,6) * t418 - qJ(4) * t342 - t334 * t408 + t336 * t411;
t322 = -mrSges(3,2) * t395 - mrSges(3,3) * t385 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t418 - qJ(3) * t335 + t323 * t412 - t328 * t409;
t321 = -pkin(1) * t327 + mrSges(2,3) * t396 - pkin(2) * t335 - mrSges(3,1) * t385 + mrSges(3,2) * t386 - t408 * t336 - t411 * t334 - pkin(3) * t419 - qJ(4) * t430 - mrSges(4,1) * t368 + mrSges(4,2) * t369 - mrSges(2,1) * t407 + (pkin(3) * t440 - Ifges(3,3) - Ifges(4,3)) * qJDD(2);
t320 = mrSges(3,1) * t395 + mrSges(3,3) * t386 + t418 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t428 + qJ(3) * t431 + t409 * t323 + t412 * t328;
t319 = mrSges(2,2) * t407 - mrSges(2,3) * t395 - pkin(5) * t327 - t320 * t415 + t322 * t417;
t1 = [-m(1) * g(1) + t433; -m(1) * g(2) + t438; -m(1) * g(3) + m(2) * t407 + t327; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t438 + t413 * t319 - t410 * t321; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t433 + t410 * t319 + t413 * t321; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t395 - mrSges(2,2) * t396 + t415 * t322 + t417 * t320 + pkin(1) * (m(3) * t395 - t428) + pkin(5) * t432;];
tauB = t1;
