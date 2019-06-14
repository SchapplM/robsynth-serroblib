% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPPPRR3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
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
% Datum: 2019-05-05 13:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPPPRR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:40:57
% EndTime: 2019-05-05 13:40:58
% DurationCPUTime: 1.24s
% Computational Cost: add. (10783->209), mult. (21784->263), div. (0->0), fcn. (12472->10), ass. (0->100)
t392 = cos(pkin(10));
t386 = t392 ^ 2;
t401 = qJD(1) ^ 2;
t390 = sin(pkin(10));
t395 = sin(qJ(5));
t398 = cos(qJ(5));
t410 = t390 * t395 - t392 * t398;
t371 = t410 * qJD(1);
t411 = -t390 * t398 - t392 * t395;
t372 = t411 * qJD(1);
t420 = t372 * qJD(5);
t357 = t410 * qJDD(1) - t420;
t428 = -pkin(1) - pkin(2);
t427 = t390 * mrSges(5,2);
t426 = t392 * mrSges(5,1);
t425 = t386 * t401;
t396 = sin(qJ(1));
t399 = cos(qJ(1));
t413 = -t399 * g(1) - t396 * g(2);
t407 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t413;
t364 = t428 * t401 + t407;
t417 = t396 * g(1) - t399 * g(2);
t406 = -t401 * qJ(2) + qJDD(2) - t417;
t367 = t428 * qJDD(1) + t406;
t391 = sin(pkin(9));
t393 = cos(pkin(9));
t346 = t393 * t364 + t391 * t367;
t343 = -t401 * pkin(3) - qJDD(1) * qJ(4) + t346;
t388 = g(3) + qJDD(3);
t419 = qJD(1) * qJD(4);
t423 = t392 * t388 + 0.2e1 * t390 * t419;
t329 = (pkin(4) * t392 * t401 + pkin(7) * qJDD(1) - t343) * t390 + t423;
t333 = t390 * t388 + (t343 - 0.2e1 * t419) * t392;
t418 = t392 * qJDD(1);
t330 = -pkin(4) * t425 - pkin(7) * t418 + t333;
t326 = t395 * t329 + t398 * t330;
t353 = -t371 * mrSges(6,1) + t372 * mrSges(6,2);
t366 = qJD(5) * mrSges(6,1) - t372 * mrSges(6,3);
t356 = -t371 * pkin(5) - t372 * pkin(8);
t400 = qJD(5) ^ 2;
t323 = -t400 * pkin(5) + qJDD(5) * pkin(8) + t371 * t356 + t326;
t385 = t390 ^ 2;
t345 = -t391 * t364 + t393 * t367;
t408 = qJDD(1) * pkin(3) + qJDD(4) - t345;
t331 = pkin(4) * t418 + (-qJ(4) + (-t385 - t386) * pkin(7)) * t401 + t408;
t421 = t371 * qJD(5);
t358 = t411 * qJDD(1) + t421;
t324 = (-t358 - t421) * pkin(8) + (-t357 + t420) * pkin(5) + t331;
t394 = sin(qJ(6));
t397 = cos(qJ(6));
t320 = -t394 * t323 + t397 * t324;
t361 = t397 * qJD(5) - t394 * t372;
t340 = t361 * qJD(6) + t394 * qJDD(5) + t397 * t358;
t362 = t394 * qJD(5) + t397 * t372;
t344 = -t361 * mrSges(7,1) + t362 * mrSges(7,2);
t369 = qJD(6) - t371;
t347 = -t369 * mrSges(7,2) + t361 * mrSges(7,3);
t355 = qJDD(6) - t357;
t318 = m(7) * t320 + t355 * mrSges(7,1) - t340 * mrSges(7,3) - t362 * t344 + t369 * t347;
t321 = t397 * t323 + t394 * t324;
t339 = -t362 * qJD(6) + t397 * qJDD(5) - t394 * t358;
t348 = t369 * mrSges(7,1) - t362 * mrSges(7,3);
t319 = m(7) * t321 - t355 * mrSges(7,2) + t339 * mrSges(7,3) + t361 * t344 - t369 * t348;
t414 = -t394 * t318 + t397 * t319;
t309 = m(6) * t326 - qJDD(5) * mrSges(6,2) + t357 * mrSges(6,3) - qJD(5) * t366 + t371 * t353 + t414;
t325 = t398 * t329 - t395 * t330;
t365 = -qJD(5) * mrSges(6,2) + t371 * mrSges(6,3);
t322 = -qJDD(5) * pkin(5) - t400 * pkin(8) + t372 * t356 - t325;
t405 = -m(7) * t322 + t339 * mrSges(7,1) - t340 * mrSges(7,2) + t361 * t347 - t362 * t348;
t314 = m(6) * t325 + qJDD(5) * mrSges(6,1) - t358 * mrSges(6,3) + qJD(5) * t365 - t372 * t353 + t405;
t424 = t395 * t309 + t398 * t314;
t310 = t397 * t318 + t394 * t319;
t332 = -t390 * t343 + t423;
t409 = mrSges(5,3) * qJDD(1) + t401 * (t426 - t427);
t303 = m(5) * t332 + t409 * t390 + t424;
t415 = t398 * t309 - t395 * t314;
t304 = m(5) * t333 - t409 * t392 + t415;
t416 = -t390 * t303 + t392 * t304;
t299 = m(4) * t346 - t401 * mrSges(4,1) + qJDD(1) * mrSges(4,2) + t416;
t342 = -t401 * qJ(4) + t408;
t404 = m(6) * t331 - t357 * mrSges(6,1) + t358 * mrSges(6,2) - t371 * t365 + t372 * t366 + t310;
t403 = -m(5) * t342 + qJDD(1) * t427 - t404 + (t385 * t401 + t425) * mrSges(5,3);
t305 = t403 + (-mrSges(4,1) - t426) * qJDD(1) + m(4) * t345 - t401 * mrSges(4,2);
t412 = t391 * t299 + t393 * t305;
t335 = Ifges(7,4) * t362 + Ifges(7,2) * t361 + Ifges(7,6) * t369;
t336 = Ifges(7,1) * t362 + Ifges(7,4) * t361 + Ifges(7,5) * t369;
t402 = mrSges(7,1) * t320 - mrSges(7,2) * t321 + Ifges(7,5) * t340 + Ifges(7,6) * t339 + Ifges(7,3) * t355 + t362 * t335 - t361 * t336;
t370 = -qJDD(1) * pkin(1) + t406;
t368 = -t401 * pkin(1) + t407;
t351 = Ifges(6,1) * t372 + Ifges(6,4) * t371 + Ifges(6,5) * qJD(5);
t350 = Ifges(6,4) * t372 + Ifges(6,2) * t371 + Ifges(6,6) * qJD(5);
t349 = Ifges(6,5) * t372 + Ifges(6,6) * t371 + Ifges(6,3) * qJD(5);
t334 = Ifges(7,5) * t362 + Ifges(7,6) * t361 + Ifges(7,3) * t369;
t312 = mrSges(7,2) * t322 - mrSges(7,3) * t320 + Ifges(7,1) * t340 + Ifges(7,4) * t339 + Ifges(7,5) * t355 + t361 * t334 - t369 * t335;
t311 = -mrSges(7,1) * t322 + mrSges(7,3) * t321 + Ifges(7,4) * t340 + Ifges(7,2) * t339 + Ifges(7,6) * t355 - t362 * t334 + t369 * t336;
t306 = mrSges(5,1) * t418 - t403;
t301 = -mrSges(6,1) * t331 + mrSges(6,3) * t326 + Ifges(6,4) * t358 + Ifges(6,2) * t357 + Ifges(6,6) * qJDD(5) - pkin(5) * t310 + qJD(5) * t351 - t372 * t349 - t402;
t300 = mrSges(6,2) * t331 - mrSges(6,3) * t325 + Ifges(6,1) * t358 + Ifges(6,4) * t357 + Ifges(6,5) * qJDD(5) - pkin(8) * t310 - qJD(5) * t350 - t394 * t311 + t397 * t312 + t371 * t349;
t298 = m(3) * t370 - qJDD(1) * mrSges(3,1) - t401 * mrSges(3,3) + t412;
t1 = [-pkin(1) * t298 + qJ(2) * (m(3) * t368 - t401 * mrSges(3,1) + t393 * t299 - t391 * t305) + mrSges(2,1) * t417 - mrSges(2,2) * t413 - pkin(2) * t412 - mrSges(3,1) * t370 + mrSges(3,3) * t368 - t390 * (mrSges(5,2) * t342 - mrSges(5,3) * t332 - pkin(7) * t424 + t398 * t300 - t395 * t301) - t392 * (-mrSges(5,1) * t342 + mrSges(5,3) * t333 - pkin(4) * t404 + pkin(7) * t415 + t395 * t300 + t398 * t301) + pkin(3) * t306 - qJ(4) * t416 - mrSges(4,1) * t345 + mrSges(4,2) * t346 + (t386 * Ifges(5,2) + qJ(2) * mrSges(3,3) + Ifges(3,2) + Ifges(2,3) + Ifges(4,3) + (Ifges(5,1) * t390 + 0.2e1 * Ifges(5,4) * t392) * t390) * qJDD(1); t298; m(4) * t388 + t392 * t303 + t390 * t304; t306; mrSges(6,1) * t325 - mrSges(6,2) * t326 + Ifges(6,5) * t358 + Ifges(6,6) * t357 + Ifges(6,3) * qJDD(5) + pkin(5) * t405 + pkin(8) * t414 + t397 * t311 + t394 * t312 + t372 * t350 - t371 * t351; t402;];
tauJ  = t1;
