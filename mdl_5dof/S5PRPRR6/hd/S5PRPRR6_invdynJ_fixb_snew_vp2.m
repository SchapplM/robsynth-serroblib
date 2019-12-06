% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRPRR6
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRPRR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR6_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR6_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR6_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR6_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR6_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:56:48
% EndTime: 2019-12-05 15:56:50
% DurationCPUTime: 1.04s
% Computational Cost: add. (6368->189), mult. (13866->247), div. (0->0), fcn. (10065->12), ass. (0->94)
t363 = sin(pkin(9));
t366 = cos(pkin(9));
t351 = t363 * g(1) - t366 * g(2);
t361 = -g(3) + qJDD(1);
t364 = sin(pkin(5));
t367 = cos(pkin(5));
t400 = t351 * t367 + t361 * t364;
t375 = qJD(2) ^ 2;
t362 = sin(pkin(10));
t365 = cos(pkin(10));
t369 = sin(qJ(4));
t372 = cos(qJ(4));
t382 = t362 * t369 - t365 * t372;
t344 = t382 * qJD(2);
t352 = -t366 * g(1) - t363 * g(2);
t370 = sin(qJ(2));
t373 = cos(qJ(2));
t327 = -t370 * t352 + t400 * t373;
t383 = t362 * t372 + t365 * t369;
t345 = t383 * qJD(2);
t390 = t345 * qJD(4);
t334 = -t382 * qJDD(2) - t390;
t360 = t365 ^ 2;
t399 = 0.2e1 * t365;
t398 = pkin(3) * t365;
t397 = mrSges(4,2) * t362;
t395 = t360 * t375;
t328 = t373 * t352 + t400 * t370;
t323 = -t375 * pkin(2) + qJDD(2) * qJ(3) + t328;
t342 = -t364 * t351 + t367 * t361;
t389 = qJD(2) * qJD(3);
t392 = t365 * t342 - 0.2e1 * t362 * t389;
t306 = (-pkin(7) * qJDD(2) + t375 * t398 - t323) * t362 + t392;
t309 = t365 * t323 + t362 * t342 + t389 * t399;
t388 = t365 * qJDD(2);
t307 = -pkin(3) * t395 + pkin(7) * t388 + t309;
t302 = t369 * t306 + t372 * t307;
t330 = mrSges(5,1) * t344 + mrSges(5,2) * t345;
t341 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t345;
t333 = pkin(4) * t344 - pkin(8) * t345;
t374 = qJD(4) ^ 2;
t300 = -t374 * pkin(4) + qJDD(4) * pkin(8) - t333 * t344 + t302;
t359 = t362 ^ 2;
t379 = qJDD(3) - t327;
t317 = (-pkin(2) - t398) * qJDD(2) + (-qJ(3) + (-t359 - t360) * pkin(7)) * t375 + t379;
t391 = t344 * qJD(4);
t335 = t383 * qJDD(2) - t391;
t303 = (-t335 + t391) * pkin(8) + (-t334 + t390) * pkin(4) + t317;
t368 = sin(qJ(5));
t371 = cos(qJ(5));
t297 = -t368 * t300 + t371 * t303;
t336 = t371 * qJD(4) - t368 * t345;
t316 = t336 * qJD(5) + t368 * qJDD(4) + t371 * t335;
t337 = t368 * qJD(4) + t371 * t345;
t318 = -mrSges(6,1) * t336 + mrSges(6,2) * t337;
t343 = qJD(5) + t344;
t321 = -mrSges(6,2) * t343 + mrSges(6,3) * t336;
t332 = qJDD(5) - t334;
t295 = m(6) * t297 + mrSges(6,1) * t332 - mrSges(6,3) * t316 - t318 * t337 + t321 * t343;
t298 = t371 * t300 + t368 * t303;
t315 = -t337 * qJD(5) + t371 * qJDD(4) - t368 * t335;
t322 = mrSges(6,1) * t343 - mrSges(6,3) * t337;
t296 = m(6) * t298 - mrSges(6,2) * t332 + mrSges(6,3) * t315 + t318 * t336 - t322 * t343;
t385 = -t368 * t295 + t371 * t296;
t285 = m(5) * t302 - qJDD(4) * mrSges(5,2) + t334 * mrSges(5,3) - qJD(4) * t341 - t344 * t330 + t385;
t301 = t372 * t306 - t369 * t307;
t340 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t344;
t299 = -qJDD(4) * pkin(4) - t374 * pkin(8) + t345 * t333 - t301;
t380 = -m(6) * t299 + t315 * mrSges(6,1) - mrSges(6,2) * t316 + t336 * t321 - t322 * t337;
t291 = m(5) * t301 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t335 + qJD(4) * t340 - t330 * t345 + t380;
t393 = t369 * t285 + t372 * t291;
t287 = t371 * t295 + t368 * t296;
t308 = -t362 * t323 + t392;
t381 = mrSges(4,3) * qJDD(2) + t375 * (-t365 * mrSges(4,1) + t397);
t281 = m(4) * t308 - t381 * t362 + t393;
t386 = t372 * t285 - t369 * t291;
t282 = m(4) * t309 + t381 * t365 + t386;
t387 = -t362 * t281 + t365 * t282;
t378 = m(5) * t317 - t334 * mrSges(5,1) + t335 * mrSges(5,2) + t344 * t340 + t345 * t341 + t287;
t320 = -qJDD(2) * pkin(2) - t375 * qJ(3) + t379;
t377 = -m(4) * t320 + mrSges(4,1) * t388 - t378 + (t359 * t375 + t395) * mrSges(4,3);
t311 = Ifges(6,4) * t337 + Ifges(6,2) * t336 + Ifges(6,6) * t343;
t312 = Ifges(6,1) * t337 + Ifges(6,4) * t336 + Ifges(6,5) * t343;
t376 = mrSges(6,1) * t297 - mrSges(6,2) * t298 + Ifges(6,5) * t316 + Ifges(6,6) * t315 + Ifges(6,3) * t332 + t311 * t337 - t312 * t336;
t326 = Ifges(5,1) * t345 - Ifges(5,4) * t344 + Ifges(5,5) * qJD(4);
t325 = Ifges(5,4) * t345 - Ifges(5,2) * t344 + Ifges(5,6) * qJD(4);
t324 = Ifges(5,5) * t345 - Ifges(5,6) * t344 + Ifges(5,3) * qJD(4);
t310 = Ifges(6,5) * t337 + Ifges(6,6) * t336 + Ifges(6,3) * t343;
t289 = mrSges(6,2) * t299 - mrSges(6,3) * t297 + Ifges(6,1) * t316 + Ifges(6,4) * t315 + Ifges(6,5) * t332 + t310 * t336 - t311 * t343;
t288 = -mrSges(6,1) * t299 + mrSges(6,3) * t298 + Ifges(6,4) * t316 + Ifges(6,2) * t315 + Ifges(6,6) * t332 - t310 * t337 + t312 * t343;
t286 = qJDD(2) * t397 - t377;
t279 = -mrSges(5,1) * t317 + mrSges(5,3) * t302 + Ifges(5,4) * t335 + Ifges(5,2) * t334 + Ifges(5,6) * qJDD(4) - pkin(4) * t287 + qJD(4) * t326 - t324 * t345 - t376;
t278 = mrSges(5,2) * t317 - mrSges(5,3) * t301 + Ifges(5,1) * t335 + Ifges(5,4) * t334 + Ifges(5,5) * qJDD(4) - pkin(8) * t287 - qJD(4) * t325 - t368 * t288 + t371 * t289 - t344 * t324;
t1 = [m(2) * t361 + t367 * (m(3) * t342 + t365 * t281 + t362 * t282) + (t370 * (m(3) * t328 - t375 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t387) + t373 * ((mrSges(3,1) - t397) * qJDD(2) + t377 - t375 * mrSges(3,2) + m(3) * t327)) * t364; mrSges(3,1) * t327 - mrSges(3,2) * t328 + t362 * (mrSges(4,2) * t320 - mrSges(4,3) * t308 - pkin(7) * t393 + t372 * t278 - t369 * t279) + t365 * (-mrSges(4,1) * t320 + mrSges(4,3) * t309 - pkin(3) * t378 + pkin(7) * t386 + t369 * t278 + t372 * t279) - pkin(2) * t286 + qJ(3) * t387 + (Ifges(4,2) * t360 + Ifges(3,3) + (Ifges(4,1) * t362 + Ifges(4,4) * t399) * t362) * qJDD(2); t286; mrSges(5,1) * t301 - mrSges(5,2) * t302 + Ifges(5,5) * t335 + Ifges(5,6) * t334 + Ifges(5,3) * qJDD(4) + pkin(4) * t380 + pkin(8) * t385 + t371 * t288 + t368 * t289 + t345 * t325 + t344 * t326; t376;];
tauJ = t1;
