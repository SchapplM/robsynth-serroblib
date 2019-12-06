% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRRRP8
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRRP8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP8_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP8_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP8_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP8_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP8_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:58:34
% EndTime: 2019-12-05 16:58:36
% DurationCPUTime: 1.06s
% Computational Cost: add. (3962->192), mult. (7446->238), div. (0->0), fcn. (4819->10), ass. (0->85)
t390 = Ifges(5,1) + Ifges(6,1);
t379 = Ifges(5,4) - Ifges(6,5);
t385 = -Ifges(5,5) - Ifges(6,4);
t389 = Ifges(5,2) + Ifges(6,3);
t377 = Ifges(5,6) - Ifges(6,6);
t351 = sin(qJ(4));
t352 = sin(qJ(3));
t370 = t352 * qJD(2);
t381 = cos(qJ(4));
t331 = -qJD(3) * t381 + t351 * t370;
t354 = cos(qJ(3));
t368 = qJD(2) * qJD(3);
t365 = t354 * t368;
t335 = t352 * qJDD(2) + t365;
t307 = -t331 * qJD(4) + t351 * qJDD(3) + t335 * t381;
t332 = t351 * qJD(3) + t370 * t381;
t311 = t331 * mrSges(6,1) - t332 * mrSges(6,3);
t347 = sin(pkin(9));
t349 = cos(pkin(9));
t338 = -t349 * g(1) - t347 * g(2);
t353 = sin(qJ(2));
t355 = cos(qJ(2));
t337 = t347 * g(1) - t349 * g(2);
t346 = -g(3) + qJDD(1);
t348 = sin(pkin(5));
t350 = cos(pkin(5));
t387 = t337 * t350 + t346 * t348;
t295 = t355 * t338 + t387 * t353;
t357 = qJD(2) ^ 2;
t293 = -t357 * pkin(2) + qJDD(2) * pkin(7) + t295;
t318 = -t348 * t337 + t350 * t346;
t289 = t354 * t293 + t352 * t318;
t334 = (-t354 * pkin(3) - t352 * pkin(8)) * qJD(2);
t356 = qJD(3) ^ 2;
t369 = t354 * qJD(2);
t285 = -t356 * pkin(3) + qJDD(3) * pkin(8) + t334 * t369 + t289;
t294 = -t353 * t338 + t387 * t355;
t292 = -qJDD(2) * pkin(2) - t357 * pkin(7) - t294;
t366 = t352 * t368;
t336 = t354 * qJDD(2) - t366;
t287 = (-t335 - t365) * pkin(8) + (-t336 + t366) * pkin(3) + t292;
t281 = -t351 * t285 + t287 * t381;
t310 = t331 * pkin(4) - t332 * qJ(5);
t328 = qJDD(4) - t336;
t343 = qJD(4) - t369;
t342 = t343 ^ 2;
t279 = -t328 * pkin(4) - t342 * qJ(5) + t332 * t310 + qJDD(5) - t281;
t317 = -t331 * mrSges(6,2) + t343 * mrSges(6,3);
t362 = -m(6) * t279 + t328 * mrSges(6,1) + t343 * t317;
t275 = t307 * mrSges(6,2) + t332 * t311 - t362;
t282 = t381 * t285 + t351 * t287;
t278 = -t342 * pkin(4) + t328 * qJ(5) + 0.2e1 * qJD(5) * t343 - t331 * t310 + t282;
t306 = t332 * qJD(4) - qJDD(3) * t381 + t351 * t335;
t316 = -t343 * mrSges(6,1) + t332 * mrSges(6,2);
t367 = m(6) * t278 + t328 * mrSges(6,3) + t343 * t316;
t373 = t389 * t331 - t379 * t332 - t377 * t343;
t382 = t379 * t331 - t390 * t332 + t385 * t343;
t384 = -Ifges(5,3) - Ifges(6,2);
t388 = -t385 * t307 - t382 * t331 - t377 * t306 - t384 * t328 + mrSges(5,1) * t281 - mrSges(6,1) * t279 - mrSges(5,2) * t282 + mrSges(6,3) * t278 - pkin(4) * t275 + qJ(5) * (-t306 * mrSges(6,2) - t331 * t311 + t367) - t373 * t332;
t380 = -mrSges(5,3) - mrSges(6,2);
t374 = t377 * t331 + t385 * t332 + t384 * t343;
t371 = -t331 * mrSges(5,1) - t332 * mrSges(5,2) - t311;
t333 = (-t354 * mrSges(4,1) + t352 * mrSges(4,2)) * qJD(2);
t339 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t370;
t315 = t343 * mrSges(5,1) - t332 * mrSges(5,3);
t272 = m(5) * t282 - t328 * mrSges(5,2) + t380 * t306 - t343 * t315 + t371 * t331 + t367;
t314 = -t343 * mrSges(5,2) - t331 * mrSges(5,3);
t273 = m(5) * t281 + t328 * mrSges(5,1) + t380 * t307 + t343 * t314 + t371 * t332 + t362;
t363 = t381 * t272 - t351 * t273;
t266 = m(4) * t289 - qJDD(3) * mrSges(4,2) + t336 * mrSges(4,3) - qJD(3) * t339 + t333 * t369 + t363;
t288 = -t352 * t293 + t354 * t318;
t340 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t369;
t284 = -qJDD(3) * pkin(3) - t356 * pkin(8) + t334 * t370 - t288;
t280 = -0.2e1 * qJD(5) * t332 + (t331 * t343 - t307) * qJ(5) + (t332 * t343 + t306) * pkin(4) + t284;
t276 = m(6) * t280 + t306 * mrSges(6,1) - t307 * mrSges(6,3) - t332 * t316 + t331 * t317;
t358 = -m(5) * t284 - t306 * mrSges(5,1) - t307 * mrSges(5,2) - t331 * t314 - t332 * t315 - t276;
t270 = m(4) * t288 + qJDD(3) * mrSges(4,1) - t335 * mrSges(4,3) + qJD(3) * t340 - t333 * t370 + t358;
t364 = t354 * t266 - t352 * t270;
t269 = t351 * t272 + t273 * t381;
t359 = -m(4) * t292 + t336 * mrSges(4,1) - t335 * mrSges(4,2) - t339 * t370 + t340 * t369 - t269;
t323 = Ifges(4,5) * qJD(3) + (t352 * Ifges(4,1) + t354 * Ifges(4,4)) * qJD(2);
t322 = Ifges(4,6) * qJD(3) + (t352 * Ifges(4,4) + t354 * Ifges(4,2)) * qJD(2);
t268 = mrSges(5,2) * t284 + mrSges(6,2) * t279 - mrSges(5,3) * t281 - mrSges(6,3) * t280 - qJ(5) * t276 - t379 * t306 + t390 * t307 - t385 * t328 + t374 * t331 + t373 * t343;
t267 = -mrSges(5,1) * t284 - mrSges(6,1) * t280 + mrSges(6,2) * t278 + mrSges(5,3) * t282 - pkin(4) * t276 - t389 * t306 + t379 * t307 + t377 * t328 + t374 * t332 - t382 * t343;
t1 = [m(2) * t346 + t350 * (m(3) * t318 + t352 * t266 + t354 * t270) + (t353 * (m(3) * t295 - t357 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t364) + t355 * (m(3) * t294 + qJDD(2) * mrSges(3,1) - t357 * mrSges(3,2) + t359)) * t348; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t294 - mrSges(3,2) * t295 + t352 * (mrSges(4,2) * t292 - mrSges(4,3) * t288 + Ifges(4,1) * t335 + Ifges(4,4) * t336 + Ifges(4,5) * qJDD(3) - pkin(8) * t269 - qJD(3) * t322 - t351 * t267 + t381 * t268) + t354 * (-mrSges(4,1) * t292 + mrSges(4,3) * t289 + Ifges(4,4) * t335 + Ifges(4,2) * t336 + Ifges(4,6) * qJDD(3) - pkin(3) * t269 + qJD(3) * t323 - t388) + pkin(2) * t359 + pkin(7) * t364; Ifges(4,5) * t335 + Ifges(4,6) * t336 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t288 - mrSges(4,2) * t289 + t351 * t268 + t381 * t267 + pkin(3) * t358 + pkin(8) * t363 + (t352 * t322 - t354 * t323) * qJD(2); t388; t275;];
tauJ = t1;
