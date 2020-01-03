% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRPP4
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
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
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRPP4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP4_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP4_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP4_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP4_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:23
% EndTime: 2019-12-31 18:14:24
% DurationCPUTime: 0.80s
% Computational Cost: add. (2984->192), mult. (6336->235), div. (0->0), fcn. (3515->6), ass. (0->77)
t379 = Ifges(5,1) + Ifges(6,1);
t373 = Ifges(5,4) - Ifges(6,5);
t372 = Ifges(5,5) + Ifges(6,4);
t378 = -Ifges(5,2) - Ifges(6,3);
t377 = (-Ifges(6,2) - Ifges(5,3));
t371 = Ifges(5,6) - Ifges(6,6);
t347 = sin(qJ(1));
t349 = cos(qJ(1));
t357 = -t349 * g(1) - t347 * g(2);
t353 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t357;
t376 = -2 * qJD(4);
t375 = -pkin(1) - pkin(6);
t374 = -mrSges(5,3) - mrSges(6,2);
t351 = qJD(1) ^ 2;
t360 = t347 * g(1) - t349 * g(2);
t352 = -t351 * qJ(2) + qJDD(2) - t360;
t319 = t375 * qJDD(1) + t352;
t346 = sin(qJ(3));
t348 = cos(qJ(3));
t308 = t346 * g(3) + t348 * t319;
t364 = qJD(1) * qJD(3);
t361 = t346 * t364;
t332 = qJDD(1) * t348 - t361;
t289 = (-t332 - t361) * qJ(4) + (-t346 * t348 * t351 + qJDD(3)) * pkin(3) + t308;
t309 = -g(3) * t348 + t346 * t319;
t331 = -qJDD(1) * t346 - t348 * t364;
t365 = qJD(1) * t348;
t334 = qJD(3) * pkin(3) - qJ(4) * t365;
t343 = t346 ^ 2;
t290 = -pkin(3) * t343 * t351 + qJ(4) * t331 - qJD(3) * t334 + t309;
t344 = sin(pkin(7));
t345 = cos(pkin(7));
t323 = (t344 * t348 + t345 * t346) * qJD(1);
t286 = t344 * t289 + t345 * t290 + t323 * t376;
t306 = -t345 * t331 + t332 * t344;
t366 = qJD(1) * t346;
t324 = -t344 * t366 + t345 * t365;
t316 = (qJD(3) * mrSges(5,1)) - mrSges(5,3) * t324;
t301 = pkin(4) * t323 - qJ(5) * t324;
t350 = qJD(3) ^ 2;
t281 = -pkin(4) * t350 + qJDD(3) * qJ(5) + (2 * qJD(5) * qJD(3)) - t301 * t323 + t286;
t317 = -(qJD(3) * mrSges(6,1)) + mrSges(6,2) * t324;
t362 = m(6) * t281 + qJDD(3) * mrSges(6,3) + qJD(3) * t317;
t302 = mrSges(6,1) * t323 - mrSges(6,3) * t324;
t367 = -mrSges(5,1) * t323 - mrSges(5,2) * t324 - t302;
t275 = m(5) * t286 - qJDD(3) * mrSges(5,2) - qJD(3) * t316 + t374 * t306 + t367 * t323 + t362;
t355 = -t345 * t289 + t344 * t290;
t285 = t324 * t376 - t355;
t307 = t331 * t344 + t332 * t345;
t315 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t323;
t282 = -qJDD(3) * pkin(4) - t350 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t301) * t324 + t355;
t318 = -mrSges(6,2) * t323 + qJD(3) * mrSges(6,3);
t358 = -m(6) * t282 + qJDD(3) * mrSges(6,1) + qJD(3) * t318;
t276 = m(5) * t285 + qJDD(3) * mrSges(5,1) + qJD(3) * t315 + t374 * t307 + t367 * t324 + t358;
t271 = t344 * t275 + t345 * t276;
t370 = (t377 * qJD(3)) + t371 * t323 - t372 * t324;
t369 = t371 * qJD(3) + t378 * t323 + t373 * t324;
t368 = t372 * qJD(3) - t373 * t323 + t379 * t324;
t359 = t345 * t275 - t276 * t344;
t330 = (mrSges(4,1) * t346 + mrSges(4,2) * t348) * qJD(1);
t333 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t366;
t335 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t365;
t356 = (m(4) * t308 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t332 + qJD(3) * t333 - t330 * t365 + t271) * t348 + (m(4) * t309 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t331 - qJD(3) * t335 - t330 * t366 + t359) * t346;
t292 = -t331 * pkin(3) + qJDD(4) + t334 * t365 + (-qJ(4) * t343 + t375) * t351 + t353;
t284 = -0.2e1 * qJD(5) * t324 + (qJD(3) * t323 - t307) * qJ(5) + (qJD(3) * t324 + t306) * pkin(4) + t292;
t278 = m(6) * t284 + t306 * mrSges(6,1) - mrSges(6,3) * t307 - t317 * t324 + t323 * t318;
t277 = m(5) * t292 + mrSges(5,1) * t306 + t307 * mrSges(5,2) + t315 * t323 + t324 * t316 + t278;
t327 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t348 - Ifges(4,4) * t346) * qJD(1);
t326 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t348 - Ifges(4,2) * t346) * qJD(1);
t322 = -qJDD(1) * pkin(1) + t352;
t320 = t351 * pkin(1) - t353;
t314 = t375 * t351 + t353;
t279 = t307 * mrSges(6,2) + t324 * t302 - t358;
t268 = mrSges(5,2) * t292 + mrSges(6,2) * t282 - mrSges(5,3) * t285 - mrSges(6,3) * t284 - qJ(5) * t278 - t369 * qJD(3) + t372 * qJDD(3) - t373 * t306 + t379 * t307 + t370 * t323;
t267 = -mrSges(5,1) * t292 - mrSges(6,1) * t284 + mrSges(6,2) * t281 + mrSges(5,3) * t286 - pkin(4) * t278 + t368 * qJD(3) + t371 * qJDD(3) + t378 * t306 + t373 * t307 + t370 * t324;
t266 = m(3) * t322 + qJDD(1) * mrSges(3,2) - (mrSges(3,3) * t351) + t356;
t1 = [mrSges(2,1) * t360 - mrSges(2,2) * t357 + mrSges(3,2) * t322 - mrSges(3,3) * t320 + t348 * (mrSges(4,2) * t314 - mrSges(4,3) * t308 + Ifges(4,1) * t332 + Ifges(4,4) * t331 + Ifges(4,5) * qJDD(3) - qJ(4) * t271 - qJD(3) * t326 - t267 * t344 + t268 * t345) - t346 * (-mrSges(4,1) * t314 + mrSges(4,3) * t309 + Ifges(4,4) * t332 + Ifges(4,2) * t331 + Ifges(4,6) * qJDD(3) - pkin(3) * t277 + qJ(4) * t359 + qJD(3) * t327 + t345 * t267 + t344 * t268) - pkin(6) * t356 - pkin(1) * t266 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-m(3) * t320 + m(4) * t314 - mrSges(4,1) * t331 + mrSges(3,2) * t351 + mrSges(4,2) * t332 + t277 + qJDD(1) * mrSges(3,3) + (t333 * t346 + t335 * t348) * qJD(1)) * qJ(2); t266; Ifges(4,5) * t332 + Ifges(4,6) * t331 + mrSges(4,1) * t308 - mrSges(4,2) * t309 + mrSges(5,1) * t285 - mrSges(5,2) * t286 - mrSges(6,1) * t282 + mrSges(6,3) * t281 - pkin(4) * t279 + qJ(5) * t362 + pkin(3) * t271 + t369 * t324 + (-qJ(5) * t302 + t368) * t323 + t372 * t307 + (-mrSges(6,2) * qJ(5) - t371) * t306 + (t326 * t348 + t327 * t346) * qJD(1) + (Ifges(4,3) - t377) * qJDD(3); t277; t279;];
tauJ = t1;
