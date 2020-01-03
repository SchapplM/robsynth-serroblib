% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4PRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4PRRP3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP3_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP3_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP3_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP3_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:47
% EndTime: 2019-12-31 16:26:48
% DurationCPUTime: 0.76s
% Computational Cost: add. (3830->166), mult. (7431->204), div. (0->0), fcn. (3674->6), ass. (0->73)
t368 = Ifges(4,1) + Ifges(5,1);
t361 = Ifges(4,4) + Ifges(5,4);
t360 = Ifges(4,5) + Ifges(5,5);
t367 = Ifges(4,2) + Ifges(5,2);
t366 = Ifges(4,6) + Ifges(5,6);
t365 = Ifges(4,3) + Ifges(5,3);
t339 = qJD(2) ^ 2;
t333 = sin(pkin(6));
t334 = cos(pkin(6));
t320 = t333 * g(1) - t334 * g(2);
t321 = -t334 * g(1) - t333 * g(2);
t336 = sin(qJ(2));
t338 = cos(qJ(2));
t298 = t338 * t320 - t336 * t321;
t341 = -qJDD(2) * pkin(2) - t298;
t296 = -t339 * pkin(5) + t341;
t335 = sin(qJ(3));
t337 = cos(qJ(3));
t351 = qJD(2) * qJD(3);
t346 = t337 * t351;
t317 = t335 * qJDD(2) + t346;
t318 = t337 * qJDD(2) - t335 * t351;
t352 = qJD(2) * t337;
t326 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t352;
t353 = qJD(2) * t335;
t322 = qJD(3) * pkin(3) - qJ(4) * t353;
t331 = t337 ^ 2;
t292 = t322 * t353 - t318 * pkin(3) + qJDD(4) + (-qJ(4) * t331 - pkin(5)) * t339 + t341;
t325 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t352;
t342 = m(5) * t292 - t318 * mrSges(5,1) - t325 * t352;
t362 = -mrSges(4,2) - mrSges(5,2);
t364 = -m(4) * t296 + t318 * mrSges(4,1) + t362 * t317 + t326 * t352 - t342;
t363 = pkin(3) * t339;
t299 = t336 * t320 + t338 * t321;
t297 = -t339 * pkin(2) + qJDD(2) * pkin(5) + t299;
t332 = -g(3) + qJDD(1);
t294 = t337 * t297 + t335 * t332;
t316 = (-mrSges(4,1) * t337 + mrSges(4,2) * t335) * qJD(2);
t350 = qJD(2) * qJD(4);
t291 = t318 * qJ(4) - qJD(3) * t322 - t331 * t363 + 0.2e1 * t337 * t350 + t294;
t315 = (-mrSges(5,1) * t337 + mrSges(5,2) * t335) * qJD(2);
t348 = m(5) * t291 + t318 * mrSges(5,3) + t315 * t352;
t323 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t353;
t354 = -qJD(3) * mrSges(4,1) + mrSges(4,3) * t353 - t323;
t286 = m(4) * t294 + t318 * mrSges(4,3) + t354 * qJD(3) + t362 * qJDD(3) + t316 * t352 + t348;
t284 = t337 * t286;
t328 = t337 * t332;
t293 = -t335 * t297 + t328;
t290 = qJDD(3) * pkin(3) + t328 + (-t317 + t346) * qJ(4) + (t337 * t363 - t297 - 0.2e1 * t350) * t335;
t349 = m(5) * t290 + qJDD(3) * mrSges(5,1) + qJD(3) * t325;
t285 = m(4) * t293 + qJDD(3) * mrSges(4,1) + qJD(3) * t326 + (-mrSges(4,3) - mrSges(5,3)) * t317 + (-t315 - t316) * t353 + t349;
t278 = m(3) * t299 - t339 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t335 * t285 + t284;
t343 = qJD(2) * t354;
t281 = m(3) * t298 + qJDD(2) * mrSges(3,1) - t339 * mrSges(3,2) + t335 * t343 + t364;
t273 = t336 * t278 + t338 * t281;
t271 = m(2) * t320 + t273;
t344 = t338 * t278 - t336 * t281;
t272 = m(2) * t321 + t344;
t358 = t334 * t271 + t333 * t272;
t279 = t337 * t285 + t335 * t286;
t357 = t365 * qJD(3) + (t360 * t335 + t366 * t337) * qJD(2);
t356 = -t366 * qJD(3) + (-t361 * t335 - t367 * t337) * qJD(2);
t355 = t360 * qJD(3) + (t368 * t335 + t361 * t337) * qJD(2);
t347 = m(3) * t332 + t279;
t345 = -t333 * t271 + t334 * t272;
t287 = -t317 * mrSges(5,3) - t315 * t353 + t349;
t275 = mrSges(4,2) * t296 + mrSges(5,2) * t292 - mrSges(4,3) * t293 - mrSges(5,3) * t290 - qJ(4) * t287 + t356 * qJD(3) + t360 * qJDD(3) + t368 * t317 + t361 * t318 + t357 * t352;
t274 = -mrSges(4,1) * t296 + mrSges(4,3) * t294 - mrSges(5,1) * t292 + mrSges(5,3) * t291 - pkin(3) * t342 + qJ(4) * t348 + t367 * t318 + (-pkin(3) * mrSges(5,2) + t361) * t317 + (-qJ(4) * mrSges(5,2) + t366) * qJDD(3) + (-qJ(4) * t323 + t355) * qJD(3) + (-pkin(3) * t323 - t357) * t353;
t267 = -mrSges(3,1) * t332 - mrSges(4,1) * t293 - mrSges(5,1) * t290 + mrSges(4,2) * t294 + mrSges(5,2) * t291 + mrSges(3,3) * t299 + t339 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t279 - pkin(3) * t287 - t366 * t318 - t360 * t317 - t365 * qJDD(3) + (t356 * t335 + t355 * t337) * qJD(2);
t266 = mrSges(3,2) * t332 - mrSges(3,3) * t298 + Ifges(3,5) * qJDD(2) - t339 * Ifges(3,6) - pkin(5) * t279 - t335 * t274 + t337 * t275;
t265 = mrSges(2,2) * t332 - mrSges(2,3) * t320 - pkin(4) * t273 + t338 * t266 - t336 * t267;
t264 = -mrSges(2,1) * t332 + mrSges(2,3) * t321 - pkin(1) * t347 + pkin(4) * t344 + t336 * t266 + t338 * t267;
t1 = [-m(1) * g(1) + t345; -m(1) * g(2) + t358; -m(1) * g(3) + m(2) * t332 + t347; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t358 - t333 * t264 + t334 * t265; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t345 + t334 * t264 + t333 * t265; pkin(1) * t273 + mrSges(2,1) * t320 - mrSges(2,2) * t321 + t337 * t274 + pkin(2) * t364 + pkin(5) * t284 + mrSges(3,1) * t298 - mrSges(3,2) * t299 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(3,3) * qJDD(2) + (pkin(2) * t343 - pkin(5) * t285 + t275) * t335;];
tauB = t1;
