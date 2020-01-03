% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4PRRR4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4PRRR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR4_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR4_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR4_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR4_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR4_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:34
% EndTime: 2019-12-31 16:32:35
% DurationCPUTime: 0.98s
% Computational Cost: add. (9459->185), mult. (18695->244), div. (0->0), fcn. (11357->8), ass. (0->78)
t334 = sin(pkin(7));
t335 = cos(pkin(7));
t322 = g(1) * t334 - g(2) * t335;
t323 = -g(1) * t335 - g(2) * t334;
t338 = sin(qJ(2));
t341 = cos(qJ(2));
t306 = t338 * t322 + t341 * t323;
t342 = qJD(2) ^ 2;
t304 = -pkin(2) * t342 + qJDD(2) * pkin(5) + t306;
t333 = -g(3) + qJDD(1);
t337 = sin(qJ(3));
t340 = cos(qJ(3));
t295 = -t337 * t304 + t340 * t333;
t352 = qJD(2) * qJD(3);
t350 = t340 * t352;
t320 = qJDD(2) * t337 + t350;
t289 = (-t320 + t350) * pkin(6) + (t337 * t340 * t342 + qJDD(3)) * pkin(3) + t295;
t296 = t340 * t304 + t337 * t333;
t321 = qJDD(2) * t340 - t337 * t352;
t354 = qJD(2) * t337;
t326 = qJD(3) * pkin(3) - pkin(6) * t354;
t332 = t340 ^ 2;
t290 = -pkin(3) * t332 * t342 + pkin(6) * t321 - qJD(3) * t326 + t296;
t336 = sin(qJ(4));
t339 = cos(qJ(4));
t287 = t289 * t339 - t290 * t336;
t312 = (-t336 * t337 + t339 * t340) * qJD(2);
t294 = qJD(4) * t312 + t320 * t339 + t321 * t336;
t313 = (t336 * t340 + t337 * t339) * qJD(2);
t302 = -mrSges(5,1) * t312 + mrSges(5,2) * t313;
t331 = qJD(3) + qJD(4);
t307 = -mrSges(5,2) * t331 + mrSges(5,3) * t312;
t330 = qJDD(3) + qJDD(4);
t285 = m(5) * t287 + mrSges(5,1) * t330 - mrSges(5,3) * t294 - t302 * t313 + t307 * t331;
t288 = t289 * t336 + t290 * t339;
t293 = -qJD(4) * t313 - t320 * t336 + t321 * t339;
t308 = mrSges(5,1) * t331 - mrSges(5,3) * t313;
t286 = m(5) * t288 - mrSges(5,2) * t330 + mrSges(5,3) * t293 + t302 * t312 - t308 * t331;
t277 = t339 * t285 + t336 * t286;
t319 = (-mrSges(4,1) * t340 + mrSges(4,2) * t337) * qJD(2);
t353 = qJD(2) * t340;
t325 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t353;
t275 = m(4) * t295 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t320 + qJD(3) * t325 - t319 * t354 + t277;
t324 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t354;
t346 = -t285 * t336 + t339 * t286;
t276 = m(4) * t296 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t321 - qJD(3) * t324 + t319 * t353 + t346;
t347 = -t275 * t337 + t340 * t276;
t270 = m(3) * t306 - mrSges(3,1) * t342 - qJDD(2) * mrSges(3,2) + t347;
t305 = t341 * t322 - t338 * t323;
t345 = -qJDD(2) * pkin(2) - t305;
t303 = -pkin(5) * t342 + t345;
t291 = t326 * t354 - t321 * pkin(3) + (-pkin(6) * t332 - pkin(5)) * t342 + t345;
t344 = m(5) * t291 - t293 * mrSges(5,1) + mrSges(5,2) * t294 - t312 * t307 + t308 * t313;
t343 = -m(4) * t303 + t321 * mrSges(4,1) - mrSges(4,2) * t320 - t324 * t354 + t325 * t353 - t344;
t281 = m(3) * t305 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t342 + t343;
t266 = t338 * t270 + t341 * t281;
t264 = m(2) * t322 + t266;
t348 = t341 * t270 - t281 * t338;
t265 = m(2) * t323 + t348;
t355 = t335 * t264 + t334 * t265;
t271 = t340 * t275 + t337 * t276;
t351 = m(3) * t333 + t271;
t349 = -t264 * t334 + t335 * t265;
t311 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t337 + Ifges(4,4) * t340) * qJD(2);
t310 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t337 + Ifges(4,2) * t340) * qJD(2);
t309 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t337 + Ifges(4,6) * t340) * qJD(2);
t299 = Ifges(5,1) * t313 + Ifges(5,4) * t312 + Ifges(5,5) * t331;
t298 = Ifges(5,4) * t313 + Ifges(5,2) * t312 + Ifges(5,6) * t331;
t297 = Ifges(5,5) * t313 + Ifges(5,6) * t312 + Ifges(5,3) * t331;
t279 = mrSges(5,2) * t291 - mrSges(5,3) * t287 + Ifges(5,1) * t294 + Ifges(5,4) * t293 + Ifges(5,5) * t330 + t297 * t312 - t298 * t331;
t278 = -mrSges(5,1) * t291 + mrSges(5,3) * t288 + Ifges(5,4) * t294 + Ifges(5,2) * t293 + Ifges(5,6) * t330 - t297 * t313 + t299 * t331;
t267 = mrSges(4,2) * t303 - mrSges(4,3) * t295 + Ifges(4,1) * t320 + Ifges(4,4) * t321 + Ifges(4,5) * qJDD(3) - pkin(6) * t277 - qJD(3) * t310 - t278 * t336 + t279 * t339 + t309 * t353;
t260 = -mrSges(4,1) * t303 + mrSges(4,3) * t296 + Ifges(4,4) * t320 + Ifges(4,2) * t321 + Ifges(4,6) * qJDD(3) - pkin(3) * t344 + pkin(6) * t346 + qJD(3) * t311 + t339 * t278 + t336 * t279 - t309 * t354;
t259 = Ifges(3,6) * qJDD(2) + t342 * Ifges(3,5) - mrSges(3,1) * t333 + mrSges(3,3) * t306 - Ifges(4,5) * t320 - Ifges(4,6) * t321 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t295 + mrSges(4,2) * t296 - Ifges(5,5) * t294 - Ifges(5,6) * t293 - Ifges(5,3) * t330 - t313 * t298 + t312 * t299 - mrSges(5,1) * t287 + mrSges(5,2) * t288 - pkin(3) * t277 - pkin(2) * t271 + (-t310 * t337 + t311 * t340) * qJD(2);
t258 = mrSges(3,2) * t333 - mrSges(3,3) * t305 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t342 - pkin(5) * t271 - t260 * t337 + t267 * t340;
t257 = mrSges(2,2) * t333 - mrSges(2,3) * t322 - pkin(4) * t266 + t258 * t341 - t259 * t338;
t256 = -mrSges(2,1) * t333 + mrSges(2,3) * t323 - pkin(1) * t351 + pkin(4) * t348 + t338 * t258 + t341 * t259;
t1 = [-m(1) * g(1) + t349; -m(1) * g(2) + t355; -m(1) * g(3) + m(2) * t333 + t351; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t355 - t334 * t256 + t335 * t257; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t349 + t335 * t256 + t334 * t257; -mrSges(1,1) * g(2) + mrSges(2,1) * t322 + mrSges(3,1) * t305 + mrSges(1,2) * g(1) - mrSges(2,2) * t323 - mrSges(3,2) * t306 + Ifges(3,3) * qJDD(2) + pkin(1) * t266 + pkin(2) * t343 + pkin(5) * t347 + t340 * t260 + t337 * t267;];
tauB = t1;
