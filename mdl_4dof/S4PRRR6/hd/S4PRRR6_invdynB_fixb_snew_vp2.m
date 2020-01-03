% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4PRRR6
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
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4PRRR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR6_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR6_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR6_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR6_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR6_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:34:45
% EndTime: 2019-12-31 16:34:46
% DurationCPUTime: 0.92s
% Computational Cost: add. (8719->186), mult. (17074->244), div. (0->0), fcn. (10311->8), ass. (0->77)
t352 = cos(pkin(7));
t332 = sin(pkin(7));
t322 = -t352 * g(1) - t332 * g(2);
t331 = -g(3) + qJDD(1);
t335 = sin(qJ(2));
t338 = cos(qJ(2));
t306 = t338 * t322 + t335 * t331;
t339 = qJD(2) ^ 2;
t302 = -t339 * pkin(2) + qJDD(2) * pkin(5) + t306;
t321 = t332 * g(1) - t352 * g(2);
t334 = sin(qJ(3));
t337 = cos(qJ(3));
t293 = -t334 * t302 - t337 * t321;
t348 = qJD(2) * qJD(3);
t347 = t337 * t348;
t319 = t334 * qJDD(2) + t347;
t287 = (-t319 + t347) * pkin(6) + (t334 * t337 * t339 + qJDD(3)) * pkin(3) + t293;
t294 = t337 * t302 - t334 * t321;
t320 = t337 * qJDD(2) - t334 * t348;
t350 = qJD(2) * t334;
t325 = qJD(3) * pkin(3) - pkin(6) * t350;
t330 = t337 ^ 2;
t288 = -t330 * t339 * pkin(3) + t320 * pkin(6) - qJD(3) * t325 + t294;
t333 = sin(qJ(4));
t336 = cos(qJ(4));
t285 = t336 * t287 - t333 * t288;
t310 = (-t333 * t334 + t336 * t337) * qJD(2);
t292 = t310 * qJD(4) + t336 * t319 + t333 * t320;
t311 = (t333 * t337 + t334 * t336) * qJD(2);
t299 = -t310 * mrSges(5,1) + t311 * mrSges(5,2);
t329 = qJD(3) + qJD(4);
t303 = -t329 * mrSges(5,2) + t310 * mrSges(5,3);
t328 = qJDD(3) + qJDD(4);
t283 = m(5) * t285 + t328 * mrSges(5,1) - t292 * mrSges(5,3) - t311 * t299 + t329 * t303;
t286 = t333 * t287 + t336 * t288;
t291 = -t311 * qJD(4) - t333 * t319 + t336 * t320;
t304 = t329 * mrSges(5,1) - t311 * mrSges(5,3);
t284 = m(5) * t286 - t328 * mrSges(5,2) + t291 * mrSges(5,3) + t310 * t299 - t329 * t304;
t275 = t336 * t283 + t333 * t284;
t318 = (-mrSges(4,1) * t337 + mrSges(4,2) * t334) * qJD(2);
t349 = qJD(2) * t337;
t324 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t349;
t273 = m(4) * t293 + qJDD(3) * mrSges(4,1) - t319 * mrSges(4,3) + qJD(3) * t324 - t318 * t350 + t275;
t323 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t350;
t343 = -t333 * t283 + t336 * t284;
t274 = m(4) * t294 - qJDD(3) * mrSges(4,2) + t320 * mrSges(4,3) - qJD(3) * t323 + t318 * t349 + t343;
t344 = -t334 * t273 + t337 * t274;
t268 = m(3) * t306 - t339 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t344;
t305 = -t335 * t322 + t338 * t331;
t342 = -qJDD(2) * pkin(2) - t305;
t301 = -t339 * pkin(5) + t342;
t289 = t325 * t350 - t320 * pkin(3) + (-pkin(6) * t330 - pkin(5)) * t339 + t342;
t341 = m(5) * t289 - t291 * mrSges(5,1) + t292 * mrSges(5,2) - t310 * t303 + t311 * t304;
t340 = -m(4) * t301 + t320 * mrSges(4,1) - t319 * mrSges(4,2) - t323 * t350 + t324 * t349 - t341;
t279 = m(3) * t305 + qJDD(2) * mrSges(3,1) - t339 * mrSges(3,2) + t340;
t345 = t338 * t268 - t335 * t279;
t263 = m(2) * t322 + t345;
t271 = t337 * t273 + t334 * t274;
t270 = (m(2) + m(3)) * t321 - t271;
t351 = t332 * t263 + t352 * t270;
t264 = t335 * t268 + t338 * t279;
t346 = t352 * t263 - t332 * t270;
t309 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t334 + Ifges(4,4) * t337) * qJD(2);
t308 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t334 + Ifges(4,2) * t337) * qJD(2);
t307 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t334 + Ifges(4,6) * t337) * qJD(2);
t297 = Ifges(5,1) * t311 + Ifges(5,4) * t310 + Ifges(5,5) * t329;
t296 = Ifges(5,4) * t311 + Ifges(5,2) * t310 + Ifges(5,6) * t329;
t295 = Ifges(5,5) * t311 + Ifges(5,6) * t310 + Ifges(5,3) * t329;
t277 = mrSges(5,2) * t289 - mrSges(5,3) * t285 + Ifges(5,1) * t292 + Ifges(5,4) * t291 + Ifges(5,5) * t328 + t310 * t295 - t329 * t296;
t276 = -mrSges(5,1) * t289 + mrSges(5,3) * t286 + Ifges(5,4) * t292 + Ifges(5,2) * t291 + Ifges(5,6) * t328 - t311 * t295 + t329 * t297;
t265 = mrSges(4,2) * t301 - mrSges(4,3) * t293 + Ifges(4,1) * t319 + Ifges(4,4) * t320 + Ifges(4,5) * qJDD(3) - pkin(6) * t275 - qJD(3) * t308 - t333 * t276 + t336 * t277 + t307 * t349;
t260 = -mrSges(4,1) * t301 + mrSges(4,3) * t294 + Ifges(4,4) * t319 + Ifges(4,2) * t320 + Ifges(4,6) * qJDD(3) - pkin(3) * t341 + pkin(6) * t343 + qJD(3) * t309 + t336 * t276 + t333 * t277 - t307 * t350;
t259 = Ifges(3,6) * qJDD(2) + t339 * Ifges(3,5) + mrSges(3,1) * t321 + mrSges(3,3) * t306 - Ifges(4,5) * t319 - Ifges(4,6) * t320 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t293 + mrSges(4,2) * t294 - Ifges(5,5) * t292 - Ifges(5,6) * t291 - Ifges(5,3) * t328 - t311 * t296 + t310 * t297 - mrSges(5,1) * t285 + mrSges(5,2) * t286 - pkin(3) * t275 - pkin(2) * t271 + (-t334 * t308 + t337 * t309) * qJD(2);
t258 = -mrSges(3,2) * t321 - mrSges(3,3) * t305 + Ifges(3,5) * qJDD(2) - t339 * Ifges(3,6) - pkin(5) * t271 - t334 * t260 + t337 * t265;
t257 = -mrSges(2,1) * t331 - mrSges(3,1) * t305 + mrSges(3,2) * t306 + mrSges(2,3) * t322 - Ifges(3,3) * qJDD(2) - pkin(1) * t264 - pkin(2) * t340 - pkin(5) * t344 - t337 * t260 - t334 * t265;
t256 = mrSges(2,2) * t331 - mrSges(2,3) * t321 - pkin(4) * t264 + t338 * t258 - t335 * t259;
t1 = [-m(1) * g(1) + t346; -m(1) * g(2) + t351; -m(1) * g(3) + m(2) * t331 + t264; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t351 + t352 * t256 - t332 * t257; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t346 + t332 * t256 + t352 * t257; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t321 - mrSges(2,2) * t322 + t335 * t258 + t338 * t259 + pkin(1) * (m(3) * t321 - t271) + pkin(4) * t345;];
tauB = t1;
