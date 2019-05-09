% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPRP8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% m [3x7]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 08:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPRP8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP8_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP8_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP8_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP8_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP8_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 08:29:21
% EndTime: 2019-05-07 08:29:39
% DurationCPUTime: 7.17s
% Computational Cost: add. (122137->376), mult. (241272->443), div. (0->0), fcn. (159910->8), ass. (0->136)
t352 = sin(qJ(1));
t355 = cos(qJ(1));
t336 = g(1) * t352 - t355 * g(2);
t357 = qJD(1) ^ 2;
t314 = -qJDD(1) * pkin(1) - pkin(7) * t357 - t336;
t351 = sin(qJ(2));
t354 = cos(qJ(2));
t379 = qJD(1) * qJD(2);
t375 = t354 * t379;
t330 = qJDD(1) * t351 + t375;
t376 = t351 * t379;
t331 = t354 * qJDD(1) - t376;
t252 = (-t330 - t375) * pkin(8) + (-t331 + t376) * pkin(2) + t314;
t337 = -g(1) * t355 - g(2) * t352;
t315 = -pkin(1) * t357 + qJDD(1) * pkin(7) + t337;
t299 = -g(3) * t351 + t354 * t315;
t329 = (-pkin(2) * t354 - pkin(8) * t351) * qJD(1);
t356 = qJD(2) ^ 2;
t380 = qJD(1) * t354;
t258 = -pkin(2) * t356 + qJDD(2) * pkin(8) + t329 * t380 + t299;
t350 = sin(qJ(3));
t388 = cos(qJ(3));
t235 = t252 * t388 - t350 * t258;
t236 = t350 * t252 + t388 * t258;
t381 = qJD(1) * t351;
t326 = -qJD(2) * t388 + t350 * t381;
t327 = t350 * qJD(2) + t381 * t388;
t340 = qJD(3) - t380;
t267 = Ifges(5,5) * t327 + Ifges(5,6) * t340 + Ifges(5,3) * t326;
t270 = Ifges(4,4) * t327 - Ifges(4,2) * t326 + Ifges(4,6) * t340;
t272 = Ifges(4,1) * t327 - Ifges(4,4) * t326 + Ifges(4,5) * t340;
t284 = qJD(3) * t327 - qJDD(2) * t388 + t330 * t350;
t285 = -t326 * qJD(3) + t350 * qJDD(2) + t330 * t388;
t292 = mrSges(5,1) * t326 - mrSges(5,3) * t327;
t325 = qJDD(3) - t331;
t291 = pkin(3) * t326 - qJ(4) * t327;
t339 = t340 ^ 2;
t214 = -t325 * pkin(3) - t339 * qJ(4) + t327 * t291 + qJDD(4) - t235;
t385 = t326 * t340;
t204 = (-t285 - t385) * pkin(9) + (t326 * t327 - t325) * pkin(4) + t214;
t389 = 2 * qJD(4);
t211 = -pkin(3) * t339 + t325 * qJ(4) - t326 * t291 + t340 * t389 + t236;
t300 = -pkin(4) * t340 - pkin(9) * t327;
t324 = t326 ^ 2;
t207 = -pkin(4) * t324 + pkin(9) * t284 + t300 * t340 + t211;
t349 = sin(qJ(5));
t353 = cos(qJ(5));
t197 = t353 * t204 - t207 * t349;
t287 = t326 * t353 - t327 * t349;
t234 = qJD(5) * t287 + t284 * t349 + t285 * t353;
t288 = t326 * t349 + t327 * t353;
t249 = -mrSges(7,1) * t287 + mrSges(7,2) * t288;
t250 = -mrSges(6,1) * t287 + mrSges(6,2) * t288;
t338 = qJD(5) - t340;
t261 = -mrSges(6,2) * t338 + mrSges(6,3) * t287;
t321 = qJDD(5) - t325;
t190 = -0.2e1 * qJD(6) * t288 + (t287 * t338 - t234) * qJ(6) + (t287 * t288 + t321) * pkin(5) + t197;
t260 = -mrSges(7,2) * t338 + mrSges(7,3) * t287;
t378 = m(7) * t190 + t321 * mrSges(7,1) + t338 * t260;
t181 = m(6) * t197 + mrSges(6,1) * t321 + t261 * t338 + (-t249 - t250) * t288 + (-mrSges(6,3) - mrSges(7,3)) * t234 + t378;
t198 = t349 * t204 + t353 * t207;
t233 = -qJD(5) * t288 + t284 * t353 - t285 * t349;
t263 = mrSges(7,1) * t338 - mrSges(7,3) * t288;
t264 = mrSges(6,1) * t338 - mrSges(6,3) * t288;
t262 = pkin(5) * t338 - qJ(6) * t288;
t286 = t287 ^ 2;
t194 = -pkin(5) * t286 + qJ(6) * t233 + 0.2e1 * qJD(6) * t287 - t262 * t338 + t198;
t377 = m(7) * t194 + t233 * mrSges(7,3) + t287 * t249;
t183 = m(6) * t198 + mrSges(6,3) * t233 + t250 * t287 + (-t263 - t264) * t338 + (-mrSges(6,2) - mrSges(7,2)) * t321 + t377;
t177 = t181 * t353 + t183 * t349;
t271 = Ifges(5,1) * t327 + Ifges(5,4) * t340 + Ifges(5,5) * t326;
t187 = -mrSges(7,3) * t234 - t249 * t288 + t378;
t242 = Ifges(6,4) * t288 + Ifges(6,2) * t287 + Ifges(6,6) * t338;
t244 = Ifges(6,1) * t288 + Ifges(6,4) * t287 + Ifges(6,5) * t338;
t241 = Ifges(7,4) * t288 + Ifges(7,2) * t287 + Ifges(7,6) * t338;
t243 = Ifges(7,1) * t288 + Ifges(7,4) * t287 + Ifges(7,5) * t338;
t370 = mrSges(7,1) * t190 - mrSges(7,2) * t194 + Ifges(7,5) * t234 + Ifges(7,6) * t233 + Ifges(7,3) * t321 + t288 * t241 - t287 * t243;
t362 = mrSges(6,1) * t197 - mrSges(6,2) * t198 + Ifges(6,5) * t234 + Ifges(6,6) * t233 + Ifges(6,3) * t321 + pkin(5) * t187 + t288 * t242 - t287 * t244 + t370;
t359 = mrSges(5,1) * t214 - mrSges(5,3) * t211 - Ifges(5,4) * t285 - Ifges(5,2) * t325 - Ifges(5,6) * t284 + pkin(4) * t177 - t326 * t271 + t362;
t297 = -mrSges(5,2) * t326 + mrSges(5,3) * t340;
t365 = -m(5) * t214 + t325 * mrSges(5,1) + t340 * t297 - t177;
t178 = -t181 * t349 + t353 * t183;
t296 = -mrSges(5,1) * t340 + mrSges(5,2) * t327;
t371 = m(5) * t211 + t325 * mrSges(5,3) + t340 * t296 + t178;
t391 = -(t267 - t270) * t327 + mrSges(4,1) * t235 - mrSges(4,2) * t236 + Ifges(4,5) * t285 - Ifges(4,6) * t284 + Ifges(4,3) * t325 + pkin(3) * (-mrSges(5,2) * t285 - t292 * t327 + t365) + qJ(4) * (-mrSges(5,2) * t284 - t292 * t326 + t371) + t326 * t272 - t359;
t298 = -t354 * g(3) - t351 * t315;
t257 = -qJDD(2) * pkin(2) - t356 * pkin(8) + t329 * t381 - t298;
t367 = t284 * pkin(3) + t257 + (-t285 + t385) * qJ(4);
t387 = pkin(3) * t340;
t208 = -pkin(4) * t284 - pkin(9) * t324 - t367 + (t300 - t387 + t389) * t327;
t201 = -pkin(5) * t233 - qJ(6) * t286 + t262 * t288 + qJDD(6) + t208;
t373 = m(7) * t201 - t233 * mrSges(7,1) + t234 * mrSges(7,2) - t287 * t260 + t288 * t263;
t185 = -m(6) * t208 + t233 * mrSges(6,1) - t234 * mrSges(6,2) + t287 * t261 - t288 * t264 - t373;
t213 = (-(2 * qJD(4)) + t387) * t327 + t367;
t184 = m(5) * t213 + t284 * mrSges(5,1) - t285 * mrSges(5,3) - t327 * t296 + t326 * t297 + t185;
t239 = Ifges(7,5) * t288 + Ifges(7,6) * t287 + Ifges(7,3) * t338;
t240 = Ifges(6,5) * t288 + Ifges(6,6) * t287 + Ifges(6,3) * t338;
t369 = -mrSges(7,1) * t201 + mrSges(7,3) * t194 + Ifges(7,4) * t234 + Ifges(7,2) * t233 + Ifges(7,6) * t321 + t338 * t243;
t169 = Ifges(6,4) * t234 + Ifges(6,2) * t233 + Ifges(6,6) * t321 + t338 * t244 - mrSges(6,1) * t208 + mrSges(6,3) * t198 - pkin(5) * t373 + qJ(6) * (-mrSges(7,2) * t321 - t263 * t338 + t377) + (-t240 - t239) * t288 + t369;
t368 = mrSges(7,2) * t201 - mrSges(7,3) * t190 + Ifges(7,1) * t234 + Ifges(7,4) * t233 + Ifges(7,5) * t321 + t287 * t239;
t176 = mrSges(6,2) * t208 - mrSges(6,3) * t197 + Ifges(6,1) * t234 + Ifges(6,4) * t233 + Ifges(6,5) * t321 - qJ(6) * t187 + t240 * t287 + (-t241 - t242) * t338 + t368;
t363 = -mrSges(5,1) * t213 + mrSges(5,2) * t211 - pkin(4) * t185 - pkin(9) * t178 - t353 * t169 - t349 * t176;
t269 = Ifges(5,4) * t327 + Ifges(5,2) * t340 + Ifges(5,6) * t326;
t383 = -Ifges(4,5) * t327 + Ifges(4,6) * t326 - Ifges(4,3) * t340 - t269;
t158 = -mrSges(4,1) * t257 + mrSges(4,3) * t236 - pkin(3) * t184 + (t272 + t271) * t340 + t383 * t327 + (Ifges(4,6) - Ifges(5,6)) * t325 + (Ifges(4,4) - Ifges(5,5)) * t285 + (-Ifges(4,2) - Ifges(5,3)) * t284 + t363;
t364 = mrSges(5,2) * t214 - mrSges(5,3) * t213 + Ifges(5,1) * t285 + Ifges(5,4) * t325 + Ifges(5,5) * t284 - pkin(9) * t177 - t169 * t349 + t353 * t176 + t340 * t267;
t159 = mrSges(4,2) * t257 - mrSges(4,3) * t235 + Ifges(4,1) * t285 - Ifges(4,4) * t284 + Ifges(4,5) * t325 - qJ(4) * t184 - t270 * t340 + t326 * t383 + t364;
t295 = mrSges(4,1) * t340 - mrSges(4,3) * t327;
t382 = -mrSges(4,1) * t326 - mrSges(4,2) * t327 - t292;
t386 = -mrSges(4,3) - mrSges(5,2);
t171 = m(4) * t236 - mrSges(4,2) * t325 + t284 * t386 - t295 * t340 + t326 * t382 + t371;
t294 = -mrSges(4,2) * t340 - mrSges(4,3) * t326;
t172 = m(4) * t235 + mrSges(4,1) * t325 + t285 * t386 + t294 * t340 + t327 * t382 + t365;
t168 = t388 * t171 - t172 * t350;
t180 = -m(4) * t257 - t284 * mrSges(4,1) - t285 * mrSges(4,2) - t326 * t294 - t327 * t295 - t184;
t312 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t351 + Ifges(3,2) * t354) * qJD(1);
t313 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t351 + Ifges(3,4) * t354) * qJD(1);
t390 = mrSges(3,1) * t298 - mrSges(3,2) * t299 + Ifges(3,5) * t330 + Ifges(3,6) * t331 + Ifges(3,3) * qJDD(2) + pkin(2) * t180 + pkin(8) * t168 + (t351 * t312 - t354 * t313) * qJD(1) + t158 * t388 + t350 * t159;
t328 = (-mrSges(3,1) * t354 + mrSges(3,2) * t351) * qJD(1);
t333 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t381;
t166 = m(3) * t299 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t331 - qJD(2) * t333 + t328 * t380 + t168;
t334 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t380;
t179 = m(3) * t298 + qJDD(2) * mrSges(3,1) - t330 * mrSges(3,3) + qJD(2) * t334 - t328 * t381 + t180;
t374 = t354 * t166 - t179 * t351;
t167 = t350 * t171 + t172 * t388;
t311 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t351 + Ifges(3,6) * t354) * qJD(1);
t155 = mrSges(3,2) * t314 - mrSges(3,3) * t298 + Ifges(3,1) * t330 + Ifges(3,4) * t331 + Ifges(3,5) * qJDD(2) - pkin(8) * t167 - qJD(2) * t312 - t350 * t158 + t159 * t388 + t311 * t380;
t157 = -mrSges(3,1) * t314 + mrSges(3,3) * t299 + Ifges(3,4) * t330 + Ifges(3,2) * t331 + Ifges(3,6) * qJDD(2) - pkin(2) * t167 + qJD(2) * t313 - t311 * t381 - t391;
t361 = -m(3) * t314 + t331 * mrSges(3,1) - t330 * mrSges(3,2) - t333 * t381 + t334 * t380 - t167;
t366 = mrSges(2,1) * t336 - mrSges(2,2) * t337 + Ifges(2,3) * qJDD(1) + pkin(1) * t361 + pkin(7) * t374 + t351 * t155 + t354 * t157;
t163 = m(2) * t336 + qJDD(1) * mrSges(2,1) - t357 * mrSges(2,2) + t361;
t162 = t166 * t351 + t179 * t354;
t160 = m(2) * t337 - mrSges(2,1) * t357 - qJDD(1) * mrSges(2,2) + t374;
t153 = mrSges(2,1) * g(3) + mrSges(2,3) * t337 + t357 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t162 - t390;
t152 = -mrSges(2,2) * g(3) - mrSges(2,3) * t336 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t357 - pkin(7) * t162 + t155 * t354 - t157 * t351;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t355 * t152 - t352 * t153 - pkin(6) * (t160 * t352 + t163 * t355), t152, t155, t159, -t269 * t326 + t364, t176, -t241 * t338 + t368; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t352 * t152 + t355 * t153 + pkin(6) * (t160 * t355 - t163 * t352), t153, t157, t158, -t327 * t267 - t359, t169, -t288 * t239 + t369; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t366, t366, t390, t391, Ifges(5,5) * t285 + Ifges(5,6) * t325 + Ifges(5,3) * t284 + t327 * t269 - t340 * t271 - t363, t362, t370;];
m_new  = t1;
