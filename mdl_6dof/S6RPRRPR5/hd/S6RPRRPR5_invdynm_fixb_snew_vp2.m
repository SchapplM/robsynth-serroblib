% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-05-05 22:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRPR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:38:34
% EndTime: 2019-05-05 22:39:02
% DurationCPUTime: 16.29s
% Computational Cost: add. (250586->365), mult. (611617->436), div. (0->0), fcn. (473082->10), ass. (0->147)
t334 = sin(qJ(1));
t337 = cos(qJ(1));
t312 = -g(1) * t337 - g(2) * t334;
t338 = qJD(1) ^ 2;
t305 = -pkin(1) * t338 + qJDD(1) * qJ(2) + t312;
t329 = sin(pkin(10));
t330 = cos(pkin(10));
t368 = qJD(1) * qJD(2);
t366 = -t330 * g(3) - 0.2e1 * t329 * t368;
t377 = pkin(2) * t330;
t273 = (-pkin(7) * qJDD(1) + t338 * t377 - t305) * t329 + t366;
t295 = -g(3) * t329 + (t305 + 0.2e1 * t368) * t330;
t367 = qJDD(1) * t330;
t325 = t330 ^ 2;
t373 = t325 * t338;
t275 = -pkin(2) * t373 + pkin(7) * t367 + t295;
t333 = sin(qJ(3));
t336 = cos(qJ(3));
t244 = t273 * t336 - t333 * t275;
t356 = t329 * t336 + t330 * t333;
t355 = -t329 * t333 + t330 * t336;
t303 = t355 * qJD(1);
t369 = t303 * qJD(3);
t293 = qJDD(1) * t356 + t369;
t304 = t356 * qJD(1);
t210 = (-t293 + t369) * pkin(8) + (t303 * t304 + qJDD(3)) * pkin(3) + t244;
t245 = t273 * t333 + t275 * t336;
t292 = -t304 * qJD(3) + qJDD(1) * t355;
t298 = qJD(3) * pkin(3) - pkin(8) * t304;
t302 = t303 ^ 2;
t218 = -pkin(3) * t302 + pkin(8) * t292 - qJD(3) * t298 + t245;
t332 = sin(qJ(4));
t378 = cos(qJ(4));
t206 = t210 * t378 - t218 * t332;
t207 = t210 * t332 + t218 * t378;
t283 = t303 * t332 + t304 * t378;
t240 = qJD(4) * t283 - t292 * t378 + t293 * t332;
t282 = -t303 * t378 + t304 * t332;
t241 = -qJD(4) * t282 + t292 * t332 + t293 * t378;
t326 = qJD(3) + qJD(4);
t251 = Ifges(5,4) * t283 - Ifges(5,2) * t282 + Ifges(5,6) * t326;
t260 = -mrSges(6,2) * t282 - mrSges(6,3) * t283;
t269 = mrSges(6,1) * t282 - mrSges(6,3) * t326;
t323 = qJDD(3) + qJDD(4);
t258 = pkin(4) * t282 - qJ(5) * t283;
t322 = t326 ^ 2;
t202 = -t323 * pkin(4) - t322 * qJ(5) + t258 * t283 + qJDD(5) - t206;
t374 = t282 * t326;
t196 = (t282 * t283 - t323) * pkin(9) + (t241 + t374) * pkin(5) + t202;
t271 = pkin(5) * t283 - pkin(9) * t326;
t278 = t282 ^ 2;
t324 = t329 ^ 2;
t311 = t334 * g(1) - g(2) * t337;
t361 = qJDD(2) - t311;
t291 = (-pkin(1) - t377) * qJDD(1) + (-qJ(2) + (-t324 - t325) * pkin(7)) * t338 + t361;
t226 = -t292 * pkin(3) - t302 * pkin(8) + t298 * t304 + t291;
t379 = -2 * qJD(5);
t341 = (-t241 + t374) * qJ(5) + t226 + (pkin(4) * t326 + t379) * t283;
t197 = t341 + (pkin(4) + pkin(9)) * t240 - t283 * t271 - t278 * pkin(5);
t331 = sin(qJ(6));
t335 = cos(qJ(6));
t193 = t196 * t335 - t197 * t331;
t263 = t282 * t335 - t326 * t331;
t215 = qJD(6) * t263 + t240 * t331 + t323 * t335;
t239 = qJDD(6) + t241;
t264 = t282 * t331 + t326 * t335;
t242 = -mrSges(7,1) * t263 + mrSges(7,2) * t264;
t277 = qJD(6) + t283;
t246 = -mrSges(7,2) * t277 + mrSges(7,3) * t263;
t190 = m(7) * t193 + mrSges(7,1) * t239 - mrSges(7,3) * t215 - t242 * t264 + t246 * t277;
t194 = t196 * t331 + t197 * t335;
t214 = -qJD(6) * t264 + t240 * t335 - t323 * t331;
t247 = mrSges(7,1) * t277 - mrSges(7,3) * t264;
t191 = m(7) * t194 - mrSges(7,2) * t239 + mrSges(7,3) * t214 + t242 * t263 - t247 * t277;
t178 = t335 * t190 + t331 * t191;
t351 = -t322 * pkin(4) + t323 * qJ(5) - t282 * t258 + t207;
t199 = -t240 * pkin(5) - t278 * pkin(9) + ((2 * qJD(5)) + t271) * t326 + t351;
t219 = Ifges(7,5) * t264 + Ifges(7,6) * t263 + Ifges(7,3) * t277;
t221 = Ifges(7,1) * t264 + Ifges(7,4) * t263 + Ifges(7,5) * t277;
t181 = -mrSges(7,1) * t199 + mrSges(7,3) * t194 + Ifges(7,4) * t215 + Ifges(7,2) * t214 + Ifges(7,6) * t239 - t219 * t264 + t221 * t277;
t220 = Ifges(7,4) * t264 + Ifges(7,2) * t263 + Ifges(7,6) * t277;
t182 = mrSges(7,2) * t199 - mrSges(7,3) * t193 + Ifges(7,1) * t215 + Ifges(7,4) * t214 + Ifges(7,5) * t239 + t219 * t263 - t220 * t277;
t200 = t326 * t379 - t351;
t248 = Ifges(6,5) * t326 - Ifges(6,6) * t283 + Ifges(6,3) * t282;
t347 = -mrSges(6,2) * t202 + mrSges(6,3) * t200 - Ifges(6,1) * t323 + Ifges(6,4) * t241 - Ifges(6,5) * t240 + pkin(9) * t178 + t331 * t181 - t182 * t335 + t248 * t283;
t195 = -m(7) * t199 + t214 * mrSges(7,1) - mrSges(7,2) * t215 + t263 * t246 - t247 * t264;
t270 = mrSges(6,1) * t283 + mrSges(6,2) * t326;
t348 = -m(6) * t200 + mrSges(6,3) * t323 + t270 * t326 - t195;
t352 = -m(6) * t202 - mrSges(6,1) * t241 - t260 * t283 - t178;
t250 = Ifges(6,4) * t326 - Ifges(6,2) * t283 + Ifges(6,6) * t282;
t371 = Ifges(5,1) * t283 - Ifges(5,4) * t282 + Ifges(5,5) * t326 - t250;
t383 = -mrSges(5,2) * t207 + pkin(4) * (-t323 * mrSges(6,2) - t326 * t269 + t352) + qJ(5) * (-t240 * mrSges(6,1) - t282 * t260 + t348) + mrSges(5,1) * t206 + t283 * t251 - Ifges(5,6) * t240 + Ifges(5,5) * t241 + Ifges(5,3) * t323 - t347 + t371 * t282;
t259 = mrSges(5,1) * t282 + mrSges(5,2) * t283;
t267 = -mrSges(5,2) * t326 - mrSges(5,3) * t282;
t174 = m(5) * t206 - t241 * mrSges(5,3) - t283 * t259 + (t267 - t269) * t326 + (mrSges(5,1) - mrSges(6,2)) * t323 + t352;
t268 = mrSges(5,1) * t326 - mrSges(5,3) * t283;
t185 = m(5) * t207 - t323 * mrSges(5,2) - t326 * t268 + (-t259 - t260) * t282 + (-mrSges(5,3) - mrSges(6,1)) * t240 + t348;
t170 = t174 * t378 + t185 * t332;
t280 = Ifges(4,4) * t304 + Ifges(4,2) * t303 + Ifges(4,6) * qJD(3);
t281 = Ifges(4,1) * t304 + Ifges(4,4) * t303 + Ifges(4,5) * qJD(3);
t382 = mrSges(4,1) * t244 - mrSges(4,2) * t245 + Ifges(4,5) * t293 + Ifges(4,6) * t292 + Ifges(4,3) * qJDD(3) + pkin(3) * t170 + t280 * t304 - t303 * t281 + t383;
t287 = -mrSges(4,1) * t303 + mrSges(4,2) * t304;
t296 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t303;
t167 = m(4) * t244 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t293 + qJD(3) * t296 - t287 * t304 + t170;
t297 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t304;
t362 = -t174 * t332 + t185 * t378;
t168 = m(4) * t245 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t292 - qJD(3) * t297 + t287 * t303 + t362;
t162 = t167 * t336 + t168 * t333;
t294 = -t329 * t305 + t366;
t359 = Ifges(3,4) * t329 + Ifges(3,2) * t330;
t360 = Ifges(3,1) * t329 + Ifges(3,4) * t330;
t381 = -mrSges(3,1) * t294 + mrSges(3,2) * t295 - pkin(2) * t162 - (t329 * t359 - t330 * t360) * t338 - t382;
t376 = Ifges(5,4) + Ifges(6,6);
t375 = mrSges(3,2) * t329;
t179 = -t190 * t331 + t191 * t335;
t252 = Ifges(6,1) * t326 - Ifges(6,4) * t283 + Ifges(6,5) * t282;
t372 = -Ifges(5,5) * t283 + Ifges(5,6) * t282 - Ifges(5,3) * t326 - t252;
t358 = Ifges(3,5) * t329 + Ifges(3,6) * t330;
t370 = t338 * t358;
t354 = mrSges(3,3) * qJDD(1) + t338 * (-mrSges(3,1) * t330 + t375);
t160 = m(3) * t294 - t329 * t354 + t162;
t363 = -t333 * t167 + t168 * t336;
t161 = m(3) * t295 + t330 * t354 + t363;
t364 = -t160 * t329 + t161 * t330;
t204 = t240 * pkin(4) + t341;
t175 = m(6) * t204 - mrSges(6,2) * t240 - mrSges(6,3) * t241 - t269 * t282 - t270 * t283 + t179;
t346 = -mrSges(6,1) * t200 + mrSges(6,2) * t204 - pkin(5) * t195 - pkin(9) * t179 - t335 * t181 - t331 * t182;
t158 = -mrSges(5,1) * t226 + mrSges(5,3) * t207 - pkin(4) * t175 + t371 * t326 + (Ifges(5,6) - Ifges(6,5)) * t323 + t372 * t283 + t376 * t241 + (-Ifges(5,2) - Ifges(6,3)) * t240 + t346;
t350 = mrSges(7,1) * t193 - mrSges(7,2) * t194 + Ifges(7,5) * t215 + Ifges(7,6) * t214 + Ifges(7,3) * t239 + t220 * t264 - t263 * t221;
t345 = mrSges(6,1) * t202 - mrSges(6,3) * t204 + pkin(5) * t178 + t350;
t163 = mrSges(5,2) * t226 - mrSges(5,3) * t206 - qJ(5) * t175 + t345 + (-t251 + t248) * t326 + (Ifges(5,5) - Ifges(6,4)) * t323 + t372 * t282 + (Ifges(5,1) + Ifges(6,2)) * t241 - t376 * t240;
t279 = Ifges(4,5) * t304 + Ifges(4,6) * t303 + Ifges(4,3) * qJD(3);
t349 = m(5) * t226 + mrSges(5,1) * t240 + mrSges(5,2) * t241 + t267 * t282 + t268 * t283 + t175;
t153 = -mrSges(4,1) * t291 + mrSges(4,3) * t245 + Ifges(4,4) * t293 + Ifges(4,2) * t292 + Ifges(4,6) * qJDD(3) - pkin(3) * t349 + pkin(8) * t362 + qJD(3) * t281 + t158 * t378 + t332 * t163 - t304 * t279;
t154 = mrSges(4,2) * t291 - mrSges(4,3) * t244 + Ifges(4,1) * t293 + Ifges(4,4) * t292 + Ifges(4,5) * qJDD(3) - pkin(8) * t170 - qJD(3) * t280 - t158 * t332 + t163 * t378 + t279 * t303;
t301 = -qJDD(1) * pkin(1) - t338 * qJ(2) + t361;
t344 = m(4) * t291 - mrSges(4,1) * t292 + t293 * mrSges(4,2) - t296 * t303 + t304 * t297 + t349;
t149 = -mrSges(3,1) * t301 + mrSges(3,3) * t295 - pkin(2) * t344 + pkin(7) * t363 + qJDD(1) * t359 + t336 * t153 + t333 * t154 - t329 * t370;
t151 = mrSges(3,2) * t301 - mrSges(3,3) * t294 - pkin(7) * t162 + qJDD(1) * t360 - t333 * t153 + t336 * t154 + t330 * t370;
t342 = -m(3) * t301 + mrSges(3,1) * t367 - t344 + (t324 * t338 + t373) * mrSges(3,3);
t353 = -mrSges(2,2) * t312 + qJ(2) * t364 + t330 * t149 + t329 * t151 + pkin(1) * (-qJDD(1) * t375 + t342) + mrSges(2,1) * t311 + Ifges(2,3) * qJDD(1);
t171 = t342 + (mrSges(2,1) - t375) * qJDD(1) - t338 * mrSges(2,2) + m(2) * t311;
t157 = t160 * t330 + t161 * t329;
t155 = m(2) * t312 - mrSges(2,1) * t338 - qJDD(1) * mrSges(2,2) + t364;
t152 = mrSges(2,1) * g(3) + (Ifges(2,6) - t358) * qJDD(1) - pkin(1) * t157 + t338 * Ifges(2,5) + mrSges(2,3) * t312 + t381;
t147 = -mrSges(2,2) * g(3) - mrSges(2,3) * t311 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t338 - qJ(2) * t157 - t149 * t329 + t151 * t330;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t337 * t147 - t334 * t152 - pkin(6) * (t155 * t334 + t171 * t337), t147, t151, t154, t163, -t282 * t250 - t347, t182; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t334 * t147 + t337 * t152 + pkin(6) * (t155 * t337 - t171 * t334), t152, t149, t153, t158, Ifges(6,4) * t323 - Ifges(6,2) * t241 + Ifges(6,6) * t240 - t326 * t248 + t282 * t252 - t345, t181; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t353, t353, qJDD(1) * t358 - t381, t382, t383, Ifges(6,5) * t323 - Ifges(6,6) * t241 + Ifges(6,3) * t240 + t326 * t250 + t283 * t252 - t346, t350;];
m_new  = t1;
