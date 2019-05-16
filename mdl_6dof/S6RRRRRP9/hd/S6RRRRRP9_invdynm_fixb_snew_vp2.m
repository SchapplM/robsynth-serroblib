% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRRP9
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 06:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRRP9_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP9_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP9_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP9_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP9_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP9_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP9_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 05:57:20
% EndTime: 2019-05-08 05:59:14
% DurationCPUTime: 43.78s
% Computational Cost: add. (794548->393), mult. (1689903->498), div. (0->0), fcn. (1355695->12), ass. (0->156)
t345 = cos(pkin(6));
t340 = qJD(1) * t345 + qJD(2);
t348 = sin(qJ(3));
t353 = cos(qJ(3));
t349 = sin(qJ(2));
t344 = sin(pkin(6));
t377 = qJD(1) * t344;
t372 = t349 * t377;
t318 = t340 * t353 - t348 * t372;
t354 = cos(qJ(2));
t375 = qJD(1) * qJD(2);
t330 = (qJDD(1) * t349 + t354 * t375) * t344;
t339 = qJDD(1) * t345 + qJDD(2);
t299 = qJD(3) * t318 + t330 * t353 + t339 * t348;
t319 = t340 * t348 + t353 * t372;
t376 = qJD(1) * t354;
t371 = t344 * t376;
t335 = qJD(3) - t371;
t347 = sin(qJ(4));
t352 = cos(qJ(4));
t306 = t319 * t352 + t335 * t347;
t331 = (-qJDD(1) * t354 + t349 * t375) * t344;
t323 = qJDD(3) + t331;
t260 = -qJD(4) * t306 - t299 * t347 + t323 * t352;
t305 = -t319 * t347 + t335 * t352;
t261 = qJD(4) * t305 + t299 * t352 + t323 * t347;
t346 = sin(qJ(5));
t351 = cos(qJ(5));
t279 = t305 * t351 - t306 * t346;
t233 = qJD(5) * t279 + t260 * t346 + t261 * t351;
t280 = t305 * t346 + t306 * t351;
t256 = -mrSges(7,1) * t279 + mrSges(7,2) * t280;
t329 = (-pkin(2) * t354 - pkin(9) * t349) * t377;
t338 = t340 ^ 2;
t350 = sin(qJ(1));
t355 = cos(qJ(1));
t336 = t350 * g(1) - g(2) * t355;
t356 = qJD(1) ^ 2;
t385 = pkin(8) * t344;
t326 = qJDD(1) * pkin(1) + t356 * t385 + t336;
t337 = -g(1) * t355 - g(2) * t350;
t327 = -pkin(1) * t356 + qJDD(1) * t385 + t337;
t381 = t345 * t349;
t378 = t326 * t381 + t354 * t327;
t276 = -pkin(2) * t338 + pkin(9) * t339 + (-g(3) * t349 + t329 * t376) * t344 + t378;
t384 = g(3) * t345;
t277 = pkin(2) * t331 - pkin(9) * t330 - t384 + (-t326 + (pkin(2) * t349 - pkin(9) * t354) * t340 * qJD(1)) * t344;
t246 = t353 * t276 + t348 * t277;
t303 = -pkin(3) * t318 - pkin(10) * t319;
t333 = t335 ^ 2;
t236 = -pkin(3) * t333 + pkin(10) * t323 + t303 * t318 + t246;
t380 = t345 * t354;
t382 = t344 * t354;
t300 = -g(3) * t382 + t326 * t380 - t349 * t327;
t275 = -pkin(2) * t339 - pkin(9) * t338 + t329 * t372 - t300;
t298 = -t319 * qJD(3) - t348 * t330 + t339 * t353;
t242 = (-t318 * t335 - t299) * pkin(10) + (t319 * t335 - t298) * pkin(3) + t275;
t217 = -t236 * t347 + t352 * t242;
t296 = qJDD(4) - t298;
t317 = qJD(4) - t318;
t214 = (t305 * t317 - t261) * pkin(11) + (t305 * t306 + t296) * pkin(4) + t217;
t218 = t352 * t236 + t347 * t242;
t285 = pkin(4) * t317 - pkin(11) * t306;
t304 = t305 ^ 2;
t216 = -pkin(4) * t304 + pkin(11) * t260 - t285 * t317 + t218;
t207 = t351 * t214 - t216 * t346;
t291 = qJDD(5) + t296;
t315 = qJD(5) + t317;
t202 = -0.2e1 * qJD(6) * t280 + (t279 * t315 - t233) * qJ(6) + (t279 * t280 + t291) * pkin(5) + t207;
t262 = -mrSges(7,2) * t315 + mrSges(7,3) * t279;
t374 = m(7) * t202 + t291 * mrSges(7,1) + t315 * t262;
t199 = -mrSges(7,3) * t233 - t256 * t280 + t374;
t208 = t346 * t214 + t351 * t216;
t232 = -qJD(5) * t280 + t260 * t351 - t261 * t346;
t250 = Ifges(6,4) * t280 + Ifges(6,2) * t279 + Ifges(6,6) * t315;
t251 = Ifges(7,1) * t280 + Ifges(7,4) * t279 + Ifges(7,5) * t315;
t252 = Ifges(6,1) * t280 + Ifges(6,4) * t279 + Ifges(6,5) * t315;
t264 = pkin(5) * t315 - qJ(6) * t280;
t278 = t279 ^ 2;
t205 = -pkin(5) * t278 + qJ(6) * t232 + 0.2e1 * qJD(6) * t279 - t264 * t315 + t208;
t249 = Ifges(7,4) * t280 + Ifges(7,2) * t279 + Ifges(7,6) * t315;
t364 = -mrSges(7,1) * t202 + mrSges(7,2) * t205 - Ifges(7,5) * t233 - Ifges(7,6) * t232 - Ifges(7,3) * t291 - t280 * t249;
t387 = mrSges(6,1) * t207 - mrSges(6,2) * t208 + Ifges(6,5) * t233 + Ifges(6,6) * t232 + Ifges(6,3) * t291 + pkin(5) * t199 + t280 * t250 - t364 + (-t252 - t251) * t279;
t257 = -mrSges(6,1) * t279 + mrSges(6,2) * t280;
t263 = -mrSges(6,2) * t315 + mrSges(6,3) * t279;
t191 = m(6) * t207 + mrSges(6,1) * t291 + t263 * t315 + (-t256 - t257) * t280 + (-mrSges(6,3) - mrSges(7,3)) * t233 + t374;
t265 = mrSges(7,1) * t315 - mrSges(7,3) * t280;
t266 = mrSges(6,1) * t315 - mrSges(6,3) * t280;
t373 = m(7) * t205 + t232 * mrSges(7,3) + t279 * t256;
t194 = m(6) * t208 + mrSges(6,3) * t232 + t257 * t279 + (-t265 - t266) * t315 + (-mrSges(6,2) - mrSges(7,2)) * t291 + t373;
t189 = t351 * t191 + t346 * t194;
t268 = Ifges(5,4) * t306 + Ifges(5,2) * t305 + Ifges(5,6) * t317;
t269 = Ifges(5,1) * t306 + Ifges(5,4) * t305 + Ifges(5,5) * t317;
t386 = mrSges(5,1) * t217 - mrSges(5,2) * t218 + Ifges(5,5) * t261 + Ifges(5,6) * t260 + Ifges(5,3) * t296 + pkin(4) * t189 + t306 * t268 - t305 * t269 + t387;
t383 = t344 * t349;
t281 = -mrSges(5,1) * t305 + mrSges(5,2) * t306;
t283 = -mrSges(5,2) * t317 + mrSges(5,3) * t305;
t186 = m(5) * t217 + mrSges(5,1) * t296 - mrSges(5,3) * t261 - t281 * t306 + t283 * t317 + t189;
t284 = mrSges(5,1) * t317 - mrSges(5,3) * t306;
t368 = -t191 * t346 + t351 * t194;
t187 = m(5) * t218 - mrSges(5,2) * t296 + mrSges(5,3) * t260 + t281 * t305 - t284 * t317 + t368;
t183 = -t186 * t347 + t352 * t187;
t302 = -mrSges(4,1) * t318 + mrSges(4,2) * t319;
t308 = mrSges(4,1) * t335 - mrSges(4,3) * t319;
t181 = m(4) * t246 - mrSges(4,2) * t323 + mrSges(4,3) * t298 + t302 * t318 - t308 * t335 + t183;
t245 = -t348 * t276 + t277 * t353;
t235 = -pkin(3) * t323 - pkin(10) * t333 + t319 * t303 - t245;
t219 = -pkin(4) * t260 - pkin(11) * t304 + t306 * t285 + t235;
t211 = -pkin(5) * t232 - qJ(6) * t278 + t264 * t280 + qJDD(6) + t219;
t367 = m(7) * t211 - t232 * mrSges(7,1) + t233 * mrSges(7,2) - t279 * t262 + t280 * t265;
t361 = m(6) * t219 - t232 * mrSges(6,1) + mrSges(6,2) * t233 - t279 * t263 + t266 * t280 + t367;
t197 = -m(5) * t235 + t260 * mrSges(5,1) - mrSges(5,2) * t261 + t305 * t283 - t284 * t306 - t361;
t307 = -mrSges(4,2) * t335 + mrSges(4,3) * t318;
t196 = m(4) * t245 + mrSges(4,1) * t323 - mrSges(4,3) * t299 - t302 * t319 + t307 * t335 + t197;
t176 = t348 * t181 + t353 * t196;
t301 = -g(3) * t383 + t378;
t324 = mrSges(3,1) * t340 - mrSges(3,3) * t372;
t328 = (-mrSges(3,1) * t354 + mrSges(3,2) * t349) * t377;
t369 = t353 * t181 - t196 * t348;
t174 = m(3) * t301 - mrSges(3,2) * t339 - mrSges(3,3) * t331 - t324 * t340 + t328 * t371 + t369;
t325 = -mrSges(3,2) * t340 + mrSges(3,3) * t371;
t182 = t186 * t352 + t187 * t347;
t360 = -m(4) * t275 + t298 * mrSges(4,1) - mrSges(4,2) * t299 + t318 * t307 - t308 * t319 - t182;
t178 = m(3) * t300 + mrSges(3,1) * t339 - mrSges(3,3) * t330 + t325 * t340 - t328 * t372 + t360;
t168 = t354 * t174 - t178 * t349;
t312 = -t326 * t344 - t384;
t175 = m(3) * t312 + mrSges(3,1) * t331 + mrSges(3,2) * t330 + (t324 * t349 - t325 * t354) * t377 + t176;
t165 = t174 * t381 - t175 * t344 + t178 * t380;
t365 = -mrSges(7,1) * t211 + mrSges(7,3) * t205 + Ifges(7,4) * t233 + Ifges(7,2) * t232 + Ifges(7,6) * t291 + t315 * t251;
t247 = Ifges(7,5) * t280 + Ifges(7,6) * t279 + Ifges(7,3) * t315;
t363 = mrSges(7,2) * t211 - mrSges(7,3) * t202 + Ifges(7,1) * t233 + Ifges(7,4) * t232 + Ifges(7,5) * t291 + t279 * t247;
t248 = Ifges(6,5) * t280 + Ifges(6,6) * t279 + Ifges(6,3) * t315;
t184 = Ifges(6,4) * t233 + Ifges(6,2) * t232 + Ifges(6,6) * t291 + t315 * t252 - mrSges(6,1) * t219 + mrSges(6,3) * t208 - pkin(5) * t367 + qJ(6) * (-t291 * mrSges(7,2) - t315 * t265 + t373) + (-t248 - t247) * t280 + t365;
t188 = mrSges(6,2) * t219 - mrSges(6,3) * t207 + Ifges(6,1) * t233 + Ifges(6,4) * t232 + Ifges(6,5) * t291 - qJ(6) * t199 + t248 * t279 + (-t249 - t250) * t315 + t363;
t267 = Ifges(5,5) * t306 + Ifges(5,6) * t305 + Ifges(5,3) * t317;
t170 = -mrSges(5,1) * t235 + mrSges(5,3) * t218 + Ifges(5,4) * t261 + Ifges(5,2) * t260 + Ifges(5,6) * t296 - pkin(4) * t361 + pkin(11) * t368 + t351 * t184 + t346 * t188 - t306 * t267 + t317 * t269;
t171 = mrSges(5,2) * t235 - mrSges(5,3) * t217 + Ifges(5,1) * t261 + Ifges(5,4) * t260 + Ifges(5,5) * t296 - pkin(11) * t189 - t184 * t346 + t188 * t351 + t267 * t305 - t268 * t317;
t292 = Ifges(4,5) * t319 + Ifges(4,6) * t318 + Ifges(4,3) * t335;
t293 = Ifges(4,4) * t319 + Ifges(4,2) * t318 + Ifges(4,6) * t335;
t161 = mrSges(4,2) * t275 - mrSges(4,3) * t245 + Ifges(4,1) * t299 + Ifges(4,4) * t298 + Ifges(4,5) * t323 - pkin(10) * t182 - t170 * t347 + t171 * t352 + t292 * t318 - t293 * t335;
t294 = Ifges(4,1) * t319 + Ifges(4,4) * t318 + Ifges(4,5) * t335;
t169 = -mrSges(4,1) * t275 + mrSges(4,3) * t246 + Ifges(4,4) * t299 + Ifges(4,2) * t298 + Ifges(4,6) * t323 - pkin(3) * t182 - t319 * t292 + t335 * t294 - t386;
t310 = Ifges(3,6) * t340 + (Ifges(3,4) * t349 + Ifges(3,2) * t354) * t377;
t311 = Ifges(3,5) * t340 + (Ifges(3,1) * t349 + Ifges(3,4) * t354) * t377;
t156 = Ifges(3,5) * t330 - Ifges(3,6) * t331 + Ifges(3,3) * t339 + mrSges(3,1) * t300 - mrSges(3,2) * t301 + t348 * t161 + t353 * t169 + pkin(2) * t360 + pkin(9) * t369 + (t310 * t349 - t311 * t354) * t377;
t309 = Ifges(3,3) * t340 + (Ifges(3,5) * t349 + Ifges(3,6) * t354) * t377;
t158 = mrSges(3,2) * t312 - mrSges(3,3) * t300 + Ifges(3,1) * t330 - Ifges(3,4) * t331 + Ifges(3,5) * t339 - pkin(9) * t176 + t161 * t353 - t169 * t348 + t309 * t371 - t310 * t340;
t358 = mrSges(4,1) * t245 - mrSges(4,2) * t246 + Ifges(4,5) * t299 + Ifges(4,6) * t298 + Ifges(4,3) * t323 + pkin(3) * t197 + pkin(10) * t183 + t352 * t170 + t347 * t171 + t319 * t293 - t318 * t294;
t160 = -mrSges(3,1) * t312 + mrSges(3,3) * t301 + Ifges(3,4) * t330 - Ifges(3,2) * t331 + Ifges(3,6) * t339 - pkin(2) * t176 - t309 * t372 + t340 * t311 - t358;
t362 = mrSges(2,1) * t336 - mrSges(2,2) * t337 + Ifges(2,3) * qJDD(1) + pkin(1) * t165 + t345 * t156 + t158 * t383 + t160 * t382 + t168 * t385;
t166 = m(2) * t337 - mrSges(2,1) * t356 - qJDD(1) * mrSges(2,2) + t168;
t164 = t175 * t345 + (t174 * t349 + t178 * t354) * t344;
t162 = m(2) * t336 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t356 + t165;
t154 = -mrSges(2,2) * g(3) - mrSges(2,3) * t336 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t356 + t158 * t354 - t160 * t349 + (-t164 * t344 - t165 * t345) * pkin(8);
t153 = mrSges(2,1) * g(3) + mrSges(2,3) * t337 + Ifges(2,5) * t356 + Ifges(2,6) * qJDD(1) - pkin(1) * t164 - t156 * t344 + (pkin(8) * t168 + t158 * t349 + t160 * t354) * t345;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t355 * t154 - t350 * t153 - pkin(7) * (t162 * t355 + t166 * t350), t154, t158, t161, t171, t188, -t249 * t315 + t363; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t350 * t154 + t355 * t153 + pkin(7) * (-t162 * t350 + t166 * t355), t153, t160, t169, t170, t184, -t280 * t247 + t365; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t362, t362, t156, t358, t386, t387, -t279 * t251 - t364;];
m_new  = t1;
