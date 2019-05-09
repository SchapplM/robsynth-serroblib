% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
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
% Datum: 2019-05-07 23:04
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRPR10_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR10_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR10_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR10_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR10_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR10_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR10_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 22:49:32
% EndTime: 2019-05-07 22:51:05
% DurationCPUTime: 33.10s
% Computational Cost: add. (591724->400), mult. (1268507->496), div. (0->0), fcn. (1016507->12), ass. (0->160)
t350 = sin(qJ(2));
t354 = cos(qJ(2));
t345 = sin(pkin(6));
t377 = qJD(1) * t345;
t328 = (-pkin(2) * t354 - pkin(9) * t350) * t377;
t346 = cos(pkin(6));
t341 = qJD(1) * t346 + qJD(2);
t339 = t341 ^ 2;
t340 = qJDD(1) * t346 + qJDD(2);
t376 = qJD(1) * t354;
t351 = sin(qJ(1));
t355 = cos(qJ(1));
t336 = t351 * g(1) - g(2) * t355;
t356 = qJD(1) ^ 2;
t388 = pkin(8) * t345;
t325 = qJDD(1) * pkin(1) + t356 * t388 + t336;
t337 = -g(1) * t355 - g(2) * t351;
t375 = qJDD(1) * t345;
t326 = -pkin(1) * t356 + pkin(8) * t375 + t337;
t382 = t346 * t350;
t378 = t325 * t382 + t354 * t326;
t275 = -pkin(2) * t339 + pkin(9) * t340 + (-g(3) * t350 + t328 * t376) * t345 + t378;
t329 = (qJD(2) * t376 + qJDD(1) * t350) * t345;
t374 = t350 * t377;
t330 = -qJD(2) * t374 + t354 * t375;
t387 = g(3) * t346;
t276 = -pkin(2) * t330 - pkin(9) * t329 - t387 + (-t325 + (pkin(2) * t350 - pkin(9) * t354) * t341 * qJD(1)) * t345;
t349 = sin(qJ(3));
t353 = cos(qJ(3));
t234 = -t275 * t349 + t353 * t276;
t317 = t341 * t353 - t349 * t374;
t294 = qJD(3) * t317 + t329 * t353 + t340 * t349;
t318 = t341 * t349 + t353 * t374;
t322 = qJDD(3) - t330;
t373 = t345 * t376;
t335 = qJD(3) - t373;
t223 = (t317 * t335 - t294) * pkin(10) + (t317 * t318 + t322) * pkin(3) + t234;
t235 = t353 * t275 + t349 * t276;
t293 = -qJD(3) * t318 - t329 * t349 + t340 * t353;
t305 = pkin(3) * t335 - pkin(10) * t318;
t316 = t317 ^ 2;
t226 = -pkin(3) * t316 + pkin(10) * t293 - t305 * t335 + t235;
t348 = sin(qJ(4));
t389 = cos(qJ(4));
t220 = t389 * t223 - t348 * t226;
t221 = t348 * t223 + t389 * t226;
t301 = t348 * t317 + t389 * t318;
t250 = qJD(4) * t301 - t389 * t293 + t294 * t348;
t300 = -t389 * t317 + t318 * t348;
t251 = -t300 * qJD(4) + t348 * t293 + t389 * t294;
t332 = -qJD(4) - t335;
t262 = Ifges(5,4) * t301 - Ifges(5,2) * t300 - Ifges(5,6) * t332;
t270 = -mrSges(6,2) * t300 - mrSges(6,3) * t301;
t280 = mrSges(6,1) * t300 + mrSges(6,3) * t332;
t321 = qJDD(4) + t322;
t268 = pkin(4) * t300 - qJ(5) * t301;
t331 = t332 ^ 2;
t216 = -t321 * pkin(4) - t331 * qJ(5) + t301 * t268 + qJDD(5) - t220;
t385 = t300 * t332;
t210 = (t300 * t301 - t321) * pkin(11) + (t251 - t385) * pkin(5) + t216;
t284 = pkin(5) * t301 + pkin(11) * t332;
t299 = t300 ^ 2;
t381 = t346 * t354;
t383 = t345 * t354;
t295 = -g(3) * t383 + t325 * t381 - t350 * t326;
t274 = -pkin(2) * t340 - pkin(9) * t339 + t328 * t374 - t295;
t233 = -pkin(3) * t293 - pkin(10) * t316 + t318 * t305 + t274;
t390 = -2 * qJD(5);
t359 = (-t251 - t385) * qJ(5) + t233 + (-t332 * pkin(4) + t390) * t301;
t213 = t359 + (pkin(4) + pkin(11)) * t250 - pkin(5) * t299 - t284 * t301;
t347 = sin(qJ(6));
t352 = cos(qJ(6));
t207 = t210 * t352 - t213 * t347;
t278 = t300 * t352 + t332 * t347;
t231 = qJD(6) * t278 + t250 * t347 + t321 * t352;
t249 = qJDD(6) + t251;
t279 = t300 * t347 - t332 * t352;
t256 = -mrSges(7,1) * t278 + mrSges(7,2) * t279;
t298 = qJD(6) + t301;
t257 = -mrSges(7,2) * t298 + mrSges(7,3) * t278;
t204 = m(7) * t207 + mrSges(7,1) * t249 - mrSges(7,3) * t231 - t256 * t279 + t257 * t298;
t208 = t210 * t347 + t213 * t352;
t230 = -qJD(6) * t279 + t250 * t352 - t321 * t347;
t258 = mrSges(7,1) * t298 - mrSges(7,3) * t279;
t205 = m(7) * t208 - mrSges(7,2) * t249 + mrSges(7,3) * t230 + t256 * t278 - t258 * t298;
t192 = t204 * t352 + t205 * t347;
t368 = -pkin(4) * t331 + qJ(5) * t321 - t268 * t300 + t221;
t212 = -pkin(5) * t250 - pkin(11) * t299 + (t390 - t284) * t332 + t368;
t236 = Ifges(7,5) * t279 + Ifges(7,6) * t278 + Ifges(7,3) * t298;
t238 = Ifges(7,1) * t279 + Ifges(7,4) * t278 + Ifges(7,5) * t298;
t195 = -mrSges(7,1) * t212 + mrSges(7,3) * t208 + Ifges(7,4) * t231 + Ifges(7,2) * t230 + Ifges(7,6) * t249 - t236 * t279 + t238 * t298;
t237 = Ifges(7,4) * t279 + Ifges(7,2) * t278 + Ifges(7,6) * t298;
t196 = mrSges(7,2) * t212 - mrSges(7,3) * t207 + Ifges(7,1) * t231 + Ifges(7,4) * t230 + Ifges(7,5) * t249 + t236 * t278 - t237 * t298;
t214 = 0.2e1 * qJD(5) * t332 - t368;
t259 = -Ifges(6,5) * t332 - Ifges(6,6) * t301 + Ifges(6,3) * t300;
t364 = -mrSges(6,2) * t216 + mrSges(6,3) * t214 - Ifges(6,1) * t321 + Ifges(6,4) * t251 - Ifges(6,5) * t250 + pkin(11) * t192 + t347 * t195 - t352 * t196 + t301 * t259;
t209 = -m(7) * t212 + mrSges(7,1) * t230 - t231 * mrSges(7,2) + t257 * t278 - t279 * t258;
t281 = mrSges(6,1) * t301 - mrSges(6,2) * t332;
t365 = -m(6) * t214 + t321 * mrSges(6,3) - t332 * t281 - t209;
t369 = -m(6) * t216 - t251 * mrSges(6,1) - t301 * t270 - t192;
t261 = -Ifges(6,4) * t332 - Ifges(6,2) * t301 + Ifges(6,6) * t300;
t379 = -Ifges(5,1) * t301 + Ifges(5,4) * t300 + Ifges(5,5) * t332 + t261;
t392 = -t379 * t300 - mrSges(5,2) * t221 + pkin(4) * (-mrSges(6,2) * t321 + t280 * t332 + t369) + qJ(5) * (-mrSges(6,1) * t250 - t270 * t300 + t365) + mrSges(5,1) * t220 - Ifges(5,6) * t250 + Ifges(5,5) * t251 + t301 * t262 + Ifges(5,3) * t321 - t364;
t269 = mrSges(5,1) * t300 + mrSges(5,2) * t301;
t282 = mrSges(5,2) * t332 - mrSges(5,3) * t300;
t188 = m(5) * t220 - mrSges(5,3) * t251 - t269 * t301 + (t280 - t282) * t332 + (mrSges(5,1) - mrSges(6,2)) * t321 + t369;
t283 = -mrSges(5,1) * t332 - mrSges(5,3) * t301;
t199 = m(5) * t221 - mrSges(5,2) * t321 + t283 * t332 + (-t269 - t270) * t300 + (-mrSges(5,3) - mrSges(6,1)) * t250 + t365;
t184 = t389 * t188 + t348 * t199;
t288 = Ifges(4,4) * t318 + Ifges(4,2) * t317 + Ifges(4,6) * t335;
t289 = Ifges(4,1) * t318 + Ifges(4,4) * t317 + Ifges(4,5) * t335;
t391 = mrSges(4,1) * t234 - mrSges(4,2) * t235 + Ifges(4,5) * t294 + Ifges(4,6) * t293 + Ifges(4,3) * t322 + pkin(3) * t184 + t318 * t288 - t317 * t289 + t392;
t386 = Ifges(5,4) + Ifges(6,6);
t384 = t345 * t350;
t302 = -mrSges(4,1) * t317 + mrSges(4,2) * t318;
t303 = -mrSges(4,2) * t335 + mrSges(4,3) * t317;
t182 = m(4) * t234 + mrSges(4,1) * t322 - mrSges(4,3) * t294 - t302 * t318 + t303 * t335 + t184;
t304 = mrSges(4,1) * t335 - mrSges(4,3) * t318;
t370 = -t188 * t348 + t389 * t199;
t183 = m(4) * t235 - mrSges(4,2) * t322 + mrSges(4,3) * t293 + t302 * t317 - t304 * t335 + t370;
t177 = t353 * t182 + t349 * t183;
t193 = -t347 * t204 + t352 * t205;
t263 = -Ifges(6,1) * t332 - Ifges(6,4) * t301 + Ifges(6,5) * t300;
t380 = -Ifges(5,5) * t301 + Ifges(5,6) * t300 + Ifges(5,3) * t332 - t263;
t296 = -g(3) * t384 + t378;
t323 = mrSges(3,1) * t341 - mrSges(3,3) * t374;
t327 = (-mrSges(3,1) * t354 + mrSges(3,2) * t350) * t377;
t371 = -t182 * t349 + t353 * t183;
t175 = m(3) * t296 - mrSges(3,2) * t340 + mrSges(3,3) * t330 - t323 * t341 + t327 * t373 + t371;
t324 = -mrSges(3,2) * t341 + mrSges(3,3) * t373;
t218 = pkin(4) * t250 + t359;
t189 = m(6) * t218 - t250 * mrSges(6,2) - t251 * mrSges(6,3) - t300 * t280 - t301 * t281 + t193;
t362 = m(5) * t233 + t250 * mrSges(5,1) + t251 * mrSges(5,2) + t300 * t282 + t301 * t283 + t189;
t358 = -m(4) * t274 + t293 * mrSges(4,1) - t294 * mrSges(4,2) + t317 * t303 - t318 * t304 - t362;
t186 = m(3) * t295 + t340 * mrSges(3,1) - t329 * mrSges(3,3) + t341 * t324 - t327 * t374 + t358;
t171 = t354 * t175 - t186 * t350;
t309 = -t325 * t345 - t387;
t176 = m(3) * t309 - mrSges(3,1) * t330 + mrSges(3,2) * t329 + (t323 * t350 - t324 * t354) * t377 + t177;
t168 = t175 * t382 - t176 * t345 + t186 * t381;
t363 = -mrSges(6,1) * t214 + mrSges(6,2) * t218 - pkin(5) * t209 - pkin(11) * t193 - t352 * t195 - t347 * t196;
t172 = -mrSges(5,1) * t233 + mrSges(5,3) * t221 - pkin(4) * t189 + t379 * t332 + (Ifges(5,6) - Ifges(6,5)) * t321 + t380 * t301 + t386 * t251 + (-Ifges(5,2) - Ifges(6,3)) * t250 + t363;
t366 = mrSges(7,1) * t207 - mrSges(7,2) * t208 + Ifges(7,5) * t231 + Ifges(7,6) * t230 + Ifges(7,3) * t249 + t279 * t237 - t278 * t238;
t361 = mrSges(6,1) * t216 - mrSges(6,3) * t218 + pkin(5) * t192 + t366;
t178 = (t262 - t259) * t332 + (Ifges(5,5) - Ifges(6,4)) * t321 + t380 * t300 + (Ifges(5,1) + Ifges(6,2)) * t251 - t386 * t250 + mrSges(5,2) * t233 - mrSges(5,3) * t220 + t361 - qJ(5) * t189;
t287 = Ifges(4,5) * t318 + Ifges(4,6) * t317 + Ifges(4,3) * t335;
t161 = -mrSges(4,1) * t274 + mrSges(4,3) * t235 + Ifges(4,4) * t294 + Ifges(4,2) * t293 + Ifges(4,6) * t322 - pkin(3) * t362 + pkin(10) * t370 + t389 * t172 + t348 * t178 - t318 * t287 + t335 * t289;
t164 = mrSges(4,2) * t274 - mrSges(4,3) * t234 + Ifges(4,1) * t294 + Ifges(4,4) * t293 + Ifges(4,5) * t322 - pkin(10) * t184 - t348 * t172 + t389 * t178 + t317 * t287 - t335 * t288;
t307 = Ifges(3,6) * t341 + (Ifges(3,4) * t350 + Ifges(3,2) * t354) * t377;
t308 = Ifges(3,5) * t341 + (Ifges(3,1) * t350 + Ifges(3,4) * t354) * t377;
t158 = Ifges(3,5) * t329 + Ifges(3,6) * t330 + Ifges(3,3) * t340 + mrSges(3,1) * t295 - mrSges(3,2) * t296 + t349 * t164 + t353 * t161 + pkin(2) * t358 + pkin(9) * t371 + (t307 * t350 - t308 * t354) * t377;
t306 = Ifges(3,3) * t341 + (Ifges(3,5) * t350 + Ifges(3,6) * t354) * t377;
t160 = mrSges(3,2) * t309 - mrSges(3,3) * t295 + Ifges(3,1) * t329 + Ifges(3,4) * t330 + Ifges(3,5) * t340 - pkin(9) * t177 - t161 * t349 + t164 * t353 + t306 * t373 - t307 * t341;
t163 = -mrSges(3,1) * t309 + mrSges(3,3) * t296 + Ifges(3,4) * t329 + Ifges(3,2) * t330 + Ifges(3,6) * t340 - pkin(2) * t177 - t306 * t374 + t341 * t308 - t391;
t367 = mrSges(2,1) * t336 - mrSges(2,2) * t337 + Ifges(2,3) * qJDD(1) + pkin(1) * t168 + t346 * t158 + t160 * t384 + t163 * t383 + t171 * t388;
t169 = m(2) * t337 - mrSges(2,1) * t356 - qJDD(1) * mrSges(2,2) + t171;
t167 = t176 * t346 + (t175 * t350 + t186 * t354) * t345;
t165 = m(2) * t336 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t356 + t168;
t156 = -mrSges(2,2) * g(3) - mrSges(2,3) * t336 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t356 + t160 * t354 - t163 * t350 + (-t167 * t345 - t168 * t346) * pkin(8);
t155 = mrSges(2,1) * g(3) + mrSges(2,3) * t337 + Ifges(2,5) * t356 + Ifges(2,6) * qJDD(1) - pkin(1) * t167 - t158 * t345 + (pkin(8) * t171 + t160 * t350 + t163 * t354) * t346;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t355 * t156 - t351 * t155 - pkin(7) * (t165 * t355 + t169 * t351), t156, t160, t164, t178, -t300 * t261 - t364, t196; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t351 * t156 + t355 * t155 + pkin(7) * (-t165 * t351 + t169 * t355), t155, t163, t161, t172, Ifges(6,4) * t321 - Ifges(6,2) * t251 + Ifges(6,6) * t250 + t332 * t259 + t300 * t263 - t361, t195; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t367, t367, t158, t391, t392, Ifges(6,5) * t321 - Ifges(6,6) * t251 + Ifges(6,3) * t250 - t332 * t261 + t301 * t263 - t363, t366;];
m_new  = t1;
