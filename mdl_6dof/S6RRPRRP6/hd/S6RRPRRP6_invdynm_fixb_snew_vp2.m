% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 18:11
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRRP6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP6_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP6_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP6_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP6_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP6_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:58:34
% EndTime: 2019-05-06 17:59:42
% DurationCPUTime: 34.08s
% Computational Cost: add. (542099->392), mult. (1415849->499), div. (0->0), fcn. (1112593->12), ass. (0->154)
t384 = -2 * qJD(3);
t338 = sin(pkin(11));
t340 = cos(pkin(11));
t344 = sin(qJ(2));
t347 = cos(qJ(2));
t339 = sin(pkin(6));
t371 = qJD(1) * t339;
t316 = (t338 * t344 - t340 * t347) * t371;
t369 = qJD(1) * qJD(2);
t325 = (qJDD(1) * t344 + t347 * t369) * t339;
t341 = cos(pkin(6));
t332 = qJDD(1) * t341 + qJDD(2);
t333 = qJD(1) * t341 + qJD(2);
t345 = sin(qJ(1));
t348 = cos(qJ(1));
t329 = t345 * g(1) - g(2) * t348;
t349 = qJD(1) ^ 2;
t381 = pkin(8) * t339;
t322 = qJDD(1) * pkin(1) + t349 * t381 + t329;
t330 = -g(1) * t348 - g(2) * t345;
t323 = -pkin(1) * t349 + qJDD(1) * t381 + t330;
t375 = t341 * t347;
t361 = t322 * t375 - t344 * t323;
t379 = t339 ^ 2 * t349;
t255 = t332 * pkin(2) - t325 * qJ(3) + (pkin(2) * t344 * t379 + (qJ(3) * qJD(1) * t333 - g(3)) * t339) * t347 + t361;
t376 = t341 * t344;
t378 = t339 * t344;
t288 = -g(3) * t378 + t322 * t376 + t347 * t323;
t366 = t344 * t371;
t319 = pkin(2) * t333 - qJ(3) * t366;
t326 = (qJDD(1) * t347 - t344 * t369) * t339;
t368 = t347 ^ 2 * t379;
t260 = -pkin(2) * t368 + qJ(3) * t326 - t319 * t333 + t288;
t317 = (t338 * t347 + t340 * t344) * t371;
t224 = t340 * t255 - t338 * t260 + t317 * t384;
t225 = t338 * t255 + t340 * t260 + t316 * t384;
t291 = pkin(3) * t316 - pkin(9) * t317;
t331 = t333 ^ 2;
t222 = -pkin(3) * t331 + pkin(9) * t332 - t291 * t316 + t225;
t307 = -t341 * g(3) - t339 * t322;
t274 = -t326 * pkin(2) - qJ(3) * t368 + t319 * t366 + qJDD(3) + t307;
t296 = -t325 * t338 + t326 * t340;
t297 = t325 * t340 + t326 * t338;
t228 = (t316 * t333 - t297) * pkin(9) + (t317 * t333 - t296) * pkin(3) + t274;
t343 = sin(qJ(4));
t346 = cos(qJ(4));
t218 = t346 * t222 + t343 * t228;
t300 = -t317 * t343 + t333 * t346;
t301 = t317 * t346 + t333 * t343;
t277 = -pkin(4) * t300 - pkin(10) * t301;
t295 = qJDD(4) - t296;
t315 = qJD(4) + t316;
t314 = t315 ^ 2;
t214 = -pkin(4) * t314 + pkin(10) * t295 + t277 * t300 + t218;
t221 = -t332 * pkin(3) - t331 * pkin(9) + t317 * t291 - t224;
t271 = -qJD(4) * t301 - t297 * t343 + t332 * t346;
t272 = qJD(4) * t300 + t297 * t346 + t332 * t343;
t216 = (-t300 * t315 - t272) * pkin(10) + (t301 * t315 - t271) * pkin(4) + t221;
t342 = sin(qJ(5));
t382 = cos(qJ(5));
t210 = -t342 * t214 + t216 * t382;
t211 = t214 * t382 + t342 * t216;
t280 = t301 * t382 + t342 * t315;
t235 = qJD(5) * t280 + t272 * t342 - t295 * t382;
t279 = t301 * t342 - t315 * t382;
t236 = -t279 * qJD(5) + t272 * t382 + t342 * t295;
t299 = qJD(5) - t300;
t238 = Ifges(7,5) * t280 + Ifges(7,6) * t299 + Ifges(7,3) * t279;
t241 = Ifges(6,4) * t280 - Ifges(6,2) * t279 + Ifges(6,6) * t299;
t243 = Ifges(6,1) * t280 - Ifges(6,4) * t279 + Ifges(6,5) * t299;
t249 = mrSges(7,1) * t279 - mrSges(7,3) * t280;
t270 = qJDD(5) - t271;
t248 = pkin(5) * t279 - qJ(6) * t280;
t298 = t299 ^ 2;
t206 = -pkin(5) * t298 + qJ(6) * t270 + 0.2e1 * qJD(6) * t299 - t248 * t279 + t211;
t208 = -t270 * pkin(5) - t298 * qJ(6) + t280 * t248 + qJDD(6) - t210;
t242 = Ifges(7,1) * t280 + Ifges(7,4) * t299 + Ifges(7,5) * t279;
t357 = mrSges(7,1) * t208 - mrSges(7,3) * t206 - Ifges(7,4) * t236 - Ifges(7,2) * t270 - Ifges(7,6) * t235 - t279 * t242;
t256 = -mrSges(7,2) * t279 + mrSges(7,3) * t299;
t360 = -m(7) * t208 + t270 * mrSges(7,1) + t299 * t256;
t259 = -mrSges(7,1) * t299 + mrSges(7,2) * t280;
t367 = m(7) * t206 + t270 * mrSges(7,3) + t299 * t259;
t383 = -(-t241 + t238) * t280 + mrSges(6,1) * t210 - mrSges(6,2) * t211 + Ifges(6,5) * t236 - Ifges(6,6) * t235 + Ifges(6,3) * t270 + pkin(5) * (-t236 * mrSges(7,2) - t280 * t249 + t360) + qJ(6) * (-t235 * mrSges(7,2) - t279 * t249 + t367) + t279 * t243 - t357;
t380 = -mrSges(6,3) - mrSges(7,2);
t377 = t339 * t347;
t289 = mrSges(4,1) * t316 + mrSges(4,2) * t317;
t303 = mrSges(4,1) * t333 - mrSges(4,3) * t317;
t258 = mrSges(6,1) * t299 - mrSges(6,3) * t280;
t372 = -mrSges(6,1) * t279 - mrSges(6,2) * t280 - t249;
t198 = m(6) * t211 - t270 * mrSges(6,2) + t235 * t380 - t299 * t258 + t279 * t372 + t367;
t257 = -mrSges(6,2) * t299 - mrSges(6,3) * t279;
t199 = m(6) * t210 + t270 * mrSges(6,1) + t236 * t380 + t299 * t257 + t280 * t372 + t360;
t194 = t198 * t382 - t199 * t342;
t276 = -mrSges(5,1) * t300 + mrSges(5,2) * t301;
t282 = mrSges(5,1) * t315 - mrSges(5,3) * t301;
t189 = m(5) * t218 - mrSges(5,2) * t295 + mrSges(5,3) * t271 + t276 * t300 - t282 * t315 + t194;
t217 = -t343 * t222 + t346 * t228;
t213 = -t295 * pkin(4) - t314 * pkin(10) + t301 * t277 - t217;
t209 = -0.2e1 * qJD(6) * t280 + (t279 * t299 - t236) * qJ(6) + (t280 * t299 + t235) * pkin(5) + t213;
t203 = m(7) * t209 + mrSges(7,1) * t235 - t236 * mrSges(7,3) + t256 * t279 - t280 * t259;
t200 = -m(6) * t213 - t235 * mrSges(6,1) - mrSges(6,2) * t236 - t279 * t257 - t258 * t280 - t203;
t281 = -mrSges(5,2) * t315 + mrSges(5,3) * t300;
t196 = m(5) * t217 + mrSges(5,1) * t295 - mrSges(5,3) * t272 - t276 * t301 + t281 * t315 + t200;
t363 = t346 * t189 - t196 * t343;
t181 = m(4) * t225 - mrSges(4,2) * t332 + mrSges(4,3) * t296 - t289 * t316 - t303 * t333 + t363;
t302 = -mrSges(4,2) * t333 - mrSges(4,3) * t316;
t193 = t342 * t198 + t199 * t382;
t352 = -m(5) * t221 + t271 * mrSges(5,1) - t272 * mrSges(5,2) + t300 * t281 - t301 * t282 - t193;
t186 = m(4) * t224 + t332 * mrSges(4,1) - t297 * mrSges(4,3) - t317 * t289 + t333 * t302 + t352;
t176 = t338 * t181 + t340 * t186;
t184 = t343 * t189 + t346 * t196;
t240 = Ifges(7,4) * t280 + Ifges(7,2) * t299 + Ifges(7,6) * t279;
t374 = -Ifges(6,5) * t280 + Ifges(6,6) * t279 - Ifges(6,3) * t299 - t240;
t365 = t347 * t371;
t287 = -g(3) * t377 + t361;
t321 = -mrSges(3,2) * t333 + mrSges(3,3) * t365;
t324 = (-mrSges(3,1) * t347 + mrSges(3,2) * t344) * t371;
t174 = m(3) * t287 + mrSges(3,1) * t332 - mrSges(3,3) * t325 + t321 * t333 - t324 * t366 + t176;
t320 = mrSges(3,1) * t333 - mrSges(3,3) * t366;
t364 = t340 * t181 - t186 * t338;
t175 = m(3) * t288 - mrSges(3,2) * t332 + mrSges(3,3) * t326 - t320 * t333 + t324 * t365 + t364;
t167 = -t174 * t344 + t347 * t175;
t354 = m(4) * t274 - t296 * mrSges(4,1) + t297 * mrSges(4,2) + t316 * t302 + t317 * t303 + t184;
t182 = m(3) * t307 - t326 * mrSges(3,1) + t325 * mrSges(3,2) + (t320 * t344 - t321 * t347) * t371 + t354;
t164 = t174 * t375 + t175 * t376 - t182 * t339;
t359 = -mrSges(7,1) * t209 + mrSges(7,2) * t206;
t356 = mrSges(7,2) * t208 - mrSges(7,3) * t209 + Ifges(7,1) * t236 + Ifges(7,4) * t270 + Ifges(7,5) * t235 + t299 * t238;
t191 = -mrSges(6,1) * t213 + mrSges(6,3) * t211 - pkin(5) * t203 + (t242 + t243) * t299 + t374 * t280 + (Ifges(6,6) - Ifges(7,6)) * t270 + (Ifges(6,4) - Ifges(7,5)) * t236 + (-Ifges(6,2) - Ifges(7,3)) * t235 + t359;
t192 = mrSges(6,2) * t213 - mrSges(6,3) * t210 + Ifges(6,1) * t236 - Ifges(6,4) * t235 + Ifges(6,5) * t270 - qJ(6) * t203 - t299 * t241 + t279 * t374 + t356;
t261 = Ifges(5,5) * t301 + Ifges(5,6) * t300 + Ifges(5,3) * t315;
t262 = Ifges(5,4) * t301 + Ifges(5,2) * t300 + Ifges(5,6) * t315;
t170 = mrSges(5,2) * t221 - mrSges(5,3) * t217 + Ifges(5,1) * t272 + Ifges(5,4) * t271 + Ifges(5,5) * t295 - pkin(10) * t193 - t342 * t191 + t192 * t382 + t300 * t261 - t315 * t262;
t263 = Ifges(5,1) * t301 + Ifges(5,4) * t300 + Ifges(5,5) * t315;
t178 = -mrSges(5,1) * t221 + mrSges(5,3) * t218 + Ifges(5,4) * t272 + Ifges(5,2) * t271 + Ifges(5,6) * t295 - pkin(4) * t193 - t301 * t261 + t315 * t263 - t383;
t283 = Ifges(4,5) * t317 - Ifges(4,6) * t316 + Ifges(4,3) * t333;
t284 = Ifges(4,4) * t317 - Ifges(4,2) * t316 + Ifges(4,6) * t333;
t160 = mrSges(4,2) * t274 - mrSges(4,3) * t224 + Ifges(4,1) * t297 + Ifges(4,4) * t296 + Ifges(4,5) * t332 - pkin(9) * t184 + t170 * t346 - t178 * t343 - t283 * t316 - t284 * t333;
t285 = Ifges(4,1) * t317 - Ifges(4,4) * t316 + Ifges(4,5) * t333;
t350 = mrSges(5,1) * t217 - mrSges(5,2) * t218 + Ifges(5,5) * t272 + Ifges(5,6) * t271 + Ifges(5,3) * t295 + pkin(4) * t200 + pkin(10) * t194 + t191 * t382 + t342 * t192 + t301 * t262 - t300 * t263;
t168 = -mrSges(4,1) * t274 + mrSges(4,3) * t225 + Ifges(4,4) * t297 + Ifges(4,2) * t296 + Ifges(4,6) * t332 - pkin(3) * t184 - t317 * t283 + t333 * t285 - t350;
t304 = Ifges(3,3) * t333 + (Ifges(3,5) * t344 + Ifges(3,6) * t347) * t371;
t306 = Ifges(3,5) * t333 + (Ifges(3,1) * t344 + Ifges(3,4) * t347) * t371;
t155 = -mrSges(3,1) * t307 + mrSges(3,3) * t288 + Ifges(3,4) * t325 + Ifges(3,2) * t326 + Ifges(3,6) * t332 - pkin(2) * t354 + qJ(3) * t364 + t338 * t160 + t340 * t168 - t304 * t366 + t333 * t306;
t305 = Ifges(3,6) * t333 + (Ifges(3,4) * t344 + Ifges(3,2) * t347) * t371;
t157 = mrSges(3,2) * t307 - mrSges(3,3) * t287 + Ifges(3,1) * t325 + Ifges(3,4) * t326 + Ifges(3,5) * t332 - qJ(3) * t176 + t160 * t340 - t168 * t338 + t304 * t365 - t305 * t333;
t353 = mrSges(4,1) * t224 - mrSges(4,2) * t225 + Ifges(4,5) * t297 + Ifges(4,6) * t296 + Ifges(4,3) * t332 + pkin(3) * t352 + pkin(9) * t363 + t343 * t170 + t346 * t178 + t317 * t284 + t316 * t285;
t159 = Ifges(3,5) * t325 + Ifges(3,6) * t326 + mrSges(3,1) * t287 - mrSges(3,2) * t288 + t353 + (t305 * t344 - t306 * t347) * t371 + Ifges(3,3) * t332 + pkin(2) * t176;
t355 = mrSges(2,1) * t329 - mrSges(2,2) * t330 + Ifges(2,3) * qJDD(1) + pkin(1) * t164 + t155 * t377 + t157 * t378 + t341 * t159 + t167 * t381;
t165 = m(2) * t330 - mrSges(2,1) * t349 - qJDD(1) * mrSges(2,2) + t167;
t163 = t341 * t182 + (t174 * t347 + t175 * t344) * t339;
t161 = m(2) * t329 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t349 + t164;
t153 = -mrSges(2,2) * g(3) - mrSges(2,3) * t329 + Ifges(2,5) * qJDD(1) - t349 * Ifges(2,6) - t344 * t155 + t347 * t157 + (-t163 * t339 - t164 * t341) * pkin(8);
t152 = mrSges(2,1) * g(3) + mrSges(2,3) * t330 + t349 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t163 - t339 * t159 + (pkin(8) * t167 + t155 * t347 + t157 * t344) * t341;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t348 * t153 - t345 * t152 - pkin(7) * (t161 * t348 + t165 * t345), t153, t157, t160, t170, t192, -t240 * t279 + t356; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t345 * t153 + t348 * t152 + pkin(7) * (-t161 * t345 + t165 * t348), t152, t155, t168, t178, t191, -t280 * t238 - t357; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t355, t355, t159, t353, t350, t383, Ifges(7,5) * t236 + Ifges(7,6) * t270 + Ifges(7,3) * t235 + t280 * t240 - t299 * t242 - t359;];
m_new  = t1;
