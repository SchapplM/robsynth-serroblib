% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-05-07 07:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPPR10_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR10_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR10_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR10_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR10_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR10_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR10_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 06:58:37
% EndTime: 2019-05-07 06:59:36
% DurationCPUTime: 29.37s
% Computational Cost: add. (534463->402), mult. (1166921->499), div. (0->0), fcn. (905081->12), ass. (0->157)
t346 = sin(pkin(6));
t351 = sin(qJ(2));
t354 = cos(qJ(2));
t376 = qJD(1) * qJD(2);
t329 = (-qJDD(1) * t354 + t351 * t376) * t346;
t379 = qJD(1) * t346;
t327 = (-pkin(2) * t354 - pkin(9) * t351) * t379;
t348 = cos(pkin(6));
t341 = qJD(1) * t348 + qJD(2);
t339 = t341 ^ 2;
t340 = qJDD(1) * t348 + qJDD(2);
t378 = qJD(1) * t354;
t352 = sin(qJ(1));
t355 = cos(qJ(1));
t337 = g(1) * t352 - g(2) * t355;
t356 = qJD(1) ^ 2;
t390 = pkin(8) * t346;
t324 = qJDD(1) * pkin(1) + t356 * t390 + t337;
t338 = -g(1) * t355 - g(2) * t352;
t325 = -pkin(1) * t356 + qJDD(1) * t390 + t338;
t384 = t348 * t351;
t380 = t324 * t384 + t325 * t354;
t252 = -t339 * pkin(2) + t340 * pkin(9) + (-g(3) * t351 + t327 * t378) * t346 + t380;
t328 = (qJDD(1) * t351 + t354 * t376) * t346;
t389 = t348 * g(3);
t253 = t329 * pkin(2) - t328 * pkin(9) - t389 + (-t324 + (pkin(2) * t351 - pkin(9) * t354) * t341 * qJD(1)) * t346;
t350 = sin(qJ(3));
t391 = cos(qJ(3));
t235 = t252 * t391 + t253 * t350;
t375 = t351 * t379;
t314 = -t341 * t391 + t350 * t375;
t315 = t341 * t350 + t375 * t391;
t289 = pkin(3) * t314 - qJ(4) * t315;
t321 = qJDD(3) + t329;
t374 = t346 * t378;
t335 = -qJD(3) + t374;
t334 = t335 ^ 2;
t224 = pkin(3) * t334 - qJ(4) * t321 + 0.2e1 * qJD(4) * t335 + t289 * t314 - t235;
t234 = -t252 * t350 + t253 * t391;
t273 = Ifges(4,4) * t315 - Ifges(4,2) * t314 - Ifges(4,6) * t335;
t285 = qJD(3) * t315 + t328 * t350 - t340 * t391;
t286 = -qJD(3) * t314 + t328 * t391 + t340 * t350;
t291 = -mrSges(5,2) * t314 - mrSges(5,3) * t315;
t302 = mrSges(5,1) * t314 + mrSges(5,3) * t335;
t301 = pkin(4) * t315 + qJ(5) * t335;
t313 = t314 ^ 2;
t222 = -pkin(4) * t285 - qJ(5) * t313 - t301 * t335 + qJDD(5) - t224;
t345 = sin(pkin(11));
t347 = cos(pkin(11));
t259 = t285 * t347 - t321 * t345;
t260 = t285 * t345 + t321 * t347;
t297 = t314 * t347 + t335 * t345;
t266 = -mrSges(6,2) * t315 + mrSges(6,3) * t297;
t298 = t314 * t345 - t335 * t347;
t267 = mrSges(6,1) * t315 - mrSges(6,3) * t298;
t268 = pkin(5) * t315 - pkin(10) * t298;
t296 = t297 ^ 2;
t216 = -pkin(5) * t259 - pkin(10) * t296 + t268 * t298 + t222;
t349 = sin(qJ(6));
t353 = cos(qJ(6));
t255 = t297 * t349 + t298 * t353;
t232 = -qJD(6) * t255 + t259 * t353 - t260 * t349;
t254 = t297 * t353 - t298 * t349;
t233 = qJD(6) * t254 + t259 * t349 + t260 * t353;
t312 = qJD(6) + t315;
t242 = -mrSges(7,2) * t312 + mrSges(7,3) * t254;
t243 = mrSges(7,1) * t312 - mrSges(7,3) * t255;
t368 = m(7) * t216 - mrSges(7,1) * t232 + mrSges(7,2) * t233 - t242 * t254 + t243 * t255;
t206 = -m(6) * t222 + mrSges(6,1) * t259 - mrSges(6,2) * t260 + t266 * t297 - t267 * t298 - t368;
t303 = mrSges(5,1) * t315 - mrSges(5,2) * t335;
t360 = -m(5) * t224 + mrSges(5,3) * t321 - t303 * t335 - t206;
t226 = -t321 * pkin(3) - t334 * qJ(4) + t289 * t315 + qJDD(4) - t234;
t387 = t314 * t335;
t219 = (t314 * t315 - t321) * qJ(5) + (t286 - t387) * pkin(4) + t226;
t383 = t348 * t354;
t385 = t346 * t354;
t287 = -g(3) * t385 + t324 * t383 - t325 * t351;
t251 = -pkin(2) * t340 - pkin(9) * t339 + t327 * t375 - t287;
t361 = (-t286 - t387) * qJ(4) + t251 + (-pkin(3) * t335 - 0.2e1 * qJD(4)) * t315;
t223 = -pkin(4) * t313 - t301 * t315 + (pkin(3) + qJ(5)) * t285 + t361;
t213 = -0.2e1 * qJD(5) * t298 + t219 * t347 - t345 * t223;
t210 = (t297 * t315 - t260) * pkin(10) + (t297 * t298 + t286) * pkin(5) + t213;
t214 = 0.2e1 * qJD(5) * t297 + t219 * t345 + t223 * t347;
t211 = -pkin(5) * t296 + pkin(10) * t259 - t268 * t315 + t214;
t209 = t210 * t349 + t211 * t353;
t236 = Ifges(7,5) * t255 + Ifges(7,6) * t254 + Ifges(7,3) * t312;
t238 = Ifges(7,1) * t255 + Ifges(7,4) * t254 + Ifges(7,5) * t312;
t282 = qJDD(6) + t286;
t195 = -mrSges(7,1) * t216 + mrSges(7,3) * t209 + Ifges(7,4) * t233 + Ifges(7,2) * t232 + Ifges(7,6) * t282 - t236 * t255 + t238 * t312;
t208 = t210 * t353 - t211 * t349;
t237 = Ifges(7,4) * t255 + Ifges(7,2) * t254 + Ifges(7,6) * t312;
t196 = mrSges(7,2) * t216 - mrSges(7,3) * t208 + Ifges(7,1) * t233 + Ifges(7,4) * t232 + Ifges(7,5) * t282 + t236 * t254 - t237 * t312;
t244 = Ifges(6,5) * t298 + Ifges(6,6) * t297 + Ifges(6,3) * t315;
t246 = Ifges(6,1) * t298 + Ifges(6,4) * t297 + Ifges(6,5) * t315;
t240 = -mrSges(7,1) * t254 + mrSges(7,2) * t255;
t203 = m(7) * t208 + mrSges(7,1) * t282 - mrSges(7,3) * t233 - t240 * t255 + t242 * t312;
t204 = m(7) * t209 - mrSges(7,2) * t282 + mrSges(7,3) * t232 + t240 * t254 - t243 * t312;
t372 = -t203 * t349 + t204 * t353;
t178 = -mrSges(6,1) * t222 + mrSges(6,3) * t214 + Ifges(6,4) * t260 + Ifges(6,2) * t259 + Ifges(6,6) * t286 - pkin(5) * t368 + pkin(10) * t372 + t353 * t195 + t349 * t196 - t298 * t244 + t315 * t246;
t194 = t203 * t353 + t204 * t349;
t245 = Ifges(6,4) * t298 + Ifges(6,2) * t297 + Ifges(6,6) * t315;
t180 = mrSges(6,2) * t222 - mrSges(6,3) * t213 + Ifges(6,1) * t260 + Ifges(6,4) * t259 + Ifges(6,5) * t286 - pkin(10) * t194 - t195 * t349 + t196 * t353 + t244 * t297 - t245 * t315;
t261 = -mrSges(6,1) * t297 + mrSges(6,2) * t298;
t191 = m(6) * t213 + mrSges(6,1) * t286 - mrSges(6,3) * t260 - t261 * t298 + t266 * t315 + t194;
t192 = m(6) * t214 - mrSges(6,2) * t286 + mrSges(6,3) * t259 + t261 * t297 - t267 * t315 + t372;
t187 = t191 * t347 + t192 * t345;
t270 = -Ifges(5,5) * t335 - Ifges(5,6) * t315 + Ifges(5,3) * t314;
t364 = -mrSges(5,2) * t226 + mrSges(5,3) * t224 - Ifges(5,1) * t321 + Ifges(5,4) * t286 - Ifges(5,5) * t285 + qJ(5) * t187 + t345 * t178 - t180 * t347 + t270 * t315;
t367 = -m(5) * t226 - mrSges(5,1) * t286 - t291 * t315 - t187;
t272 = -Ifges(5,4) * t335 - Ifges(5,2) * t315 + Ifges(5,6) * t314;
t381 = -Ifges(4,1) * t315 + Ifges(4,4) * t314 + Ifges(4,5) * t335 + t272;
t392 = -t314 * t381 + mrSges(4,1) * t234 - mrSges(4,2) * t235 + Ifges(4,5) * t286 - Ifges(4,6) * t285 + Ifges(4,3) * t321 + pkin(3) * (-mrSges(5,2) * t321 + t302 * t335 + t367) + qJ(4) * (-mrSges(5,1) * t285 - t291 * t314 + t360) + t315 * t273 - t364;
t388 = -Ifges(5,6) - Ifges(4,4);
t386 = t346 * t351;
t188 = -t191 * t345 + t192 * t347;
t290 = mrSges(4,1) * t314 + mrSges(4,2) * t315;
t299 = mrSges(4,2) * t335 - mrSges(4,3) * t314;
t184 = m(4) * t234 - mrSges(4,3) * t286 - t290 * t315 + (-t299 + t302) * t335 + (mrSges(4,1) - mrSges(5,2)) * t321 + t367;
t300 = -mrSges(4,1) * t335 - mrSges(4,3) * t315;
t199 = (-t290 - t291) * t314 + (-mrSges(4,3) - mrSges(5,1)) * t285 + m(4) * t235 - mrSges(4,2) * t321 + t300 * t335 + t360;
t177 = t184 * t391 + t199 * t350;
t274 = -Ifges(5,1) * t335 - Ifges(5,4) * t315 + Ifges(5,5) * t314;
t382 = -Ifges(4,5) * t315 + Ifges(4,6) * t314 + Ifges(4,3) * t335 - t274;
t288 = -g(3) * t386 + t380;
t322 = mrSges(3,1) * t341 - mrSges(3,3) * t375;
t326 = (-mrSges(3,1) * t354 + mrSges(3,2) * t351) * t379;
t373 = -t184 * t350 + t199 * t391;
t175 = m(3) * t288 - mrSges(3,2) * t340 - mrSges(3,3) * t329 - t322 * t341 + t326 * t374 + t373;
t323 = -mrSges(3,2) * t341 + mrSges(3,3) * t374;
t227 = pkin(3) * t285 + t361;
t371 = -m(5) * t227 + mrSges(5,2) * t285 + t302 * t314 - t188;
t362 = -m(4) * t251 - t285 * mrSges(4,1) - t314 * t299 + (-t300 + t303) * t315 + (-mrSges(4,2) + mrSges(5,3)) * t286 + t371;
t182 = m(3) * t287 + mrSges(3,1) * t340 - mrSges(3,3) * t328 + t323 * t341 - t326 * t375 + t362;
t172 = t175 * t354 - t182 * t351;
t308 = -t346 * t324 - t389;
t176 = m(3) * t308 + mrSges(3,1) * t329 + mrSges(3,2) * t328 + (t322 * t351 - t323 * t354) * t379 + t177;
t168 = t175 * t384 - t176 * t346 + t182 * t383;
t185 = -mrSges(5,3) * t286 - t303 * t315 - t371;
t363 = -mrSges(5,1) * t224 + mrSges(5,2) * t227 - pkin(4) * t206 - qJ(5) * t188 - t347 * t178 - t345 * t180;
t164 = -mrSges(4,1) * t251 + mrSges(4,3) * t235 - pkin(3) * t185 + t381 * t335 + (Ifges(4,6) - Ifges(5,5)) * t321 + t382 * t315 - t388 * t286 + (-Ifges(4,2) - Ifges(5,3)) * t285 + t363;
t365 = -mrSges(7,1) * t208 + mrSges(7,2) * t209 - Ifges(7,5) * t233 - Ifges(7,6) * t232 - Ifges(7,3) * t282 - t237 * t255 + t254 * t238;
t359 = -mrSges(6,1) * t213 + mrSges(6,2) * t214 - Ifges(6,5) * t260 - Ifges(6,6) * t259 - Ifges(6,3) * t286 - pkin(5) * t194 - t245 * t298 + t297 * t246 + t365;
t358 = -mrSges(5,1) * t226 + mrSges(5,3) * t227 - pkin(4) * t187 + t359;
t169 = -qJ(4) * t185 + (-t270 + t273) * t335 + (Ifges(4,5) - Ifges(5,4)) * t321 + t382 * t314 + (Ifges(4,1) + Ifges(5,2)) * t286 + t388 * t285 - t358 - mrSges(4,3) * t234 + mrSges(4,2) * t251;
t306 = Ifges(3,6) * t341 + (Ifges(3,4) * t351 + Ifges(3,2) * t354) * t379;
t307 = Ifges(3,5) * t341 + (Ifges(3,1) * t351 + Ifges(3,4) * t354) * t379;
t159 = Ifges(3,5) * t328 - Ifges(3,6) * t329 + Ifges(3,3) * t340 + mrSges(3,1) * t287 - mrSges(3,2) * t288 + t350 * t169 + t391 * t164 + pkin(2) * t362 + pkin(9) * t373 + (t306 * t351 - t307 * t354) * t379;
t305 = Ifges(3,3) * t341 + (Ifges(3,5) * t351 + Ifges(3,6) * t354) * t379;
t161 = mrSges(3,2) * t308 - mrSges(3,3) * t287 + Ifges(3,1) * t328 - Ifges(3,4) * t329 + Ifges(3,5) * t340 - pkin(9) * t177 - t164 * t350 + t169 * t391 + t305 * t374 - t306 * t341;
t163 = -mrSges(3,1) * t308 + mrSges(3,3) * t288 + Ifges(3,4) * t328 - Ifges(3,2) * t329 + Ifges(3,6) * t340 - pkin(2) * t177 - t305 * t375 + t341 * t307 - t392;
t366 = mrSges(2,1) * t337 - mrSges(2,2) * t338 + Ifges(2,3) * qJDD(1) + pkin(1) * t168 + t159 * t348 + t161 * t386 + t163 * t385 + t172 * t390;
t170 = m(2) * t338 - mrSges(2,1) * t356 - qJDD(1) * mrSges(2,2) + t172;
t167 = t348 * t176 + (t175 * t351 + t182 * t354) * t346;
t165 = m(2) * t337 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t356 + t168;
t157 = -mrSges(2,2) * g(3) - mrSges(2,3) * t337 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t356 + t161 * t354 - t163 * t351 + (-t167 * t346 - t168 * t348) * pkin(8);
t156 = mrSges(2,1) * g(3) + mrSges(2,3) * t338 + Ifges(2,5) * t356 + Ifges(2,6) * qJDD(1) - pkin(1) * t167 - t159 * t346 + (pkin(8) * t172 + t161 * t351 + t163 * t354) * t348;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t355 * t157 - t352 * t156 - pkin(7) * (t165 * t355 + t170 * t352), t157, t161, t169, -t314 * t272 - t364, t180, t196; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t352 * t157 + t355 * t156 + pkin(7) * (-t165 * t352 + t170 * t355), t156, t163, t164, Ifges(5,4) * t321 - Ifges(5,2) * t286 + Ifges(5,6) * t285 + t335 * t270 + t314 * t274 + t358, t178, t195; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t366, t366, t159, t392, Ifges(5,5) * t321 - Ifges(5,6) * t286 + Ifges(5,3) * t285 - t335 * t272 + t315 * t274 - t363, -t359, -t365;];
m_new  = t1;
