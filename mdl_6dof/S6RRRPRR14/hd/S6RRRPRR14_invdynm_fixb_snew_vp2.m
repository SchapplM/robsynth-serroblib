% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-05-07 16:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPRR14_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR14_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR14_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR14_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR14_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR14_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR14_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 16:33:25
% EndTime: 2019-05-07 16:34:44
% DurationCPUTime: 31.48s
% Computational Cost: add. (558603->403), mult. (1195721->497), div. (0->0), fcn. (933170->12), ass. (0->159)
t345 = sin(pkin(6));
t350 = sin(qJ(2));
t354 = cos(qJ(2));
t376 = qJD(1) * qJD(2);
t329 = (-qJDD(1) * t354 + t350 * t376) * t345;
t379 = qJD(1) * t345;
t327 = (-pkin(2) * t354 - pkin(9) * t350) * t379;
t346 = cos(pkin(6));
t341 = t346 * qJD(1) + qJD(2);
t339 = t341 ^ 2;
t340 = t346 * qJDD(1) + qJDD(2);
t378 = qJD(1) * t354;
t351 = sin(qJ(1));
t355 = cos(qJ(1));
t337 = t351 * g(1) - t355 * g(2);
t356 = qJD(1) ^ 2;
t390 = pkin(8) * t345;
t324 = qJDD(1) * pkin(1) + t356 * t390 + t337;
t338 = -t355 * g(1) - t351 * g(2);
t325 = -t356 * pkin(1) + qJDD(1) * t390 + t338;
t384 = t346 * t350;
t380 = t324 * t384 + t354 * t325;
t257 = -t339 * pkin(2) + t340 * pkin(9) + (-g(3) * t350 + t327 * t378) * t345 + t380;
t328 = (qJDD(1) * t350 + t354 * t376) * t345;
t389 = t346 * g(3);
t258 = t329 * pkin(2) - t328 * pkin(9) - t389 + (-t324 + (pkin(2) * t350 - pkin(9) * t354) * t341 * qJD(1)) * t345;
t349 = sin(qJ(3));
t391 = cos(qJ(3));
t235 = t391 * t257 + t349 * t258;
t375 = t350 * t379;
t314 = -t391 * t341 + t349 * t375;
t315 = t349 * t341 + t391 * t375;
t290 = t314 * pkin(3) - t315 * qJ(4);
t321 = qJDD(3) + t329;
t374 = t345 * t378;
t335 = -qJD(3) + t374;
t334 = t335 ^ 2;
t229 = pkin(3) * t334 - t321 * qJ(4) + 0.2e1 * qJD(4) * t335 + t314 * t290 - t235;
t234 = -t349 * t257 + t391 * t258;
t275 = Ifges(4,4) * t315 - Ifges(4,2) * t314 - Ifges(4,6) * t335;
t286 = t315 * qJD(3) + t349 * t328 - t391 * t340;
t287 = -t314 * qJD(3) + t391 * t328 + t349 * t340;
t292 = -t314 * mrSges(5,2) - t315 * mrSges(5,3);
t300 = t314 * mrSges(5,1) + t335 * mrSges(5,3);
t302 = t315 * pkin(4) + t335 * pkin(10);
t313 = t314 ^ 2;
t222 = -pkin(4) * t286 - pkin(10) * t313 - t335 * t302 - t229;
t348 = sin(qJ(5));
t353 = cos(qJ(5));
t297 = t348 * t314 - t353 * t335;
t245 = -t297 * qJD(5) + t353 * t286 - t348 * t321;
t296 = t353 * t314 + t348 * t335;
t246 = t296 * qJD(5) + t348 * t286 + t353 * t321;
t312 = qJD(5) + t315;
t266 = -t312 * mrSges(6,2) + t296 * mrSges(6,3);
t267 = t312 * mrSges(6,1) - t297 * mrSges(6,3);
t268 = t312 * pkin(5) - t297 * pkin(11);
t295 = t296 ^ 2;
t216 = -pkin(5) * t245 - pkin(11) * t295 + t268 * t297 + t222;
t347 = sin(qJ(6));
t352 = cos(qJ(6));
t260 = t347 * t296 + t352 * t297;
t227 = -qJD(6) * t260 + t245 * t352 - t246 * t347;
t259 = t352 * t296 - t347 * t297;
t228 = qJD(6) * t259 + t245 * t347 + t246 * t352;
t310 = qJD(6) + t312;
t247 = -t310 * mrSges(7,2) + t259 * mrSges(7,3);
t248 = t310 * mrSges(7,1) - t260 * mrSges(7,3);
t368 = m(7) * t216 - t227 * mrSges(7,1) + t228 * mrSges(7,2) - t259 * t247 + t260 * t248;
t206 = -m(6) * t222 + t245 * mrSges(6,1) - t246 * mrSges(6,2) + t296 * t266 - t297 * t267 - t368;
t301 = t315 * mrSges(5,1) - t335 * mrSges(5,2);
t360 = -m(5) * t229 + t321 * mrSges(5,3) - t335 * t301 - t206;
t231 = -t321 * pkin(3) - t334 * qJ(4) + t315 * t290 + qJDD(4) - t234;
t387 = t314 * t335;
t219 = (t314 * t315 - t321) * pkin(10) + (t287 - t387) * pkin(4) + t231;
t383 = t346 * t354;
t385 = t345 * t354;
t288 = -g(3) * t385 + t324 * t383 - t350 * t325;
t256 = -t340 * pkin(2) - t339 * pkin(9) + t327 * t375 - t288;
t361 = (-t287 - t387) * qJ(4) + t256 + (-pkin(3) * t335 - 0.2e1 * qJD(4)) * t315;
t223 = -t313 * pkin(4) - t315 * t302 + (pkin(3) + pkin(10)) * t286 + t361;
t213 = t353 * t219 - t348 * t223;
t283 = qJDD(5) + t287;
t210 = (t296 * t312 - t246) * pkin(11) + (t296 * t297 + t283) * pkin(5) + t213;
t214 = t348 * t219 + t353 * t223;
t211 = -pkin(5) * t295 + pkin(11) * t245 - t268 * t312 + t214;
t209 = t210 * t347 + t211 * t352;
t236 = Ifges(7,5) * t260 + Ifges(7,6) * t259 + Ifges(7,3) * t310;
t238 = Ifges(7,1) * t260 + Ifges(7,4) * t259 + Ifges(7,5) * t310;
t271 = qJDD(6) + t283;
t195 = -mrSges(7,1) * t216 + mrSges(7,3) * t209 + Ifges(7,4) * t228 + Ifges(7,2) * t227 + Ifges(7,6) * t271 - t236 * t260 + t238 * t310;
t208 = t210 * t352 - t211 * t347;
t237 = Ifges(7,4) * t260 + Ifges(7,2) * t259 + Ifges(7,6) * t310;
t196 = mrSges(7,2) * t216 - mrSges(7,3) * t208 + Ifges(7,1) * t228 + Ifges(7,4) * t227 + Ifges(7,5) * t271 + t236 * t259 - t237 * t310;
t249 = Ifges(6,5) * t297 + Ifges(6,6) * t296 + Ifges(6,3) * t312;
t251 = Ifges(6,1) * t297 + Ifges(6,4) * t296 + Ifges(6,5) * t312;
t240 = -t259 * mrSges(7,1) + t260 * mrSges(7,2);
t204 = m(7) * t208 + mrSges(7,1) * t271 - t228 * mrSges(7,3) - t240 * t260 + t247 * t310;
t205 = m(7) * t209 - mrSges(7,2) * t271 + t227 * mrSges(7,3) + t240 * t259 - t248 * t310;
t372 = -t204 * t347 + t352 * t205;
t178 = -mrSges(6,1) * t222 + mrSges(6,3) * t214 + Ifges(6,4) * t246 + Ifges(6,2) * t245 + Ifges(6,6) * t283 - pkin(5) * t368 + pkin(11) * t372 + t352 * t195 + t347 * t196 - t297 * t249 + t312 * t251;
t194 = t352 * t204 + t347 * t205;
t250 = Ifges(6,4) * t297 + Ifges(6,2) * t296 + Ifges(6,6) * t312;
t180 = mrSges(6,2) * t222 - mrSges(6,3) * t213 + Ifges(6,1) * t246 + Ifges(6,4) * t245 + Ifges(6,5) * t283 - pkin(11) * t194 - t195 * t347 + t196 * t352 + t249 * t296 - t250 * t312;
t261 = -t296 * mrSges(6,1) + t297 * mrSges(6,2);
t191 = m(6) * t213 + mrSges(6,1) * t283 - mrSges(6,3) * t246 - t261 * t297 + t266 * t312 + t194;
t192 = m(6) * t214 - mrSges(6,2) * t283 + mrSges(6,3) * t245 + t261 * t296 - t267 * t312 + t372;
t187 = t353 * t191 + t348 * t192;
t272 = -Ifges(5,5) * t335 - Ifges(5,6) * t315 + Ifges(5,3) * t314;
t364 = -mrSges(5,2) * t231 + mrSges(5,3) * t229 - Ifges(5,1) * t321 + Ifges(5,4) * t287 - Ifges(5,5) * t286 + pkin(10) * t187 + t348 * t178 - t353 * t180 + t315 * t272;
t367 = -m(5) * t231 - t287 * mrSges(5,1) - t315 * t292 - t187;
t274 = -Ifges(5,4) * t335 - Ifges(5,2) * t315 + Ifges(5,6) * t314;
t381 = -Ifges(4,1) * t315 + Ifges(4,4) * t314 + Ifges(4,5) * t335 + t274;
t392 = -t381 * t314 + mrSges(4,1) * t234 - mrSges(4,2) * t235 + Ifges(4,5) * t287 - Ifges(4,6) * t286 + Ifges(4,3) * t321 + pkin(3) * (-t321 * mrSges(5,2) + t335 * t300 + t367) + qJ(4) * (-t286 * mrSges(5,1) - t314 * t292 + t360) + t315 * t275 - t364;
t388 = -Ifges(5,6) - Ifges(4,4);
t386 = t345 * t350;
t188 = -t348 * t191 + t353 * t192;
t291 = t314 * mrSges(4,1) + t315 * mrSges(4,2);
t298 = t335 * mrSges(4,2) - t314 * mrSges(4,3);
t184 = m(4) * t234 - t287 * mrSges(4,3) - t315 * t291 + (-t298 + t300) * t335 + (mrSges(4,1) - mrSges(5,2)) * t321 + t367;
t299 = -t335 * mrSges(4,1) - t315 * mrSges(4,3);
t199 = m(4) * t235 + (-t291 - t292) * t314 + (-mrSges(4,3) - mrSges(5,1)) * t286 + t360 - t321 * mrSges(4,2) + t335 * t299;
t177 = t391 * t184 + t349 * t199;
t276 = -Ifges(5,1) * t335 - Ifges(5,4) * t315 + Ifges(5,5) * t314;
t382 = -Ifges(4,5) * t315 + Ifges(4,6) * t314 + Ifges(4,3) * t335 - t276;
t289 = -g(3) * t386 + t380;
t322 = t341 * mrSges(3,1) - mrSges(3,3) * t375;
t326 = (-mrSges(3,1) * t354 + mrSges(3,2) * t350) * t379;
t373 = -t184 * t349 + t391 * t199;
t175 = m(3) * t289 - mrSges(3,2) * t340 - mrSges(3,3) * t329 - t322 * t341 + t326 * t374 + t373;
t323 = -t341 * mrSges(3,2) + mrSges(3,3) * t374;
t232 = t286 * pkin(3) + t361;
t371 = -m(5) * t232 + t286 * mrSges(5,2) + t314 * t300 - t188;
t362 = -m(4) * t256 - t286 * mrSges(4,1) - t314 * t298 + (-t299 + t301) * t315 + (-mrSges(4,2) + mrSges(5,3)) * t287 + t371;
t182 = m(3) * t288 + t340 * mrSges(3,1) - t328 * mrSges(3,3) + t341 * t323 - t326 * t375 + t362;
t172 = t354 * t175 - t182 * t350;
t307 = -t345 * t324 - t389;
t176 = m(3) * t307 + t329 * mrSges(3,1) + t328 * mrSges(3,2) + (t322 * t350 - t323 * t354) * t379 + t177;
t168 = t175 * t384 - t176 * t345 + t182 * t383;
t185 = -t287 * mrSges(5,3) - t315 * t301 - t371;
t363 = -mrSges(5,1) * t229 + mrSges(5,2) * t232 - pkin(4) * t206 - pkin(10) * t188 - t353 * t178 - t348 * t180;
t164 = -mrSges(4,1) * t256 + mrSges(4,3) * t235 - pkin(3) * t185 + t381 * t335 + (Ifges(4,6) - Ifges(5,5)) * t321 + t382 * t315 - t388 * t287 + (-Ifges(4,2) - Ifges(5,3)) * t286 + t363;
t365 = -mrSges(7,1) * t208 + mrSges(7,2) * t209 - Ifges(7,5) * t228 - Ifges(7,6) * t227 - Ifges(7,3) * t271 - t260 * t237 + t259 * t238;
t359 = -mrSges(6,1) * t213 + mrSges(6,2) * t214 - Ifges(6,5) * t246 - Ifges(6,6) * t245 - Ifges(6,3) * t283 - pkin(5) * t194 - t297 * t250 + t296 * t251 + t365;
t358 = -mrSges(5,1) * t231 + mrSges(5,3) * t232 - pkin(4) * t187 + t359;
t169 = -qJ(4) * t185 - mrSges(4,3) * t234 + (-t272 + t275) * t335 + (Ifges(4,5) - Ifges(5,4)) * t321 + t382 * t314 + (Ifges(4,1) + Ifges(5,2)) * t287 + t388 * t286 - t358 + mrSges(4,2) * t256;
t305 = Ifges(3,6) * t341 + (Ifges(3,4) * t350 + Ifges(3,2) * t354) * t379;
t306 = Ifges(3,5) * t341 + (Ifges(3,1) * t350 + Ifges(3,4) * t354) * t379;
t159 = Ifges(3,5) * t328 - Ifges(3,6) * t329 + Ifges(3,3) * t340 + mrSges(3,1) * t288 - mrSges(3,2) * t289 + t349 * t169 + t391 * t164 + pkin(2) * t362 + pkin(9) * t373 + (t305 * t350 - t306 * t354) * t379;
t304 = Ifges(3,3) * t341 + (Ifges(3,5) * t350 + Ifges(3,6) * t354) * t379;
t161 = mrSges(3,2) * t307 - mrSges(3,3) * t288 + Ifges(3,1) * t328 - Ifges(3,4) * t329 + Ifges(3,5) * t340 - pkin(9) * t177 - t349 * t164 + t391 * t169 + t304 * t374 - t341 * t305;
t163 = -mrSges(3,1) * t307 + mrSges(3,3) * t289 + Ifges(3,4) * t328 - Ifges(3,2) * t329 + Ifges(3,6) * t340 - pkin(2) * t177 - t304 * t375 + t341 * t306 - t392;
t366 = mrSges(2,1) * t337 - mrSges(2,2) * t338 + Ifges(2,3) * qJDD(1) + pkin(1) * t168 + t346 * t159 + t161 * t386 + t163 * t385 + t172 * t390;
t170 = m(2) * t338 - mrSges(2,1) * t356 - qJDD(1) * mrSges(2,2) + t172;
t167 = t346 * t176 + (t175 * t350 + t182 * t354) * t345;
t165 = m(2) * t337 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t356 + t168;
t157 = -mrSges(2,2) * g(3) - mrSges(2,3) * t337 + Ifges(2,5) * qJDD(1) - t356 * Ifges(2,6) + t354 * t161 - t350 * t163 + (-t167 * t345 - t168 * t346) * pkin(8);
t156 = mrSges(2,1) * g(3) + mrSges(2,3) * t338 + t356 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t167 - t345 * t159 + (pkin(8) * t172 + t161 * t350 + t163 * t354) * t346;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t355 * t157 - t351 * t156 - pkin(7) * (t165 * t355 + t170 * t351), t157, t161, t169, -t314 * t274 - t364, t180, t196; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t351 * t157 + t355 * t156 + pkin(7) * (-t165 * t351 + t170 * t355), t156, t163, t164, Ifges(5,4) * t321 - Ifges(5,2) * t287 + Ifges(5,6) * t286 + t335 * t272 + t314 * t276 + t358, t178, t195; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t366, t366, t159, t392, Ifges(5,5) * t321 - Ifges(5,6) * t287 + Ifges(5,3) * t286 - t335 * t274 + t315 * t276 - t363, -t359, -t365;];
m_new  = t1;
