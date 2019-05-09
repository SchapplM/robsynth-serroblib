% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-05-07 09:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPRP11_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP11_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP11_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP11_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP11_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP11_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP11_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 09:07:47
% EndTime: 2019-05-07 09:08:23
% DurationCPUTime: 14.86s
% Computational Cost: add. (247345->399), mult. (530251->479), div. (0->0), fcn. (403672->10), ass. (0->148)
t347 = cos(pkin(6));
t342 = qJD(1) * t347 + qJD(2);
t349 = sin(qJ(3));
t350 = sin(qJ(2));
t346 = sin(pkin(6));
t383 = qJD(1) * t346;
t376 = t350 * t383;
t398 = cos(qJ(3));
t316 = t349 * t342 + t398 * t376;
t353 = cos(qJ(2));
t380 = qJD(1) * qJD(2);
t329 = (qJDD(1) * t350 + t353 * t380) * t346;
t341 = qJDD(1) * t347 + qJDD(2);
t288 = t316 * qJD(3) + t349 * t329 - t398 * t341;
t315 = -t398 * t342 + t349 * t376;
t382 = qJD(1) * t353;
t375 = t346 * t382;
t336 = -qJD(3) + t375;
t348 = sin(qJ(5));
t352 = cos(qJ(5));
t298 = t315 * t352 + t336 * t348;
t330 = (-qJDD(1) * t353 + t350 * t380) * t346;
t322 = qJDD(3) + t330;
t241 = qJD(5) * t298 + t288 * t348 + t322 * t352;
t299 = t315 * t348 - t336 * t352;
t258 = -mrSges(7,1) * t298 + mrSges(7,2) * t299;
t328 = (-pkin(2) * t353 - pkin(9) * t350) * t383;
t340 = t342 ^ 2;
t351 = sin(qJ(1));
t354 = cos(qJ(1));
t338 = t351 * g(1) - g(2) * t354;
t355 = qJD(1) ^ 2;
t397 = pkin(8) * t346;
t325 = qJDD(1) * pkin(1) + t355 * t397 + t338;
t339 = -g(1) * t354 - g(2) * t351;
t326 = -pkin(1) * t355 + qJDD(1) * t397 + t339;
t390 = t347 * t350;
t384 = t325 * t390 + t353 * t326;
t255 = -t340 * pkin(2) + t341 * pkin(9) + (-g(3) * t350 + t328 * t382) * t346 + t384;
t396 = t347 * g(3);
t256 = t330 * pkin(2) - t329 * pkin(9) - t396 + (-t325 + (pkin(2) * t350 - pkin(9) * t353) * t342 * qJD(1)) * t346;
t223 = -t349 * t255 + t398 * t256;
t292 = pkin(3) * t315 - qJ(4) * t316;
t335 = t336 ^ 2;
t221 = -t322 * pkin(3) - t335 * qJ(4) + t316 * t292 + qJDD(4) - t223;
t289 = -t315 * qJD(3) + t398 * t329 + t349 * t341;
t393 = t315 * t336;
t214 = (t315 * t316 - t322) * pkin(10) + (t289 - t393) * pkin(4) + t221;
t304 = pkin(4) * t316 + pkin(10) * t336;
t314 = t315 ^ 2;
t389 = t347 * t353;
t391 = t346 * t353;
t290 = -g(3) * t391 + t325 * t389 - t350 * t326;
t254 = -t341 * pkin(2) - t340 * pkin(9) + t328 * t376 - t290;
t358 = (-t289 - t393) * qJ(4) + t254 + (-pkin(3) * t336 - (2 * qJD(4))) * t316;
t218 = -t314 * pkin(4) - t316 * t304 + (pkin(3) + pkin(10)) * t288 + t358;
t207 = t352 * t214 - t348 * t218;
t285 = qJDD(5) + t289;
t313 = qJD(5) + t316;
t202 = -0.2e1 * qJD(6) * t299 + (t298 * t313 - t241) * qJ(6) + (t298 * t299 + t285) * pkin(5) + t207;
t264 = -mrSges(7,2) * t313 + mrSges(7,3) * t298;
t379 = m(7) * t202 + t285 * mrSges(7,1) + t313 * t264;
t199 = -t241 * mrSges(7,3) - t299 * t258 + t379;
t208 = t348 * t214 + t352 * t218;
t240 = -qJD(5) * t299 + t288 * t352 - t322 * t348;
t246 = Ifges(6,4) * t299 + Ifges(6,2) * t298 + Ifges(6,6) * t313;
t247 = Ifges(7,1) * t299 + Ifges(7,4) * t298 + Ifges(7,5) * t313;
t248 = Ifges(6,1) * t299 + Ifges(6,4) * t298 + Ifges(6,5) * t313;
t266 = pkin(5) * t313 - qJ(6) * t299;
t297 = t298 ^ 2;
t205 = -pkin(5) * t297 + qJ(6) * t240 + 0.2e1 * qJD(6) * t298 - t266 * t313 + t208;
t245 = Ifges(7,4) * t299 + Ifges(7,2) * t298 + Ifges(7,6) * t313;
t367 = -mrSges(7,1) * t202 + mrSges(7,2) * t205 - Ifges(7,5) * t241 - Ifges(7,6) * t240 - Ifges(7,3) * t285 - t299 * t245;
t402 = mrSges(6,1) * t207 - mrSges(6,2) * t208 + Ifges(6,5) * t241 + Ifges(6,6) * t240 + Ifges(6,3) * t285 + pkin(5) * t199 + t299 * t246 - t367 - (t248 + t247) * t298;
t259 = -mrSges(6,1) * t298 + mrSges(6,2) * t299;
t265 = -mrSges(6,2) * t313 + mrSges(6,3) * t298;
t193 = m(6) * t207 + t285 * mrSges(6,1) + t313 * t265 + (-t258 - t259) * t299 + (-mrSges(6,3) - mrSges(7,3)) * t241 + t379;
t267 = mrSges(7,1) * t313 - mrSges(7,3) * t299;
t268 = mrSges(6,1) * t313 - mrSges(6,3) * t299;
t378 = m(7) * t205 + t240 * mrSges(7,3) + t298 * t258;
t195 = m(6) * t208 + t240 * mrSges(6,3) + t298 * t259 + (-t267 - t268) * t313 + (-mrSges(6,2) - mrSges(7,2)) * t285 + t378;
t187 = t352 * t193 + t348 * t195;
t222 = t288 * pkin(3) + t358;
t401 = mrSges(5,1) * t221 - mrSges(5,3) * t222 + pkin(4) * t187 + t402;
t224 = t398 * t255 + t349 * t256;
t219 = pkin(3) * t335 - t322 * qJ(4) + 0.2e1 * qJD(4) * t336 + t315 * t292 - t224;
t277 = Ifges(4,4) * t316 - Ifges(4,2) * t315 - Ifges(4,6) * t336;
t294 = -mrSges(5,2) * t315 - mrSges(5,3) * t316;
t302 = mrSges(5,1) * t315 + mrSges(5,3) * t336;
t303 = mrSges(5,1) * t316 - mrSges(5,2) * t336;
t217 = -pkin(4) * t288 - pkin(10) * t314 - t336 * t304 - t219;
t211 = -pkin(5) * t240 - qJ(6) * t297 + t266 * t299 + qJDD(6) + t217;
t377 = m(7) * t211 + t241 * mrSges(7,2) + t299 * t267;
t399 = -m(6) * t217 - t241 * mrSges(6,2) + (mrSges(6,1) + mrSges(7,1)) * t240 - t299 * t268 + (t264 + t265) * t298 - t377;
t360 = -m(5) * t219 + t322 * mrSges(5,3) - t336 * t303 - t399;
t243 = Ifges(7,5) * t299 + Ifges(7,6) * t298 + Ifges(7,3) * t313;
t244 = Ifges(6,5) * t299 + Ifges(6,6) * t298 + Ifges(6,3) * t313;
t368 = -mrSges(7,1) * t211 + mrSges(7,3) * t205 + Ifges(7,4) * t241 + Ifges(7,2) * t240 + Ifges(7,6) * t285 + t313 * t247;
t178 = Ifges(6,4) * t241 + Ifges(6,2) * t240 + Ifges(6,6) * t285 + t313 * t248 - mrSges(6,1) * t217 + mrSges(6,3) * t208 - pkin(5) * (-t240 * mrSges(7,1) - t298 * t264 + t377) + qJ(6) * (-t285 * mrSges(7,2) - t313 * t267 + t378) + (-t244 - t243) * t299 + t368;
t366 = mrSges(7,2) * t211 - mrSges(7,3) * t202 + Ifges(7,1) * t241 + Ifges(7,4) * t240 + Ifges(7,5) * t285 + t298 * t243;
t186 = mrSges(6,2) * t217 - mrSges(6,3) * t207 + Ifges(6,1) * t241 + Ifges(6,4) * t240 + Ifges(6,5) * t285 - qJ(6) * t199 + t298 * t244 + (-t245 - t246) * t313 + t366;
t274 = -Ifges(5,5) * t336 - Ifges(5,6) * t316 + Ifges(5,3) * t315;
t363 = -mrSges(5,2) * t221 + mrSges(5,3) * t219 - Ifges(5,1) * t322 + Ifges(5,4) * t289 - Ifges(5,5) * t288 + pkin(10) * t187 + t348 * t178 - t352 * t186 + t316 * t274;
t365 = -m(5) * t221 - t289 * mrSges(5,1) - t316 * t294 - t187;
t276 = -Ifges(5,4) * t336 - Ifges(5,2) * t316 + Ifges(5,6) * t315;
t385 = -Ifges(4,1) * t316 + Ifges(4,4) * t315 + Ifges(4,5) * t336 + t276;
t400 = -t385 * t315 + mrSges(4,1) * t223 - mrSges(4,2) * t224 + Ifges(4,5) * t289 - Ifges(4,6) * t288 + Ifges(4,3) * t322 + pkin(3) * (-t322 * mrSges(5,2) + t336 * t302 + t365) + qJ(4) * (-t288 * mrSges(5,1) - t315 * t294 + t360) + t316 * t277 - t363;
t394 = -Ifges(5,6) - Ifges(4,4);
t392 = t346 * t350;
t293 = mrSges(4,1) * t315 + mrSges(4,2) * t316;
t300 = mrSges(4,2) * t336 - mrSges(4,3) * t315;
t182 = m(4) * t223 - t289 * mrSges(4,3) - t316 * t293 + (-t300 + t302) * t336 + (mrSges(4,1) - mrSges(5,2)) * t322 + t365;
t301 = -mrSges(4,1) * t336 - mrSges(4,3) * t316;
t191 = (-t293 - t294) * t315 + (-mrSges(4,3) - mrSges(5,1)) * t288 + t360 + t336 * t301 - t322 * mrSges(4,2) + m(4) * t224;
t177 = t398 * t182 + t349 * t191;
t188 = -t348 * t193 + t352 * t195;
t278 = -Ifges(5,1) * t336 - Ifges(5,4) * t316 + Ifges(5,5) * t315;
t386 = -Ifges(4,5) * t316 + Ifges(4,6) * t315 + Ifges(4,3) * t336 - t278;
t291 = -g(3) * t392 + t384;
t323 = mrSges(3,1) * t342 - mrSges(3,3) * t376;
t327 = (-mrSges(3,1) * t353 + mrSges(3,2) * t350) * t383;
t373 = -t182 * t349 + t398 * t191;
t175 = m(3) * t291 - mrSges(3,2) * t341 - mrSges(3,3) * t330 - t323 * t342 + t327 * t375 + t373;
t324 = -mrSges(3,2) * t342 + mrSges(3,3) * t375;
t371 = -m(5) * t222 + t288 * mrSges(5,2) + t315 * t302 - t188;
t359 = -m(4) * t254 - t288 * mrSges(4,1) - t315 * t300 + (-t301 + t303) * t316 + (-mrSges(4,2) + mrSges(5,3)) * t289 + t371;
t180 = m(3) * t290 + t341 * mrSges(3,1) - t329 * mrSges(3,3) + t342 * t324 - t327 * t376 + t359;
t172 = t353 * t175 - t180 * t350;
t309 = -t346 * t325 - t396;
t176 = m(3) * t309 + t330 * mrSges(3,1) + t329 * mrSges(3,2) + (t323 * t350 - t324 * t353) * t383 + t177;
t168 = t175 * t390 - t176 * t346 + t180 * t389;
t183 = -t289 * mrSges(5,3) - t316 * t303 - t371;
t362 = -mrSges(5,1) * t219 + mrSges(5,2) * t222 - pkin(4) * t399 - pkin(10) * t188 - t352 * t178 - t348 * t186;
t164 = -mrSges(4,1) * t254 + mrSges(4,3) * t224 - pkin(3) * t183 + t385 * t336 + (Ifges(4,6) - Ifges(5,5)) * t322 + t386 * t316 - t394 * t289 + (-Ifges(4,2) - Ifges(5,3)) * t288 + t362;
t169 = (-t274 + t277) * t336 + (Ifges(4,5) - Ifges(5,4)) * t322 + t386 * t315 + (Ifges(4,1) + Ifges(5,2)) * t289 + t394 * t288 + mrSges(4,2) * t254 - mrSges(4,3) * t223 - qJ(4) * t183 + t401;
t307 = Ifges(3,6) * t342 + (Ifges(3,4) * t350 + Ifges(3,2) * t353) * t383;
t308 = Ifges(3,5) * t342 + (Ifges(3,1) * t350 + Ifges(3,4) * t353) * t383;
t159 = Ifges(3,5) * t329 - Ifges(3,6) * t330 + Ifges(3,3) * t341 + mrSges(3,1) * t290 - mrSges(3,2) * t291 + t349 * t169 + t398 * t164 + pkin(2) * t359 + pkin(9) * t373 + (t307 * t350 - t308 * t353) * t383;
t306 = Ifges(3,3) * t342 + (Ifges(3,5) * t350 + Ifges(3,6) * t353) * t383;
t161 = mrSges(3,2) * t309 - mrSges(3,3) * t290 + Ifges(3,1) * t329 - Ifges(3,4) * t330 + Ifges(3,5) * t341 - pkin(9) * t177 - t349 * t164 + t398 * t169 + t306 * t375 - t342 * t307;
t163 = -mrSges(3,1) * t309 + mrSges(3,3) * t291 + Ifges(3,4) * t329 - Ifges(3,2) * t330 + Ifges(3,6) * t341 - pkin(2) * t177 - t306 * t376 + t342 * t308 - t400;
t364 = mrSges(2,1) * t338 - mrSges(2,2) * t339 + Ifges(2,3) * qJDD(1) + pkin(1) * t168 + t347 * t159 + t161 * t392 + t163 * t391 + t172 * t397;
t170 = m(2) * t339 - mrSges(2,1) * t355 - qJDD(1) * mrSges(2,2) + t172;
t167 = t347 * t176 + (t175 * t350 + t180 * t353) * t346;
t165 = m(2) * t338 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t355 + t168;
t157 = -mrSges(2,2) * g(3) - mrSges(2,3) * t338 + Ifges(2,5) * qJDD(1) - t355 * Ifges(2,6) + t353 * t161 - t350 * t163 + (-t167 * t346 - t168 * t347) * pkin(8);
t156 = mrSges(2,1) * g(3) + mrSges(2,3) * t339 + t355 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t167 - t346 * t159 + (pkin(8) * t172 + t161 * t350 + t163 * t353) * t347;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t354 * t157 - t351 * t156 - pkin(7) * (t165 * t354 + t170 * t351), t157, t161, t169, -t315 * t276 - t363, t186, -t245 * t313 + t366; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t351 * t157 + t354 * t156 + pkin(7) * (-t165 * t351 + t170 * t354), t156, t163, t164, Ifges(5,4) * t322 - Ifges(5,2) * t289 + Ifges(5,6) * t288 + t336 * t274 + t315 * t278 - t401, t178, -t299 * t243 + t368; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t364, t364, t159, t400, Ifges(5,5) * t322 - Ifges(5,6) * t289 + Ifges(5,3) * t288 - t336 * t276 + t316 * t278 - t362, t402, -t298 * t247 - t367;];
m_new  = t1;
