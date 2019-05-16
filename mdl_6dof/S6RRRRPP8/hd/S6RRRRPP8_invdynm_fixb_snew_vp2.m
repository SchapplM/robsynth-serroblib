% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRPP8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-05-07 19:10
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRPP8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP8_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP8_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP8_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP8_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP8_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:57:28
% EndTime: 2019-05-07 18:58:02
% DurationCPUTime: 15.33s
% Computational Cost: add. (271799->396), mult. (576593->479), div. (0->0), fcn. (446166->10), ass. (0->151)
t346 = sin(pkin(6));
t350 = sin(qJ(2));
t352 = cos(qJ(2));
t374 = qJD(1) * qJD(2);
t332 = (-qJDD(1) * t352 + t350 * t374) * t346;
t376 = qJD(1) * t346;
t330 = (-pkin(2) * t352 - pkin(9) * t350) * t376;
t347 = cos(pkin(6));
t342 = qJD(1) * t347 + qJD(2);
t341 = t342 ^ 2;
t368 = qJDD(1) * t347 + qJDD(2);
t375 = qJD(1) * t352;
t351 = sin(qJ(1));
t353 = cos(qJ(1));
t339 = t351 * g(1) - g(2) * t353;
t354 = qJD(1) ^ 2;
t389 = pkin(8) * t346;
t327 = qJDD(1) * pkin(1) + t354 * t389 + t339;
t340 = -g(1) * t353 - g(2) * t351;
t328 = -pkin(1) * t354 + qJDD(1) * t389 + t340;
t382 = t347 * t350;
t377 = t327 * t382 + t352 * t328;
t269 = t368 * pkin(9) - t341 * pkin(2) + (-g(3) * t350 + t330 * t375) * t346 + t377;
t331 = (qJDD(1) * t350 + t352 * t374) * t346;
t388 = t347 * g(3);
t270 = t332 * pkin(2) - t331 * pkin(9) - t388 + (-t327 + (pkin(2) * t350 - pkin(9) * t352) * t342 * qJD(1)) * t346;
t349 = sin(qJ(3));
t392 = cos(qJ(3));
t222 = t392 * t269 + t349 * t270;
t372 = t350 * t376;
t320 = t392 * t342 - t349 * t372;
t321 = t349 * t342 + t372 * t392;
t304 = -pkin(3) * t320 - pkin(10) * t321;
t324 = qJDD(3) + t332;
t371 = t346 * t375;
t338 = qJD(3) - t371;
t337 = t338 ^ 2;
t218 = -pkin(3) * t337 + pkin(10) * t324 + t304 * t320 + t222;
t381 = t347 * t352;
t383 = t346 * t352;
t301 = -g(3) * t383 + t327 * t381 - t350 * t328;
t268 = -pkin(2) * t368 - t341 * pkin(9) + t330 * t372 - t301;
t299 = -qJD(3) * t321 - t331 * t349 + t392 * t368;
t300 = t320 * qJD(3) + t331 * t392 + t349 * t368;
t220 = (-t320 * t338 - t300) * pkin(10) + (t321 * t338 - t299) * pkin(3) + t268;
t348 = sin(qJ(4));
t391 = cos(qJ(4));
t214 = t391 * t218 + t348 * t220;
t306 = t348 * t321 - t338 * t391;
t307 = t321 * t391 + t348 * t338;
t273 = pkin(4) * t306 - qJ(5) * t307;
t297 = qJDD(4) - t299;
t318 = qJD(4) - t320;
t317 = t318 ^ 2;
t393 = 2 * qJD(5);
t209 = -pkin(4) * t317 + t297 * qJ(5) - t306 * t273 + t318 * t393 + t214;
t284 = -mrSges(6,1) * t318 + mrSges(6,2) * t307;
t398 = m(6) * t209 + t297 * mrSges(6,3) + t318 * t284;
t213 = -t348 * t218 + t220 * t391;
t211 = -t297 * pkin(4) - t317 * qJ(5) + t307 * t273 + qJDD(5) - t213;
t278 = -mrSges(6,2) * t306 + mrSges(6,3) * t318;
t397 = -m(6) * t211 + t297 * mrSges(6,1) + t318 * t278;
t246 = -t306 * qJD(4) + t300 * t391 + t348 * t324;
t221 = -t349 * t269 + t392 * t270;
t362 = t324 * pkin(3) + t337 * pkin(10) - t321 * t304 + t221;
t386 = t306 * t318;
t396 = (-t246 + t386) * qJ(5) - t362;
t279 = mrSges(7,2) * t318 + mrSges(7,3) * t306;
t394 = -0.2e1 * t307;
t199 = qJD(6) * t394 + (-t246 - t386) * qJ(6) + (t306 * t307 - t297) * pkin(5) + t211;
t275 = -mrSges(7,1) * t306 + mrSges(7,2) * t307;
t367 = -m(7) * t199 + t246 * mrSges(7,3) + t307 * t275;
t195 = -t297 * mrSges(7,1) - t318 * t279 - t367;
t282 = -mrSges(7,1) * t318 - mrSges(7,3) * t307;
t245 = t307 * qJD(4) + t348 * t300 - t324 * t391;
t281 = -pkin(5) * t318 - qJ(6) * t307;
t305 = t306 ^ 2;
t203 = -pkin(5) * t305 + qJ(6) * t245 + 0.2e1 * qJD(6) * t306 + t281 * t318 + t209;
t373 = m(7) * t203 + t245 * mrSges(7,3) + t306 * t275;
t196 = t297 * mrSges(7,2) + t318 * t282 + t373;
t252 = Ifges(6,5) * t307 + Ifges(6,6) * t318 + Ifges(6,3) * t306;
t256 = Ifges(5,4) * t307 - Ifges(5,2) * t306 + Ifges(5,6) * t318;
t259 = Ifges(5,1) * t307 - Ifges(5,4) * t306 + Ifges(5,5) * t318;
t274 = mrSges(6,1) * t306 - mrSges(6,3) * t307;
t258 = Ifges(6,1) * t307 + Ifges(6,4) * t318 + Ifges(6,5) * t306;
t254 = Ifges(7,4) * t307 + Ifges(7,2) * t306 - Ifges(7,6) * t318;
t257 = Ifges(7,1) * t307 + Ifges(7,4) * t306 - Ifges(7,5) * t318;
t365 = mrSges(7,1) * t199 - mrSges(7,2) * t203 + Ifges(7,5) * t246 + Ifges(7,6) * t245 - Ifges(7,3) * t297 + t307 * t254 - t306 * t257;
t359 = mrSges(6,1) * t211 - mrSges(6,3) * t209 - Ifges(6,4) * t246 - Ifges(6,2) * t297 - Ifges(6,6) * t245 + pkin(5) * t195 - t306 * t258 + t365;
t395 = (t256 - t252) * t307 + mrSges(5,1) * t213 - mrSges(5,2) * t214 + Ifges(5,5) * t246 - Ifges(5,6) * t245 + Ifges(5,3) * t297 + pkin(4) * (-t246 * mrSges(6,2) - t307 * t274 - t195 + t397) + qJ(5) * (-t245 * mrSges(6,2) - t306 * t274 + t196 + t398) + t306 * t259 - t359;
t387 = -mrSges(5,3) - mrSges(6,2);
t385 = t318 * t254;
t384 = t346 * t350;
t280 = -mrSges(5,2) * t318 - mrSges(5,3) * t306;
t378 = -mrSges(5,1) * t306 - mrSges(5,2) * t307 - t274;
t188 = m(5) * t213 + (t279 + t280) * t318 + t378 * t307 + (mrSges(5,1) + mrSges(7,1)) * t297 + t387 * t246 + t367 + t397;
t283 = mrSges(5,1) * t318 - mrSges(5,3) * t307;
t190 = m(5) * t214 + (t282 - t283) * t318 + t378 * t306 + (-mrSges(5,2) + mrSges(7,2)) * t297 + t387 * t245 + t373 + t398;
t185 = -t188 * t348 + t391 * t190;
t303 = -mrSges(4,1) * t320 + mrSges(4,2) * t321;
t309 = mrSges(4,1) * t338 - mrSges(4,3) * t321;
t183 = m(4) * t222 - mrSges(4,2) * t324 + mrSges(4,3) * t299 + t303 * t320 - t309 * t338 + t185;
t206 = -t305 * qJ(6) + qJDD(6) + (-pkin(4) - pkin(5)) * t245 + (-pkin(4) * t318 + t281 + t393) * t307 - t396;
t197 = -m(7) * t206 + t245 * mrSges(7,1) - t246 * mrSges(7,2) + t306 * t279 - t307 * t282;
t212 = qJD(5) * t394 + (t307 * t318 + t245) * pkin(4) + t396;
t194 = m(6) * t212 + mrSges(6,1) * t245 - t246 * mrSges(6,3) + t278 * t306 - t307 * t284 + t197;
t191 = m(5) * t362 - t245 * mrSges(5,1) - mrSges(5,2) * t246 - t306 * t280 - t283 * t307 - t194;
t308 = -mrSges(4,2) * t338 + mrSges(4,3) * t320;
t187 = m(4) * t221 + mrSges(4,1) * t324 - mrSges(4,3) * t300 - t303 * t321 + t308 * t338 + t191;
t177 = t349 * t183 + t392 * t187;
t255 = Ifges(6,4) * t307 + Ifges(6,2) * t318 + Ifges(6,6) * t306;
t380 = -Ifges(5,5) * t307 + Ifges(5,6) * t306 - Ifges(5,3) * t318 - t255;
t302 = -g(3) * t384 + t377;
t325 = mrSges(3,1) * t342 - mrSges(3,3) * t372;
t329 = (-mrSges(3,1) * t352 + mrSges(3,2) * t350) * t376;
t370 = t392 * t183 - t349 * t187;
t175 = m(3) * t302 - mrSges(3,2) * t368 - t332 * mrSges(3,3) - t342 * t325 + t329 * t371 + t370;
t326 = -mrSges(3,2) * t342 + mrSges(3,3) * t371;
t184 = t188 * t391 + t348 * t190;
t360 = -m(4) * t268 + t299 * mrSges(4,1) - t300 * mrSges(4,2) + t320 * t308 - t321 * t309 - t184;
t180 = m(3) * t301 + mrSges(3,1) * t368 - t331 * mrSges(3,3) + t342 * t326 - t329 * t372 + t360;
t171 = t352 * t175 - t180 * t350;
t313 = -t346 * t327 - t388;
t176 = m(3) * t313 + t332 * mrSges(3,1) + t331 * mrSges(3,2) + (t325 * t350 - t326 * t352) * t376 + t177;
t166 = t175 * t382 - t176 * t346 + t180 * t381;
t364 = mrSges(7,1) * t206 - mrSges(7,3) * t203 - Ifges(7,4) * t246 - Ifges(7,2) * t245 + Ifges(7,6) * t297 + t318 * t257;
t251 = Ifges(7,5) * t307 + Ifges(7,6) * t306 - Ifges(7,3) * t318;
t363 = mrSges(7,2) * t206 - mrSges(7,3) * t199 + Ifges(7,1) * t246 + Ifges(7,4) * t245 - Ifges(7,5) * t297 + t306 * t251;
t358 = mrSges(6,1) * t212 - mrSges(6,2) * t209 + pkin(5) * t197 + qJ(6) * t196 - t364;
t172 = (t251 + t380) * t307 + (Ifges(5,6) - Ifges(6,6)) * t297 + (Ifges(5,4) - Ifges(6,5)) * t246 + mrSges(5,1) * t362 + mrSges(5,3) * t214 - pkin(4) * t194 + (-Ifges(5,2) - Ifges(6,3)) * t245 + (t259 + t258) * t318 - t358;
t357 = mrSges(6,2) * t211 - mrSges(6,3) * t212 + Ifges(6,1) * t246 + Ifges(6,4) * t297 + Ifges(6,5) * t245 - qJ(6) * t195 + t318 * t252 + t363;
t178 = Ifges(5,5) * t297 - Ifges(5,4) * t245 + Ifges(5,1) * t246 - mrSges(5,2) * t362 - mrSges(5,3) * t213 + t357 - qJ(5) * t194 + (-t256 + t254) * t318 + t380 * t306;
t293 = Ifges(4,5) * t321 + Ifges(4,6) * t320 + Ifges(4,3) * t338;
t294 = Ifges(4,4) * t321 + Ifges(4,2) * t320 + Ifges(4,6) * t338;
t167 = mrSges(4,2) * t268 - mrSges(4,3) * t221 + Ifges(4,1) * t300 + Ifges(4,4) * t299 + Ifges(4,5) * t324 - pkin(10) * t184 - t348 * t172 + t178 * t391 + t320 * t293 - t338 * t294;
t295 = Ifges(4,1) * t321 + Ifges(4,4) * t320 + Ifges(4,5) * t338;
t168 = -mrSges(4,1) * t268 + mrSges(4,3) * t222 + Ifges(4,4) * t300 + Ifges(4,2) * t299 + Ifges(4,6) * t324 - pkin(3) * t184 - t321 * t293 + t338 * t295 - t395;
t311 = Ifges(3,6) * t342 + (Ifges(3,4) * t350 + Ifges(3,2) * t352) * t376;
t312 = Ifges(3,5) * t342 + (Ifges(3,1) * t350 + Ifges(3,4) * t352) * t376;
t158 = Ifges(3,5) * t331 - Ifges(3,6) * t332 + Ifges(3,3) * t368 + mrSges(3,1) * t301 - mrSges(3,2) * t302 + t349 * t167 + t392 * t168 + pkin(2) * t360 + pkin(9) * t370 + (t311 * t350 - t312 * t352) * t376;
t310 = Ifges(3,3) * t342 + (Ifges(3,5) * t350 + Ifges(3,6) * t352) * t376;
t160 = mrSges(3,2) * t313 - mrSges(3,3) * t301 + Ifges(3,1) * t331 - Ifges(3,4) * t332 + Ifges(3,5) * t368 - pkin(9) * t177 + t167 * t392 - t349 * t168 + t310 * t371 - t342 * t311;
t356 = mrSges(4,1) * t221 - mrSges(4,2) * t222 + Ifges(4,5) * t300 + Ifges(4,6) * t299 + Ifges(4,3) * t324 + pkin(3) * t191 + pkin(10) * t185 + t172 * t391 + t348 * t178 + t321 * t294 - t320 * t295;
t162 = -mrSges(3,1) * t313 + mrSges(3,3) * t302 + Ifges(3,4) * t331 - Ifges(3,2) * t332 + Ifges(3,6) * t368 - pkin(2) * t177 - t310 * t372 + t342 * t312 - t356;
t361 = mrSges(2,1) * t339 - mrSges(2,2) * t340 + Ifges(2,3) * qJDD(1) + pkin(1) * t166 + t347 * t158 + t160 * t384 + t162 * t383 + t171 * t389;
t169 = m(2) * t340 - mrSges(2,1) * t354 - qJDD(1) * mrSges(2,2) + t171;
t165 = t347 * t176 + (t175 * t350 + t180 * t352) * t346;
t163 = m(2) * t339 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t354 + t166;
t156 = -mrSges(2,2) * g(3) - mrSges(2,3) * t339 + Ifges(2,5) * qJDD(1) - t354 * Ifges(2,6) + t352 * t160 - t350 * t162 + (-t165 * t346 - t166 * t347) * pkin(8);
t155 = mrSges(2,1) * g(3) + mrSges(2,3) * t340 + t354 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t165 - t346 * t158 + (pkin(8) * t171 + t160 * t350 + t162 * t352) * t347;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t353 * t156 - t351 * t155 - pkin(7) * (t163 * t353 + t169 * t351), t156, t160, t167, t178, -t306 * t255 + t357 + t385, t363 + t385; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t351 * t156 + t353 * t155 + pkin(7) * (-t163 * t351 + t169 * t353), t155, t162, t168, t172, -t307 * t252 - t359, -t307 * t251 - t364; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t361, t361, t158, t356, t395, -t318 * t258 + Ifges(6,6) * t297 + Ifges(6,5) * t246 + Ifges(6,3) * t245 + (-t251 + t255) * t307 + t358, t365;];
m_new  = t1;
