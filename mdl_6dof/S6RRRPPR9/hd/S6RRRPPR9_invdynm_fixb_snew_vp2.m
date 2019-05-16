% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-05-07 06:43
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPPR9_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR9_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR9_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR9_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR9_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR9_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR9_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 06:32:05
% EndTime: 2019-05-07 06:33:00
% DurationCPUTime: 27.91s
% Computational Cost: add. (500510->397), mult. (1094947->499), div. (0->0), fcn. (858669->12), ass. (0->158)
t398 = -2 * qJD(4);
t350 = sin(pkin(6));
t354 = sin(qJ(2));
t357 = cos(qJ(2));
t377 = qJD(1) * qJD(2);
t334 = (-qJDD(1) * t357 + t354 * t377) * t350;
t380 = qJD(1) * t350;
t332 = (-pkin(2) * t357 - pkin(9) * t354) * t380;
t351 = cos(pkin(6));
t345 = qJD(1) * t351 + qJD(2);
t343 = t345 ^ 2;
t344 = qJDD(1) * t351 + qJDD(2);
t379 = qJD(1) * t357;
t355 = sin(qJ(1));
t358 = cos(qJ(1));
t341 = t355 * g(1) - g(2) * t358;
t359 = qJD(1) ^ 2;
t393 = pkin(8) * t350;
t329 = qJDD(1) * pkin(1) + t359 * t393 + t341;
t342 = -g(1) * t358 - g(2) * t355;
t330 = -pkin(1) * t359 + qJDD(1) * t393 + t342;
t386 = t351 * t354;
t381 = t329 * t386 + t357 * t330;
t264 = -t343 * pkin(2) + t344 * pkin(9) + (-g(3) * t354 + t332 * t379) * t350 + t381;
t333 = (qJDD(1) * t354 + t357 * t377) * t350;
t392 = t351 * g(3);
t265 = t334 * pkin(2) - t333 * pkin(9) - t392 + (-t329 + (pkin(2) * t354 - pkin(9) * t357) * t345 * qJD(1)) * t350;
t353 = sin(qJ(3));
t394 = cos(qJ(3));
t238 = t394 * t264 + t353 * t265;
t376 = t354 * t380;
t322 = -t394 * t345 + t353 * t376;
t323 = t353 * t345 + t376 * t394;
t303 = pkin(3) * t322 - qJ(4) * t323;
t326 = qJDD(3) + t334;
t375 = t350 * t379;
t340 = qJD(3) - t375;
t339 = t340 ^ 2;
t226 = -pkin(3) * t339 + qJ(4) * t326 - t303 * t322 + t238;
t385 = t351 * t357;
t387 = t350 * t357;
t301 = -g(3) * t387 + t329 * t385 - t354 * t330;
t263 = -t344 * pkin(2) - t343 * pkin(9) + t332 * t376 - t301;
t299 = qJD(3) * t323 + t333 * t353 - t394 * t344;
t300 = -t322 * qJD(3) + t333 * t394 + t353 * t344;
t230 = (t322 * t340 - t300) * qJ(4) + (t323 * t340 + t299) * pkin(3) + t263;
t349 = sin(pkin(11));
t390 = cos(pkin(11));
t309 = t323 * t390 + t349 * t340;
t221 = -t349 * t226 + t230 * t390 + t309 * t398;
t277 = t300 * t390 + t349 * t326;
t237 = -t353 * t264 + t394 * t265;
t368 = t326 * pkin(3) + t339 * qJ(4) - t323 * t303 - qJDD(4) + t237;
t308 = t323 * t349 - t340 * t390;
t389 = t308 * t322;
t397 = (-t277 + t389) * qJ(5) - t368;
t222 = t390 * t226 + t349 * t230 + t308 * t398;
t251 = Ifges(6,5) * t309 + Ifges(6,6) * t322 + Ifges(6,3) * t308;
t254 = Ifges(5,4) * t309 - Ifges(5,2) * t308 + Ifges(5,6) * t322;
t256 = Ifges(5,1) * t309 - Ifges(5,4) * t308 + Ifges(5,5) * t322;
t276 = t300 * t349 - t326 * t390;
t279 = mrSges(6,1) * t308 - mrSges(6,3) * t309;
t278 = pkin(4) * t308 - qJ(5) * t309;
t321 = t322 ^ 2;
t219 = -t299 * pkin(4) - t321 * qJ(5) + t309 * t278 + qJDD(5) - t221;
t211 = (-t277 - t389) * pkin(10) + (t308 * t309 - t299) * pkin(5) + t219;
t395 = 2 * qJD(5);
t217 = -pkin(4) * t321 + t299 * qJ(5) - t308 * t278 + t322 * t395 + t222;
t286 = -pkin(5) * t322 - pkin(10) * t309;
t307 = t308 ^ 2;
t212 = -pkin(5) * t307 + pkin(10) * t276 + t286 * t322 + t217;
t352 = sin(qJ(6));
t356 = cos(qJ(6));
t208 = t211 * t356 - t212 * t352;
t268 = t308 * t356 - t309 * t352;
t236 = qJD(6) * t268 + t276 * t352 + t277 * t356;
t269 = t308 * t352 + t309 * t356;
t244 = -mrSges(7,1) * t268 + mrSges(7,2) * t269;
t319 = qJD(6) - t322;
t249 = -mrSges(7,2) * t319 + mrSges(7,3) * t268;
t297 = qJDD(6) - t299;
t203 = m(7) * t208 + mrSges(7,1) * t297 - mrSges(7,3) * t236 - t244 * t269 + t249 * t319;
t209 = t211 * t352 + t212 * t356;
t235 = -qJD(6) * t269 + t276 * t356 - t277 * t352;
t250 = mrSges(7,1) * t319 - mrSges(7,3) * t269;
t204 = m(7) * t209 - mrSges(7,2) * t297 + mrSges(7,3) * t235 + t244 * t268 - t250 * t319;
t194 = t356 * t203 + t352 * t204;
t255 = Ifges(6,1) * t309 + Ifges(6,4) * t322 + Ifges(6,5) * t308;
t240 = Ifges(7,4) * t269 + Ifges(7,2) * t268 + Ifges(7,6) * t319;
t241 = Ifges(7,1) * t269 + Ifges(7,4) * t268 + Ifges(7,5) * t319;
t369 = mrSges(7,1) * t208 - mrSges(7,2) * t209 + Ifges(7,5) * t236 + Ifges(7,6) * t235 + Ifges(7,3) * t297 + t269 * t240 - t268 * t241;
t362 = mrSges(6,1) * t219 - mrSges(6,3) * t217 - Ifges(6,4) * t277 - Ifges(6,2) * t299 - Ifges(6,6) * t276 + pkin(5) * t194 - t308 * t255 + t369;
t282 = -mrSges(6,2) * t308 + mrSges(6,3) * t322;
t367 = -m(6) * t219 + t299 * mrSges(6,1) + t322 * t282 - t194;
t195 = -t352 * t203 + t356 * t204;
t285 = -mrSges(6,1) * t322 + mrSges(6,2) * t309;
t371 = m(6) * t217 + t299 * mrSges(6,3) + t322 * t285 + t195;
t396 = -(t254 - t251) * t309 - mrSges(5,1) * t221 + mrSges(5,2) * t222 - Ifges(5,5) * t277 + Ifges(5,6) * t276 - pkin(4) * (-t277 * mrSges(6,2) - t309 * t279 + t367) - qJ(5) * (-t276 * mrSges(6,2) - t308 * t279 + t371) - t308 * t256 + t362;
t391 = -mrSges(5,3) - mrSges(6,2);
t388 = t350 * t354;
t284 = mrSges(5,1) * t322 - mrSges(5,3) * t309;
t382 = -mrSges(5,1) * t308 - mrSges(5,2) * t309 - t279;
t190 = m(5) * t222 - t299 * mrSges(5,2) + t276 * t391 - t322 * t284 + t308 * t382 + t371;
t283 = -mrSges(5,2) * t322 - mrSges(5,3) * t308;
t191 = m(5) * t221 + t299 * mrSges(5,1) + t277 * t391 + t322 * t283 + t309 * t382 + t367;
t188 = t390 * t190 - t191 * t349;
t304 = mrSges(4,1) * t322 + mrSges(4,2) * t323;
t311 = mrSges(4,1) * t340 - mrSges(4,3) * t323;
t186 = m(4) * t238 - mrSges(4,2) * t326 - mrSges(4,3) * t299 - t304 * t322 - t311 * t340 + t188;
t214 = -t307 * pkin(10) + (-pkin(4) - pkin(5)) * t276 + (-pkin(4) * t322 + t286 + t395) * t309 - t397;
t210 = -m(7) * t214 + t235 * mrSges(7,1) - t236 * mrSges(7,2) + t268 * t249 - t269 * t250;
t220 = -0.2e1 * qJD(5) * t309 + (t309 * t322 + t276) * pkin(4) + t397;
t205 = m(6) * t220 + mrSges(6,1) * t276 - t277 * mrSges(6,3) + t282 * t308 - t309 * t285 + t210;
t201 = m(5) * t368 - t276 * mrSges(5,1) - mrSges(5,2) * t277 - t308 * t283 - t284 * t309 - t205;
t310 = -mrSges(4,2) * t340 - mrSges(4,3) * t322;
t200 = m(4) * t237 + mrSges(4,1) * t326 - mrSges(4,3) * t300 - t304 * t323 + t310 * t340 + t201;
t181 = t353 * t186 + t394 * t200;
t253 = Ifges(6,4) * t309 + Ifges(6,2) * t322 + Ifges(6,6) * t308;
t384 = -Ifges(5,5) * t309 + Ifges(5,6) * t308 - Ifges(5,3) * t322 - t253;
t302 = -g(3) * t388 + t381;
t327 = mrSges(3,1) * t345 - mrSges(3,3) * t376;
t331 = (-mrSges(3,1) * t357 + mrSges(3,2) * t354) * t380;
t374 = t394 * t186 - t200 * t353;
t178 = m(3) * t302 - mrSges(3,2) * t344 - mrSges(3,3) * t334 - t327 * t345 + t331 * t375 + t374;
t328 = -mrSges(3,2) * t345 + mrSges(3,3) * t375;
t187 = t349 * t190 + t191 * t390;
t363 = -m(4) * t263 - t299 * mrSges(4,1) - t300 * mrSges(4,2) - t322 * t310 - t323 * t311 - t187;
t183 = m(3) * t301 + t344 * mrSges(3,1) - t333 * mrSges(3,3) + t345 * t328 - t331 * t376 + t363;
t174 = t357 * t178 - t183 * t354;
t315 = -t350 * t329 - t392;
t179 = m(3) * t315 + t334 * mrSges(3,1) + t333 * mrSges(3,2) + (t327 * t354 - t328 * t357) * t380 + t181;
t170 = t178 * t386 - t179 * t350 + t183 * t385;
t239 = Ifges(7,5) * t269 + Ifges(7,6) * t268 + Ifges(7,3) * t319;
t197 = -mrSges(7,1) * t214 + mrSges(7,3) * t209 + Ifges(7,4) * t236 + Ifges(7,2) * t235 + Ifges(7,6) * t297 - t239 * t269 + t241 * t319;
t198 = mrSges(7,2) * t214 - mrSges(7,3) * t208 + Ifges(7,1) * t236 + Ifges(7,4) * t235 + Ifges(7,5) * t297 + t239 * t268 - t240 * t319;
t364 = -mrSges(6,1) * t220 + mrSges(6,2) * t217 - pkin(5) * t210 - pkin(10) * t195 - t356 * t197 - t352 * t198;
t175 = mrSges(5,1) * t368 + mrSges(5,3) * t222 - pkin(4) * t205 + (t256 + t255) * t322 + t384 * t309 + (Ifges(5,6) - Ifges(6,6)) * t299 + (Ifges(5,4) - Ifges(6,5)) * t277 + (-Ifges(5,2) - Ifges(6,3)) * t276 + t364;
t365 = mrSges(6,2) * t219 - mrSges(6,3) * t220 + Ifges(6,1) * t277 + Ifges(6,4) * t299 + Ifges(6,5) * t276 - pkin(10) * t194 - t352 * t197 + t356 * t198 + t322 * t251;
t180 = -mrSges(5,2) * t368 - mrSges(5,3) * t221 + Ifges(5,1) * t277 - Ifges(5,4) * t276 + Ifges(5,5) * t299 - qJ(5) * t205 - t322 * t254 + t308 * t384 + t365;
t288 = Ifges(4,5) * t323 - Ifges(4,6) * t322 + Ifges(4,3) * t340;
t289 = Ifges(4,4) * t323 - Ifges(4,2) * t322 + Ifges(4,6) * t340;
t166 = mrSges(4,2) * t263 - mrSges(4,3) * t237 + Ifges(4,1) * t300 - Ifges(4,4) * t299 + Ifges(4,5) * t326 - qJ(4) * t187 - t349 * t175 + t180 * t390 - t322 * t288 - t340 * t289;
t290 = Ifges(4,1) * t323 - Ifges(4,4) * t322 + Ifges(4,5) * t340;
t171 = -pkin(3) * t187 + (-Ifges(4,2) - Ifges(5,3)) * t299 + t340 * t290 + Ifges(4,6) * t326 - t323 * t288 + Ifges(4,4) * t300 - mrSges(4,1) * t263 + mrSges(4,3) * t238 + t396;
t313 = Ifges(3,6) * t345 + (Ifges(3,4) * t354 + Ifges(3,2) * t357) * t380;
t314 = Ifges(3,5) * t345 + (Ifges(3,1) * t354 + Ifges(3,4) * t357) * t380;
t161 = Ifges(3,5) * t333 - Ifges(3,6) * t334 + Ifges(3,3) * t344 + mrSges(3,1) * t301 - mrSges(3,2) * t302 + t353 * t166 + t394 * t171 + pkin(2) * t363 + pkin(9) * t374 + (t313 * t354 - t314 * t357) * t380;
t312 = Ifges(3,3) * t345 + (Ifges(3,5) * t354 + Ifges(3,6) * t357) * t380;
t163 = mrSges(3,2) * t315 - mrSges(3,3) * t301 + Ifges(3,1) * t333 - Ifges(3,4) * t334 + Ifges(3,5) * t344 - pkin(9) * t181 + t166 * t394 - t353 * t171 + t312 * t375 - t345 * t313;
t361 = mrSges(4,1) * t237 - mrSges(4,2) * t238 + Ifges(4,5) * t300 - Ifges(4,6) * t299 + Ifges(4,3) * t326 + pkin(3) * t201 + qJ(4) * t188 + t175 * t390 + t349 * t180 + t323 * t289 + t322 * t290;
t165 = -mrSges(3,1) * t315 + mrSges(3,3) * t302 + Ifges(3,4) * t333 - Ifges(3,2) * t334 + Ifges(3,6) * t344 - pkin(2) * t181 - t312 * t376 + t345 * t314 - t361;
t366 = mrSges(2,1) * t341 - mrSges(2,2) * t342 + Ifges(2,3) * qJDD(1) + pkin(1) * t170 + t351 * t161 + t163 * t388 + t165 * t387 + t174 * t393;
t172 = m(2) * t342 - mrSges(2,1) * t359 - qJDD(1) * mrSges(2,2) + t174;
t169 = t351 * t179 + (t178 * t354 + t183 * t357) * t350;
t167 = m(2) * t341 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t359 + t170;
t159 = -mrSges(2,2) * g(3) - mrSges(2,3) * t341 + Ifges(2,5) * qJDD(1) - t359 * Ifges(2,6) + t357 * t163 - t354 * t165 + (-t169 * t350 - t170 * t351) * pkin(8);
t158 = mrSges(2,1) * g(3) + mrSges(2,3) * t342 + t359 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t169 - t350 * t161 + (pkin(8) * t174 + t163 * t354 + t165 * t357) * t351;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t358 * t159 - t355 * t158 - pkin(7) * (t167 * t358 + t172 * t355), t159, t163, t166, t180, -t253 * t308 + t365, t198; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t355 * t159 + t358 * t158 + pkin(7) * (-t167 * t355 + t172 * t358), t158, t165, t171, t175, -t309 * t251 - t362, t197; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t366, t366, t161, t361, Ifges(5,3) * t299 - t396, Ifges(6,5) * t277 + Ifges(6,6) * t299 + Ifges(6,3) * t276 + t309 * t253 - t322 * t255 - t364, t369;];
m_new  = t1;
