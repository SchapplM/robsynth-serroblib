% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRRP12
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 07:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRRP12_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP12_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP12_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP12_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP12_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP12_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP12_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 07:18:23
% EndTime: 2019-05-08 07:22:02
% DurationCPUTime: 98.81s
% Computational Cost: add. (1751126->412), mult. (4316653->537), div. (0->0), fcn. (3628284->14), ass. (0->171)
t345 = cos(pkin(6));
t339 = qJD(1) * t345 + qJD(2);
t342 = sin(pkin(7));
t344 = cos(pkin(7));
t343 = sin(pkin(6));
t353 = cos(qJ(2));
t376 = qJD(1) * t353;
t371 = t343 * t376;
t323 = (t339 * t342 + t344 * t371) * pkin(10);
t349 = sin(qJ(2));
t378 = qJD(1) * t343;
t394 = pkin(10) * t342;
t327 = (-pkin(2) * t353 - t349 * t394) * t378;
t375 = qJD(1) * qJD(2);
t333 = (qJDD(1) * t349 + t353 * t375) * t343;
t338 = qJDD(1) * t345 + qJDD(2);
t350 = sin(qJ(1));
t354 = cos(qJ(1));
t336 = t350 * g(1) - g(2) * t354;
t355 = qJD(1) ^ 2;
t395 = pkin(9) * t343;
t330 = qJDD(1) * pkin(1) + t355 * t395 + t336;
t337 = -g(1) * t354 - g(2) * t350;
t331 = -pkin(1) * t355 + qJDD(1) * t395 + t337;
t383 = t345 * t353;
t369 = t330 * t383 - t349 * t331;
t377 = qJD(1) * t349;
t393 = pkin(10) * t344;
t281 = -t333 * t393 + t338 * pkin(2) + t339 * t323 + (-g(3) * t353 - t327 * t377) * t343 + t369;
t372 = t343 * t377;
t326 = pkin(2) * t339 - t372 * t393;
t334 = (qJDD(1) * t353 - t349 * t375) * t343;
t365 = t334 * t344 + t338 * t342;
t384 = t345 * t349;
t379 = t330 * t384 + t353 * t331;
t282 = -t339 * t326 + (-g(3) * t349 + t327 * t376) * t343 + t365 * pkin(10) + t379;
t392 = t345 * g(3);
t287 = -t333 * t394 - t334 * pkin(2) - t392 + (-t330 + (-t323 * t353 + t326 * t349) * qJD(1)) * t343;
t348 = sin(qJ(3));
t352 = cos(qJ(3));
t243 = -t348 * t282 + (t281 * t344 + t287 * t342) * t352;
t386 = t344 * t348;
t390 = t342 * t348;
t244 = t281 * t386 + t352 * t282 + t287 * t390;
t385 = t344 * t353;
t389 = t342 * t352;
t313 = (-t348 * t349 + t352 * t385) * t378 + t339 * t389;
t314 = t339 * t390 + (t348 * t385 + t349 * t352) * t378;
t301 = -pkin(3) * t313 - pkin(11) * t314;
t315 = -t334 * t342 + t338 * t344 + qJDD(3);
t324 = t339 * t344 - t342 * t371 + qJD(3);
t322 = t324 ^ 2;
t232 = -pkin(3) * t322 + pkin(11) * t315 + t301 * t313 + t244;
t255 = -t342 * t281 + t344 * t287;
t298 = -t314 * qJD(3) - t348 * t333 + t365 * t352;
t299 = t313 * qJD(3) + t352 * t333 + t365 * t348;
t234 = (-t313 * t324 - t299) * pkin(11) + (t314 * t324 - t298) * pkin(3) + t255;
t347 = sin(qJ(4));
t351 = cos(qJ(4));
t228 = t351 * t232 + t347 * t234;
t304 = -t314 * t347 + t324 * t351;
t305 = t314 * t351 + t324 * t347;
t284 = -pkin(4) * t304 - pkin(12) * t305;
t297 = qJDD(4) - t298;
t312 = qJD(4) - t313;
t311 = t312 ^ 2;
t224 = -pkin(4) * t311 + pkin(12) * t297 + t284 * t304 + t228;
t231 = -t315 * pkin(3) - t322 * pkin(11) + t314 * t301 - t243;
t267 = -qJD(4) * t305 - t299 * t347 + t315 * t351;
t268 = qJD(4) * t304 + t299 * t351 + t315 * t347;
t226 = (-t304 * t312 - t268) * pkin(12) + (t305 * t312 - t267) * pkin(4) + t231;
t346 = sin(qJ(5));
t396 = cos(qJ(5));
t220 = -t346 * t224 + t396 * t226;
t221 = t396 * t224 + t346 * t226;
t290 = t396 * t305 + t346 * t312;
t241 = t290 * qJD(5) + t346 * t268 - t396 * t297;
t289 = t346 * t305 - t396 * t312;
t242 = -t289 * qJD(5) + t396 * t268 + t346 * t297;
t303 = qJD(5) - t304;
t247 = Ifges(7,5) * t290 + Ifges(7,6) * t303 + Ifges(7,3) * t289;
t250 = Ifges(6,4) * t290 - Ifges(6,2) * t289 + Ifges(6,6) * t303;
t252 = Ifges(6,1) * t290 - Ifges(6,4) * t289 + Ifges(6,5) * t303;
t259 = mrSges(7,1) * t289 - mrSges(7,3) * t290;
t266 = qJDD(5) - t267;
t258 = pkin(5) * t289 - qJ(6) * t290;
t302 = t303 ^ 2;
t216 = -pkin(5) * t302 + qJ(6) * t266 + 0.2e1 * qJD(6) * t303 - t258 * t289 + t221;
t218 = -t266 * pkin(5) - t302 * qJ(6) + t290 * t258 + qJDD(6) - t220;
t251 = Ifges(7,1) * t290 + Ifges(7,4) * t303 + Ifges(7,5) * t289;
t361 = mrSges(7,1) * t218 - mrSges(7,3) * t216 - Ifges(7,4) * t242 - Ifges(7,2) * t266 - Ifges(7,6) * t241 - t289 * t251;
t270 = -mrSges(7,2) * t289 + mrSges(7,3) * t303;
t368 = -m(7) * t218 + t266 * mrSges(7,1) + t303 * t270;
t273 = -mrSges(7,1) * t303 + mrSges(7,2) * t290;
t373 = m(7) * t216 + t266 * mrSges(7,3) + t303 * t273;
t397 = -(-t250 + t247) * t290 + mrSges(6,1) * t220 - mrSges(6,2) * t221 + Ifges(6,5) * t242 - Ifges(6,6) * t241 + Ifges(6,3) * t266 + pkin(5) * (-t242 * mrSges(7,2) - t290 * t259 + t368) + qJ(6) * (-t241 * mrSges(7,2) - t289 * t259 + t373) + t289 * t252 - t361;
t391 = -mrSges(6,3) - mrSges(7,2);
t388 = t343 * t349;
t387 = t343 * t353;
t272 = mrSges(6,1) * t303 - mrSges(6,3) * t290;
t380 = -mrSges(6,1) * t289 - mrSges(6,2) * t290 - t259;
t208 = m(6) * t221 - t266 * mrSges(6,2) + t391 * t241 - t303 * t272 + t380 * t289 + t373;
t271 = -mrSges(6,2) * t303 - mrSges(6,3) * t289;
t209 = m(6) * t220 + t266 * mrSges(6,1) + t391 * t242 + t303 * t271 + t380 * t290 + t368;
t204 = t396 * t208 - t209 * t346;
t283 = -mrSges(5,1) * t304 + mrSges(5,2) * t305;
t292 = mrSges(5,1) * t312 - mrSges(5,3) * t305;
t200 = m(5) * t228 - mrSges(5,2) * t297 + mrSges(5,3) * t267 + t283 * t304 - t292 * t312 + t204;
t227 = -t347 * t232 + t351 * t234;
t223 = -t297 * pkin(4) - t311 * pkin(12) + t305 * t284 - t227;
t219 = -0.2e1 * qJD(6) * t290 + (t289 * t303 - t242) * qJ(6) + (t290 * t303 + t241) * pkin(5) + t223;
t213 = m(7) * t219 + mrSges(7,1) * t241 - t242 * mrSges(7,3) + t270 * t289 - t290 * t273;
t210 = -m(6) * t223 - t241 * mrSges(6,1) - mrSges(6,2) * t242 - t289 * t271 - t272 * t290 - t213;
t291 = -mrSges(5,2) * t312 + mrSges(5,3) * t304;
t206 = m(5) * t227 + mrSges(5,1) * t297 - mrSges(5,3) * t268 - t283 * t305 + t291 * t312 + t210;
t194 = t347 * t200 + t351 * t206;
t249 = Ifges(7,4) * t290 + Ifges(7,2) * t303 + Ifges(7,6) * t289;
t382 = -Ifges(6,5) * t290 + Ifges(6,6) * t289 - Ifges(6,3) * t303 - t249;
t300 = -mrSges(4,1) * t313 + mrSges(4,2) * t314;
t307 = mrSges(4,1) * t324 - mrSges(4,3) * t314;
t370 = t351 * t200 - t206 * t347;
t191 = m(4) * t244 - mrSges(4,2) * t315 + mrSges(4,3) * t298 + t300 * t313 - t307 * t324 + t370;
t306 = -mrSges(4,2) * t324 + mrSges(4,3) * t313;
t193 = m(4) * t255 - mrSges(4,1) * t298 + mrSges(4,2) * t299 - t306 * t313 + t307 * t314 + t194;
t203 = t346 * t208 + t396 * t209;
t358 = -m(5) * t231 + t267 * mrSges(5,1) - t268 * mrSges(5,2) + t304 * t291 - t305 * t292 - t203;
t197 = m(4) * t243 + t315 * mrSges(4,1) - t299 * mrSges(4,3) - t314 * t300 + t324 * t306 + t358;
t180 = t191 * t390 + t344 * t193 + t197 * t389;
t181 = t344 * t352 * t197 + t191 * t386 - t193 * t342;
t308 = -g(3) * t387 + t369;
t329 = -mrSges(3,2) * t339 + mrSges(3,3) * t371;
t332 = (-mrSges(3,1) * t353 + mrSges(3,2) * t349) * t378;
t178 = m(3) * t308 + mrSges(3,1) * t338 - mrSges(3,3) * t333 + t329 * t339 - t332 * t372 + t181;
t186 = t352 * t191 - t197 * t348;
t309 = -g(3) * t388 + t379;
t328 = mrSges(3,1) * t339 - mrSges(3,3) * t372;
t185 = m(3) * t309 - mrSges(3,2) * t338 + mrSges(3,3) * t334 - t328 * t339 + t332 * t371 + t186;
t175 = -t178 * t349 + t353 * t185;
t319 = -t343 * t330 - t392;
t179 = m(3) * t319 - t334 * mrSges(3,1) + t333 * mrSges(3,2) + (t328 * t349 - t329 * t353) * t378 + t180;
t170 = t178 * t383 - t179 * t343 + t185 * t384;
t367 = -mrSges(7,1) * t219 + mrSges(7,2) * t216;
t201 = -mrSges(6,1) * t223 + mrSges(6,3) * t221 - pkin(5) * t213 + (t251 + t252) * t303 + t382 * t290 + (Ifges(6,6) - Ifges(7,6)) * t266 + (Ifges(6,4) - Ifges(7,5)) * t242 + (-Ifges(6,2) - Ifges(7,3)) * t241 + t367;
t360 = mrSges(7,2) * t218 - mrSges(7,3) * t219 + Ifges(7,1) * t242 + Ifges(7,4) * t266 + Ifges(7,5) * t241 + t303 * t247;
t202 = mrSges(6,2) * t223 - mrSges(6,3) * t220 + Ifges(6,1) * t242 - Ifges(6,4) * t241 + Ifges(6,5) * t266 - qJ(6) * t213 - t303 * t250 + t382 * t289 + t360;
t274 = Ifges(5,5) * t305 + Ifges(5,6) * t304 + Ifges(5,3) * t312;
t275 = Ifges(5,4) * t305 + Ifges(5,2) * t304 + Ifges(5,6) * t312;
t182 = mrSges(5,2) * t231 - mrSges(5,3) * t227 + Ifges(5,1) * t268 + Ifges(5,4) * t267 + Ifges(5,5) * t297 - pkin(12) * t203 - t346 * t201 + t396 * t202 + t304 * t274 - t312 * t275;
t276 = Ifges(5,1) * t305 + Ifges(5,4) * t304 + Ifges(5,5) * t312;
t187 = -mrSges(5,1) * t231 + mrSges(5,3) * t228 + Ifges(5,4) * t268 + Ifges(5,2) * t267 + Ifges(5,6) * t297 - pkin(4) * t203 - t305 * t274 + t312 * t276 - t397;
t293 = Ifges(4,5) * t314 + Ifges(4,6) * t313 + Ifges(4,3) * t324;
t294 = Ifges(4,4) * t314 + Ifges(4,2) * t313 + Ifges(4,6) * t324;
t172 = mrSges(4,2) * t255 - mrSges(4,3) * t243 + Ifges(4,1) * t299 + Ifges(4,4) * t298 + Ifges(4,5) * t315 - pkin(11) * t194 + t182 * t351 - t187 * t347 + t293 * t313 - t294 * t324;
t295 = Ifges(4,1) * t314 + Ifges(4,4) * t313 + Ifges(4,5) * t324;
t356 = mrSges(5,1) * t227 - mrSges(5,2) * t228 + Ifges(5,5) * t268 + Ifges(5,6) * t267 + Ifges(5,3) * t297 + pkin(4) * t210 + pkin(12) * t204 + t396 * t201 + t346 * t202 + t305 * t275 - t304 * t276;
t176 = -mrSges(4,1) * t255 + mrSges(4,3) * t244 + Ifges(4,4) * t299 + Ifges(4,2) * t298 + Ifges(4,6) * t315 - pkin(3) * t194 - t314 * t293 + t324 * t295 - t356;
t362 = pkin(10) * t186 + t172 * t348 + t176 * t352;
t171 = mrSges(4,1) * t243 - mrSges(4,2) * t244 + Ifges(4,5) * t299 + Ifges(4,6) * t298 + Ifges(4,3) * t315 + pkin(3) * t358 + pkin(11) * t370 + t347 * t182 + t351 * t187 + t314 * t294 - t313 * t295;
t317 = Ifges(3,6) * t339 + (Ifges(3,4) * t349 + Ifges(3,2) * t353) * t378;
t318 = Ifges(3,5) * t339 + (Ifges(3,1) * t349 + Ifges(3,4) * t353) * t378;
t162 = mrSges(3,1) * t308 - mrSges(3,2) * t309 + Ifges(3,5) * t333 + Ifges(3,6) * t334 + Ifges(3,3) * t338 + pkin(2) * t181 + t344 * t171 + (t317 * t349 - t318 * t353) * t378 + t362 * t342;
t316 = Ifges(3,3) * t339 + (Ifges(3,5) * t349 + Ifges(3,6) * t353) * t378;
t164 = -mrSges(3,1) * t319 + mrSges(3,3) * t309 + Ifges(3,4) * t333 + Ifges(3,2) * t334 + Ifges(3,6) * t338 - pkin(2) * t180 - t342 * t171 - t316 * t372 + t339 * t318 + t362 * t344;
t166 = t316 * t371 + mrSges(3,2) * t319 - mrSges(3,3) * t308 + Ifges(3,1) * t333 + Ifges(3,4) * t334 + Ifges(3,5) * t338 + t352 * t172 - t348 * t176 - t339 * t317 + (-t180 * t342 - t181 * t344) * pkin(10);
t359 = mrSges(2,1) * t336 - mrSges(2,2) * t337 + Ifges(2,3) * qJDD(1) + pkin(1) * t170 + t345 * t162 + t164 * t387 + t166 * t388 + t175 * t395;
t173 = m(2) * t337 - mrSges(2,1) * t355 - qJDD(1) * mrSges(2,2) + t175;
t169 = t345 * t179 + (t178 * t353 + t185 * t349) * t343;
t167 = m(2) * t336 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t355 + t170;
t160 = -mrSges(2,2) * g(3) - mrSges(2,3) * t336 + Ifges(2,5) * qJDD(1) - t355 * Ifges(2,6) - t349 * t164 + t353 * t166 + (-t169 * t343 - t170 * t345) * pkin(9);
t159 = mrSges(2,1) * g(3) + mrSges(2,3) * t337 + t355 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t169 - t343 * t162 + (pkin(9) * t175 + t164 * t353 + t166 * t349) * t345;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t354 * t160 - t350 * t159 - pkin(8) * (t167 * t354 + t173 * t350), t160, t166, t172, t182, t202, -t249 * t289 + t360; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t350 * t160 + t354 * t159 + pkin(8) * (-t167 * t350 + t173 * t354), t159, t164, t176, t187, t201, -t290 * t247 - t361; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t359, t359, t162, t171, t356, t397, Ifges(7,5) * t242 + Ifges(7,6) * t266 + Ifges(7,3) * t241 + t290 * t249 - t303 * t251 - t367;];
m_new  = t1;
