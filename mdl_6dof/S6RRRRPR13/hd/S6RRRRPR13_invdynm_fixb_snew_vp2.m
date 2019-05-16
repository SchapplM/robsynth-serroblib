% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRPR13
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
% Datum: 2019-05-08 01:41
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRPR13_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR13_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR13_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR13_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR13_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR13_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR13_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 01:19:23
% EndTime: 2019-05-08 01:20:43
% DurationCPUTime: 29.22s
% Computational Cost: add. (536271->397), mult. (1134387->498), div. (0->0), fcn. (896370->12), ass. (0->159)
t348 = sin(pkin(6));
t353 = sin(qJ(2));
t357 = cos(qJ(2));
t376 = qJD(1) * qJD(2);
t333 = (-qJDD(1) * t357 + t353 * t376) * t348;
t349 = cos(pkin(6));
t344 = qJD(1) * t349 + qJD(2);
t352 = sin(qJ(3));
t356 = cos(qJ(3));
t378 = qJD(1) * t348;
t375 = t353 * t378;
t321 = t356 * t344 - t352 * t375;
t332 = (qJDD(1) * t353 + t357 * t376) * t348;
t343 = qJDD(1) * t349 + qJDD(2);
t300 = t321 * qJD(3) + t332 * t356 + t343 * t352;
t322 = t344 * t352 + t356 * t375;
t377 = qJD(1) * t357;
t374 = t348 * t377;
t339 = qJD(3) - t374;
t351 = sin(qJ(4));
t391 = cos(qJ(4));
t306 = t322 * t351 - t391 * t339;
t325 = qJDD(3) + t333;
t253 = -t306 * qJD(4) + t391 * t300 + t351 * t325;
t331 = (-pkin(2) * t357 - pkin(9) * t353) * t378;
t342 = t344 ^ 2;
t354 = sin(qJ(1));
t358 = cos(qJ(1));
t340 = t354 * g(1) - g(2) * t358;
t359 = qJD(1) ^ 2;
t390 = pkin(8) * t348;
t328 = qJDD(1) * pkin(1) + t359 * t390 + t340;
t341 = -g(1) * t358 - g(2) * t354;
t329 = -t359 * pkin(1) + qJDD(1) * t390 + t341;
t384 = t349 * t353;
t379 = t328 * t384 + t357 * t329;
t271 = -pkin(2) * t342 + pkin(9) * t343 + (-g(3) * t353 + t331 * t377) * t348 + t379;
t389 = g(3) * t349;
t272 = pkin(2) * t333 - pkin(9) * t332 - t389 + (-t328 + (pkin(2) * t353 - pkin(9) * t357) * t344 * qJD(1)) * t348;
t236 = -t352 * t271 + t356 * t272;
t304 = -t321 * pkin(3) - pkin(10) * t322;
t338 = t339 ^ 2;
t368 = pkin(3) * t325 + pkin(10) * t338 - t322 * t304 + t236;
t319 = qJD(4) - t321;
t387 = t306 * t319;
t394 = (-t253 + t387) * qJ(5) - t368;
t237 = t356 * t271 + t352 * t272;
t231 = -pkin(3) * t338 + pkin(10) * t325 + t321 * t304 + t237;
t383 = t349 * t357;
t385 = t348 * t357;
t301 = -g(3) * t385 + t328 * t383 - t353 * t329;
t270 = -pkin(2) * t343 - pkin(9) * t342 + t331 * t375 - t301;
t299 = -qJD(3) * t322 - t332 * t352 + t356 * t343;
t235 = (-t321 * t339 - t300) * pkin(10) + (t322 * t339 - t299) * pkin(3) + t270;
t220 = -t351 * t231 + t391 * t235;
t221 = t391 * t231 + t351 * t235;
t307 = t391 * t322 + t351 * t339;
t252 = t307 * qJD(4) + t300 * t351 - t391 * t325;
t258 = Ifges(6,5) * t307 + Ifges(6,6) * t319 + Ifges(6,3) * t306;
t261 = Ifges(5,4) * t307 - Ifges(5,2) * t306 + Ifges(5,6) * t319;
t263 = Ifges(5,1) * t307 - Ifges(5,4) * t306 + Ifges(5,5) * t319;
t278 = mrSges(6,1) * t306 - mrSges(6,3) * t307;
t297 = qJDD(4) - t299;
t277 = pkin(4) * t306 - qJ(5) * t307;
t318 = t319 ^ 2;
t218 = -t297 * pkin(4) - t318 * qJ(5) + t307 * t277 + qJDD(5) - t220;
t210 = (-t253 - t387) * pkin(11) + (t306 * t307 - t297) * pkin(5) + t218;
t392 = 2 * qJD(5);
t216 = -pkin(4) * t318 + t297 * qJ(5) - t306 * t277 + t319 * t392 + t221;
t285 = -pkin(5) * t319 - pkin(11) * t307;
t305 = t306 ^ 2;
t211 = -pkin(5) * t305 + pkin(11) * t252 + t285 * t319 + t216;
t350 = sin(qJ(6));
t355 = cos(qJ(6));
t207 = t210 * t355 - t211 * t350;
t275 = t306 * t355 - t307 * t350;
t227 = t275 * qJD(6) + t252 * t350 + t253 * t355;
t276 = t306 * t350 + t307 * t355;
t243 = -mrSges(7,1) * t275 + mrSges(7,2) * t276;
t317 = qJD(6) - t319;
t256 = -mrSges(7,2) * t317 + mrSges(7,3) * t275;
t292 = qJDD(6) - t297;
t203 = m(7) * t207 + mrSges(7,1) * t292 - mrSges(7,3) * t227 - t243 * t276 + t256 * t317;
t208 = t210 * t350 + t211 * t355;
t226 = -t276 * qJD(6) + t252 * t355 - t253 * t350;
t257 = mrSges(7,1) * t317 - mrSges(7,3) * t276;
t204 = m(7) * t208 - mrSges(7,2) * t292 + mrSges(7,3) * t226 + t243 * t275 - t257 * t317;
t193 = t203 * t355 + t204 * t350;
t262 = Ifges(6,1) * t307 + Ifges(6,4) * t319 + Ifges(6,5) * t306;
t239 = Ifges(7,4) * t276 + Ifges(7,2) * t275 + Ifges(7,6) * t317;
t240 = Ifges(7,1) * t276 + Ifges(7,4) * t275 + Ifges(7,5) * t317;
t369 = mrSges(7,1) * t207 - mrSges(7,2) * t208 + Ifges(7,5) * t227 + Ifges(7,6) * t226 + Ifges(7,3) * t292 + t276 * t239 - t275 * t240;
t362 = mrSges(6,1) * t218 - mrSges(6,3) * t216 - Ifges(6,4) * t253 - Ifges(6,2) * t297 - Ifges(6,6) * t252 + pkin(5) * t193 - t306 * t262 + t369;
t281 = -mrSges(6,2) * t306 + mrSges(6,3) * t319;
t367 = -m(6) * t218 + t297 * mrSges(6,1) + t319 * t281 - t193;
t194 = -t203 * t350 + t355 * t204;
t284 = -mrSges(6,1) * t319 + mrSges(6,2) * t307;
t371 = m(6) * t216 + t297 * mrSges(6,3) + t319 * t284 + t194;
t393 = (t261 - t258) * t307 + mrSges(5,1) * t220 - mrSges(5,2) * t221 + Ifges(5,5) * t253 - Ifges(5,6) * t252 + Ifges(5,3) * t297 + pkin(4) * (-t253 * mrSges(6,2) - t307 * t278 + t367) + qJ(5) * (-t252 * mrSges(6,2) - t306 * t278 + t371) + t306 * t263 - t362;
t388 = -mrSges(5,3) - mrSges(6,2);
t386 = t348 * t353;
t283 = mrSges(5,1) * t319 - mrSges(5,3) * t307;
t380 = -mrSges(5,1) * t306 - mrSges(5,2) * t307 - t278;
t189 = m(5) * t221 - t297 * mrSges(5,2) + t388 * t252 - t319 * t283 + t380 * t306 + t371;
t282 = -mrSges(5,2) * t319 - mrSges(5,3) * t306;
t190 = m(5) * t220 + t297 * mrSges(5,1) + t388 * t253 + t319 * t282 + t380 * t307 + t367;
t187 = t391 * t189 - t190 * t351;
t303 = -t321 * mrSges(4,1) + mrSges(4,2) * t322;
t309 = mrSges(4,1) * t339 - mrSges(4,3) * t322;
t185 = m(4) * t237 - mrSges(4,2) * t325 + t299 * mrSges(4,3) + t321 * t303 - t309 * t339 + t187;
t213 = -t305 * pkin(11) + (-pkin(4) - pkin(5)) * t252 + (-pkin(4) * t319 + t285 + t392) * t307 - t394;
t209 = -m(7) * t213 + t226 * mrSges(7,1) - t227 * mrSges(7,2) + t275 * t256 - t276 * t257;
t219 = -0.2e1 * qJD(5) * t307 + (t307 * t319 + t252) * pkin(4) + t394;
t201 = m(6) * t219 + t252 * mrSges(6,1) - t253 * mrSges(6,3) + t306 * t281 - t307 * t284 + t209;
t200 = m(5) * t368 - t252 * mrSges(5,1) - t253 * mrSges(5,2) - t306 * t282 - t307 * t283 - t201;
t308 = -mrSges(4,2) * t339 + t321 * mrSges(4,3);
t199 = m(4) * t236 + mrSges(4,1) * t325 - t300 * mrSges(4,3) - t303 * t322 + t308 * t339 + t200;
t180 = t352 * t185 + t356 * t199;
t260 = Ifges(6,4) * t307 + Ifges(6,2) * t319 + Ifges(6,6) * t306;
t382 = -Ifges(5,5) * t307 + Ifges(5,6) * t306 - Ifges(5,3) * t319 - t260;
t302 = -g(3) * t386 + t379;
t326 = mrSges(3,1) * t344 - mrSges(3,3) * t375;
t330 = (-mrSges(3,1) * t357 + mrSges(3,2) * t353) * t378;
t373 = t356 * t185 - t199 * t352;
t178 = m(3) * t302 - mrSges(3,2) * t343 - mrSges(3,3) * t333 - t326 * t344 + t330 * t374 + t373;
t327 = -mrSges(3,2) * t344 + mrSges(3,3) * t374;
t186 = t351 * t189 + t391 * t190;
t363 = -m(4) * t270 + t299 * mrSges(4,1) - t300 * mrSges(4,2) + t321 * t308 - t322 * t309 - t186;
t182 = m(3) * t301 + t343 * mrSges(3,1) - t332 * mrSges(3,3) + t344 * t327 - t330 * t375 + t363;
t173 = t357 * t178 - t182 * t353;
t313 = -t328 * t348 - t389;
t179 = m(3) * t313 + mrSges(3,1) * t333 + mrSges(3,2) * t332 + (t326 * t353 - t327 * t357) * t378 + t180;
t169 = t178 * t384 - t179 * t348 + t182 * t383;
t238 = Ifges(7,5) * t276 + Ifges(7,6) * t275 + Ifges(7,3) * t317;
t196 = -mrSges(7,1) * t213 + mrSges(7,3) * t208 + Ifges(7,4) * t227 + Ifges(7,2) * t226 + Ifges(7,6) * t292 - t238 * t276 + t240 * t317;
t197 = mrSges(7,2) * t213 - mrSges(7,3) * t207 + Ifges(7,1) * t227 + Ifges(7,4) * t226 + Ifges(7,5) * t292 + t238 * t275 - t239 * t317;
t364 = -mrSges(6,1) * t219 + mrSges(6,2) * t216 - pkin(5) * t209 - pkin(11) * t194 - t355 * t196 - t350 * t197;
t174 = mrSges(5,1) * t368 + mrSges(5,3) * t221 - pkin(4) * t201 + (t263 + t262) * t319 + t382 * t307 + (Ifges(5,6) - Ifges(6,6)) * t297 + (Ifges(5,4) - Ifges(6,5)) * t253 + (-Ifges(5,2) - Ifges(6,3)) * t252 + t364;
t365 = mrSges(6,2) * t218 - mrSges(6,3) * t219 + Ifges(6,1) * t253 + Ifges(6,4) * t297 + Ifges(6,5) * t252 - pkin(11) * t193 - t196 * t350 + t355 * t197 + t319 * t258;
t175 = -mrSges(5,2) * t368 - mrSges(5,3) * t220 + Ifges(5,1) * t253 - Ifges(5,4) * t252 + Ifges(5,5) * t297 - qJ(5) * t201 - t319 * t261 + t382 * t306 + t365;
t293 = Ifges(4,5) * t322 + Ifges(4,6) * t321 + Ifges(4,3) * t339;
t294 = Ifges(4,4) * t322 + Ifges(4,2) * t321 + Ifges(4,6) * t339;
t165 = mrSges(4,2) * t270 - mrSges(4,3) * t236 + Ifges(4,1) * t300 + Ifges(4,4) * t299 + Ifges(4,5) * t325 - pkin(10) * t186 - t351 * t174 + t391 * t175 + t321 * t293 - t339 * t294;
t295 = Ifges(4,1) * t322 + Ifges(4,4) * t321 + Ifges(4,5) * t339;
t170 = -mrSges(4,1) * t270 + mrSges(4,3) * t237 + Ifges(4,4) * t300 + Ifges(4,2) * t299 + Ifges(4,6) * t325 - pkin(3) * t186 - t322 * t293 + t339 * t295 - t393;
t311 = Ifges(3,6) * t344 + (Ifges(3,4) * t353 + Ifges(3,2) * t357) * t378;
t312 = Ifges(3,5) * t344 + (Ifges(3,1) * t353 + Ifges(3,4) * t357) * t378;
t160 = Ifges(3,5) * t332 - Ifges(3,6) * t333 + Ifges(3,3) * t343 + mrSges(3,1) * t301 - mrSges(3,2) * t302 + t352 * t165 + t356 * t170 + pkin(2) * t363 + pkin(9) * t373 + (t311 * t353 - t312 * t357) * t378;
t310 = Ifges(3,3) * t344 + (Ifges(3,5) * t353 + Ifges(3,6) * t357) * t378;
t162 = mrSges(3,2) * t313 - mrSges(3,3) * t301 + Ifges(3,1) * t332 - Ifges(3,4) * t333 + Ifges(3,5) * t343 - pkin(9) * t180 + t165 * t356 - t170 * t352 + t310 * t374 - t311 * t344;
t361 = mrSges(4,1) * t236 - mrSges(4,2) * t237 + Ifges(4,5) * t300 + Ifges(4,6) * t299 + Ifges(4,3) * t325 + pkin(3) * t200 + pkin(10) * t187 + t391 * t174 + t351 * t175 + t322 * t294 - t321 * t295;
t164 = -mrSges(3,1) * t313 + mrSges(3,3) * t302 + Ifges(3,4) * t332 - Ifges(3,2) * t333 + Ifges(3,6) * t343 - pkin(2) * t180 - t310 * t375 + t344 * t312 - t361;
t366 = mrSges(2,1) * t340 - mrSges(2,2) * t341 + Ifges(2,3) * qJDD(1) + pkin(1) * t169 + t349 * t160 + t162 * t386 + t164 * t385 + t173 * t390;
t171 = m(2) * t341 - t359 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t173;
t168 = t179 * t349 + (t178 * t353 + t182 * t357) * t348;
t166 = m(2) * t340 + qJDD(1) * mrSges(2,1) - t359 * mrSges(2,2) + t169;
t158 = -mrSges(2,2) * g(3) - mrSges(2,3) * t340 + Ifges(2,5) * qJDD(1) - t359 * Ifges(2,6) + t162 * t357 - t164 * t353 + (-t168 * t348 - t169 * t349) * pkin(8);
t157 = mrSges(2,1) * g(3) + mrSges(2,3) * t341 + t359 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t168 - t160 * t348 + (pkin(8) * t173 + t162 * t353 + t164 * t357) * t349;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t358 * t158 - t354 * t157 - pkin(7) * (t166 * t358 + t171 * t354), t158, t162, t165, t175, -t306 * t260 + t365, t197; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t354 * t158 + t358 * t157 + pkin(7) * (-t354 * t166 + t358 * t171), t157, t164, t170, t174, -t307 * t258 - t362, t196; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t366, t366, t160, t361, t393, Ifges(6,5) * t253 + Ifges(6,6) * t297 + Ifges(6,3) * t252 + t307 * t260 - t319 * t262 - t364, t369;];
m_new  = t1;
