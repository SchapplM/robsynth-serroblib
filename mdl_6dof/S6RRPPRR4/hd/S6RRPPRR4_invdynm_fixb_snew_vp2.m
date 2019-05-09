% Calculate vector of cutting torques with Newton-Euler for
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-05-06 10:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPPRR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 10:19:51
% EndTime: 2019-05-06 10:20:34
% DurationCPUTime: 24.35s
% Computational Cost: add. (386308->399), mult. (1039485->498), div. (0->0), fcn. (795624->12), ass. (0->162)
t401 = -2 * qJD(4);
t348 = sin(pkin(11));
t357 = cos(qJ(2));
t349 = sin(pkin(6));
t385 = qJD(1) * t349;
t378 = t357 * t385;
t353 = sin(qJ(2));
t379 = t353 * t385;
t396 = cos(pkin(11));
t320 = t348 * t379 - t396 * t378;
t321 = (t348 * t357 + t396 * t353) * t385;
t282 = pkin(3) * t320 - qJ(4) * t321;
t350 = cos(pkin(6));
t343 = qJD(1) * t350 + qJD(2);
t341 = t343 ^ 2;
t342 = qJDD(1) * t350 + qJDD(2);
t381 = qJD(1) * qJD(2);
t330 = (qJDD(1) * t353 + t357 * t381) * t349;
t354 = sin(qJ(1));
t358 = cos(qJ(1));
t339 = t354 * g(1) - g(2) * t358;
t359 = qJD(1) ^ 2;
t399 = pkin(8) * t349;
t327 = qJDD(1) * pkin(1) + t359 * t399 + t339;
t340 = -g(1) * t358 - g(2) * t354;
t328 = -pkin(1) * t359 + qJDD(1) * t399 + t340;
t390 = t350 * t357;
t375 = t327 * t390 - t353 * t328;
t394 = t349 ^ 2 * t359;
t245 = t342 * pkin(2) - t330 * qJ(3) + (pkin(2) * t353 * t394 + (qJ(3) * qJD(1) * t343 - g(3)) * t349) * t357 + t375;
t391 = t350 * t353;
t393 = t349 * t353;
t281 = -g(3) * t393 + t327 * t391 + t357 * t328;
t324 = pkin(2) * t343 - qJ(3) * t379;
t331 = (qJDD(1) * t357 - t353 * t381) * t349;
t380 = t357 ^ 2 * t394;
t248 = -pkin(2) * t380 + qJ(3) * t331 - t324 * t343 + t281;
t389 = t348 * t245 + t396 * t248;
t400 = t341 * pkin(3) - t342 * qJ(4) + t320 * t282 + t343 * t401 - t389;
t231 = -0.2e1 * qJD(3) * t321 + t396 * t245 - t348 * t248;
t398 = mrSges(4,1) - mrSges(5,2);
t397 = -Ifges(5,6) - Ifges(4,4);
t395 = t320 * t343;
t392 = t349 * t357;
t283 = mrSges(4,1) * t320 + mrSges(4,2) * t321;
t295 = t396 * t330 + t348 * t331;
t227 = -t342 * pkin(3) - t341 * qJ(4) + t321 * t282 + qJDD(4) - t231;
t220 = (t320 * t321 - t342) * pkin(9) + (t295 + t395) * pkin(4) + t227;
t294 = t348 * t330 - t396 * t331;
t305 = pkin(4) * t321 - pkin(9) * t343;
t319 = t320 ^ 2;
t310 = -t350 * g(3) - t349 * t327;
t260 = -t331 * pkin(2) - qJ(3) * t380 + t324 * t379 + qJDD(3) + t310;
t361 = (-t295 + t395) * qJ(4) + t260 + (pkin(3) * t343 + t401) * t321;
t224 = -t319 * pkin(4) - t321 * t305 + (pkin(3) + pkin(9)) * t294 + t361;
t352 = sin(qJ(5));
t356 = cos(qJ(5));
t217 = t352 * t220 + t356 * t224;
t300 = t320 * t352 + t343 * t356;
t257 = -qJD(5) * t300 + t294 * t356 - t342 * t352;
t299 = t320 * t356 - t343 * t352;
t263 = -mrSges(6,1) * t299 + mrSges(6,2) * t300;
t318 = qJD(5) + t321;
t271 = mrSges(6,1) * t318 - mrSges(6,3) * t300;
t293 = qJDD(5) + t295;
t264 = -pkin(5) * t299 - pkin(10) * t300;
t317 = t318 ^ 2;
t214 = -pkin(5) * t317 + pkin(10) * t293 + t264 * t299 + t217;
t384 = qJD(3) * t320;
t313 = -0.2e1 * t384;
t222 = -t294 * pkin(4) - t319 * pkin(9) + t343 * t305 + t313 - t400;
t258 = qJD(5) * t299 + t294 * t352 + t342 * t356;
t218 = (-t299 * t318 - t258) * pkin(10) + (t300 * t318 - t257) * pkin(5) + t222;
t351 = sin(qJ(6));
t355 = cos(qJ(6));
t211 = -t214 * t351 + t218 * t355;
t268 = -t300 * t351 + t318 * t355;
t235 = qJD(6) * t268 + t258 * t355 + t293 * t351;
t269 = t300 * t355 + t318 * t351;
t241 = -mrSges(7,1) * t268 + mrSges(7,2) * t269;
t298 = qJD(6) - t299;
t246 = -mrSges(7,2) * t298 + mrSges(7,3) * t268;
t256 = qJDD(6) - t257;
t207 = m(7) * t211 + mrSges(7,1) * t256 - mrSges(7,3) * t235 - t241 * t269 + t246 * t298;
t212 = t214 * t355 + t218 * t351;
t234 = -qJD(6) * t269 - t258 * t351 + t293 * t355;
t247 = mrSges(7,1) * t298 - mrSges(7,3) * t269;
t208 = m(7) * t212 - mrSges(7,2) * t256 + mrSges(7,3) * t234 + t241 * t268 - t247 * t298;
t376 = -t207 * t351 + t355 * t208;
t194 = m(6) * t217 - mrSges(6,2) * t293 + mrSges(6,3) * t257 + t263 * t299 - t271 * t318 + t376;
t216 = t220 * t356 - t224 * t352;
t270 = -mrSges(6,2) * t318 + mrSges(6,3) * t299;
t213 = -pkin(5) * t293 - pkin(10) * t317 + t264 * t300 - t216;
t371 = -m(7) * t213 + t234 * mrSges(7,1) - mrSges(7,2) * t235 + t268 * t246 - t247 * t269;
t203 = m(6) * t216 + mrSges(6,1) * t293 - mrSges(6,3) * t258 - t263 * t300 + t270 * t318 + t371;
t187 = t352 * t194 + t356 * t203;
t284 = -mrSges(5,2) * t320 - mrSges(5,3) * t321;
t370 = -m(5) * t227 - t295 * mrSges(5,1) - t321 * t284 - t187;
t303 = mrSges(5,1) * t320 - mrSges(5,3) * t343;
t386 = -mrSges(4,2) * t343 - mrSges(4,3) * t320 - t303;
t183 = m(4) * t231 - t295 * mrSges(4,3) - t321 * t283 + t398 * t342 + t386 * t343 + t370;
t232 = t313 + t389;
t302 = mrSges(4,1) * t343 - mrSges(4,3) * t321;
t197 = t355 * t207 + t351 * t208;
t195 = -m(6) * t222 + t257 * mrSges(6,1) - t258 * mrSges(6,2) + t299 * t270 - t300 * t271 - t197;
t225 = 0.2e1 * t384 + t400;
t304 = mrSges(5,1) * t321 + mrSges(5,2) * t343;
t366 = -m(5) * t225 + t342 * mrSges(5,3) + t343 * t304 - t195;
t191 = (-mrSges(4,3) - mrSges(5,1)) * t294 + (-t283 - t284) * t320 + t366 - t343 * t302 - t342 * mrSges(4,2) + m(4) * t232;
t178 = t396 * t183 + t348 * t191;
t276 = Ifges(5,1) * t343 - Ifges(5,4) * t321 + Ifges(5,5) * t320;
t388 = -Ifges(4,5) * t321 + Ifges(4,6) * t320 - Ifges(4,3) * t343 - t276;
t274 = Ifges(5,4) * t343 - Ifges(5,2) * t321 + Ifges(5,6) * t320;
t387 = Ifges(4,1) * t321 - Ifges(4,4) * t320 + Ifges(4,5) * t343 - t274;
t280 = -g(3) * t392 + t375;
t326 = -mrSges(3,2) * t343 + mrSges(3,3) * t378;
t329 = (-mrSges(3,1) * t357 + mrSges(3,2) * t353) * t385;
t176 = m(3) * t280 + mrSges(3,1) * t342 - mrSges(3,3) * t330 + t326 * t343 - t329 * t379 + t178;
t325 = mrSges(3,1) * t343 - mrSges(3,3) * t379;
t377 = -t183 * t348 + t396 * t191;
t177 = m(3) * t281 - mrSges(3,2) * t342 + mrSges(3,3) * t331 - t325 * t343 + t329 * t378 + t377;
t172 = -t176 * t353 + t357 * t177;
t188 = t356 * t194 - t352 * t203;
t229 = t294 * pkin(3) + t361;
t372 = m(5) * t229 - t295 * mrSges(5,3) - t321 * t304 + t188;
t363 = m(4) * t260 + t295 * mrSges(4,2) + t398 * t294 + t321 * t302 + t386 * t320 + t372;
t184 = (t325 * t353 - t326 * t357) * t385 + t330 * mrSges(3,2) - t331 * mrSges(3,1) + m(3) * t310 + t363;
t168 = t176 * t390 + t177 * t391 - t184 * t349;
t186 = -t294 * mrSges(5,2) - t320 * t303 + t372;
t236 = Ifges(7,5) * t269 + Ifges(7,6) * t268 + Ifges(7,3) * t298;
t238 = Ifges(7,1) * t269 + Ifges(7,4) * t268 + Ifges(7,5) * t298;
t201 = -mrSges(7,1) * t213 + mrSges(7,3) * t212 + Ifges(7,4) * t235 + Ifges(7,2) * t234 + Ifges(7,6) * t256 - t236 * t269 + t238 * t298;
t237 = Ifges(7,4) * t269 + Ifges(7,2) * t268 + Ifges(7,6) * t298;
t202 = mrSges(7,2) * t213 - mrSges(7,3) * t211 + Ifges(7,1) * t235 + Ifges(7,4) * t234 + Ifges(7,5) * t256 + t236 * t268 - t237 * t298;
t249 = Ifges(6,5) * t300 + Ifges(6,6) * t299 + Ifges(6,3) * t318;
t250 = Ifges(6,4) * t300 + Ifges(6,2) * t299 + Ifges(6,6) * t318;
t180 = mrSges(6,2) * t222 - mrSges(6,3) * t216 + Ifges(6,1) * t258 + Ifges(6,4) * t257 + Ifges(6,5) * t293 - pkin(10) * t197 - t201 * t351 + t202 * t355 + t249 * t299 - t250 * t318;
t251 = Ifges(6,1) * t300 + Ifges(6,4) * t299 + Ifges(6,5) * t318;
t364 = mrSges(7,1) * t211 - mrSges(7,2) * t212 + Ifges(7,5) * t235 + Ifges(7,6) * t234 + Ifges(7,3) * t256 + t237 * t269 - t238 * t268;
t181 = -mrSges(6,1) * t222 + mrSges(6,3) * t217 + Ifges(6,4) * t258 + Ifges(6,2) * t257 + Ifges(6,6) * t293 - pkin(5) * t197 - t249 * t300 + t251 * t318 - t364;
t365 = -mrSges(5,1) * t225 + mrSges(5,2) * t229 - pkin(4) * t195 - pkin(9) * t188 - t352 * t180 - t356 * t181;
t164 = -mrSges(4,1) * t260 + mrSges(4,3) * t232 - pkin(3) * t186 + t387 * t343 + (Ifges(4,6) - Ifges(5,5)) * t342 + t388 * t321 - t397 * t295 + (-Ifges(4,2) - Ifges(5,3)) * t294 + t365;
t272 = Ifges(5,5) * t343 - Ifges(5,6) * t321 + Ifges(5,3) * t320;
t275 = Ifges(4,4) * t321 - Ifges(4,2) * t320 + Ifges(4,6) * t343;
t367 = mrSges(6,1) * t216 - mrSges(6,2) * t217 + Ifges(6,5) * t258 + Ifges(6,6) * t257 + Ifges(6,3) * t293 + pkin(5) * t371 + pkin(10) * t376 + t355 * t201 + t351 * t202 + t300 * t250 - t299 * t251;
t362 = mrSges(5,1) * t227 - mrSges(5,3) * t229 + pkin(4) * t187 + t367;
t169 = -qJ(4) * t186 + (-t275 + t272) * t343 + (Ifges(4,5) - Ifges(5,4)) * t342 + t388 * t320 + (Ifges(4,1) + Ifges(5,2)) * t295 + t397 * t294 + mrSges(4,2) * t260 - mrSges(4,3) * t231 + t362;
t307 = Ifges(3,3) * t343 + (Ifges(3,5) * t353 + Ifges(3,6) * t357) * t385;
t309 = Ifges(3,5) * t343 + (Ifges(3,1) * t353 + Ifges(3,4) * t357) * t385;
t159 = -mrSges(3,1) * t310 + mrSges(3,3) * t281 + Ifges(3,4) * t330 + Ifges(3,2) * t331 + Ifges(3,6) * t342 - pkin(2) * t363 + qJ(3) * t377 + t396 * t164 + t348 * t169 - t307 * t379 + t343 * t309;
t308 = Ifges(3,6) * t343 + (Ifges(3,4) * t353 + Ifges(3,2) * t357) * t385;
t161 = mrSges(3,2) * t310 - mrSges(3,3) * t280 + Ifges(3,1) * t330 + Ifges(3,4) * t331 + Ifges(3,5) * t342 - qJ(3) * t178 - t348 * t164 + t396 * t169 + t307 * t378 - t343 * t308;
t368 = mrSges(5,2) * t227 - mrSges(5,3) * t225 + Ifges(5,1) * t342 - Ifges(5,4) * t295 + Ifges(5,5) * t294 - pkin(9) * t187 + t356 * t180 - t352 * t181 - t321 * t272;
t360 = -mrSges(4,2) * t232 + t387 * t320 + pkin(3) * (-t342 * mrSges(5,2) - t343 * t303 + t370) + qJ(4) * (-t294 * mrSges(5,1) - t320 * t284 + t366) + mrSges(4,1) * t231 + t321 * t275 - Ifges(4,6) * t294 + Ifges(4,5) * t295 + Ifges(4,3) * t342 + t368;
t163 = (t308 * t353 - t309 * t357) * t385 + pkin(2) * t178 + t360 + Ifges(3,3) * t342 + Ifges(3,5) * t330 + Ifges(3,6) * t331 + mrSges(3,1) * t280 - mrSges(3,2) * t281;
t369 = mrSges(2,1) * t339 - mrSges(2,2) * t340 + Ifges(2,3) * qJDD(1) + pkin(1) * t168 + t159 * t392 + t161 * t393 + t350 * t163 + t172 * t399;
t170 = m(2) * t340 - mrSges(2,1) * t359 - qJDD(1) * mrSges(2,2) + t172;
t167 = t350 * t184 + (t176 * t357 + t177 * t353) * t349;
t165 = m(2) * t339 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t359 + t168;
t157 = -mrSges(2,2) * g(3) - mrSges(2,3) * t339 + Ifges(2,5) * qJDD(1) - t359 * Ifges(2,6) - t353 * t159 + t357 * t161 + (-t167 * t349 - t168 * t350) * pkin(8);
t156 = mrSges(2,1) * g(3) + mrSges(2,3) * t340 + t359 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t167 - t349 * t163 + (pkin(8) * t172 + t159 * t357 + t161 * t353) * t350;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t358 * t157 - t354 * t156 - pkin(7) * (t165 * t358 + t170 * t354), t157, t161, t169, -t320 * t274 + t368, t180, t202; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t354 * t157 + t358 * t156 + pkin(7) * (-t165 * t354 + t170 * t358), t156, t159, t164, Ifges(5,4) * t342 - Ifges(5,2) * t295 + Ifges(5,6) * t294 - t343 * t272 + t320 * t276 - t362, t181, t201; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t369, t369, t163, t360, Ifges(5,5) * t342 - Ifges(5,6) * t295 + Ifges(5,3) * t294 + t343 * t274 + t321 * t276 - t365, t367, t364;];
m_new  = t1;
