% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRPP2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
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
% Datum: 2019-05-06 12:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRPP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:23:01
% EndTime: 2019-05-06 12:23:17
% DurationCPUTime: 7.47s
% Computational Cost: add. (113567->388), mult. (256430->451), div. (0->0), fcn. (173719->8), ass. (0->144)
t344 = sin(qJ(2));
t346 = cos(qJ(2));
t374 = qJD(1) * qJD(2);
t328 = qJDD(1) * t344 + t346 * t374;
t345 = sin(qJ(1));
t347 = cos(qJ(1));
t334 = -g(1) * t347 - g(2) * t345;
t349 = qJD(1) ^ 2;
t323 = -pkin(1) * t349 + qJDD(1) * pkin(7) + t334;
t383 = t344 * t323;
t390 = pkin(2) * t349;
t267 = qJDD(2) * pkin(2) - t328 * qJ(3) - t383 + (qJ(3) * t374 + t344 * t390 - g(3)) * t346;
t306 = -t344 * g(3) + t323 * t346;
t377 = qJD(1) * t344;
t330 = qJD(2) * pkin(2) - qJ(3) * t377;
t341 = t346 ^ 2;
t373 = t346 * qJDD(1);
t364 = -t344 * t374 + t373;
t268 = qJ(3) * t364 - qJD(2) * t330 - t341 * t390 + t306;
t342 = sin(pkin(9));
t376 = qJD(1) * t346;
t387 = cos(pkin(9));
t316 = -t342 * t377 + t376 * t387;
t218 = 0.2e1 * qJD(3) * t316 + t267 * t342 + t268 * t387;
t317 = (t342 * t346 + t344 * t387) * qJD(1);
t286 = -mrSges(4,1) * t316 + mrSges(4,2) * t317;
t300 = -t328 * t342 + t364 * t387;
t308 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t317;
t287 = -pkin(3) * t316 - pkin(8) * t317;
t348 = qJD(2) ^ 2;
t213 = -pkin(3) * t348 + qJDD(2) * pkin(8) + t287 * t316 + t218;
t333 = t345 * g(1) - g(2) * t347;
t365 = -qJDD(1) * pkin(1) - t333;
t274 = t330 * t377 - t364 * pkin(2) + qJDD(3) + (-qJ(3) * t341 - pkin(7)) * t349 + t365;
t301 = t328 * t387 + t342 * t364;
t215 = (-qJD(2) * t316 - t301) * pkin(8) + (qJD(2) * t317 - t300) * pkin(3) + t274;
t343 = sin(qJ(4));
t392 = cos(qJ(4));
t208 = -t213 * t343 + t215 * t392;
t303 = -qJD(2) * t392 + t317 * t343;
t259 = -qJD(4) * t303 + qJDD(2) * t343 + t301 * t392;
t315 = qJD(4) - t316;
t276 = mrSges(7,2) * t315 + mrSges(7,3) * t303;
t277 = -mrSges(5,2) * t315 - mrSges(5,3) * t303;
t299 = qJDD(4) - t300;
t304 = qJD(2) * t343 + t317 * t392;
t269 = pkin(4) * t303 - qJ(5) * t304;
t314 = t315 ^ 2;
t206 = -t299 * pkin(4) - t314 * qJ(5) + t269 * t304 + qJDD(5) - t208;
t385 = t303 * t315;
t394 = -0.2e1 * t304;
t194 = qJD(6) * t394 + (-t259 - t385) * qJ(6) + (t303 * t304 - t299) * pkin(5) + t206;
t271 = -mrSges(7,1) * t303 + mrSges(7,2) * t304;
t366 = -m(7) * t194 + mrSges(7,3) * t259 + t271 * t304;
t270 = mrSges(6,1) * t303 - mrSges(6,3) * t304;
t378 = -mrSges(5,1) * t303 - mrSges(5,2) * t304 - t270;
t389 = -mrSges(5,3) - mrSges(6,2);
t275 = -mrSges(6,2) * t303 + mrSges(6,3) * t315;
t397 = -m(6) * t206 + mrSges(6,1) * t299 + t275 * t315;
t183 = m(5) * t208 + (t276 + t277) * t315 + t378 * t304 + (mrSges(5,1) + mrSges(7,1)) * t299 + t389 * t259 + t366 + t397;
t209 = t213 * t392 + t215 * t343;
t258 = qJD(4) * t304 - qJDD(2) * t392 + t301 * t343;
t279 = -mrSges(7,1) * t315 - mrSges(7,3) * t304;
t280 = mrSges(5,1) * t315 - mrSges(5,3) * t304;
t393 = 2 * qJD(5);
t204 = -pkin(4) * t314 + qJ(5) * t299 - t269 * t303 + t315 * t393 + t209;
t278 = -pkin(5) * t315 - qJ(6) * t304;
t302 = t303 ^ 2;
t198 = -pkin(5) * t302 + qJ(6) * t258 + 0.2e1 * qJD(6) * t303 + t278 * t315 + t204;
t372 = m(7) * t198 + mrSges(7,3) * t258 + t271 * t303;
t281 = -mrSges(6,1) * t315 + mrSges(6,2) * t304;
t398 = m(6) * t204 + mrSges(6,3) * t299 + t281 * t315;
t184 = m(5) * t209 + (t279 - t280) * t315 + t378 * t303 + (-mrSges(5,2) + mrSges(7,2)) * t299 + t389 * t258 + t372 + t398;
t368 = -t183 * t343 + t184 * t392;
t174 = m(4) * t218 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t300 - qJD(2) * t308 + t286 * t316 + t368;
t375 = qJD(3) * t317;
t311 = -0.2e1 * t375;
t379 = t267 * t387 - t268 * t342;
t217 = t311 + t379;
t307 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t316;
t360 = qJDD(2) * pkin(3) + t348 * pkin(8) - t287 * t317 + t379;
t399 = (-t259 + t385) * qJ(5);
t201 = -t302 * qJ(6) + qJDD(6) + t311 + (-pkin(4) - pkin(5)) * t258 - t399 + (-pkin(4) * t315 + t278 + t393) * t304 + t360;
t192 = -m(7) * t201 + mrSges(7,1) * t258 - mrSges(7,2) * t259 + t276 * t303 - t279 * t304;
t212 = 0.2e1 * t375 - t360;
t207 = qJD(5) * t394 + t399 + (t304 * t315 + t258) * pkin(4) + t212;
t189 = m(6) * t207 + mrSges(6,1) * t258 - mrSges(6,3) * t259 + t275 * t303 - t281 * t304 + t192;
t352 = -m(5) * t212 - mrSges(5,1) * t258 - mrSges(5,2) * t259 - t277 * t303 - t280 * t304 - t189;
t179 = m(4) * t217 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t301 + qJD(2) * t307 - t286 * t317 + t352;
t167 = t174 * t342 + t179 * t387;
t305 = -t346 * g(3) - t383;
t386 = Ifges(3,6) * qJD(2);
t319 = t386 + (Ifges(3,4) * t344 + Ifges(3,2) * t346) * qJD(1);
t320 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t344 + Ifges(3,4) * t346) * qJD(1);
t229 = Ifges(7,5) * t304 + Ifges(7,6) * t303 - Ifges(7,3) * t315;
t236 = Ifges(6,1) * t304 + Ifges(6,4) * t315 + Ifges(6,5) * t303;
t237 = Ifges(5,1) * t304 - Ifges(5,4) * t303 + Ifges(5,5) * t315;
t191 = t299 * mrSges(7,2) + t315 * t279 + t372;
t235 = Ifges(7,1) * t304 + Ifges(7,4) * t303 - Ifges(7,5) * t315;
t362 = mrSges(7,1) * t201 - mrSges(7,3) * t198 - Ifges(7,4) * t259 - Ifges(7,2) * t258 + Ifges(7,6) * t299 + t235 * t315;
t355 = mrSges(6,1) * t207 - mrSges(6,2) * t204 + pkin(5) * t192 + qJ(6) * t191 - t362;
t233 = Ifges(6,4) * t304 + Ifges(6,2) * t315 + Ifges(6,6) * t303;
t381 = -Ifges(5,5) * t304 + Ifges(5,6) * t303 - Ifges(5,3) * t315 - t233;
t163 = (t237 + t236) * t315 + (t229 + t381) * t304 + (Ifges(5,6) - Ifges(6,6)) * t299 + (Ifges(5,4) - Ifges(6,5)) * t259 + (-Ifges(5,2) - Ifges(6,3)) * t258 + mrSges(5,3) * t209 - mrSges(5,1) * t212 - pkin(4) * t189 - t355;
t232 = Ifges(7,4) * t304 + Ifges(7,2) * t303 - Ifges(7,6) * t315;
t234 = Ifges(5,4) * t304 - Ifges(5,2) * t303 + Ifges(5,6) * t315;
t190 = -t299 * mrSges(7,1) - t315 * t276 - t366;
t230 = Ifges(6,5) * t304 + Ifges(6,6) * t315 + Ifges(6,3) * t303;
t361 = mrSges(7,2) * t201 - mrSges(7,3) * t194 + Ifges(7,1) * t259 + Ifges(7,4) * t258 - Ifges(7,5) * t299 + t229 * t303;
t354 = mrSges(6,2) * t206 - mrSges(6,3) * t207 + Ifges(6,1) * t259 + Ifges(6,4) * t299 + Ifges(6,5) * t258 - qJ(6) * t190 + t230 * t315 + t361;
t169 = (-t234 + t232) * t315 + t381 * t303 + Ifges(5,5) * t299 - Ifges(5,4) * t258 + Ifges(5,1) * t259 - mrSges(5,3) * t208 + mrSges(5,2) * t212 - qJ(5) * t189 + t354;
t283 = Ifges(4,4) * t317 + Ifges(4,2) * t316 + Ifges(4,6) * qJD(2);
t284 = Ifges(4,1) * t317 + Ifges(4,4) * t316 + Ifges(4,5) * qJD(2);
t357 = -mrSges(4,1) * t217 + mrSges(4,2) * t218 - Ifges(4,5) * t301 - Ifges(4,6) * t300 - Ifges(4,3) * qJDD(2) - pkin(3) * t352 - pkin(8) * t368 - t163 * t392 - t169 * t343 - t283 * t317 + t316 * t284;
t400 = -mrSges(3,1) * t305 + mrSges(3,2) * t306 - Ifges(3,5) * t328 - Ifges(3,3) * qJDD(2) - pkin(2) * t167 + (t344 * (-t319 + t386) + t346 * t320) * qJD(1) + t357;
t363 = mrSges(7,1) * t194 - mrSges(7,2) * t198 + Ifges(7,5) * t259 + Ifges(7,6) * t258 - Ifges(7,3) * t299 + t232 * t304 - t235 * t303;
t356 = mrSges(6,1) * t206 - mrSges(6,3) * t204 - Ifges(6,4) * t259 - Ifges(6,2) * t299 - Ifges(6,6) * t258 + pkin(5) * t190 - t236 * t303 + t363;
t396 = (t234 - t230) * t304 + mrSges(5,1) * t208 - mrSges(5,2) * t209 + Ifges(5,5) * t259 - Ifges(5,6) * t258 + Ifges(5,3) * t299 + pkin(4) * (-t259 * mrSges(6,2) - t304 * t270 - t190 + t397) + qJ(5) * (-t258 * mrSges(6,2) - t303 * t270 + t191 + t398) + t303 * t237 - t356;
t388 = Ifges(3,6) * t346;
t384 = t315 * t232;
t176 = t183 * t392 + t184 * t343;
t327 = (-mrSges(3,1) * t346 + mrSges(3,2) * t344) * qJD(1);
t332 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t376;
t165 = m(3) * t305 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t328 + qJD(2) * t332 - t327 * t377 + t167;
t331 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t377;
t369 = t174 * t387 - t342 * t179;
t166 = m(3) * t306 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t364 - qJD(2) * t331 + t327 * t376 + t369;
t370 = -t165 * t344 + t166 * t346;
t282 = Ifges(4,5) * t317 + Ifges(4,6) * t316 + Ifges(4,3) * qJD(2);
t157 = mrSges(4,2) * t274 - mrSges(4,3) * t217 + Ifges(4,1) * t301 + Ifges(4,4) * t300 + Ifges(4,5) * qJDD(2) - pkin(8) * t176 - qJD(2) * t283 - t163 * t343 + t169 * t392 + t282 * t316;
t161 = -mrSges(4,1) * t274 + mrSges(4,3) * t218 + Ifges(4,4) * t301 + Ifges(4,2) * t300 + Ifges(4,6) * qJDD(2) - pkin(3) * t176 + qJD(2) * t284 - t317 * t282 - t396;
t318 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t344 + t388) * qJD(1);
t322 = -t349 * pkin(7) + t365;
t358 = m(4) * t274 - mrSges(4,1) * t300 + mrSges(4,2) * t301 - t307 * t316 + t308 * t317 + t176;
t153 = -mrSges(3,1) * t322 + mrSges(3,3) * t306 + Ifges(3,4) * t328 + Ifges(3,2) * t364 + Ifges(3,6) * qJDD(2) - pkin(2) * t358 + qJ(3) * t369 + qJD(2) * t320 + t342 * t157 + t161 * t387 - t318 * t377;
t156 = mrSges(3,2) * t322 - mrSges(3,3) * t305 + Ifges(3,1) * t328 + Ifges(3,4) * t364 + Ifges(3,5) * qJDD(2) - qJ(3) * t167 - qJD(2) * t319 + t157 * t387 - t342 * t161 + t318 * t376;
t353 = -m(3) * t322 + mrSges(3,1) * t364 - mrSges(3,2) * t328 - t331 * t377 + t332 * t376 - t358;
t359 = mrSges(2,1) * t333 - mrSges(2,2) * t334 + Ifges(2,3) * qJDD(1) + pkin(1) * t353 + pkin(7) * t370 + t153 * t346 + t156 * t344;
t170 = m(2) * t333 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t349 + t353;
t160 = t165 * t346 + t166 * t344;
t158 = m(2) * t334 - mrSges(2,1) * t349 - qJDD(1) * mrSges(2,2) + t370;
t154 = mrSges(2,1) * g(3) + t349 * Ifges(2,5) + mrSges(2,3) * t334 - pkin(1) * t160 + (Ifges(2,6) - t388) * qJDD(1) + t400;
t151 = -mrSges(2,2) * g(3) - mrSges(2,3) * t333 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t349 - pkin(7) * t160 - t153 * t344 + t156 * t346;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t347 * t151 - t345 * t154 - pkin(6) * (t158 * t345 + t170 * t347), t151, t156, t157, t169, -t303 * t233 + t354 + t384, t361 + t384; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t345 * t151 + t347 * t154 + pkin(6) * (t158 * t347 - t170 * t345), t154, t153, t161, t163, -t304 * t230 - t356, -t304 * t229 - t362; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t359, t359, Ifges(3,6) * t373 - t400, -t357, t396, -t315 * t236 + Ifges(6,6) * t299 + Ifges(6,3) * t258 + Ifges(6,5) * t259 + t355 + (-t229 + t233) * t304, t363;];
m_new  = t1;
