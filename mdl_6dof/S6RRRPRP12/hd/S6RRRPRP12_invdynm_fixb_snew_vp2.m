% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPRP12
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
% Datum: 2019-05-07 09:38
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPRP12_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP12_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP12_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP12_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP12_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP12_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP12_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 09:23:34
% EndTime: 2019-05-07 09:24:16
% DurationCPUTime: 14.01s
% Computational Cost: add. (240783->399), mult. (515186->479), div. (0->0), fcn. (391337->10), ass. (0->148)
t341 = sin(pkin(6));
t345 = sin(qJ(2));
t347 = cos(qJ(2));
t372 = qJD(1) * qJD(2);
t325 = (-qJDD(1) * t347 + t345 * t372) * t341;
t375 = qJD(1) * t341;
t323 = (-pkin(2) * t347 - pkin(9) * t345) * t375;
t342 = cos(pkin(6));
t337 = t342 * qJD(1) + qJD(2);
t335 = t337 ^ 2;
t336 = t342 * qJDD(1) + qJDD(2);
t374 = qJD(1) * t347;
t346 = sin(qJ(1));
t348 = cos(qJ(1));
t333 = t346 * g(1) - t348 * g(2);
t349 = qJD(1) ^ 2;
t389 = pkin(8) * t341;
t320 = qJDD(1) * pkin(1) + t349 * t389 + t333;
t334 = -t348 * g(1) - t346 * g(2);
t321 = -t349 * pkin(1) + qJDD(1) * t389 + t334;
t382 = t342 * t345;
t376 = t320 * t382 + t347 * t321;
t247 = -t335 * pkin(2) + t336 * pkin(9) + (-g(3) * t345 + t323 * t374) * t341 + t376;
t324 = (qJDD(1) * t345 + t347 * t372) * t341;
t388 = t342 * g(3);
t248 = t325 * pkin(2) - t324 * pkin(9) - t388 + (-t320 + (pkin(2) * t345 - pkin(9) * t347) * t337 * qJD(1)) * t341;
t344 = sin(qJ(3));
t391 = cos(qJ(3));
t221 = t391 * t247 + t344 * t248;
t370 = t345 * t375;
t310 = -t337 * t391 + t344 * t370;
t311 = t344 * t337 + t370 * t391;
t286 = t310 * pkin(3) - t311 * qJ(4);
t317 = qJDD(3) + t325;
t369 = t341 * t374;
t331 = -qJD(3) + t369;
t330 = t331 ^ 2;
t216 = t330 * pkin(3) - t317 * qJ(4) + 0.2e1 * qJD(4) * t331 + t310 * t286 - t221;
t220 = -t344 * t247 + t248 * t391;
t271 = Ifges(4,4) * t311 - Ifges(4,2) * t310 - Ifges(4,6) * t331;
t282 = t311 * qJD(3) + t344 * t324 - t336 * t391;
t283 = -t310 * qJD(3) + t324 * t391 + t344 * t336;
t288 = -t310 * mrSges(5,2) - t311 * mrSges(5,3);
t295 = t310 * mrSges(5,1) + t331 * mrSges(5,3);
t296 = t311 * mrSges(5,1) - t331 * mrSges(5,2);
t297 = t311 * pkin(4) + t331 * pkin(10);
t309 = t310 ^ 2;
t214 = -t282 * pkin(4) - t309 * pkin(10) - t331 * t297 - t216;
t343 = sin(qJ(5));
t390 = cos(qJ(5));
t292 = t343 * t310 - t331 * t390;
t233 = t292 * qJD(5) - t282 * t390 + t343 * t317;
t291 = -t310 * t390 - t343 * t331;
t234 = -t291 * qJD(5) + t343 * t282 + t317 * t390;
t308 = qJD(5) + t311;
t209 = -0.2e1 * qJD(6) * t292 + (t291 * t308 - t234) * qJ(6) + (t292 * t308 + t233) * pkin(5) + t214;
t258 = -t291 * mrSges(7,2) + t308 * mrSges(7,3);
t261 = -t308 * mrSges(7,1) + t292 * mrSges(7,2);
t199 = m(7) * t209 + t233 * mrSges(7,1) - t234 * mrSges(7,3) + t291 * t258 - t292 * t261;
t259 = -t308 * mrSges(6,2) - t291 * mrSges(6,3);
t260 = t308 * mrSges(6,1) - t292 * mrSges(6,3);
t356 = m(6) * t214 + t233 * mrSges(6,1) + t234 * mrSges(6,2) + t291 * t259 + t292 * t260 + t199;
t353 = -m(5) * t216 + t317 * mrSges(5,3) - t331 * t296 + t356;
t218 = -t317 * pkin(3) - t330 * qJ(4) + t311 * t286 + qJDD(4) - t220;
t385 = t310 * t331;
t211 = (t310 * t311 - t317) * pkin(10) + (t283 - t385) * pkin(4) + t218;
t381 = t342 * t347;
t383 = t341 * t347;
t284 = -g(3) * t383 + t320 * t381 - t345 * t321;
t246 = -t336 * pkin(2) - t335 * pkin(9) + t323 * t370 - t284;
t354 = (-t283 - t385) * qJ(4) + t246 + (-t331 * pkin(3) - 0.2e1 * qJD(4)) * t311;
t215 = -t309 * pkin(4) - t311 * t297 + (pkin(3) + pkin(10)) * t282 + t354;
t207 = t343 * t211 + t390 * t215;
t239 = Ifges(7,1) * t292 + Ifges(7,4) * t308 + Ifges(7,5) * t291;
t240 = Ifges(6,1) * t292 - Ifges(6,4) * t291 + Ifges(6,5) * t308;
t279 = qJDD(5) + t283;
t251 = t291 * pkin(5) - t292 * qJ(6);
t307 = t308 ^ 2;
t202 = -t307 * pkin(5) + t279 * qJ(6) + 0.2e1 * qJD(6) * t308 - t291 * t251 + t207;
t366 = -mrSges(7,1) * t209 + mrSges(7,2) * t202;
t237 = Ifges(7,4) * t292 + Ifges(7,2) * t308 + Ifges(7,6) * t291;
t380 = -Ifges(6,5) * t292 + Ifges(6,6) * t291 - Ifges(6,3) * t308 - t237;
t181 = -mrSges(6,1) * t214 + mrSges(6,3) * t207 - pkin(5) * t199 + (t239 + t240) * t308 + t380 * t292 + (Ifges(6,6) - Ifges(7,6)) * t279 + (Ifges(6,4) - Ifges(7,5)) * t234 + (-Ifges(6,2) - Ifges(7,3)) * t233 + t366;
t206 = t211 * t390 - t343 * t215;
t238 = Ifges(6,4) * t292 - Ifges(6,2) * t291 + Ifges(6,6) * t308;
t204 = -t279 * pkin(5) - t307 * qJ(6) + t292 * t251 + qJDD(6) - t206;
t235 = Ifges(7,5) * t292 + Ifges(7,6) * t308 + Ifges(7,3) * t291;
t362 = mrSges(7,2) * t204 - mrSges(7,3) * t209 + Ifges(7,1) * t234 + Ifges(7,4) * t279 + Ifges(7,5) * t233 + t308 * t235;
t183 = mrSges(6,2) * t214 - mrSges(6,3) * t206 + Ifges(6,1) * t234 - Ifges(6,4) * t233 + Ifges(6,5) * t279 - qJ(6) * t199 - t308 * t238 + t291 * t380 + t362;
t371 = m(7) * t202 + t279 * mrSges(7,3) + t308 * t261;
t252 = t291 * mrSges(7,1) - t292 * mrSges(7,3);
t379 = -t291 * mrSges(6,1) - t292 * mrSges(6,2) - t252;
t387 = -mrSges(6,3) - mrSges(7,2);
t190 = m(6) * t207 - t279 * mrSges(6,2) + t233 * t387 - t308 * t260 + t291 * t379 + t371;
t367 = -m(7) * t204 + t279 * mrSges(7,1) + t308 * t258;
t192 = m(6) * t206 + t279 * mrSges(6,1) + t234 * t387 + t308 * t259 + t292 * t379 + t367;
t184 = t343 * t190 + t192 * t390;
t268 = -Ifges(5,5) * t331 - Ifges(5,6) * t311 + Ifges(5,3) * t310;
t358 = -mrSges(5,2) * t218 + mrSges(5,3) * t216 - Ifges(5,1) * t317 + Ifges(5,4) * t283 - Ifges(5,5) * t282 + pkin(10) * t184 + t343 * t181 - t390 * t183 + t311 * t268;
t359 = -m(5) * t218 - t283 * mrSges(5,1) - t311 * t288 - t184;
t270 = -Ifges(5,4) * t331 - Ifges(5,2) * t311 + Ifges(5,6) * t310;
t377 = -Ifges(4,1) * t311 + Ifges(4,4) * t310 + Ifges(4,5) * t331 + t270;
t392 = -t310 * t377 + mrSges(4,1) * t220 - mrSges(4,2) * t221 + Ifges(4,5) * t283 - Ifges(4,6) * t282 + Ifges(4,3) * t317 + pkin(3) * (-t317 * mrSges(5,2) + t331 * t295 + t359) + qJ(4) * (-t282 * mrSges(5,1) - t310 * t288 + t353) + t311 * t271 - t358;
t386 = -Ifges(5,6) - Ifges(4,4);
t384 = t341 * t345;
t287 = t310 * mrSges(4,1) + t311 * mrSges(4,2);
t293 = t331 * mrSges(4,2) - t310 * mrSges(4,3);
t178 = m(4) * t220 - t283 * mrSges(4,3) - t311 * t287 + (-t293 + t295) * t331 + (mrSges(4,1) - mrSges(5,2)) * t317 + t359;
t294 = -t331 * mrSges(4,1) - t311 * mrSges(4,3);
t188 = (-t287 - t288) * t310 + (-mrSges(4,3) - mrSges(5,1)) * t282 + t331 * t294 - t317 * mrSges(4,2) + m(4) * t221 + t353;
t174 = t178 * t391 + t344 * t188;
t185 = t390 * t190 - t343 * t192;
t272 = -Ifges(5,1) * t331 - Ifges(5,4) * t311 + Ifges(5,5) * t310;
t378 = -Ifges(4,5) * t311 + Ifges(4,6) * t310 + Ifges(4,3) * t331 - t272;
t285 = -g(3) * t384 + t376;
t318 = t337 * mrSges(3,1) - mrSges(3,3) * t370;
t322 = (-mrSges(3,1) * t347 + mrSges(3,2) * t345) * t375;
t368 = -t344 * t178 + t188 * t391;
t172 = m(3) * t285 - t336 * mrSges(3,2) - t325 * mrSges(3,3) - t337 * t318 + t322 * t369 + t368;
t319 = -t337 * mrSges(3,2) + mrSges(3,3) * t369;
t219 = t282 * pkin(3) + t354;
t365 = -m(5) * t219 + t282 * mrSges(5,2) + t310 * t295 - t185;
t355 = -m(4) * t246 - t282 * mrSges(4,1) - t310 * t293 + (-t294 + t296) * t311 + (-mrSges(4,2) + mrSges(5,3)) * t283 + t365;
t176 = m(3) * t284 + t336 * mrSges(3,1) - t324 * mrSges(3,3) + t337 * t319 - t322 * t370 + t355;
t169 = t347 * t172 - t345 * t176;
t302 = -t341 * t320 - t388;
t173 = m(3) * t302 + t325 * mrSges(3,1) + t324 * mrSges(3,2) + (t318 * t345 - t319 * t347) * t375 + t174;
t165 = t172 * t382 - t341 * t173 + t176 * t381;
t179 = -t283 * mrSges(5,3) - t311 * t296 - t365;
t357 = mrSges(5,1) * t216 - mrSges(5,2) * t219 - pkin(4) * t356 + pkin(10) * t185 + t181 * t390 + t343 * t183;
t161 = -mrSges(4,1) * t246 + mrSges(4,3) * t221 - pkin(3) * t179 + t377 * t331 + (Ifges(4,6) - Ifges(5,5)) * t317 + t378 * t311 - t386 * t283 + (-Ifges(4,2) - Ifges(5,3)) * t282 - t357;
t360 = mrSges(7,1) * t204 - mrSges(7,3) * t202 - Ifges(7,4) * t234 - Ifges(7,2) * t279 - Ifges(7,6) * t233 + t292 * t235 - t291 * t239;
t352 = mrSges(6,2) * t207 - t291 * t240 - qJ(6) * (-t233 * mrSges(7,2) - t291 * t252 + t371) - pkin(5) * (-t234 * mrSges(7,2) - t292 * t252 + t367) - mrSges(6,1) * t206 - t292 * t238 + Ifges(6,6) * t233 - Ifges(6,5) * t234 - Ifges(6,3) * t279 + t360;
t351 = -mrSges(5,1) * t218 + mrSges(5,3) * t219 - pkin(4) * t184 + t352;
t166 = -t351 + mrSges(4,2) * t246 - mrSges(4,3) * t220 - qJ(4) * t179 + (-t268 + t271) * t331 + (-Ifges(5,4) + Ifges(4,5)) * t317 + t378 * t310 + (Ifges(4,1) + Ifges(5,2)) * t283 + t386 * t282;
t300 = Ifges(3,6) * t337 + (Ifges(3,4) * t345 + Ifges(3,2) * t347) * t375;
t301 = Ifges(3,5) * t337 + (Ifges(3,1) * t345 + Ifges(3,4) * t347) * t375;
t156 = Ifges(3,5) * t324 - Ifges(3,6) * t325 + Ifges(3,3) * t336 + mrSges(3,1) * t284 - mrSges(3,2) * t285 + t344 * t166 + t391 * t161 + pkin(2) * t355 + pkin(9) * t368 + (t300 * t345 - t301 * t347) * t375;
t299 = Ifges(3,3) * t337 + (Ifges(3,5) * t345 + Ifges(3,6) * t347) * t375;
t158 = mrSges(3,2) * t302 - mrSges(3,3) * t284 + Ifges(3,1) * t324 - Ifges(3,4) * t325 + Ifges(3,5) * t336 - pkin(9) * t174 - t344 * t161 + t166 * t391 + t299 * t369 - t337 * t300;
t160 = -mrSges(3,1) * t302 + mrSges(3,3) * t285 + Ifges(3,4) * t324 - Ifges(3,2) * t325 + Ifges(3,6) * t336 - pkin(2) * t174 - t299 * t370 + t337 * t301 - t392;
t361 = mrSges(2,1) * t333 - mrSges(2,2) * t334 + Ifges(2,3) * qJDD(1) + pkin(1) * t165 + t342 * t156 + t158 * t384 + t160 * t383 + t169 * t389;
t167 = m(2) * t334 - t349 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t169;
t164 = t342 * t173 + (t172 * t345 + t176 * t347) * t341;
t162 = m(2) * t333 + qJDD(1) * mrSges(2,1) - t349 * mrSges(2,2) + t165;
t154 = -mrSges(2,2) * g(3) - mrSges(2,3) * t333 + Ifges(2,5) * qJDD(1) - t349 * Ifges(2,6) + t347 * t158 - t345 * t160 + (-t164 * t341 - t165 * t342) * pkin(8);
t153 = mrSges(2,1) * g(3) + mrSges(2,3) * t334 + t349 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t164 - t341 * t156 + (pkin(8) * t169 + t158 * t345 + t160 * t347) * t342;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t348 * t154 - t346 * t153 - pkin(7) * (t348 * t162 + t346 * t167), t154, t158, t166, -t310 * t270 - t358, t183, -t291 * t237 + t362; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t346 * t154 + t348 * t153 + pkin(7) * (-t346 * t162 + t348 * t167), t153, t160, t161, Ifges(5,4) * t317 - Ifges(5,2) * t283 + Ifges(5,6) * t282 + t331 * t268 + t310 * t272 + t351, t181, -t360; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t361, t361, t156, t392, Ifges(5,5) * t317 - Ifges(5,6) * t283 + Ifges(5,3) * t282 - t331 * t270 + t311 * t272 + t357, -t352, Ifges(7,5) * t234 + Ifges(7,6) * t279 + Ifges(7,3) * t233 + t292 * t237 - t308 * t239 - t366;];
m_new  = t1;
