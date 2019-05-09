% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPRP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-05-07 08:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPRP7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP7_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP7_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP7_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP7_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP7_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 08:13:08
% EndTime: 2019-05-07 08:14:44
% DurationCPUTime: 37.77s
% Computational Cost: add. (673492->394), mult. (1483947->500), div. (0->0), fcn. (1184914->12), ass. (0->154)
t339 = sin(pkin(6));
t344 = sin(qJ(2));
t347 = cos(qJ(2));
t369 = qJD(1) * qJD(2);
t325 = (-qJDD(1) * t347 + t344 * t369) * t339;
t372 = qJD(1) * t339;
t323 = (-pkin(2) * t347 - pkin(9) * t344) * t372;
t341 = cos(pkin(6));
t334 = qJD(1) * t341 + qJD(2);
t332 = t334 ^ 2;
t333 = qJDD(1) * t341 + qJDD(2);
t371 = qJD(1) * t347;
t345 = sin(qJ(1));
t348 = cos(qJ(1));
t330 = t345 * g(1) - g(2) * t348;
t349 = qJD(1) ^ 2;
t383 = pkin(8) * t339;
t320 = qJDD(1) * pkin(1) + t349 * t383 + t330;
t331 = -g(1) * t348 - g(2) * t345;
t321 = -pkin(1) * t349 + qJDD(1) * t383 + t331;
t378 = t341 * t344;
t373 = t320 * t378 + t347 * t321;
t278 = -t332 * pkin(2) + t333 * pkin(9) + (-g(3) * t344 + t323 * t371) * t339 + t373;
t324 = (qJDD(1) * t344 + t347 * t369) * t339;
t382 = t341 * g(3);
t279 = t325 * pkin(2) - t324 * pkin(9) - t382 + (-t320 + (pkin(2) * t344 - pkin(9) * t347) * t334 * qJD(1)) * t339;
t343 = sin(qJ(3));
t346 = cos(qJ(3));
t237 = -t343 * t278 + t346 * t279;
t367 = t344 * t372;
t312 = t334 * t346 - t343 * t367;
t291 = qJD(3) * t312 + t324 * t346 + t333 * t343;
t313 = t334 * t343 + t346 * t367;
t317 = qJDD(3) + t325;
t366 = t339 * t371;
t329 = qJD(3) - t366;
t222 = (t312 * t329 - t291) * qJ(4) + (t312 * t313 + t317) * pkin(3) + t237;
t238 = t346 * t278 + t343 * t279;
t290 = -qJD(3) * t313 - t324 * t343 + t333 * t346;
t302 = pkin(3) * t329 - qJ(4) * t313;
t311 = t312 ^ 2;
t225 = -pkin(3) * t311 + qJ(4) * t290 - t302 * t329 + t238;
t338 = sin(pkin(11));
t340 = cos(pkin(11));
t299 = t312 * t338 + t313 * t340;
t217 = -0.2e1 * qJD(4) * t299 + t340 * t222 - t338 * t225;
t298 = t312 * t340 - t313 * t338;
t218 = 0.2e1 * qJD(4) * t298 + t338 * t222 + t340 * t225;
t273 = -pkin(4) * t298 - pkin(10) * t299;
t328 = t329 ^ 2;
t215 = -pkin(4) * t328 + pkin(10) * t317 + t273 * t298 + t218;
t377 = t341 * t347;
t379 = t339 * t347;
t292 = -g(3) * t379 + t320 * t377 - t344 * t321;
t277 = -t333 * pkin(2) - t332 * pkin(9) + t323 * t367 - t292;
t228 = -t290 * pkin(3) - t311 * qJ(4) + t313 * t302 + qJDD(4) + t277;
t266 = t290 * t340 - t291 * t338;
t267 = t290 * t338 + t291 * t340;
t220 = (-t298 * t329 - t267) * pkin(10) + (t299 * t329 - t266) * pkin(4) + t228;
t342 = sin(qJ(5));
t384 = cos(qJ(5));
t211 = -t342 * t215 + t220 * t384;
t212 = t215 * t384 + t342 * t220;
t281 = t299 * t384 + t342 * t329;
t235 = qJD(5) * t281 + t267 * t342 - t317 * t384;
t280 = t299 * t342 - t329 * t384;
t236 = -t280 * qJD(5) + t267 * t384 + t342 * t317;
t297 = qJD(5) - t298;
t239 = Ifges(7,5) * t281 + Ifges(7,6) * t297 + Ifges(7,3) * t280;
t242 = Ifges(6,4) * t281 - Ifges(6,2) * t280 + Ifges(6,6) * t297;
t244 = Ifges(6,1) * t281 - Ifges(6,4) * t280 + Ifges(6,5) * t297;
t252 = mrSges(7,1) * t280 - mrSges(7,3) * t281;
t265 = qJDD(5) - t266;
t251 = pkin(5) * t280 - qJ(6) * t281;
t296 = t297 ^ 2;
t207 = -pkin(5) * t296 + qJ(6) * t265 + 0.2e1 * qJD(6) * t297 - t251 * t280 + t212;
t209 = -t265 * pkin(5) - t296 * qJ(6) + t281 * t251 + qJDD(6) - t211;
t243 = Ifges(7,1) * t281 + Ifges(7,4) * t297 + Ifges(7,5) * t280;
t358 = mrSges(7,1) * t209 - mrSges(7,3) * t207 - Ifges(7,4) * t236 - Ifges(7,2) * t265 - Ifges(7,6) * t235 - t280 * t243;
t254 = -mrSges(7,2) * t280 + mrSges(7,3) * t297;
t361 = -m(7) * t209 + t265 * mrSges(7,1) + t297 * t254;
t257 = -mrSges(7,1) * t297 + mrSges(7,2) * t281;
t368 = m(7) * t207 + t265 * mrSges(7,3) + t297 * t257;
t385 = -(-t242 + t239) * t281 + mrSges(6,1) * t211 - mrSges(6,2) * t212 + Ifges(6,5) * t236 - Ifges(6,6) * t235 + Ifges(6,3) * t265 + pkin(5) * (-t236 * mrSges(7,2) - t281 * t252 + t361) + qJ(6) * (-t235 * mrSges(7,2) - t280 * t252 + t368) + t280 * t244 - t358;
t381 = -mrSges(6,3) - mrSges(7,2);
t380 = t339 * t344;
t272 = -mrSges(5,1) * t298 + mrSges(5,2) * t299;
t283 = mrSges(5,1) * t329 - mrSges(5,3) * t299;
t256 = mrSges(6,1) * t297 - mrSges(6,3) * t281;
t374 = -mrSges(6,1) * t280 - mrSges(6,2) * t281 - t252;
t197 = m(6) * t212 - t265 * mrSges(6,2) + t235 * t381 - t297 * t256 + t280 * t374 + t368;
t255 = -mrSges(6,2) * t297 - mrSges(6,3) * t280;
t199 = m(6) * t211 + t265 * mrSges(6,1) + t236 * t381 + t297 * t255 + t281 * t374 + t361;
t363 = t197 * t384 - t199 * t342;
t185 = m(5) * t218 - mrSges(5,2) * t317 + mrSges(5,3) * t266 + t272 * t298 - t283 * t329 + t363;
t282 = -mrSges(5,2) * t329 + mrSges(5,3) * t298;
t214 = -t317 * pkin(4) - t328 * pkin(10) + t299 * t273 - t217;
t210 = -0.2e1 * qJD(6) * t281 + (t280 * t297 - t236) * qJ(6) + (t281 * t297 + t235) * pkin(5) + t214;
t204 = m(7) * t210 + mrSges(7,1) * t235 - t236 * mrSges(7,3) + t254 * t280 - t281 * t257;
t353 = -m(6) * t214 - t235 * mrSges(6,1) - mrSges(6,2) * t236 - t280 * t255 - t256 * t281 - t204;
t194 = m(5) * t217 + mrSges(5,1) * t317 - mrSges(5,3) * t267 - t272 * t299 + t282 * t329 + t353;
t180 = t338 * t185 + t340 * t194;
t300 = -mrSges(4,1) * t312 + mrSges(4,2) * t313;
t301 = -mrSges(4,2) * t329 + mrSges(4,3) * t312;
t178 = m(4) * t237 + mrSges(4,1) * t317 - mrSges(4,3) * t291 - t300 * t313 + t301 * t329 + t180;
t303 = mrSges(4,1) * t329 - mrSges(4,3) * t313;
t364 = t340 * t185 - t194 * t338;
t179 = m(4) * t238 - mrSges(4,2) * t317 + mrSges(4,3) * t290 + t300 * t312 - t303 * t329 + t364;
t172 = t346 * t178 + t343 * t179;
t191 = t342 * t197 + t199 * t384;
t241 = Ifges(7,4) * t281 + Ifges(7,2) * t297 + Ifges(7,6) * t280;
t376 = -Ifges(6,5) * t281 + Ifges(6,6) * t280 - Ifges(6,3) * t297 - t241;
t293 = -g(3) * t380 + t373;
t318 = mrSges(3,1) * t334 - mrSges(3,3) * t367;
t322 = (-mrSges(3,1) * t347 + mrSges(3,2) * t344) * t372;
t365 = -t178 * t343 + t346 * t179;
t170 = m(3) * t293 - mrSges(3,2) * t333 - mrSges(3,3) * t325 - t318 * t334 + t322 * t366 + t365;
t319 = -mrSges(3,2) * t334 + mrSges(3,3) * t366;
t355 = m(5) * t228 - t266 * mrSges(5,1) + mrSges(5,2) * t267 - t298 * t282 + t283 * t299 + t191;
t352 = -m(4) * t277 + t290 * mrSges(4,1) - mrSges(4,2) * t291 + t312 * t301 - t303 * t313 - t355;
t182 = m(3) * t292 + mrSges(3,1) * t333 - mrSges(3,3) * t324 + t319 * t334 - t322 * t367 + t352;
t167 = t347 * t170 - t182 * t344;
t307 = -t339 * t320 - t382;
t171 = m(3) * t307 + t325 * mrSges(3,1) + t324 * mrSges(3,2) + (t318 * t344 - t319 * t347) * t372 + t172;
t163 = t170 * t378 - t171 * t339 + t182 * t377;
t360 = -mrSges(7,1) * t210 + mrSges(7,2) * t207;
t357 = mrSges(7,2) * t209 - mrSges(7,3) * t210 + Ifges(7,1) * t236 + Ifges(7,4) * t265 + Ifges(7,5) * t235 + t297 * t239;
t187 = -mrSges(6,1) * t214 + mrSges(6,3) * t212 - pkin(5) * t204 + (t243 + t244) * t297 + t376 * t281 + (Ifges(6,6) - Ifges(7,6)) * t265 + (Ifges(6,4) - Ifges(7,5)) * t236 + (-Ifges(6,2) - Ifges(7,3)) * t235 + t360;
t189 = mrSges(6,2) * t214 - mrSges(6,3) * t211 + Ifges(6,1) * t236 - Ifges(6,4) * t235 + Ifges(6,5) * t265 - qJ(6) * t204 - t297 * t242 + t280 * t376 + t357;
t268 = Ifges(5,5) * t299 + Ifges(5,6) * t298 + Ifges(5,3) * t329;
t269 = Ifges(5,4) * t299 + Ifges(5,2) * t298 + Ifges(5,6) * t329;
t173 = mrSges(5,2) * t228 - mrSges(5,3) * t217 + Ifges(5,1) * t267 + Ifges(5,4) * t266 + Ifges(5,5) * t317 - pkin(10) * t191 - t342 * t187 + t189 * t384 + t298 * t268 - t329 * t269;
t270 = Ifges(5,1) * t299 + Ifges(5,4) * t298 + Ifges(5,5) * t329;
t174 = -mrSges(5,1) * t228 + mrSges(5,3) * t218 + Ifges(5,4) * t267 + Ifges(5,2) * t266 + Ifges(5,6) * t317 - pkin(4) * t191 - t299 * t268 + t329 * t270 - t385;
t284 = Ifges(4,5) * t313 + Ifges(4,6) * t312 + Ifges(4,3) * t329;
t286 = Ifges(4,1) * t313 + Ifges(4,4) * t312 + Ifges(4,5) * t329;
t159 = -mrSges(4,1) * t277 + mrSges(4,3) * t238 + Ifges(4,4) * t291 + Ifges(4,2) * t290 + Ifges(4,6) * t317 - pkin(3) * t355 + qJ(4) * t364 + t338 * t173 + t340 * t174 - t313 * t284 + t329 * t286;
t285 = Ifges(4,4) * t313 + Ifges(4,2) * t312 + Ifges(4,6) * t329;
t164 = mrSges(4,2) * t277 - mrSges(4,3) * t237 + Ifges(4,1) * t291 + Ifges(4,4) * t290 + Ifges(4,5) * t317 - qJ(4) * t180 + t173 * t340 - t174 * t338 + t284 * t312 - t285 * t329;
t305 = Ifges(3,6) * t334 + (Ifges(3,4) * t344 + Ifges(3,2) * t347) * t372;
t306 = Ifges(3,5) * t334 + (Ifges(3,1) * t344 + Ifges(3,4) * t347) * t372;
t154 = Ifges(3,5) * t324 - Ifges(3,6) * t325 + Ifges(3,3) * t333 + mrSges(3,1) * t292 - mrSges(3,2) * t293 + t343 * t164 + t346 * t159 + pkin(2) * t352 + pkin(9) * t365 + (t305 * t344 - t306 * t347) * t372;
t304 = Ifges(3,3) * t334 + (Ifges(3,5) * t344 + Ifges(3,6) * t347) * t372;
t156 = mrSges(3,2) * t307 - mrSges(3,3) * t292 + Ifges(3,1) * t324 - Ifges(3,4) * t325 + Ifges(3,5) * t333 - pkin(9) * t172 - t159 * t343 + t164 * t346 + t304 * t366 - t305 * t334;
t354 = -mrSges(5,1) * t217 + mrSges(5,2) * t218 - Ifges(5,5) * t267 - Ifges(5,6) * t266 - Ifges(5,3) * t317 - pkin(4) * t353 - pkin(10) * t363 - t187 * t384 - t342 * t189 - t299 * t269 + t298 * t270;
t350 = mrSges(4,1) * t237 - mrSges(4,2) * t238 + Ifges(4,5) * t291 + Ifges(4,6) * t290 + Ifges(4,3) * t317 + pkin(3) * t180 + t313 * t285 - t312 * t286 - t354;
t158 = -mrSges(3,1) * t307 + mrSges(3,3) * t293 + Ifges(3,4) * t324 - Ifges(3,2) * t325 + Ifges(3,6) * t333 - pkin(2) * t172 - t304 * t367 + t334 * t306 - t350;
t356 = mrSges(2,1) * t330 - mrSges(2,2) * t331 + Ifges(2,3) * qJDD(1) + pkin(1) * t163 + t341 * t154 + t156 * t380 + t158 * t379 + t167 * t383;
t165 = m(2) * t331 - mrSges(2,1) * t349 - qJDD(1) * mrSges(2,2) + t167;
t162 = t341 * t171 + (t170 * t344 + t182 * t347) * t339;
t160 = m(2) * t330 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t349 + t163;
t152 = -mrSges(2,2) * g(3) - mrSges(2,3) * t330 + Ifges(2,5) * qJDD(1) - t349 * Ifges(2,6) + t347 * t156 - t344 * t158 + (-t162 * t339 - t163 * t341) * pkin(8);
t151 = mrSges(2,1) * g(3) + mrSges(2,3) * t331 + t349 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t162 - t339 * t154 + (pkin(8) * t167 + t156 * t344 + t158 * t347) * t341;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t348 * t152 - t345 * t151 - pkin(7) * (t160 * t348 + t165 * t345), t152, t156, t164, t173, t189, -t241 * t280 + t357; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t345 * t152 + t348 * t151 + pkin(7) * (-t160 * t345 + t165 * t348), t151, t158, t159, t174, t187, -t281 * t239 - t358; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t356, t356, t154, t350, -t354, t385, Ifges(7,5) * t236 + Ifges(7,6) * t265 + Ifges(7,3) * t235 + t281 * t241 - t297 * t243 - t360;];
m_new  = t1;
