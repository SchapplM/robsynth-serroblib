% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-05-06 14:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRPR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR6_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR6_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR6_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR6_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR6_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 14:20:37
% EndTime: 2019-05-06 14:21:32
% DurationCPUTime: 27.16s
% Computational Cost: add. (426050->400), mult. (1119805->498), div. (0->0), fcn. (872130->12), ass. (0->158)
t345 = sin(pkin(6));
t350 = sin(qJ(2));
t353 = cos(qJ(2));
t377 = qJD(1) * qJD(2);
t328 = (qJDD(1) * t350 + t353 * t377) * t345;
t347 = cos(pkin(6));
t338 = qJDD(1) * t347 + qJDD(2);
t339 = qJD(1) * t347 + qJD(2);
t351 = sin(qJ(1));
t354 = cos(qJ(1));
t335 = t351 * g(1) - g(2) * t354;
t355 = qJD(1) ^ 2;
t389 = pkin(8) * t345;
t325 = qJDD(1) * pkin(1) + t355 * t389 + t335;
t336 = -g(1) * t354 - g(2) * t351;
t326 = -pkin(1) * t355 + qJDD(1) * t389 + t336;
t382 = t347 * t353;
t370 = t325 * t382 - t350 * t326;
t386 = t345 ^ 2 * t355;
t248 = t338 * pkin(2) - t328 * qJ(3) + (pkin(2) * t350 * t386 + (qJ(3) * qJD(1) * t339 - g(3)) * t345) * t353 + t370;
t383 = t347 * t350;
t385 = t345 * t350;
t290 = -g(3) * t385 + t325 * t383 + t353 * t326;
t379 = qJD(1) * t345;
t375 = t350 * t379;
t322 = pkin(2) * t339 - qJ(3) * t375;
t329 = (qJDD(1) * t353 - t350 * t377) * t345;
t376 = t353 ^ 2 * t386;
t251 = -pkin(2) * t376 + qJ(3) * t329 - t322 * t339 + t290;
t344 = sin(pkin(11));
t346 = cos(pkin(11));
t320 = (t344 * t353 + t346 * t350) * t379;
t227 = -0.2e1 * qJD(3) * t320 + t346 * t248 - t344 * t251;
t374 = t353 * t379;
t319 = -t344 * t375 + t346 * t374;
t228 = 0.2e1 * qJD(3) * t319 + t344 * t248 + t346 * t251;
t292 = -pkin(3) * t319 - pkin(9) * t320;
t337 = t339 ^ 2;
t225 = -pkin(3) * t337 + pkin(9) * t338 + t292 * t319 + t228;
t311 = -t347 * g(3) - t345 * t325;
t271 = -t329 * pkin(2) - qJ(3) * t376 + t322 * t375 + qJDD(3) + t311;
t299 = -t328 * t344 + t329 * t346;
t300 = t328 * t346 + t329 * t344;
t231 = (-t319 * t339 - t300) * pkin(9) + (t320 * t339 - t299) * pkin(3) + t271;
t349 = sin(qJ(4));
t390 = cos(qJ(4));
t220 = -t349 * t225 + t390 * t231;
t221 = t390 * t225 + t349 * t231;
t304 = t349 * t320 - t390 * t339;
t305 = t390 * t320 + t349 * t339;
t318 = qJD(4) - t319;
t255 = Ifges(5,4) * t305 - Ifges(5,2) * t304 + Ifges(5,6) * t318;
t267 = t305 * qJD(4) + t349 * t300 - t390 * t338;
t268 = -t304 * qJD(4) + t390 * t300 + t349 * t338;
t275 = -mrSges(6,2) * t304 - mrSges(6,3) * t305;
t279 = mrSges(6,1) * t304 - mrSges(6,3) * t318;
t298 = qJDD(4) - t299;
t273 = pkin(4) * t304 - qJ(5) * t305;
t317 = t318 ^ 2;
t218 = -t298 * pkin(4) - t317 * qJ(5) + t305 * t273 + qJDD(5) - t220;
t387 = t304 * t318;
t212 = (t304 * t305 - t298) * pkin(10) + (t268 + t387) * pkin(5) + t218;
t283 = pkin(5) * t305 - pkin(10) * t318;
t303 = t304 ^ 2;
t224 = -t338 * pkin(3) - t337 * pkin(9) + t320 * t292 - t227;
t391 = -2 * qJD(5);
t358 = (-t268 + t387) * qJ(5) + t224 + (pkin(4) * t318 + t391) * t305;
t215 = -t303 * pkin(5) - t305 * t283 + (pkin(4) + pkin(10)) * t267 + t358;
t348 = sin(qJ(6));
t352 = cos(qJ(6));
t209 = t212 * t352 - t215 * t348;
t277 = t304 * t352 - t318 * t348;
t236 = qJD(6) * t277 + t267 * t348 + t298 * t352;
t278 = t304 * t348 + t318 * t352;
t242 = -mrSges(7,1) * t277 + mrSges(7,2) * t278;
t302 = qJD(6) + t305;
t249 = -mrSges(7,2) * t302 + mrSges(7,3) * t277;
t266 = qJDD(6) + t268;
t206 = m(7) * t209 + mrSges(7,1) * t266 - mrSges(7,3) * t236 - t242 * t278 + t249 * t302;
t210 = t212 * t348 + t215 * t352;
t235 = -qJD(6) * t278 + t267 * t352 - t298 * t348;
t250 = mrSges(7,1) * t302 - mrSges(7,3) * t278;
t207 = m(7) * t210 - mrSges(7,2) * t266 + mrSges(7,3) * t235 + t242 * t277 - t250 * t302;
t195 = t352 * t206 + t348 * t207;
t367 = -t317 * pkin(4) + t298 * qJ(5) - t304 * t273 + t221;
t214 = -t267 * pkin(5) - t303 * pkin(10) + ((2 * qJD(5)) + t283) * t318 + t367;
t237 = Ifges(7,5) * t278 + Ifges(7,6) * t277 + Ifges(7,3) * t302;
t239 = Ifges(7,1) * t278 + Ifges(7,4) * t277 + Ifges(7,5) * t302;
t198 = -mrSges(7,1) * t214 + mrSges(7,3) * t210 + Ifges(7,4) * t236 + Ifges(7,2) * t235 + Ifges(7,6) * t266 - t237 * t278 + t239 * t302;
t238 = Ifges(7,4) * t278 + Ifges(7,2) * t277 + Ifges(7,6) * t302;
t199 = mrSges(7,2) * t214 - mrSges(7,3) * t209 + Ifges(7,1) * t236 + Ifges(7,4) * t235 + Ifges(7,5) * t266 + t237 * t277 - t238 * t302;
t216 = t318 * t391 - t367;
t252 = Ifges(6,5) * t318 - Ifges(6,6) * t305 + Ifges(6,3) * t304;
t362 = -mrSges(6,2) * t218 + mrSges(6,3) * t216 - Ifges(6,1) * t298 + Ifges(6,4) * t268 - Ifges(6,5) * t267 + pkin(10) * t195 + t348 * t198 - t352 * t199 + t305 * t252;
t211 = -m(7) * t214 + t235 * mrSges(7,1) - t236 * mrSges(7,2) + t277 * t249 - t278 * t250;
t280 = mrSges(6,1) * t305 + mrSges(6,2) * t318;
t363 = -m(6) * t216 + t298 * mrSges(6,3) + t318 * t280 - t211;
t368 = -m(6) * t218 - t268 * mrSges(6,1) - t305 * t275 - t195;
t254 = Ifges(6,4) * t318 - Ifges(6,2) * t305 + Ifges(6,6) * t304;
t380 = Ifges(5,1) * t305 - Ifges(5,4) * t304 + Ifges(5,5) * t318 - t254;
t392 = t380 * t304 + mrSges(5,1) * t220 - mrSges(5,2) * t221 + Ifges(5,5) * t268 - Ifges(5,6) * t267 + Ifges(5,3) * t298 + pkin(4) * (-t298 * mrSges(6,2) - t318 * t279 + t368) + qJ(5) * (-t267 * mrSges(6,1) - t304 * t275 + t363) + t305 * t255 - t362;
t388 = Ifges(5,4) + Ifges(6,6);
t384 = t345 * t353;
t291 = -mrSges(4,1) * t319 + mrSges(4,2) * t320;
t307 = mrSges(4,1) * t339 - mrSges(4,3) * t320;
t274 = mrSges(5,1) * t304 + mrSges(5,2) * t305;
t281 = -mrSges(5,2) * t318 - mrSges(5,3) * t304;
t191 = m(5) * t220 - t268 * mrSges(5,3) - t305 * t274 + (-t279 + t281) * t318 + (mrSges(5,1) - mrSges(6,2)) * t298 + t368;
t282 = mrSges(5,1) * t318 - mrSges(5,3) * t305;
t202 = m(5) * t221 - t298 * mrSges(5,2) - t318 * t282 + (-t274 - t275) * t304 + (-mrSges(5,3) - mrSges(6,1)) * t267 + t363;
t372 = -t191 * t349 + t390 * t202;
t184 = m(4) * t228 - mrSges(4,2) * t338 + mrSges(4,3) * t299 + t291 * t319 - t307 * t339 + t372;
t306 = -mrSges(4,2) * t339 + mrSges(4,3) * t319;
t196 = -t348 * t206 + t352 * t207;
t219 = t267 * pkin(4) + t358;
t369 = -m(6) * t219 + t267 * mrSges(6,2) + t304 * t279 - t196;
t357 = -m(5) * t224 - t304 * t281 - t267 * mrSges(5,1) + (t280 - t282) * t305 + (-mrSges(5,2) + mrSges(6,3)) * t268 + t369;
t189 = m(4) * t227 + t338 * mrSges(4,1) - t300 * mrSges(4,3) - t320 * t291 + t339 * t306 + t357;
t181 = t344 * t184 + t346 * t189;
t187 = t390 * t191 + t349 * t202;
t256 = Ifges(6,1) * t318 - Ifges(6,4) * t305 + Ifges(6,5) * t304;
t381 = -Ifges(5,5) * t305 + Ifges(5,6) * t304 - Ifges(5,3) * t318 - t256;
t289 = -g(3) * t384 + t370;
t324 = -mrSges(3,2) * t339 + mrSges(3,3) * t374;
t327 = (-mrSges(3,1) * t353 + mrSges(3,2) * t350) * t379;
t179 = m(3) * t289 + mrSges(3,1) * t338 - mrSges(3,3) * t328 + t324 * t339 - t327 * t375 + t181;
t323 = mrSges(3,1) * t339 - mrSges(3,3) * t375;
t373 = t346 * t184 - t189 * t344;
t180 = m(3) * t290 - mrSges(3,2) * t338 + mrSges(3,3) * t329 - t323 * t339 + t327 * t374 + t373;
t170 = -t179 * t350 + t353 * t180;
t365 = m(4) * t271 - t299 * mrSges(4,1) + t300 * mrSges(4,2) - t319 * t306 + t320 * t307 + t187;
t185 = m(3) * t311 - t329 * mrSges(3,1) + t328 * mrSges(3,2) + (t323 * t350 - t324 * t353) * t379 + t365;
t167 = t179 * t382 + t180 * t383 - t185 * t345;
t194 = -t268 * mrSges(6,3) - t305 * t280 - t369;
t360 = -mrSges(6,1) * t216 + mrSges(6,2) * t219 - pkin(5) * t211 - pkin(10) * t196 - t352 * t198 - t348 * t199;
t173 = -mrSges(5,1) * t224 + mrSges(5,3) * t221 - pkin(4) * t194 + t380 * t318 + t381 * t305 + (Ifges(5,6) - Ifges(6,5)) * t298 + t388 * t268 + (-Ifges(5,2) - Ifges(6,3)) * t267 + t360;
t364 = mrSges(7,1) * t209 - mrSges(7,2) * t210 + Ifges(7,5) * t236 + Ifges(7,6) * t235 + Ifges(7,3) * t266 + t278 * t238 - t277 * t239;
t359 = mrSges(6,1) * t218 - mrSges(6,3) * t219 + pkin(5) * t195 + t364;
t175 = (-t255 + t252) * t318 + t381 * t304 + (Ifges(5,5) - Ifges(6,4)) * t298 + (Ifges(5,1) + Ifges(6,2)) * t268 - t388 * t267 + t359 + mrSges(5,2) * t224 - mrSges(5,3) * t220 - qJ(5) * t194;
t285 = Ifges(4,5) * t320 + Ifges(4,6) * t319 + Ifges(4,3) * t339;
t286 = Ifges(4,4) * t320 + Ifges(4,2) * t319 + Ifges(4,6) * t339;
t163 = mrSges(4,2) * t271 - mrSges(4,3) * t227 + Ifges(4,1) * t300 + Ifges(4,4) * t299 + Ifges(4,5) * t338 - pkin(9) * t187 - t349 * t173 + t390 * t175 + t319 * t285 - t339 * t286;
t287 = Ifges(4,1) * t320 + Ifges(4,4) * t319 + Ifges(4,5) * t339;
t171 = -mrSges(4,1) * t271 + mrSges(4,3) * t228 + Ifges(4,4) * t300 + Ifges(4,2) * t299 + Ifges(4,6) * t338 - pkin(3) * t187 - t320 * t285 + t339 * t287 - t392;
t308 = Ifges(3,3) * t339 + (Ifges(3,5) * t350 + Ifges(3,6) * t353) * t379;
t310 = Ifges(3,5) * t339 + (Ifges(3,1) * t350 + Ifges(3,4) * t353) * t379;
t158 = -mrSges(3,1) * t311 + mrSges(3,3) * t290 + Ifges(3,4) * t328 + Ifges(3,2) * t329 + Ifges(3,6) * t338 - pkin(2) * t365 + qJ(3) * t373 + t344 * t163 + t346 * t171 - t308 * t375 + t339 * t310;
t309 = Ifges(3,6) * t339 + (Ifges(3,4) * t350 + Ifges(3,2) * t353) * t379;
t160 = mrSges(3,2) * t311 - mrSges(3,3) * t289 + Ifges(3,1) * t328 + Ifges(3,4) * t329 + Ifges(3,5) * t338 - qJ(3) * t181 + t163 * t346 - t171 * t344 + t308 * t374 - t309 * t339;
t361 = mrSges(4,1) * t227 - mrSges(4,2) * t228 + Ifges(4,5) * t300 + Ifges(4,6) * t299 + Ifges(4,3) * t338 + pkin(3) * t357 + pkin(9) * t372 + t390 * t173 + t349 * t175 + t320 * t286 - t319 * t287;
t162 = t361 + (t309 * t350 - t310 * t353) * t379 + Ifges(3,3) * t338 + Ifges(3,5) * t328 + Ifges(3,6) * t329 + mrSges(3,1) * t289 - mrSges(3,2) * t290 + pkin(2) * t181;
t366 = mrSges(2,1) * t335 - mrSges(2,2) * t336 + Ifges(2,3) * qJDD(1) + pkin(1) * t167 + t158 * t384 + t160 * t385 + t347 * t162 + t170 * t389;
t168 = m(2) * t336 - mrSges(2,1) * t355 - qJDD(1) * mrSges(2,2) + t170;
t166 = t347 * t185 + (t179 * t353 + t180 * t350) * t345;
t164 = m(2) * t335 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t355 + t167;
t156 = -mrSges(2,2) * g(3) - mrSges(2,3) * t335 + Ifges(2,5) * qJDD(1) - t355 * Ifges(2,6) - t350 * t158 + t353 * t160 + (-t166 * t345 - t167 * t347) * pkin(8);
t155 = mrSges(2,1) * g(3) + mrSges(2,3) * t336 + t355 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t166 - t345 * t162 + (pkin(8) * t170 + t158 * t353 + t160 * t350) * t347;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t354 * t156 - t351 * t155 - pkin(7) * (t164 * t354 + t168 * t351), t156, t160, t163, t175, -t304 * t254 - t362, t199; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t351 * t156 + t354 * t155 + pkin(7) * (-t164 * t351 + t168 * t354), t155, t158, t171, t173, Ifges(6,4) * t298 - Ifges(6,2) * t268 + Ifges(6,6) * t267 - t318 * t252 + t304 * t256 - t359, t198; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t366, t366, t162, t361, t392, Ifges(6,5) * t298 - Ifges(6,6) * t268 + Ifges(6,3) * t267 + t318 * t254 + t305 * t256 - t360, t364;];
m_new  = t1;
