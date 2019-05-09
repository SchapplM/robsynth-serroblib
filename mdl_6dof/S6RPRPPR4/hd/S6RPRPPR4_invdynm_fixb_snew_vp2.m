% Calculate vector of cutting torques with Newton-Euler for
% S6RPRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-05-05 16:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRPPR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:47:51
% EndTime: 2019-05-05 16:48:15
% DurationCPUTime: 12.72s
% Computational Cost: add. (189749->364), mult. (469884->441), div. (0->0), fcn. (345189->10), ass. (0->149)
t400 = -2 * qJD(4);
t349 = qJD(1) ^ 2;
t345 = sin(qJ(1));
t347 = cos(qJ(1));
t325 = -t347 * g(1) - t345 * g(2);
t318 = -t349 * pkin(1) + qJDD(1) * qJ(2) + t325;
t341 = sin(pkin(9));
t342 = cos(pkin(9));
t380 = qJD(1) * qJD(2);
t376 = -t342 * g(3) - 0.2e1 * t341 * t380;
t393 = pkin(2) * t342;
t268 = (-pkin(7) * qJDD(1) + t349 * t393 - t318) * t341 + t376;
t303 = -t341 * g(3) + (t318 + 0.2e1 * t380) * t342;
t378 = qJDD(1) * t342;
t337 = t342 ^ 2;
t388 = t337 * t349;
t282 = -pkin(2) * t388 + pkin(7) * t378 + t303;
t344 = sin(qJ(3));
t394 = cos(qJ(3));
t241 = t344 * t268 + t394 * t282;
t377 = t342 * t394;
t384 = qJD(1) * t341;
t316 = -qJD(1) * t377 + t344 * t384;
t364 = t341 * t394 + t342 * t344;
t317 = t364 * qJD(1);
t289 = t316 * pkin(3) - t317 * qJ(4);
t348 = qJD(3) ^ 2;
t229 = -t348 * pkin(3) + qJDD(3) * qJ(4) - t316 * t289 + t241;
t336 = t341 ^ 2;
t324 = t345 * g(1) - t347 * g(2);
t371 = qJDD(2) - t324;
t299 = (-pkin(1) - t393) * qJDD(1) + (-qJ(2) + (-t336 - t337) * pkin(7)) * t349 + t371;
t379 = qJDD(1) * t341;
t381 = t317 * qJD(3);
t300 = -qJDD(1) * t377 + t344 * t379 + t381;
t382 = t316 * qJD(3);
t301 = qJDD(1) * t364 - t382;
t231 = (-t301 + t382) * qJ(4) + (t300 + t381) * pkin(3) + t299;
t340 = sin(pkin(10));
t390 = cos(pkin(10));
t308 = t340 * qJD(3) + t317 * t390;
t216 = -t340 * t229 + t231 * t390 + t308 * t400;
t399 = qJD(1) * t342;
t281 = t340 * qJDD(3) + t301 * t390;
t240 = t394 * t268 - t344 * t282;
t359 = qJDD(3) * pkin(3) + t348 * qJ(4) - t317 * t289 - qJDD(4) + t240;
t307 = -qJD(3) * t390 + t340 * t317;
t389 = t307 * t316;
t398 = (-t281 + t389) * qJ(5) - t359;
t217 = t390 * t229 + t340 * t231 + t307 * t400;
t247 = Ifges(6,5) * t308 + Ifges(6,6) * t316 + Ifges(6,3) * t307;
t250 = Ifges(5,4) * t308 - Ifges(5,2) * t307 + Ifges(5,6) * t316;
t252 = Ifges(5,1) * t308 - Ifges(5,4) * t307 + Ifges(5,5) * t316;
t261 = t307 * mrSges(6,1) - t308 * mrSges(6,3);
t280 = -qJDD(3) * t390 + t340 * t301;
t260 = t307 * pkin(4) - t308 * qJ(5);
t315 = t316 ^ 2;
t214 = -t300 * pkin(4) - t315 * qJ(5) + t308 * t260 + qJDD(5) - t216;
t206 = (-t281 - t389) * pkin(8) + (t307 * t308 - t300) * pkin(5) + t214;
t395 = 2 * qJD(5);
t212 = -t315 * pkin(4) + t300 * qJ(5) - t307 * t260 + t316 * t395 + t217;
t279 = -t316 * pkin(5) - t308 * pkin(8);
t306 = t307 ^ 2;
t207 = -t306 * pkin(5) + t280 * pkin(8) + t316 * t279 + t212;
t343 = sin(qJ(6));
t346 = cos(qJ(6));
t203 = t346 * t206 - t343 * t207;
t257 = t346 * t307 - t343 * t308;
t228 = t257 * qJD(6) + t343 * t280 + t346 * t281;
t258 = t343 * t307 + t346 * t308;
t237 = -t257 * mrSges(7,1) + t258 * mrSges(7,2);
t313 = qJD(6) - t316;
t245 = -t313 * mrSges(7,2) + t257 * mrSges(7,3);
t298 = qJDD(6) - t300;
t198 = m(7) * t203 + t298 * mrSges(7,1) - t228 * mrSges(7,3) - t258 * t237 + t313 * t245;
t204 = t343 * t206 + t346 * t207;
t227 = -t258 * qJD(6) + t346 * t280 - t343 * t281;
t246 = t313 * mrSges(7,1) - t258 * mrSges(7,3);
t199 = m(7) * t204 - t298 * mrSges(7,2) + t227 * mrSges(7,3) + t257 * t237 - t313 * t246;
t188 = t346 * t198 + t343 * t199;
t251 = Ifges(6,1) * t308 + Ifges(6,4) * t316 + Ifges(6,5) * t307;
t233 = Ifges(7,4) * t258 + Ifges(7,2) * t257 + Ifges(7,6) * t313;
t234 = Ifges(7,1) * t258 + Ifges(7,4) * t257 + Ifges(7,5) * t313;
t362 = mrSges(7,1) * t203 - mrSges(7,2) * t204 + Ifges(7,5) * t228 + Ifges(7,6) * t227 + Ifges(7,3) * t298 + t258 * t233 - t257 * t234;
t353 = mrSges(6,1) * t214 - mrSges(6,3) * t212 - Ifges(6,4) * t281 - Ifges(6,2) * t300 - Ifges(6,6) * t280 + pkin(5) * t188 - t307 * t251 + t362;
t275 = -t307 * mrSges(6,2) + t316 * mrSges(6,3);
t360 = -m(6) * t214 + t300 * mrSges(6,1) + t316 * t275 - t188;
t189 = -t343 * t198 + t346 * t199;
t278 = -t316 * mrSges(6,1) + t308 * mrSges(6,2);
t363 = m(6) * t212 + t300 * mrSges(6,3) + t316 * t278 + t189;
t397 = -(t250 - t247) * t308 - mrSges(5,1) * t216 + mrSges(5,2) * t217 - Ifges(5,5) * t281 + Ifges(5,6) * t280 - pkin(4) * (-t281 * mrSges(6,2) - t308 * t261 + t360) - qJ(5) * (-t280 * mrSges(6,2) - t307 * t261 + t363) - t307 * t252 + t353;
t290 = t316 * mrSges(4,1) + t317 * mrSges(4,2);
t310 = qJD(3) * mrSges(4,1) - t317 * mrSges(4,3);
t277 = t316 * mrSges(5,1) - t308 * mrSges(5,3);
t385 = -t307 * mrSges(5,1) - t308 * mrSges(5,2) - t261;
t392 = -mrSges(5,3) - mrSges(6,2);
t183 = m(5) * t217 - t300 * mrSges(5,2) - t316 * t277 + t280 * t392 + t307 * t385 + t363;
t276 = -t316 * mrSges(5,2) - t307 * mrSges(5,3);
t185 = m(5) * t216 + t300 * mrSges(5,1) + t316 * t276 + t281 * t392 + t308 * t385 + t360;
t373 = t390 * t183 - t340 * t185;
t177 = m(4) * t241 - qJDD(3) * mrSges(4,2) - t300 * mrSges(4,3) - qJD(3) * t310 - t316 * t290 + t373;
t309 = -qJD(3) * mrSges(4,2) - t316 * mrSges(4,3);
t209 = -t306 * pkin(8) + (-pkin(4) - pkin(5)) * t280 + (-pkin(4) * t316 + t279 + t395) * t308 - t398;
t205 = -m(7) * t209 + t227 * mrSges(7,1) - t228 * mrSges(7,2) + t257 * t245 - t258 * t246;
t215 = -0.2e1 * qJD(5) * t308 + (t308 * t316 + t280) * pkin(4) + t398;
t200 = m(6) * t215 + t280 * mrSges(6,1) - t281 * mrSges(6,3) + t307 * t275 - t308 * t278 + t205;
t351 = m(5) * t359 - t280 * mrSges(5,1) - t281 * mrSges(5,2) - t307 * t276 - t308 * t277 - t200;
t194 = m(4) * t240 + qJDD(3) * mrSges(4,1) - t301 * mrSges(4,3) + qJD(3) * t309 - t317 * t290 + t351;
t172 = t344 * t177 + t394 * t194;
t302 = -t341 * t318 + t376;
t232 = Ifges(7,5) * t258 + Ifges(7,6) * t257 + Ifges(7,3) * t313;
t191 = -mrSges(7,1) * t209 + mrSges(7,3) * t204 + Ifges(7,4) * t228 + Ifges(7,2) * t227 + Ifges(7,6) * t298 - t258 * t232 + t313 * t234;
t192 = mrSges(7,2) * t209 - mrSges(7,3) * t203 + Ifges(7,1) * t228 + Ifges(7,4) * t227 + Ifges(7,5) * t298 + t257 * t232 - t313 * t233;
t355 = -mrSges(6,1) * t215 + mrSges(6,2) * t212 - pkin(5) * t205 - pkin(8) * t189 - t346 * t191 - t343 * t192;
t249 = Ifges(6,4) * t308 + Ifges(6,2) * t316 + Ifges(6,6) * t307;
t387 = -Ifges(5,5) * t308 + Ifges(5,6) * t307 - Ifges(5,3) * t316 - t249;
t166 = mrSges(5,1) * t359 + mrSges(5,3) * t217 - pkin(4) * t200 + (t252 + t251) * t316 + t387 * t308 + (Ifges(5,6) - Ifges(6,6)) * t300 + (Ifges(5,4) - Ifges(6,5)) * t281 + (-Ifges(5,2) - Ifges(6,3)) * t280 + t355;
t357 = mrSges(6,2) * t214 - mrSges(6,3) * t215 + Ifges(6,1) * t281 + Ifges(6,4) * t300 + Ifges(6,5) * t280 - pkin(8) * t188 - t343 * t191 + t346 * t192 + t316 * t247;
t168 = -mrSges(5,2) * t359 - mrSges(5,3) * t216 + Ifges(5,1) * t281 - Ifges(5,4) * t280 + Ifges(5,5) * t300 - qJ(5) * t200 - t316 * t250 + t307 * t387 + t357;
t284 = Ifges(4,4) * t317 - Ifges(4,2) * t316 + Ifges(4,6) * qJD(3);
t285 = Ifges(4,1) * t317 - Ifges(4,4) * t316 + Ifges(4,5) * qJD(3);
t356 = -mrSges(4,1) * t240 + mrSges(4,2) * t241 - Ifges(4,5) * t301 + Ifges(4,6) * t300 - Ifges(4,3) * qJDD(3) - pkin(3) * t351 - qJ(4) * t373 - t390 * t166 - t340 * t168 - t317 * t284 - t316 * t285;
t369 = Ifges(3,4) * t341 + Ifges(3,2) * t342;
t370 = Ifges(3,1) * t341 + Ifges(3,4) * t342;
t396 = -mrSges(3,1) * t302 + mrSges(3,2) * t303 - pkin(2) * t172 - (t369 * t384 - t370 * t399) * qJD(1) + t356;
t391 = mrSges(3,2) * t341;
t179 = t340 * t183 + t390 * t185;
t365 = mrSges(3,3) * qJDD(1) + t349 * (-mrSges(3,1) * t342 + t391);
t170 = m(3) * t302 - t341 * t365 + t172;
t374 = t394 * t177 - t344 * t194;
t171 = m(3) * t303 + t342 * t365 + t374;
t375 = -t341 * t170 + t342 * t171;
t368 = Ifges(3,5) * t341 + Ifges(3,6) * t342;
t283 = Ifges(4,5) * t317 - Ifges(4,6) * t316 + Ifges(4,3) * qJD(3);
t160 = mrSges(4,2) * t299 - mrSges(4,3) * t240 + Ifges(4,1) * t301 - Ifges(4,4) * t300 + Ifges(4,5) * qJDD(3) - qJ(4) * t179 - qJD(3) * t284 - t340 * t166 + t168 * t390 - t316 * t283;
t161 = Ifges(4,6) * qJDD(3) + (-Ifges(4,2) - Ifges(5,3)) * t300 - t317 * t283 + Ifges(4,4) * t301 - mrSges(4,1) * t299 + qJD(3) * t285 + mrSges(4,3) * t241 - pkin(3) * t179 + t397;
t314 = -qJDD(1) * pkin(1) - t349 * qJ(2) + t371;
t320 = t368 * qJD(1);
t358 = m(4) * t299 + t300 * mrSges(4,1) + t301 * mrSges(4,2) + t316 * t309 + t317 * t310 + t179;
t156 = -mrSges(3,1) * t314 + mrSges(3,3) * t303 - pkin(2) * t358 + pkin(7) * t374 + qJDD(1) * t369 + t344 * t160 + t161 * t394 - t320 * t384;
t158 = mrSges(3,2) * t314 - mrSges(3,3) * t302 - pkin(7) * t172 + qJDD(1) * t370 + t160 * t394 - t344 * t161 + t320 * t399;
t354 = -m(3) * t314 + mrSges(3,1) * t378 - t358 + (t336 * t349 + t388) * mrSges(3,3);
t361 = -mrSges(2,2) * t325 + qJ(2) * t375 + t342 * t156 + t341 * t158 + pkin(1) * (-mrSges(3,2) * t379 + t354) + mrSges(2,1) * t324 + Ifges(2,3) * qJDD(1);
t173 = t354 - t349 * mrSges(2,2) + m(2) * t324 + (mrSges(2,1) - t391) * qJDD(1);
t164 = t342 * t170 + t341 * t171;
t162 = m(2) * t325 - t349 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t375;
t159 = -pkin(1) * t164 + mrSges(2,1) * g(3) + (Ifges(2,6) - t368) * qJDD(1) + t349 * Ifges(2,5) + mrSges(2,3) * t325 + t396;
t154 = -mrSges(2,2) * g(3) - mrSges(2,3) * t324 + Ifges(2,5) * qJDD(1) - t349 * Ifges(2,6) - qJ(2) * t164 - t341 * t156 + t342 * t158;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t347 * t154 - t345 * t159 - pkin(6) * (t345 * t162 + t347 * t173), t154, t158, t160, t168, -t307 * t249 + t357, t192; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t345 * t154 + t347 * t159 + pkin(6) * (t347 * t162 - t345 * t173), t159, t156, t161, t166, -t308 * t247 - t353, t191; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t361, t361, qJDD(1) * t368 - t396, -t356, Ifges(5,3) * t300 - t397, Ifges(6,5) * t281 + Ifges(6,6) * t300 + Ifges(6,3) * t280 + t308 * t249 - t316 * t251 - t355, t362;];
m_new  = t1;
