% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
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
% Datum: 2019-05-07 04:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPPP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPP1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPP1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPP1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPP1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPP1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 03:57:25
% EndTime: 2019-05-07 03:58:08
% DurationCPUTime: 19.90s
% Computational Cost: add. (349824->410), mult. (738729->481), div. (0->0), fcn. (539757->10), ass. (0->160)
t398 = -2 * qJD(4);
t337 = sin(qJ(1));
t339 = cos(qJ(1));
t328 = t337 * g(1) - t339 * g(2);
t341 = qJD(1) ^ 2;
t312 = -qJDD(1) * pkin(1) - t341 * pkin(8) - t328;
t336 = sin(qJ(2));
t392 = cos(qJ(2));
t373 = t392 * qJD(1);
t365 = qJD(2) * t373;
t322 = t336 * qJDD(1) + t365;
t378 = qJD(1) * t336;
t372 = qJD(2) * t378;
t323 = t392 * qJDD(1) - t372;
t275 = (-t365 - t322) * pkin(9) + (-t323 + t372) * pkin(2) + t312;
t329 = -g(1) * t339 - g(2) * t337;
t313 = -pkin(1) * t341 + qJDD(1) * pkin(8) + t329;
t303 = -t336 * g(3) + t392 * t313;
t321 = (-t392 * pkin(2) - pkin(9) * t336) * qJD(1);
t340 = qJD(2) ^ 2;
t284 = -t340 * pkin(2) + qJDD(2) * pkin(9) + t321 * t373 + t303;
t335 = sin(qJ(3));
t338 = cos(qJ(3));
t238 = t338 * t275 - t335 * t284;
t318 = qJD(2) * t338 - t335 * t378;
t292 = qJD(3) * t318 + qJDD(2) * t335 + t322 * t338;
t319 = qJD(2) * t335 + t338 * t378;
t334 = sin(pkin(6));
t385 = qJ(4) * t334;
t293 = -pkin(3) * t318 - t319 * t385;
t363 = t373 - qJD(3);
t361 = t363 * t334;
t387 = cos(pkin(6));
t350 = t318 * t387 - t361;
t294 = t350 * qJ(4);
t317 = qJDD(3) - t323;
t371 = qJ(4) * t387;
t213 = t317 * pkin(3) - t292 * t371 - t319 * t293 - t363 * t294 + t238;
t239 = t335 * t275 + t338 * t284;
t299 = -t363 * pkin(3) - t319 * t371;
t291 = -qJD(3) * t319 + qJDD(2) * t338 - t322 * t335;
t360 = t387 * t291 + t317 * t334;
t214 = t360 * qJ(4) + t318 * t293 + t363 * t299 + t239;
t302 = -t392 * g(3) - t336 * t313;
t283 = -qJDD(2) * pkin(2) - pkin(9) * t340 + t321 * t378 - t302;
t217 = -pkin(3) * t291 - t292 * t385 - t294 * t318 + t299 * t319 + t283;
t333 = sin(pkin(10));
t386 = cos(pkin(10));
t281 = t386 * t319 + t350 * t333;
t364 = t387 * t386;
t369 = t334 * t386;
t207 = t213 * t364 - t333 * t214 + t217 * t369 + t281 * t398;
t255 = -t291 * t364 + t333 * t292 - t317 * t369;
t209 = -t334 * t213 + t387 * t217 + qJDD(4);
t256 = t386 * t292 + t360 * t333;
t297 = t334 * t318 + t387 * t363;
t280 = -t318 * t364 + t333 * t319 + t386 * t361;
t384 = t280 * t297;
t394 = -2 * qJD(5);
t348 = (-t256 - t384) * qJ(5) + t209 + (-pkin(4) * t297 + t394) * t281;
t206 = t255 * pkin(4) + t348;
t263 = mrSges(6,1) * t281 - mrSges(6,2) * t297;
t397 = m(6) * t206 - t256 * mrSges(6,3) - t281 * t263;
t277 = t280 * t398;
t370 = t333 * t387;
t375 = t334 * t333 * t217 + t213 * t370 + t386 * t214;
t208 = t277 + t375;
t227 = -Ifges(7,4) * t297 + Ifges(7,2) * t280 + Ifges(7,6) * t281;
t229 = Ifges(5,4) * t281 - Ifges(5,2) * t280 - Ifges(5,6) * t297;
t240 = -mrSges(7,2) * t281 + mrSges(7,3) * t280;
t261 = mrSges(6,1) * t280 + mrSges(6,3) * t297;
t274 = -t334 * t291 + t387 * t317;
t241 = pkin(4) * t280 - qJ(5) * t281;
t296 = t297 ^ 2;
t204 = -t274 * pkin(4) - t296 * qJ(5) + t281 * t241 + qJDD(5) - t207;
t393 = 2 * qJD(6);
t195 = t297 * t393 + (t280 * t281 - t274) * qJ(6) + (t256 - t384) * pkin(5) + t204;
t262 = -mrSges(7,1) * t280 - mrSges(7,2) * t297;
t366 = -m(7) * t195 + t274 * mrSges(7,3) - t297 * t262;
t388 = t256 * mrSges(7,1);
t191 = t281 * t240 - t366 + t388;
t354 = t296 * pkin(4) - t274 * qJ(5) - t375;
t202 = 0.2e1 * qJD(5) * t297 + ((2 * qJD(4)) + t241) * t280 + t354;
t225 = -Ifges(6,5) * t297 - Ifges(6,6) * t281 + Ifges(6,3) * t280;
t259 = pkin(5) * t281 + qJ(6) * t297;
t279 = t280 ^ 2;
t198 = -t255 * pkin(5) - t279 * qJ(6) - t280 * t241 + qJDD(6) + t277 + (t394 - t259) * t297 - t354;
t224 = -Ifges(7,5) * t297 + Ifges(7,6) * t280 + Ifges(7,3) * t281;
t355 = mrSges(7,2) * t198 - mrSges(7,3) * t195 + Ifges(7,1) * t274 + Ifges(7,4) * t255 + Ifges(7,5) * t256 + t280 * t224;
t344 = mrSges(6,2) * t204 - mrSges(6,3) * t202 + Ifges(6,1) * t274 - Ifges(6,4) * t256 + Ifges(6,5) * t255 - qJ(6) * t191 - t281 * t225 + t355;
t243 = -mrSges(6,2) * t280 - mrSges(6,3) * t281;
t352 = -m(6) * t204 - t256 * mrSges(6,1) - t281 * t243 + t366;
t260 = mrSges(7,1) * t281 + mrSges(7,3) * t297;
t374 = -m(7) * t198 - t274 * mrSges(7,2) + t297 * t260;
t358 = -m(6) * t202 + t274 * mrSges(6,3) - t297 * t263 - t374;
t380 = -t240 - t243;
t228 = -Ifges(6,4) * t297 - Ifges(6,2) * t281 + Ifges(6,6) * t280;
t382 = -Ifges(5,1) * t281 + Ifges(5,4) * t280 + Ifges(5,5) * t297 + t228;
t173 = qJ(5) * t358 + pkin(4) * (-t274 * mrSges(6,2) + t297 * t261 + t352 - t388) + Ifges(5,3) * t274 + Ifges(5,5) * t256 - mrSges(5,2) * t208 + mrSges(5,1) * t207 + t344 + (-pkin(4) * t240 - t227 + t229) * t281 + (qJ(5) * t380 - t382) * t280 + (-Ifges(5,6) + qJ(5) * (-mrSges(6,1) - mrSges(7,1))) * t255;
t201 = -t279 * pkin(5) + t280 * t393 - t281 * t259 + (pkin(4) + qJ(6)) * t255 + t348;
t376 = m(7) * t201 + t255 * mrSges(7,3) + t280 * t262;
t193 = -t256 * mrSges(7,2) - t281 * t260 + t376;
t190 = -t255 * mrSges(6,2) - t280 * t261 + t193 + t397;
t226 = Ifges(5,5) * t281 - Ifges(5,6) * t280 - Ifges(5,3) * t297;
t231 = -Ifges(6,1) * t297 - Ifges(6,4) * t281 + Ifges(6,5) * t280;
t230 = -Ifges(7,1) * t297 + Ifges(7,4) * t280 + Ifges(7,5) * t281;
t357 = mrSges(7,1) * t198 - mrSges(7,3) * t201 - Ifges(7,4) * t274 - Ifges(7,2) * t255 - Ifges(7,6) * t256 - t281 * t230;
t345 = mrSges(6,1) * t202 - mrSges(6,2) * t206 + pkin(5) * (t255 * mrSges(7,1) + t280 * t240 + t374) + qJ(6) * t193 - t357;
t389 = Ifges(5,4) + Ifges(6,6);
t174 = mrSges(5,3) * t208 - mrSges(5,1) * t209 - pkin(4) * t190 - t345 + (-Ifges(5,2) - Ifges(6,3)) * t255 + t389 * t256 + (Ifges(5,6) - Ifges(6,5)) * t274 + (-t226 - t231) * t281 + (-t224 + t382) * t297;
t242 = mrSges(5,1) * t280 + mrSges(5,2) * t281;
t379 = mrSges(5,2) * t297 - mrSges(5,3) * t280 - t261;
t390 = -mrSges(7,1) - mrSges(5,3);
t391 = mrSges(5,1) - mrSges(6,2);
t185 = m(5) * t207 - t379 * t297 + (-t240 - t242) * t281 + t391 * t274 + t390 * t256 + t352;
t258 = -mrSges(5,1) * t297 - mrSges(5,3) * t281;
t188 = m(5) * t208 - t274 * mrSges(5,2) + t297 * t258 + (-t242 + t380) * t280 + (-mrSges(6,1) + t390) * t255 + t358;
t189 = m(5) * t209 + (t258 - t260) * t281 + t379 * t280 + (mrSges(5,2) - mrSges(7,2)) * t256 + t391 * t255 + t376 + t397;
t179 = t185 * t364 + t188 * t370 - t334 * t189;
t356 = -mrSges(7,1) * t195 + mrSges(7,2) * t201 - Ifges(7,5) * t274 - Ifges(7,6) * t255 - Ifges(7,3) * t256 + t297 * t227;
t347 = -mrSges(6,1) * t204 + mrSges(6,3) * t206 - pkin(5) * t191 + t356;
t381 = t230 + t231;
t180 = (t229 - t225) * t297 + (-t226 - t381) * t280 + (Ifges(5,5) - Ifges(6,4)) * t274 + (Ifges(5,1) + Ifges(6,2)) * t256 - t389 * t255 - t347 + mrSges(5,2) * t209 - mrSges(5,3) * t207 - qJ(5) * t190;
t183 = -t333 * t185 + t386 * t188;
t287 = Ifges(4,4) * t319 + Ifges(4,2) * t318 - Ifges(4,6) * t363;
t288 = Ifges(4,1) * t319 + Ifges(4,4) * t318 - Ifges(4,5) * t363;
t396 = mrSges(4,1) * t238 - mrSges(4,2) * t239 + Ifges(4,5) * t292 + Ifges(4,6) * t291 + Ifges(4,3) * t317 + pkin(3) * t179 + t387 * t173 + t319 * t287 - t318 * t288 + (qJ(4) * t183 + t386 * t174 + t180 * t333) * t334;
t178 = (t386 * t185 + t188 * t333) * t334 + t387 * t189;
t286 = Ifges(4,5) * t319 + Ifges(4,6) * t318 - Ifges(4,3) * t363;
t162 = -mrSges(4,1) * t283 + mrSges(4,3) * t239 + Ifges(4,4) * t292 + Ifges(4,2) * t291 + Ifges(4,6) * t317 - pkin(3) * t178 - t334 * t173 + t174 * t364 + t180 * t370 + t183 * t371 - t319 * t286 - t363 * t288;
t163 = Ifges(4,1) * t292 + Ifges(4,4) * t291 + Ifges(4,5) * t317 + t318 * t286 + t363 * t287 + mrSges(4,2) * t283 - mrSges(4,3) * t238 + t386 * t180 - t333 * t174 + (-t334 * t178 - t387 * t179) * qJ(4);
t295 = -mrSges(4,1) * t318 + mrSges(4,2) * t319;
t300 = t363 * mrSges(4,2) + t318 * mrSges(4,3);
t176 = m(4) * t238 + t317 * mrSges(4,1) - t292 * mrSges(4,3) - t319 * t295 - t363 * t300 + t179;
t301 = -t363 * mrSges(4,1) - t319 * mrSges(4,3);
t182 = m(4) * t239 - t317 * mrSges(4,2) + t291 * mrSges(4,3) + t318 * t295 + t363 * t301 + t183;
t172 = -t335 * t176 + t338 * t182;
t177 = -m(4) * t283 + t291 * mrSges(4,1) - t292 * mrSges(4,2) + t318 * t300 - t319 * t301 - t178;
t310 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t336 + t392 * Ifges(3,2)) * qJD(1);
t311 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t336 + t392 * Ifges(3,4)) * qJD(1);
t395 = mrSges(3,1) * t302 - mrSges(3,2) * t303 + Ifges(3,5) * t322 + Ifges(3,6) * t323 + Ifges(3,3) * qJDD(2) + pkin(2) * t177 + pkin(9) * t172 - (-t310 * t336 + t392 * t311) * qJD(1) + t338 * t162 + t335 * t163;
t383 = t281 * t227;
t320 = (-mrSges(3,1) * t392 + mrSges(3,2) * t336) * qJD(1);
t326 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t378;
t170 = m(3) * t303 - qJDD(2) * mrSges(3,2) + t323 * mrSges(3,3) - qJD(2) * t326 + t320 * t373 + t172;
t327 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t373;
t175 = m(3) * t302 + qJDD(2) * mrSges(3,1) - t322 * mrSges(3,3) + qJD(2) * t327 - t320 * t378 + t177;
t367 = t392 * t170 - t175 * t336;
t171 = t176 * t338 + t182 * t335;
t309 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t336 + t392 * Ifges(3,6)) * qJD(1);
t159 = mrSges(3,2) * t312 - mrSges(3,3) * t302 + Ifges(3,1) * t322 + Ifges(3,4) * t323 + Ifges(3,5) * qJDD(2) - pkin(9) * t171 - qJD(2) * t310 - t335 * t162 + t338 * t163 + t309 * t373;
t161 = -mrSges(3,1) * t312 + mrSges(3,3) * t303 + Ifges(3,4) * t322 + Ifges(3,2) * t323 + Ifges(3,6) * qJDD(2) - pkin(2) * t171 + qJD(2) * t311 - t309 * t378 - t396;
t346 = -m(3) * t312 + t323 * mrSges(3,1) - mrSges(3,2) * t322 - t326 * t378 + t327 * t373 - t171;
t353 = mrSges(2,1) * t328 - mrSges(2,2) * t329 + Ifges(2,3) * qJDD(1) + pkin(1) * t346 + pkin(8) * t367 + t336 * t159 + t392 * t161;
t167 = m(2) * t328 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t341 + t346;
t166 = t336 * t170 + t392 * t175;
t164 = m(2) * t329 - mrSges(2,1) * t341 - qJDD(1) * mrSges(2,2) + t367;
t157 = mrSges(2,1) * g(3) + mrSges(2,3) * t329 + t341 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t166 - t395;
t156 = -mrSges(2,2) * g(3) - mrSges(2,3) * t328 + Ifges(2,5) * qJDD(1) - t341 * Ifges(2,6) - pkin(8) * t166 + t392 * t159 - t336 * t161;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t339 * t156 - t337 * t157 - pkin(7) * (t164 * t337 + t167 * t339), t156, t159, t163, t180, -t280 * t228 + t344 - t383, t355 - t383; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t337 * t156 + t339 * t157 + pkin(7) * (t164 * t339 - t167 * t337), t157, t161, t162, t174, Ifges(6,4) * t274 - Ifges(6,2) * t256 + Ifges(6,6) * t255 + t297 * t225 + t381 * t280 + t347, t297 * t224 - t357; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t353, t353, t395, t396, t173, t281 * t231 + Ifges(6,5) * t274 + Ifges(6,3) * t255 - Ifges(6,6) * t256 + (t224 - t228) * t297 + t345, -t280 * t230 - t356;];
m_new  = t1;
