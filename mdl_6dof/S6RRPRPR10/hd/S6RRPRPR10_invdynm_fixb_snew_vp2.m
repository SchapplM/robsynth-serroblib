% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRPR10
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
% Datum: 2019-05-06 15:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRPR10_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR10_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR10_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR10_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR10_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR10_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR10_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 15:37:24
% EndTime: 2019-05-06 15:38:43
% DurationCPUTime: 31.24s
% Computational Cost: add. (514935->400), mult. (1173867->498), div. (0->0), fcn. (934177->12), ass. (0->158)
t351 = sin(qJ(2));
t354 = cos(qJ(2));
t346 = sin(pkin(6));
t377 = qJD(1) * t346;
t328 = (-pkin(2) * t354 - qJ(3) * t351) * t377;
t348 = cos(pkin(6));
t341 = qJD(1) * t348 + qJD(2);
t339 = t341 ^ 2;
t340 = qJDD(1) * t348 + qJDD(2);
t376 = qJD(1) * t354;
t352 = sin(qJ(1));
t355 = cos(qJ(1));
t336 = t352 * g(1) - g(2) * t355;
t356 = qJD(1) ^ 2;
t388 = pkin(8) * t346;
t326 = qJDD(1) * pkin(1) + t356 * t388 + t336;
t337 = -g(1) * t355 - g(2) * t352;
t375 = qJDD(1) * t346;
t327 = -pkin(1) * t356 + pkin(8) * t375 + t337;
t382 = t348 * t351;
t378 = t326 * t382 + t354 * t327;
t275 = -pkin(2) * t339 + qJ(3) * t340 + (-g(3) * t351 + t328 * t376) * t346 + t378;
t330 = (qJD(2) * t376 + qJDD(1) * t351) * t346;
t374 = t351 * t377;
t331 = -qJD(2) * t374 + t354 * t375;
t387 = g(3) * t348;
t276 = -pkin(2) * t331 - t387 - qJ(3) * t330 + (-t326 + (pkin(2) * t351 - qJ(3) * t354) * t341 * qJD(1)) * t346;
t345 = sin(pkin(11));
t347 = cos(pkin(11));
t317 = t341 * t345 + t347 * t374;
t228 = -0.2e1 * qJD(3) * t317 - t275 * t345 + t347 * t276;
t304 = t330 * t347 + t340 * t345;
t316 = t341 * t347 - t345 * t374;
t373 = t346 * t376;
t223 = (-t316 * t373 - t304) * pkin(9) + (t316 * t317 - t331) * pkin(3) + t228;
t229 = 0.2e1 * qJD(3) * t316 + t347 * t275 + t345 * t276;
t303 = -t330 * t345 + t340 * t347;
t305 = -pkin(3) * t373 - pkin(9) * t317;
t315 = t316 ^ 2;
t226 = -pkin(3) * t315 + pkin(9) * t303 + t305 * t373 + t229;
t350 = sin(qJ(4));
t389 = cos(qJ(4));
t218 = t223 * t389 - t350 * t226;
t219 = t350 * t223 + t389 * t226;
t297 = t350 * t316 + t317 * t389;
t255 = qJD(4) * t297 - t303 * t389 + t304 * t350;
t296 = -t316 * t389 + t317 * t350;
t256 = -t296 * qJD(4) + t350 * t303 + t304 * t389;
t334 = -qJD(4) + t373;
t262 = Ifges(5,4) * t297 - Ifges(5,2) * t296 - Ifges(5,6) * t334;
t270 = -mrSges(6,2) * t296 - mrSges(6,3) * t297;
t280 = mrSges(6,1) * t296 + mrSges(6,3) * t334;
t323 = qJDD(4) - t331;
t268 = pkin(4) * t296 - qJ(5) * t297;
t333 = t334 ^ 2;
t216 = -t323 * pkin(4) - t333 * qJ(5) + t297 * t268 + qJDD(5) - t218;
t385 = t296 * t334;
t210 = (t296 * t297 - t323) * pkin(10) + (t256 - t385) * pkin(5) + t216;
t284 = pkin(5) * t297 + pkin(10) * t334;
t295 = t296 ^ 2;
t381 = t348 * t354;
t383 = t346 * t354;
t292 = -g(3) * t383 + t326 * t381 - t351 * t327;
t274 = -pkin(2) * t340 - qJ(3) * t339 + t328 * t374 + qJDD(3) - t292;
t235 = -pkin(3) * t303 - pkin(9) * t315 + t317 * t305 + t274;
t390 = -2 * qJD(5);
t358 = (-t256 - t385) * qJ(5) + t235 + (-t334 * pkin(4) + t390) * t297;
t213 = t358 + (pkin(4) + pkin(10)) * t255 - pkin(5) * t295 - t284 * t297;
t349 = sin(qJ(6));
t353 = cos(qJ(6));
t207 = t210 * t353 - t213 * t349;
t278 = t296 * t353 + t334 * t349;
t234 = qJD(6) * t278 + t255 * t349 + t323 * t353;
t279 = t296 * t349 - t334 * t353;
t244 = -mrSges(7,1) * t278 + mrSges(7,2) * t279;
t254 = qJDD(6) + t256;
t294 = qJD(6) + t297;
t257 = -mrSges(7,2) * t294 + mrSges(7,3) * t278;
t204 = m(7) * t207 + mrSges(7,1) * t254 - mrSges(7,3) * t234 - t244 * t279 + t257 * t294;
t208 = t210 * t349 + t213 * t353;
t233 = -qJD(6) * t279 + t255 * t353 - t323 * t349;
t258 = mrSges(7,1) * t294 - mrSges(7,3) * t279;
t205 = m(7) * t208 - mrSges(7,2) * t254 + mrSges(7,3) * t233 + t244 * t278 - t258 * t294;
t192 = t204 * t353 + t205 * t349;
t368 = -pkin(4) * t333 + qJ(5) * t323 - t268 * t296 + t219;
t212 = -pkin(5) * t255 - pkin(10) * t295 + (t390 - t284) * t334 + t368;
t236 = Ifges(7,5) * t279 + Ifges(7,6) * t278 + Ifges(7,3) * t294;
t238 = Ifges(7,1) * t279 + Ifges(7,4) * t278 + Ifges(7,5) * t294;
t195 = -mrSges(7,1) * t212 + mrSges(7,3) * t208 + Ifges(7,4) * t234 + Ifges(7,2) * t233 + Ifges(7,6) * t254 - t236 * t279 + t238 * t294;
t237 = Ifges(7,4) * t279 + Ifges(7,2) * t278 + Ifges(7,6) * t294;
t196 = mrSges(7,2) * t212 - mrSges(7,3) * t207 + Ifges(7,1) * t234 + Ifges(7,4) * t233 + Ifges(7,5) * t254 + t236 * t278 - t237 * t294;
t214 = 0.2e1 * qJD(5) * t334 - t368;
t259 = -Ifges(6,5) * t334 - Ifges(6,6) * t297 + Ifges(6,3) * t296;
t364 = -mrSges(6,2) * t216 + mrSges(6,3) * t214 - Ifges(6,1) * t323 + Ifges(6,4) * t256 - Ifges(6,5) * t255 + pkin(10) * t192 + t349 * t195 - t353 * t196 + t297 * t259;
t209 = -m(7) * t212 + mrSges(7,1) * t233 - t234 * mrSges(7,2) + t257 * t278 - t279 * t258;
t281 = mrSges(6,1) * t297 - mrSges(6,2) * t334;
t365 = -m(6) * t214 + t323 * mrSges(6,3) - t334 * t281 - t209;
t369 = -m(6) * t216 - t256 * mrSges(6,1) - t297 * t270 - t192;
t261 = -Ifges(6,4) * t334 - Ifges(6,2) * t297 + Ifges(6,6) * t296;
t379 = -Ifges(5,1) * t297 + Ifges(5,4) * t296 + Ifges(5,5) * t334 + t261;
t392 = -mrSges(5,2) * t219 + pkin(4) * (-mrSges(6,2) * t323 + t280 * t334 + t369) + qJ(5) * (-mrSges(6,1) * t255 - t270 * t296 + t365) + mrSges(5,1) * t218 + t297 * t262 - Ifges(5,6) * t255 + Ifges(5,5) * t256 + Ifges(5,3) * t323 - t364 - t379 * t296;
t269 = mrSges(5,1) * t296 + mrSges(5,2) * t297;
t282 = mrSges(5,2) * t334 - mrSges(5,3) * t296;
t186 = m(5) * t218 - mrSges(5,3) * t256 - t269 * t297 + (t280 - t282) * t334 + (mrSges(5,1) - mrSges(6,2)) * t323 + t369;
t283 = -mrSges(5,1) * t334 - mrSges(5,3) * t297;
t199 = m(5) * t219 - mrSges(5,2) * t323 + t283 * t334 + (-t269 - t270) * t296 + (-mrSges(5,3) - mrSges(6,1)) * t255 + t365;
t184 = t389 * t186 + t350 * t199;
t287 = Ifges(4,4) * t317 + Ifges(4,2) * t316 - Ifges(4,6) * t373;
t288 = Ifges(4,1) * t317 + Ifges(4,4) * t316 - Ifges(4,5) * t373;
t391 = -mrSges(4,1) * t228 + mrSges(4,2) * t229 - Ifges(4,5) * t304 - Ifges(4,6) * t303 - pkin(3) * t184 - t317 * t287 + t316 * t288 - t392;
t386 = Ifges(5,4) + Ifges(6,6);
t384 = t346 * t351;
t298 = -mrSges(4,1) * t316 + mrSges(4,2) * t317;
t301 = mrSges(4,2) * t373 + mrSges(4,3) * t316;
t182 = m(4) * t228 - mrSges(4,1) * t331 - mrSges(4,3) * t304 - t298 * t317 - t301 * t373 + t184;
t302 = -mrSges(4,1) * t373 - mrSges(4,3) * t317;
t370 = -t186 * t350 + t389 * t199;
t183 = m(4) * t229 + mrSges(4,2) * t331 + mrSges(4,3) * t303 + t298 * t316 + t302 * t373 + t370;
t177 = t347 * t182 + t345 * t183;
t193 = -t349 * t204 + t353 * t205;
t263 = -Ifges(6,1) * t334 - Ifges(6,4) * t297 + Ifges(6,5) * t296;
t380 = -Ifges(5,5) * t297 + Ifges(5,6) * t296 + Ifges(5,3) * t334 - t263;
t293 = -g(3) * t384 + t378;
t324 = mrSges(3,1) * t341 - mrSges(3,3) * t374;
t329 = (-mrSges(3,1) * t354 + mrSges(3,2) * t351) * t377;
t371 = -t182 * t345 + t347 * t183;
t175 = m(3) * t293 - mrSges(3,2) * t340 + mrSges(3,3) * t331 - t324 * t341 + t329 * t373 + t371;
t325 = -mrSges(3,2) * t341 + mrSges(3,3) * t373;
t221 = pkin(4) * t255 + t358;
t191 = m(6) * t221 - t255 * mrSges(6,2) - t256 * mrSges(6,3) - t296 * t280 - t297 * t281 + t193;
t362 = m(5) * t235 + t255 * mrSges(5,1) + t256 * mrSges(5,2) + t296 * t282 + t297 * t283 + t191;
t359 = -m(4) * t274 + t303 * mrSges(4,1) - t304 * mrSges(4,2) + t316 * t301 - t317 * t302 - t362;
t188 = m(3) * t292 + t340 * mrSges(3,1) - t330 * mrSges(3,3) + t341 * t325 - t329 * t374 + t359;
t171 = t354 * t175 - t188 * t351;
t309 = -t326 * t346 - t387;
t176 = m(3) * t309 - mrSges(3,1) * t331 + mrSges(3,2) * t330 + (t324 * t351 - t325 * t354) * t377 + t177;
t168 = t175 * t382 - t176 * t346 + t188 * t381;
t363 = -mrSges(6,1) * t214 + mrSges(6,2) * t221 - pkin(5) * t209 - pkin(10) * t193 - t353 * t195 - t349 * t196;
t172 = -mrSges(5,1) * t235 + mrSges(5,3) * t219 - pkin(4) * t191 + t379 * t334 + (Ifges(5,6) - Ifges(6,5)) * t323 + t380 * t297 + t386 * t256 + (-Ifges(5,2) - Ifges(6,3)) * t255 + t363;
t366 = mrSges(7,1) * t207 - mrSges(7,2) * t208 + Ifges(7,5) * t234 + Ifges(7,6) * t233 + Ifges(7,3) * t254 + t279 * t237 - t278 * t238;
t361 = mrSges(6,1) * t216 - mrSges(6,3) * t221 + pkin(5) * t192 + t366;
t178 = mrSges(5,2) * t235 - mrSges(5,3) * t218 - qJ(5) * t191 + (t262 - t259) * t334 + (Ifges(5,5) - Ifges(6,4)) * t323 + t361 + t380 * t296 + (Ifges(5,1) + Ifges(6,2)) * t256 - t386 * t255;
t286 = Ifges(4,5) * t317 + Ifges(4,6) * t316 - Ifges(4,3) * t373;
t161 = -mrSges(4,1) * t274 + mrSges(4,3) * t229 + Ifges(4,4) * t304 + Ifges(4,2) * t303 - Ifges(4,6) * t331 - pkin(3) * t362 + pkin(9) * t370 + t172 * t389 + t350 * t178 - t317 * t286 - t288 * t373;
t164 = mrSges(4,2) * t274 - mrSges(4,3) * t228 + Ifges(4,1) * t304 + Ifges(4,4) * t303 - Ifges(4,5) * t331 - pkin(9) * t184 - t350 * t172 + t178 * t389 + t316 * t286 + t287 * t373;
t307 = Ifges(3,6) * t341 + (Ifges(3,4) * t351 + Ifges(3,2) * t354) * t377;
t308 = Ifges(3,5) * t341 + (Ifges(3,1) * t351 + Ifges(3,4) * t354) * t377;
t158 = Ifges(3,5) * t330 + Ifges(3,6) * t331 + Ifges(3,3) * t340 + mrSges(3,1) * t292 - mrSges(3,2) * t293 + t345 * t164 + t347 * t161 + pkin(2) * t359 + qJ(3) * t371 + (t307 * t351 - t308 * t354) * t377;
t306 = Ifges(3,3) * t341 + (Ifges(3,5) * t351 + Ifges(3,6) * t354) * t377;
t160 = mrSges(3,2) * t309 - mrSges(3,3) * t292 + Ifges(3,1) * t330 + Ifges(3,4) * t331 + Ifges(3,5) * t340 - qJ(3) * t177 - t161 * t345 + t164 * t347 + t306 * t373 - t307 * t341;
t163 = t341 * t308 + Ifges(3,6) * t340 + Ifges(3,4) * t330 - mrSges(3,1) * t309 + mrSges(3,3) * t293 - t306 * t374 - pkin(2) * t177 + (Ifges(3,2) + Ifges(4,3)) * t331 + t391;
t367 = mrSges(2,1) * t336 - mrSges(2,2) * t337 + Ifges(2,3) * qJDD(1) + pkin(1) * t168 + t348 * t158 + t160 * t384 + t163 * t383 + t171 * t388;
t169 = m(2) * t337 - mrSges(2,1) * t356 - qJDD(1) * mrSges(2,2) + t171;
t167 = t176 * t348 + (t175 * t351 + t188 * t354) * t346;
t165 = m(2) * t336 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t356 + t168;
t156 = -mrSges(2,2) * g(3) - mrSges(2,3) * t336 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t356 + t160 * t354 - t163 * t351 + (-t167 * t346 - t168 * t348) * pkin(8);
t155 = mrSges(2,1) * g(3) + mrSges(2,3) * t337 + Ifges(2,5) * t356 + Ifges(2,6) * qJDD(1) - pkin(1) * t167 - t158 * t346 + (pkin(8) * t171 + t160 * t351 + t163 * t354) * t348;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t355 * t156 - t352 * t155 - pkin(7) * (t165 * t355 + t169 * t352), t156, t160, t164, t178, -t296 * t261 - t364, t196; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t352 * t156 + t355 * t155 + pkin(7) * (-t165 * t352 + t169 * t355), t155, t163, t161, t172, Ifges(6,4) * t323 - Ifges(6,2) * t256 + Ifges(6,6) * t255 + t334 * t259 + t296 * t263 - t361, t195; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t367, t367, t158, -Ifges(4,3) * t331 - t391, t392, Ifges(6,5) * t323 - Ifges(6,6) * t256 + Ifges(6,3) * t255 - t334 * t261 + t297 * t263 - t363, t366;];
m_new  = t1;
