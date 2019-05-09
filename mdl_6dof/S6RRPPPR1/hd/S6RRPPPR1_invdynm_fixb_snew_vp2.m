% Calculate vector of cutting torques with Newton-Euler for
% S6RRPPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
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
% Datum: 2019-05-06 08:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPPPR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:10:32
% EndTime: 2019-05-06 08:10:58
% DurationCPUTime: 13.15s
% Computational Cost: add. (206071->387), mult. (489808->471), div. (0->0), fcn. (339511->10), ass. (0->148)
t395 = -2 * qJD(4);
t348 = sin(qJ(2));
t351 = cos(qJ(2));
t376 = qJD(1) * qJD(2);
t330 = qJDD(1) * t348 + t351 * t376;
t331 = qJDD(1) * t351 - t348 * t376;
t346 = sin(pkin(9));
t388 = cos(pkin(9));
t301 = t388 * t330 + t346 * t331;
t345 = sin(pkin(10));
t387 = cos(pkin(10));
t283 = t345 * qJDD(2) + t387 * t301;
t319 = (t346 * t351 + t388 * t348) * qJD(1);
t305 = -t387 * qJD(2) + t345 * t319;
t379 = qJD(1) * t351;
t380 = qJD(1) * t348;
t318 = t346 * t380 - t388 * t379;
t386 = t305 * t318;
t394 = (-t283 + t386) * qJ(5);
t349 = sin(qJ(1));
t352 = cos(qJ(1));
t337 = -g(1) * t352 - g(2) * t349;
t354 = qJD(1) ^ 2;
t325 = -pkin(1) * t354 + qJDD(1) * pkin(7) + t337;
t385 = t348 * t325;
t390 = pkin(2) * t354;
t264 = qJDD(2) * pkin(2) - t330 * qJ(3) - t385 + (qJ(3) * t376 + t348 * t390 - g(3)) * t351;
t308 = -g(3) * t348 + t351 * t325;
t333 = qJD(2) * pkin(2) - qJ(3) * t380;
t344 = t351 ^ 2;
t268 = qJ(3) * t331 - qJD(2) * t333 - t344 * t390 + t308;
t239 = -0.2e1 * qJD(3) * t318 + t346 * t264 + t388 * t268;
t288 = pkin(3) * t318 - qJ(4) * t319;
t353 = qJD(2) ^ 2;
t223 = -pkin(3) * t353 + qJDD(2) * qJ(4) - t288 * t318 + t239;
t336 = t349 * g(1) - t352 * g(2);
t369 = -qJDD(1) * pkin(1) - t336;
t270 = -t331 * pkin(2) + qJDD(3) + t333 * t380 + (-qJ(3) * t344 - pkin(7)) * t354 + t369;
t300 = t330 * t346 - t388 * t331;
t225 = (qJD(2) * t318 - t301) * qJ(4) + (qJD(2) * t319 + t300) * pkin(3) + t270;
t306 = t345 * qJD(2) + t387 * t319;
t216 = -t345 * t223 + t387 * t225 + t306 * t395;
t217 = t387 * t223 + t345 * t225 + t305 * t395;
t247 = Ifges(6,5) * t306 + Ifges(6,6) * t318 + Ifges(6,3) * t305;
t250 = Ifges(5,4) * t306 - Ifges(5,2) * t305 + Ifges(5,6) * t318;
t252 = Ifges(5,1) * t306 - Ifges(5,4) * t305 + Ifges(5,5) * t318;
t266 = mrSges(6,1) * t305 - mrSges(6,3) * t306;
t282 = -t387 * qJDD(2) + t345 * t301;
t265 = pkin(4) * t305 - qJ(5) * t306;
t317 = t318 ^ 2;
t214 = -t300 * pkin(4) - t317 * qJ(5) + t306 * t265 + qJDD(5) - t216;
t206 = (-t283 - t386) * pkin(8) + (t305 * t306 - t300) * pkin(5) + t214;
t391 = 2 * qJD(5);
t212 = -pkin(4) * t317 + t300 * qJ(5) - t305 * t265 + t318 * t391 + t217;
t275 = -pkin(5) * t318 - pkin(8) * t306;
t304 = t305 ^ 2;
t207 = -pkin(5) * t304 + pkin(8) * t282 + t275 * t318 + t212;
t347 = sin(qJ(6));
t350 = cos(qJ(6));
t203 = t206 * t350 - t207 * t347;
t259 = t305 * t350 - t306 * t347;
t231 = qJD(6) * t259 + t282 * t347 + t283 * t350;
t260 = t305 * t347 + t306 * t350;
t240 = -mrSges(7,1) * t259 + mrSges(7,2) * t260;
t316 = qJD(6) - t318;
t245 = -mrSges(7,2) * t316 + mrSges(7,3) * t259;
t299 = qJDD(6) - t300;
t198 = m(7) * t203 + mrSges(7,1) * t299 - mrSges(7,3) * t231 - t240 * t260 + t245 * t316;
t204 = t206 * t347 + t207 * t350;
t230 = -qJD(6) * t260 + t282 * t350 - t283 * t347;
t246 = mrSges(7,1) * t316 - mrSges(7,3) * t260;
t199 = m(7) * t204 - mrSges(7,2) * t299 + mrSges(7,3) * t230 + t240 * t259 - t246 * t316;
t188 = t350 * t198 + t347 * t199;
t251 = Ifges(6,1) * t306 + Ifges(6,4) * t318 + Ifges(6,5) * t305;
t233 = Ifges(7,4) * t260 + Ifges(7,2) * t259 + Ifges(7,6) * t316;
t234 = Ifges(7,1) * t260 + Ifges(7,4) * t259 + Ifges(7,5) * t316;
t367 = mrSges(7,1) * t203 - mrSges(7,2) * t204 + Ifges(7,5) * t231 + Ifges(7,6) * t230 + Ifges(7,3) * t299 + t260 * t233 - t259 * t234;
t359 = mrSges(6,1) * t214 - mrSges(6,3) * t212 - Ifges(6,4) * t283 - Ifges(6,2) * t300 - Ifges(6,6) * t282 + pkin(5) * t188 - t305 * t251 + t367;
t271 = -mrSges(6,2) * t305 + mrSges(6,3) * t318;
t365 = -m(6) * t214 + t300 * mrSges(6,1) + t318 * t271 - t188;
t189 = -t347 * t198 + t350 * t199;
t274 = -mrSges(6,1) * t318 + mrSges(6,2) * t306;
t368 = m(6) * t212 + t300 * mrSges(6,3) + t318 * t274 + t189;
t393 = -(t250 - t247) * t306 - mrSges(5,1) * t216 + mrSges(5,2) * t217 - Ifges(5,5) * t283 + Ifges(5,6) * t282 - pkin(4) * (-t283 * mrSges(6,2) - t306 * t266 + t365) - qJ(5) * (-t282 * mrSges(6,2) - t305 * t266 + t368) - t305 * t252 + t359;
t289 = mrSges(4,1) * t318 + mrSges(4,2) * t319;
t310 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t319;
t273 = mrSges(5,1) * t318 - mrSges(5,3) * t306;
t381 = -mrSges(5,1) * t305 - mrSges(5,2) * t306 - t266;
t389 = -mrSges(5,3) - mrSges(6,2);
t183 = m(5) * t217 - t300 * mrSges(5,2) - t318 * t273 + t389 * t282 + t381 * t305 + t368;
t272 = -mrSges(5,2) * t318 - mrSges(5,3) * t305;
t185 = m(5) * t216 + t300 * mrSges(5,1) + t318 * t272 + t389 * t283 + t381 * t306 + t365;
t373 = t387 * t183 - t185 * t345;
t177 = m(4) * t239 - qJDD(2) * mrSges(4,2) - mrSges(4,3) * t300 - qJD(2) * t310 - t289 * t318 + t373;
t378 = qJD(3) * t319;
t313 = -0.2e1 * t378;
t382 = t388 * t264 - t346 * t268;
t238 = t313 + t382;
t309 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t318;
t364 = qJDD(2) * pkin(3) + t353 * qJ(4) - t319 * t288 - qJDD(4) + t382;
t211 = -t304 * pkin(8) + t313 + (-pkin(4) - pkin(5)) * t282 - t394 + (-pkin(4) * t318 + t275 + t391) * t306 + t364;
t205 = -m(7) * t211 + t230 * mrSges(7,1) - t231 * mrSges(7,2) + t259 * t245 - t260 * t246;
t222 = 0.2e1 * t378 - t364;
t215 = -0.2e1 * qJD(5) * t306 + t394 + (t306 * t318 + t282) * pkin(4) + t222;
t200 = m(6) * t215 + mrSges(6,1) * t282 - t283 * mrSges(6,3) + t271 * t305 - t306 * t274 + t205;
t357 = -m(5) * t222 - t282 * mrSges(5,1) - mrSges(5,2) * t283 - t305 * t272 - t273 * t306 - t200;
t194 = m(4) * t238 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t301 + qJD(2) * t309 - t289 * t319 + t357;
t172 = t346 * t177 + t388 * t194;
t307 = -t351 * g(3) - t385;
t321 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t348 + Ifges(3,2) * t351) * qJD(1);
t322 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t348 + Ifges(3,4) * t351) * qJD(1);
t232 = Ifges(7,5) * t260 + Ifges(7,6) * t259 + Ifges(7,3) * t316;
t191 = -mrSges(7,1) * t211 + mrSges(7,3) * t204 + Ifges(7,4) * t231 + Ifges(7,2) * t230 + Ifges(7,6) * t299 - t232 * t260 + t234 * t316;
t192 = mrSges(7,2) * t211 - mrSges(7,3) * t203 + Ifges(7,1) * t231 + Ifges(7,4) * t230 + Ifges(7,5) * t299 + t232 * t259 - t233 * t316;
t360 = -mrSges(6,1) * t215 + mrSges(6,2) * t212 - pkin(5) * t205 - pkin(8) * t189 - t350 * t191 - t347 * t192;
t249 = Ifges(6,4) * t306 + Ifges(6,2) * t318 + Ifges(6,6) * t305;
t384 = -Ifges(5,5) * t306 + Ifges(5,6) * t305 - Ifges(5,3) * t318 - t249;
t166 = -mrSges(5,1) * t222 + mrSges(5,3) * t217 - pkin(4) * t200 + (t252 + t251) * t318 + t384 * t306 + (Ifges(5,6) - Ifges(6,6)) * t300 + (Ifges(5,4) - Ifges(6,5)) * t283 + (-Ifges(5,2) - Ifges(6,3)) * t282 + t360;
t362 = mrSges(6,2) * t214 - mrSges(6,3) * t215 + Ifges(6,1) * t283 + Ifges(6,4) * t300 + Ifges(6,5) * t282 - pkin(8) * t188 - t347 * t191 + t350 * t192 + t318 * t247;
t171 = mrSges(5,2) * t222 - mrSges(5,3) * t216 + Ifges(5,1) * t283 - Ifges(5,4) * t282 + Ifges(5,5) * t300 - qJ(5) * t200 - t318 * t250 + t384 * t305 + t362;
t285 = Ifges(4,4) * t319 - Ifges(4,2) * t318 + Ifges(4,6) * qJD(2);
t286 = Ifges(4,1) * t319 - Ifges(4,4) * t318 + Ifges(4,5) * qJD(2);
t361 = -mrSges(4,1) * t238 + mrSges(4,2) * t239 - Ifges(4,5) * t301 + Ifges(4,6) * t300 - Ifges(4,3) * qJDD(2) - pkin(3) * t357 - qJ(4) * t373 - t387 * t166 - t345 * t171 - t319 * t285 - t318 * t286;
t392 = mrSges(3,1) * t307 - mrSges(3,2) * t308 + Ifges(3,5) * t330 + Ifges(3,6) * t331 + Ifges(3,3) * qJDD(2) + pkin(2) * t172 + (t321 * t348 - t322 * t351) * qJD(1) - t361;
t179 = t345 * t183 + t387 * t185;
t329 = (-mrSges(3,1) * t351 + mrSges(3,2) * t348) * qJD(1);
t335 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t379;
t168 = m(3) * t307 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t330 + qJD(2) * t335 - t329 * t380 + t172;
t334 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t380;
t374 = t388 * t177 - t194 * t346;
t169 = m(3) * t308 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t331 - qJD(2) * t334 + t329 * t379 + t374;
t375 = -t168 * t348 + t351 * t169;
t284 = Ifges(4,5) * t319 - Ifges(4,6) * t318 + Ifges(4,3) * qJD(2);
t160 = mrSges(4,2) * t270 - mrSges(4,3) * t238 + Ifges(4,1) * t301 - Ifges(4,4) * t300 + Ifges(4,5) * qJDD(2) - qJ(4) * t179 - qJD(2) * t285 - t345 * t166 + t387 * t171 - t318 * t284;
t161 = Ifges(4,6) * qJDD(2) + (-Ifges(5,3) - Ifges(4,2)) * t300 - t319 * t284 + Ifges(4,4) * t301 + qJD(2) * t286 - mrSges(4,1) * t270 + mrSges(4,3) * t239 - pkin(3) * t179 + t393;
t320 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t348 + Ifges(3,6) * t351) * qJD(1);
t324 = -t354 * pkin(7) + t369;
t363 = m(4) * t270 + t300 * mrSges(4,1) + mrSges(4,2) * t301 + t318 * t309 + t310 * t319 + t179;
t156 = -mrSges(3,1) * t324 + mrSges(3,3) * t308 + Ifges(3,4) * t330 + Ifges(3,2) * t331 + Ifges(3,6) * qJDD(2) - pkin(2) * t363 + qJ(3) * t374 + qJD(2) * t322 + t346 * t160 + t388 * t161 - t320 * t380;
t158 = mrSges(3,2) * t324 - mrSges(3,3) * t307 + Ifges(3,1) * t330 + Ifges(3,4) * t331 + Ifges(3,5) * qJDD(2) - qJ(3) * t172 - qJD(2) * t321 + t388 * t160 - t346 * t161 + t320 * t379;
t358 = -m(3) * t324 + t331 * mrSges(3,1) - mrSges(3,2) * t330 - t334 * t380 + t335 * t379 - t363;
t366 = mrSges(2,1) * t336 - mrSges(2,2) * t337 + Ifges(2,3) * qJDD(1) + pkin(1) * t358 + pkin(7) * t375 + t351 * t156 + t348 * t158;
t173 = m(2) * t336 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t354 + t358;
t164 = t168 * t351 + t169 * t348;
t162 = m(2) * t337 - mrSges(2,1) * t354 - qJDD(1) * mrSges(2,2) + t375;
t159 = mrSges(2,1) * g(3) + mrSges(2,3) * t337 + t354 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t164 - t392;
t154 = -mrSges(2,2) * g(3) - mrSges(2,3) * t336 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t354 - pkin(7) * t164 - t156 * t348 + t158 * t351;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t352 * t154 - t349 * t159 - pkin(6) * (t162 * t349 + t173 * t352), t154, t158, t160, t171, -t249 * t305 + t362, t192; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t349 * t154 + t352 * t159 + pkin(6) * (t162 * t352 - t173 * t349), t159, t156, t161, t166, -t306 * t247 - t359, t191; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t366, t366, t392, -t361, Ifges(5,3) * t300 - t393, Ifges(6,5) * t283 + Ifges(6,6) * t300 + Ifges(6,3) * t282 + t306 * t249 - t318 * t251 - t360, t367;];
m_new  = t1;
