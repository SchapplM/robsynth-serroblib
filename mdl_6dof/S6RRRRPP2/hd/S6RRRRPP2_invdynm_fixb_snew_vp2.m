% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-05-07 18:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRPP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:01:20
% EndTime: 2019-05-07 18:01:40
% DurationCPUTime: 8.32s
% Computational Cost: add. (138099->386), mult. (276535->448), div. (0->0), fcn. (190787->8), ass. (0->142)
t342 = sin(qJ(2));
t344 = cos(qJ(2));
t371 = qJD(1) * qJD(2);
t323 = qJDD(1) * t342 + t344 * t371;
t343 = sin(qJ(1));
t345 = cos(qJ(1));
t329 = -g(1) * t345 - g(2) * t343;
t346 = qJD(1) ^ 2;
t317 = -pkin(1) * t346 + qJDD(1) * pkin(7) + t329;
t378 = t342 * t317;
t384 = pkin(2) * t346;
t269 = qJDD(2) * pkin(2) - t323 * pkin(8) - t378 + (pkin(8) * t371 + t342 * t384 - g(3)) * t344;
t304 = -t342 * g(3) + t344 * t317;
t373 = qJD(1) * t342;
t327 = qJD(2) * pkin(2) - pkin(8) * t373;
t339 = t344 ^ 2;
t370 = qJDD(1) * t344;
t361 = -t342 * t371 + t370;
t270 = pkin(8) * t361 - qJD(2) * t327 - t339 * t384 + t304;
t341 = sin(qJ(3));
t387 = cos(qJ(3));
t216 = t341 * t269 + t387 * t270;
t315 = (t341 * t344 + t342 * t387) * qJD(1);
t284 = -qJD(3) * t315 - t323 * t341 + t387 * t361;
t372 = qJD(1) * t344;
t314 = -t341 * t373 + t387 * t372;
t298 = -mrSges(4,1) * t314 + mrSges(4,2) * t315;
t337 = qJD(2) + qJD(3);
t306 = mrSges(4,1) * t337 - mrSges(4,3) * t315;
t336 = qJDD(2) + qJDD(3);
t285 = t314 * qJD(3) + t323 * t387 + t341 * t361;
t328 = t343 * g(1) - t345 * g(2);
t362 = -qJDD(1) * pkin(1) - t328;
t286 = -t361 * pkin(2) + t327 * t373 + (-pkin(8) * t339 - pkin(7)) * t346 + t362;
t209 = (-t314 * t337 - t285) * pkin(9) + (t315 * t337 - t284) * pkin(3) + t286;
t299 = -pkin(3) * t314 - pkin(9) * t315;
t335 = t337 ^ 2;
t213 = -pkin(3) * t335 + pkin(9) * t336 + t299 * t314 + t216;
t340 = sin(qJ(4));
t386 = cos(qJ(4));
t206 = t209 * t386 - t340 * t213;
t301 = t315 * t340 - t337 * t386;
t240 = -t301 * qJD(4) + t285 * t386 + t340 * t336;
t283 = qJDD(4) - t284;
t310 = qJD(4) - t314;
t288 = mrSges(7,2) * t310 + mrSges(7,3) * t301;
t289 = -mrSges(5,2) * t310 - mrSges(5,3) * t301;
t302 = t315 * t386 + t340 * t337;
t264 = pkin(4) * t301 - qJ(5) * t302;
t309 = t310 ^ 2;
t204 = -t283 * pkin(4) - t309 * qJ(5) + t302 * t264 + qJDD(5) - t206;
t380 = t301 * t310;
t389 = -0.2e1 * t302;
t192 = qJD(6) * t389 + (-t240 - t380) * qJ(6) + (t301 * t302 - t283) * pkin(5) + t204;
t266 = -mrSges(7,1) * t301 + mrSges(7,2) * t302;
t363 = -m(7) * t192 + t240 * mrSges(7,3) + t302 * t266;
t265 = mrSges(6,1) * t301 - mrSges(6,3) * t302;
t374 = -mrSges(5,1) * t301 - mrSges(5,2) * t302 - t265;
t383 = -mrSges(5,3) - mrSges(6,2);
t287 = -mrSges(6,2) * t301 + mrSges(6,3) * t310;
t393 = -m(6) * t204 + t283 * mrSges(6,1) + t310 * t287;
t177 = m(5) * t206 + (t288 + t289) * t310 + t374 * t302 + (mrSges(5,1) + mrSges(7,1)) * t283 + t383 * t240 + t363 + t393;
t207 = t340 * t209 + t386 * t213;
t239 = qJD(4) * t302 + t285 * t340 - t336 * t386;
t291 = -mrSges(7,1) * t310 - mrSges(7,3) * t302;
t292 = mrSges(5,1) * t310 - mrSges(5,3) * t302;
t388 = 2 * qJD(5);
t202 = -pkin(4) * t309 + t283 * qJ(5) - t301 * t264 + t310 * t388 + t207;
t290 = -pkin(5) * t310 - qJ(6) * t302;
t300 = t301 ^ 2;
t196 = -pkin(5) * t300 + qJ(6) * t239 + 0.2e1 * qJD(6) * t301 + t290 * t310 + t202;
t369 = m(7) * t196 + t239 * mrSges(7,3) + t301 * t266;
t293 = -mrSges(6,1) * t310 + mrSges(6,2) * t302;
t394 = m(6) * t202 + t283 * mrSges(6,3) + t310 * t293;
t182 = m(5) * t207 + (t291 - t292) * t310 + t374 * t301 + (-mrSges(5,2) + mrSges(7,2)) * t283 + t383 * t239 + t369 + t394;
t365 = -t177 * t340 + t386 * t182;
t172 = m(4) * t216 - mrSges(4,2) * t336 + mrSges(4,3) * t284 + t298 * t314 - t306 * t337 + t365;
t215 = t387 * t269 - t341 * t270;
t305 = -mrSges(4,2) * t337 + mrSges(4,3) * t314;
t357 = t336 * pkin(3) + t335 * pkin(9) - t315 * t299 + t215;
t392 = (-t240 + t380) * qJ(5) - t357;
t199 = -t300 * qJ(6) + qJDD(6) + (-pkin(4) - pkin(5)) * t239 + (-pkin(4) * t310 + t290 + t388) * t302 - t392;
t190 = -m(7) * t199 + t239 * mrSges(7,1) - t240 * mrSges(7,2) + t301 * t288 - t302 * t291;
t205 = qJD(5) * t389 + (t302 * t310 + t239) * pkin(4) + t392;
t187 = m(6) * t205 + mrSges(6,1) * t239 - t240 * mrSges(6,3) + t287 * t301 - t302 * t293 + t190;
t349 = m(5) * t357 - t239 * mrSges(5,1) - mrSges(5,2) * t240 - t301 * t289 - t292 * t302 - t187;
t179 = m(4) * t215 + mrSges(4,1) * t336 - mrSges(4,3) * t285 - t298 * t315 + t305 * t337 + t349;
t165 = t341 * t172 + t387 * t179;
t303 = -t344 * g(3) - t378;
t381 = Ifges(3,6) * qJD(2);
t312 = t381 + (Ifges(3,4) * t342 + Ifges(3,2) * t344) * qJD(1);
t313 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t342 + Ifges(3,4) * t344) * qJD(1);
t245 = Ifges(7,5) * t302 + Ifges(7,6) * t301 - Ifges(7,3) * t310;
t252 = Ifges(6,1) * t302 + Ifges(6,4) * t310 + Ifges(6,5) * t301;
t253 = Ifges(5,1) * t302 - Ifges(5,4) * t301 + Ifges(5,5) * t310;
t189 = t283 * mrSges(7,2) + t310 * t291 + t369;
t251 = Ifges(7,1) * t302 + Ifges(7,4) * t301 - Ifges(7,5) * t310;
t359 = mrSges(7,1) * t199 - mrSges(7,3) * t196 - Ifges(7,4) * t240 - Ifges(7,2) * t239 + Ifges(7,6) * t283 + t310 * t251;
t352 = mrSges(6,1) * t205 - mrSges(6,2) * t202 + pkin(5) * t190 + qJ(6) * t189 - t359;
t249 = Ifges(6,4) * t302 + Ifges(6,2) * t310 + Ifges(6,6) * t301;
t376 = -Ifges(5,5) * t302 + Ifges(5,6) * t301 - Ifges(5,3) * t310 - t249;
t161 = -t352 + (t253 + t252) * t310 + (t245 + t376) * t302 + (Ifges(5,6) - Ifges(6,6)) * t283 + (Ifges(5,4) - Ifges(6,5)) * t240 + (-Ifges(5,2) - Ifges(6,3)) * t239 - pkin(4) * t187 + mrSges(5,3) * t207 + mrSges(5,1) * t357;
t248 = Ifges(7,4) * t302 + Ifges(7,2) * t301 - Ifges(7,6) * t310;
t250 = Ifges(5,4) * t302 - Ifges(5,2) * t301 + Ifges(5,6) * t310;
t188 = -t283 * mrSges(7,1) - t310 * t288 - t363;
t246 = Ifges(6,5) * t302 + Ifges(6,6) * t310 + Ifges(6,3) * t301;
t358 = mrSges(7,2) * t199 - mrSges(7,3) * t192 + Ifges(7,1) * t240 + Ifges(7,4) * t239 - Ifges(7,5) * t283 + t301 * t245;
t351 = mrSges(6,2) * t204 - mrSges(6,3) * t205 + Ifges(6,1) * t240 + Ifges(6,4) * t283 + Ifges(6,5) * t239 - qJ(6) * t188 + t310 * t246 + t358;
t167 = t351 + (-t250 + t248) * t310 + t376 * t301 - qJ(5) * t187 - mrSges(5,3) * t206 - mrSges(5,2) * t357 - Ifges(5,4) * t239 + Ifges(5,1) * t240 + Ifges(5,5) * t283;
t295 = Ifges(4,4) * t315 + Ifges(4,2) * t314 + Ifges(4,6) * t337;
t296 = Ifges(4,1) * t315 + Ifges(4,4) * t314 + Ifges(4,5) * t337;
t354 = -mrSges(4,1) * t215 + mrSges(4,2) * t216 - Ifges(4,5) * t285 - Ifges(4,6) * t284 - Ifges(4,3) * t336 - pkin(3) * t349 - pkin(9) * t365 - t386 * t161 - t340 * t167 - t315 * t295 + t314 * t296;
t395 = -mrSges(3,1) * t303 + mrSges(3,2) * t304 - Ifges(3,5) * t323 - Ifges(3,3) * qJDD(2) - pkin(2) * t165 + (t342 * (-t312 + t381) + t344 * t313) * qJD(1) + t354;
t360 = mrSges(7,1) * t192 - mrSges(7,2) * t196 + Ifges(7,5) * t240 + Ifges(7,6) * t239 - Ifges(7,3) * t283 + t302 * t248 - t301 * t251;
t353 = mrSges(6,1) * t204 - mrSges(6,3) * t202 - Ifges(6,4) * t240 - Ifges(6,2) * t283 - Ifges(6,6) * t239 + pkin(5) * t188 - t301 * t252 + t360;
t391 = (t250 - t246) * t302 + mrSges(5,1) * t206 - mrSges(5,2) * t207 + Ifges(5,5) * t240 - Ifges(5,6) * t239 + Ifges(5,3) * t283 + pkin(4) * (-t240 * mrSges(6,2) - t302 * t265 - t188 + t393) + qJ(5) * (-t239 * mrSges(6,2) - t301 * t265 + t189 + t394) + t301 * t253 - t353;
t382 = Ifges(3,6) * t344;
t379 = t310 * t248;
t174 = t386 * t177 + t340 * t182;
t322 = (-mrSges(3,1) * t344 + mrSges(3,2) * t342) * qJD(1);
t326 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t372;
t163 = m(3) * t303 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t323 + qJD(2) * t326 - t322 * t373 + t165;
t325 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t373;
t366 = t387 * t172 - t341 * t179;
t164 = m(3) * t304 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t361 - qJD(2) * t325 + t322 * t372 + t366;
t367 = -t163 * t342 + t344 * t164;
t294 = Ifges(4,5) * t315 + Ifges(4,6) * t314 + Ifges(4,3) * t337;
t155 = mrSges(4,2) * t286 - mrSges(4,3) * t215 + Ifges(4,1) * t285 + Ifges(4,4) * t284 + Ifges(4,5) * t336 - pkin(9) * t174 - t340 * t161 + t167 * t386 + t314 * t294 - t337 * t295;
t159 = -mrSges(4,1) * t286 + mrSges(4,3) * t216 + Ifges(4,4) * t285 + Ifges(4,2) * t284 + Ifges(4,6) * t336 - pkin(3) * t174 - t315 * t294 + t337 * t296 - t391;
t311 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t342 + t382) * qJD(1);
t316 = -t346 * pkin(7) + t362;
t355 = m(4) * t286 - t284 * mrSges(4,1) + mrSges(4,2) * t285 - t314 * t305 + t306 * t315 + t174;
t151 = -mrSges(3,1) * t316 + mrSges(3,3) * t304 + Ifges(3,4) * t323 + Ifges(3,2) * t361 + Ifges(3,6) * qJDD(2) - pkin(2) * t355 + pkin(8) * t366 + qJD(2) * t313 + t341 * t155 + t159 * t387 - t311 * t373;
t154 = mrSges(3,2) * t316 - mrSges(3,3) * t303 + Ifges(3,1) * t323 + Ifges(3,4) * t361 + Ifges(3,5) * qJDD(2) - pkin(8) * t165 - qJD(2) * t312 + t155 * t387 - t341 * t159 + t311 * t372;
t350 = -m(3) * t316 + t361 * mrSges(3,1) - mrSges(3,2) * t323 - t325 * t373 + t326 * t372 - t355;
t356 = mrSges(2,1) * t328 - mrSges(2,2) * t329 + Ifges(2,3) * qJDD(1) + pkin(1) * t350 + pkin(7) * t367 + t344 * t151 + t342 * t154;
t168 = m(2) * t328 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t346 + t350;
t158 = t163 * t344 + t164 * t342;
t156 = m(2) * t329 - mrSges(2,1) * t346 - qJDD(1) * mrSges(2,2) + t367;
t152 = -pkin(1) * t158 + mrSges(2,1) * g(3) + (Ifges(2,6) - t382) * qJDD(1) + mrSges(2,3) * t329 + t346 * Ifges(2,5) + t395;
t149 = -mrSges(2,2) * g(3) - mrSges(2,3) * t328 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t346 - pkin(7) * t158 - t151 * t342 + t154 * t344;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t345 * t149 - t343 * t152 - pkin(6) * (t156 * t343 + t168 * t345), t149, t154, t155, t167, -t301 * t249 + t351 + t379, t358 + t379; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t343 * t149 + t345 * t152 + pkin(6) * (t156 * t345 - t168 * t343), t152, t151, t159, t161, -t302 * t246 - t353, -t302 * t245 - t359; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t356, t356, Ifges(3,6) * t370 - t395, -t354, t391, (-t245 + t249) * t302 + t352 + Ifges(6,3) * t239 + Ifges(6,5) * t240 + Ifges(6,6) * t283 - t310 * t252, t360;];
m_new  = t1;
