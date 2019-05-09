% Calculate vector of cutting torques with Newton-Euler for
% S6RRPPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
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
% Datum: 2019-05-06 08:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPPPR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:22:14
% EndTime: 2019-05-06 08:22:40
% DurationCPUTime: 14.15s
% Computational Cost: add. (221618->388), mult. (528119->470), div. (0->0), fcn. (359413->10), ass. (0->147)
t342 = sin(qJ(2));
t345 = cos(qJ(2));
t372 = qJD(1) * qJD(2);
t321 = qJDD(1) * t342 + t345 * t372;
t343 = sin(qJ(1));
t346 = cos(qJ(1));
t328 = -g(1) * t346 - g(2) * t343;
t348 = qJD(1) ^ 2;
t316 = -pkin(1) * t348 + qJDD(1) * pkin(7) + t328;
t381 = t342 * t316;
t384 = pkin(2) * t348;
t248 = qJDD(2) * pkin(2) - t321 * qJ(3) - t381 + (qJ(3) * t372 + t342 * t384 - g(3)) * t345;
t296 = -g(3) * t342 + t316 * t345;
t322 = qJDD(1) * t345 - t342 * t372;
t377 = qJD(1) * t342;
t324 = qJD(2) * pkin(2) - qJ(3) * t377;
t337 = t345 ^ 2;
t250 = qJ(3) * t322 - qJD(2) * t324 - t337 * t384 + t296;
t339 = sin(pkin(9));
t382 = cos(pkin(9));
t310 = (t339 * t345 + t342 * t382) * qJD(1);
t231 = -0.2e1 * qJD(3) * t310 + t248 * t382 - t250 * t339;
t376 = qJD(1) * t345;
t309 = t339 * t377 - t376 * t382;
t374 = qJD(3) * t309;
t304 = -0.2e1 * t374;
t380 = t248 * t339 + t250 * t382;
t232 = t304 + t380;
t266 = Ifges(4,4) * t310 - Ifges(4,2) * t309 + Ifges(4,6) * qJD(2);
t274 = -mrSges(5,2) * t309 - mrSges(5,3) * t310;
t286 = t321 * t339 - t322 * t382;
t287 = t321 * t382 + t322 * t339;
t300 = mrSges(5,1) * t309 - qJD(2) * mrSges(5,3);
t299 = pkin(4) * t310 - qJD(2) * qJ(5);
t308 = t309 ^ 2;
t272 = pkin(3) * t309 - qJ(4) * t310;
t347 = qJD(2) ^ 2;
t388 = -2 * qJD(4);
t386 = pkin(3) * t347 - qJDD(2) * qJ(4) + qJD(2) * t388 + t272 * t309 - t380;
t212 = -pkin(4) * t286 - qJ(5) * t308 + qJD(2) * t299 + qJDD(5) + t304 - t386;
t338 = sin(pkin(10));
t340 = cos(pkin(10));
t291 = -qJD(2) * t338 + t309 * t340;
t255 = -mrSges(6,2) * t310 + mrSges(6,3) * t291;
t292 = qJD(2) * t340 + t309 * t338;
t256 = mrSges(6,1) * t310 - mrSges(6,3) * t292;
t261 = -qJDD(2) * t338 + t286 * t340;
t262 = qJDD(2) * t340 + t286 * t338;
t257 = pkin(5) * t310 - pkin(8) * t292;
t290 = t291 ^ 2;
t207 = -pkin(5) * t261 - pkin(8) * t290 + t257 * t292 + t212;
t341 = sin(qJ(6));
t344 = cos(qJ(6));
t243 = t291 * t341 + t292 * t344;
t224 = -qJD(6) * t243 + t261 * t344 - t262 * t341;
t242 = t291 * t344 - t292 * t341;
t225 = qJD(6) * t242 + t261 * t341 + t262 * t344;
t307 = qJD(6) + t310;
t235 = -mrSges(7,2) * t307 + mrSges(7,3) * t242;
t236 = mrSges(7,1) * t307 - mrSges(7,3) * t243;
t362 = m(7) * t207 - t224 * mrSges(7,1) + mrSges(7,2) * t225 - t242 * t235 + t236 * t243;
t197 = -m(6) * t212 + t261 * mrSges(6,1) - mrSges(6,2) * t262 + t291 * t255 - t256 * t292 - t362;
t215 = 0.2e1 * t374 + t386;
t301 = mrSges(5,1) * t310 + qJD(2) * mrSges(5,2);
t355 = -m(5) * t215 + qJDD(2) * mrSges(5,3) + qJD(2) * t301 - t197;
t218 = -qJDD(2) * pkin(3) - t347 * qJ(4) + t272 * t310 + qJDD(4) - t231;
t375 = qJD(2) * t309;
t210 = (t309 * t310 - qJDD(2)) * qJ(5) + (t287 + t375) * pkin(4) + t218;
t327 = t343 * g(1) - g(2) * t346;
t364 = -qJDD(1) * pkin(1) - t327;
t254 = -t322 * pkin(2) + qJDD(3) + t324 * t377 + (-qJ(3) * t337 - pkin(7)) * t348 + t364;
t352 = (-t287 + t375) * qJ(4) + t254 + (pkin(3) * qJD(2) + t388) * t310;
t214 = -t308 * pkin(4) - t310 * t299 + (pkin(3) + qJ(5)) * t286 + t352;
t204 = -0.2e1 * qJD(5) * t292 + t210 * t340 - t338 * t214;
t201 = (t291 * t310 - t262) * pkin(8) + (t291 * t292 + t287) * pkin(5) + t204;
t205 = 0.2e1 * qJD(5) * t291 + t210 * t338 + t214 * t340;
t202 = -pkin(5) * t290 + pkin(8) * t261 - t257 * t310 + t205;
t200 = t201 * t341 + t202 * t344;
t226 = Ifges(7,5) * t243 + Ifges(7,6) * t242 + Ifges(7,3) * t307;
t228 = Ifges(7,1) * t243 + Ifges(7,4) * t242 + Ifges(7,5) * t307;
t285 = qJDD(6) + t287;
t185 = -mrSges(7,1) * t207 + mrSges(7,3) * t200 + Ifges(7,4) * t225 + Ifges(7,2) * t224 + Ifges(7,6) * t285 - t226 * t243 + t228 * t307;
t199 = t201 * t344 - t202 * t341;
t227 = Ifges(7,4) * t243 + Ifges(7,2) * t242 + Ifges(7,6) * t307;
t186 = mrSges(7,2) * t207 - mrSges(7,3) * t199 + Ifges(7,1) * t225 + Ifges(7,4) * t224 + Ifges(7,5) * t285 + t226 * t242 - t227 * t307;
t237 = Ifges(6,5) * t292 + Ifges(6,6) * t291 + Ifges(6,3) * t310;
t239 = Ifges(6,1) * t292 + Ifges(6,4) * t291 + Ifges(6,5) * t310;
t233 = -mrSges(7,1) * t242 + mrSges(7,2) * t243;
t193 = m(7) * t199 + mrSges(7,1) * t285 - mrSges(7,3) * t225 - t233 * t243 + t235 * t307;
t194 = m(7) * t200 - mrSges(7,2) * t285 + mrSges(7,3) * t224 + t233 * t242 - t236 * t307;
t367 = -t193 * t341 + t194 * t344;
t167 = -mrSges(6,1) * t212 + mrSges(6,3) * t205 + Ifges(6,4) * t262 + Ifges(6,2) * t261 + Ifges(6,6) * t287 - pkin(5) * t362 + pkin(8) * t367 + t344 * t185 + t341 * t186 - t292 * t237 + t310 * t239;
t184 = t193 * t344 + t194 * t341;
t238 = Ifges(6,4) * t292 + Ifges(6,2) * t291 + Ifges(6,6) * t310;
t169 = mrSges(6,2) * t212 - mrSges(6,3) * t204 + Ifges(6,1) * t262 + Ifges(6,4) * t261 + Ifges(6,5) * t287 - pkin(8) * t184 - t185 * t341 + t186 * t344 + t237 * t291 - t238 * t310;
t249 = -mrSges(6,1) * t291 + mrSges(6,2) * t292;
t181 = m(6) * t204 + mrSges(6,1) * t287 - mrSges(6,3) * t262 - t249 * t292 + t255 * t310 + t184;
t182 = m(6) * t205 - mrSges(6,2) * t287 + mrSges(6,3) * t261 + t249 * t291 - t256 * t310 + t367;
t177 = t340 * t181 + t338 * t182;
t263 = Ifges(5,5) * qJD(2) - Ifges(5,6) * t310 + Ifges(5,3) * t309;
t358 = -mrSges(5,2) * t218 + mrSges(5,3) * t215 - Ifges(5,1) * qJDD(2) + Ifges(5,4) * t287 - Ifges(5,5) * t286 + qJ(5) * t177 + t338 * t167 - t169 * t340 + t263 * t310;
t360 = -m(5) * t218 - mrSges(5,1) * t287 - t274 * t310 - t177;
t265 = Ifges(5,4) * qJD(2) - Ifges(5,2) * t310 + Ifges(5,6) * t309;
t378 = Ifges(4,1) * t310 - Ifges(4,4) * t309 + Ifges(4,5) * qJD(2) - t265;
t389 = -mrSges(4,2) * t232 + pkin(3) * (-qJDD(2) * mrSges(5,2) - qJD(2) * t300 + t360) + qJ(4) * (-t286 * mrSges(5,1) - t309 * t274 + t355) + mrSges(4,1) * t231 + t310 * t266 - Ifges(4,6) * t286 + Ifges(4,5) * t287 + Ifges(4,3) * qJDD(2) - t358 + t378 * t309;
t273 = mrSges(4,1) * t309 + mrSges(4,2) * t310;
t297 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t309;
t173 = m(4) * t231 - t287 * mrSges(4,3) - t310 * t273 + (mrSges(4,1) - mrSges(5,2)) * qJDD(2) + (t297 - t300) * qJD(2) + t360;
t298 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t310;
t189 = -qJDD(2) * mrSges(4,2) - qJD(2) * t298 + t355 + m(4) * t232 + (-t273 - t274) * t309 + (-mrSges(4,3) - mrSges(5,1)) * t286;
t166 = t173 * t382 + t189 * t339;
t295 = -t345 * g(3) - t381;
t312 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t342 + Ifges(3,2) * t345) * qJD(1);
t313 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t342 + Ifges(3,4) * t345) * qJD(1);
t387 = mrSges(3,1) * t295 - mrSges(3,2) * t296 + Ifges(3,5) * t321 + Ifges(3,6) * t322 + Ifges(3,3) * qJDD(2) + pkin(2) * t166 + (t312 * t342 - t313 * t345) * qJD(1) + t389;
t383 = -Ifges(5,6) - Ifges(4,4);
t178 = -t181 * t338 + t182 * t340;
t267 = Ifges(5,1) * qJD(2) - Ifges(5,4) * t310 + Ifges(5,5) * t309;
t379 = -Ifges(4,5) * t310 + Ifges(4,6) * t309 - Ifges(4,3) * qJD(2) - t267;
t320 = (-mrSges(3,1) * t345 + mrSges(3,2) * t342) * qJD(1);
t326 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t376;
t164 = m(3) * t295 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t321 + qJD(2) * t326 - t320 * t377 + t166;
t325 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t377;
t368 = -t173 * t339 + t189 * t382;
t165 = m(3) * t296 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t322 - qJD(2) * t325 + t320 * t376 + t368;
t369 = -t164 * t342 + t165 * t345;
t220 = t286 * pkin(3) + t352;
t174 = m(5) * t220 - mrSges(5,2) * t286 - mrSges(5,3) * t287 - t300 * t309 - t301 * t310 + t178;
t357 = -mrSges(5,1) * t215 + mrSges(5,2) * t220 - pkin(4) * t197 - qJ(5) * t178 - t340 * t167 - t338 * t169;
t158 = -mrSges(4,1) * t254 + mrSges(4,3) * t232 - pkin(3) * t174 + t379 * t310 - t383 * t287 + (-Ifges(4,2) - Ifges(5,3)) * t286 + (Ifges(4,6) - Ifges(5,5)) * qJDD(2) + t378 * qJD(2) + t357;
t359 = -mrSges(7,1) * t199 + mrSges(7,2) * t200 - Ifges(7,5) * t225 - Ifges(7,6) * t224 - Ifges(7,3) * t285 - t227 * t243 + t242 * t228;
t354 = -mrSges(6,1) * t204 + mrSges(6,2) * t205 - Ifges(6,5) * t262 - Ifges(6,6) * t261 - Ifges(6,3) * t287 - pkin(5) * t184 - t238 * t292 + t291 * t239 + t359;
t350 = -mrSges(5,1) * t218 + mrSges(5,3) * t220 - pkin(4) * t177 + t354;
t159 = -t350 + mrSges(4,2) * t254 - mrSges(4,3) * t231 - qJ(4) * t174 + t379 * t309 + (Ifges(5,2) + Ifges(4,1)) * t287 + t383 * t286 + (Ifges(4,5) - Ifges(5,4)) * qJDD(2) + (t263 - t266) * qJD(2);
t311 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t342 + Ifges(3,6) * t345) * qJD(1);
t315 = -t348 * pkin(7) + t364;
t356 = m(4) * t254 + mrSges(4,1) * t286 + mrSges(4,2) * t287 + t297 * t309 + t298 * t310 + t174;
t154 = -mrSges(3,1) * t315 + mrSges(3,3) * t296 + Ifges(3,4) * t321 + Ifges(3,2) * t322 + Ifges(3,6) * qJDD(2) - pkin(2) * t356 + qJ(3) * t368 + qJD(2) * t313 + t158 * t382 + t339 * t159 - t311 * t377;
t156 = mrSges(3,2) * t315 - mrSges(3,3) * t295 + Ifges(3,1) * t321 + Ifges(3,4) * t322 + Ifges(3,5) * qJDD(2) - qJ(3) * t166 - qJD(2) * t312 - t158 * t339 + t159 * t382 + t311 * t376;
t351 = -m(3) * t315 + mrSges(3,1) * t322 - mrSges(3,2) * t321 - t325 * t377 + t326 * t376 - t356;
t361 = mrSges(2,1) * t327 - mrSges(2,2) * t328 + Ifges(2,3) * qJDD(1) + pkin(1) * t351 + pkin(7) * t369 + t154 * t345 + t156 * t342;
t170 = m(2) * t327 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t348 + t351;
t162 = t164 * t345 + t165 * t342;
t160 = m(2) * t328 - mrSges(2,1) * t348 - qJDD(1) * mrSges(2,2) + t369;
t157 = mrSges(2,1) * g(3) + mrSges(2,3) * t328 + t348 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t162 - t387;
t152 = -mrSges(2,2) * g(3) - mrSges(2,3) * t327 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t348 - pkin(7) * t162 - t154 * t342 + t156 * t345;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t346 * t152 - t343 * t157 - pkin(6) * (t160 * t343 + t170 * t346), t152, t156, t159, -t309 * t265 - t358, t169, t186; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t343 * t152 + t346 * t157 + pkin(6) * (t160 * t346 - t170 * t343), t157, t154, t158, Ifges(5,4) * qJDD(2) - Ifges(5,2) * t287 + Ifges(5,6) * t286 - qJD(2) * t263 + t309 * t267 + t350, t167, t185; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t361, t361, t387, t389, Ifges(5,5) * qJDD(2) - Ifges(5,6) * t287 + Ifges(5,3) * t286 + qJD(2) * t265 + t310 * t267 - t357, -t354, -t359;];
m_new  = t1;
