% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 17:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRRP3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:29:32
% EndTime: 2019-05-06 17:30:13
% DurationCPUTime: 19.64s
% Computational Cost: add. (326690->380), mult. (745220->469), div. (0->0), fcn. (536410->10), ass. (0->142)
t341 = sin(qJ(2));
t345 = cos(qJ(2));
t372 = qJD(1) * qJD(2);
t324 = qJDD(1) * t341 + t345 * t372;
t325 = qJDD(1) * t345 - t341 * t372;
t337 = sin(pkin(10));
t338 = cos(pkin(10));
t299 = t324 * t338 + t325 * t337;
t313 = (t337 * t345 + t338 * t341) * qJD(1);
t340 = sin(qJ(4));
t344 = cos(qJ(4));
t302 = qJD(2) * t340 + t313 * t344;
t265 = -qJD(4) * t302 + qJDD(2) * t344 - t299 * t340;
t301 = qJD(2) * t344 - t313 * t340;
t266 = qJD(4) * t301 + qJDD(2) * t340 + t299 * t344;
t339 = sin(qJ(5));
t343 = cos(qJ(5));
t271 = t301 * t343 - t302 * t339;
t227 = qJD(5) * t271 + t265 * t339 + t266 * t343;
t272 = t301 * t339 + t302 * t343;
t251 = -mrSges(7,1) * t271 + mrSges(7,2) * t272;
t342 = sin(qJ(1));
t346 = cos(qJ(1));
t331 = -g(1) * t346 - g(2) * t342;
t348 = qJD(1) ^ 2;
t319 = -pkin(1) * t348 + qJDD(1) * pkin(7) + t331;
t377 = t319 * t341;
t378 = pkin(2) * t348;
t274 = qJDD(2) * pkin(2) - qJ(3) * t324 - t377 + (qJ(3) * t372 + t341 * t378 - g(3)) * t345;
t304 = -g(3) * t341 + t345 * t319;
t375 = qJD(1) * t341;
t327 = qJD(2) * pkin(2) - qJ(3) * t375;
t336 = t345 ^ 2;
t275 = qJ(3) * t325 - qJD(2) * t327 - t336 * t378 + t304;
t374 = qJD(1) * t345;
t312 = -t337 * t375 + t338 * t374;
t250 = 0.2e1 * qJD(3) * t312 + t337 * t274 + t338 * t275;
t292 = -pkin(3) * t312 - pkin(8) * t313;
t347 = qJD(2) ^ 2;
t233 = -pkin(3) * t347 + qJDD(2) * pkin(8) + t292 * t312 + t250;
t330 = g(1) * t342 - t346 * g(2);
t361 = -qJDD(1) * pkin(1) - t330;
t278 = -pkin(2) * t325 + qJDD(3) + t327 * t375 + (-qJ(3) * t336 - pkin(7)) * t348 + t361;
t298 = -t337 * t324 + t325 * t338;
t238 = (-qJD(2) * t312 - t299) * pkin(8) + (qJD(2) * t313 - t298) * pkin(3) + t278;
t212 = -t233 * t340 + t344 * t238;
t297 = qJDD(4) - t298;
t311 = qJD(4) - t312;
t208 = (t301 * t311 - t266) * pkin(9) + (t301 * t302 + t297) * pkin(4) + t212;
t213 = t344 * t233 + t340 * t238;
t281 = pkin(4) * t311 - pkin(9) * t302;
t300 = t301 ^ 2;
t210 = -pkin(4) * t300 + pkin(9) * t265 - t281 * t311 + t213;
t201 = t343 * t208 - t210 * t339;
t293 = qJDD(5) + t297;
t307 = qJD(5) + t311;
t196 = -0.2e1 * qJD(6) * t272 + (t271 * t307 - t227) * qJ(6) + (t271 * t272 + t293) * pkin(5) + t201;
t254 = -mrSges(7,2) * t307 + mrSges(7,3) * t271;
t371 = m(7) * t196 + t293 * mrSges(7,1) + t307 * t254;
t193 = -mrSges(7,3) * t227 - t251 * t272 + t371;
t202 = t339 * t208 + t343 * t210;
t226 = -qJD(5) * t272 + t265 * t343 - t266 * t339;
t242 = Ifges(6,4) * t272 + Ifges(6,2) * t271 + Ifges(6,6) * t307;
t243 = Ifges(7,1) * t272 + Ifges(7,4) * t271 + Ifges(7,5) * t307;
t244 = Ifges(6,1) * t272 + Ifges(6,4) * t271 + Ifges(6,5) * t307;
t256 = pkin(5) * t307 - qJ(6) * t272;
t270 = t271 ^ 2;
t199 = -pkin(5) * t270 + qJ(6) * t226 + 0.2e1 * qJD(6) * t271 - t256 * t307 + t202;
t241 = Ifges(7,4) * t272 + Ifges(7,2) * t271 + Ifges(7,6) * t307;
t359 = -mrSges(7,1) * t196 + mrSges(7,2) * t199 - Ifges(7,5) * t227 - Ifges(7,6) * t226 - Ifges(7,3) * t293 - t272 * t241;
t381 = mrSges(6,1) * t201 - mrSges(6,2) * t202 + Ifges(6,5) * t227 + Ifges(6,6) * t226 + Ifges(6,3) * t293 + pkin(5) * t193 + t272 * t242 - t359 + (-t244 - t243) * t271;
t252 = -mrSges(6,1) * t271 + mrSges(6,2) * t272;
t255 = -mrSges(6,2) * t307 + mrSges(6,3) * t271;
t184 = m(6) * t201 + mrSges(6,1) * t293 + t255 * t307 + (-t251 - t252) * t272 + (-mrSges(6,3) - mrSges(7,3)) * t227 + t371;
t257 = mrSges(7,1) * t307 - mrSges(7,3) * t272;
t258 = mrSges(6,1) * t307 - mrSges(6,3) * t272;
t370 = m(7) * t199 + t226 * mrSges(7,3) + t271 * t251;
t187 = m(6) * t202 + mrSges(6,3) * t226 + t252 * t271 + (-t257 - t258) * t307 + (-mrSges(6,2) - mrSges(7,2)) * t293 + t370;
t182 = t343 * t184 + t339 * t187;
t260 = Ifges(5,4) * t302 + Ifges(5,2) * t301 + Ifges(5,6) * t311;
t261 = Ifges(5,1) * t302 + Ifges(5,4) * t301 + Ifges(5,5) * t311;
t380 = mrSges(5,1) * t212 - mrSges(5,2) * t213 + Ifges(5,5) * t266 + Ifges(5,6) * t265 + Ifges(5,3) * t297 + pkin(4) * t182 + t302 * t260 - t301 * t261 + t381;
t249 = -0.2e1 * qJD(3) * t313 + t274 * t338 - t337 * t275;
t286 = -mrSges(4,1) * t312 + mrSges(4,2) * t313;
t306 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t313;
t276 = -mrSges(5,1) * t301 + mrSges(5,2) * t302;
t279 = -mrSges(5,2) * t311 + mrSges(5,3) * t301;
t179 = m(5) * t212 + mrSges(5,1) * t297 - mrSges(5,3) * t266 - t276 * t302 + t279 * t311 + t182;
t280 = mrSges(5,1) * t311 - mrSges(5,3) * t302;
t365 = -t184 * t339 + t343 * t187;
t180 = m(5) * t213 - mrSges(5,2) * t297 + mrSges(5,3) * t265 + t276 * t301 - t280 * t311 + t365;
t366 = -t179 * t340 + t344 * t180;
t171 = m(4) * t250 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t298 - qJD(2) * t306 + t286 * t312 + t366;
t305 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t312;
t232 = -qJDD(2) * pkin(3) - pkin(8) * t347 + t313 * t292 - t249;
t211 = -pkin(4) * t265 - pkin(9) * t300 + t302 * t281 + t232;
t205 = -pkin(5) * t226 - qJ(6) * t270 + t256 * t272 + qJDD(6) + t211;
t363 = m(7) * t205 - t226 * mrSges(7,1) + t227 * mrSges(7,2) - t271 * t254 + t272 * t257;
t354 = m(6) * t211 - t226 * mrSges(6,1) + mrSges(6,2) * t227 - t271 * t255 + t258 * t272 + t363;
t351 = -m(5) * t232 + t265 * mrSges(5,1) - mrSges(5,2) * t266 + t301 * t279 - t280 * t302 - t354;
t189 = m(4) * t249 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t299 + qJD(2) * t305 - t286 * t313 + t351;
t166 = t337 * t171 + t338 * t189;
t303 = -g(3) * t345 - t377;
t315 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t341 + Ifges(3,2) * t345) * qJD(1);
t316 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t341 + Ifges(3,4) * t345) * qJD(1);
t239 = Ifges(7,5) * t272 + Ifges(7,6) * t271 + Ifges(7,3) * t307;
t240 = Ifges(6,5) * t272 + Ifges(6,6) * t271 + Ifges(6,3) * t307;
t360 = -mrSges(7,1) * t205 + mrSges(7,3) * t199 + Ifges(7,4) * t227 + Ifges(7,2) * t226 + Ifges(7,6) * t293 + t307 * t243;
t175 = Ifges(6,4) * t227 + Ifges(6,2) * t226 + Ifges(6,6) * t293 + t307 * t244 - mrSges(6,1) * t211 + mrSges(6,3) * t202 - pkin(5) * t363 + qJ(6) * (-mrSges(7,2) * t293 - t257 * t307 + t370) + (-t240 - t239) * t272 + t360;
t358 = mrSges(7,2) * t205 - mrSges(7,3) * t196 + Ifges(7,1) * t227 + Ifges(7,4) * t226 + Ifges(7,5) * t293 + t271 * t239;
t181 = mrSges(6,2) * t211 - mrSges(6,3) * t201 + Ifges(6,1) * t227 + Ifges(6,4) * t226 + Ifges(6,5) * t293 - qJ(6) * t193 + t240 * t271 + (-t241 - t242) * t307 + t358;
t259 = Ifges(5,5) * t302 + Ifges(5,6) * t301 + Ifges(5,3) * t311;
t160 = -mrSges(5,1) * t232 + mrSges(5,3) * t213 + Ifges(5,4) * t266 + Ifges(5,2) * t265 + Ifges(5,6) * t297 - pkin(4) * t354 + pkin(9) * t365 + t343 * t175 + t339 * t181 - t302 * t259 + t311 * t261;
t162 = mrSges(5,2) * t232 - mrSges(5,3) * t212 + Ifges(5,1) * t266 + Ifges(5,4) * t265 + Ifges(5,5) * t297 - pkin(9) * t182 - t175 * t339 + t181 * t343 + t259 * t301 - t260 * t311;
t283 = Ifges(4,4) * t313 + Ifges(4,2) * t312 + Ifges(4,6) * qJD(2);
t284 = Ifges(4,1) * t313 + Ifges(4,4) * t312 + Ifges(4,5) * qJD(2);
t355 = -mrSges(4,1) * t249 + mrSges(4,2) * t250 - Ifges(4,5) * t299 - Ifges(4,6) * t298 - Ifges(4,3) * qJDD(2) - pkin(3) * t351 - pkin(8) * t366 - t344 * t160 - t340 * t162 - t313 * t283 + t312 * t284;
t379 = mrSges(3,1) * t303 - mrSges(3,2) * t304 + Ifges(3,5) * t324 + Ifges(3,6) * t325 + Ifges(3,3) * qJDD(2) + pkin(2) * t166 + (t341 * t315 - t345 * t316) * qJD(1) - t355;
t173 = t344 * t179 + t340 * t180;
t323 = (-mrSges(3,1) * t345 + mrSges(3,2) * t341) * qJD(1);
t329 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t374;
t164 = m(3) * t303 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t324 + qJD(2) * t329 - t323 * t375 + t166;
t328 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t375;
t367 = t338 * t171 - t337 * t189;
t165 = m(3) * t304 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t325 - qJD(2) * t328 + t323 * t374 + t367;
t368 = -t164 * t341 + t345 * t165;
t282 = Ifges(4,5) * t313 + Ifges(4,6) * t312 + Ifges(4,3) * qJD(2);
t154 = mrSges(4,2) * t278 - mrSges(4,3) * t249 + Ifges(4,1) * t299 + Ifges(4,4) * t298 + Ifges(4,5) * qJDD(2) - pkin(8) * t173 - qJD(2) * t283 - t160 * t340 + t162 * t344 + t282 * t312;
t158 = -mrSges(4,1) * t278 + mrSges(4,3) * t250 + Ifges(4,4) * t299 + Ifges(4,2) * t298 + Ifges(4,6) * qJDD(2) - pkin(3) * t173 + qJD(2) * t284 - t313 * t282 - t380;
t314 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t341 + Ifges(3,6) * t345) * qJD(1);
t318 = -pkin(7) * t348 + t361;
t356 = m(4) * t278 - t298 * mrSges(4,1) + mrSges(4,2) * t299 - t312 * t305 + t306 * t313 + t173;
t150 = -mrSges(3,1) * t318 + mrSges(3,3) * t304 + Ifges(3,4) * t324 + Ifges(3,2) * t325 + Ifges(3,6) * qJDD(2) - pkin(2) * t356 + qJ(3) * t367 + qJD(2) * t316 + t337 * t154 + t338 * t158 - t314 * t375;
t153 = mrSges(3,2) * t318 - mrSges(3,3) * t303 + Ifges(3,1) * t324 + Ifges(3,4) * t325 + Ifges(3,5) * qJDD(2) - qJ(3) * t166 - qJD(2) * t315 + t154 * t338 - t158 * t337 + t314 * t374;
t352 = -m(3) * t318 + t325 * mrSges(3,1) - mrSges(3,2) * t324 - t328 * t375 + t329 * t374 - t356;
t357 = mrSges(2,1) * t330 - mrSges(2,2) * t331 + Ifges(2,3) * qJDD(1) + pkin(1) * t352 + pkin(7) * t368 + t345 * t150 + t341 * t153;
t167 = m(2) * t330 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t348 + t352;
t157 = t164 * t345 + t165 * t341;
t155 = m(2) * t331 - mrSges(2,1) * t348 - qJDD(1) * mrSges(2,2) + t368;
t151 = mrSges(2,1) * g(3) + mrSges(2,3) * t331 + t348 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t157 - t379;
t148 = -mrSges(2,2) * g(3) - mrSges(2,3) * t330 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t348 - pkin(7) * t157 - t150 * t341 + t153 * t345;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t346 * t148 - t342 * t151 - pkin(6) * (t155 * t342 + t167 * t346), t148, t153, t154, t162, t181, -t241 * t307 + t358; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t342 * t148 + t346 * t151 + pkin(6) * (t155 * t346 - t342 * t167), t151, t150, t158, t160, t175, -t272 * t239 + t360; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t357, t357, t379, -t355, t380, t381, -t271 * t243 - t359;];
m_new  = t1;
