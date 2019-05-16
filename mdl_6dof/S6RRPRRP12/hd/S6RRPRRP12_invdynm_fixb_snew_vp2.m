% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRRP12
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-05-06 19:04
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRRP12_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP12_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP12_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP12_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP12_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP12_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP12_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:59:23
% EndTime: 2019-05-06 18:59:42
% DurationCPUTime: 7.94s
% Computational Cost: add. (125249->387), mult. (255442->452), div. (0->0), fcn. (156292->8), ass. (0->137)
t379 = -2 * qJD(3);
t335 = sin(qJ(1));
t338 = cos(qJ(1));
t316 = -g(1) * t338 - g(2) * t335;
t340 = qJD(1) ^ 2;
t288 = -pkin(1) * t340 + qJDD(1) * pkin(7) + t316;
t334 = sin(qJ(2));
t337 = cos(qJ(2));
t266 = -t337 * g(3) - t334 * t288;
t267 = -t334 * g(3) + t337 * t288;
t283 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t334 + Ifges(3,4) * t337) * qJD(1);
t303 = (mrSges(4,2) * t337 - mrSges(4,3) * t334) * qJD(1);
t366 = qJD(1) * qJD(2);
t363 = t337 * t366;
t305 = qJDD(1) * t334 + t363;
t362 = t334 * t366;
t306 = qJDD(1) * t337 - t362;
t367 = qJD(1) * t337;
t312 = -mrSges(4,1) * t367 - qJD(2) * mrSges(4,3);
t323 = t334 * qJD(1);
t302 = (-pkin(2) * t337 - qJ(3) * t334) * qJD(1);
t339 = qJD(2) ^ 2;
t243 = t339 * pkin(2) - qJDD(2) * qJ(3) + qJD(2) * t379 - t302 * t367 - t267;
t314 = pkin(3) * t323 - qJD(2) * pkin(8);
t331 = t337 ^ 2;
t226 = -t331 * t340 * pkin(8) + t306 * pkin(3) + qJD(2) * t314 - t243;
t333 = sin(qJ(4));
t336 = cos(qJ(4));
t301 = qJD(2) * t336 - t333 * t367;
t259 = -qJD(4) * t301 - qJDD(2) * t333 - t306 * t336;
t300 = -qJD(2) * t333 - t336 * t367;
t260 = qJD(4) * t300 + qJDD(2) * t336 - t306 * t333;
t320 = t323 + qJD(4);
t264 = -mrSges(5,2) * t320 + mrSges(5,3) * t300;
t265 = mrSges(5,1) * t320 - mrSges(5,3) * t301;
t268 = pkin(4) * t320 - pkin(9) * t301;
t298 = t300 ^ 2;
t207 = -t259 * pkin(4) - t298 * pkin(9) + t301 * t268 + t226;
t332 = sin(qJ(5));
t376 = cos(qJ(5));
t262 = t332 * t300 + t376 * t301;
t218 = t262 * qJD(5) - t376 * t259 + t332 * t260;
t261 = -t376 * t300 + t332 * t301;
t219 = -t261 * qJD(5) + t332 * t259 + t376 * t260;
t318 = qJD(5) + t320;
t198 = (t262 * t318 + t218) * pkin(5) + (t261 * t318 - t219) * qJ(6) + t207 - 0.2e1 * qJD(6) * t262;
t246 = -mrSges(7,2) * t261 + mrSges(7,3) * t318;
t249 = -mrSges(7,1) * t318 + mrSges(7,2) * t262;
t188 = m(7) * t198 + t218 * mrSges(7,1) - t219 * mrSges(7,3) + t261 * t246 - t262 * t249;
t247 = -mrSges(6,2) * t318 - mrSges(6,3) * t261;
t248 = mrSges(6,1) * t318 - mrSges(6,3) * t262;
t347 = m(6) * t207 + t218 * mrSges(6,1) + t219 * mrSges(6,2) + t261 * t247 + t262 * t248 + t188;
t183 = -m(5) * t226 + t259 * mrSges(5,1) - t260 * mrSges(5,2) + t300 * t264 - t301 * t265 - t347;
t313 = mrSges(4,1) * t323 + qJD(2) * mrSges(4,2);
t344 = -m(4) * t243 + qJDD(2) * mrSges(4,3) + qJD(2) * t313 + t303 * t367 - t183;
t315 = t335 * g(1) - t338 * g(2);
t357 = -qJDD(1) * pkin(1) - t315;
t348 = pkin(2) * t362 + t323 * t379 + (-t305 - t363) * qJ(3) + t357;
t221 = -t314 * t323 + (-pkin(3) * t331 - pkin(7)) * t340 + (-pkin(2) - pkin(8)) * t306 + t348;
t245 = -qJDD(2) * pkin(2) - t339 * qJ(3) + t302 * t323 + qJDD(3) - t266;
t227 = (-t334 * t337 * t340 - qJDD(2)) * pkin(8) + (t305 - t363) * pkin(3) + t245;
t204 = -t333 * t221 + t336 * t227;
t299 = qJDD(4) + t305;
t200 = (t300 * t320 - t260) * pkin(9) + (t300 * t301 + t299) * pkin(4) + t204;
t205 = t336 * t221 + t333 * t227;
t202 = -pkin(4) * t298 + pkin(9) * t259 - t268 * t320 + t205;
t196 = t332 * t200 + t376 * t202;
t232 = Ifges(7,1) * t262 + Ifges(7,4) * t318 + Ifges(7,5) * t261;
t233 = Ifges(6,1) * t262 - Ifges(6,4) * t261 + Ifges(6,5) * t318;
t290 = qJDD(5) + t299;
t236 = pkin(5) * t261 - qJ(6) * t262;
t317 = t318 ^ 2;
t191 = -pkin(5) * t317 + qJ(6) * t290 + 0.2e1 * qJD(6) * t318 - t236 * t261 + t196;
t358 = -mrSges(7,1) * t198 + mrSges(7,2) * t191;
t230 = Ifges(7,4) * t262 + Ifges(7,2) * t318 + Ifges(7,6) * t261;
t371 = -Ifges(6,5) * t262 + Ifges(6,6) * t261 - Ifges(6,3) * t318 - t230;
t171 = -mrSges(6,1) * t207 + mrSges(6,3) * t196 - pkin(5) * t188 + (t232 + t233) * t318 + (Ifges(6,6) - Ifges(7,6)) * t290 + t371 * t262 + (Ifges(6,4) - Ifges(7,5)) * t219 + (-Ifges(6,2) - Ifges(7,3)) * t218 + t358;
t195 = t376 * t200 - t332 * t202;
t231 = Ifges(6,4) * t262 - Ifges(6,2) * t261 + Ifges(6,6) * t318;
t193 = -t290 * pkin(5) - t317 * qJ(6) + t262 * t236 + qJDD(6) - t195;
t228 = Ifges(7,5) * t262 + Ifges(7,6) * t318 + Ifges(7,3) * t261;
t355 = mrSges(7,2) * t193 - mrSges(7,3) * t198 + Ifges(7,1) * t219 + Ifges(7,4) * t290 + Ifges(7,5) * t218 + t318 * t228;
t172 = mrSges(6,2) * t207 - mrSges(6,3) * t195 + Ifges(6,1) * t219 - Ifges(6,4) * t218 + Ifges(6,5) * t290 - qJ(6) * t188 - t318 * t231 + t371 * t261 + t355;
t251 = Ifges(5,5) * t301 + Ifges(5,6) * t300 + Ifges(5,3) * t320;
t253 = Ifges(5,1) * t301 + Ifges(5,4) * t300 + Ifges(5,5) * t320;
t364 = m(7) * t191 + t290 * mrSges(7,3) + t318 * t249;
t237 = mrSges(7,1) * t261 - mrSges(7,3) * t262;
t370 = -mrSges(6,1) * t261 - mrSges(6,2) * t262 - t237;
t374 = -mrSges(6,3) - mrSges(7,2);
t179 = m(6) * t196 - t290 * mrSges(6,2) + t374 * t218 - t318 * t248 + t370 * t261 + t364;
t359 = -m(7) * t193 + t290 * mrSges(7,1) + t318 * t246;
t181 = m(6) * t195 + t290 * mrSges(6,1) + t374 * t219 + t318 * t247 + t370 * t262 + t359;
t360 = t376 * t179 - t181 * t332;
t154 = -mrSges(5,1) * t226 + mrSges(5,3) * t205 + Ifges(5,4) * t260 + Ifges(5,2) * t259 + Ifges(5,6) * t299 - pkin(4) * t347 + pkin(9) * t360 + t376 * t171 + t332 * t172 - t301 * t251 + t320 * t253;
t174 = t332 * t179 + t376 * t181;
t252 = Ifges(5,4) * t301 + Ifges(5,2) * t300 + Ifges(5,6) * t320;
t156 = mrSges(5,2) * t226 - mrSges(5,3) * t204 + Ifges(5,1) * t260 + Ifges(5,4) * t259 + Ifges(5,5) * t299 - pkin(9) * t174 - t332 * t171 + t376 * t172 + t300 * t251 - t320 * t252;
t263 = -mrSges(5,1) * t300 + mrSges(5,2) * t301;
t169 = m(5) * t204 + mrSges(5,1) * t299 - mrSges(5,3) * t260 - t263 * t301 + t264 * t320 + t174;
t170 = m(5) * t205 - mrSges(5,2) * t299 + mrSges(5,3) * t259 + t263 * t300 - t265 * t320 + t360;
t165 = t336 * t169 + t333 * t170;
t285 = Ifges(4,4) * qJD(2) + (-Ifges(4,2) * t334 - Ifges(4,6) * t337) * qJD(1);
t350 = -mrSges(4,2) * t245 + mrSges(4,3) * t243 - Ifges(4,1) * qJDD(2) + Ifges(4,4) * t305 + Ifges(4,5) * t306 + pkin(8) * t165 + t333 * t154 - t336 * t156 - t285 * t367;
t353 = -m(4) * t245 - t305 * mrSges(4,1) - t165;
t284 = Ifges(4,5) * qJD(2) + (-Ifges(4,6) * t334 - Ifges(4,3) * t337) * qJD(1);
t368 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t334 + Ifges(3,2) * t337) * qJD(1) - t284;
t378 = qJD(1) * (-t337 * t283 + t368 * t334) + mrSges(3,1) * t266 - mrSges(3,2) * t267 + Ifges(3,5) * t305 + Ifges(3,6) * t306 + Ifges(3,3) * qJDD(2) + pkin(2) * (-qJDD(2) * mrSges(4,2) - qJD(2) * t312 - t303 * t323 + t353) + qJ(3) * (t306 * mrSges(4,1) + t344) - t350;
t375 = t340 * pkin(7);
t373 = Ifges(3,4) + Ifges(4,6);
t166 = -t333 * t169 + t336 * t170;
t286 = Ifges(4,1) * qJD(2) + (-Ifges(4,4) * t334 - Ifges(4,5) * t337) * qJD(1);
t369 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t334 + Ifges(3,6) * t337) * qJD(1) + t286;
t304 = (-mrSges(3,1) * t337 + mrSges(3,2) * t334) * qJD(1);
t311 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t367;
t162 = m(3) * t266 - t305 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t311 - t312) * qJD(2) + (-t303 - t304) * t323 + t353;
t310 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t323;
t176 = -qJDD(2) * mrSges(3,2) + (mrSges(4,1) + mrSges(3,3)) * t306 - qJD(2) * t310 + m(3) * t267 + t344 + t304 * t367;
t361 = -t162 * t334 + t337 * t176;
t239 = -t306 * pkin(2) + t348 - t375;
t356 = -m(4) * t239 - t306 * mrSges(4,2) + t313 * t323 - t166;
t163 = -t305 * mrSges(4,3) + t312 * t367 - t356;
t287 = t357 - t375;
t349 = -mrSges(4,1) * t243 + mrSges(4,2) * t239 - pkin(3) * t183 - pkin(8) * t166 - t336 * t154 - t333 * t156;
t151 = -mrSges(3,1) * t287 + mrSges(3,3) * t267 - pkin(2) * t163 + (Ifges(3,2) + Ifges(4,3)) * t306 + t373 * t305 + (Ifges(3,6) - Ifges(4,5)) * qJDD(2) + (t283 - t285) * qJD(2) - t369 * t323 + t349;
t351 = mrSges(7,1) * t193 - mrSges(7,3) * t191 - Ifges(7,4) * t219 - Ifges(7,2) * t290 - Ifges(7,6) * t218 + t262 * t228 - t261 * t232;
t345 = mrSges(6,2) * t196 - t261 * t233 - qJ(6) * (-t218 * mrSges(7,2) - t261 * t237 + t364) - pkin(5) * (-t219 * mrSges(7,2) - t262 * t237 + t359) - mrSges(6,1) * t195 + Ifges(6,6) * t218 - Ifges(6,5) * t219 - t262 * t231 - Ifges(6,3) * t290 + t351;
t342 = -mrSges(5,1) * t204 + mrSges(5,2) * t205 - Ifges(5,5) * t260 - Ifges(5,6) * t259 - Ifges(5,3) * t299 - pkin(4) * t174 - t301 * t252 + t300 * t253 + t345;
t341 = -mrSges(4,1) * t245 + mrSges(4,3) * t239 - pkin(3) * t165 + t342;
t153 = -t341 + t373 * t306 + (Ifges(3,1) + Ifges(4,2)) * t305 + (Ifges(3,5) - Ifges(4,4)) * qJDD(2) - t368 * qJD(2) + mrSges(3,2) * t287 - mrSges(3,3) * t266 - qJ(3) * t163 + t369 * t367;
t346 = -m(3) * t287 + t311 * t367 + t306 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t305 + (-t310 * t334 - t312 * t337) * qJD(1) + t356;
t352 = mrSges(2,1) * t315 - mrSges(2,2) * t316 + Ifges(2,3) * qJDD(1) + pkin(1) * t346 + pkin(7) * t361 + t337 * t151 + t334 * t153;
t160 = m(2) * t315 + qJDD(1) * mrSges(2,1) - t340 * mrSges(2,2) + t346;
t159 = t162 * t337 + t176 * t334;
t157 = m(2) * t316 - mrSges(2,1) * t340 - qJDD(1) * mrSges(2,2) + t361;
t149 = mrSges(2,1) * g(3) + mrSges(2,3) * t316 + t340 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t159 - t378;
t148 = -mrSges(2,2) * g(3) - mrSges(2,3) * t315 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t340 - pkin(7) * t159 - t151 * t334 + t153 * t337;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t338 * t148 - t335 * t149 - pkin(6) * (t157 * t335 + t160 * t338), t148, t153, -t284 * t323 - t350, t156, t172, -t230 * t261 + t355; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t335 * t148 + t338 * t149 + pkin(6) * (t157 * t338 - t160 * t335), t149, t151, Ifges(4,4) * qJDD(2) - Ifges(4,2) * t305 - Ifges(4,6) * t306 - qJD(2) * t284 - t286 * t367 + t341, t154, t171, -t351; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t352, t352, t378, Ifges(4,5) * qJDD(2) - Ifges(4,6) * t305 - Ifges(4,3) * t306 + qJD(2) * t285 + t286 * t323 - t349, -t342, -t345, Ifges(7,5) * t219 + Ifges(7,6) * t290 + Ifges(7,3) * t218 + t262 * t230 - t318 * t232 - t358;];
m_new  = t1;
