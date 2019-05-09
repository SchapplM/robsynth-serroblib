% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
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
% Datum: 2019-05-06 12:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRPP4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:35:34
% EndTime: 2019-05-06 12:35:51
% DurationCPUTime: 7.51s
% Computational Cost: add. (119351->387), mult. (250358->454), div. (0->0), fcn. (151156->8), ass. (0->137)
t380 = -2 * qJD(3);
t333 = sin(qJ(1));
t336 = cos(qJ(1));
t316 = -g(1) * t336 - g(2) * t333;
t338 = qJD(1) ^ 2;
t284 = -pkin(1) * t338 + qJDD(1) * pkin(7) + t316;
t332 = sin(qJ(2));
t335 = cos(qJ(2));
t268 = -t335 * g(3) - t332 * t284;
t269 = -t332 * g(3) + t335 * t284;
t279 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t332 + Ifges(3,4) * t335) * qJD(1);
t303 = (mrSges(4,2) * t335 - mrSges(4,3) * t332) * qJD(1);
t365 = qJD(1) * qJD(2);
t362 = t335 * t365;
t305 = qJDD(1) * t332 + t362;
t361 = t332 * t365;
t306 = qJDD(1) * t335 - t361;
t366 = qJD(1) * t335;
t312 = -mrSges(4,1) * t366 - qJD(2) * mrSges(4,3);
t302 = (-pkin(2) * t335 - qJ(3) * t332) * qJD(1);
t337 = qJD(2) ^ 2;
t243 = t337 * pkin(2) - qJDD(2) * qJ(3) + qJD(2) * t380 - t302 * t366 - t269;
t367 = qJD(1) * t332;
t314 = pkin(3) * t367 - qJD(2) * pkin(8);
t329 = t335 ^ 2;
t224 = -t329 * t338 * pkin(8) + t306 * pkin(3) + qJD(2) * t314 - t243;
t331 = sin(qJ(4));
t334 = cos(qJ(4));
t301 = qJD(2) * t334 - t331 * t366;
t260 = -qJD(4) * t301 - qJDD(2) * t331 - t306 * t334;
t300 = -qJD(2) * t331 - t334 * t366;
t261 = qJD(4) * t300 + qJDD(2) * t334 - t306 * t331;
t319 = qJD(4) + t367;
t265 = -mrSges(5,2) * t319 + mrSges(5,3) * t300;
t267 = mrSges(5,1) * t319 - mrSges(5,3) * t301;
t266 = pkin(4) * t319 - qJ(5) * t301;
t298 = t300 ^ 2;
t207 = -t260 * pkin(4) - t298 * qJ(5) + t301 * t266 + qJDD(5) + t224;
t330 = sin(pkin(9));
t373 = cos(pkin(9));
t232 = -t373 * t260 + t330 * t261;
t233 = t330 * t260 + t373 * t261;
t262 = -t373 * t300 + t330 * t301;
t263 = t330 * t300 + t373 * t301;
t198 = (t262 * t319 - t233) * qJ(6) + (t263 * t319 + t232) * pkin(5) + t207 - 0.2e1 * qJD(6) * t263;
t247 = -mrSges(7,2) * t262 + mrSges(7,3) * t319;
t249 = -mrSges(7,1) * t319 + mrSges(7,2) * t263;
t188 = m(7) * t198 + t232 * mrSges(7,1) - t233 * mrSges(7,3) + t262 * t247 - t263 * t249;
t246 = -mrSges(6,2) * t319 - mrSges(6,3) * t262;
t248 = mrSges(6,1) * t319 - mrSges(6,3) * t263;
t345 = m(6) * t207 + t232 * mrSges(6,1) + t233 * mrSges(6,2) + t262 * t246 + t263 * t248 + t188;
t183 = -m(5) * t224 + t260 * mrSges(5,1) - t261 * mrSges(5,2) + t300 * t265 - t301 * t267 - t345;
t313 = mrSges(4,1) * t367 + qJD(2) * mrSges(4,2);
t342 = -m(4) * t243 + qJDD(2) * mrSges(4,3) + qJD(2) * t313 + t303 * t366 - t183;
t315 = t333 * g(1) - t336 * g(2);
t356 = -qJDD(1) * pkin(1) - t315;
t346 = pkin(2) * t361 + t367 * t380 + (-t305 - t362) * qJ(3) + t356;
t211 = -t314 * t367 + (-pkin(3) * t329 - pkin(7)) * t338 + (-pkin(2) - pkin(8)) * t306 + t346;
t245 = -qJDD(2) * pkin(2) - t337 * qJ(3) + t302 * t367 + qJDD(3) - t268;
t225 = (-t332 * t335 * t338 - qJDD(2)) * pkin(8) + (t305 - t362) * pkin(3) + t245;
t204 = -t331 * t211 + t334 * t225;
t299 = qJDD(4) + t305;
t200 = (t300 * t319 - t261) * qJ(5) + (t300 * t301 + t299) * pkin(4) + t204;
t205 = t334 * t211 + t331 * t225;
t202 = -pkin(4) * t298 + qJ(5) * t260 - t266 * t319 + t205;
t377 = -2 * qJD(5);
t196 = t330 * t200 + t373 * t202 + t262 * t377;
t230 = Ifges(7,1) * t263 + Ifges(7,4) * t319 + Ifges(7,5) * t262;
t231 = Ifges(6,1) * t263 - Ifges(6,4) * t262 + Ifges(6,5) * t319;
t236 = pkin(5) * t262 - qJ(6) * t263;
t317 = t319 ^ 2;
t191 = -pkin(5) * t317 + qJ(6) * t299 + 0.2e1 * qJD(6) * t319 - t236 * t262 + t196;
t357 = -mrSges(7,1) * t198 + mrSges(7,2) * t191;
t228 = Ifges(7,4) * t263 + Ifges(7,2) * t319 + Ifges(7,6) * t262;
t371 = -Ifges(6,5) * t263 + Ifges(6,6) * t262 - Ifges(6,3) * t319 - t228;
t173 = -mrSges(6,1) * t207 + mrSges(6,3) * t196 - pkin(5) * t188 + (t230 + t231) * t319 + (Ifges(6,6) - Ifges(7,6)) * t299 + t371 * t263 + (Ifges(6,4) - Ifges(7,5)) * t233 + (-Ifges(6,2) - Ifges(7,3)) * t232 + t357;
t354 = t373 * t200 - t330 * t202;
t195 = t263 * t377 + t354;
t229 = Ifges(6,4) * t263 - Ifges(6,2) * t262 + Ifges(6,6) * t319;
t193 = -t299 * pkin(5) - t317 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t236) * t263 - t354;
t226 = Ifges(7,5) * t263 + Ifges(7,6) * t319 + Ifges(7,3) * t262;
t353 = mrSges(7,2) * t193 - mrSges(7,3) * t198 + Ifges(7,1) * t233 + Ifges(7,4) * t299 + Ifges(7,5) * t232 + t319 * t226;
t174 = mrSges(6,2) * t207 - mrSges(6,3) * t195 + Ifges(6,1) * t233 - Ifges(6,4) * t232 + Ifges(6,5) * t299 - qJ(6) * t188 - t319 * t229 + t371 * t262 + t353;
t251 = Ifges(5,5) * t301 + Ifges(5,6) * t300 + Ifges(5,3) * t319;
t253 = Ifges(5,1) * t301 + Ifges(5,4) * t300 + Ifges(5,5) * t319;
t363 = m(7) * t191 + t299 * mrSges(7,3) + t319 * t249;
t237 = mrSges(7,1) * t262 - mrSges(7,3) * t263;
t370 = -mrSges(6,1) * t262 - mrSges(6,2) * t263 - t237;
t375 = -mrSges(6,3) - mrSges(7,2);
t177 = m(6) * t196 - t299 * mrSges(6,2) + t375 * t232 - t319 * t248 + t370 * t262 + t363;
t358 = -m(7) * t193 + t299 * mrSges(7,1) + t319 * t247;
t179 = m(6) * t195 + t299 * mrSges(6,1) + t375 * t233 + t319 * t246 + t370 * t263 + t358;
t359 = t373 * t177 - t179 * t330;
t154 = -mrSges(5,1) * t224 + mrSges(5,3) * t205 + Ifges(5,4) * t261 + Ifges(5,2) * t260 + Ifges(5,6) * t299 - pkin(4) * t345 + qJ(5) * t359 + t373 * t173 + t330 * t174 - t301 * t251 + t319 * t253;
t172 = t330 * t177 + t373 * t179;
t252 = Ifges(5,4) * t301 + Ifges(5,2) * t300 + Ifges(5,6) * t319;
t156 = mrSges(5,2) * t224 - mrSges(5,3) * t204 + Ifges(5,1) * t261 + Ifges(5,4) * t260 + Ifges(5,5) * t299 - qJ(5) * t172 - t330 * t173 + t373 * t174 + t300 * t251 - t319 * t252;
t264 = -mrSges(5,1) * t300 + mrSges(5,2) * t301;
t169 = m(5) * t204 + mrSges(5,1) * t299 - mrSges(5,3) * t261 - t264 * t301 + t265 * t319 + t172;
t170 = m(5) * t205 - mrSges(5,2) * t299 + mrSges(5,3) * t260 + t264 * t300 - t267 * t319 + t359;
t165 = t334 * t169 + t331 * t170;
t281 = Ifges(4,4) * qJD(2) + (-Ifges(4,2) * t332 - Ifges(4,6) * t335) * qJD(1);
t348 = -mrSges(4,2) * t245 + mrSges(4,3) * t243 - Ifges(4,1) * qJDD(2) + Ifges(4,4) * t305 + Ifges(4,5) * t306 + pkin(8) * t165 + t331 * t154 - t334 * t156 - t281 * t366;
t351 = -m(4) * t245 - t305 * mrSges(4,1) - t165;
t280 = Ifges(4,5) * qJD(2) + (-Ifges(4,6) * t332 - Ifges(4,3) * t335) * qJD(1);
t368 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t332 + Ifges(3,2) * t335) * qJD(1) - t280;
t379 = qJD(1) * (-t335 * t279 + t368 * t332) + mrSges(3,1) * t268 - mrSges(3,2) * t269 + Ifges(3,5) * t305 + Ifges(3,6) * t306 + Ifges(3,3) * qJDD(2) + pkin(2) * (-qJDD(2) * mrSges(4,2) - qJD(2) * t312 - t303 * t367 + t351) + qJ(3) * (t306 * mrSges(4,1) + t342) - t348;
t376 = t338 * pkin(7);
t374 = Ifges(3,4) + Ifges(4,6);
t166 = -t331 * t169 + t334 * t170;
t282 = Ifges(4,1) * qJD(2) + (-Ifges(4,4) * t332 - Ifges(4,5) * t335) * qJD(1);
t369 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t332 + Ifges(3,6) * t335) * qJD(1) + t282;
t304 = (-mrSges(3,1) * t335 + mrSges(3,2) * t332) * qJD(1);
t311 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t366;
t162 = m(3) * t268 - t305 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t311 - t312) * qJD(2) + (-t303 - t304) * t367 + t351;
t310 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t367;
t181 = -qJDD(2) * mrSges(3,2) + t342 + t304 * t366 + (mrSges(3,3) + mrSges(4,1)) * t306 + m(3) * t269 - qJD(2) * t310;
t360 = -t162 * t332 + t335 * t181;
t239 = -t306 * pkin(2) + t346 - t376;
t355 = -m(4) * t239 - t306 * mrSges(4,2) + t313 * t367 - t166;
t163 = -t305 * mrSges(4,3) + t312 * t366 - t355;
t283 = t356 - t376;
t347 = -mrSges(4,1) * t243 + mrSges(4,2) * t239 - pkin(3) * t183 - pkin(8) * t166 - t334 * t154 - t331 * t156;
t151 = -mrSges(3,1) * t283 + mrSges(3,3) * t269 - pkin(2) * t163 + (Ifges(3,2) + Ifges(4,3)) * t306 + t374 * t305 + (Ifges(3,6) - Ifges(4,5)) * qJDD(2) + (t279 - t281) * qJD(2) - t369 * t367 + t347;
t349 = mrSges(7,1) * t193 - mrSges(7,3) * t191 - Ifges(7,4) * t233 - Ifges(7,2) * t299 - Ifges(7,6) * t232 + t263 * t226 - t262 * t230;
t343 = mrSges(6,2) * t196 - t262 * t231 - qJ(6) * (-t232 * mrSges(7,2) - t262 * t237 + t363) - pkin(5) * (-t233 * mrSges(7,2) - t263 * t237 + t358) - mrSges(6,1) * t195 - t263 * t229 + Ifges(6,6) * t232 - Ifges(6,5) * t233 - Ifges(6,3) * t299 + t349;
t340 = -mrSges(5,1) * t204 + mrSges(5,2) * t205 - Ifges(5,5) * t261 - Ifges(5,6) * t260 - Ifges(5,3) * t299 - pkin(4) * t172 - t301 * t252 + t300 * t253 + t343;
t339 = -mrSges(4,1) * t245 + mrSges(4,3) * t239 - pkin(3) * t165 + t340;
t153 = -qJ(3) * t163 - t339 + t369 * t366 + t374 * t306 + (Ifges(3,1) + Ifges(4,2)) * t305 + (Ifges(3,5) - Ifges(4,4)) * qJDD(2) - t368 * qJD(2) - mrSges(3,3) * t268 + mrSges(3,2) * t283;
t344 = -m(3) * t283 + t311 * t366 + t306 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t305 + (-t310 * t332 - t312 * t335) * qJD(1) + t355;
t350 = mrSges(2,1) * t315 - mrSges(2,2) * t316 + Ifges(2,3) * qJDD(1) + pkin(1) * t344 + pkin(7) * t360 + t335 * t151 + t332 * t153;
t160 = m(2) * t315 + qJDD(1) * mrSges(2,1) - t338 * mrSges(2,2) + t344;
t159 = t162 * t335 + t181 * t332;
t157 = m(2) * t316 - mrSges(2,1) * t338 - qJDD(1) * mrSges(2,2) + t360;
t149 = mrSges(2,1) * g(3) + mrSges(2,3) * t316 + t338 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t159 - t379;
t148 = -mrSges(2,2) * g(3) - mrSges(2,3) * t315 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t338 - pkin(7) * t159 - t151 * t332 + t153 * t335;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t336 * t148 - t333 * t149 - pkin(6) * (t157 * t333 + t160 * t336), t148, t153, -t280 * t367 - t348, t156, t174, -t228 * t262 + t353; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t333 * t148 + t336 * t149 + pkin(6) * (t157 * t336 - t160 * t333), t149, t151, Ifges(4,4) * qJDD(2) - Ifges(4,2) * t305 - Ifges(4,6) * t306 - qJD(2) * t280 - t282 * t366 + t339, t154, t173, -t349; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t350, t350, t379, Ifges(4,5) * qJDD(2) - Ifges(4,6) * t305 - Ifges(4,3) * t306 + qJD(2) * t281 + t282 * t367 - t347, -t340, -t343, Ifges(7,5) * t233 + Ifges(7,6) * t299 + Ifges(7,3) * t232 + t263 * t228 - t319 * t230 - t357;];
m_new  = t1;
