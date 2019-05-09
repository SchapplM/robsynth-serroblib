% Calculate vector of cutting torques with Newton-Euler for
% S6RRPPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
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
% Datum: 2019-05-06 09:29
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPPRP5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:25:12
% EndTime: 2019-05-06 09:25:27
% DurationCPUTime: 7.20s
% Computational Cost: add. (112550->386), mult. (245082->454), div. (0->0), fcn. (147376->8), ass. (0->135)
t379 = -2 * qJD(3);
t335 = sin(qJ(1));
t337 = cos(qJ(1));
t317 = -g(1) * t337 - g(2) * t335;
t339 = qJD(1) ^ 2;
t285 = -pkin(1) * t339 + qJDD(1) * pkin(7) + t317;
t334 = sin(qJ(2));
t336 = cos(qJ(2));
t259 = -t336 * g(3) - t334 * t285;
t260 = -t334 * g(3) + t336 * t285;
t280 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t334 + Ifges(3,4) * t336) * qJD(1);
t304 = (mrSges(4,2) * t336 - mrSges(4,3) * t334) * qJD(1);
t365 = qJD(1) * qJD(2);
t362 = t336 * t365;
t306 = qJDD(1) * t334 + t362;
t361 = t334 * t365;
t307 = qJDD(1) * t336 - t361;
t366 = qJD(1) * t336;
t314 = -mrSges(4,1) * t366 - qJD(2) * mrSges(4,3);
t303 = (-pkin(2) * t336 - qJ(3) * t334) * qJD(1);
t338 = qJD(2) ^ 2;
t243 = t338 * pkin(2) - qJDD(2) * qJ(3) + qJD(2) * t379 - t303 * t366 - t260;
t367 = qJD(1) * t334;
t313 = pkin(3) * t367 - qJD(2) * qJ(4);
t330 = t336 ^ 2;
t226 = -t330 * t339 * qJ(4) + t307 * pkin(3) + qJD(2) * t313 + qJDD(4) - t243;
t331 = sin(pkin(9));
t332 = cos(pkin(9));
t294 = -qJD(2) * t331 - t332 * t366;
t264 = -mrSges(5,2) * t367 + mrSges(5,3) * t294;
t295 = qJD(2) * t332 - t331 * t366;
t265 = mrSges(5,1) * t367 - mrSges(5,3) * t295;
t266 = -qJDD(2) * t331 - t307 * t332;
t267 = qJDD(2) * t332 - t307 * t331;
t268 = pkin(4) * t367 - pkin(8) * t295;
t293 = t294 ^ 2;
t207 = -t266 * pkin(4) - t293 * pkin(8) + t295 * t268 + t226;
t333 = sin(qJ(5));
t376 = cos(qJ(5));
t256 = t333 * t294 + t376 * t295;
t221 = qJD(5) * t256 - t376 * t266 + t267 * t333;
t255 = -t376 * t294 + t295 * t333;
t222 = -t255 * qJD(5) + t333 * t266 + t376 * t267;
t320 = qJD(5) + t367;
t198 = -0.2e1 * qJD(6) * t256 + t207 + (t255 * t320 - t222) * qJ(6) + (t256 * t320 + t221) * pkin(5);
t247 = -mrSges(7,2) * t255 + mrSges(7,3) * t320;
t249 = -mrSges(7,1) * t320 + mrSges(7,2) * t256;
t188 = m(7) * t198 + t221 * mrSges(7,1) - t222 * mrSges(7,3) + t255 * t247 - t256 * t249;
t246 = -mrSges(6,2) * t320 - mrSges(6,3) * t255;
t248 = mrSges(6,1) * t320 - mrSges(6,3) * t256;
t346 = m(6) * t207 + t221 * mrSges(6,1) + t222 * mrSges(6,2) + t255 * t246 + t256 * t248 + t188;
t183 = -m(5) * t226 + t266 * mrSges(5,1) - t267 * mrSges(5,2) + t294 * t264 - t295 * t265 - t346;
t315 = mrSges(4,1) * t367 + qJD(2) * mrSges(4,2);
t343 = -m(4) * t243 + qJDD(2) * mrSges(4,3) + qJD(2) * t315 + t304 * t366 - t183;
t316 = t335 * g(1) - t337 * g(2);
t356 = -qJDD(1) * pkin(1) - t316;
t348 = pkin(2) * t361 + t367 * t379 + (-t306 - t362) * qJ(3) + t356;
t211 = -t313 * t367 + (-pkin(3) * t330 - pkin(7)) * t339 + (-pkin(2) - qJ(4)) * t307 + t348;
t245 = -qJDD(2) * pkin(2) - t338 * qJ(3) + t303 * t367 + qJDD(3) - t259;
t233 = (-t334 * t336 * t339 - qJDD(2)) * qJ(4) + (t306 - t362) * pkin(3) + t245;
t204 = -0.2e1 * qJD(4) * t295 - t331 * t211 + t332 * t233;
t200 = (t294 * t367 - t267) * pkin(8) + (t294 * t295 + t306) * pkin(4) + t204;
t205 = 0.2e1 * qJD(4) * t294 + t332 * t211 + t331 * t233;
t202 = -pkin(4) * t293 + pkin(8) * t266 - t268 * t367 + t205;
t196 = t333 * t200 + t376 * t202;
t231 = Ifges(7,1) * t256 + Ifges(7,4) * t320 + Ifges(7,5) * t255;
t232 = Ifges(6,1) * t256 - Ifges(6,4) * t255 + Ifges(6,5) * t320;
t302 = qJDD(5) + t306;
t236 = pkin(5) * t255 - qJ(6) * t256;
t318 = t320 ^ 2;
t191 = -pkin(5) * t318 + qJ(6) * t302 + 0.2e1 * qJD(6) * t320 - t236 * t255 + t196;
t357 = -mrSges(7,1) * t198 + mrSges(7,2) * t191;
t229 = Ifges(7,4) * t256 + Ifges(7,2) * t320 + Ifges(7,6) * t255;
t371 = -Ifges(6,5) * t256 + Ifges(6,6) * t255 - Ifges(6,3) * t320 - t229;
t173 = -mrSges(6,1) * t207 + mrSges(6,3) * t196 - pkin(5) * t188 + (t231 + t232) * t320 + (Ifges(6,6) - Ifges(7,6)) * t302 + t371 * t256 + (Ifges(6,4) - Ifges(7,5)) * t222 + (-Ifges(6,2) - Ifges(7,3)) * t221 + t357;
t195 = t376 * t200 - t333 * t202;
t230 = Ifges(6,4) * t256 - Ifges(6,2) * t255 + Ifges(6,6) * t320;
t193 = -t302 * pkin(5) - t318 * qJ(6) + t256 * t236 + qJDD(6) - t195;
t227 = Ifges(7,5) * t256 + Ifges(7,6) * t320 + Ifges(7,3) * t255;
t354 = mrSges(7,2) * t193 - mrSges(7,3) * t198 + Ifges(7,1) * t222 + Ifges(7,4) * t302 + Ifges(7,5) * t221 + t320 * t227;
t174 = mrSges(6,2) * t207 - mrSges(6,3) * t195 + Ifges(6,1) * t222 - Ifges(6,4) * t221 + Ifges(6,5) * t302 - qJ(6) * t188 - t320 * t230 + t371 * t255 + t354;
t250 = Ifges(5,5) * t295 + Ifges(5,6) * t294 + Ifges(5,3) * t367;
t252 = Ifges(5,1) * t295 + Ifges(5,4) * t294 + Ifges(5,5) * t367;
t363 = m(7) * t191 + t302 * mrSges(7,3) + t320 * t249;
t237 = mrSges(7,1) * t255 - mrSges(7,3) * t256;
t370 = -mrSges(6,1) * t255 - mrSges(6,2) * t256 - t237;
t374 = -mrSges(6,3) - mrSges(7,2);
t177 = m(6) * t196 - t302 * mrSges(6,2) + t374 * t221 - t320 * t248 + t370 * t255 + t363;
t358 = -m(7) * t193 + t302 * mrSges(7,1) + t320 * t247;
t179 = m(6) * t195 + t302 * mrSges(6,1) + t374 * t222 + t320 * t246 + t370 * t256 + t358;
t359 = t376 * t177 - t179 * t333;
t154 = -mrSges(5,1) * t226 + mrSges(5,3) * t205 + Ifges(5,4) * t267 + Ifges(5,2) * t266 + Ifges(5,6) * t306 - pkin(4) * t346 + pkin(8) * t359 + t376 * t173 + t333 * t174 - t295 * t250 + t252 * t367;
t172 = t333 * t177 + t376 * t179;
t251 = Ifges(5,4) * t295 + Ifges(5,2) * t294 + Ifges(5,6) * t367;
t156 = mrSges(5,2) * t226 - mrSges(5,3) * t204 + Ifges(5,1) * t267 + Ifges(5,4) * t266 + Ifges(5,5) * t306 - pkin(8) * t172 - t333 * t173 + t376 * t174 + t294 * t250 - t251 * t367;
t257 = -mrSges(5,1) * t294 + mrSges(5,2) * t295;
t169 = m(5) * t204 + mrSges(5,1) * t306 - mrSges(5,3) * t267 - t257 * t295 + t264 * t367 + t172;
t170 = m(5) * t205 - mrSges(5,2) * t306 + mrSges(5,3) * t266 + t257 * t294 - t265 * t367 + t359;
t165 = t332 * t169 + t331 * t170;
t282 = Ifges(4,4) * qJD(2) + (-Ifges(4,2) * t334 - Ifges(4,6) * t336) * qJD(1);
t349 = -mrSges(4,2) * t245 + mrSges(4,3) * t243 - Ifges(4,1) * qJDD(2) + Ifges(4,4) * t306 + Ifges(4,5) * t307 + qJ(4) * t165 + t331 * t154 - t332 * t156 - t282 * t366;
t352 = -m(4) * t245 - t306 * mrSges(4,1) - t165;
t281 = Ifges(4,5) * qJD(2) + (-Ifges(4,6) * t334 - Ifges(4,3) * t336) * qJD(1);
t368 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t334 + Ifges(3,2) * t336) * qJD(1) - t281;
t378 = (-t336 * t280 + t368 * t334) * qJD(1) + mrSges(3,1) * t259 - mrSges(3,2) * t260 + Ifges(3,5) * t306 + Ifges(3,6) * t307 + Ifges(3,3) * qJDD(2) + pkin(2) * (-qJDD(2) * mrSges(4,2) - qJD(2) * t314 - t304 * t367 + t352) + qJ(3) * (t307 * mrSges(4,1) + t343) - t349;
t375 = t339 * pkin(7);
t373 = Ifges(3,4) + Ifges(4,6);
t166 = -t331 * t169 + t332 * t170;
t283 = Ifges(4,1) * qJD(2) + (-Ifges(4,4) * t334 - Ifges(4,5) * t336) * qJD(1);
t369 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t334 + Ifges(3,6) * t336) * qJD(1) + t283;
t305 = (-mrSges(3,1) * t336 + mrSges(3,2) * t334) * qJD(1);
t312 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t366;
t162 = m(3) * t259 - t306 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t312 - t314) * qJD(2) + (-t304 - t305) * t367 + t352;
t311 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t367;
t181 = (mrSges(3,3) + mrSges(4,1)) * t307 - qJDD(2) * mrSges(3,2) + t343 - qJD(2) * t311 + m(3) * t260 + t305 * t366;
t360 = -t162 * t334 + t336 * t181;
t239 = -t307 * pkin(2) + t348 - t375;
t355 = -m(4) * t239 - t307 * mrSges(4,2) + t315 * t367 - t166;
t163 = -t306 * mrSges(4,3) + t314 * t366 - t355;
t284 = t356 - t375;
t347 = -mrSges(4,1) * t243 + mrSges(4,2) * t239 - pkin(3) * t183 - qJ(4) * t166 - t332 * t154 - t331 * t156;
t151 = -mrSges(3,1) * t284 + mrSges(3,3) * t260 - pkin(2) * t163 + (Ifges(3,2) + Ifges(4,3)) * t307 + t373 * t306 + (Ifges(3,6) - Ifges(4,5)) * qJDD(2) + (t280 - t282) * qJD(2) - t369 * t367 + t347;
t350 = mrSges(7,1) * t193 - mrSges(7,3) * t191 - Ifges(7,4) * t222 - Ifges(7,2) * t302 - Ifges(7,6) * t221 + t256 * t227 - t255 * t231;
t344 = mrSges(6,2) * t196 - t255 * t232 - qJ(6) * (-t221 * mrSges(7,2) - t255 * t237 + t363) - pkin(5) * (-t222 * mrSges(7,2) - t256 * t237 + t358) - mrSges(6,1) * t195 - t256 * t230 + Ifges(6,6) * t221 - Ifges(6,5) * t222 - Ifges(6,3) * t302 + t350;
t341 = -mrSges(5,1) * t204 + mrSges(5,2) * t205 - Ifges(5,5) * t267 - Ifges(5,6) * t266 - Ifges(5,3) * t306 - pkin(4) * t172 - t295 * t251 + t294 * t252 + t344;
t340 = -mrSges(4,1) * t245 + mrSges(4,3) * t239 - pkin(3) * t165 + t341;
t153 = -t340 + t373 * t307 + (Ifges(3,1) + Ifges(4,2)) * t306 + (Ifges(3,5) - Ifges(4,4)) * qJDD(2) - t368 * qJD(2) + mrSges(3,2) * t284 - mrSges(3,3) * t259 - qJ(3) * t163 + t369 * t366;
t345 = -m(3) * t284 + t312 * t366 + t307 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t306 + (-t311 * t334 - t314 * t336) * qJD(1) + t355;
t351 = mrSges(2,1) * t316 - mrSges(2,2) * t317 + Ifges(2,3) * qJDD(1) + pkin(1) * t345 + pkin(7) * t360 + t336 * t151 + t334 * t153;
t160 = m(2) * t316 + qJDD(1) * mrSges(2,1) - t339 * mrSges(2,2) + t345;
t159 = t162 * t336 + t181 * t334;
t157 = m(2) * t317 - mrSges(2,1) * t339 - qJDD(1) * mrSges(2,2) + t360;
t149 = mrSges(2,1) * g(3) + mrSges(2,3) * t317 + t339 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t159 - t378;
t148 = -mrSges(2,2) * g(3) - mrSges(2,3) * t316 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t339 - pkin(7) * t159 - t151 * t334 + t153 * t336;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t337 * t148 - t335 * t149 - pkin(6) * (t157 * t335 + t160 * t337), t148, t153, -t281 * t367 - t349, t156, t174, -t229 * t255 + t354; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t335 * t148 + t337 * t149 + pkin(6) * (t157 * t337 - t160 * t335), t149, t151, Ifges(4,4) * qJDD(2) - Ifges(4,2) * t306 - Ifges(4,6) * t307 - qJD(2) * t281 - t283 * t366 + t340, t154, t173, -t350; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t351, t351, t378, Ifges(4,5) * qJDD(2) - Ifges(4,6) * t306 - Ifges(4,3) * t307 + qJD(2) * t282 + t283 * t367 - t347, -t341, -t344, Ifges(7,5) * t222 + Ifges(7,6) * t302 + Ifges(7,3) * t221 + t256 * t229 - t320 * t231 - t357;];
m_new  = t1;
