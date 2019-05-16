% Calculate vector of cutting torques with Newton-Euler for
% S6PRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2019-05-05 04:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRPRP4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:00:40
% EndTime: 2019-05-05 04:00:55
% DurationCPUTime: 7.12s
% Computational Cost: add. (89578->344), mult. (175687->407), div. (0->0), fcn. (104362->10), ass. (0->134)
t315 = sin(qJ(5));
t318 = cos(qJ(5));
t319 = cos(qJ(3));
t351 = qJD(2) * t319;
t283 = -qJD(3) * t315 - t318 * t351;
t316 = sin(qJ(3));
t350 = qJD(2) * qJD(3);
t344 = t316 * t350;
t289 = qJDD(2) * t319 - t344;
t245 = qJD(5) * t283 + qJDD(3) * t318 - t289 * t315;
t284 = qJD(3) * t318 - t315 * t351;
t247 = -mrSges(7,1) * t283 + mrSges(7,2) * t284;
t311 = sin(pkin(10));
t313 = cos(pkin(10));
t293 = g(1) * t311 - g(2) * t313;
t310 = -g(3) + qJDD(1);
t312 = sin(pkin(6));
t314 = cos(pkin(6));
t255 = -t293 * t312 + t310 * t314;
t345 = t319 * t350;
t288 = qJDD(2) * t316 + t345;
t294 = -g(1) * t313 - g(2) * t311;
t320 = cos(qJ(2));
t317 = sin(qJ(2));
t358 = t314 * t317;
t359 = t312 * t317;
t224 = t293 * t358 + t320 * t294 + t310 * t359;
t322 = qJD(2) ^ 2;
t221 = -pkin(2) * t322 + qJDD(2) * pkin(8) + t224;
t217 = t316 * t221;
t285 = (-pkin(3) * t319 - qJ(4) * t316) * qJD(2);
t321 = qJD(3) ^ 2;
t352 = qJD(2) * t316;
t340 = -qJ(4) * t321 + t285 * t352 + qJDD(4) + t217;
t366 = pkin(9) * t322;
t368 = -pkin(3) - pkin(9);
t205 = pkin(4) * t288 + t368 * qJDD(3) + (-pkin(4) * t350 - t316 * t366 - t255) * t319 + t340;
t300 = pkin(4) * t352 - qJD(3) * pkin(9);
t309 = t319 ^ 2;
t223 = -t317 * t294 + (t293 * t314 + t310 * t312) * t320;
t331 = -qJDD(2) * pkin(2) - t223;
t374 = -2 * qJD(4);
t326 = pkin(3) * t344 + t352 * t374 + (-t288 - t345) * qJ(4) + t331;
t207 = -t300 * t352 + (-pkin(4) * t309 - pkin(8)) * t322 + t368 * t289 + t326;
t196 = t318 * t205 - t315 * t207;
t280 = qJDD(5) + t288;
t302 = qJD(5) + t352;
t191 = -0.2e1 * qJD(6) * t284 + (t283 * t302 - t245) * qJ(6) + (t283 * t284 + t280) * pkin(5) + t196;
t250 = -mrSges(7,2) * t302 + mrSges(7,3) * t283;
t348 = m(7) * t191 + t280 * mrSges(7,1) + t302 * t250;
t188 = -mrSges(7,3) * t245 - t247 * t284 + t348;
t197 = t315 * t205 + t318 * t207;
t228 = Ifges(6,4) * t284 + Ifges(6,2) * t283 + Ifges(6,6) * t302;
t229 = Ifges(7,1) * t284 + Ifges(7,4) * t283 + Ifges(7,5) * t302;
t230 = Ifges(6,1) * t284 + Ifges(6,4) * t283 + Ifges(6,5) * t302;
t244 = -qJD(5) * t284 - qJDD(3) * t315 - t289 * t318;
t252 = pkin(5) * t302 - qJ(6) * t284;
t279 = t283 ^ 2;
t194 = -pkin(5) * t279 + qJ(6) * t244 + 0.2e1 * qJD(6) * t283 - t252 * t302 + t197;
t227 = Ifges(7,4) * t284 + Ifges(7,2) * t283 + Ifges(7,6) * t302;
t335 = -mrSges(7,1) * t191 + mrSges(7,2) * t194 - Ifges(7,5) * t245 - Ifges(7,6) * t244 - Ifges(7,3) * t280 - t284 * t227;
t375 = mrSges(6,1) * t196 - mrSges(6,2) * t197 + Ifges(6,5) * t245 + Ifges(6,6) * t244 + Ifges(6,3) * t280 + pkin(5) * t188 + t284 * t228 - t335 - (t230 + t229) * t283;
t248 = -mrSges(6,1) * t283 + mrSges(6,2) * t284;
t251 = -mrSges(6,2) * t302 + mrSges(6,3) * t283;
t182 = m(6) * t196 + mrSges(6,1) * t280 + t251 * t302 + (-t247 - t248) * t284 + (-mrSges(6,3) - mrSges(7,3)) * t245 + t348;
t253 = mrSges(7,1) * t302 - mrSges(7,3) * t284;
t254 = mrSges(6,1) * t302 - mrSges(6,3) * t284;
t347 = m(7) * t194 + t244 * mrSges(7,3) + t283 * t247;
t184 = m(6) * t197 + mrSges(6,3) * t244 + t248 * t283 + (-t253 - t254) * t302 + (-mrSges(6,2) - mrSges(7,2)) * t280 + t347;
t176 = t182 * t318 + t184 * t315;
t360 = t255 * t319;
t210 = -qJDD(3) * pkin(3) + t340 - t360;
t365 = t322 * pkin(8);
t211 = -pkin(3) * t289 + t326 - t365;
t373 = mrSges(5,1) * t210 - mrSges(5,3) * t211 + pkin(4) * t176 + t375;
t212 = -t217 + t360;
t286 = (mrSges(5,2) * t319 - mrSges(5,3) * t316) * qJD(2);
t287 = (-mrSges(4,1) * t319 + mrSges(4,2) * t316) * qJD(2);
t296 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t351;
t297 = -mrSges(5,1) * t351 - qJD(3) * mrSges(5,3);
t333 = -m(5) * t210 - t288 * mrSges(5,1) - t176;
t171 = m(4) * t212 - mrSges(4,3) * t288 + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t296 - t297) * qJD(3) + (-t286 - t287) * t352 + t333;
t213 = t319 * t221 + t316 * t255;
t295 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t352;
t208 = pkin(3) * t321 - qJDD(3) * qJ(4) + qJD(3) * t374 - t285 * t351 - t213;
t298 = mrSges(5,1) * t352 + qJD(3) * mrSges(5,2);
t204 = pkin(4) * t289 + qJD(3) * t300 - t309 * t366 - t208;
t200 = -pkin(5) * t244 - qJ(6) * t279 + t252 * t284 + qJDD(6) + t204;
t346 = m(7) * t200 + t245 * mrSges(7,2) + t284 * t253;
t369 = -m(6) * t204 - t245 * mrSges(6,2) + (mrSges(6,1) + mrSges(7,1)) * t244 - t284 * t254 + (t250 + t251) * t283 - t346;
t327 = -m(5) * t208 + qJDD(3) * mrSges(5,3) + qJD(3) * t298 + t286 * t351 - t369;
t180 = (mrSges(4,3) + mrSges(5,1)) * t289 - qJD(3) * t295 + t287 * t351 + m(4) * t213 - qJDD(3) * mrSges(4,2) + t327;
t166 = t319 * t171 + t316 * t180;
t263 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t316 + Ifges(4,4) * t319) * qJD(2);
t225 = Ifges(7,5) * t284 + Ifges(7,6) * t283 + Ifges(7,3) * t302;
t226 = Ifges(6,5) * t284 + Ifges(6,6) * t283 + Ifges(6,3) * t302;
t336 = -mrSges(7,1) * t200 + mrSges(7,3) * t194 + Ifges(7,4) * t245 + Ifges(7,2) * t244 + Ifges(7,6) * t280 + t302 * t229;
t167 = Ifges(6,4) * t245 + Ifges(6,2) * t244 + Ifges(6,6) * t280 + t302 * t230 - mrSges(6,1) * t204 + mrSges(6,3) * t197 - pkin(5) * (-mrSges(7,1) * t244 - t250 * t283 + t346) + qJ(6) * (-mrSges(7,2) * t280 - t253 * t302 + t347) + (-t226 - t225) * t284 + t336;
t334 = mrSges(7,2) * t200 - mrSges(7,3) * t191 + Ifges(7,1) * t245 + Ifges(7,4) * t244 + Ifges(7,5) * t280 + t283 * t225;
t175 = mrSges(6,2) * t204 - mrSges(6,3) * t196 + Ifges(6,1) * t245 + Ifges(6,4) * t244 + Ifges(6,5) * t280 - qJ(6) * t188 + t226 * t283 + (-t227 - t228) * t302 + t334;
t265 = Ifges(5,4) * qJD(3) + (-Ifges(5,2) * t316 - Ifges(5,6) * t319) * qJD(2);
t330 = -mrSges(5,2) * t210 + mrSges(5,3) * t208 - Ifges(5,1) * qJDD(3) + Ifges(5,4) * t288 + Ifges(5,5) * t289 + pkin(9) * t176 + t315 * t167 - t318 * t175 - t265 * t351;
t264 = Ifges(5,5) * qJD(3) + (-Ifges(5,6) * t316 - Ifges(5,3) * t319) * qJD(2);
t353 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t316 + Ifges(4,2) * t319) * qJD(2) - t264;
t371 = qJD(2) * (-t319 * t263 + t316 * t353) + mrSges(4,1) * t212 - mrSges(4,2) * t213 + Ifges(4,5) * t288 + Ifges(4,6) * t289 + Ifges(4,3) * qJDD(3) + pkin(3) * (-qJDD(3) * mrSges(5,2) - qJD(3) * t297 - t286 * t352 + t333) + qJ(4) * (mrSges(5,1) * t289 + t327) - t330;
t152 = -mrSges(3,1) * t255 + mrSges(3,3) * t224 + t322 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t166 - t371;
t342 = -t171 * t316 + t319 * t180;
t164 = m(3) * t224 - mrSges(3,1) * t322 - qJDD(2) * mrSges(3,2) + t342;
t220 = t331 - t365;
t177 = -t315 * t182 + t318 * t184;
t338 = -m(5) * t211 - t289 * mrSges(5,2) + t298 * t352 - t177;
t325 = -m(4) * t220 + t296 * t351 + t289 * mrSges(4,1) + (-mrSges(4,2) + mrSges(5,3)) * t288 + (-t295 * t316 - t297 * t319) * qJD(2) + t338;
t169 = m(3) * t223 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t322 + t325;
t161 = t320 * t164 - t169 * t317;
t372 = pkin(7) * t161 + t152 * t320;
t363 = Ifges(4,4) + Ifges(5,6);
t361 = t169 * t320;
t266 = Ifges(5,1) * qJD(3) + (-Ifges(5,4) * t316 - Ifges(5,5) * t319) * qJD(2);
t354 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t316 + Ifges(4,6) * t319) * qJD(2) + t266;
t165 = m(3) * t255 + t166;
t157 = t164 * t358 - t165 * t312 + t314 * t361;
t172 = -mrSges(5,3) * t288 + t297 * t351 - t338;
t329 = -mrSges(5,1) * t208 + mrSges(5,2) * t211 - pkin(4) * t369 - pkin(9) * t177 - t318 * t167 - t315 * t175;
t153 = -mrSges(4,1) * t220 + mrSges(4,3) * t213 - pkin(3) * t172 + (Ifges(4,2) + Ifges(5,3)) * t289 + t363 * t288 + (Ifges(4,6) - Ifges(5,5)) * qJDD(3) + (t263 - t265) * qJD(3) - t354 * t352 + t329;
t158 = t363 * t289 + (Ifges(4,1) + Ifges(5,2)) * t288 + (Ifges(4,5) - Ifges(5,4)) * qJDD(3) - t353 * qJD(3) + t354 * t351 + mrSges(4,2) * t220 - mrSges(4,3) * t212 - qJ(4) * t172 + t373;
t148 = mrSges(3,1) * t223 - mrSges(3,2) * t224 + Ifges(3,3) * qJDD(2) + pkin(2) * t325 + pkin(8) * t342 + t319 * t153 + t316 * t158;
t150 = mrSges(3,2) * t255 - mrSges(3,3) * t223 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t322 - pkin(8) * t166 - t153 * t316 + t158 * t319;
t332 = mrSges(2,1) * t293 - mrSges(2,2) * t294 + pkin(1) * t157 + t314 * t148 + t150 * t359 + t312 * t372;
t159 = m(2) * t294 + t161;
t156 = t314 * t165 + (t164 * t317 + t361) * t312;
t154 = m(2) * t293 + t157;
t146 = mrSges(2,2) * t310 - mrSges(2,3) * t293 + t320 * t150 - t317 * t152 + (-t156 * t312 - t157 * t314) * pkin(7);
t145 = -mrSges(2,1) * t310 + mrSges(2,3) * t294 - pkin(1) * t156 - t312 * t148 + (t150 * t317 + t372) * t314;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t313 * t146 - t311 * t145 - qJ(1) * (t154 * t313 + t159 * t311), t146, t150, t158, -t264 * t352 - t330, t175, -t227 * t302 + t334; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t311 * t146 + t313 * t145 + qJ(1) * (-t154 * t311 + t159 * t313), t145, t152, t153, Ifges(5,4) * qJDD(3) - Ifges(5,2) * t288 - Ifges(5,6) * t289 - qJD(3) * t264 - t266 * t351 - t373, t167, -t284 * t225 + t336; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t332, t332, t148, t371, Ifges(5,5) * qJDD(3) - Ifges(5,6) * t288 - Ifges(5,3) * t289 + qJD(3) * t265 + t266 * t352 - t329, t375, -t283 * t229 - t335;];
m_new  = t1;
