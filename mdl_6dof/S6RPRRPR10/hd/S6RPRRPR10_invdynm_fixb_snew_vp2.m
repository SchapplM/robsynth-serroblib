% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
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
% Datum: 2019-05-05 23:57
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRPR10_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR10_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR10_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR10_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR10_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR10_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR10_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 23:53:09
% EndTime: 2019-05-05 23:53:23
% DurationCPUTime: 6.05s
% Computational Cost: add. (84205->339), mult. (160992->397), div. (0->0), fcn. (97785->8), ass. (0->131)
t310 = sin(qJ(1));
t313 = cos(qJ(1));
t291 = -t313 * g(1) - t310 * g(2);
t355 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t291;
t308 = sin(qJ(4));
t312 = cos(qJ(3));
t340 = qJD(1) * t312;
t350 = cos(qJ(4));
t280 = -qJD(3) * t350 + t308 * t340;
t309 = sin(qJ(3));
t339 = qJD(1) * qJD(3);
t337 = t309 * t339;
t285 = qJDD(1) * t312 - t337;
t236 = -qJD(4) * t280 + qJDD(3) * t308 + t285 * t350;
t290 = t310 * g(1) - g(2) * t313;
t315 = qJD(1) ^ 2;
t330 = -t315 * qJ(2) + qJDD(2) - t290;
t351 = -pkin(1) - pkin(7);
t256 = qJDD(1) * t351 + t330;
t245 = g(3) * t309 + t256 * t312;
t283 = (pkin(3) * t309 - pkin(8) * t312) * qJD(1);
t314 = qJD(3) ^ 2;
t329 = qJDD(3) * pkin(3) + t314 * pkin(8) - t283 * t340 + t245;
t341 = qJD(1) * t309;
t294 = qJD(4) + t341;
t345 = t280 * t294;
t354 = (-t236 + t345) * qJ(5) - t329;
t255 = t315 * t351 - t355;
t336 = t312 * t339;
t284 = -qJDD(1) * t309 - t336;
t210 = (-t285 + t337) * pkin(8) + (-t284 + t336) * pkin(3) + t255;
t246 = -g(3) * t312 + t256 * t309;
t215 = -pkin(3) * t314 + qJDD(3) * pkin(8) - t283 * t341 + t246;
t193 = t210 * t350 - t215 * t308;
t194 = t210 * t308 + t215 * t350;
t281 = qJD(3) * t308 + t340 * t350;
t221 = Ifges(6,5) * t281 + Ifges(6,6) * t294 + Ifges(6,3) * t280;
t224 = Ifges(5,4) * t281 - Ifges(5,2) * t280 + Ifges(5,6) * t294;
t226 = Ifges(5,1) * t281 - Ifges(5,4) * t280 + Ifges(5,5) * t294;
t235 = qJD(4) * t281 - qJDD(3) * t350 + t285 * t308;
t243 = mrSges(6,1) * t280 - mrSges(6,3) * t281;
t279 = qJDD(4) - t284;
t242 = pkin(4) * t280 - qJ(5) * t281;
t293 = t294 ^ 2;
t191 = -t279 * pkin(4) - t293 * qJ(5) + t242 * t281 + qJDD(5) - t193;
t183 = (-t236 - t345) * pkin(9) + (t280 * t281 - t279) * pkin(5) + t191;
t352 = 2 * qJD(5);
t189 = -pkin(4) * t293 + qJ(5) * t279 - t242 * t280 + t294 * t352 + t194;
t254 = -pkin(5) * t294 - pkin(9) * t281;
t278 = t280 ^ 2;
t184 = -pkin(5) * t278 + pkin(9) * t235 + t254 * t294 + t189;
t307 = sin(qJ(6));
t311 = cos(qJ(6));
t181 = t183 * t311 - t184 * t307;
t237 = t280 * t311 - t281 * t307;
t202 = qJD(6) * t237 + t235 * t307 + t236 * t311;
t238 = t280 * t307 + t281 * t311;
t208 = -mrSges(7,1) * t237 + mrSges(7,2) * t238;
t292 = qJD(6) - t294;
t217 = -mrSges(7,2) * t292 + mrSges(7,3) * t237;
t272 = qJDD(6) - t279;
t176 = m(7) * t181 + mrSges(7,1) * t272 - mrSges(7,3) * t202 - t208 * t238 + t217 * t292;
t182 = t183 * t307 + t184 * t311;
t201 = -qJD(6) * t238 + t235 * t311 - t236 * t307;
t218 = mrSges(7,1) * t292 - mrSges(7,3) * t238;
t177 = m(7) * t182 - mrSges(7,2) * t272 + mrSges(7,3) * t201 + t208 * t237 - t218 * t292;
t166 = t311 * t176 + t307 * t177;
t225 = Ifges(6,1) * t281 + Ifges(6,4) * t294 + Ifges(6,5) * t280;
t204 = Ifges(7,4) * t238 + Ifges(7,2) * t237 + Ifges(7,6) * t292;
t205 = Ifges(7,1) * t238 + Ifges(7,4) * t237 + Ifges(7,5) * t292;
t331 = mrSges(7,1) * t181 - mrSges(7,2) * t182 + Ifges(7,5) * t202 + Ifges(7,6) * t201 + Ifges(7,3) * t272 + t204 * t238 - t205 * t237;
t319 = mrSges(6,1) * t191 - mrSges(6,3) * t189 - Ifges(6,4) * t236 - Ifges(6,2) * t279 - Ifges(6,6) * t235 + pkin(5) * t166 - t225 * t280 + t331;
t250 = -mrSges(6,2) * t280 + mrSges(6,3) * t294;
t327 = -m(6) * t191 + mrSges(6,1) * t279 + t250 * t294 - t166;
t167 = -t307 * t176 + t177 * t311;
t249 = -mrSges(6,1) * t294 + mrSges(6,2) * t281;
t332 = m(6) * t189 + mrSges(6,3) * t279 + t249 * t294 + t167;
t353 = (t224 - t221) * t281 + mrSges(5,1) * t193 - mrSges(5,2) * t194 + Ifges(5,5) * t236 - Ifges(5,6) * t235 + Ifges(5,3) * t279 + pkin(4) * (-t236 * mrSges(6,2) - t281 * t243 + t327) + qJ(5) * (-t235 * mrSges(6,2) - t280 * t243 + t332) + t280 * t226 - t319;
t349 = mrSges(2,1) - mrSges(3,2);
t348 = -mrSges(5,3) - mrSges(6,2);
t347 = -Ifges(3,4) + Ifges(2,5);
t346 = (Ifges(3,5) - Ifges(2,6));
t248 = mrSges(5,1) * t294 - mrSges(5,3) * t281;
t342 = -mrSges(5,1) * t280 - mrSges(5,2) * t281 - t243;
t161 = m(5) * t194 - t279 * mrSges(5,2) + t235 * t348 - t294 * t248 + t280 * t342 + t332;
t247 = -mrSges(5,2) * t294 - mrSges(5,3) * t280;
t163 = m(5) * t193 + t279 * mrSges(5,1) + t236 * t348 + t294 * t247 + t281 * t342 + t327;
t157 = t161 * t308 + t163 * t350;
t223 = Ifges(6,4) * t281 + Ifges(6,2) * t294 + Ifges(6,6) * t280;
t344 = -Ifges(5,5) * t281 + Ifges(5,6) * t280 - Ifges(5,3) * t294 - t223;
t282 = (mrSges(4,1) * t309 + mrSges(4,2) * t312) * qJD(1);
t288 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t340;
t335 = t161 * t350 - t163 * t308;
t155 = m(4) * t246 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t284 - qJD(3) * t288 - t282 * t341 + t335;
t287 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t341;
t186 = -t278 * pkin(9) + (-pkin(4) - pkin(5)) * t235 + (-pkin(4) * t294 + t254 + t352) * t281 - t354;
t178 = -m(7) * t186 + mrSges(7,1) * t201 - mrSges(7,2) * t202 + t217 * t237 - t218 * t238;
t192 = -0.2e1 * qJD(5) * t281 + (t281 * t294 + t235) * pkin(4) + t354;
t174 = m(6) * t192 + mrSges(6,1) * t235 - mrSges(6,3) * t236 - t249 * t281 + t250 * t280 + t178;
t317 = m(5) * t329 - mrSges(5,1) * t235 - mrSges(5,2) * t236 - t247 * t280 - t248 * t281 - t174;
t171 = m(4) * t245 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t285 + qJD(3) * t287 - t282 * t340 + t317;
t150 = t155 * t312 - t171 * t309;
t149 = t309 * t155 + t312 * t171;
t261 = -qJDD(1) * pkin(1) + t330;
t328 = -m(3) * t261 + (mrSges(3,3) * t315) - t149;
t153 = -m(4) * t255 + mrSges(4,1) * t284 - mrSges(4,2) * t285 - t287 * t341 - t288 * t340 - t157;
t203 = Ifges(7,5) * t238 + Ifges(7,6) * t237 + Ifges(7,3) * t292;
t169 = -mrSges(7,1) * t186 + mrSges(7,3) * t182 + Ifges(7,4) * t202 + Ifges(7,2) * t201 + Ifges(7,6) * t272 - t203 * t238 + t205 * t292;
t170 = mrSges(7,2) * t186 - mrSges(7,3) * t181 + Ifges(7,1) * t202 + Ifges(7,4) * t201 + Ifges(7,5) * t272 + t203 * t237 - t204 * t292;
t321 = -mrSges(6,1) * t192 + mrSges(6,2) * t189 - pkin(5) * t178 - pkin(9) * t167 - t311 * t169 - t307 * t170;
t143 = mrSges(5,1) * t329 + mrSges(5,3) * t194 - pkin(4) * t174 + (t226 + t225) * t294 + t344 * t281 + (Ifges(5,6) - Ifges(6,6)) * t279 + (Ifges(5,4) - Ifges(6,5)) * t236 + (-Ifges(5,2) - Ifges(6,3)) * t235 + t321;
t322 = mrSges(6,2) * t191 - mrSges(6,3) * t192 + Ifges(6,1) * t236 + Ifges(6,4) * t279 + Ifges(6,5) * t235 - pkin(9) * t166 - t169 * t307 + t170 * t311 + t221 * t294;
t145 = -mrSges(5,2) * t329 - mrSges(5,3) * t193 + Ifges(5,1) * t236 - Ifges(5,4) * t235 + Ifges(5,5) * t279 - qJ(5) * t174 - t224 * t294 + t280 * t344 + t322;
t264 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t312 - Ifges(4,6) * t309) * qJD(1);
t265 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t312 - Ifges(4,2) * t309) * qJD(1);
t140 = mrSges(4,2) * t255 - mrSges(4,3) * t245 + Ifges(4,1) * t285 + Ifges(4,4) * t284 + Ifges(4,5) * qJDD(3) - pkin(8) * t157 - qJD(3) * t265 - t143 * t308 + t145 * t350 - t264 * t341;
t266 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t312 - Ifges(4,4) * t309) * qJD(1);
t141 = -mrSges(4,1) * t255 + mrSges(4,3) * t246 + Ifges(4,4) * t285 + Ifges(4,2) * t284 + Ifges(4,6) * qJDD(3) - pkin(3) * t157 + qJD(3) * t266 - t264 * t340 - t353;
t259 = t315 * pkin(1) + t355;
t326 = mrSges(3,2) * t261 - mrSges(3,3) * t259 + Ifges(3,1) * qJDD(1) - pkin(7) * t149 + t140 * t312 - t141 * t309;
t325 = -mrSges(3,1) * t259 - pkin(2) * t153 - pkin(7) * t150 - t309 * t140 - t312 * t141;
t324 = mrSges(4,1) * t245 - mrSges(4,2) * t246 + Ifges(4,5) * t285 + Ifges(4,6) * t284 + Ifges(4,3) * qJDD(3) + pkin(3) * t317 + pkin(8) * t335 + t143 * t350 + t145 * t308 + t265 * t340 + t266 * t341;
t323 = -m(3) * t259 + mrSges(3,2) * t315 + qJDD(1) * mrSges(3,3) - t153;
t320 = -mrSges(2,2) * t291 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t328) + qJ(2) * t323 + mrSges(2,1) * t290 + Ifges(2,3) * qJDD(1) + t326;
t318 = mrSges(3,1) * t261 + pkin(2) * t149 + t324;
t151 = m(2) * t291 - mrSges(2,1) * t315 - qJDD(1) * mrSges(2,2) + t323;
t148 = -m(3) * g(3) + t150;
t146 = m(2) * t290 - t315 * mrSges(2,2) + qJDD(1) * t349 + t328;
t138 = (t346 * t315) + t347 * qJDD(1) - qJ(2) * t148 + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t290 + t318;
t137 = mrSges(2,3) * t291 - pkin(1) * t148 + g(3) * t349 - qJDD(1) * t346 + t315 * t347 + t325;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t313 * t138 - t310 * t137 - pkin(6) * (t146 * t313 + t151 * t310), t138, t326, t140, t145, -t223 * t280 + t322, t170; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t310 * t138 + t313 * t137 + pkin(6) * (-t146 * t310 + t151 * t313), t137, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - (Ifges(3,5) * t315) - t318, t141, t143, -t281 * t221 - t319, t169; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t320, t320, mrSges(3,2) * g(3) + t315 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t325, t324, t353, Ifges(6,5) * t236 + Ifges(6,6) * t279 + Ifges(6,3) * t235 + t281 * t223 - t294 * t225 - t321, t331;];
m_new  = t1;
