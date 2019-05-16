% Calculate vector of cutting torques with Newton-Euler for
% S6PRRPRP5
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
% Datum: 2019-05-05 04:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRPRP5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:06:59
% EndTime: 2019-05-05 04:07:12
% DurationCPUTime: 6.80s
% Computational Cost: add. (87506->344), mult. (170478->407), div. (0->0), fcn. (100770->10), ass. (0->134)
t363 = -2 * qJD(4);
t304 = sin(pkin(10));
t306 = cos(pkin(10));
t285 = g(1) * t304 - g(2) * t306;
t286 = -g(1) * t306 - g(2) * t304;
t303 = -g(3) + qJDD(1);
t313 = cos(qJ(2));
t307 = cos(pkin(6));
t310 = sin(qJ(2));
t350 = t307 * t310;
t305 = sin(pkin(6));
t351 = t305 * t310;
t218 = t285 * t350 + t313 * t286 + t303 * t351;
t315 = qJD(2) ^ 2;
t216 = -pkin(2) * t315 + qJDD(2) * pkin(8) + t218;
t309 = sin(qJ(3));
t212 = t309 * t216;
t247 = -t285 * t305 + t303 * t307;
t312 = cos(qJ(3));
t349 = t312 * t247;
t208 = -t212 + t349;
t278 = (mrSges(5,2) * t312 - mrSges(5,3) * t309) * qJD(2);
t279 = (-mrSges(4,1) * t312 + mrSges(4,2) * t309) * qJD(2);
t341 = qJD(2) * qJD(3);
t338 = t312 * t341;
t280 = qJDD(2) * t309 + t338;
t342 = qJD(2) * t312;
t288 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t342;
t289 = -mrSges(5,1) * t342 - qJD(3) * mrSges(5,3);
t277 = (-pkin(3) * t312 - qJ(4) * t309) * qJD(2);
t314 = qJD(3) ^ 2;
t343 = qJD(2) * t309;
t332 = -t314 * qJ(4) + t277 * t343 + qJDD(4) + t212;
t357 = pkin(9) * t315;
t359 = -pkin(3) - pkin(9);
t201 = t280 * pkin(4) + t359 * qJDD(3) + (-pkin(4) * t341 - t309 * t357 - t247) * t312 + t332;
t337 = t309 * t341;
t281 = qJDD(2) * t312 - t337;
t292 = pkin(4) * t343 - qJD(3) * pkin(9);
t302 = t312 ^ 2;
t217 = -t310 * t286 + (t285 * t307 + t303 * t305) * t313;
t325 = -qJDD(2) * pkin(2) - t217;
t321 = pkin(3) * t337 + t343 * t363 + (-t280 - t338) * qJ(4) + t325;
t203 = -t292 * t343 + (-pkin(4) * t302 - pkin(8)) * t315 + t359 * t281 + t321;
t308 = sin(qJ(5));
t311 = cos(qJ(5));
t196 = t308 * t201 + t311 * t203;
t276 = qJD(3) * t311 - t308 * t342;
t235 = qJD(5) * t276 + qJDD(3) * t308 + t311 * t281;
t295 = qJD(5) + t343;
t245 = mrSges(6,1) * t295 - mrSges(6,3) * t276;
t272 = qJDD(5) + t280;
t275 = qJD(3) * t308 + t311 * t342;
t239 = pkin(5) * t275 - qJ(6) * t276;
t293 = t295 ^ 2;
t190 = -pkin(5) * t293 + qJ(6) * t272 + 0.2e1 * qJD(6) * t295 - t239 * t275 + t196;
t246 = -mrSges(7,1) * t295 + mrSges(7,2) * t276;
t339 = m(7) * t190 + t272 * mrSges(7,3) + t295 * t246;
t240 = mrSges(7,1) * t275 - mrSges(7,3) * t276;
t346 = -mrSges(6,1) * t275 - mrSges(6,2) * t276 - t240;
t355 = -mrSges(6,3) - mrSges(7,2);
t178 = m(6) * t196 - t272 * mrSges(6,2) + t355 * t235 - t295 * t245 + t346 * t275 + t339;
t195 = t201 * t311 - t203 * t308;
t236 = -qJD(5) * t275 + qJDD(3) * t311 - t281 * t308;
t243 = -mrSges(6,2) * t295 - mrSges(6,3) * t275;
t192 = -pkin(5) * t272 - qJ(6) * t293 + t239 * t276 + qJDD(6) - t195;
t244 = -mrSges(7,2) * t275 + mrSges(7,3) * t295;
t335 = -m(7) * t192 + t272 * mrSges(7,1) + t295 * t244;
t180 = m(6) * t195 + t272 * mrSges(6,1) + t355 * t236 + t295 * t243 + t346 * t276 + t335;
t172 = t308 * t178 + t311 * t180;
t206 = -qJDD(3) * pkin(3) + t332 - t349;
t328 = -m(5) * t206 - t280 * mrSges(5,1) - t172;
t166 = m(4) * t208 - t280 * mrSges(4,3) + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t288 - t289) * qJD(3) + (-t278 - t279) * t343 + t328;
t209 = t312 * t216 + t309 * t247;
t287 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t343;
t204 = t314 * pkin(3) - qJDD(3) * qJ(4) + qJD(3) * t363 - t277 * t342 - t209;
t290 = mrSges(5,1) * t343 + qJD(3) * mrSges(5,2);
t200 = t281 * pkin(4) + qJD(3) * t292 - t302 * t357 - t204;
t197 = -0.2e1 * qJD(6) * t276 + (t275 * t295 - t236) * qJ(6) + (t276 * t295 + t235) * pkin(5) + t200;
t187 = m(7) * t197 + t235 * mrSges(7,1) - mrSges(7,3) * t236 + t275 * t244 - t246 * t276;
t322 = m(6) * t200 + mrSges(6,1) * t235 + t236 * mrSges(6,2) + t243 * t275 + t276 * t245 + t187;
t320 = -m(5) * t204 + qJDD(3) * mrSges(5,3) + qJD(3) * t290 + t278 * t342 + t322;
t176 = t279 * t342 + (mrSges(4,3) + mrSges(5,1)) * t281 + m(4) * t209 - qJD(3) * t287 - qJDD(3) * mrSges(4,2) + t320;
t162 = t312 * t166 + t309 * t176;
t255 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t309 + Ifges(4,4) * t312) * qJD(2);
t223 = Ifges(7,1) * t276 + Ifges(7,4) * t295 + Ifges(7,5) * t275;
t224 = Ifges(6,1) * t276 - Ifges(6,4) * t275 + Ifges(6,5) * t295;
t334 = -mrSges(7,1) * t197 + mrSges(7,2) * t190;
t221 = Ifges(7,4) * t276 + Ifges(7,2) * t295 + Ifges(7,6) * t275;
t347 = -Ifges(6,5) * t276 + Ifges(6,6) * t275 - Ifges(6,3) * t295 - t221;
t169 = -mrSges(6,1) * t200 + mrSges(6,3) * t196 - pkin(5) * t187 + (t223 + t224) * t295 + t347 * t276 + (Ifges(6,6) - Ifges(7,6)) * t272 + (Ifges(6,4) - Ifges(7,5)) * t236 + (-Ifges(6,2) - Ifges(7,3)) * t235 + t334;
t222 = Ifges(6,4) * t276 - Ifges(6,2) * t275 + Ifges(6,6) * t295;
t219 = Ifges(7,5) * t276 + Ifges(7,6) * t295 + Ifges(7,3) * t275;
t329 = mrSges(7,2) * t192 - mrSges(7,3) * t197 + Ifges(7,1) * t236 + Ifges(7,4) * t272 + Ifges(7,5) * t235 + t295 * t219;
t171 = mrSges(6,2) * t200 - mrSges(6,3) * t195 + Ifges(6,1) * t236 - Ifges(6,4) * t235 + Ifges(6,5) * t272 - qJ(6) * t187 - t295 * t222 + t347 * t275 + t329;
t257 = Ifges(5,4) * qJD(3) + (-Ifges(5,2) * t309 - Ifges(5,6) * t312) * qJD(2);
t324 = -mrSges(5,2) * t206 + mrSges(5,3) * t204 - Ifges(5,1) * qJDD(3) + Ifges(5,4) * t280 + Ifges(5,5) * t281 + pkin(9) * t172 + t308 * t169 - t311 * t171 - t257 * t342;
t256 = Ifges(5,5) * qJD(3) + (-Ifges(5,6) * t309 - Ifges(5,3) * t312) * qJD(2);
t344 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t309 + Ifges(4,2) * t312) * qJD(2) - t256;
t361 = qJD(2) * (-t312 * t255 + t344 * t309) + mrSges(4,1) * t208 - mrSges(4,2) * t209 + Ifges(4,5) * t280 + Ifges(4,6) * t281 + Ifges(4,3) * qJDD(3) + pkin(3) * (-qJDD(3) * mrSges(5,2) - qJD(3) * t289 - t278 * t343 + t328) + qJ(4) * (mrSges(5,1) * t281 + t320) - t324;
t148 = -mrSges(3,1) * t247 + mrSges(3,3) * t218 + t315 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t162 - t361;
t336 = -t166 * t309 + t312 * t176;
t160 = m(3) * t218 - mrSges(3,1) * t315 - qJDD(2) * mrSges(3,2) + t336;
t356 = t315 * pkin(8);
t215 = t325 - t356;
t173 = t311 * t178 - t308 * t180;
t207 = -t281 * pkin(3) + t321 - t356;
t331 = -m(5) * t207 - t281 * mrSges(5,2) + t290 * t343 - t173;
t319 = -m(4) * t215 + t288 * t342 + t281 * mrSges(4,1) + (-mrSges(4,2) + mrSges(5,3)) * t280 + (-t287 * t309 - t289 * t312) * qJD(2) + t331;
t164 = m(3) * t217 + qJDD(2) * mrSges(3,1) - t315 * mrSges(3,2) + t319;
t157 = t313 * t160 - t164 * t310;
t362 = pkin(7) * t157 + t148 * t313;
t354 = Ifges(4,4) + Ifges(5,6);
t352 = t164 * t313;
t258 = Ifges(5,1) * qJD(3) + (-Ifges(5,4) * t309 - Ifges(5,5) * t312) * qJD(2);
t345 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t309 + Ifges(4,6) * t312) * qJD(2) + t258;
t161 = m(3) * t247 + t162;
t153 = t160 * t350 - t161 * t305 + t307 * t352;
t167 = -t280 * mrSges(5,3) + t289 * t342 - t331;
t323 = -mrSges(5,1) * t204 + mrSges(5,2) * t207 + pkin(4) * t322 - pkin(9) * t173 - t311 * t169 - t308 * t171;
t149 = -mrSges(4,1) * t215 + mrSges(4,3) * t209 - pkin(3) * t167 + (Ifges(4,2) + Ifges(5,3)) * t281 + t354 * t280 + (Ifges(4,6) - Ifges(5,5)) * qJDD(3) + (t255 - t257) * qJD(3) - t345 * t343 + t323;
t326 = mrSges(7,1) * t192 - mrSges(7,3) * t190 - Ifges(7,4) * t236 - Ifges(7,2) * t272 - Ifges(7,6) * t235 + t276 * t219 - t275 * t223;
t318 = mrSges(6,2) * t196 - t275 * t224 - qJ(6) * (-t235 * mrSges(7,2) - t275 * t240 + t339) - pkin(5) * (-t236 * mrSges(7,2) - t276 * t240 + t335) - mrSges(6,1) * t195 - t276 * t222 + Ifges(6,6) * t235 - Ifges(6,5) * t236 - Ifges(6,3) * t272 + t326;
t317 = -mrSges(5,1) * t206 + mrSges(5,3) * t207 - pkin(4) * t172 + t318;
t154 = t345 * t342 + t354 * t281 + (Ifges(4,1) + Ifges(5,2)) * t280 + (Ifges(4,5) - Ifges(5,4)) * qJDD(3) - t344 * qJD(3) - t317 + mrSges(4,2) * t215 - mrSges(4,3) * t208 - qJ(4) * t167;
t144 = mrSges(3,1) * t217 - mrSges(3,2) * t218 + Ifges(3,3) * qJDD(2) + pkin(2) * t319 + pkin(8) * t336 + t312 * t149 + t309 * t154;
t146 = mrSges(3,2) * t247 - mrSges(3,3) * t217 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t315 - pkin(8) * t162 - t149 * t309 + t154 * t312;
t327 = mrSges(2,1) * t285 - mrSges(2,2) * t286 + pkin(1) * t153 + t307 * t144 + t146 * t351 + t305 * t362;
t155 = m(2) * t286 + t157;
t152 = t307 * t161 + (t160 * t310 + t352) * t305;
t150 = m(2) * t285 + t153;
t142 = mrSges(2,2) * t303 - mrSges(2,3) * t285 + t313 * t146 - t310 * t148 + (-t152 * t305 - t153 * t307) * pkin(7);
t141 = -mrSges(2,1) * t303 + mrSges(2,3) * t286 - pkin(1) * t152 - t305 * t144 + (t146 * t310 + t362) * t307;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t306 * t142 - t304 * t141 - qJ(1) * (t150 * t306 + t155 * t304), t142, t146, t154, -t256 * t343 - t324, t171, -t221 * t275 + t329; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t304 * t142 + t306 * t141 + qJ(1) * (-t150 * t304 + t155 * t306), t141, t148, t149, Ifges(5,4) * qJDD(3) - Ifges(5,2) * t280 - Ifges(5,6) * t281 - qJD(3) * t256 - t258 * t342 + t317, t169, -t326; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t327, t327, t144, t361, Ifges(5,5) * qJDD(3) - Ifges(5,6) * t280 - Ifges(5,3) * t281 + qJD(3) * t257 + t258 * t343 - t323, -t318, Ifges(7,5) * t236 + Ifges(7,6) * t272 + Ifges(7,3) * t235 + t276 * t221 - t295 * t223 - t334;];
m_new  = t1;
