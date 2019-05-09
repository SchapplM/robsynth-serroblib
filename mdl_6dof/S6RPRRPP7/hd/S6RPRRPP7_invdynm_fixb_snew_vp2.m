% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRPP7
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% Datum: 2019-05-05 21:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRPP7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP7_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP7_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP7_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP7_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP7_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:46:23
% EndTime: 2019-05-05 21:46:29
% DurationCPUTime: 3.38s
% Computational Cost: add. (41322->338), mult. (78102->378), div. (0->0), fcn. (44014->6), ass. (0->123)
t308 = sin(qJ(1));
t310 = cos(qJ(1));
t291 = -t310 * g(1) - t308 * g(2);
t359 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t291;
t312 = qJD(1) ^ 2;
t352 = (-pkin(1) - pkin(7));
t254 = (t352 * t312) - t359;
t307 = sin(qJ(3));
t309 = cos(qJ(3));
t338 = qJD(1) * qJD(3);
t334 = t309 * t338;
t284 = -t307 * qJDD(1) - t334;
t335 = t307 * t338;
t285 = t309 * qJDD(1) - t335;
t189 = (-t285 + t335) * pkin(8) + (-t284 + t334) * pkin(3) + t254;
t290 = t308 * g(1) - t310 * g(2);
t328 = -t312 * qJ(2) + qJDD(2) - t290;
t255 = t352 * qJDD(1) + t328;
t243 = -t309 * g(3) + t307 * t255;
t283 = (pkin(3) * t307 - pkin(8) * t309) * qJD(1);
t311 = qJD(3) ^ 2;
t339 = t307 * qJD(1);
t197 = -t311 * pkin(3) + qJDD(3) * pkin(8) - t283 * t339 + t243;
t306 = sin(qJ(4));
t351 = cos(qJ(4));
t187 = t306 * t189 + t351 * t197;
t340 = qJD(1) * t309;
t280 = -t351 * qJD(3) + t306 * t340;
t281 = t306 * qJD(3) + t351 * t340;
t238 = t280 * pkin(4) - t281 * qJ(5);
t279 = qJDD(4) - t284;
t293 = qJD(4) + t339;
t292 = t293 ^ 2;
t353 = 2 * qJD(5);
t182 = -t292 * pkin(4) + t279 * qJ(5) - t280 * t238 + t293 * t353 + t187;
t249 = -t293 * mrSges(6,1) + t281 * mrSges(6,2);
t358 = m(6) * t182 + t279 * mrSges(6,3) + t293 * t249;
t186 = t351 * t189 - t306 * t197;
t184 = -t279 * pkin(4) - t292 * qJ(5) + t281 * t238 + qJDD(5) - t186;
t250 = -t280 * mrSges(6,2) + t293 * mrSges(6,3);
t357 = -m(6) * t184 + t279 * mrSges(6,1) + t293 * t250;
t234 = -t280 * qJD(4) + t306 * qJDD(3) + t351 * t285;
t242 = t307 * g(3) + t309 * t255;
t325 = qJDD(3) * pkin(3) + t311 * pkin(8) - t283 * t340 + t242;
t345 = t280 * t293;
t356 = (-t234 + t345) * qJ(5) - t325;
t244 = t293 * mrSges(7,2) + t280 * mrSges(7,3);
t354 = -0.2e1 * t281;
t172 = qJD(6) * t354 + (-t234 - t345) * qJ(6) + (t280 * t281 - t279) * pkin(5) + t184;
t240 = -t280 * mrSges(7,1) + t281 * mrSges(7,2);
t331 = -m(7) * t172 + t234 * mrSges(7,3) + t281 * t240;
t168 = -t279 * mrSges(7,1) - t293 * t244 - t331;
t247 = -t293 * mrSges(7,1) - t281 * mrSges(7,3);
t233 = t281 * qJD(4) - t351 * qJDD(3) + t306 * t285;
t246 = -t293 * pkin(5) - t281 * qJ(6);
t278 = t280 ^ 2;
t176 = -t278 * pkin(5) + t233 * qJ(6) + 0.2e1 * qJD(6) * t280 + t293 * t246 + t182;
t336 = m(7) * t176 + t233 * mrSges(7,3) + t280 * t240;
t170 = t279 * mrSges(7,2) + t293 * t247 + t336;
t205 = Ifges(6,5) * t281 + Ifges(6,6) * t293 + Ifges(6,3) * t280;
t209 = Ifges(5,4) * t281 - Ifges(5,2) * t280 + Ifges(5,6) * t293;
t212 = Ifges(5,1) * t281 - Ifges(5,4) * t280 + Ifges(5,5) * t293;
t239 = t280 * mrSges(6,1) - t281 * mrSges(6,3);
t211 = Ifges(6,1) * t281 + Ifges(6,4) * t293 + Ifges(6,5) * t280;
t207 = Ifges(7,4) * t281 + Ifges(7,2) * t280 - Ifges(7,6) * t293;
t210 = Ifges(7,1) * t281 + Ifges(7,4) * t280 - Ifges(7,5) * t293;
t329 = mrSges(7,1) * t172 - mrSges(7,2) * t176 + Ifges(7,5) * t234 + Ifges(7,6) * t233 - Ifges(7,3) * t279 + t281 * t207 - t280 * t210;
t318 = mrSges(6,1) * t184 - mrSges(6,3) * t182 - Ifges(6,4) * t234 - Ifges(6,2) * t279 - Ifges(6,6) * t233 + pkin(5) * t168 - t280 * t211 + t329;
t355 = (t209 - t205) * t281 + mrSges(5,1) * t186 - mrSges(5,2) * t187 + Ifges(5,5) * t234 - Ifges(5,6) * t233 + Ifges(5,3) * t279 + pkin(4) * (-t234 * mrSges(6,2) - t281 * t239 - t168 + t357) + qJ(5) * (-t233 * mrSges(6,2) - t280 * t239 + t170 + t358) + t280 * t212 - t318;
t349 = mrSges(2,1) - mrSges(3,2);
t348 = -mrSges(5,3) - mrSges(6,2);
t347 = -Ifges(3,4) + Ifges(2,5);
t346 = (Ifges(3,5) - Ifges(2,6));
t344 = t293 * t207;
t245 = -t293 * mrSges(5,2) - t280 * mrSges(5,3);
t341 = -t280 * mrSges(5,1) - t281 * mrSges(5,2) - t239;
t161 = m(5) * t186 + (t244 + t245) * t293 + t341 * t281 + (mrSges(5,1) + mrSges(7,1)) * t279 + t348 * t234 + t331 + t357;
t248 = t293 * mrSges(5,1) - t281 * mrSges(5,3);
t162 = m(5) * t187 + (t247 - t248) * t293 + t341 * t280 + (-mrSges(5,2) + mrSges(7,2)) * t279 + t348 * t233 + t336 + t358;
t155 = t351 * t161 + t306 * t162;
t208 = Ifges(6,4) * t281 + Ifges(6,2) * t293 + Ifges(6,6) * t280;
t343 = -Ifges(5,5) * t281 + Ifges(5,6) * t280 - Ifges(5,3) * t293 - t208;
t282 = (mrSges(4,1) * t307 + mrSges(4,2) * t309) * qJD(1);
t288 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t340;
t333 = -t306 * t161 + t351 * t162;
t153 = m(4) * t243 - qJDD(3) * mrSges(4,2) + t284 * mrSges(4,3) - qJD(3) * t288 - t282 * t339 + t333;
t287 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t339;
t179 = -t278 * qJ(6) + qJDD(6) + (-pkin(4) - pkin(5)) * t233 + (-pkin(4) * t293 + t246 + t353) * t281 - t356;
t169 = -m(7) * t179 + t233 * mrSges(7,1) - t234 * mrSges(7,2) + t280 * t244 - t281 * t247;
t185 = qJD(5) * t354 + (t281 * t293 + t233) * pkin(4) + t356;
t165 = m(6) * t185 + t233 * mrSges(6,1) - t234 * mrSges(6,3) - t281 * t249 + t280 * t250 + t169;
t314 = m(5) * t325 - t233 * mrSges(5,1) - t234 * mrSges(5,2) - t280 * t245 - t281 * t248 - t165;
t157 = m(4) * t242 + qJDD(3) * mrSges(4,1) - t285 * mrSges(4,3) + qJD(3) * t287 - t282 * t340 + t314;
t146 = t309 * t153 - t307 * t157;
t145 = t307 * t153 + t309 * t157;
t327 = mrSges(7,1) * t179 - mrSges(7,3) * t176 - Ifges(7,4) * t234 - Ifges(7,2) * t233 + Ifges(7,6) * t279 + t293 * t210;
t204 = Ifges(7,5) * t281 + Ifges(7,6) * t280 - Ifges(7,3) * t293;
t326 = mrSges(7,2) * t179 - mrSges(7,3) * t172 + Ifges(7,1) * t234 + Ifges(7,4) * t233 - Ifges(7,5) * t279 + t280 * t204;
t260 = -qJDD(1) * pkin(1) + t328;
t324 = -m(3) * t260 + (t312 * mrSges(3,3)) - t145;
t151 = -m(4) * t254 + t284 * mrSges(4,1) - t285 * mrSges(4,2) - t287 * t339 - t288 * t340 - t155;
t317 = mrSges(6,1) * t185 - mrSges(6,2) * t182 + pkin(5) * t169 + qJ(6) * t170 - t327;
t141 = (t212 + t211) * t293 + (t204 + t343) * t281 + (Ifges(5,6) - Ifges(6,6)) * t279 + (Ifges(5,4) - Ifges(6,5)) * t234 + (-Ifges(5,2) - Ifges(6,3)) * t233 + mrSges(5,1) * t325 + mrSges(5,3) * t187 - pkin(4) * t165 - t317;
t315 = mrSges(6,2) * t184 - mrSges(6,3) * t185 + Ifges(6,1) * t234 + Ifges(6,4) * t279 + Ifges(6,5) * t233 - qJ(6) * t168 + t293 * t205 + t326;
t148 = (-t209 + t207) * t293 + t343 * t280 + Ifges(5,5) * t279 - Ifges(5,4) * t233 + Ifges(5,1) * t234 - mrSges(5,2) * t325 - mrSges(5,3) * t186 - qJ(5) * t165 + t315;
t262 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t309 - Ifges(4,6) * t307) * qJD(1);
t263 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t309 - Ifges(4,2) * t307) * qJD(1);
t138 = mrSges(4,2) * t254 - mrSges(4,3) * t242 + Ifges(4,1) * t285 + Ifges(4,4) * t284 + Ifges(4,5) * qJDD(3) - pkin(8) * t155 - qJD(3) * t263 - t306 * t141 + t351 * t148 - t262 * t339;
t264 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t309 - Ifges(4,4) * t307) * qJD(1);
t139 = -mrSges(4,1) * t254 + mrSges(4,3) * t243 + Ifges(4,4) * t285 + Ifges(4,2) * t284 + Ifges(4,6) * qJDD(3) - pkin(3) * t155 + qJD(3) * t264 - t262 * t340 - t355;
t258 = t312 * pkin(1) + t359;
t323 = mrSges(3,2) * t260 - mrSges(3,3) * t258 + Ifges(3,1) * qJDD(1) - pkin(7) * t145 + t309 * t138 - t307 * t139;
t322 = -mrSges(3,1) * t258 - pkin(2) * t151 - pkin(7) * t146 - t307 * t138 - t309 * t139;
t321 = mrSges(4,1) * t242 - mrSges(4,2) * t243 + Ifges(4,5) * t285 + Ifges(4,6) * t284 + Ifges(4,3) * qJDD(3) + pkin(3) * t314 + pkin(8) * t333 + t351 * t141 + t306 * t148 + t263 * t340 + t264 * t339;
t320 = -m(3) * t258 + t312 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t151;
t319 = -mrSges(2,2) * t291 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t324) + qJ(2) * t320 + mrSges(2,1) * t290 + Ifges(2,3) * qJDD(1) + t323;
t316 = mrSges(3,1) * t260 + pkin(2) * t145 + t321;
t149 = m(2) * t291 - t312 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t320;
t144 = -m(3) * g(3) + t146;
t142 = m(2) * t290 - t312 * mrSges(2,2) + t349 * qJDD(1) + t324;
t136 = -mrSges(2,3) * t290 - qJ(2) * t144 + t316 + (t346 * t312) + t347 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t135 = mrSges(2,3) * t291 - pkin(1) * t144 + t349 * g(3) - t346 * qJDD(1) + t347 * t312 + t322;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t310 * t136 - t308 * t135 - pkin(6) * (t310 * t142 + t308 * t149), t136, t323, t138, t148, -t280 * t208 + t315 + t344, t326 + t344; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t308 * t136 + t310 * t135 + pkin(6) * (-t308 * t142 + t310 * t149), t135, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - (t312 * Ifges(3,5)) - t316, t139, t141, -t281 * t205 - t318, -t281 * t204 - t327; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t319, t319, mrSges(3,2) * g(3) + t312 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t322, t321, t355, -t293 * t211 + Ifges(6,6) * t279 + Ifges(6,5) * t234 + Ifges(6,3) * t233 + t317 + (-t204 + t208) * t281, t329;];
m_new  = t1;
