% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-05-06 01:54
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRRP9_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP9_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP9_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP9_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP9_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP9_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP9_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:50:53
% EndTime: 2019-05-06 01:51:09
% DurationCPUTime: 7.18s
% Computational Cost: add. (124388->335), mult. (242764->397), div. (0->0), fcn. (156800->8), ass. (0->128)
t304 = sin(qJ(4));
t308 = cos(qJ(4));
t309 = cos(qJ(3));
t338 = qJD(1) * t309;
t279 = t304 * qJD(3) + t308 * t338;
t305 = sin(qJ(3));
t337 = qJD(1) * qJD(3);
t333 = t305 * t337;
t283 = t309 * qJDD(1) - t333;
t240 = -t279 * qJD(4) + t308 * qJDD(3) - t304 * t283;
t278 = t308 * qJD(3) - t304 * t338;
t241 = t278 * qJD(4) + t304 * qJDD(3) + t308 * t283;
t303 = sin(qJ(5));
t307 = cos(qJ(5));
t243 = t307 * t278 - t303 * t279;
t208 = t243 * qJD(5) + t303 * t240 + t307 * t241;
t244 = t303 * t278 + t307 * t279;
t220 = -t243 * mrSges(7,1) + t244 * mrSges(7,2);
t312 = qJD(1) ^ 2;
t343 = -pkin(1) - pkin(7);
t306 = sin(qJ(1));
t310 = cos(qJ(1));
t288 = -t310 * g(1) - t306 * g(2);
t344 = -qJDD(1) * qJ(2) - 0.2e1 * qJD(2) * qJD(1) - t288;
t254 = t343 * t312 - t344;
t291 = t309 * t337;
t282 = -t305 * qJDD(1) - t291;
t224 = (-t283 + t333) * pkin(8) + (-t282 + t291) * pkin(3) + t254;
t287 = t306 * g(1) - t310 * g(2);
t327 = -t312 * qJ(2) + qJDD(2) - t287;
t255 = t343 * qJDD(1) + t327;
t248 = -t309 * g(3) + t305 * t255;
t281 = (pkin(3) * t305 - pkin(8) * t309) * qJD(1);
t293 = t305 * qJD(1);
t311 = qJD(3) ^ 2;
t228 = -t311 * pkin(3) + qJDD(3) * pkin(8) - t281 * t293 + t248;
t190 = t308 * t224 - t304 * t228;
t277 = qJDD(4) - t282;
t290 = t293 + qJD(4);
t186 = (t278 * t290 - t241) * pkin(9) + (t278 * t279 + t277) * pkin(4) + t190;
t191 = t304 * t224 + t308 * t228;
t253 = t290 * pkin(4) - t279 * pkin(9);
t276 = t278 ^ 2;
t188 = -t276 * pkin(4) + t240 * pkin(9) - t290 * t253 + t191;
t179 = t307 * t186 - t303 * t188;
t270 = qJDD(5) + t277;
t289 = qJD(5) + t290;
t174 = -0.2e1 * qJD(6) * t244 + (t243 * t289 - t208) * qJ(6) + (t243 * t244 + t270) * pkin(5) + t179;
t229 = -t289 * mrSges(7,2) + t243 * mrSges(7,3);
t335 = m(7) * t174 + t270 * mrSges(7,1) + t289 * t229;
t171 = -t208 * mrSges(7,3) - t244 * t220 + t335;
t180 = t303 * t186 + t307 * t188;
t207 = -t244 * qJD(5) + t307 * t240 - t303 * t241;
t214 = Ifges(6,4) * t244 + Ifges(6,2) * t243 + Ifges(6,6) * t289;
t215 = Ifges(7,1) * t244 + Ifges(7,4) * t243 + Ifges(7,5) * t289;
t216 = Ifges(6,1) * t244 + Ifges(6,4) * t243 + Ifges(6,5) * t289;
t231 = t289 * pkin(5) - t244 * qJ(6);
t242 = t243 ^ 2;
t177 = -t242 * pkin(5) + t207 * qJ(6) + 0.2e1 * qJD(6) * t243 - t289 * t231 + t180;
t213 = Ifges(7,4) * t244 + Ifges(7,2) * t243 + Ifges(7,6) * t289;
t325 = -mrSges(7,1) * t174 + mrSges(7,2) * t177 - Ifges(7,5) * t208 - Ifges(7,6) * t207 - Ifges(7,3) * t270 - t244 * t213;
t346 = mrSges(6,1) * t179 - mrSges(6,2) * t180 + Ifges(6,5) * t208 + Ifges(6,6) * t207 + Ifges(6,3) * t270 + pkin(5) * t171 + t244 * t214 - t325 + (-t216 - t215) * t243;
t221 = -t243 * mrSges(6,1) + t244 * mrSges(6,2);
t230 = -t289 * mrSges(6,2) + t243 * mrSges(6,3);
t163 = m(6) * t179 + t270 * mrSges(6,1) + t289 * t230 + (-t220 - t221) * t244 + (-mrSges(6,3) - mrSges(7,3)) * t208 + t335;
t232 = t289 * mrSges(7,1) - t244 * mrSges(7,3);
t233 = t289 * mrSges(6,1) - t244 * mrSges(6,3);
t334 = m(7) * t177 + t207 * mrSges(7,3) + t243 * t220;
t166 = m(6) * t180 + t207 * mrSges(6,3) + t243 * t221 + (-t232 - t233) * t289 + (-mrSges(6,2) - mrSges(7,2)) * t270 + t334;
t161 = t307 * t163 + t303 * t166;
t235 = Ifges(5,4) * t279 + Ifges(5,2) * t278 + Ifges(5,6) * t290;
t236 = Ifges(5,1) * t279 + Ifges(5,4) * t278 + Ifges(5,5) * t290;
t345 = mrSges(5,1) * t190 - mrSges(5,2) * t191 + Ifges(5,5) * t241 + Ifges(5,6) * t240 + Ifges(5,3) * t277 + pkin(4) * t161 + t279 * t235 - t278 * t236 + t346;
t342 = mrSges(2,1) - mrSges(3,2);
t341 = -Ifges(3,4) + Ifges(2,5);
t340 = Ifges(3,5) - Ifges(2,6);
t246 = -t278 * mrSges(5,1) + t279 * mrSges(5,2);
t249 = -t290 * mrSges(5,2) + t278 * mrSges(5,3);
t158 = m(5) * t190 + t277 * mrSges(5,1) - t241 * mrSges(5,3) - t279 * t246 + t290 * t249 + t161;
t250 = t290 * mrSges(5,1) - t279 * mrSges(5,3);
t330 = -t303 * t163 + t307 * t166;
t159 = m(5) * t191 - t277 * mrSges(5,2) + t240 * mrSges(5,3) + t278 * t246 - t290 * t250 + t330;
t152 = t308 * t158 + t304 * t159;
t280 = (mrSges(4,1) * t305 + mrSges(4,2) * t309) * qJD(1);
t286 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t338;
t331 = -t304 * t158 + t308 * t159;
t150 = m(4) * t248 - qJDD(3) * mrSges(4,2) + t282 * mrSges(4,3) - qJD(3) * t286 - t280 * t293 + t331;
t247 = t305 * g(3) + t309 * t255;
t285 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t293;
t227 = -qJDD(3) * pkin(3) - t311 * pkin(8) + t281 * t338 - t247;
t189 = -t240 * pkin(4) - t276 * pkin(9) + t279 * t253 + t227;
t183 = -t207 * pkin(5) - t242 * qJ(6) + t244 * t231 + qJDD(6) + t189;
t329 = m(7) * t183 - t207 * mrSges(7,1) + t208 * mrSges(7,2) - t243 * t229 + t244 * t232;
t318 = m(6) * t189 - t207 * mrSges(6,1) + t208 * mrSges(6,2) - t243 * t230 + t244 * t233 + t329;
t314 = -m(5) * t227 + t240 * mrSges(5,1) - t241 * mrSges(5,2) + t278 * t249 - t279 * t250 - t318;
t167 = m(4) * t247 + qJDD(3) * mrSges(4,1) - t283 * mrSges(4,3) + qJD(3) * t285 - t280 * t338 + t314;
t145 = t309 * t150 - t305 * t167;
t144 = t305 * t150 + t309 * t167;
t326 = -mrSges(7,1) * t183 + mrSges(7,3) * t177 + Ifges(7,4) * t208 + Ifges(7,2) * t207 + Ifges(7,6) * t270 + t289 * t215;
t211 = Ifges(7,5) * t244 + Ifges(7,6) * t243 + Ifges(7,3) * t289;
t324 = mrSges(7,2) * t183 - mrSges(7,3) * t174 + Ifges(7,1) * t208 + Ifges(7,4) * t207 + Ifges(7,5) * t270 + t243 * t211;
t260 = -qJDD(1) * pkin(1) + t327;
t323 = -m(3) * t260 + t312 * mrSges(3,3) - t144;
t148 = -m(4) * t254 + t282 * mrSges(4,1) - t283 * mrSges(4,2) - t285 * t293 - t286 * t338 - t152;
t212 = Ifges(6,5) * t244 + Ifges(6,6) * t243 + Ifges(6,3) * t289;
t154 = Ifges(6,4) * t208 + Ifges(6,2) * t207 + Ifges(6,6) * t270 + t289 * t216 - mrSges(6,1) * t189 + mrSges(6,3) * t180 - pkin(5) * t329 + qJ(6) * (-t270 * mrSges(7,2) - t289 * t232 + t334) + (-t212 - t211) * t244 + t326;
t160 = mrSges(6,2) * t189 - mrSges(6,3) * t179 + Ifges(6,1) * t208 + Ifges(6,4) * t207 + Ifges(6,5) * t270 - qJ(6) * t171 + t243 * t212 + (-t213 - t214) * t289 + t324;
t234 = Ifges(5,5) * t279 + Ifges(5,6) * t278 + Ifges(5,3) * t290;
t138 = -mrSges(5,1) * t227 + mrSges(5,3) * t191 + Ifges(5,4) * t241 + Ifges(5,2) * t240 + Ifges(5,6) * t277 - pkin(4) * t318 + pkin(9) * t330 + t307 * t154 + t303 * t160 - t279 * t234 + t290 * t236;
t140 = mrSges(5,2) * t227 - mrSges(5,3) * t190 + Ifges(5,1) * t241 + Ifges(5,4) * t240 + Ifges(5,5) * t277 - pkin(9) * t161 - t303 * t154 + t307 * t160 + t278 * t234 - t290 * t235;
t267 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t309 - Ifges(4,6) * t305) * qJD(1);
t268 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t309 - Ifges(4,2) * t305) * qJD(1);
t135 = mrSges(4,2) * t254 - mrSges(4,3) * t247 + Ifges(4,1) * t283 + Ifges(4,4) * t282 + Ifges(4,5) * qJDD(3) - pkin(8) * t152 - qJD(3) * t268 - t304 * t138 + t308 * t140 - t267 * t293;
t269 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t309 - Ifges(4,4) * t305) * qJD(1);
t136 = -mrSges(4,1) * t254 + mrSges(4,3) * t248 + Ifges(4,4) * t283 + Ifges(4,2) * t282 + Ifges(4,6) * qJDD(3) - pkin(3) * t152 + qJD(3) * t269 - t267 * t338 - t345;
t258 = t312 * pkin(1) + t344;
t322 = mrSges(3,2) * t260 - mrSges(3,3) * t258 + Ifges(3,1) * qJDD(1) - pkin(7) * t144 + t309 * t135 - t305 * t136;
t321 = -mrSges(3,1) * t258 - pkin(2) * t148 - pkin(7) * t145 - t305 * t135 - t309 * t136;
t320 = mrSges(4,1) * t247 - mrSges(4,2) * t248 + Ifges(4,5) * t283 + Ifges(4,6) * t282 + Ifges(4,3) * qJDD(3) + pkin(3) * t314 + pkin(8) * t331 + t308 * t138 + t304 * t140 + t268 * t338 + t269 * t293;
t319 = -m(3) * t258 + t312 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t148;
t316 = -mrSges(2,2) * t288 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t323) + qJ(2) * t319 + mrSges(2,1) * t287 + Ifges(2,3) * qJDD(1) + t322;
t315 = mrSges(3,1) * t260 + pkin(2) * t144 + t320;
t146 = m(2) * t288 - t312 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t319;
t143 = -m(3) * g(3) + t145;
t141 = m(2) * t287 - t312 * mrSges(2,2) + t342 * qJDD(1) + t323;
t133 = -mrSges(2,3) * t287 + t315 - qJ(2) * t143 + t340 * t312 + t341 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t132 = mrSges(2,3) * t288 - pkin(1) * t143 + t342 * g(3) - t340 * qJDD(1) + t341 * t312 + t321;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t310 * t133 - t306 * t132 - pkin(6) * (t310 * t141 + t306 * t146), t133, t322, t135, t140, t160, -t289 * t213 + t324; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t306 * t133 + t310 * t132 + pkin(6) * (-t306 * t141 + t310 * t146), t132, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - t312 * Ifges(3,5) - t315, t136, t138, t154, -t244 * t211 + t326; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t316, t316, mrSges(3,2) * g(3) + t312 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t321, t320, t345, t346, -t243 * t215 - t325;];
m_new  = t1;
