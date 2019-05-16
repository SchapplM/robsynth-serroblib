% Calculate vector of cutting torques with Newton-Euler for
% S6RPRPRR2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-05-05 18:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRPRR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:22:31
% EndTime: 2019-05-05 18:22:58
% DurationCPUTime: 18.64s
% Computational Cost: add. (343872->341), mult. (746529->432), div. (0->0), fcn. (508237->12), ass. (0->137)
t307 = sin(qJ(1));
t311 = cos(qJ(1));
t289 = t307 * g(1) - t311 * g(2);
t280 = qJDD(1) * pkin(1) + t289;
t290 = -t311 * g(1) - t307 * g(2);
t313 = qJD(1) ^ 2;
t282 = -t313 * pkin(1) + t290;
t301 = sin(pkin(10));
t303 = cos(pkin(10));
t256 = t301 * t280 + t303 * t282;
t244 = -t313 * pkin(2) + qJDD(1) * pkin(7) + t256;
t299 = -g(3) + qJDD(2);
t306 = sin(qJ(3));
t310 = cos(qJ(3));
t234 = -t306 * t244 + t310 * t299;
t333 = qJD(1) * qJD(3);
t332 = t310 * t333;
t283 = t306 * qJDD(1) + t332;
t217 = (-t283 + t332) * qJ(4) + (t306 * t310 * t313 + qJDD(3)) * pkin(3) + t234;
t235 = t310 * t244 + t306 * t299;
t284 = t310 * qJDD(1) - t306 * t333;
t336 = qJD(1) * t306;
t286 = qJD(3) * pkin(3) - qJ(4) * t336;
t298 = t310 ^ 2;
t220 = -t298 * t313 * pkin(3) + t284 * qJ(4) - qJD(3) * t286 + t235;
t300 = sin(pkin(11));
t302 = cos(pkin(11));
t269 = (t300 * t310 + t302 * t306) * qJD(1);
t204 = -0.2e1 * qJD(4) * t269 + t302 * t217 - t300 * t220;
t335 = qJD(1) * t310;
t268 = -t300 * t336 + t302 * t335;
t205 = 0.2e1 * qJD(4) * t268 + t300 * t217 + t302 * t220;
t246 = -t268 * mrSges(5,1) + t269 * mrSges(5,2);
t257 = -t300 * t283 + t302 * t284;
t263 = qJD(3) * mrSges(5,1) - t269 * mrSges(5,3);
t248 = -t268 * pkin(4) - t269 * pkin(8);
t312 = qJD(3) ^ 2;
t196 = -t312 * pkin(4) + qJDD(3) * pkin(8) + t268 * t248 + t205;
t255 = t303 * t280 - t301 * t282;
t324 = -qJDD(1) * pkin(2) - t255;
t221 = -t284 * pkin(3) + qJDD(4) + t286 * t336 + (-qJ(4) * t298 - pkin(7)) * t313 + t324;
t258 = t302 * t283 + t300 * t284;
t208 = (-qJD(3) * t268 - t258) * pkin(8) + (qJD(3) * t269 - t257) * pkin(4) + t221;
t305 = sin(qJ(5));
t309 = cos(qJ(5));
t191 = -t305 * t196 + t309 * t208;
t260 = t309 * qJD(3) - t305 * t269;
t228 = t260 * qJD(5) + t305 * qJDD(3) + t309 * t258;
t254 = qJDD(5) - t257;
t261 = t305 * qJD(3) + t309 * t269;
t267 = qJD(5) - t268;
t189 = (t260 * t267 - t228) * pkin(9) + (t260 * t261 + t254) * pkin(5) + t191;
t192 = t309 * t196 + t305 * t208;
t227 = -t261 * qJD(5) + t309 * qJDD(3) - t305 * t258;
t238 = t267 * pkin(5) - t261 * pkin(9);
t259 = t260 ^ 2;
t190 = -t259 * pkin(5) + t227 * pkin(9) - t267 * t238 + t192;
t304 = sin(qJ(6));
t308 = cos(qJ(6));
t187 = t308 * t189 - t304 * t190;
t229 = t308 * t260 - t304 * t261;
t201 = t229 * qJD(6) + t304 * t227 + t308 * t228;
t230 = t304 * t260 + t308 * t261;
t213 = -t229 * mrSges(7,1) + t230 * mrSges(7,2);
t264 = qJD(6) + t267;
t218 = -t264 * mrSges(7,2) + t229 * mrSges(7,3);
t249 = qJDD(6) + t254;
t182 = m(7) * t187 + t249 * mrSges(7,1) - t201 * mrSges(7,3) - t230 * t213 + t264 * t218;
t188 = t304 * t189 + t308 * t190;
t200 = -t230 * qJD(6) + t308 * t227 - t304 * t228;
t219 = t264 * mrSges(7,1) - t230 * mrSges(7,3);
t183 = m(7) * t188 - t249 * mrSges(7,2) + t200 * mrSges(7,3) + t229 * t213 - t264 * t219;
t174 = t308 * t182 + t304 * t183;
t232 = -t260 * mrSges(6,1) + t261 * mrSges(6,2);
t236 = -t267 * mrSges(6,2) + t260 * mrSges(6,3);
t172 = m(6) * t191 + t254 * mrSges(6,1) - t228 * mrSges(6,3) - t261 * t232 + t267 * t236 + t174;
t237 = t267 * mrSges(6,1) - t261 * mrSges(6,3);
t327 = -t304 * t182 + t308 * t183;
t173 = m(6) * t192 - t254 * mrSges(6,2) + t227 * mrSges(6,3) + t260 * t232 - t267 * t237 + t327;
t328 = -t305 * t172 + t309 * t173;
t165 = m(5) * t205 - qJDD(3) * mrSges(5,2) + t257 * mrSges(5,3) - qJD(3) * t263 + t268 * t246 + t328;
t262 = -qJD(3) * mrSges(5,2) + t268 * mrSges(5,3);
t195 = -qJDD(3) * pkin(4) - t312 * pkin(8) + t269 * t248 - t204;
t193 = -t227 * pkin(5) - t259 * pkin(9) + t261 * t238 + t195;
t322 = m(7) * t193 - t200 * mrSges(7,1) + t201 * mrSges(7,2) - t229 * t218 + t230 * t219;
t317 = -m(6) * t195 + t227 * mrSges(6,1) - t228 * mrSges(6,2) + t260 * t236 - t261 * t237 - t322;
t178 = m(5) * t204 + qJDD(3) * mrSges(5,1) - t258 * mrSges(5,3) + qJD(3) * t262 - t269 * t246 + t317;
t155 = t300 * t165 + t302 * t178;
t274 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t306 + Ifges(4,2) * t310) * qJD(1);
t275 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t306 + Ifges(4,4) * t310) * qJD(1);
t209 = Ifges(7,5) * t230 + Ifges(7,6) * t229 + Ifges(7,3) * t264;
t211 = Ifges(7,1) * t230 + Ifges(7,4) * t229 + Ifges(7,5) * t264;
t175 = -mrSges(7,1) * t193 + mrSges(7,3) * t188 + Ifges(7,4) * t201 + Ifges(7,2) * t200 + Ifges(7,6) * t249 - t230 * t209 + t264 * t211;
t210 = Ifges(7,4) * t230 + Ifges(7,2) * t229 + Ifges(7,6) * t264;
t176 = mrSges(7,2) * t193 - mrSges(7,3) * t187 + Ifges(7,1) * t201 + Ifges(7,4) * t200 + Ifges(7,5) * t249 + t229 * t209 - t264 * t210;
t222 = Ifges(6,5) * t261 + Ifges(6,6) * t260 + Ifges(6,3) * t267;
t224 = Ifges(6,1) * t261 + Ifges(6,4) * t260 + Ifges(6,5) * t267;
t157 = -mrSges(6,1) * t195 + mrSges(6,3) * t192 + Ifges(6,4) * t228 + Ifges(6,2) * t227 + Ifges(6,6) * t254 - pkin(5) * t322 + pkin(9) * t327 + t308 * t175 + t304 * t176 - t261 * t222 + t267 * t224;
t223 = Ifges(6,4) * t261 + Ifges(6,2) * t260 + Ifges(6,6) * t267;
t159 = mrSges(6,2) * t195 - mrSges(6,3) * t191 + Ifges(6,1) * t228 + Ifges(6,4) * t227 + Ifges(6,5) * t254 - pkin(9) * t174 - t304 * t175 + t308 * t176 + t260 * t222 - t267 * t223;
t241 = Ifges(5,4) * t269 + Ifges(5,2) * t268 + Ifges(5,6) * qJD(3);
t242 = Ifges(5,1) * t269 + Ifges(5,4) * t268 + Ifges(5,5) * qJD(3);
t318 = -mrSges(5,1) * t204 + mrSges(5,2) * t205 - Ifges(5,5) * t258 - Ifges(5,6) * t257 - Ifges(5,3) * qJDD(3) - pkin(4) * t317 - pkin(8) * t328 - t309 * t157 - t305 * t159 - t269 * t241 + t268 * t242;
t337 = mrSges(4,1) * t234 - mrSges(4,2) * t235 + Ifges(4,5) * t283 + Ifges(4,6) * t284 + Ifges(4,3) * qJDD(3) + pkin(3) * t155 + (t306 * t274 - t310 * t275) * qJD(1) - t318;
t281 = (-mrSges(4,1) * t310 + mrSges(4,2) * t306) * qJD(1);
t288 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t335;
t153 = m(4) * t234 + qJDD(3) * mrSges(4,1) - t283 * mrSges(4,3) + qJD(3) * t288 - t281 * t336 + t155;
t287 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t336;
t329 = t302 * t165 - t300 * t178;
t154 = m(4) * t235 - qJDD(3) * mrSges(4,2) + t284 * mrSges(4,3) - qJD(3) * t287 + t281 * t335 + t329;
t330 = -t306 * t153 + t310 * t154;
t146 = m(3) * t256 - t313 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t330;
t243 = -t313 * pkin(7) + t324;
t167 = t309 * t172 + t305 * t173;
t320 = m(5) * t221 - t257 * mrSges(5,1) + t258 * mrSges(5,2) - t268 * t262 + t269 * t263 + t167;
t316 = -m(4) * t243 + t284 * mrSges(4,1) - t283 * mrSges(4,2) - t287 * t336 + t288 * t335 - t320;
t161 = m(3) * t255 + qJDD(1) * mrSges(3,1) - t313 * mrSges(3,2) + t316;
t142 = t301 * t146 + t303 * t161;
t148 = t310 * t153 + t306 * t154;
t331 = t303 * t146 - t301 * t161;
t240 = Ifges(5,5) * t269 + Ifges(5,6) * t268 + Ifges(5,3) * qJD(3);
t143 = mrSges(5,2) * t221 - mrSges(5,3) * t204 + Ifges(5,1) * t258 + Ifges(5,4) * t257 + Ifges(5,5) * qJDD(3) - pkin(8) * t167 - qJD(3) * t241 - t305 * t157 + t309 * t159 + t268 * t240;
t321 = -mrSges(7,1) * t187 + mrSges(7,2) * t188 - Ifges(7,5) * t201 - Ifges(7,6) * t200 - Ifges(7,3) * t249 - t230 * t210 + t229 * t211;
t314 = mrSges(6,1) * t191 - mrSges(6,2) * t192 + Ifges(6,5) * t228 + Ifges(6,6) * t227 + Ifges(6,3) * t254 + pkin(5) * t174 + t261 * t223 - t260 * t224 - t321;
t149 = -mrSges(5,1) * t221 + mrSges(5,3) * t205 + Ifges(5,4) * t258 + Ifges(5,2) * t257 + Ifges(5,6) * qJDD(3) - pkin(4) * t167 + qJD(3) * t242 - t269 * t240 - t314;
t273 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t306 + Ifges(4,6) * t310) * qJD(1);
t135 = -mrSges(4,1) * t243 + mrSges(4,3) * t235 + Ifges(4,4) * t283 + Ifges(4,2) * t284 + Ifges(4,6) * qJDD(3) - pkin(3) * t320 + qJ(4) * t329 + qJD(3) * t275 + t300 * t143 + t302 * t149 - t273 * t336;
t138 = mrSges(4,2) * t243 - mrSges(4,3) * t234 + Ifges(4,1) * t283 + Ifges(4,4) * t284 + Ifges(4,5) * qJDD(3) - qJ(4) * t155 - qJD(3) * t274 + t302 * t143 - t300 * t149 + t273 * t335;
t323 = mrSges(3,1) * t255 - mrSges(3,2) * t256 + Ifges(3,3) * qJDD(1) + pkin(2) * t316 + pkin(7) * t330 + t310 * t135 + t306 * t138;
t319 = mrSges(2,1) * t289 - mrSges(2,2) * t290 + Ifges(2,3) * qJDD(1) + pkin(1) * t142 + t323;
t140 = m(2) * t290 - t313 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t331;
t139 = m(2) * t289 + qJDD(1) * mrSges(2,1) - t313 * mrSges(2,2) + t142;
t136 = -mrSges(3,1) * t299 + mrSges(3,3) * t256 + t313 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t148 - t337;
t133 = mrSges(3,2) * t299 - mrSges(3,3) * t255 + Ifges(3,5) * qJDD(1) - t313 * Ifges(3,6) - pkin(7) * t148 - t306 * t135 + t310 * t138;
t132 = -mrSges(2,2) * g(3) - mrSges(2,3) * t289 + Ifges(2,5) * qJDD(1) - t313 * Ifges(2,6) - qJ(2) * t142 + t303 * t133 - t301 * t136;
t131 = Ifges(2,6) * qJDD(1) + t313 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t290 + t301 * t133 + t303 * t136 - pkin(1) * (m(3) * t299 + t148) + qJ(2) * t331;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t311 * t132 - t307 * t131 - pkin(6) * (t311 * t139 + t307 * t140), t132, t133, t138, t143, t159, t176; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t307 * t132 + t311 * t131 + pkin(6) * (-t307 * t139 + t311 * t140), t131, t136, t135, t149, t157, t175; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t319, t319, t323, t337, -t318, t314, -t321;];
m_new  = t1;
