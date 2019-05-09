% Calculate vector of cutting torques with Newton-Euler for
% S6RPRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-05-05 18:46
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRPRR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:42:47
% EndTime: 2019-05-05 18:43:01
% DurationCPUTime: 7.63s
% Computational Cost: add. (133889->349), mult. (265811->417), div. (0->0), fcn. (153623->10), ass. (0->138)
t351 = -2 * qJD(4);
t309 = sin(qJ(1));
t313 = cos(qJ(1));
t286 = t309 * g(1) - g(2) * t313;
t271 = qJDD(1) * pkin(1) + t286;
t287 = -g(1) * t313 - g(2) * t309;
t315 = qJD(1) ^ 2;
t275 = -pkin(1) * t315 + t287;
t304 = sin(pkin(10));
t305 = cos(pkin(10));
t239 = t304 * t271 + t305 * t275;
t225 = -pkin(2) * t315 + qJDD(1) * pkin(7) + t239;
t308 = sin(qJ(3));
t222 = t308 * t225;
t303 = -g(3) + qJDD(2);
t312 = cos(qJ(3));
t343 = t312 * t303;
t218 = -t222 + t343;
t219 = t312 * t225 + t308 * t303;
t254 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t308 + Ifges(4,4) * t312) * qJD(1);
t273 = (mrSges(5,2) * t312 - mrSges(5,3) * t308) * qJD(1);
t339 = qJD(1) * qJD(3);
t337 = t312 * t339;
t276 = qJDD(1) * t308 + t337;
t336 = t308 * t339;
t277 = qJDD(1) * t312 - t336;
t340 = qJD(1) * t312;
t283 = -mrSges(5,1) * t340 - qJD(3) * mrSges(5,3);
t294 = t308 * qJD(1);
t272 = (-pkin(3) * t312 - qJ(4) * t308) * qJD(1);
t314 = qJD(3) ^ 2;
t213 = pkin(3) * t314 - qJDD(3) * qJ(4) + qJD(3) * t351 - t272 * t340 - t219;
t285 = pkin(4) * t294 - qJD(3) * pkin(8);
t302 = t312 ^ 2;
t347 = pkin(8) * t315;
t206 = pkin(4) * t277 + qJD(3) * t285 - t302 * t347 - t213;
t307 = sin(qJ(5));
t311 = cos(qJ(5));
t270 = qJD(3) * t311 - t307 * t340;
t233 = -qJD(5) * t270 - qJDD(3) * t307 - t277 * t311;
t269 = -qJD(3) * t307 - t311 * t340;
t234 = qJD(5) * t269 + qJDD(3) * t311 - t277 * t307;
t290 = t294 + qJD(5);
t241 = -mrSges(6,2) * t290 + mrSges(6,3) * t269;
t242 = mrSges(6,1) * t290 - mrSges(6,3) * t270;
t243 = pkin(5) * t290 - pkin(9) * t270;
t267 = t269 ^ 2;
t194 = -pkin(5) * t233 - pkin(9) * t267 + t243 * t270 + t206;
t306 = sin(qJ(6));
t310 = cos(qJ(6));
t237 = t269 * t306 + t270 * t310;
t201 = -qJD(6) * t237 + t233 * t310 - t234 * t306;
t236 = t269 * t310 - t270 * t306;
t202 = qJD(6) * t236 + t233 * t306 + t234 * t310;
t288 = qJD(6) + t290;
t220 = -mrSges(7,2) * t288 + mrSges(7,3) * t236;
t221 = mrSges(7,1) * t288 - mrSges(7,3) * t237;
t328 = m(7) * t194 - t201 * mrSges(7,1) + t202 * mrSges(7,2) - t220 * t236 + t237 * t221;
t184 = -m(6) * t206 + mrSges(6,1) * t233 - t234 * mrSges(6,2) + t241 * t269 - t270 * t242 - t328;
t284 = mrSges(5,1) * t294 + qJD(3) * mrSges(5,2);
t319 = -m(5) * t213 + qJDD(3) * mrSges(5,3) + qJD(3) * t284 + t273 * t340 - t184;
t238 = t305 * t271 - t304 * t275;
t331 = -qJDD(1) * pkin(2) - t238;
t321 = pkin(3) * t336 + t294 * t351 + (-t276 - t337) * qJ(4) + t331;
t348 = -pkin(3) - pkin(8);
t196 = -t285 * t294 + (-pkin(4) * t302 - pkin(7)) * t315 + t348 * t277 + t321;
t332 = -t314 * qJ(4) + t272 * t294 + qJDD(4) + t222;
t207 = t276 * pkin(4) + t348 * qJDD(3) + (-pkin(4) * t339 - t308 * t347 - t303) * t312 + t332;
t191 = -t307 * t196 + t311 * t207;
t268 = qJDD(5) + t276;
t188 = (t269 * t290 - t234) * pkin(9) + (t269 * t270 + t268) * pkin(5) + t191;
t192 = t311 * t196 + t307 * t207;
t189 = -pkin(5) * t267 + pkin(9) * t233 - t243 * t290 + t192;
t187 = t188 * t306 + t189 * t310;
t209 = Ifges(7,5) * t237 + Ifges(7,6) * t236 + Ifges(7,3) * t288;
t211 = Ifges(7,1) * t237 + Ifges(7,4) * t236 + Ifges(7,5) * t288;
t259 = qJDD(6) + t268;
t173 = -mrSges(7,1) * t194 + mrSges(7,3) * t187 + Ifges(7,4) * t202 + Ifges(7,2) * t201 + Ifges(7,6) * t259 - t209 * t237 + t211 * t288;
t186 = t188 * t310 - t189 * t306;
t210 = Ifges(7,4) * t237 + Ifges(7,2) * t236 + Ifges(7,6) * t288;
t174 = mrSges(7,2) * t194 - mrSges(7,3) * t186 + Ifges(7,1) * t202 + Ifges(7,4) * t201 + Ifges(7,5) * t259 + t209 * t236 - t210 * t288;
t226 = Ifges(6,5) * t270 + Ifges(6,6) * t269 + Ifges(6,3) * t290;
t228 = Ifges(6,1) * t270 + Ifges(6,4) * t269 + Ifges(6,5) * t290;
t216 = -mrSges(7,1) * t236 + mrSges(7,2) * t237;
t182 = m(7) * t186 + mrSges(7,1) * t259 - t202 * mrSges(7,3) - t216 * t237 + t220 * t288;
t183 = m(7) * t187 - mrSges(7,2) * t259 + t201 * mrSges(7,3) + t216 * t236 - t221 * t288;
t333 = -t182 * t306 + t310 * t183;
t155 = -mrSges(6,1) * t206 + mrSges(6,3) * t192 + Ifges(6,4) * t234 + Ifges(6,2) * t233 + Ifges(6,6) * t268 - pkin(5) * t328 + pkin(9) * t333 + t310 * t173 + t306 * t174 - t270 * t226 + t290 * t228;
t172 = t310 * t182 + t306 * t183;
t227 = Ifges(6,4) * t270 + Ifges(6,2) * t269 + Ifges(6,6) * t290;
t157 = mrSges(6,2) * t206 - mrSges(6,3) * t191 + Ifges(6,1) * t234 + Ifges(6,4) * t233 + Ifges(6,5) * t268 - pkin(9) * t172 - t173 * t306 + t174 * t310 + t226 * t269 - t227 * t290;
t240 = -mrSges(6,1) * t269 + mrSges(6,2) * t270;
t169 = m(6) * t191 + mrSges(6,1) * t268 - mrSges(6,3) * t234 - t240 * t270 + t241 * t290 + t172;
t170 = m(6) * t192 - mrSges(6,2) * t268 + mrSges(6,3) * t233 + t240 * t269 - t242 * t290 + t333;
t165 = t169 * t311 + t170 * t307;
t215 = -qJDD(3) * pkin(3) + t332 - t343;
t256 = Ifges(5,4) * qJD(3) + (-Ifges(5,2) * t308 - Ifges(5,6) * t312) * qJD(1);
t324 = -mrSges(5,2) * t215 + mrSges(5,3) * t213 - Ifges(5,1) * qJDD(3) + Ifges(5,4) * t276 + Ifges(5,5) * t277 + pkin(8) * t165 + t307 * t155 - t311 * t157 - t256 * t340;
t327 = -m(5) * t215 - t276 * mrSges(5,1) - t165;
t255 = Ifges(5,5) * qJD(3) + (-Ifges(5,6) * t308 - Ifges(5,3) * t312) * qJD(1);
t341 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t308 + Ifges(4,2) * t312) * qJD(1) - t255;
t350 = (-t312 * t254 + t341 * t308) * qJD(1) + mrSges(4,1) * t218 - mrSges(4,2) * t219 + Ifges(4,5) * t276 + Ifges(4,6) * t277 + Ifges(4,3) * qJDD(3) + pkin(3) * (-qJDD(3) * mrSges(5,2) - qJD(3) * t283 - t273 * t294 + t327) + qJ(4) * (mrSges(5,1) * t277 + t319) - t324;
t346 = t315 * pkin(7);
t345 = Ifges(4,4) + Ifges(5,6);
t274 = (-mrSges(4,1) * t312 + mrSges(4,2) * t308) * qJD(1);
t282 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t340;
t162 = m(4) * t218 - mrSges(4,3) * t276 + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t282 - t283) * qJD(3) + (-t273 - t274) * t294 + t327;
t281 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t294;
t177 = t274 * t340 + m(4) * t219 - qJD(3) * t281 + (mrSges(4,3) + mrSges(5,1)) * t277 - qJDD(3) * mrSges(4,2) + t319;
t334 = -t162 * t308 + t312 * t177;
t152 = m(3) * t239 - mrSges(3,1) * t315 - qJDD(1) * mrSges(3,2) + t334;
t224 = t331 - t346;
t166 = -t307 * t169 + t311 * t170;
t208 = -t277 * pkin(3) + t321 - t346;
t330 = -m(5) * t208 - t277 * mrSges(5,2) + t284 * t294 - t166;
t318 = -m(4) * t224 + t282 * t340 + t277 * mrSges(4,1) + (-mrSges(4,2) + mrSges(5,3)) * t276 + (-t281 * t308 - t283 * t312) * qJD(1) + t330;
t159 = m(3) * t238 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t315 + t318;
t149 = t304 * t152 + t305 * t159;
t154 = t312 * t162 + t308 * t177;
t257 = Ifges(5,1) * qJD(3) + (-Ifges(5,4) * t308 - Ifges(5,5) * t312) * qJD(1);
t342 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t308 + Ifges(4,6) * t312) * qJD(1) + t257;
t335 = t305 * t152 - t159 * t304;
t163 = -mrSges(5,3) * t276 + t283 * t340 - t330;
t322 = -mrSges(5,1) * t213 + mrSges(5,2) * t208 - pkin(4) * t184 - pkin(8) * t166 - t311 * t155 - t307 * t157;
t143 = -mrSges(4,1) * t224 + mrSges(4,3) * t219 - pkin(3) * t163 + (Ifges(4,2) + Ifges(5,3)) * t277 + t345 * t276 + (Ifges(4,6) - Ifges(5,5)) * qJDD(3) + (t254 - t256) * qJD(3) - t342 * t294 + t322;
t325 = -mrSges(7,1) * t186 + mrSges(7,2) * t187 - Ifges(7,5) * t202 - Ifges(7,6) * t201 - Ifges(7,3) * t259 - t237 * t210 + t236 * t211;
t320 = -mrSges(6,1) * t191 + mrSges(6,2) * t192 - Ifges(6,5) * t234 - Ifges(6,6) * t233 - Ifges(6,3) * t268 - pkin(5) * t172 - t270 * t227 + t269 * t228 + t325;
t317 = -mrSges(5,1) * t215 + mrSges(5,3) * t208 - pkin(4) * t165 + t320;
t145 = t345 * t277 + (Ifges(4,1) + Ifges(5,2)) * t276 + (Ifges(4,5) - Ifges(5,4)) * qJDD(3) - t341 * qJD(3) - t317 + t342 * t340 - mrSges(4,3) * t218 + mrSges(4,2) * t224 - qJ(4) * t163;
t326 = mrSges(3,1) * t238 - mrSges(3,2) * t239 + Ifges(3,3) * qJDD(1) + pkin(2) * t318 + pkin(7) * t334 + t312 * t143 + t308 * t145;
t323 = mrSges(2,1) * t286 - mrSges(2,2) * t287 + Ifges(2,3) * qJDD(1) + pkin(1) * t149 + t326;
t147 = m(2) * t287 - mrSges(2,1) * t315 - qJDD(1) * mrSges(2,2) + t335;
t146 = m(2) * t286 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t315 + t149;
t141 = -mrSges(3,1) * t303 + mrSges(3,3) * t239 + t315 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t154 - t350;
t140 = mrSges(3,2) * t303 - mrSges(3,3) * t238 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t315 - pkin(7) * t154 - t143 * t308 + t145 * t312;
t139 = -mrSges(2,2) * g(3) - mrSges(2,3) * t286 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t315 - qJ(2) * t149 + t140 * t305 - t141 * t304;
t138 = Ifges(2,6) * qJDD(1) + t315 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t287 + t304 * t140 + t305 * t141 - pkin(1) * (m(3) * t303 + t154) + qJ(2) * t335;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t313 * t139 - t309 * t138 - pkin(6) * (t146 * t313 + t147 * t309), t139, t140, t145, -t255 * t294 - t324, t157, t174; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t309 * t139 + t313 * t138 + pkin(6) * (-t146 * t309 + t147 * t313), t138, t141, t143, Ifges(5,4) * qJDD(3) - Ifges(5,2) * t276 - Ifges(5,6) * t277 - qJD(3) * t255 - t257 * t340 + t317, t155, t173; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t323, t323, t326, t350, Ifges(5,5) * qJDD(3) - Ifges(5,6) * t276 - Ifges(5,3) * t277 + qJD(3) * t256 + t257 * t294 - t322, -t320, -t325;];
m_new  = t1;
