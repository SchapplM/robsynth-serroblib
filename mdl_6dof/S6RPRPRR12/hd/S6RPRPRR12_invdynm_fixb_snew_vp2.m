% Calculate vector of cutting torques with Newton-Euler for
% S6RPRPRR12
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
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
% Datum: 2019-05-05 20:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRPRR12_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR12_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR12_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR12_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR12_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR12_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR12_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 20:44:46
% EndTime: 2019-05-05 20:45:00
% DurationCPUTime: 6.35s
% Computational Cost: add. (86159->347), mult. (172175->408), div. (0->0), fcn. (95987->8), ass. (0->138)
t360 = -2 * qJD(4);
t301 = sin(qJ(1));
t305 = cos(qJ(1));
t275 = t301 * g(1) - t305 * g(2);
t307 = qJD(1) ^ 2;
t327 = -t307 * qJ(2) + qJDD(2) - t275;
t354 = -pkin(1) - pkin(7);
t233 = t354 * qJDD(1) + t327;
t304 = cos(qJ(3));
t343 = t304 * t233;
t300 = sin(qJ(3));
t350 = t300 * g(3);
t225 = t343 + t350;
t341 = qJD(1) * qJD(3);
t282 = t300 * t341;
t266 = t304 * qJDD(1) - t282;
t342 = qJD(1) * t300;
t270 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t342;
t272 = mrSges(5,1) * t342 - qJD(3) * mrSges(5,3);
t337 = t304 * t341;
t265 = t300 * qJDD(1) + t337;
t283 = t304 * qJD(1);
t274 = pkin(4) * t283 - qJD(3) * pkin(8);
t295 = t300 ^ 2;
t276 = -t305 * g(1) - t301 * g(2);
t329 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t276;
t318 = pkin(3) * t337 + t283 * t360 + t329 + (-t266 + t282) * qJ(4);
t353 = pkin(3) + pkin(8);
t193 = -t274 * t283 + t353 * t265 + (-pkin(4) * t295 + t354) * t307 + t318;
t262 = (pkin(3) * t300 - qJ(4) * t304) * qJD(1);
t306 = qJD(3) ^ 2;
t326 = -t306 * qJ(4) + t262 * t283 + qJDD(4) - t343;
t351 = pkin(8) * t307;
t198 = t266 * pkin(4) - t353 * qJDD(3) + (pkin(4) * t341 + t304 * t351 - g(3)) * t300 + t326;
t299 = sin(qJ(5));
t303 = cos(qJ(5));
t182 = -t299 * t193 + t303 * t198;
t260 = -t299 * qJD(3) + t303 * t342;
t220 = t260 * qJD(5) + t303 * qJDD(3) + t299 * t265;
t259 = qJDD(5) + t266;
t261 = t303 * qJD(3) + t299 * t342;
t280 = t283 + qJD(5);
t179 = (t260 * t280 - t220) * pkin(9) + (t260 * t261 + t259) * pkin(5) + t182;
t183 = t303 * t193 + t299 * t198;
t219 = -t261 * qJD(5) - t299 * qJDD(3) + t303 * t265;
t231 = t280 * pkin(5) - t261 * pkin(9);
t258 = t260 ^ 2;
t180 = -t258 * pkin(5) + t219 * pkin(9) - t280 * t231 + t183;
t298 = sin(qJ(6));
t302 = cos(qJ(6));
t177 = t302 * t179 - t298 * t180;
t221 = t302 * t260 - t298 * t261;
t191 = t221 * qJD(6) + t298 * t219 + t302 * t220;
t222 = t298 * t260 + t302 * t261;
t204 = -t221 * mrSges(7,1) + t222 * mrSges(7,2);
t277 = qJD(6) + t280;
t210 = -t277 * mrSges(7,2) + t221 * mrSges(7,3);
t248 = qJDD(6) + t259;
t173 = m(7) * t177 + t248 * mrSges(7,1) - t191 * mrSges(7,3) - t222 * t204 + t277 * t210;
t178 = t298 * t179 + t302 * t180;
t190 = -t222 * qJD(6) + t302 * t219 - t298 * t220;
t211 = t277 * mrSges(7,1) - t222 * mrSges(7,3);
t174 = m(7) * t178 - t248 * mrSges(7,2) + t190 * mrSges(7,3) + t221 * t204 - t277 * t211;
t163 = t302 * t173 + t298 * t174;
t224 = -t260 * mrSges(6,1) + t261 * mrSges(6,2);
t227 = -t280 * mrSges(6,2) + t260 * mrSges(6,3);
t160 = m(6) * t182 + t259 * mrSges(6,1) - t220 * mrSges(6,3) - t261 * t224 + t280 * t227 + t163;
t228 = t280 * mrSges(6,1) - t261 * mrSges(6,3);
t336 = -t298 * t173 + t302 * t174;
t161 = m(6) * t183 - t259 * mrSges(6,2) + t219 * mrSges(6,3) + t260 * t224 - t280 * t228 + t336;
t157 = t303 * t160 + t299 * t161;
t208 = -qJDD(3) * pkin(3) + t326 - t350;
t323 = -m(5) * t208 - t266 * mrSges(5,1) - t157;
t263 = (-mrSges(5,2) * t300 - mrSges(5,3) * t304) * qJD(1);
t334 = qJD(1) * (-t263 - (mrSges(4,1) * t300 + mrSges(4,2) * t304) * qJD(1));
t348 = mrSges(4,1) - mrSges(5,2);
t153 = m(4) * t225 - t266 * mrSges(4,3) + t348 * qJDD(3) + (t270 - t272) * qJD(3) + t304 * t334 + t323;
t226 = -t304 * g(3) + t300 * t233;
t271 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t283;
t206 = t306 * pkin(3) - qJDD(3) * qJ(4) + qJD(3) * t360 + t262 * t342 - t226;
t197 = -t265 * pkin(4) + qJD(3) * t274 - t295 * t351 - t206;
t185 = -t219 * pkin(5) - t258 * pkin(9) + t261 * t231 + t197;
t325 = m(7) * t185 - t190 * mrSges(7,1) + t191 * mrSges(7,2) - t221 * t210 + t222 * t211;
t175 = -m(6) * t197 + t219 * mrSges(6,1) - t220 * mrSges(6,2) + t260 * t227 - t261 * t228 - t325;
t273 = mrSges(5,1) * t283 + qJD(3) * mrSges(5,2);
t313 = -m(5) * t206 + qJDD(3) * mrSges(5,3) + qJD(3) * t273 - t175;
t167 = -qJDD(3) * mrSges(4,2) + t300 * t334 + t313 - qJD(3) * t271 + m(4) * t226 + (-mrSges(4,3) - mrSges(5,1)) * t265;
t145 = t304 * t153 + t300 * t167;
t238 = -qJDD(1) * pkin(1) + t327;
t242 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t304 - Ifges(4,2) * t300) * qJD(1);
t243 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t304 - Ifges(4,4) * t300) * qJD(1);
t199 = Ifges(7,5) * t222 + Ifges(7,6) * t221 + Ifges(7,3) * t277;
t201 = Ifges(7,1) * t222 + Ifges(7,4) * t221 + Ifges(7,5) * t277;
t164 = -mrSges(7,1) * t185 + mrSges(7,3) * t178 + Ifges(7,4) * t191 + Ifges(7,2) * t190 + Ifges(7,6) * t248 - t222 * t199 + t277 * t201;
t200 = Ifges(7,4) * t222 + Ifges(7,2) * t221 + Ifges(7,6) * t277;
t165 = mrSges(7,2) * t185 - mrSges(7,3) * t177 + Ifges(7,1) * t191 + Ifges(7,4) * t190 + Ifges(7,5) * t248 + t221 * t199 - t277 * t200;
t212 = Ifges(6,5) * t261 + Ifges(6,6) * t260 + Ifges(6,3) * t280;
t214 = Ifges(6,1) * t261 + Ifges(6,4) * t260 + Ifges(6,5) * t280;
t147 = -mrSges(6,1) * t197 + mrSges(6,3) * t183 + Ifges(6,4) * t220 + Ifges(6,2) * t219 + Ifges(6,6) * t259 - pkin(5) * t325 + pkin(9) * t336 + t302 * t164 + t298 * t165 - t261 * t212 + t280 * t214;
t213 = Ifges(6,4) * t261 + Ifges(6,2) * t260 + Ifges(6,6) * t280;
t149 = mrSges(6,2) * t197 - mrSges(6,3) * t182 + Ifges(6,1) * t220 + Ifges(6,4) * t219 + Ifges(6,5) * t259 - pkin(9) * t163 - t298 * t164 + t302 * t165 + t260 * t212 - t280 * t213;
t244 = Ifges(5,5) * qJD(3) + (-Ifges(5,6) * t304 + Ifges(5,3) * t300) * qJD(1);
t245 = Ifges(5,4) * qJD(3) + (-Ifges(5,2) * t304 + Ifges(5,6) * t300) * qJD(1);
t357 = mrSges(5,2) * t208 - mrSges(5,3) * t206 + Ifges(5,1) * qJDD(3) - Ifges(5,4) * t266 + Ifges(5,5) * t265 - pkin(8) * t157 - t299 * t147 + t303 * t149 - (t304 * t244 + t300 * t245) * qJD(1);
t308 = -mrSges(4,2) * t226 + pkin(3) * (-qJDD(3) * mrSges(5,2) - qJD(3) * t272 - t263 * t283 + t323) + qJ(4) * (-t265 * mrSges(5,1) - t263 * t342 + t313) + mrSges(4,1) * t225 + t243 * t342 + t242 * t283 - Ifges(4,6) * t265 + Ifges(4,5) * t266 + Ifges(4,3) * qJDD(3) + t357;
t359 = mrSges(3,1) * t238 + pkin(2) * t145 + t308;
t338 = t354 * t307;
t232 = t338 + t329;
t158 = -t299 * t160 + t303 * t161;
t205 = t265 * pkin(3) + t318 + t338;
t316 = m(5) * t205 - t266 * mrSges(5,3) - (t272 * t300 + t273 * t304) * qJD(1) + t158;
t358 = -m(4) * t232 - t266 * mrSges(4,2) - t348 * t265 - t270 * t342 - t271 * t283 - t316;
t349 = mrSges(2,1) - mrSges(3,2);
t347 = Ifges(2,5) - Ifges(3,4);
t346 = -Ifges(2,6) + Ifges(3,5);
t345 = -Ifges(5,6) - Ifges(4,4);
t146 = -t300 * t153 + t304 * t167;
t246 = Ifges(5,1) * qJD(3) + (-Ifges(5,4) * t304 + Ifges(5,5) * t300) * qJD(1);
t335 = qJD(1) * (-Ifges(4,3) * qJD(3) - (Ifges(4,5) * t304 - Ifges(4,6) * t300) * qJD(1) - t246);
t324 = -m(3) * t238 + t307 * mrSges(3,3) - t145;
t322 = -mrSges(7,1) * t177 + mrSges(7,2) * t178 - Ifges(7,5) * t191 - Ifges(7,6) * t190 - Ifges(7,3) * t248 - t222 * t200 + t221 * t201;
t154 = -t265 * mrSges(5,2) + t316;
t315 = -mrSges(5,1) * t206 + mrSges(5,2) * t205 - pkin(4) * t175 - pkin(8) * t158 - t303 * t147 - t299 * t149;
t139 = -mrSges(4,1) * t232 + mrSges(4,3) * t226 - pkin(3) * t154 - t345 * t266 + (-Ifges(4,2) - Ifges(5,3)) * t265 + (Ifges(4,6) - Ifges(5,5)) * qJDD(3) + (t243 - t245) * qJD(3) + t304 * t335 + t315;
t312 = -mrSges(6,1) * t182 + mrSges(6,2) * t183 - Ifges(6,5) * t220 - Ifges(6,6) * t219 - Ifges(6,3) * t259 - pkin(5) * t163 - t261 * t213 + t260 * t214 + t322;
t309 = -mrSges(5,1) * t208 + mrSges(5,3) * t205 - pkin(4) * t157 + t312;
t141 = (Ifges(4,1) + Ifges(5,2)) * t266 + t345 * t265 + (Ifges(4,5) - Ifges(5,4)) * qJDD(3) + (-t242 + t244) * qJD(3) - mrSges(4,3) * t225 + mrSges(4,2) * t232 - qJ(4) * t154 - t309 + t300 * t335;
t236 = t307 * pkin(1) - t329;
t320 = mrSges(3,2) * t238 - mrSges(3,3) * t236 + Ifges(3,1) * qJDD(1) - pkin(7) * t145 - t300 * t139 + t304 * t141;
t319 = -mrSges(3,1) * t236 - pkin(2) * t358 - pkin(7) * t146 - t304 * t139 - t300 * t141;
t310 = -m(3) * t236 + t307 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t358;
t314 = -mrSges(2,2) * t276 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t324) + qJ(2) * t310 + mrSges(2,1) * t275 + Ifges(2,3) * qJDD(1) + t320;
t150 = m(2) * t276 - t307 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t310;
t144 = -m(3) * g(3) + t146;
t142 = m(2) * t275 - t307 * mrSges(2,2) + t349 * qJDD(1) + t324;
t138 = t346 * t307 + t347 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t275 - qJ(2) * t144 + t359;
t137 = mrSges(2,3) * t276 - pkin(1) * t144 + t349 * g(3) - t346 * qJDD(1) + t347 * t307 + t319;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t305 * t138 - t301 * t137 - pkin(6) * (t305 * t142 + t301 * t150), t138, t320, t141, t357, t149, t165; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t301 * t138 + t305 * t137 + pkin(6) * (-t301 * t142 + t305 * t150), t137, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - t307 * Ifges(3,5) - t359, t139, Ifges(5,4) * qJDD(3) - Ifges(5,2) * t266 + Ifges(5,6) * t265 - qJD(3) * t244 + t246 * t342 + t309, t147, t164; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t314, t314, mrSges(3,2) * g(3) + t307 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t319, t308, Ifges(5,5) * qJDD(3) - Ifges(5,6) * t266 + Ifges(5,3) * t265 + qJD(3) * t245 + t246 * t283 - t315, -t312, -t322;];
m_new  = t1;
