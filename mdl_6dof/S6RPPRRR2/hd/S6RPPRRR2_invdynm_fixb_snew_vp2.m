% Calculate vector of cutting torques with Newton-Euler for
% S6RPPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Datum: 2019-05-05 15:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPPRRR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:20:27
% EndTime: 2019-05-05 15:20:50
% DurationCPUTime: 17.50s
% Computational Cost: add. (318475->322), mult. (712740->403), div. (0->0), fcn. (511323->12), ass. (0->142)
t308 = qJD(1) ^ 2;
t295 = sin(pkin(11));
t336 = qJD(1) * t295;
t297 = cos(pkin(11));
t335 = qJD(1) * t297;
t302 = sin(qJ(1));
t306 = cos(qJ(1));
t277 = t302 * g(1) - t306 * g(2);
t274 = qJDD(1) * pkin(1) + t277;
t278 = -t306 * g(1) - t302 * g(2);
t275 = -t308 * pkin(1) + t278;
t296 = sin(pkin(10));
t298 = cos(pkin(10));
t256 = t296 * t274 + t298 * t275;
t243 = -t308 * pkin(2) + qJDD(1) * qJ(3) + t256;
t294 = -g(3) + qJDD(2);
t333 = qJD(1) * qJD(3);
t337 = t297 * t294 - 0.2e1 * t295 * t333;
t340 = pkin(3) * t297;
t222 = (-pkin(7) * qJDD(1) + t308 * t340 - t243) * t295 + t337;
t230 = t295 * t294 + (t243 + 0.2e1 * t333) * t297;
t331 = qJDD(1) * t297;
t290 = t297 ^ 2;
t338 = t290 * t308;
t225 = -pkin(3) * t338 + pkin(7) * t331 + t230;
t301 = sin(qJ(4));
t305 = cos(qJ(4));
t206 = t301 * t222 + t305 * t225;
t265 = -t301 * t336 + t305 * t335;
t320 = t295 * t305 + t297 * t301;
t266 = t320 * qJD(1);
t246 = -t265 * mrSges(5,1) + t266 * mrSges(5,2);
t263 = t266 * qJD(4);
t332 = qJDD(1) * t295;
t252 = -t301 * t332 + t305 * t331 - t263;
t261 = qJD(4) * mrSges(5,1) - t266 * mrSges(5,3);
t251 = -t265 * pkin(4) - t266 * pkin(8);
t307 = qJD(4) ^ 2;
t199 = -t307 * pkin(4) + qJDD(4) * pkin(8) + t265 * t251 + t206;
t289 = t295 ^ 2;
t255 = t298 * t274 - t296 * t275;
t322 = qJDD(3) - t255;
t226 = (-pkin(2) - t340) * qJDD(1) + (-qJ(3) + (-t289 - t290) * pkin(7)) * t308 + t322;
t334 = t265 * qJD(4);
t253 = t320 * qJDD(1) + t334;
t203 = (-t253 - t334) * pkin(8) + (-t252 + t263) * pkin(4) + t226;
t300 = sin(qJ(5));
t304 = cos(qJ(5));
t189 = -t300 * t199 + t304 * t203;
t258 = t304 * qJD(4) - t300 * t266;
t224 = t258 * qJD(5) + t300 * qJDD(4) + t304 * t253;
t250 = qJDD(5) - t252;
t259 = t300 * qJD(4) + t304 * t266;
t264 = qJD(5) - t265;
t187 = (t258 * t264 - t224) * pkin(9) + (t258 * t259 + t250) * pkin(5) + t189;
t190 = t304 * t199 + t300 * t203;
t223 = -t259 * qJD(5) + t304 * qJDD(4) - t300 * t253;
t236 = t264 * pkin(5) - t259 * pkin(9);
t257 = t258 ^ 2;
t188 = -t257 * pkin(5) + t223 * pkin(9) - t264 * t236 + t190;
t299 = sin(qJ(6));
t303 = cos(qJ(6));
t185 = t303 * t187 - t299 * t188;
t227 = t303 * t258 - t299 * t259;
t197 = t227 * qJD(6) + t299 * t223 + t303 * t224;
t228 = t299 * t258 + t303 * t259;
t211 = -t227 * mrSges(7,1) + t228 * mrSges(7,2);
t262 = qJD(6) + t264;
t212 = -t262 * mrSges(7,2) + t227 * mrSges(7,3);
t245 = qJDD(6) + t250;
t180 = m(7) * t185 + t245 * mrSges(7,1) - t197 * mrSges(7,3) - t228 * t211 + t262 * t212;
t186 = t299 * t187 + t303 * t188;
t196 = -t228 * qJD(6) + t303 * t223 - t299 * t224;
t213 = t262 * mrSges(7,1) - t228 * mrSges(7,3);
t181 = m(7) * t186 - t245 * mrSges(7,2) + t196 * mrSges(7,3) + t227 * t211 - t262 * t213;
t172 = t303 * t180 + t299 * t181;
t232 = -t258 * mrSges(6,1) + t259 * mrSges(6,2);
t234 = -t264 * mrSges(6,2) + t258 * mrSges(6,3);
t170 = m(6) * t189 + t250 * mrSges(6,1) - t224 * mrSges(6,3) - t259 * t232 + t264 * t234 + t172;
t235 = t264 * mrSges(6,1) - t259 * mrSges(6,3);
t326 = -t299 * t180 + t303 * t181;
t171 = m(6) * t190 - t250 * mrSges(6,2) + t223 * mrSges(6,3) + t258 * t232 - t264 * t235 + t326;
t327 = -t300 * t170 + t304 * t171;
t163 = m(5) * t206 - qJDD(4) * mrSges(5,2) + t252 * mrSges(5,3) - qJD(4) * t261 + t265 * t246 + t327;
t205 = t305 * t222 - t301 * t225;
t260 = -qJD(4) * mrSges(5,2) + t265 * mrSges(5,3);
t198 = -qJDD(4) * pkin(4) - t307 * pkin(8) + t266 * t251 - t205;
t191 = -t223 * pkin(5) - t257 * pkin(9) + t259 * t236 + t198;
t317 = m(7) * t191 - t196 * mrSges(7,1) + t197 * mrSges(7,2) - t227 * t212 + t228 * t213;
t311 = -m(6) * t198 + t223 * mrSges(6,1) - t224 * mrSges(6,2) + t258 * t234 - t259 * t235 - t317;
t176 = m(5) * t205 + qJDD(4) * mrSges(5,1) - t253 * mrSges(5,3) + qJD(4) * t260 - t266 * t246 + t311;
t153 = t301 * t163 + t305 * t176;
t229 = -t295 * t243 + t337;
t207 = Ifges(7,5) * t228 + Ifges(7,6) * t227 + Ifges(7,3) * t262;
t209 = Ifges(7,1) * t228 + Ifges(7,4) * t227 + Ifges(7,5) * t262;
t173 = -mrSges(7,1) * t191 + mrSges(7,3) * t186 + Ifges(7,4) * t197 + Ifges(7,2) * t196 + Ifges(7,6) * t245 - t228 * t207 + t262 * t209;
t208 = Ifges(7,4) * t228 + Ifges(7,2) * t227 + Ifges(7,6) * t262;
t174 = mrSges(7,2) * t191 - mrSges(7,3) * t185 + Ifges(7,1) * t197 + Ifges(7,4) * t196 + Ifges(7,5) * t245 + t227 * t207 - t262 * t208;
t215 = Ifges(6,5) * t259 + Ifges(6,6) * t258 + Ifges(6,3) * t264;
t217 = Ifges(6,1) * t259 + Ifges(6,4) * t258 + Ifges(6,5) * t264;
t155 = -mrSges(6,1) * t198 + mrSges(6,3) * t190 + Ifges(6,4) * t224 + Ifges(6,2) * t223 + Ifges(6,6) * t250 - pkin(5) * t317 + pkin(9) * t326 + t303 * t173 + t299 * t174 - t259 * t215 + t264 * t217;
t216 = Ifges(6,4) * t259 + Ifges(6,2) * t258 + Ifges(6,6) * t264;
t157 = mrSges(6,2) * t198 - mrSges(6,3) * t189 + Ifges(6,1) * t224 + Ifges(6,4) * t223 + Ifges(6,5) * t250 - pkin(9) * t172 - t299 * t173 + t303 * t174 + t258 * t215 - t264 * t216;
t240 = Ifges(5,4) * t266 + Ifges(5,2) * t265 + Ifges(5,6) * qJD(4);
t241 = Ifges(5,1) * t266 + Ifges(5,4) * t265 + Ifges(5,5) * qJD(4);
t313 = -mrSges(5,1) * t205 + mrSges(5,2) * t206 - Ifges(5,5) * t253 - Ifges(5,6) * t252 - Ifges(5,3) * qJDD(4) - pkin(4) * t311 - pkin(8) * t327 - t304 * t155 - t300 * t157 - t266 * t240 + t265 * t241;
t324 = Ifges(4,4) * t295 + Ifges(4,2) * t297;
t325 = Ifges(4,1) * t295 + Ifges(4,4) * t297;
t341 = -mrSges(4,1) * t229 + mrSges(4,2) * t230 - pkin(3) * t153 - (t324 * t336 - t325 * t335) * qJD(1) + t313;
t339 = mrSges(4,2) * t295;
t319 = mrSges(4,3) * qJDD(1) + t308 * (-mrSges(4,1) * t297 + t339);
t151 = m(4) * t229 - t319 * t295 + t153;
t328 = t305 * t163 - t301 * t176;
t152 = m(4) * t230 + t319 * t297 + t328;
t329 = -t295 * t151 + t297 * t152;
t144 = m(3) * t256 - t308 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t329;
t238 = -qJDD(1) * pkin(2) - t308 * qJ(3) + t322;
t165 = t304 * t170 + t300 * t171;
t315 = m(5) * t226 - t252 * mrSges(5,1) + t253 * mrSges(5,2) - t265 * t260 + t266 * t261 + t165;
t312 = -m(4) * t238 + mrSges(4,1) * t331 - t315 + (t289 * t308 + t338) * mrSges(4,3);
t159 = t312 + (mrSges(3,1) - t339) * qJDD(1) - t308 * mrSges(3,2) + m(3) * t255;
t140 = t296 * t144 + t298 * t159;
t146 = t297 * t151 + t295 * t152;
t330 = t298 * t144 - t296 * t159;
t323 = Ifges(4,5) * t295 + Ifges(4,6) * t297;
t239 = Ifges(5,5) * t266 + Ifges(5,6) * t265 + Ifges(5,3) * qJD(4);
t141 = mrSges(5,2) * t226 - mrSges(5,3) * t205 + Ifges(5,1) * t253 + Ifges(5,4) * t252 + Ifges(5,5) * qJDD(4) - pkin(8) * t165 - qJD(4) * t240 - t300 * t155 + t304 * t157 + t265 * t239;
t316 = -mrSges(7,1) * t185 + mrSges(7,2) * t186 - Ifges(7,5) * t197 - Ifges(7,6) * t196 - Ifges(7,3) * t245 - t228 * t208 + t227 * t209;
t309 = mrSges(6,1) * t189 - mrSges(6,2) * t190 + Ifges(6,5) * t224 + Ifges(6,6) * t223 + Ifges(6,3) * t250 + pkin(5) * t172 + t259 * t216 - t258 * t217 - t316;
t147 = -mrSges(5,1) * t226 + mrSges(5,3) * t206 + Ifges(5,4) * t253 + Ifges(5,2) * t252 + Ifges(5,6) * qJDD(4) - pkin(4) * t165 + qJD(4) * t241 - t266 * t239 - t309;
t271 = t323 * qJD(1);
t133 = -mrSges(4,1) * t238 + mrSges(4,3) * t230 - pkin(3) * t315 + pkin(7) * t328 + t324 * qJDD(1) + t301 * t141 + t305 * t147 - t271 * t336;
t136 = mrSges(4,2) * t238 - mrSges(4,3) * t229 - pkin(7) * t153 + t325 * qJDD(1) + t305 * t141 - t301 * t147 + t271 * t335;
t318 = -mrSges(3,2) * t256 + qJ(3) * t329 + t297 * t133 + t295 * t136 + pkin(2) * (-mrSges(4,2) * t332 + t312) + mrSges(3,1) * t255 + Ifges(3,3) * qJDD(1);
t314 = mrSges(2,1) * t277 - mrSges(2,2) * t278 + Ifges(2,3) * qJDD(1) + pkin(1) * t140 + t318;
t138 = m(2) * t278 - t308 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t330;
t137 = m(2) * t277 + qJDD(1) * mrSges(2,1) - t308 * mrSges(2,2) + t140;
t134 = -pkin(2) * t146 + (Ifges(3,6) - t323) * qJDD(1) + t308 * Ifges(3,5) - mrSges(3,1) * t294 + mrSges(3,3) * t256 + t341;
t131 = mrSges(3,2) * t294 - mrSges(3,3) * t255 + Ifges(3,5) * qJDD(1) - t308 * Ifges(3,6) - qJ(3) * t146 - t295 * t133 + t297 * t136;
t130 = -mrSges(2,2) * g(3) - mrSges(2,3) * t277 + Ifges(2,5) * qJDD(1) - t308 * Ifges(2,6) - qJ(2) * t140 + t298 * t131 - t296 * t134;
t129 = Ifges(2,6) * qJDD(1) + t308 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t278 + t296 * t131 + t298 * t134 - pkin(1) * (m(3) * t294 + t146) + qJ(2) * t330;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t306 * t130 - t302 * t129 - pkin(6) * (t306 * t137 + t302 * t138), t130, t131, t136, t141, t157, t174; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t302 * t130 + t306 * t129 + pkin(6) * (-t302 * t137 + t306 * t138), t129, t134, t133, t147, t155, t173; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t314, t314, t318, t323 * qJDD(1) - t341, -t313, t309, -t316;];
m_new  = t1;
