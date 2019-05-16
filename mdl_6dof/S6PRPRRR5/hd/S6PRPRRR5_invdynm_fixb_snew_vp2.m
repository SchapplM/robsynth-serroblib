% Calculate vector of cutting torques with Newton-Euler for
% S6PRPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 01:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRPRRR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 01:15:43
% EndTime: 2019-05-05 01:15:59
% DurationCPUTime: 9.36s
% Computational Cost: add. (172974->299), mult. (327496->375), div. (0->0), fcn. (222750->12), ass. (0->131)
t273 = sin(pkin(11));
t275 = cos(pkin(11));
t254 = g(1) * t273 - g(2) * t275;
t255 = -g(1) * t275 - g(2) * t273;
t270 = -g(3) + qJDD(1);
t274 = sin(pkin(6));
t276 = cos(pkin(6));
t280 = sin(qJ(2));
t284 = cos(qJ(2));
t221 = -t280 * t255 + (t254 * t276 + t270 * t274) * t284;
t285 = qJD(2) ^ 2;
t292 = -t285 * qJ(3) + qJDD(3) - t221;
t317 = -pkin(2) - pkin(8);
t206 = t317 * qJDD(2) + t292;
t283 = cos(qJ(4));
t203 = t283 * t206;
t232 = -t254 * t274 + t270 * t276;
t279 = sin(qJ(4));
t305 = qJD(2) * qJD(4);
t252 = qJDD(2) * t283 - t279 * t305;
t316 = pkin(4) * t285;
t186 = qJDD(4) * pkin(4) - t252 * pkin(9) + t203 + (-pkin(9) * t305 - t283 * t316 - t232) * t279;
t195 = t279 * t206 + t283 * t232;
t251 = -qJDD(2) * t279 - t283 * t305;
t306 = qJD(2) * t283;
t259 = qJD(4) * pkin(4) - pkin(9) * t306;
t269 = t279 ^ 2;
t187 = pkin(9) * t251 - qJD(4) * t259 - t269 * t316 + t195;
t278 = sin(qJ(5));
t282 = cos(qJ(5));
t182 = t278 * t186 + t282 * t187;
t241 = (-t278 * t279 + t282 * t283) * qJD(2);
t213 = -qJD(5) * t241 + t251 * t282 - t252 * t278;
t240 = (t278 * t283 + t279 * t282) * qJD(2);
t224 = mrSges(6,1) * t240 + mrSges(6,2) * t241;
t264 = qJD(4) + qJD(5);
t230 = mrSges(6,1) * t264 - mrSges(6,3) * t241;
t263 = qJDD(4) + qJDD(5);
t225 = pkin(5) * t240 - pkin(10) * t241;
t262 = t264 ^ 2;
t179 = -pkin(5) * t262 + pkin(10) * t263 - t225 * t240 + t182;
t308 = t276 * t280;
t309 = t274 * t280;
t222 = t254 * t308 + t284 * t255 + t270 * t309;
t300 = qJDD(2) * qJ(3) + 0.2e1 * qJD(3) * qJD(2) + t222;
t192 = -t251 * pkin(4) + t259 * t306 + (-pkin(9) * t269 + t317) * t285 + t300;
t214 = -qJD(5) * t240 + t251 * t278 + t252 * t282;
t183 = (t240 * t264 - t214) * pkin(10) + (t241 * t264 - t213) * pkin(5) + t192;
t277 = sin(qJ(6));
t281 = cos(qJ(6));
t176 = -t179 * t277 + t183 * t281;
t226 = -t241 * t277 + t264 * t281;
t190 = qJD(6) * t226 + t214 * t281 + t263 * t277;
t227 = t241 * t281 + t264 * t277;
t200 = -mrSges(7,1) * t226 + mrSges(7,2) * t227;
t210 = qJDD(6) - t213;
t235 = qJD(6) + t240;
t216 = -mrSges(7,2) * t235 + mrSges(7,3) * t226;
t172 = m(7) * t176 + mrSges(7,1) * t210 - t190 * mrSges(7,3) - t200 * t227 + t216 * t235;
t177 = t179 * t281 + t183 * t277;
t189 = -qJD(6) * t227 - t214 * t277 + t263 * t281;
t217 = mrSges(7,1) * t235 - mrSges(7,3) * t227;
t173 = m(7) * t177 - mrSges(7,2) * t210 + t189 * mrSges(7,3) + t200 * t226 - t217 * t235;
t302 = -t172 * t277 + t281 * t173;
t159 = m(6) * t182 - mrSges(6,2) * t263 + mrSges(6,3) * t213 - t224 * t240 - t230 * t264 + t302;
t181 = t186 * t282 - t187 * t278;
t229 = -mrSges(6,2) * t264 - mrSges(6,3) * t240;
t178 = -pkin(5) * t263 - pkin(10) * t262 + t225 * t241 - t181;
t295 = -m(7) * t178 + t189 * mrSges(7,1) - t190 * mrSges(7,2) + t226 * t216 - t217 * t227;
t168 = m(6) * t181 + mrSges(6,1) * t263 - mrSges(6,3) * t214 - t224 * t241 + t229 * t264 + t295;
t152 = t278 * t159 + t282 * t168;
t194 = -t279 * t232 + t203;
t250 = (mrSges(5,1) * t279 + mrSges(5,2) * t283) * qJD(2);
t307 = qJD(2) * t279;
t256 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t307;
t149 = m(5) * t194 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t252 + qJD(4) * t256 - t250 * t306 + t152;
t257 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t306;
t303 = t282 * t159 - t168 * t278;
t150 = m(5) * t195 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t251 - qJD(4) * t257 - t250 * t307 + t303;
t145 = -t149 * t279 + t283 * t150;
t143 = m(4) * t232 + t145;
t161 = t281 * t172 + t277 * t173;
t196 = Ifges(7,5) * t227 + Ifges(7,6) * t226 + Ifges(7,3) * t235;
t198 = Ifges(7,1) * t227 + Ifges(7,4) * t226 + Ifges(7,5) * t235;
t165 = -mrSges(7,1) * t178 + mrSges(7,3) * t177 + Ifges(7,4) * t190 + Ifges(7,2) * t189 + Ifges(7,6) * t210 - t196 * t227 + t198 * t235;
t197 = Ifges(7,4) * t227 + Ifges(7,2) * t226 + Ifges(7,6) * t235;
t166 = mrSges(7,2) * t178 - mrSges(7,3) * t176 + Ifges(7,1) * t190 + Ifges(7,4) * t189 + Ifges(7,5) * t210 + t196 * t226 - t197 * t235;
t218 = Ifges(6,5) * t241 - Ifges(6,6) * t240 + Ifges(6,3) * t264;
t219 = Ifges(6,4) * t241 - Ifges(6,2) * t240 + Ifges(6,6) * t264;
t146 = mrSges(6,2) * t192 - mrSges(6,3) * t181 + Ifges(6,1) * t214 + Ifges(6,4) * t213 + Ifges(6,5) * t263 - pkin(10) * t161 - t165 * t277 + t166 * t281 - t218 * t240 - t219 * t264;
t220 = Ifges(6,1) * t241 - Ifges(6,4) * t240 + Ifges(6,5) * t264;
t289 = mrSges(7,1) * t176 - mrSges(7,2) * t177 + Ifges(7,5) * t190 + Ifges(7,6) * t189 + Ifges(7,3) * t210 + t197 * t227 - t198 * t226;
t147 = -mrSges(6,1) * t192 + mrSges(6,3) * t182 + Ifges(6,4) * t214 + Ifges(6,2) * t213 + Ifges(6,6) * t263 - pkin(5) * t161 - t218 * t241 + t220 * t264 - t289;
t205 = t317 * t285 + t300;
t237 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t283 - Ifges(5,6) * t279) * qJD(2);
t239 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t283 - Ifges(5,4) * t279) * qJD(2);
t294 = m(6) * t192 - mrSges(6,1) * t213 + t214 * mrSges(6,2) + t229 * t240 + t241 * t230 + t161;
t134 = -mrSges(5,1) * t205 + mrSges(5,3) * t195 + Ifges(5,4) * t252 + Ifges(5,2) * t251 + Ifges(5,6) * qJDD(4) - pkin(4) * t294 + pkin(9) * t303 + qJD(4) * t239 + t278 * t146 + t282 * t147 - t237 * t306;
t238 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t283 - Ifges(5,2) * t279) * qJD(2);
t136 = mrSges(5,2) * t205 - mrSges(5,3) * t194 + Ifges(5,1) * t252 + Ifges(5,4) * t251 + Ifges(5,5) * qJDD(4) - pkin(9) * t152 - qJD(4) * t238 + t146 * t282 - t147 * t278 - t237 * t307;
t156 = -m(5) * t205 + mrSges(5,1) * t251 - t252 * mrSges(5,2) - t256 * t307 - t257 * t306 - t294;
t211 = t285 * pkin(2) - t300;
t291 = -mrSges(4,1) * t211 - pkin(3) * t156 - pkin(8) * t145 - t283 * t134 - t279 * t136;
t312 = -Ifges(3,6) + Ifges(4,5);
t313 = Ifges(3,5) - Ifges(4,4);
t314 = mrSges(3,1) - mrSges(4,2);
t127 = mrSges(3,3) * t222 - pkin(2) * t143 - t312 * qJDD(2) - t314 * t232 + t313 * t285 + t291;
t144 = t283 * t149 + t279 * t150;
t215 = -qJDD(2) * pkin(2) + t292;
t297 = -m(4) * t215 + t285 * mrSges(4,3) - t144;
t141 = m(3) * t221 - t285 * mrSges(3,2) + t314 * qJDD(2) + t297;
t288 = -m(4) * t211 + t285 * mrSges(4,2) + qJDD(2) * mrSges(4,3) - t156;
t155 = m(3) * t222 - mrSges(3,1) * t285 - qJDD(2) * mrSges(3,2) + t288;
t139 = -t141 * t280 + t284 * t155;
t318 = pkin(7) * t139 + t127 * t284;
t310 = t141 * t284;
t142 = m(3) * t232 + t143;
t133 = -t142 * t274 + t155 * t308 + t276 * t310;
t293 = mrSges(4,2) * t215 - mrSges(4,3) * t211 + Ifges(4,1) * qJDD(2) - pkin(8) * t144 - t279 * t134 + t283 * t136;
t125 = Ifges(3,3) * qJDD(2) + mrSges(3,1) * t221 - mrSges(3,2) * t222 + pkin(2) * (-qJDD(2) * mrSges(4,2) + t297) + qJ(3) * t288 + t293;
t290 = -mrSges(6,1) * t181 + mrSges(6,2) * t182 - Ifges(6,5) * t214 - Ifges(6,6) * t213 - Ifges(6,3) * t263 - pkin(5) * t295 - pkin(10) * t302 - t281 * t165 - t277 * t166 - t241 * t219 - t240 * t220;
t287 = -mrSges(5,1) * t194 + mrSges(5,2) * t195 - Ifges(5,5) * t252 - Ifges(5,6) * t251 - Ifges(5,3) * qJDD(4) - pkin(4) * t152 - t238 * t306 - t239 * t307 + t290;
t286 = -mrSges(4,1) * t215 - pkin(3) * t144 + t287;
t129 = -t286 + t312 * t285 + (mrSges(3,2) - mrSges(4,3)) * t232 + t313 * qJDD(2) - mrSges(3,3) * t221 - qJ(3) * t143;
t296 = mrSges(2,1) * t254 - mrSges(2,2) * t255 + pkin(1) * t133 + t276 * t125 + t129 * t309 + t318 * t274;
t137 = m(2) * t255 + t139;
t132 = t276 * t142 + (t155 * t280 + t310) * t274;
t130 = m(2) * t254 + t133;
t123 = mrSges(2,2) * t270 - mrSges(2,3) * t254 - t280 * t127 + t284 * t129 + (-t132 * t274 - t133 * t276) * pkin(7);
t122 = -mrSges(2,1) * t270 + mrSges(2,3) * t255 - pkin(1) * t132 - t274 * t125 + (t129 * t280 + t318) * t276;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t275 * t123 - t273 * t122 - qJ(1) * (t130 * t275 + t137 * t273), t123, t129, t293, t136, t146, t166; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t273 * t123 + t275 * t122 + qJ(1) * (-t130 * t273 + t137 * t275), t122, t127, mrSges(4,3) * t232 + Ifges(4,4) * qJDD(2) - t285 * Ifges(4,5) + t286, t134, t147, t165; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t296, t296, t125, -mrSges(4,2) * t232 + t285 * Ifges(4,4) + Ifges(4,5) * qJDD(2) - t291, -t287, -t290, t289;];
m_new  = t1;
