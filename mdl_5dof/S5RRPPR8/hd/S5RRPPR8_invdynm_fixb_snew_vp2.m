% Calculate vector of cutting torques with Newton-Euler for
% S5RRPPR8
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x6]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPPR8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR8_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR8_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR8_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR8_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR8_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:38:13
% EndTime: 2019-12-31 19:38:22
% DurationCPUTime: 4.97s
% Computational Cost: add. (53576->306), mult. (121257->383), div. (0->0), fcn. (71608->8), ass. (0->114)
t297 = sin(qJ(1));
t300 = cos(qJ(1));
t272 = -t300 * g(1) - t297 * g(2);
t302 = qJD(1) ^ 2;
t245 = -t302 * pkin(1) + qJDD(1) * pkin(6) + t272;
t296 = sin(qJ(2));
t299 = cos(qJ(2));
t221 = -t299 * g(3) - t296 * t245;
t222 = -t296 * g(3) + t299 * t245;
t238 = Ifges(4,6) * qJD(2) + (Ifges(4,5) * t296 - Ifges(4,3) * t299) * qJD(1);
t241 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t296 + Ifges(3,2) * t299) * qJD(1);
t259 = (-mrSges(4,1) * t299 - mrSges(4,3) * t296) * qJD(1);
t319 = qJD(1) * qJD(2);
t317 = t299 * t319;
t261 = t296 * qJDD(1) + t317;
t318 = t296 * t319;
t262 = t299 * qJDD(1) - t318;
t258 = (-pkin(2) * t299 - qJ(3) * t296) * qJD(1);
t301 = qJD(2) ^ 2;
t320 = qJD(1) * t299;
t326 = 2 * qJD(3);
t202 = -t301 * pkin(2) + qJDD(2) * qJ(3) + qJD(2) * t326 + t258 * t320 + t222;
t321 = qJD(1) * t296;
t266 = -qJD(2) * pkin(3) - qJ(4) * t321;
t324 = t299 ^ 2 * t302;
t195 = -pkin(3) * t324 - t262 * qJ(4) + qJD(2) * t266 + t202;
t204 = -qJDD(2) * pkin(2) - t301 * qJ(3) + t258 * t321 + qJDD(3) - t221;
t196 = (-t261 + t317) * qJ(4) + (-t296 * t299 * t302 - qJDD(2)) * pkin(3) + t204;
t292 = sin(pkin(8));
t293 = cos(pkin(8));
t237 = (-t292 * t299 + t293 * t296) * qJD(1);
t172 = -0.2e1 * qJD(4) * t237 - t292 * t195 + t293 * t196;
t220 = t293 * t261 - t292 * t262;
t236 = (-t292 * t296 - t293 * t299) * qJD(1);
t168 = (-qJD(2) * t236 - t220) * pkin(7) + (t236 * t237 - qJDD(2)) * pkin(4) + t172;
t173 = 0.2e1 * qJD(4) * t236 + t293 * t195 + t292 * t196;
t219 = -t292 * t261 - t293 * t262;
t225 = -qJD(2) * pkin(4) - t237 * pkin(7);
t235 = t236 ^ 2;
t169 = -t235 * pkin(4) + t219 * pkin(7) + qJD(2) * t225 + t173;
t295 = sin(qJ(5));
t298 = cos(qJ(5));
t166 = t298 * t168 - t295 * t169;
t210 = t298 * t236 - t295 * t237;
t183 = t210 * qJD(5) + t295 * t219 + t298 * t220;
t211 = t295 * t236 + t298 * t237;
t194 = -t210 * mrSges(6,1) + t211 * mrSges(6,2);
t282 = -qJD(2) + qJD(5);
t205 = -t282 * mrSges(6,2) + t210 * mrSges(6,3);
t281 = -qJDD(2) + qJDD(5);
t161 = m(6) * t166 + t281 * mrSges(6,1) - t183 * mrSges(6,3) - t211 * t194 + t282 * t205;
t167 = t295 * t168 + t298 * t169;
t182 = -t211 * qJD(5) + t298 * t219 - t295 * t220;
t206 = t282 * mrSges(6,1) - t211 * mrSges(6,3);
t162 = m(6) * t167 - t281 * mrSges(6,2) + t182 * mrSges(6,3) + t210 * t194 - t282 * t206;
t152 = t298 * t161 + t295 * t162;
t214 = -t236 * mrSges(5,1) + t237 * mrSges(5,2);
t223 = qJD(2) * mrSges(5,2) + t236 * mrSges(5,3);
t149 = m(5) * t172 - qJDD(2) * mrSges(5,1) - t220 * mrSges(5,3) - qJD(2) * t223 - t237 * t214 + t152;
t224 = -qJD(2) * mrSges(5,1) - t237 * mrSges(5,3);
t315 = -t295 * t161 + t298 * t162;
t150 = m(5) * t173 + qJDD(2) * mrSges(5,2) + t219 * mrSges(5,3) + qJD(2) * t224 + t236 * t214 + t315;
t146 = t293 * t149 + t292 * t150;
t208 = Ifges(5,4) * t237 + Ifges(5,2) * t236 - Ifges(5,6) * qJD(2);
t209 = Ifges(5,1) * t237 + Ifges(5,4) * t236 - Ifges(5,5) * qJD(2);
t187 = Ifges(6,4) * t211 + Ifges(6,2) * t210 + Ifges(6,6) * t282;
t188 = Ifges(6,1) * t211 + Ifges(6,4) * t210 + Ifges(6,5) * t282;
t313 = mrSges(6,1) * t166 - mrSges(6,2) * t167 + Ifges(6,5) * t183 + Ifges(6,6) * t182 + Ifges(6,3) * t281 + t211 * t187 - t210 * t188;
t306 = mrSges(5,1) * t172 - mrSges(5,2) * t173 + Ifges(5,5) * t220 + Ifges(5,6) * t219 - Ifges(5,3) * qJDD(2) + pkin(4) * t152 + t237 * t208 - t236 * t209 + t313;
t305 = -mrSges(4,1) * t204 + mrSges(4,3) * t202 + Ifges(4,4) * t261 + Ifges(4,2) * qJDD(2) - Ifges(4,6) * t262 - pkin(3) * t146 - t306;
t270 = mrSges(4,2) * t320 + qJD(2) * mrSges(4,3);
t309 = -m(4) * t204 + qJDD(2) * mrSges(4,1) + qJD(2) * t270 - t146;
t147 = -t292 * t149 + t293 * t150;
t268 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t321;
t311 = m(4) * t202 + qJDD(2) * mrSges(4,3) + qJD(2) * t268 + t259 * t320 + t147;
t242 = Ifges(4,4) * qJD(2) + (Ifges(4,1) * t296 - Ifges(4,5) * t299) * qJD(1);
t322 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t296 + Ifges(3,4) * t299) * qJD(1) + t242;
t328 = -((t238 - t241) * t296 + t299 * t322) * qJD(1) + mrSges(3,1) * t221 - mrSges(3,2) * t222 + Ifges(3,5) * t261 + Ifges(3,6) * t262 + Ifges(3,3) * qJDD(2) + pkin(2) * (-t261 * mrSges(4,2) - t259 * t321 + t309) + qJ(3) * (t262 * mrSges(4,2) + t311) + t305;
t325 = mrSges(3,3) + mrSges(4,2);
t271 = t297 * g(1) - t300 * g(2);
t260 = (-mrSges(3,1) * t299 + mrSges(3,2) * t296) * qJD(1);
t267 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t321;
t142 = m(3) * t222 - qJDD(2) * mrSges(3,2) - qJD(2) * t267 + t260 * t320 + t325 * t262 + t311;
t269 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t320;
t143 = m(3) * t221 + qJDD(2) * mrSges(3,1) + qJD(2) * t269 - t325 * t261 + (-t259 - t260) * t321 + t309;
t316 = t299 * t142 - t296 * t143;
t244 = -qJDD(1) * pkin(1) - t302 * pkin(6) - t271;
t312 = -t262 * pkin(2) + t244 + (-t261 - t317) * qJ(3);
t185 = -pkin(2) * t318 + t262 * pkin(3) - qJ(4) * t324 + qJDD(4) - t312 + (t266 + t326) * t321;
t175 = -t219 * pkin(4) - t235 * pkin(7) + t237 * t225 + t185;
t314 = m(6) * t175 - t182 * mrSges(6,1) + t183 * mrSges(6,2) - t210 * t205 + t211 * t206;
t163 = -m(5) * t185 + t219 * mrSges(5,1) - t220 * mrSges(5,2) + t236 * t223 - t237 * t224 - t314;
t197 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t321 + t312;
t157 = m(4) * t197 - t262 * mrSges(4,1) - t261 * mrSges(4,3) - t268 * t321 - t270 * t320 + t163;
t239 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t296 + Ifges(3,6) * t299) * qJD(1);
t240 = Ifges(4,2) * qJD(2) + (Ifges(4,4) * t296 - Ifges(4,6) * t299) * qJD(1);
t186 = Ifges(6,5) * t211 + Ifges(6,6) * t210 + Ifges(6,3) * t282;
t153 = -mrSges(6,1) * t175 + mrSges(6,3) * t167 + Ifges(6,4) * t183 + Ifges(6,2) * t182 + Ifges(6,6) * t281 - t211 * t186 + t282 * t188;
t154 = mrSges(6,2) * t175 - mrSges(6,3) * t166 + Ifges(6,1) * t183 + Ifges(6,4) * t182 + Ifges(6,5) * t281 + t210 * t186 - t282 * t187;
t207 = Ifges(5,5) * t237 + Ifges(5,6) * t236 - Ifges(5,3) * qJD(2);
t138 = -mrSges(5,1) * t185 + mrSges(5,3) * t173 + Ifges(5,4) * t220 + Ifges(5,2) * t219 - Ifges(5,6) * qJDD(2) - pkin(4) * t314 + pkin(7) * t315 - qJD(2) * t209 + t298 * t153 + t295 * t154 - t237 * t207;
t140 = mrSges(5,2) * t185 - mrSges(5,3) * t172 + Ifges(5,1) * t220 + Ifges(5,4) * t219 - Ifges(5,5) * qJDD(2) - pkin(7) * t152 + qJD(2) * t208 - t295 * t153 + t298 * t154 + t236 * t207;
t307 = -mrSges(4,1) * t197 + mrSges(4,2) * t202 - pkin(3) * t163 - qJ(4) * t147 - t293 * t138 - t292 * t140;
t132 = -mrSges(3,1) * t244 + mrSges(3,3) * t222 - pkin(2) * t157 + (Ifges(3,2) + Ifges(4,3)) * t262 + (Ifges(3,4) - Ifges(4,5)) * t261 + (Ifges(3,6) - Ifges(4,6)) * qJDD(2) + t322 * qJD(2) + (-t239 - t240) * t321 + t307;
t308 = mrSges(4,2) * t204 - mrSges(4,3) * t197 + Ifges(4,1) * t261 + Ifges(4,4) * qJDD(2) - Ifges(4,5) * t262 - qJ(4) * t146 + qJD(2) * t238 - t292 * t138 + t293 * t140 + t240 * t320;
t134 = mrSges(3,2) * t244 - mrSges(3,3) * t221 + Ifges(3,1) * t261 + Ifges(3,4) * t262 + Ifges(3,5) * qJDD(2) - qJ(3) * t157 - qJD(2) * t241 + t239 * t320 + t308;
t304 = -m(3) * t244 + t262 * mrSges(3,1) - t261 * mrSges(3,2) - t267 * t321 + t269 * t320 - t157;
t310 = mrSges(2,1) * t271 - mrSges(2,2) * t272 + Ifges(2,3) * qJDD(1) + pkin(1) * t304 + pkin(6) * t316 + t299 * t132 + t296 * t134;
t155 = m(2) * t271 + qJDD(1) * mrSges(2,1) - t302 * mrSges(2,2) + t304;
t137 = t296 * t142 + t299 * t143;
t135 = m(2) * t272 - t302 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t316;
t130 = mrSges(2,1) * g(3) + mrSges(2,3) * t272 + t302 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t137 - t328;
t129 = -mrSges(2,2) * g(3) - mrSges(2,3) * t271 + Ifges(2,5) * qJDD(1) - t302 * Ifges(2,6) - pkin(6) * t137 - t296 * t132 + t299 * t134;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t300 * t129 - t297 * t130 - pkin(5) * (t297 * t135 + t300 * t155), t129, t134, t308, t140, t154; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t297 * t129 + t300 * t130 + pkin(5) * (t300 * t135 - t297 * t155), t130, t132, t305 + (-t296 * t238 - t299 * t242) * qJD(1), t138, t153; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t310, t310, t328, Ifges(4,5) * t261 + Ifges(4,6) * qJDD(2) - Ifges(4,3) * t262 - qJD(2) * t242 + t240 * t321 - t307, t306, t313;];
m_new = t1;
