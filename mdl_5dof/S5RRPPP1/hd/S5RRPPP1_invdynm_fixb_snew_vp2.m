% Calculate vector of cutting torques with Newton-Euler for
% S5RRPPP1
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
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
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
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPPP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPP1_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPP1_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPP1_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPP1_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPP1_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:23:27
% EndTime: 2019-12-31 19:23:40
% DurationCPUTime: 6.44s
% Computational Cost: add. (68354->338), mult. (171896->400), div. (0->0), fcn. (119601->8), ass. (0->131)
t330 = -2 * qJD(3);
t275 = cos(qJ(2));
t320 = cos(pkin(5));
t299 = qJD(1) * t320;
t272 = sin(pkin(5));
t309 = qJD(2) * t272;
t244 = (t275 * t299 + t309) * qJ(3);
t274 = sin(qJ(1));
t276 = cos(qJ(1));
t264 = t274 * g(1) - t276 * g(2);
t277 = qJD(1) ^ 2;
t248 = -qJDD(1) * pkin(1) - t277 * pkin(7) - t264;
t273 = sin(qJ(2));
t318 = qJ(3) * t273;
t256 = qJD(2) * pkin(2) - t299 * t318;
t307 = qJD(1) * qJD(2);
t258 = t273 * qJDD(1) + t275 * t307;
t259 = t275 * qJDD(1) - t273 * t307;
t184 = -t258 * t272 * qJ(3) - t259 * pkin(2) + (-t244 * t275 + t256 * t273) * qJD(1) + t248;
t265 = -t276 * g(1) - t274 * g(2);
t249 = -t277 * pkin(1) + qJDD(1) * pkin(7) + t265;
t250 = (-pkin(2) * t275 - t272 * t318) * qJD(1);
t303 = qJ(3) * t320;
t325 = t275 * g(3);
t185 = -t258 * t303 + qJDD(2) * pkin(2) - t325 + qJD(2) * t244 + (-qJD(1) * t250 - t249) * t273;
t233 = -t273 * g(3) + t275 * t249;
t293 = qJDD(2) * t272 + t320 * t259;
t310 = qJD(1) * t275;
t186 = t293 * qJ(3) - qJD(2) * t256 + t250 * t310 + t233;
t271 = sin(pkin(8));
t302 = t271 * t320;
t319 = cos(pkin(8));
t231 = t271 * t309 + (t319 * t273 + t275 * t302) * qJD(1);
t296 = t320 * t319;
t301 = t272 * t319;
t170 = t184 * t301 + t185 * t296 - t271 * t186 + t231 * t330;
t223 = -qJDD(2) * t301 + t271 * t258 - t259 * t296;
t172 = t320 * t184 - t272 * t185 + qJDD(3);
t224 = t319 * t258 + t293 * t271;
t253 = -t320 * qJD(2) + t272 * t310;
t311 = qJD(1) * t273;
t230 = -qJD(2) * t301 + t271 * t311 - t296 * t310;
t317 = t230 * t253;
t327 = -2 * qJD(4);
t283 = (-t224 - t317) * qJ(4) + t172 + (-t253 * pkin(3) + t327) * t231;
t169 = t223 * pkin(3) + t283;
t211 = t231 * mrSges(5,1) - t253 * mrSges(5,2);
t329 = m(5) * t169 - t224 * mrSges(5,3) - t231 * t211;
t201 = -t231 * mrSges(6,2) + t230 * mrSges(6,3);
t203 = t230 * mrSges(4,1) + t231 * mrSges(4,2);
t241 = t320 * qJDD(2) - t272 * t259;
t202 = t230 * pkin(3) - t231 * qJ(4);
t252 = t253 ^ 2;
t167 = -t241 * pkin(3) - t252 * qJ(4) + t231 * t202 + qJDD(4) - t170;
t204 = -t230 * mrSges(5,2) - t231 * mrSges(5,3);
t326 = 2 * qJD(5);
t158 = t253 * t326 + (t230 * t231 - t241) * qJ(5) + (t224 - t317) * pkin(4) + t167;
t210 = -t230 * mrSges(6,1) - t253 * mrSges(6,2);
t297 = -m(6) * t158 + t241 * mrSges(6,3) - t253 * t210;
t286 = -m(5) * t167 - t224 * mrSges(5,1) - t231 * t204 + t297;
t209 = t230 * mrSges(5,1) + t253 * mrSges(5,3);
t312 = t253 * mrSges(4,2) - t230 * mrSges(4,3) - t209;
t323 = -mrSges(6,1) - mrSges(4,3);
t324 = mrSges(4,1) - mrSges(5,2);
t148 = m(4) * t170 - t312 * t253 + t324 * t241 + (-t201 - t203) * t231 + t323 * t224 + t286;
t227 = t230 * t330;
t305 = t272 * t271 * t184 + t185 * t302 + t319 * t186;
t171 = t227 + t305;
t206 = -t253 * mrSges(4,1) - t231 * mrSges(4,3);
t288 = t252 * pkin(3) - t241 * qJ(4) - t305;
t165 = 0.2e1 * qJD(4) * t253 + ((2 * qJD(3)) + t202) * t230 + t288;
t207 = t231 * pkin(4) + t253 * qJ(5);
t229 = t230 ^ 2;
t161 = -t223 * pkin(4) - t229 * qJ(5) - t230 * t202 + qJDD(5) + t227 + (t327 - t207) * t253 - t288;
t208 = t231 * mrSges(6,1) + t253 * mrSges(6,3);
t304 = -m(6) * t161 - t241 * mrSges(6,2) + t253 * t208;
t292 = -m(5) * t165 + t241 * mrSges(5,3) - t253 * t211 - t304;
t313 = -t201 - t204;
t151 = m(4) * t171 - t241 * mrSges(4,2) + t253 * t206 + (-t203 + t313) * t230 + (-mrSges(5,1) + t323) * t223 + t292;
t164 = -t229 * pkin(4) + t230 * t326 - t231 * t207 + (pkin(3) + qJ(5)) * t223 + t283;
t306 = m(6) * t164 + t223 * mrSges(6,3) + t230 * t210;
t152 = m(4) * t172 + (t206 - t208) * t231 + t312 * t230 + (mrSges(4,2) - mrSges(6,2)) * t224 + t324 * t223 + t306 + t329;
t141 = (t319 * t148 + t151 * t271) * t272 + t320 * t152;
t190 = -Ifges(6,4) * t253 + Ifges(6,2) * t230 + Ifges(6,6) * t231;
t192 = Ifges(4,4) * t231 - Ifges(4,2) * t230 - Ifges(4,6) * t253;
t321 = t224 * mrSges(6,1);
t154 = t231 * t201 - t297 + t321;
t188 = -Ifges(5,5) * t253 - Ifges(5,6) * t231 + Ifges(5,3) * t230;
t187 = -Ifges(6,5) * t253 + Ifges(6,6) * t230 + Ifges(6,3) * t231;
t289 = mrSges(6,2) * t161 - mrSges(6,3) * t158 + Ifges(6,1) * t241 + Ifges(6,4) * t223 + Ifges(6,5) * t224 + t230 * t187;
t280 = mrSges(5,2) * t167 - mrSges(5,3) * t165 + Ifges(5,1) * t241 - Ifges(5,4) * t224 + Ifges(5,5) * t223 - qJ(5) * t154 - t231 * t188 + t289;
t191 = -Ifges(5,4) * t253 - Ifges(5,2) * t231 + Ifges(5,6) * t230;
t315 = -Ifges(4,1) * t231 + Ifges(4,4) * t230 + Ifges(4,5) * t253 + t191;
t136 = t280 + qJ(4) * t292 + pkin(3) * (-t241 * mrSges(5,2) + t253 * t209 + t286 - t321) + Ifges(4,3) * t241 + Ifges(4,5) * t224 + mrSges(4,1) * t170 - mrSges(4,2) * t171 + (-pkin(3) * t201 - t190 + t192) * t231 + (qJ(4) * t313 - t315) * t230 + (-Ifges(4,6) + qJ(4) * (-mrSges(5,1) - mrSges(6,1))) * t223;
t155 = -t224 * mrSges(6,2) - t231 * t208 + t306;
t153 = -t223 * mrSges(5,2) - t230 * t209 + t155 + t329;
t189 = Ifges(4,5) * t231 - Ifges(4,6) * t230 - Ifges(4,3) * t253;
t194 = -Ifges(5,1) * t253 - Ifges(5,4) * t231 + Ifges(5,5) * t230;
t193 = -Ifges(6,1) * t253 + Ifges(6,4) * t230 + Ifges(6,5) * t231;
t291 = mrSges(6,1) * t161 - mrSges(6,3) * t164 - Ifges(6,4) * t241 - Ifges(6,2) * t223 - Ifges(6,6) * t224 - t231 * t193;
t281 = mrSges(5,1) * t165 - mrSges(5,2) * t169 + pkin(4) * (t223 * mrSges(6,1) + t230 * t201 + t304) + qJ(5) * t155 - t291;
t322 = Ifges(4,4) + Ifges(5,6);
t137 = -t281 + mrSges(4,3) * t171 - mrSges(4,1) * t172 - pkin(3) * t153 + (-t187 + t315) * t253 + (Ifges(4,6) - Ifges(5,5)) * t241 + (-t189 - t194) * t231 + t322 * t224 + (-Ifges(4,2) - Ifges(5,3)) * t223;
t142 = t148 * t296 + t151 * t302 - t272 * t152;
t290 = -mrSges(6,1) * t158 + mrSges(6,2) * t164 - Ifges(6,5) * t241 - Ifges(6,6) * t223 - Ifges(6,3) * t224 + t253 * t190;
t282 = -mrSges(5,1) * t167 + mrSges(5,3) * t169 - pkin(4) * t154 + t290;
t314 = t193 + t194;
t143 = -t282 - mrSges(4,3) * t170 + mrSges(4,2) * t172 - qJ(4) * t153 - t322 * t223 + (Ifges(4,1) + Ifges(5,2)) * t224 + (-t189 - t314) * t230 + (Ifges(4,5) - Ifges(5,4)) * t241 + (t192 - t188) * t253;
t146 = -t271 * t148 + t319 * t151;
t232 = -t273 * t249 - t325;
t246 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t273 + Ifges(3,2) * t275) * qJD(1);
t247 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t273 + Ifges(3,4) * t275) * qJD(1);
t328 = mrSges(3,1) * t232 - mrSges(3,2) * t233 + Ifges(3,5) * t258 + Ifges(3,6) * t259 + Ifges(3,3) * qJDD(2) + pkin(2) * t142 + (t246 * t273 - t247 * t275) * qJD(1) + t320 * t136 + (qJ(3) * t146 + t319 * t137 + t143 * t271) * t272;
t316 = t231 * t190;
t257 = (-mrSges(3,1) * t275 + mrSges(3,2) * t273) * qJD(1);
t263 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t310;
t140 = m(3) * t232 + qJDD(2) * mrSges(3,1) - t258 * mrSges(3,3) + qJD(2) * t263 - t257 * t311 + t142;
t262 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t311;
t145 = m(3) * t233 - qJDD(2) * mrSges(3,2) + t259 * mrSges(3,3) - qJD(2) * t262 + t257 * t310 + t146;
t298 = -t273 * t140 + t275 * t145;
t245 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t273 + Ifges(3,6) * t275) * qJD(1);
t130 = -mrSges(3,1) * t248 + mrSges(3,3) * t233 + Ifges(3,4) * t258 + Ifges(3,2) * t259 + Ifges(3,6) * qJDD(2) - pkin(2) * t141 + qJD(2) * t247 - t272 * t136 + t137 * t296 + t143 * t302 + t146 * t303 - t245 * t311;
t132 = t245 * t310 + mrSges(3,2) * t248 - mrSges(3,3) * t232 + t319 * t143 + Ifges(3,1) * t258 + Ifges(3,4) * t259 + Ifges(3,5) * qJDD(2) - qJD(2) * t246 - t271 * t137 + (-t141 * t272 - t320 * t142) * qJ(3);
t278 = -m(3) * t248 + t259 * mrSges(3,1) - t258 * mrSges(3,2) - t262 * t311 + t263 * t310 - t141;
t287 = mrSges(2,1) * t264 - mrSges(2,2) * t265 + Ifges(2,3) * qJDD(1) + pkin(1) * t278 + pkin(7) * t298 + t275 * t130 + t273 * t132;
t138 = m(2) * t264 + qJDD(1) * mrSges(2,1) - t277 * mrSges(2,2) + t278;
t135 = t275 * t140 + t273 * t145;
t133 = m(2) * t265 - t277 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t298;
t128 = mrSges(2,1) * g(3) + mrSges(2,3) * t265 + t277 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t135 - t328;
t127 = -mrSges(2,2) * g(3) - mrSges(2,3) * t264 + Ifges(2,5) * qJDD(1) - t277 * Ifges(2,6) - pkin(7) * t135 - t273 * t130 + t275 * t132;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t276 * t127 - t274 * t128 - pkin(6) * (t274 * t133 + t276 * t138), t127, t132, t143, -t230 * t191 + t280 - t316, t289 - t316; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t274 * t127 + t276 * t128 + pkin(6) * (t276 * t133 - t274 * t138), t128, t130, t137, Ifges(5,4) * t241 - Ifges(5,2) * t224 + Ifges(5,6) * t223 + t253 * t188 + t314 * t230 + t282, t253 * t187 - t291; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t287, t287, t328, t136, t281 + (t187 - t191) * t253 + Ifges(5,5) * t241 + t231 * t194 + Ifges(5,3) * t223 - Ifges(5,6) * t224, -t230 * t193 - t290;];
m_new = t1;
