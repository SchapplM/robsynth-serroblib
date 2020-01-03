% Calculate vector of cutting torques with Newton-Euler for
% S5RRPPR6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPPR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR6_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR6_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR6_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR6_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:31:45
% EndTime: 2019-12-31 19:32:01
% DurationCPUTime: 8.85s
% Computational Cost: add. (135344->311), mult. (318634->398), div. (0->0), fcn. (217597->10), ass. (0->123)
t298 = -2 * qJD(3);
t267 = sin(qJ(2));
t270 = cos(qJ(2));
t290 = qJD(1) * qJD(2);
t250 = t267 * qJDD(1) + t270 * t290;
t268 = sin(qJ(1));
t271 = cos(qJ(1));
t257 = -t271 * g(1) - t268 * g(2);
t273 = qJD(1) ^ 2;
t245 = -t273 * pkin(1) + qJDD(1) * pkin(6) + t257;
t294 = t267 * t245;
t296 = pkin(2) * t273;
t201 = qJDD(2) * pkin(2) - t250 * qJ(3) - t294 + (qJ(3) * t290 + t267 * t296 - g(3)) * t270;
t231 = -t267 * g(3) + t270 * t245;
t251 = t270 * qJDD(1) - t267 * t290;
t293 = qJD(1) * t267;
t253 = qJD(2) * pkin(2) - qJ(3) * t293;
t262 = t270 ^ 2;
t203 = t251 * qJ(3) - qJD(2) * t253 - t262 * t296 + t231;
t264 = sin(pkin(8));
t295 = cos(pkin(8));
t239 = (t264 * t270 + t295 * t267) * qJD(1);
t186 = t295 * t201 - t264 * t203 + t239 * t298;
t292 = qJD(1) * t270;
t238 = t264 * t293 - t295 * t292;
t187 = t264 * t201 + t295 * t203 + t238 * t298;
t217 = t238 * mrSges(4,1) + t239 * mrSges(4,2);
t223 = t264 * t250 - t295 * t251;
t233 = qJD(2) * mrSges(4,1) - t239 * mrSges(4,3);
t216 = t238 * pkin(3) - t239 * qJ(4);
t272 = qJD(2) ^ 2;
t172 = -t272 * pkin(3) + qJDD(2) * qJ(4) - t238 * t216 + t187;
t256 = t268 * g(1) - t271 * g(2);
t283 = -qJDD(1) * pkin(1) - t256;
t205 = -t251 * pkin(2) + qJDD(3) + t253 * t293 + (-qJ(3) * t262 - pkin(6)) * t273 + t283;
t224 = t295 * t250 + t264 * t251;
t175 = (qJD(2) * t238 - t224) * qJ(4) + (qJD(2) * t239 + t223) * pkin(3) + t205;
t263 = sin(pkin(9));
t265 = cos(pkin(9));
t229 = t263 * qJD(2) + t265 * t239;
t166 = -0.2e1 * qJD(4) * t229 - t263 * t172 + t265 * t175;
t211 = t263 * qJDD(2) + t265 * t224;
t228 = t265 * qJD(2) - t263 * t239;
t164 = (t228 * t238 - t211) * pkin(7) + (t228 * t229 + t223) * pkin(4) + t166;
t167 = 0.2e1 * qJD(4) * t228 + t265 * t172 + t263 * t175;
t208 = t238 * pkin(4) - t229 * pkin(7);
t210 = t265 * qJDD(2) - t263 * t224;
t227 = t228 ^ 2;
t165 = -t227 * pkin(4) + t210 * pkin(7) - t238 * t208 + t167;
t266 = sin(qJ(5));
t269 = cos(qJ(5));
t162 = t269 * t164 - t266 * t165;
t196 = t269 * t228 - t266 * t229;
t180 = t196 * qJD(5) + t266 * t210 + t269 * t211;
t197 = t266 * t228 + t269 * t229;
t188 = -t196 * mrSges(6,1) + t197 * mrSges(6,2);
t237 = qJD(5) + t238;
t189 = -t237 * mrSges(6,2) + t196 * mrSges(6,3);
t222 = qJDD(5) + t223;
t157 = m(6) * t162 + t222 * mrSges(6,1) - t180 * mrSges(6,3) - t197 * t188 + t237 * t189;
t163 = t266 * t164 + t269 * t165;
t179 = -t197 * qJD(5) + t269 * t210 - t266 * t211;
t190 = t237 * mrSges(6,1) - t197 * mrSges(6,3);
t158 = m(6) * t163 - t222 * mrSges(6,2) + t179 * mrSges(6,3) + t196 * t188 - t237 * t190;
t149 = t269 * t157 + t266 * t158;
t202 = -t228 * mrSges(5,1) + t229 * mrSges(5,2);
t206 = -t238 * mrSges(5,2) + t228 * mrSges(5,3);
t147 = m(5) * t166 + t223 * mrSges(5,1) - t211 * mrSges(5,3) - t229 * t202 + t238 * t206 + t149;
t207 = t238 * mrSges(5,1) - t229 * mrSges(5,3);
t286 = -t266 * t157 + t269 * t158;
t148 = m(5) * t167 - t223 * mrSges(5,2) + t210 * mrSges(5,3) + t228 * t202 - t238 * t207 + t286;
t287 = -t263 * t147 + t265 * t148;
t140 = m(4) * t187 - qJDD(2) * mrSges(4,2) - t223 * mrSges(4,3) - qJD(2) * t233 - t238 * t217 + t287;
t232 = -qJD(2) * mrSges(4,2) - t238 * mrSges(4,3);
t171 = -qJDD(2) * pkin(3) - t272 * qJ(4) + t239 * t216 + qJDD(4) - t186;
t168 = -t210 * pkin(4) - t227 * pkin(7) + t229 * t208 + t171;
t281 = m(6) * t168 - t179 * mrSges(6,1) + t180 * mrSges(6,2) - t196 * t189 + t197 * t190;
t277 = -m(5) * t171 + t210 * mrSges(5,1) - t211 * mrSges(5,2) + t228 * t206 - t229 * t207 - t281;
t153 = m(4) * t186 + qJDD(2) * mrSges(4,1) - t224 * mrSges(4,3) + qJD(2) * t232 - t239 * t217 + t277;
t131 = t264 * t140 + t295 * t153;
t230 = -t270 * g(3) - t294;
t241 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t267 + Ifges(3,2) * t270) * qJD(1);
t242 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t267 + Ifges(3,4) * t270) * qJD(1);
t181 = Ifges(6,5) * t197 + Ifges(6,6) * t196 + Ifges(6,3) * t237;
t183 = Ifges(6,1) * t197 + Ifges(6,4) * t196 + Ifges(6,5) * t237;
t150 = -mrSges(6,1) * t168 + mrSges(6,3) * t163 + Ifges(6,4) * t180 + Ifges(6,2) * t179 + Ifges(6,6) * t222 - t197 * t181 + t237 * t183;
t182 = Ifges(6,4) * t197 + Ifges(6,2) * t196 + Ifges(6,6) * t237;
t151 = mrSges(6,2) * t168 - mrSges(6,3) * t162 + Ifges(6,1) * t180 + Ifges(6,4) * t179 + Ifges(6,5) * t222 + t196 * t181 - t237 * t182;
t191 = Ifges(5,5) * t229 + Ifges(5,6) * t228 + Ifges(5,3) * t238;
t193 = Ifges(5,1) * t229 + Ifges(5,4) * t228 + Ifges(5,5) * t238;
t133 = -mrSges(5,1) * t171 + mrSges(5,3) * t167 + Ifges(5,4) * t211 + Ifges(5,2) * t210 + Ifges(5,6) * t223 - pkin(4) * t281 + pkin(7) * t286 + t269 * t150 + t266 * t151 - t229 * t191 + t238 * t193;
t192 = Ifges(5,4) * t229 + Ifges(5,2) * t228 + Ifges(5,6) * t238;
t135 = mrSges(5,2) * t171 - mrSges(5,3) * t166 + Ifges(5,1) * t211 + Ifges(5,4) * t210 + Ifges(5,5) * t223 - pkin(7) * t149 - t266 * t150 + t269 * t151 + t228 * t191 - t238 * t192;
t213 = Ifges(4,4) * t239 - Ifges(4,2) * t238 + Ifges(4,6) * qJD(2);
t214 = Ifges(4,1) * t239 - Ifges(4,4) * t238 + Ifges(4,5) * qJD(2);
t278 = -mrSges(4,1) * t186 + mrSges(4,2) * t187 - Ifges(4,5) * t224 + Ifges(4,6) * t223 - Ifges(4,3) * qJDD(2) - pkin(3) * t277 - qJ(4) * t287 - t265 * t133 - t263 * t135 - t239 * t213 - t238 * t214;
t297 = mrSges(3,1) * t230 - mrSges(3,2) * t231 + Ifges(3,5) * t250 + Ifges(3,6) * t251 + Ifges(3,3) * qJDD(2) + pkin(2) * t131 + (t267 * t241 - t270 * t242) * qJD(1) - t278;
t142 = t265 * t147 + t263 * t148;
t249 = (-mrSges(3,1) * t270 + mrSges(3,2) * t267) * qJD(1);
t255 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t292;
t129 = m(3) * t230 + qJDD(2) * mrSges(3,1) - t250 * mrSges(3,3) + qJD(2) * t255 - t249 * t293 + t131;
t254 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t293;
t288 = t295 * t140 - t264 * t153;
t130 = m(3) * t231 - qJDD(2) * mrSges(3,2) + t251 * mrSges(3,3) - qJD(2) * t254 + t249 * t292 + t288;
t289 = -t267 * t129 + t270 * t130;
t212 = Ifges(4,5) * t239 - Ifges(4,6) * t238 + Ifges(4,3) * qJD(2);
t123 = mrSges(4,2) * t205 - mrSges(4,3) * t186 + Ifges(4,1) * t224 - Ifges(4,4) * t223 + Ifges(4,5) * qJDD(2) - qJ(4) * t142 - qJD(2) * t213 - t263 * t133 + t265 * t135 - t238 * t212;
t280 = -mrSges(6,1) * t162 + mrSges(6,2) * t163 - Ifges(6,5) * t180 - Ifges(6,6) * t179 - Ifges(6,3) * t222 - t197 * t182 + t196 * t183;
t275 = -mrSges(5,1) * t166 + mrSges(5,2) * t167 - Ifges(5,5) * t211 - Ifges(5,6) * t210 - pkin(4) * t149 - t229 * t192 + t228 * t193 + t280;
t127 = t275 + Ifges(4,6) * qJDD(2) + (-Ifges(4,2) - Ifges(5,3)) * t223 - t239 * t212 + qJD(2) * t214 + Ifges(4,4) * t224 - mrSges(4,1) * t205 + mrSges(4,3) * t187 - pkin(3) * t142;
t240 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t267 + Ifges(3,6) * t270) * qJD(1);
t244 = -t273 * pkin(6) + t283;
t279 = m(4) * t205 + t223 * mrSges(4,1) + t224 * mrSges(4,2) + t238 * t232 + t239 * t233 + t142;
t119 = -mrSges(3,1) * t244 + mrSges(3,3) * t231 + Ifges(3,4) * t250 + Ifges(3,2) * t251 + Ifges(3,6) * qJDD(2) - pkin(2) * t279 + qJ(3) * t288 + qJD(2) * t242 + t264 * t123 + t295 * t127 - t240 * t293;
t122 = mrSges(3,2) * t244 - mrSges(3,3) * t230 + Ifges(3,1) * t250 + Ifges(3,4) * t251 + Ifges(3,5) * qJDD(2) - qJ(3) * t131 - qJD(2) * t241 + t295 * t123 - t264 * t127 + t240 * t292;
t276 = -m(3) * t244 + t251 * mrSges(3,1) - t250 * mrSges(3,2) - t254 * t293 + t255 * t292 - t279;
t282 = mrSges(2,1) * t256 - mrSges(2,2) * t257 + Ifges(2,3) * qJDD(1) + pkin(1) * t276 + pkin(6) * t289 + t270 * t119 + t267 * t122;
t136 = m(2) * t256 + qJDD(1) * mrSges(2,1) - t273 * mrSges(2,2) + t276;
t126 = t270 * t129 + t267 * t130;
t124 = m(2) * t257 - t273 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t289;
t120 = mrSges(2,1) * g(3) + mrSges(2,3) * t257 + t273 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t126 - t297;
t117 = -mrSges(2,2) * g(3) - mrSges(2,3) * t256 + Ifges(2,5) * qJDD(1) - t273 * Ifges(2,6) - pkin(6) * t126 - t267 * t119 + t270 * t122;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t271 * t117 - t268 * t120 - pkin(5) * (t268 * t124 + t271 * t136), t117, t122, t123, t135, t151; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t268 * t117 + t271 * t120 + pkin(5) * (t271 * t124 - t268 * t136), t120, t119, t127, t133, t150; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t282, t282, t297, -t278, Ifges(5,3) * t223 - t275, -t280;];
m_new = t1;
