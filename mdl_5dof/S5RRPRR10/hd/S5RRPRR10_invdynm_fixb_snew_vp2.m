% Calculate vector of cutting torques with Newton-Euler for
% S5RRPRR10
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPRR10_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR10_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR10_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR10_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR10_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR10_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR10_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:24:06
% EndTime: 2019-12-31 20:24:34
% DurationCPUTime: 16.62s
% Computational Cost: add. (257378->321), mult. (674936->427), div. (0->0), fcn. (517798->12), ass. (0->136)
t308 = -2 * qJD(3);
t272 = sin(pkin(10));
t274 = cos(pkin(10));
t278 = sin(qJ(2));
t282 = cos(qJ(2));
t273 = sin(pkin(5));
t301 = qJD(1) * t273;
t250 = (t272 * t278 - t274 * t282) * t301;
t299 = qJD(1) * qJD(2);
t259 = (qJDD(1) * t278 + t282 * t299) * t273;
t275 = cos(pkin(5));
t266 = t275 * qJDD(1) + qJDD(2);
t267 = t275 * qJD(1) + qJD(2);
t279 = sin(qJ(1));
t283 = cos(qJ(1));
t263 = t279 * g(1) - t283 * g(2);
t284 = qJD(1) ^ 2;
t307 = pkin(7) * t273;
t256 = qJDD(1) * pkin(1) + t284 * t307 + t263;
t264 = -t283 * g(1) - t279 * g(2);
t257 = -t284 * pkin(1) + qJDD(1) * t307 + t264;
t302 = t275 * t282;
t292 = t256 * t302 - t278 * t257;
t306 = t273 ^ 2 * t284;
t198 = t266 * pkin(2) - t259 * qJ(3) + (pkin(2) * t278 * t306 + (qJ(3) * qJD(1) * t267 - g(3)) * t273) * t282 + t292;
t303 = t275 * t278;
t305 = t273 * t278;
t225 = -g(3) * t305 + t256 * t303 + t282 * t257;
t297 = t278 * t301;
t253 = t267 * pkin(2) - qJ(3) * t297;
t260 = (qJDD(1) * t282 - t278 * t299) * t273;
t298 = t282 ^ 2 * t306;
t201 = -pkin(2) * t298 + t260 * qJ(3) - t267 * t253 + t225;
t251 = (t272 * t282 + t274 * t278) * t301;
t183 = t274 * t198 - t272 * t201 + t251 * t308;
t304 = t273 * t282;
t184 = t272 * t198 + t274 * t201 + t250 * t308;
t226 = t250 * mrSges(4,1) + t251 * mrSges(4,2);
t232 = -t272 * t259 + t274 * t260;
t238 = t267 * mrSges(4,1) - t251 * mrSges(4,3);
t227 = t250 * pkin(3) - t251 * pkin(8);
t265 = t267 ^ 2;
t181 = -t265 * pkin(3) + t266 * pkin(8) - t250 * t227 + t184;
t242 = -t275 * g(3) - t273 * t256;
t211 = -t260 * pkin(2) - qJ(3) * t298 + t253 * t297 + qJDD(3) + t242;
t233 = t274 * t259 + t272 * t260;
t186 = (t250 * t267 - t233) * pkin(8) + (t251 * t267 - t232) * pkin(3) + t211;
t277 = sin(qJ(4));
t281 = cos(qJ(4));
t178 = t281 * t181 + t277 * t186;
t235 = -t277 * t251 + t281 * t267;
t236 = t281 * t251 + t277 * t267;
t214 = -t235 * pkin(4) - t236 * pkin(9);
t231 = qJDD(4) - t232;
t249 = qJD(4) + t250;
t248 = t249 ^ 2;
t175 = -t248 * pkin(4) + t231 * pkin(9) + t235 * t214 + t178;
t180 = -t266 * pkin(3) - t265 * pkin(8) + t251 * t227 - t183;
t208 = -t236 * qJD(4) - t277 * t233 + t281 * t266;
t209 = t235 * qJD(4) + t281 * t233 + t277 * t266;
t176 = (-t235 * t249 - t209) * pkin(9) + (t236 * t249 - t208) * pkin(4) + t180;
t276 = sin(qJ(5));
t280 = cos(qJ(5));
t172 = -t276 * t175 + t280 * t176;
t216 = -t276 * t236 + t280 * t249;
t189 = t216 * qJD(5) + t280 * t209 + t276 * t231;
t217 = t280 * t236 + t276 * t249;
t194 = -t216 * mrSges(6,1) + t217 * mrSges(6,2);
t234 = qJD(5) - t235;
t199 = -t234 * mrSges(6,2) + t216 * mrSges(6,3);
t207 = qJDD(5) - t208;
t169 = m(6) * t172 + t207 * mrSges(6,1) - t189 * mrSges(6,3) - t217 * t194 + t234 * t199;
t173 = t280 * t175 + t276 * t176;
t188 = -t217 * qJD(5) - t276 * t209 + t280 * t231;
t200 = t234 * mrSges(6,1) - t217 * mrSges(6,3);
t170 = m(6) * t173 - t207 * mrSges(6,2) + t188 * mrSges(6,3) + t216 * t194 - t234 * t200;
t163 = -t276 * t169 + t280 * t170;
t213 = -t235 * mrSges(5,1) + t236 * mrSges(5,2);
t219 = t249 * mrSges(5,1) - t236 * mrSges(5,3);
t160 = m(5) * t178 - t231 * mrSges(5,2) + t208 * mrSges(5,3) + t235 * t213 - t249 * t219 + t163;
t177 = -t277 * t181 + t281 * t186;
t174 = -t231 * pkin(4) - t248 * pkin(9) + t236 * t214 - t177;
t171 = -m(6) * t174 + t188 * mrSges(6,1) - t189 * mrSges(6,2) + t216 * t199 - t217 * t200;
t218 = -t249 * mrSges(5,2) + t235 * mrSges(5,3);
t167 = m(5) * t177 + t231 * mrSges(5,1) - t209 * mrSges(5,3) - t236 * t213 + t249 * t218 + t171;
t294 = t281 * t160 - t277 * t167;
t152 = m(4) * t184 - t266 * mrSges(4,2) + t232 * mrSges(4,3) - t250 * t226 - t267 * t238 + t294;
t237 = -t267 * mrSges(4,2) - t250 * mrSges(4,3);
t162 = t280 * t169 + t276 * t170;
t287 = -m(5) * t180 + t208 * mrSges(5,1) - t209 * mrSges(5,2) + t235 * t218 - t236 * t219 - t162;
t157 = m(4) * t183 + t266 * mrSges(4,1) - t233 * mrSges(4,3) - t251 * t226 + t267 * t237 + t287;
t145 = t272 * t152 + t274 * t157;
t155 = t277 * t160 + t281 * t167;
t296 = t282 * t301;
t224 = -g(3) * t304 + t292;
t255 = -t267 * mrSges(3,2) + mrSges(3,3) * t296;
t258 = (-mrSges(3,1) * t282 + mrSges(3,2) * t278) * t301;
t143 = m(3) * t224 + t266 * mrSges(3,1) - t259 * mrSges(3,3) + t267 * t255 - t258 * t297 + t145;
t254 = t267 * mrSges(3,1) - mrSges(3,3) * t297;
t295 = t274 * t152 - t272 * t157;
t144 = m(3) * t225 - t266 * mrSges(3,2) + t260 * mrSges(3,3) - t267 * t254 + t258 * t296 + t295;
t138 = -t278 * t143 + t282 * t144;
t289 = m(4) * t211 - t232 * mrSges(4,1) + t233 * mrSges(4,2) + t250 * t237 + t251 * t238 + t155;
t153 = m(3) * t242 - t260 * mrSges(3,1) + t259 * mrSges(3,2) + (t254 * t278 - t255 * t282) * t301 + t289;
t134 = t143 * t302 + t144 * t303 - t273 * t153;
t190 = Ifges(6,5) * t217 + Ifges(6,6) * t216 + Ifges(6,3) * t234;
t192 = Ifges(6,1) * t217 + Ifges(6,4) * t216 + Ifges(6,5) * t234;
t164 = -mrSges(6,1) * t174 + mrSges(6,3) * t173 + Ifges(6,4) * t189 + Ifges(6,2) * t188 + Ifges(6,6) * t207 - t217 * t190 + t234 * t192;
t191 = Ifges(6,4) * t217 + Ifges(6,2) * t216 + Ifges(6,6) * t234;
t165 = mrSges(6,2) * t174 - mrSges(6,3) * t172 + Ifges(6,1) * t189 + Ifges(6,4) * t188 + Ifges(6,5) * t207 + t216 * t190 - t234 * t191;
t202 = Ifges(5,5) * t236 + Ifges(5,6) * t235 + Ifges(5,3) * t249;
t203 = Ifges(5,4) * t236 + Ifges(5,2) * t235 + Ifges(5,6) * t249;
t147 = mrSges(5,2) * t180 - mrSges(5,3) * t177 + Ifges(5,1) * t209 + Ifges(5,4) * t208 + Ifges(5,5) * t231 - pkin(9) * t162 - t276 * t164 + t280 * t165 + t235 * t202 - t249 * t203;
t204 = Ifges(5,1) * t236 + Ifges(5,4) * t235 + Ifges(5,5) * t249;
t286 = mrSges(6,1) * t172 - mrSges(6,2) * t173 + Ifges(6,5) * t189 + Ifges(6,6) * t188 + Ifges(6,3) * t207 + t217 * t191 - t216 * t192;
t149 = -mrSges(5,1) * t180 + mrSges(5,3) * t178 + Ifges(5,4) * t209 + Ifges(5,2) * t208 + Ifges(5,6) * t231 - pkin(4) * t162 - t236 * t202 + t249 * t204 - t286;
t220 = Ifges(4,5) * t251 - Ifges(4,6) * t250 + Ifges(4,3) * t267;
t221 = Ifges(4,4) * t251 - Ifges(4,2) * t250 + Ifges(4,6) * t267;
t135 = mrSges(4,2) * t211 - mrSges(4,3) * t183 + Ifges(4,1) * t233 + Ifges(4,4) * t232 + Ifges(4,5) * t266 - pkin(8) * t155 + t281 * t147 - t277 * t149 - t250 * t220 - t267 * t221;
t222 = Ifges(4,1) * t251 - Ifges(4,4) * t250 + Ifges(4,5) * t267;
t285 = mrSges(5,1) * t177 - mrSges(5,2) * t178 + Ifges(5,5) * t209 + Ifges(5,6) * t208 + Ifges(5,3) * t231 + pkin(4) * t171 + pkin(9) * t163 + t280 * t164 + t276 * t165 + t236 * t203 - t235 * t204;
t139 = -mrSges(4,1) * t211 + mrSges(4,3) * t184 + Ifges(4,4) * t233 + Ifges(4,2) * t232 + Ifges(4,6) * t266 - pkin(3) * t155 - t251 * t220 + t267 * t222 - t285;
t239 = Ifges(3,3) * t267 + (Ifges(3,5) * t278 + Ifges(3,6) * t282) * t301;
t241 = Ifges(3,5) * t267 + (Ifges(3,1) * t278 + Ifges(3,4) * t282) * t301;
t126 = -mrSges(3,1) * t242 + mrSges(3,3) * t225 + Ifges(3,4) * t259 + Ifges(3,2) * t260 + Ifges(3,6) * t266 - pkin(2) * t289 + qJ(3) * t295 + t272 * t135 + t274 * t139 - t239 * t297 + t267 * t241;
t240 = Ifges(3,6) * t267 + (Ifges(3,4) * t278 + Ifges(3,2) * t282) * t301;
t128 = mrSges(3,2) * t242 - mrSges(3,3) * t224 + Ifges(3,1) * t259 + Ifges(3,4) * t260 + Ifges(3,5) * t266 - qJ(3) * t145 + t274 * t135 - t272 * t139 + t239 * t296 - t267 * t240;
t288 = mrSges(4,1) * t183 - mrSges(4,2) * t184 + Ifges(4,5) * t233 + Ifges(4,6) * t232 + Ifges(4,3) * t266 + pkin(3) * t287 + pkin(8) * t294 + t277 * t147 + t281 * t149 + t251 * t221 + t250 * t222;
t130 = t288 + pkin(2) * t145 + (t240 * t278 - t241 * t282) * t301 + mrSges(3,1) * t224 - mrSges(3,2) * t225 + Ifges(3,5) * t259 + Ifges(3,6) * t260 + Ifges(3,3) * t266;
t290 = mrSges(2,1) * t263 - mrSges(2,2) * t264 + Ifges(2,3) * qJDD(1) + pkin(1) * t134 + t126 * t304 + t128 * t305 + t275 * t130 + t138 * t307;
t136 = m(2) * t264 - t284 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t138;
t133 = t275 * t153 + (t143 * t282 + t144 * t278) * t273;
t131 = m(2) * t263 + qJDD(1) * mrSges(2,1) - t284 * mrSges(2,2) + t134;
t124 = -mrSges(2,2) * g(3) - mrSges(2,3) * t263 + Ifges(2,5) * qJDD(1) - t284 * Ifges(2,6) - t278 * t126 + t282 * t128 + (-t133 * t273 - t134 * t275) * pkin(7);
t123 = mrSges(2,1) * g(3) + mrSges(2,3) * t264 + t284 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t133 - t273 * t130 + (pkin(7) * t138 + t126 * t282 + t128 * t278) * t275;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t283 * t124 - t279 * t123 - pkin(6) * (t283 * t131 + t279 * t136), t124, t128, t135, t147, t165; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t279 * t124 + t283 * t123 + pkin(6) * (-t279 * t131 + t283 * t136), t123, t126, t139, t149, t164; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t290, t290, t130, t288, t285, t286;];
m_new = t1;
