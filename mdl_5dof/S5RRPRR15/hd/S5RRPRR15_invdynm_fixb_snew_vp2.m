% Calculate vector of cutting torques with Newton-Euler for
% S5RRPRR15
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPRR15_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR15_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR15_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR15_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR15_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR15_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR15_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:41:17
% EndTime: 2019-12-31 20:41:27
% DurationCPUTime: 4.92s
% Computational Cost: add. (57520->316), mult. (118149->380), div. (0->0), fcn. (68479->8), ass. (0->119)
t306 = -2 * qJD(3);
t270 = sin(qJ(1));
t274 = cos(qJ(1));
t252 = -t274 * g(1) - t270 * g(2);
t276 = qJD(1) ^ 2;
t224 = -t276 * pkin(1) + qJDD(1) * pkin(6) + t252;
t269 = sin(qJ(2));
t273 = cos(qJ(2));
t207 = -t273 * g(3) - t269 * t224;
t208 = -t269 * g(3) + t273 * t224;
t219 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t269 + Ifges(3,4) * t273) * qJD(1);
t239 = (mrSges(4,2) * t273 - mrSges(4,3) * t269) * qJD(1);
t297 = qJD(1) * qJD(2);
t294 = t273 * t297;
t241 = t269 * qJDD(1) + t294;
t295 = t269 * t297;
t242 = t273 * qJDD(1) - t295;
t298 = qJD(1) * t273;
t248 = -mrSges(4,1) * t298 - qJD(2) * mrSges(4,3);
t258 = t269 * qJD(1);
t238 = (-pkin(2) * t273 - qJ(3) * t269) * qJD(1);
t275 = qJD(2) ^ 2;
t188 = t275 * pkin(2) - qJDD(2) * qJ(3) + qJD(2) * t306 - t238 * t298 - t208;
t250 = pkin(3) * t258 - qJD(2) * pkin(7);
t266 = t273 ^ 2;
t179 = -t266 * t276 * pkin(7) + t242 * pkin(3) + qJD(2) * t250 - t188;
t268 = sin(qJ(4));
t272 = cos(qJ(4));
t237 = t272 * qJD(2) - t268 * t298;
t200 = -t237 * qJD(4) - t268 * qJDD(2) - t272 * t242;
t236 = -t268 * qJD(2) - t272 * t298;
t201 = t236 * qJD(4) + t272 * qJDD(2) - t268 * t242;
t255 = t258 + qJD(4);
t205 = -t255 * mrSges(5,2) + t236 * mrSges(5,3);
t206 = t255 * mrSges(5,1) - t237 * mrSges(5,3);
t209 = t255 * pkin(4) - t237 * pkin(8);
t234 = t236 ^ 2;
t167 = -t200 * pkin(4) - t234 * pkin(8) + t237 * t209 + t179;
t267 = sin(qJ(5));
t271 = cos(qJ(5));
t203 = t267 * t236 + t271 * t237;
t172 = -t203 * qJD(5) + t271 * t200 - t267 * t201;
t202 = t271 * t236 - t267 * t237;
t173 = t202 * qJD(5) + t267 * t200 + t271 * t201;
t253 = qJD(5) + t255;
t191 = -t253 * mrSges(6,2) + t202 * mrSges(6,3);
t192 = t253 * mrSges(6,1) - t203 * mrSges(6,3);
t288 = m(6) * t167 - t172 * mrSges(6,1) + t173 * mrSges(6,2) - t202 * t191 + t203 * t192;
t157 = -m(5) * t179 + t200 * mrSges(5,1) - t201 * mrSges(5,2) + t236 * t205 - t237 * t206 - t288;
t249 = mrSges(4,1) * t258 + qJD(2) * mrSges(4,2);
t280 = -m(4) * t188 + qJDD(2) * mrSges(4,3) + qJD(2) * t249 + t239 * t298 - t157;
t251 = t270 * g(1) - t274 * g(2);
t291 = -qJDD(1) * pkin(1) - t251;
t282 = pkin(2) * t295 + t258 * t306 + (-t241 - t294) * qJ(3) + t291;
t175 = -t250 * t258 + (-pkin(3) * t266 - pkin(6)) * t276 + (-pkin(2) - pkin(7)) * t242 + t282;
t190 = -qJDD(2) * pkin(2) - t275 * qJ(3) + t238 * t258 + qJDD(3) - t207;
t180 = (-t269 * t273 * t276 - qJDD(2)) * pkin(7) + (t241 - t294) * pkin(3) + t190;
t164 = -t268 * t175 + t272 * t180;
t235 = qJDD(4) + t241;
t161 = (t236 * t255 - t201) * pkin(8) + (t236 * t237 + t235) * pkin(4) + t164;
t165 = t272 * t175 + t268 * t180;
t162 = -t234 * pkin(4) + t200 * pkin(8) - t255 * t209 + t165;
t160 = t267 * t161 + t271 * t162;
t181 = Ifges(6,5) * t203 + Ifges(6,6) * t202 + Ifges(6,3) * t253;
t183 = Ifges(6,1) * t203 + Ifges(6,4) * t202 + Ifges(6,5) * t253;
t226 = qJDD(5) + t235;
t147 = -mrSges(6,1) * t167 + mrSges(6,3) * t160 + Ifges(6,4) * t173 + Ifges(6,2) * t172 + Ifges(6,6) * t226 - t203 * t181 + t253 * t183;
t159 = t271 * t161 - t267 * t162;
t182 = Ifges(6,4) * t203 + Ifges(6,2) * t202 + Ifges(6,6) * t253;
t148 = mrSges(6,2) * t167 - mrSges(6,3) * t159 + Ifges(6,1) * t173 + Ifges(6,4) * t172 + Ifges(6,5) * t226 + t202 * t181 - t253 * t182;
t193 = Ifges(5,5) * t237 + Ifges(5,6) * t236 + Ifges(5,3) * t255;
t195 = Ifges(5,1) * t237 + Ifges(5,4) * t236 + Ifges(5,5) * t255;
t185 = -t202 * mrSges(6,1) + t203 * mrSges(6,2);
t155 = m(6) * t159 + t226 * mrSges(6,1) - t173 * mrSges(6,3) - t203 * t185 + t253 * t191;
t156 = m(6) * t160 - t226 * mrSges(6,2) + t172 * mrSges(6,3) + t202 * t185 - t253 * t192;
t292 = -t267 * t155 + t271 * t156;
t131 = -mrSges(5,1) * t179 + mrSges(5,3) * t165 + Ifges(5,4) * t201 + Ifges(5,2) * t200 + Ifges(5,6) * t235 - pkin(4) * t288 + pkin(8) * t292 + t271 * t147 + t267 * t148 - t237 * t193 + t255 * t195;
t146 = t271 * t155 + t267 * t156;
t194 = Ifges(5,4) * t237 + Ifges(5,2) * t236 + Ifges(5,6) * t255;
t133 = mrSges(5,2) * t179 - mrSges(5,3) * t164 + Ifges(5,1) * t201 + Ifges(5,4) * t200 + Ifges(5,5) * t235 - pkin(8) * t146 - t267 * t147 + t271 * t148 + t236 * t193 - t255 * t194;
t204 = -t236 * mrSges(5,1) + t237 * mrSges(5,2);
t143 = m(5) * t164 + t235 * mrSges(5,1) - t201 * mrSges(5,3) - t237 * t204 + t255 * t205 + t146;
t144 = m(5) * t165 - t235 * mrSges(5,2) + t200 * mrSges(5,3) + t236 * t204 - t255 * t206 + t292;
t139 = t272 * t143 + t268 * t144;
t221 = Ifges(4,4) * qJD(2) + (-Ifges(4,2) * t269 - Ifges(4,6) * t273) * qJD(1);
t284 = -mrSges(4,2) * t190 + mrSges(4,3) * t188 - Ifges(4,1) * qJDD(2) + Ifges(4,4) * t241 + Ifges(4,5) * t242 + pkin(7) * t139 + t268 * t131 - t272 * t133 - t221 * t298;
t287 = -m(4) * t190 - t241 * mrSges(4,1) - t139;
t220 = Ifges(4,5) * qJD(2) + (-Ifges(4,6) * t269 - Ifges(4,3) * t273) * qJD(1);
t299 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t269 + Ifges(3,2) * t273) * qJD(1) - t220;
t305 = (-t273 * t219 + t299 * t269) * qJD(1) + mrSges(3,1) * t207 - mrSges(3,2) * t208 + Ifges(3,5) * t241 + Ifges(3,6) * t242 + Ifges(3,3) * qJDD(2) + pkin(2) * (-qJDD(2) * mrSges(4,2) - qJD(2) * t248 - t239 * t258 + t287) + qJ(3) * (t242 * mrSges(4,1) + t280) - t284;
t303 = t276 * pkin(6);
t302 = Ifges(3,4) + Ifges(4,6);
t140 = -t268 * t143 + t272 * t144;
t222 = Ifges(4,1) * qJD(2) + (-Ifges(4,4) * t269 - Ifges(4,5) * t273) * qJD(1);
t300 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t269 + Ifges(3,6) * t273) * qJD(1) + t222;
t240 = (-mrSges(3,1) * t273 + mrSges(3,2) * t269) * qJD(1);
t247 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t298;
t136 = m(3) * t207 - t241 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t247 - t248) * qJD(2) + (-t239 - t240) * t258 + t287;
t246 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t258;
t150 = -qJD(2) * t246 + m(3) * t208 - qJDD(2) * mrSges(3,2) + (mrSges(3,3) + mrSges(4,1)) * t242 + t280 + t240 * t298;
t293 = -t269 * t136 + t273 * t150;
t186 = -t242 * pkin(2) + t282 - t303;
t290 = -m(4) * t186 - t242 * mrSges(4,2) + t249 * t258 - t140;
t137 = -t241 * mrSges(4,3) + t248 * t298 - t290;
t223 = t291 - t303;
t283 = -mrSges(4,1) * t188 + mrSges(4,2) * t186 - pkin(3) * t157 - pkin(7) * t140 - t272 * t131 - t268 * t133;
t125 = -mrSges(3,1) * t223 + mrSges(3,3) * t208 - pkin(2) * t137 + (Ifges(3,2) + Ifges(4,3)) * t242 + t302 * t241 + (Ifges(3,6) - Ifges(4,5)) * qJDD(2) + (t219 - t221) * qJD(2) - t300 * t258 + t283;
t285 = -mrSges(6,1) * t159 + mrSges(6,2) * t160 - Ifges(6,5) * t173 - Ifges(6,6) * t172 - Ifges(6,3) * t226 - t203 * t182 + t202 * t183;
t281 = -mrSges(5,1) * t164 + mrSges(5,2) * t165 - Ifges(5,5) * t201 - Ifges(5,6) * t200 - Ifges(5,3) * t235 - pkin(4) * t146 - t237 * t194 + t236 * t195 + t285;
t278 = -mrSges(4,1) * t190 + mrSges(4,3) * t186 - pkin(3) * t139 + t281;
t127 = mrSges(3,2) * t223 - mrSges(3,3) * t207 - qJ(3) * t137 - t278 + t302 * t242 + (Ifges(3,1) + Ifges(4,2)) * t241 + (Ifges(3,5) - Ifges(4,4)) * qJDD(2) - t299 * qJD(2) + t300 * t298;
t279 = -m(3) * t223 + t247 * t298 + t242 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t241 + (-t246 * t269 - t248 * t273) * qJD(1) + t290;
t286 = mrSges(2,1) * t251 - mrSges(2,2) * t252 + Ifges(2,3) * qJDD(1) + pkin(1) * t279 + pkin(6) * t293 + t273 * t125 + t269 * t127;
t134 = m(2) * t251 + qJDD(1) * mrSges(2,1) - t276 * mrSges(2,2) + t279;
t130 = t273 * t136 + t269 * t150;
t128 = m(2) * t252 - t276 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t293;
t123 = mrSges(2,1) * g(3) + mrSges(2,3) * t252 + t276 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t130 - t305;
t122 = -mrSges(2,2) * g(3) - mrSges(2,3) * t251 + Ifges(2,5) * qJDD(1) - t276 * Ifges(2,6) - pkin(6) * t130 - t269 * t125 + t273 * t127;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t274 * t122 - t270 * t123 - pkin(5) * (t270 * t128 + t274 * t134), t122, t127, -t220 * t258 - t284, t133, t148; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t270 * t122 + t274 * t123 + pkin(5) * (t274 * t128 - t270 * t134), t123, t125, Ifges(4,4) * qJDD(2) - Ifges(4,2) * t241 - Ifges(4,6) * t242 - qJD(2) * t220 - t222 * t298 + t278, t131, t147; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t286, t286, t305, Ifges(4,5) * qJDD(2) - Ifges(4,6) * t241 - Ifges(4,3) * t242 + qJD(2) * t221 + t222 * t258 - t283, -t281, -t285;];
m_new = t1;
