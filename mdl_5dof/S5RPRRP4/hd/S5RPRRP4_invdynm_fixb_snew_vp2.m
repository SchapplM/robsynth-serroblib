% Calculate vector of cutting torques with Newton-Euler for
% S5RPRRP4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% m [6x1]
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
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRRP4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP4_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP4_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP4_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:31:56
% EndTime: 2022-01-23 09:32:04
% DurationCPUTime: 4.90s
% Computational Cost: add. (49974->283), mult. (120788->358), div. (0->0), fcn. (79109->8), ass. (0->116)
t266 = sin(qJ(1));
t269 = cos(qJ(1));
t248 = -t269 * g(1) - t266 * g(2);
t270 = qJD(1) ^ 2;
t311 = -t270 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t248;
t264 = sin(qJ(4));
t265 = sin(qJ(3));
t267 = cos(qJ(4));
t268 = cos(qJ(3));
t262 = sin(pkin(8));
t303 = qJD(1) * t262;
t221 = (-t264 * t268 - t265 * t267) * t303;
t299 = qJD(1) * qJD(3);
t230 = (-qJDD(1) * t265 - t268 * t299) * t262;
t231 = (qJDD(1) * t268 - t265 * t299) * t262;
t190 = t221 * qJD(4) + t264 * t230 + t267 * t231;
t222 = (-t264 * t265 + t267 * t268) * t303;
t203 = -t221 * mrSges(6,1) + t222 * mrSges(6,2);
t263 = cos(pkin(8));
t214 = -t262 * g(3) + t311 * t263;
t288 = -pkin(2) * t263 - pkin(6) * t262;
t239 = t288 * qJD(1);
t301 = t263 * qJD(1);
t202 = t239 * t301 + t214;
t247 = t266 * g(1) - t269 * g(2);
t281 = -t270 * qJ(2) + qJDD(2) - t247;
t215 = (-pkin(1) + t288) * qJDD(1) + t281;
t212 = t268 * t215;
t298 = t263 * qJDD(1);
t250 = qJDD(3) - t298;
t251 = qJD(3) - t301;
t306 = t262 ^ 2 * t270;
t167 = t250 * pkin(3) - t231 * pkin(7) + t212 + (-pkin(3) * t268 * t306 - pkin(7) * t251 * t303 - t202) * t265;
t171 = t268 * t202 + t265 * t215;
t293 = t268 * t303;
t229 = t251 * pkin(3) - pkin(7) * t293;
t297 = t265 ^ 2 * t306;
t168 = -pkin(3) * t297 + t230 * pkin(7) - t251 * t229 + t171;
t159 = t267 * t167 - t264 * t168;
t246 = qJDD(4) + t250;
t249 = qJD(4) + t251;
t154 = -0.2e1 * qJD(5) * t222 + (t221 * t249 - t190) * qJ(5) + (t221 * t222 + t246) * pkin(4) + t159;
t206 = -t249 * mrSges(6,2) + t221 * mrSges(6,3);
t296 = m(6) * t154 + t246 * mrSges(6,1) + t249 * t206;
t151 = -t190 * mrSges(6,3) - t222 * t203 + t296;
t160 = t264 * t167 + t267 * t168;
t189 = -t222 * qJD(4) + t267 * t230 - t264 * t231;
t194 = Ifges(5,4) * t222 + Ifges(5,2) * t221 + Ifges(5,6) * t249;
t195 = Ifges(6,1) * t222 + Ifges(6,4) * t221 + Ifges(6,5) * t249;
t196 = Ifges(5,1) * t222 + Ifges(5,4) * t221 + Ifges(5,5) * t249;
t208 = t249 * pkin(4) - t222 * qJ(5);
t220 = t221 ^ 2;
t157 = -t220 * pkin(4) + t189 * qJ(5) + 0.2e1 * qJD(5) * t221 - t249 * t208 + t160;
t193 = Ifges(6,4) * t222 + Ifges(6,2) * t221 + Ifges(6,6) * t249;
t279 = -mrSges(6,1) * t154 + mrSges(6,2) * t157 - Ifges(6,5) * t190 - Ifges(6,6) * t189 - Ifges(6,3) * t246 - t222 * t193;
t310 = mrSges(5,1) * t159 - mrSges(5,2) * t160 + Ifges(5,5) * t190 + Ifges(5,6) * t189 + Ifges(5,3) * t246 + pkin(4) * t151 + t222 * t194 - t279 + (-t196 - t195) * t221;
t204 = -t221 * mrSges(5,1) + t222 * mrSges(5,2);
t207 = -t249 * mrSges(5,2) + t221 * mrSges(5,3);
t144 = m(5) * t159 + t246 * mrSges(5,1) + t249 * t207 + (-t203 - t204) * t222 + (-mrSges(5,3) - mrSges(6,3)) * t190 + t296;
t209 = t249 * mrSges(6,1) - t222 * mrSges(6,3);
t210 = t249 * mrSges(5,1) - t222 * mrSges(5,3);
t295 = m(6) * t157 + t189 * mrSges(6,3) + t221 * t203;
t147 = m(5) * t160 + t189 * mrSges(5,3) + t221 * t204 + (-t209 - t210) * t249 + (-mrSges(5,2) - mrSges(6,2)) * t246 + t295;
t142 = t267 * t144 + t264 * t147;
t170 = -t265 * t202 + t212;
t309 = -mrSges(4,1) * t170 + mrSges(4,2) * t171 - Ifges(4,5) * t231 - Ifges(4,6) * t230 - Ifges(4,3) * t250 - pkin(3) * t142 - t310;
t213 = -t263 * g(3) - t311 * t262;
t201 = t239 * t303 - t213;
t169 = -t230 * pkin(3) - pkin(7) * t297 + t229 * t293 + t201;
t191 = Ifges(6,5) * t222 + Ifges(6,6) * t221 + Ifges(6,3) * t249;
t192 = Ifges(5,5) * t222 + Ifges(5,6) * t221 + Ifges(5,3) * t249;
t163 = -t189 * pkin(4) - t220 * qJ(5) + t222 * t208 + qJDD(5) + t169;
t280 = -mrSges(6,1) * t163 + mrSges(6,3) * t157 + Ifges(6,4) * t190 + Ifges(6,2) * t189 + Ifges(6,6) * t246 + t249 * t195;
t285 = m(6) * t163 - t189 * mrSges(6,1) + t190 * mrSges(6,2) - t221 * t206 + t222 * t209;
t137 = Ifges(5,4) * t190 + Ifges(5,2) * t189 + Ifges(5,6) * t246 + t249 * t196 - mrSges(5,1) * t169 + mrSges(5,3) * t160 - pkin(4) * t285 + qJ(5) * (-t246 * mrSges(6,2) - t249 * t209 + t295) + (-t192 - t191) * t222 + t280;
t278 = mrSges(6,2) * t163 - mrSges(6,3) * t154 + Ifges(6,1) * t190 + Ifges(6,4) * t189 + Ifges(6,5) * t246 + t221 * t191;
t141 = mrSges(5,2) * t169 - mrSges(5,3) * t159 + Ifges(5,1) * t190 + Ifges(5,4) * t189 + Ifges(5,5) * t246 - qJ(5) * t151 + t221 * t192 + (-t193 - t194) * t249 + t278;
t216 = Ifges(4,3) * t251 + (Ifges(4,5) * t268 - Ifges(4,6) * t265) * t303;
t218 = Ifges(4,5) * t251 + (Ifges(4,1) * t268 - Ifges(4,4) * t265) * t303;
t275 = m(5) * t169 - t189 * mrSges(5,1) + t190 * mrSges(5,2) - t221 * t207 + t222 * t210 + t285;
t289 = -t264 * t144 + t267 * t147;
t126 = -mrSges(4,1) * t201 + mrSges(4,3) * t171 + Ifges(4,4) * t231 + Ifges(4,2) * t230 + Ifges(4,6) * t250 - pkin(3) * t275 + pkin(7) * t289 + t267 * t137 + t264 * t141 - t216 * t293 + t251 * t218;
t217 = Ifges(4,6) * t251 + (Ifges(4,4) * t268 - Ifges(4,2) * t265) * t303;
t294 = t265 * t303;
t127 = mrSges(4,2) * t201 - mrSges(4,3) * t170 + Ifges(4,1) * t231 + Ifges(4,4) * t230 + Ifges(4,5) * t250 - pkin(7) * t142 - t264 * t137 + t267 * t141 - t216 * t294 - t251 * t217;
t226 = -t251 * mrSges(4,2) - mrSges(4,3) * t294;
t228 = (mrSges(4,1) * t265 + mrSges(4,2) * t268) * t303;
t139 = m(4) * t170 + t250 * mrSges(4,1) - t231 * mrSges(4,3) + t251 * t226 - t228 * t293 + t142;
t227 = t251 * mrSges(4,1) - mrSges(4,3) * t293;
t140 = m(4) * t171 - t250 * mrSges(4,2) + t230 * mrSges(4,3) - t251 * t227 - t228 * t294 + t289;
t136 = -t265 * t139 + t268 * t140;
t272 = -m(4) * t201 + t230 * mrSges(4,1) - t231 * mrSges(4,2) - t275;
t283 = -t226 * t265 - t227 * t268;
t287 = Ifges(3,1) * t262 + Ifges(3,4) * t263;
t308 = -((Ifges(3,4) * t262 + Ifges(3,2) * t263) * t303 - t287 * t301) * qJD(1) - mrSges(3,1) * t213 + mrSges(3,2) * t214 - pkin(2) * (t283 * t303 + t272) - pkin(6) * t136 - t268 * t126 - t265 * t127;
t307 = mrSges(3,2) * t262;
t302 = qJDD(1) * mrSges(3,3);
t235 = (-mrSges(3,1) * t263 + t307) * qJD(1);
t133 = m(3) * t214 + (qJD(1) * t235 + t302) * t263 + t136;
t148 = t272 + (-t302 + (-t235 + t283) * qJD(1)) * t262 + m(3) * t213;
t290 = t263 * t133 - t262 * t148;
t286 = Ifges(3,5) * t262 + Ifges(3,6) * t263;
t135 = t268 * t139 + t265 * t140;
t284 = t217 * t268 + t218 * t265;
t233 = -qJDD(1) * pkin(1) + t281;
t236 = t286 * qJD(1);
t123 = mrSges(3,2) * t233 - mrSges(3,3) * t213 - pkin(6) * t135 + t287 * qJDD(1) - t265 * t126 + t268 * t127 + t236 * t301;
t125 = -pkin(2) * t135 + (Ifges(3,4) * qJDD(1) + (-t236 - t284) * qJD(1)) * t262 + Ifges(3,2) * t298 + mrSges(3,3) * t214 - mrSges(3,1) * t233 + t309;
t276 = -m(3) * t233 + mrSges(3,1) * t298 - t135 + (t263 ^ 2 * t270 + t306) * mrSges(3,3);
t277 = -mrSges(2,2) * t248 + qJ(2) * t290 + t262 * t123 + t263 * t125 + pkin(1) * (-qJDD(1) * t307 + t276) + mrSges(2,1) * t247 + Ifges(2,3) * qJDD(1);
t131 = m(2) * t247 - t270 * mrSges(2,2) + (mrSges(2,1) - t307) * qJDD(1) + t276;
t130 = t262 * t133 + t263 * t148;
t128 = m(2) * t248 - t270 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t290;
t121 = mrSges(2,1) * g(3) + mrSges(2,3) * t248 + t270 * Ifges(2,5) - pkin(1) * t130 + (Ifges(2,6) - t286) * qJDD(1) + t308;
t120 = -mrSges(2,2) * g(3) - mrSges(2,3) * t247 + Ifges(2,5) * qJDD(1) - t270 * Ifges(2,6) - qJ(2) * t130 + t263 * t123 - t262 * t125;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t269 * t120 - t266 * t121 - pkin(5) * (t266 * t128 + t269 * t131), t120, t123, t127, t141, -t249 * t193 + t278; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t266 * t120 + t269 * t121 + pkin(5) * (t269 * t128 - t266 * t131), t121, t125, t126, t137, -t222 * t191 + t280; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t277, t277, t286 * qJDD(1) - t308, t284 * t303 - t309, t310, -t221 * t195 - t279;];
m_new = t1;
