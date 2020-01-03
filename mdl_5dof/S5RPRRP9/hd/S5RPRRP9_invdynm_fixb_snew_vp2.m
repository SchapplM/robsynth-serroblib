% Calculate vector of cutting torques with Newton-Euler for
% S5RPRRP9
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
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRRP9_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP9_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP9_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP9_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP9_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP9_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP9_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:48:42
% EndTime: 2019-12-31 18:48:52
% DurationCPUTime: 6.35s
% Computational Cost: add. (70339->287), mult. (170742->347), div. (0->0), fcn. (123288->8), ass. (0->114)
t264 = qJD(1) ^ 2;
t261 = sin(qJ(1));
t263 = cos(qJ(1));
t238 = -t263 * g(1) - t261 * g(2);
t231 = -t264 * pkin(1) + qJDD(1) * qJ(2) + t238;
t257 = sin(pkin(8));
t258 = cos(pkin(8));
t290 = qJD(1) * qJD(2);
t287 = -t258 * g(3) - 0.2e1 * t257 * t290;
t298 = pkin(2) * t258;
t203 = (-pkin(6) * qJDD(1) + t264 * t298 - t231) * t257 + t287;
t221 = -t257 * g(3) + (t231 + 0.2e1 * t290) * t258;
t289 = qJDD(1) * t258;
t253 = t258 ^ 2;
t295 = t253 * t264;
t204 = -pkin(2) * t295 + pkin(6) * t289 + t221;
t260 = sin(qJ(3));
t262 = cos(qJ(3));
t177 = t262 * t203 - t260 * t204;
t276 = t257 * t262 + t258 * t260;
t275 = -t257 * t260 + t258 * t262;
t229 = t275 * qJD(1);
t291 = t229 * qJD(3);
t219 = t276 * qJDD(1) + t291;
t230 = t276 * qJD(1);
t157 = (-t219 + t291) * pkin(7) + (t229 * t230 + qJDD(3)) * pkin(3) + t177;
t178 = t260 * t203 + t262 * t204;
t218 = -t230 * qJD(3) + t275 * qJDD(1);
t224 = qJD(3) * pkin(3) - t230 * pkin(7);
t228 = t229 ^ 2;
t159 = -t228 * pkin(3) + t218 * pkin(7) - qJD(3) * t224 + t178;
t259 = sin(qJ(4));
t299 = cos(qJ(4));
t155 = t259 * t157 + t299 * t159;
t209 = t259 * t229 + t299 * t230;
t174 = t209 * qJD(4) - t299 * t218 + t259 * t219;
t254 = qJD(3) + qJD(4);
t199 = t254 * mrSges(5,1) - t209 * mrSges(5,3);
t208 = -t299 * t229 + t259 * t230;
t251 = qJDD(3) + qJDD(4);
t190 = t208 * pkin(4) - t209 * qJ(5);
t250 = t254 ^ 2;
t148 = -t250 * pkin(4) + t251 * qJ(5) + 0.2e1 * qJD(5) * t254 - t208 * t190 + t155;
t200 = -t254 * mrSges(6,1) + t209 * mrSges(6,2);
t288 = m(6) * t148 + t251 * mrSges(6,3) + t254 * t200;
t191 = t208 * mrSges(6,1) - t209 * mrSges(6,3);
t293 = -t208 * mrSges(5,1) - t209 * mrSges(5,2) - t191;
t297 = -mrSges(5,3) - mrSges(6,2);
t138 = m(5) * t155 - t251 * mrSges(5,2) + t297 * t174 - t254 * t199 + t293 * t208 + t288;
t154 = t299 * t157 - t259 * t159;
t175 = -t208 * qJD(4) + t259 * t218 + t299 * t219;
t198 = -t254 * mrSges(5,2) - t208 * mrSges(5,3);
t150 = -t251 * pkin(4) - t250 * qJ(5) + t209 * t190 + qJDD(5) - t154;
t201 = -t208 * mrSges(6,2) + t254 * mrSges(6,3);
t283 = -m(6) * t150 + t251 * mrSges(6,1) + t254 * t201;
t140 = m(5) * t154 + t251 * mrSges(5,1) + t297 * t175 + t254 * t198 + t293 * t209 + t283;
t133 = t259 * t138 + t299 * t140;
t213 = -t229 * mrSges(4,1) + t230 * mrSges(4,2);
t222 = -qJD(3) * mrSges(4,2) + t229 * mrSges(4,3);
t130 = m(4) * t177 + qJDD(3) * mrSges(4,1) - t219 * mrSges(4,3) + qJD(3) * t222 - t230 * t213 + t133;
t223 = qJD(3) * mrSges(4,1) - t230 * mrSges(4,3);
t284 = t299 * t138 - t259 * t140;
t131 = m(4) * t178 - qJDD(3) * mrSges(4,2) + t218 * mrSges(4,3) - qJD(3) * t223 + t229 * t213 + t284;
t124 = t262 * t130 + t260 * t131;
t220 = -t257 * t231 + t287;
t206 = Ifges(4,4) * t230 + Ifges(4,2) * t229 + Ifges(4,6) * qJD(3);
t207 = Ifges(4,1) * t230 + Ifges(4,4) * t229 + Ifges(4,5) * qJD(3);
t183 = Ifges(5,4) * t209 - Ifges(5,2) * t208 + Ifges(5,6) * t254;
t185 = Ifges(5,1) * t209 - Ifges(5,4) * t208 + Ifges(5,5) * t254;
t180 = Ifges(6,5) * t209 + Ifges(6,6) * t254 + Ifges(6,3) * t208;
t184 = Ifges(6,1) * t209 + Ifges(6,4) * t254 + Ifges(6,5) * t208;
t270 = mrSges(6,1) * t150 - mrSges(6,3) * t148 - Ifges(6,4) * t175 - Ifges(6,2) * t251 - Ifges(6,6) * t174 + t209 * t180 - t208 * t184;
t268 = mrSges(5,2) * t155 - t208 * t185 - qJ(5) * (-t174 * mrSges(6,2) - t208 * t191 + t288) - pkin(4) * (-t175 * mrSges(6,2) - t209 * t191 + t283) - mrSges(5,1) * t154 - t209 * t183 + Ifges(5,6) * t174 - Ifges(5,5) * t175 - Ifges(5,3) * t251 + t270;
t266 = -mrSges(4,1) * t177 + mrSges(4,2) * t178 - Ifges(4,5) * t219 - Ifges(4,6) * t218 - Ifges(4,3) * qJDD(3) - pkin(3) * t133 - t230 * t206 + t229 * t207 + t268;
t279 = Ifges(3,4) * t257 + Ifges(3,2) * t258;
t280 = Ifges(3,1) * t257 + Ifges(3,4) * t258;
t300 = -mrSges(3,1) * t220 + mrSges(3,2) * t221 - pkin(2) * t124 - (t257 * t279 - t258 * t280) * t264 + t266;
t296 = mrSges(3,2) * t257;
t182 = Ifges(6,4) * t209 + Ifges(6,2) * t254 + Ifges(6,6) * t208;
t294 = -Ifges(5,5) * t209 + Ifges(5,6) * t208 - Ifges(5,3) * t254 - t182;
t278 = Ifges(3,5) * t257 + Ifges(3,6) * t258;
t292 = t264 * t278;
t237 = t261 * g(1) - t263 * g(2);
t274 = mrSges(3,3) * qJDD(1) + t264 * (-mrSges(3,1) * t258 + t296);
t122 = m(3) * t220 - t274 * t257 + t124;
t285 = -t260 * t130 + t262 * t131;
t123 = m(3) * t221 + t274 * t258 + t285;
t286 = -t257 * t122 + t258 * t123;
t282 = qJDD(2) - t237;
t252 = t257 ^ 2;
t217 = (-pkin(1) - t298) * qJDD(1) + (-qJ(2) + (-t252 - t253) * pkin(6)) * t264 + t282;
t163 = -t218 * pkin(3) - t228 * pkin(7) + t230 * t224 + t217;
t152 = -0.2e1 * qJD(5) * t209 + (t208 * t254 - t175) * qJ(5) + (t209 * t254 + t174) * pkin(4) + t163;
t281 = -mrSges(6,1) * t152 + mrSges(6,2) * t148;
t141 = m(6) * t152 + t174 * mrSges(6,1) - t175 * mrSges(6,3) - t209 * t200 + t208 * t201;
t273 = mrSges(6,2) * t150 - mrSges(6,3) * t152 + Ifges(6,1) * t175 + Ifges(6,4) * t251 + Ifges(6,5) * t174 + t254 * t180;
t125 = -mrSges(5,1) * t163 + mrSges(5,3) * t155 - pkin(4) * t141 + (t184 + t185) * t254 + (Ifges(5,6) - Ifges(6,6)) * t251 + t294 * t209 + (Ifges(5,4) - Ifges(6,5)) * t175 + (-Ifges(5,2) - Ifges(6,3)) * t174 + t281;
t126 = mrSges(5,2) * t163 - mrSges(5,3) * t154 + Ifges(5,1) * t175 - Ifges(5,4) * t174 + Ifges(5,5) * t251 - qJ(5) * t141 - t254 * t183 + t294 * t208 + t273;
t205 = Ifges(4,5) * t230 + Ifges(4,6) * t229 + Ifges(4,3) * qJD(3);
t271 = m(5) * t163 + t174 * mrSges(5,1) + t175 * mrSges(5,2) + t208 * t198 + t209 * t199 + t141;
t119 = -mrSges(4,1) * t217 + mrSges(4,3) * t178 + Ifges(4,4) * t219 + Ifges(4,2) * t218 + Ifges(4,6) * qJDD(3) - pkin(3) * t271 + pkin(7) * t284 + qJD(3) * t207 + t299 * t125 + t259 * t126 - t230 * t205;
t120 = mrSges(4,2) * t217 - mrSges(4,3) * t177 + Ifges(4,1) * t219 + Ifges(4,4) * t218 + Ifges(4,5) * qJDD(3) - pkin(7) * t133 - qJD(3) * t206 - t259 * t125 + t299 * t126 + t229 * t205;
t227 = -qJDD(1) * pkin(1) - t264 * qJ(2) + t282;
t269 = m(4) * t217 - t218 * mrSges(4,1) + t219 * mrSges(4,2) - t229 * t222 + t230 * t223 + t271;
t112 = -mrSges(3,1) * t227 + mrSges(3,3) * t221 - pkin(2) * t269 + pkin(6) * t285 + t279 * qJDD(1) + t262 * t119 + t260 * t120 - t257 * t292;
t114 = mrSges(3,2) * t227 - mrSges(3,3) * t220 - pkin(6) * t124 + t280 * qJDD(1) - t260 * t119 + t262 * t120 + t258 * t292;
t267 = -m(3) * t227 + mrSges(3,1) * t289 - t269 + (t252 * t264 + t295) * mrSges(3,3);
t272 = -mrSges(2,2) * t238 + qJ(2) * t286 + t258 * t112 + t257 * t114 + pkin(1) * (-qJDD(1) * t296 + t267) + mrSges(2,1) * t237 + Ifges(2,3) * qJDD(1);
t134 = (mrSges(2,1) - t296) * qJDD(1) + t267 - t264 * mrSges(2,2) + m(2) * t237;
t118 = t258 * t122 + t257 * t123;
t116 = m(2) * t238 - t264 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t286;
t115 = mrSges(2,1) * g(3) + (Ifges(2,6) - t278) * qJDD(1) + t264 * Ifges(2,5) + mrSges(2,3) * t238 - pkin(1) * t118 + t300;
t110 = -mrSges(2,2) * g(3) - mrSges(2,3) * t237 + Ifges(2,5) * qJDD(1) - t264 * Ifges(2,6) - qJ(2) * t118 - t257 * t112 + t258 * t114;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t263 * t110 - t261 * t115 - pkin(5) * (t261 * t116 + t263 * t134), t110, t114, t120, t126, -t208 * t182 + t273; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t261 * t110 + t263 * t115 + pkin(5) * (t263 * t116 - t261 * t134), t115, t112, t119, t125, -t270; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t272, t272, t278 * qJDD(1) - t300, -t266, -t268, Ifges(6,5) * t175 + Ifges(6,6) * t251 + Ifges(6,3) * t174 + t209 * t182 - t254 * t184 - t281;];
m_new = t1;
