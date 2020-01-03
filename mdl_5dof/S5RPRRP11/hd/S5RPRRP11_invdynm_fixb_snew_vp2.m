% Calculate vector of cutting torques with Newton-Euler for
% S5RPRRP11
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
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRRP11_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP11_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP11_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP11_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP11_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP11_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP11_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:53:32
% EndTime: 2019-12-31 18:53:40
% DurationCPUTime: 4.99s
% Computational Cost: add. (55271->286), mult. (129389->346), div. (0->0), fcn. (90328->8), ass. (0->115)
t258 = qJD(1) ^ 2;
t250 = sin(pkin(8));
t251 = cos(pkin(8));
t253 = sin(qJ(3));
t255 = cos(qJ(3));
t269 = t250 * t253 - t251 * t255;
t228 = t269 * qJD(1);
t254 = sin(qJ(1));
t256 = cos(qJ(1));
t237 = -t256 * g(1) - t254 * g(2);
t230 = -t258 * pkin(1) + qJDD(1) * qJ(2) + t237;
t284 = qJD(1) * qJD(2);
t281 = -t251 * g(3) - 0.2e1 * t250 * t284;
t294 = pkin(2) * t251;
t199 = (-pkin(6) * qJDD(1) + t258 * t294 - t230) * t250 + t281;
t219 = -t250 * g(3) + (t230 + 0.2e1 * t284) * t251;
t283 = qJDD(1) * t251;
t247 = t251 ^ 2;
t291 = t247 * t258;
t200 = -pkin(2) * t291 + pkin(6) * t283 + t219;
t167 = t253 * t199 + t255 * t200;
t270 = t250 * t255 + t251 * t253;
t229 = t270 * qJD(1);
t214 = t228 * pkin(3) - t229 * pkin(7);
t257 = qJD(3) ^ 2;
t161 = -t257 * pkin(3) + qJDD(3) * pkin(7) - t228 * t214 + t167;
t246 = t250 ^ 2;
t236 = t254 * g(1) - t256 * g(2);
t276 = qJDD(2) - t236;
t215 = (-pkin(1) - t294) * qJDD(1) + (-qJ(2) + (-t246 - t247) * pkin(6)) * t258 + t276;
t285 = t229 * qJD(3);
t216 = -t269 * qJDD(1) - t285;
t286 = t228 * qJD(3);
t217 = t270 * qJDD(1) - t286;
t163 = (-t217 + t286) * pkin(7) + (-t216 + t285) * pkin(3) + t215;
t252 = sin(qJ(4));
t295 = cos(qJ(4));
t157 = -t252 * t161 + t295 * t163;
t158 = t295 * t161 + t252 * t163;
t220 = -t295 * qJD(3) + t252 * t229;
t221 = t252 * qJD(3) + t295 * t229;
t226 = qJD(4) + t228;
t169 = Ifges(6,5) * t221 + Ifges(6,6) * t226 + Ifges(6,3) * t220;
t172 = Ifges(5,4) * t221 - Ifges(5,2) * t220 + Ifges(5,6) * t226;
t174 = Ifges(5,1) * t221 - Ifges(5,4) * t220 + Ifges(5,5) * t226;
t183 = t221 * qJD(4) - t295 * qJDD(3) + t252 * t217;
t184 = -t220 * qJD(4) + t252 * qJDD(3) + t295 * t217;
t189 = t220 * mrSges(6,1) - t221 * mrSges(6,3);
t213 = qJDD(4) - t216;
t188 = t220 * pkin(4) - t221 * qJ(5);
t225 = t226 ^ 2;
t153 = -t225 * pkin(4) + t213 * qJ(5) + 0.2e1 * qJD(5) * t226 - t220 * t188 + t158;
t155 = -t213 * pkin(4) - t225 * qJ(5) + t221 * t188 + qJDD(5) - t157;
t173 = Ifges(6,1) * t221 + Ifges(6,4) * t226 + Ifges(6,5) * t220;
t267 = mrSges(6,1) * t155 - mrSges(6,3) * t153 - Ifges(6,4) * t184 - Ifges(6,2) * t213 - Ifges(6,6) * t183 - t220 * t173;
t193 = -t220 * mrSges(6,2) + t226 * mrSges(6,3);
t277 = -m(6) * t155 + t213 * mrSges(6,1) + t226 * t193;
t196 = -t226 * mrSges(6,1) + t221 * mrSges(6,2);
t282 = m(6) * t153 + t213 * mrSges(6,3) + t226 * t196;
t297 = -(-t172 + t169) * t221 + mrSges(5,1) * t157 - mrSges(5,2) * t158 + Ifges(5,5) * t184 - Ifges(5,6) * t183 + Ifges(5,3) * t213 + pkin(4) * (-t184 * mrSges(6,2) - t221 * t189 + t277) + qJ(5) * (-t183 * mrSges(6,2) - t220 * t189 + t282) + t220 * t174 - t267;
t209 = t228 * mrSges(4,1) + t229 * mrSges(4,2);
t223 = qJD(3) * mrSges(4,1) - t229 * mrSges(4,3);
t195 = t226 * mrSges(5,1) - t221 * mrSges(5,3);
t288 = -t220 * mrSges(5,1) - t221 * mrSges(5,2) - t189;
t293 = -mrSges(5,3) - mrSges(6,2);
t143 = m(5) * t158 - t213 * mrSges(5,2) + t293 * t183 - t226 * t195 + t288 * t220 + t282;
t194 = -t226 * mrSges(5,2) - t220 * mrSges(5,3);
t145 = m(5) * t157 + t213 * mrSges(5,1) + t293 * t184 + t226 * t194 + t288 * t221 + t277;
t278 = t295 * t143 - t252 * t145;
t135 = m(4) * t167 - qJDD(3) * mrSges(4,2) + t216 * mrSges(4,3) - qJD(3) * t223 - t228 * t209 + t278;
t166 = t255 * t199 - t253 * t200;
t222 = -qJD(3) * mrSges(4,2) - t228 * mrSges(4,3);
t160 = -qJDD(3) * pkin(3) - t257 * pkin(7) + t229 * t214 - t166;
t156 = -0.2e1 * qJD(5) * t221 + (t220 * t226 - t184) * qJ(5) + (t221 * t226 + t183) * pkin(4) + t160;
t149 = m(6) * t156 + t183 * mrSges(6,1) - t184 * mrSges(6,3) + t220 * t193 - t221 * t196;
t261 = -m(5) * t160 - t183 * mrSges(5,1) - t184 * mrSges(5,2) - t220 * t194 - t221 * t195 - t149;
t140 = m(4) * t166 + qJDD(3) * mrSges(4,1) - t217 * mrSges(4,3) + qJD(3) * t222 - t229 * t209 + t261;
t126 = t253 * t135 + t255 * t140;
t218 = -t250 * t230 + t281;
t275 = -mrSges(6,1) * t156 + mrSges(6,2) * t153;
t171 = Ifges(6,4) * t221 + Ifges(6,2) * t226 + Ifges(6,6) * t220;
t290 = -Ifges(5,5) * t221 + Ifges(5,6) * t220 - Ifges(5,3) * t226 - t171;
t129 = -mrSges(5,1) * t160 + mrSges(5,3) * t158 - pkin(4) * t149 + (t173 + t174) * t226 + t290 * t221 + (Ifges(5,6) - Ifges(6,6)) * t213 + (Ifges(5,4) - Ifges(6,5)) * t184 + (-Ifges(5,2) - Ifges(6,3)) * t183 + t275;
t266 = mrSges(6,2) * t155 - mrSges(6,3) * t156 + Ifges(6,1) * t184 + Ifges(6,4) * t213 + Ifges(6,5) * t183 + t226 * t169;
t132 = mrSges(5,2) * t160 - mrSges(5,3) * t157 + Ifges(5,1) * t184 - Ifges(5,4) * t183 + Ifges(5,5) * t213 - qJ(5) * t149 - t226 * t172 + t290 * t220 + t266;
t202 = Ifges(4,4) * t229 - Ifges(4,2) * t228 + Ifges(4,6) * qJD(3);
t203 = Ifges(4,1) * t229 - Ifges(4,4) * t228 + Ifges(4,5) * qJD(3);
t263 = -mrSges(4,1) * t166 + mrSges(4,2) * t167 - Ifges(4,5) * t217 - Ifges(4,6) * t216 - Ifges(4,3) * qJDD(3) - pkin(3) * t261 - pkin(7) * t278 - t295 * t129 - t252 * t132 - t229 * t202 - t228 * t203;
t273 = Ifges(3,4) * t250 + Ifges(3,2) * t251;
t274 = Ifges(3,1) * t250 + Ifges(3,4) * t251;
t296 = -mrSges(3,1) * t218 + mrSges(3,2) * t219 - pkin(2) * t126 - (t250 * t273 - t251 * t274) * t258 + t263;
t292 = mrSges(3,2) * t250;
t137 = t252 * t143 + t295 * t145;
t272 = Ifges(3,5) * t250 + Ifges(3,6) * t251;
t287 = t258 * t272;
t268 = mrSges(3,3) * qJDD(1) + t258 * (-mrSges(3,1) * t251 + t292);
t124 = m(3) * t218 - t268 * t250 + t126;
t279 = t255 * t135 - t253 * t140;
t125 = m(3) * t219 + t268 * t251 + t279;
t280 = -t250 * t124 + t251 * t125;
t201 = Ifges(4,5) * t229 - Ifges(4,6) * t228 + Ifges(4,3) * qJD(3);
t118 = mrSges(4,2) * t215 - mrSges(4,3) * t166 + Ifges(4,1) * t217 + Ifges(4,4) * t216 + Ifges(4,5) * qJDD(3) - pkin(7) * t137 - qJD(3) * t202 - t252 * t129 + t295 * t132 - t228 * t201;
t122 = -mrSges(4,1) * t215 + mrSges(4,3) * t167 + Ifges(4,4) * t217 + Ifges(4,2) * t216 + Ifges(4,6) * qJDD(3) - pkin(3) * t137 + qJD(3) * t203 - t229 * t201 - t297;
t227 = -qJDD(1) * pkin(1) - t258 * qJ(2) + t276;
t264 = m(4) * t215 - t216 * mrSges(4,1) + t217 * mrSges(4,2) + t228 * t222 + t229 * t223 + t137;
t115 = -mrSges(3,1) * t227 + mrSges(3,3) * t219 - pkin(2) * t264 + pkin(6) * t279 + t273 * qJDD(1) + t253 * t118 + t255 * t122 - t250 * t287;
t117 = mrSges(3,2) * t227 - mrSges(3,3) * t218 - pkin(6) * t126 + t274 * qJDD(1) + t255 * t118 - t253 * t122 + t251 * t287;
t262 = -m(3) * t227 + mrSges(3,1) * t283 - t264 + (t246 * t258 + t291) * mrSges(3,3);
t265 = -mrSges(2,2) * t237 + qJ(2) * t280 + t251 * t115 + t250 * t117 + pkin(1) * (-qJDD(1) * t292 + t262) + mrSges(2,1) * t236 + Ifges(2,3) * qJDD(1);
t127 = t262 - t258 * mrSges(2,2) + m(2) * t236 + (mrSges(2,1) - t292) * qJDD(1);
t121 = t251 * t124 + t250 * t125;
t119 = m(2) * t237 - t258 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t280;
t113 = (Ifges(2,6) - t272) * qJDD(1) + t258 * Ifges(2,5) + mrSges(2,3) * t237 + mrSges(2,1) * g(3) - pkin(1) * t121 + t296;
t112 = -mrSges(2,2) * g(3) - mrSges(2,3) * t236 + Ifges(2,5) * qJDD(1) - t258 * Ifges(2,6) - qJ(2) * t121 - t250 * t115 + t251 * t117;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t256 * t112 - t254 * t113 - pkin(5) * (t254 * t119 + t256 * t127), t112, t117, t118, t132, -t220 * t171 + t266; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t254 * t112 + t256 * t113 + pkin(5) * (t256 * t119 - t254 * t127), t113, t115, t122, t129, -t221 * t169 - t267; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t265, t265, t272 * qJDD(1) - t296, -t263, t297, Ifges(6,5) * t184 + Ifges(6,6) * t213 + Ifges(6,3) * t183 + t221 * t171 - t226 * t173 - t275;];
m_new = t1;
