% Calculate vector of cutting torques with Newton-Euler for
% S6PRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
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
% Datum: 2019-05-04 21:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRPPRR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPPRR2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 21:50:16
% EndTime: 2019-05-04 21:50:26
% DurationCPUTime: 6.44s
% Computational Cost: add. (97549->254), mult. (171464->318), div. (0->0), fcn. (112963->12), ass. (0->119)
t244 = sin(pkin(10));
t247 = cos(pkin(10));
t226 = t244 * g(1) - t247 * g(2);
t227 = -t247 * g(1) - t244 * g(2);
t240 = -g(3) + qJDD(1);
t251 = sin(qJ(2));
t248 = cos(pkin(6));
t254 = cos(qJ(2));
t279 = t248 * t254;
t245 = sin(pkin(6));
t281 = t245 * t254;
t188 = t226 * t279 - t251 * t227 + t240 * t281;
t186 = qJDD(2) * pkin(2) + t188;
t280 = t248 * t251;
t282 = t245 * t251;
t189 = t226 * t280 + t254 * t227 + t240 * t282;
t256 = qJD(2) ^ 2;
t187 = -t256 * pkin(2) + t189;
t243 = sin(pkin(11));
t246 = cos(pkin(11));
t182 = t243 * t186 + t246 * t187;
t288 = -qJDD(2) * qJ(4) - (2 * qJD(4) * qJD(2)) - t182;
t287 = -pkin(3) - pkin(8);
t181 = t246 * t186 - t243 * t187;
t267 = -t256 * qJ(4) + qJDD(4) - t181;
t176 = t287 * qJDD(2) + t267;
t204 = -t245 * t226 + t248 * t240;
t203 = qJDD(3) + t204;
t250 = sin(qJ(5));
t253 = cos(qJ(5));
t172 = t250 * t176 + t253 * t203;
t221 = (mrSges(6,1) * t250 + mrSges(6,2) * t253) * qJD(2);
t275 = qJD(2) * qJD(5);
t271 = t253 * t275;
t223 = -t250 * qJDD(2) - t271;
t277 = qJD(2) * t253;
t229 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t277;
t222 = (pkin(5) * t250 - pkin(9) * t253) * qJD(2);
t255 = qJD(5) ^ 2;
t276 = t250 * qJD(2);
t168 = -t255 * pkin(5) + qJDD(5) * pkin(9) - t222 * t276 + t172;
t175 = t287 * t256 - t288;
t272 = t250 * t275;
t224 = t253 * qJDD(2) - t272;
t169 = (-t224 + t272) * pkin(9) + (-t223 + t271) * pkin(5) + t175;
t249 = sin(qJ(6));
t252 = cos(qJ(6));
t163 = -t249 * t168 + t252 * t169;
t219 = t252 * qJD(5) - t249 * t277;
t196 = t219 * qJD(6) + t249 * qJDD(5) + t252 * t224;
t220 = t249 * qJD(5) + t252 * t277;
t197 = -t219 * mrSges(7,1) + t220 * mrSges(7,2);
t232 = qJD(6) + t276;
t201 = -t232 * mrSges(7,2) + t219 * mrSges(7,3);
t217 = qJDD(6) - t223;
t161 = m(7) * t163 + t217 * mrSges(7,1) - t196 * mrSges(7,3) - t220 * t197 + t232 * t201;
t164 = t252 * t168 + t249 * t169;
t195 = -t220 * qJD(6) + t252 * qJDD(5) - t249 * t224;
t202 = t232 * mrSges(7,1) - t220 * mrSges(7,3);
t162 = m(7) * t164 - t217 * mrSges(7,2) + t195 * mrSges(7,3) + t219 * t197 - t232 * t202;
t269 = -t249 * t161 + t252 * t162;
t148 = m(6) * t172 - qJDD(5) * mrSges(6,2) + t223 * mrSges(6,3) - qJD(5) * t229 - t221 * t276 + t269;
t278 = t250 * t203;
t171 = t253 * t176 - t278;
t228 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t276;
t167 = -qJDD(5) * pkin(5) - t255 * pkin(9) + t278 + (qJD(2) * t222 - t176) * t253;
t264 = -m(7) * t167 + t195 * mrSges(7,1) - t196 * mrSges(7,2) + t219 * t201 - t220 * t202;
t157 = m(6) * t171 + qJDD(5) * mrSges(6,1) - t224 * mrSges(6,3) + qJD(5) * t228 - t221 * t277 + t264;
t141 = t250 * t148 + t253 * t157;
t179 = -qJDD(2) * pkin(3) + t267;
t266 = -m(5) * t179 + t256 * mrSges(5,3) - t141;
t285 = mrSges(4,1) - mrSges(5,2);
t137 = m(4) * t181 - t256 * mrSges(4,2) + t285 * qJDD(2) + t266;
t151 = t252 * t161 + t249 * t162;
t149 = -m(6) * t175 + t223 * mrSges(6,1) - t224 * mrSges(6,2) - t228 * t276 - t229 * t277 - t151;
t177 = t256 * pkin(3) + t288;
t260 = -m(5) * t177 + t256 * mrSges(5,2) + qJDD(2) * mrSges(5,3) - t149;
t145 = m(4) * t182 - t256 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t260;
t132 = t246 * t137 + t243 * t145;
t130 = m(3) * t188 + qJDD(2) * mrSges(3,1) - t256 * mrSges(3,2) + t132;
t270 = -t243 * t137 + t246 * t145;
t131 = m(3) * t189 - t256 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t270;
t126 = -t251 * t130 + t254 * t131;
t286 = pkin(7) * t126;
t284 = -Ifges(5,4) + Ifges(4,5);
t283 = Ifges(5,5) - Ifges(4,6);
t142 = t253 * t148 - t250 * t157;
t140 = m(5) * t203 + t142;
t268 = m(4) * t203 + t140;
t139 = m(3) * t204 + t268;
t122 = t130 * t279 + t131 * t280 - t245 * t139;
t190 = Ifges(7,5) * t220 + Ifges(7,6) * t219 + Ifges(7,3) * t232;
t192 = Ifges(7,1) * t220 + Ifges(7,4) * t219 + Ifges(7,5) * t232;
t155 = -mrSges(7,1) * t167 + mrSges(7,3) * t164 + Ifges(7,4) * t196 + Ifges(7,2) * t195 + Ifges(7,6) * t217 - t220 * t190 + t232 * t192;
t191 = Ifges(7,4) * t220 + Ifges(7,2) * t219 + Ifges(7,6) * t232;
t156 = mrSges(7,2) * t167 - mrSges(7,3) * t163 + Ifges(7,1) * t196 + Ifges(7,4) * t195 + Ifges(7,5) * t217 + t219 * t190 - t232 * t191;
t209 = (Ifges(6,3) * qJD(5)) + (Ifges(6,5) * t253 - Ifges(6,6) * t250) * qJD(2);
t210 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t253 - Ifges(6,2) * t250) * qJD(2);
t134 = mrSges(6,2) * t175 - mrSges(6,3) * t171 + Ifges(6,1) * t224 + Ifges(6,4) * t223 + Ifges(6,5) * qJDD(5) - pkin(9) * t151 - qJD(5) * t210 - t249 * t155 + t252 * t156 - t209 * t276;
t211 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t253 - Ifges(6,4) * t250) * qJD(2);
t257 = mrSges(7,1) * t163 - mrSges(7,2) * t164 + Ifges(7,5) * t196 + Ifges(7,6) * t195 + Ifges(7,3) * t217 + t220 * t191 - t219 * t192;
t135 = -mrSges(6,1) * t175 + mrSges(6,3) * t172 + Ifges(6,4) * t224 + Ifges(6,2) * t223 + Ifges(6,6) * qJDD(5) - pkin(5) * t151 + qJD(5) * t211 - t209 * t277 - t257;
t262 = -mrSges(5,1) * t177 - pkin(4) * t149 - pkin(8) * t142 - t250 * t134 - t253 * t135;
t118 = mrSges(4,3) * t182 - pkin(3) * t140 - t283 * qJDD(2) - t285 * t203 + t284 * t256 + t262;
t261 = mrSges(6,1) * t171 - mrSges(6,2) * t172 + Ifges(6,5) * t224 + Ifges(6,6) * t223 + Ifges(6,3) * qJDD(5) + pkin(5) * t264 + pkin(9) * t269 + t252 * t155 + t249 * t156 + t210 * t277 + t211 * t276;
t258 = mrSges(5,1) * t179 + pkin(4) * t141 + t261;
t123 = t258 + t283 * t256 + (mrSges(4,2) - mrSges(5,3)) * t203 + t284 * qJDD(2) - mrSges(4,3) * t181 - qJ(4) * t140;
t113 = -mrSges(3,1) * t204 + mrSges(3,3) * t189 + t256 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t268 + qJ(3) * t270 + t246 * t118 + t243 * t123;
t115 = mrSges(3,2) * t204 - mrSges(3,3) * t188 + Ifges(3,5) * qJDD(2) - t256 * Ifges(3,6) - qJ(3) * t132 - t243 * t118 + t246 * t123;
t263 = mrSges(5,2) * t179 - mrSges(5,3) * t177 + Ifges(5,1) * qJDD(2) - pkin(8) * t141 + t253 * t134 - t250 * t135;
t259 = -mrSges(4,2) * t182 + pkin(3) * (-qJDD(2) * mrSges(5,2) + t266) + qJ(4) * t260 + mrSges(4,1) * t181 + Ifges(4,3) * qJDD(2) + t263;
t117 = mrSges(3,1) * t188 - mrSges(3,2) * t189 + Ifges(3,3) * qJDD(2) + pkin(2) * t132 + t259;
t265 = mrSges(2,1) * t226 - mrSges(2,2) * t227 + pkin(1) * t122 + t113 * t281 + t115 * t282 + t248 * t117 + t245 * t286;
t124 = m(2) * t227 + t126;
t121 = t248 * t139 + (t130 * t254 + t131 * t251) * t245;
t119 = m(2) * t226 + t122;
t111 = mrSges(2,2) * t240 - mrSges(2,3) * t226 - t251 * t113 + t254 * t115 + (-t121 * t245 - t122 * t248) * pkin(7);
t110 = -mrSges(2,1) * t240 + mrSges(2,3) * t227 - pkin(1) * t121 - t245 * t117 + (t113 * t254 + t115 * t251 + t286) * t248;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t247 * t111 - t244 * t110 - qJ(1) * (t247 * t119 + t244 * t124), t111, t115, t123, t263, t134, t156; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t244 * t111 + t247 * t110 + qJ(1) * (-t244 * t119 + t247 * t124), t110, t113, t118, mrSges(5,3) * t203 + Ifges(5,4) * qJDD(2) - t256 * Ifges(5,5) - t258, t135, t155; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t265, t265, t117, t259, -mrSges(5,2) * t203 + t256 * Ifges(5,4) + Ifges(5,5) * qJDD(2) - t262, t261, t257;];
m_new  = t1;
