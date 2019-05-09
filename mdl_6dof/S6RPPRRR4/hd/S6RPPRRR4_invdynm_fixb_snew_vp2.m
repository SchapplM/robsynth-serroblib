% Calculate vector of cutting torques with Newton-Euler for
% S6RPPRRR4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-05-05 15:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPPRRR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:40:16
% EndTime: 2019-05-05 15:40:30
% DurationCPUTime: 6.73s
% Computational Cost: add. (141515->298), mult. (257485->364), div. (0->0), fcn. (141485->10), ass. (0->123)
t256 = qJD(1) ^ 2;
t249 = sin(qJ(1));
t253 = cos(qJ(1));
t225 = -t253 * g(1) - t249 * g(2);
t268 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t225;
t277 = -pkin(1) - pkin(2);
t199 = t277 * t256 + t268;
t224 = t249 * g(1) - t253 * g(2);
t267 = -t256 * qJ(2) + qJDD(2) - t224;
t202 = t277 * qJDD(1) + t267;
t244 = sin(pkin(10));
t245 = cos(pkin(10));
t179 = -t244 * t199 + t245 * t202;
t176 = qJDD(1) * pkin(3) - t256 * pkin(7) - t179;
t248 = sin(qJ(4));
t252 = cos(qJ(4));
t274 = qJD(1) * qJD(4);
t272 = t252 * t274;
t219 = -t248 * qJDD(1) - t272;
t273 = t248 * t274;
t220 = -t252 * qJDD(1) + t273;
t163 = (-t219 + t272) * pkin(8) + (-t220 - t273) * pkin(4) + t176;
t180 = t245 * t199 + t244 * t202;
t177 = -t256 * pkin(3) - qJDD(1) * pkin(7) + t180;
t240 = g(3) + qJDD(3);
t173 = t252 * t177 + t248 * t240;
t218 = (pkin(4) * t252 + pkin(8) * t248) * qJD(1);
t230 = t252 * qJD(1);
t255 = qJD(4) ^ 2;
t166 = -t255 * pkin(4) + qJDD(4) * pkin(8) - t218 * t230 + t173;
t247 = sin(qJ(5));
t251 = cos(qJ(5));
t152 = t251 * t163 - t247 * t166;
t275 = qJD(1) * t248;
t215 = t251 * qJD(4) + t247 * t275;
t189 = t215 * qJD(5) + t247 * qJDD(4) + t251 * t219;
t214 = qJDD(5) - t220;
t216 = t247 * qJD(4) - t251 * t275;
t227 = t230 + qJD(5);
t150 = (t215 * t227 - t189) * pkin(9) + (t215 * t216 + t214) * pkin(5) + t152;
t153 = t247 * t163 + t251 * t166;
t188 = -t216 * qJD(5) + t251 * qJDD(4) - t247 * t219;
t198 = t227 * pkin(5) - t216 * pkin(9);
t213 = t215 ^ 2;
t151 = -t213 * pkin(5) + t188 * pkin(9) - t227 * t198 + t153;
t246 = sin(qJ(6));
t250 = cos(qJ(6));
t149 = t246 * t150 + t250 * t151;
t172 = -t248 * t177 + t252 * t240;
t165 = -qJDD(4) * pkin(4) - t255 * pkin(8) - t218 * t275 - t172;
t154 = -t188 * pkin(5) - t213 * pkin(9) + t216 * t198 + t165;
t191 = t246 * t215 + t250 * t216;
t159 = -t191 * qJD(6) + t250 * t188 - t246 * t189;
t190 = t250 * t215 - t246 * t216;
t160 = t190 * qJD(6) + t246 * t188 + t250 * t189;
t226 = qJD(6) + t227;
t167 = Ifges(7,5) * t191 + Ifges(7,6) * t190 + Ifges(7,3) * t226;
t169 = Ifges(7,1) * t191 + Ifges(7,4) * t190 + Ifges(7,5) * t226;
t210 = qJDD(6) + t214;
t138 = -mrSges(7,1) * t154 + mrSges(7,3) * t149 + Ifges(7,4) * t160 + Ifges(7,2) * t159 + Ifges(7,6) * t210 - t191 * t167 + t226 * t169;
t148 = t250 * t150 - t246 * t151;
t168 = Ifges(7,4) * t191 + Ifges(7,2) * t190 + Ifges(7,6) * t226;
t139 = mrSges(7,2) * t154 - mrSges(7,3) * t148 + Ifges(7,1) * t160 + Ifges(7,4) * t159 + Ifges(7,5) * t210 + t190 * t167 - t226 * t168;
t183 = Ifges(6,5) * t216 + Ifges(6,6) * t215 + Ifges(6,3) * t227;
t185 = Ifges(6,1) * t216 + Ifges(6,4) * t215 + Ifges(6,5) * t227;
t181 = -t226 * mrSges(7,2) + t190 * mrSges(7,3);
t182 = t226 * mrSges(7,1) - t191 * mrSges(7,3);
t266 = m(7) * t154 - t159 * mrSges(7,1) + t160 * mrSges(7,2) - t190 * t181 + t191 * t182;
t171 = -t190 * mrSges(7,1) + t191 * mrSges(7,2);
t144 = m(7) * t148 + t210 * mrSges(7,1) - t160 * mrSges(7,3) - t191 * t171 + t226 * t181;
t145 = m(7) * t149 - t210 * mrSges(7,2) + t159 * mrSges(7,3) + t190 * t171 - t226 * t182;
t271 = -t246 * t144 + t250 * t145;
t122 = -mrSges(6,1) * t165 + mrSges(6,3) * t153 + Ifges(6,4) * t189 + Ifges(6,2) * t188 + Ifges(6,6) * t214 - pkin(5) * t266 + pkin(9) * t271 + t250 * t138 + t246 * t139 - t216 * t183 + t227 * t185;
t137 = t250 * t144 + t246 * t145;
t184 = Ifges(6,4) * t216 + Ifges(6,2) * t215 + Ifges(6,6) * t227;
t127 = mrSges(6,2) * t165 - mrSges(6,3) * t152 + Ifges(6,1) * t189 + Ifges(6,4) * t188 + Ifges(6,5) * t214 - pkin(9) * t137 - t246 * t138 + t250 * t139 + t215 * t183 - t227 * t184;
t192 = -t215 * mrSges(6,1) + t216 * mrSges(6,2);
t196 = -t227 * mrSges(6,2) + t215 * mrSges(6,3);
t135 = m(6) * t152 + t214 * mrSges(6,1) - t189 * mrSges(6,3) - t216 * t192 + t227 * t196 + t137;
t197 = t227 * mrSges(6,1) - t216 * mrSges(6,3);
t136 = m(6) * t153 - t214 * mrSges(6,2) + t188 * mrSges(6,3) + t215 * t192 - t227 * t197 + t271;
t133 = -t247 * t135 + t251 * t136;
t146 = -m(6) * t165 + t188 * mrSges(6,1) - t189 * mrSges(6,2) + t215 * t196 - t216 * t197 - t266;
t208 = (Ifges(5,6) * qJD(4)) + (-Ifges(5,4) * t248 - Ifges(5,2) * t252) * qJD(1);
t209 = (Ifges(5,5) * qJD(4)) + (-Ifges(5,1) * t248 - Ifges(5,4) * t252) * qJD(1);
t278 = mrSges(5,1) * t172 - mrSges(5,2) * t173 + Ifges(5,5) * t219 + Ifges(5,6) * t220 + Ifges(5,3) * qJDD(4) + pkin(4) * t146 + pkin(8) * t133 - (t248 * t208 - t252 * t209) * qJD(1) + t251 * t122 + t247 * t127;
t276 = mrSges(2,1) + mrSges(3,1);
t217 = (mrSges(5,1) * t252 - mrSges(5,2) * t248) * qJD(1);
t222 = qJD(4) * mrSges(5,1) + mrSges(5,3) * t275;
t130 = m(5) * t173 - qJDD(4) * mrSges(5,2) + t220 * mrSges(5,3) - qJD(4) * t222 - t217 * t230 + t133;
t223 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t230;
t140 = m(5) * t172 + qJDD(4) * mrSges(5,1) - t219 * mrSges(5,3) + qJD(4) * t223 + t217 * t275 + t146;
t126 = t252 * t130 - t248 * t140;
t121 = m(4) * t180 - t256 * mrSges(4,1) + qJDD(1) * mrSges(4,2) + t126;
t132 = t251 * t135 + t247 * t136;
t131 = -m(5) * t176 + t220 * mrSges(5,1) - t219 * mrSges(5,2) + t222 * t275 - t223 * t230 - t132;
t128 = m(4) * t179 - qJDD(1) * mrSges(4,1) - t256 * mrSges(4,2) + t131;
t118 = t245 * t121 - t244 * t128;
t117 = t244 * t121 + t245 * t128;
t125 = t248 * t130 + t252 * t140;
t203 = -t256 * pkin(1) + t268;
t269 = m(3) * t203 + qJDD(1) * mrSges(3,3) + t118;
t124 = -m(4) * t240 - t125;
t205 = -qJDD(1) * pkin(1) + t267;
t265 = -m(3) * t205 + qJDD(1) * mrSges(3,1) + t256 * mrSges(3,3) - t117;
t264 = -mrSges(7,1) * t148 + mrSges(7,2) * t149 - Ifges(7,5) * t160 - Ifges(7,6) * t159 - Ifges(7,3) * t210 - t191 * t168 + t190 * t169;
t207 = Ifges(5,3) * qJD(4) + (-Ifges(5,5) * t248 - Ifges(5,6) * t252) * qJD(1);
t112 = mrSges(5,2) * t176 - mrSges(5,3) * t172 + Ifges(5,1) * t219 + Ifges(5,4) * t220 + Ifges(5,5) * qJDD(4) - pkin(8) * t132 - qJD(4) * t208 - t247 * t122 + t251 * t127 - t207 * t230;
t257 = mrSges(6,1) * t152 - mrSges(6,2) * t153 + Ifges(6,5) * t189 + Ifges(6,6) * t188 + Ifges(6,3) * t214 + pkin(5) * t137 + t216 * t184 - t215 * t185 - t264;
t119 = -mrSges(5,1) * t176 + mrSges(5,3) * t173 + Ifges(5,4) * t219 + Ifges(5,2) * t220 + Ifges(5,6) * qJDD(4) - pkin(4) * t132 + qJD(4) * t209 + t207 * t275 - t257;
t110 = mrSges(4,2) * t240 - mrSges(4,3) * t179 - Ifges(4,5) * qJDD(1) - t256 * Ifges(4,6) - pkin(7) * t125 + t252 * t112 - t248 * t119;
t111 = -mrSges(4,1) * t240 + mrSges(4,3) * t180 + t256 * Ifges(4,5) - Ifges(4,6) * qJDD(1) - pkin(3) * t125 - t278;
t263 = mrSges(3,2) * t205 + mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) + t256 * Ifges(3,6) - qJ(3) * t117 + t245 * t110 - t244 * t111;
t262 = mrSges(3,2) * t203 - pkin(2) * t124 - qJ(3) * t118 - t244 * t110 - t245 * t111;
t261 = mrSges(4,1) * t179 - mrSges(4,2) * t180 - Ifges(4,3) * qJDD(1) + pkin(3) * t131 + pkin(7) * t126 + t248 * t112 + t252 * t119;
t260 = -mrSges(3,1) * t205 + mrSges(3,3) * t203 + Ifges(3,2) * qJDD(1) - pkin(2) * t117 - t261;
t258 = -mrSges(2,2) * t225 + qJ(2) * (-t256 * mrSges(3,1) + t269) + pkin(1) * t265 + mrSges(2,1) * t224 + Ifges(2,3) * qJDD(1) + t260;
t123 = -m(3) * g(3) + t124;
t114 = m(2) * t224 + qJDD(1) * mrSges(2,1) - t256 * mrSges(2,2) + t265;
t113 = m(2) * t225 - qJDD(1) * mrSges(2,2) - t276 * t256 + t269;
t108 = -mrSges(2,2) * g(3) - mrSges(2,3) * t224 + Ifges(2,5) * qJDD(1) - t256 * Ifges(2,6) - qJ(2) * t123 + t263;
t107 = mrSges(2,3) * t225 - pkin(1) * t123 + (Ifges(3,4) + Ifges(2,5)) * t256 + (Ifges(2,6) - Ifges(3,6)) * qJDD(1) + t276 * g(3) + t262;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t253 * t108 - t249 * t107 - pkin(6) * (t249 * t113 + t253 * t114), t108, t263, t110, t112, t127, t139; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t249 * t108 + t253 * t107 + pkin(6) * (t253 * t113 - t249 * t114), t107, t260, t111, t119, t122, t138; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t258, t258, -mrSges(3,1) * g(3) - t256 * Ifges(3,4) + Ifges(3,6) * qJDD(1) - t262, t261, t278, t257, -t264;];
m_new  = t1;
