% Calculate vector of inverse dynamics joint torques for
% S6RPPPRR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPPRR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:35:18
% EndTime: 2019-03-09 01:35:26
% DurationCPUTime: 6.24s
% Computational Cost: add. (3154->428), mult. (5274->563), div. (0->0), fcn. (2757->8), ass. (0->204)
t237 = m(6) + m(5);
t111 = sin(qJ(6));
t234 = t111 / 0.2e1;
t255 = m(7) + t237;
t275 = m(4) + t255;
t113 = cos(qJ(6));
t191 = qJD(5) * t113;
t114 = cos(qJ(5));
t194 = qJD(1) * t114;
t73 = t111 * t194 + t191;
t274 = t73 / 0.2e1;
t184 = t111 * qJD(5);
t74 = t113 * t194 - t184;
t239 = -t74 / 0.2e1;
t112 = sin(qJ(5));
t195 = qJD(1) * t112;
t90 = qJD(6) - t195;
t273 = t90 / 0.2e1;
t272 = -qJD(1) / 0.2e1;
t176 = mrSges(6,3) * t194;
t222 = -qJD(5) * mrSges(6,1) - mrSges(7,1) * t73 - mrSges(7,2) * t74 - t176;
t271 = t114 * t222;
t182 = qJD(1) * qJD(5);
t78 = -qJDD(1) * t114 + t112 * t182;
t30 = qJD(6) * t73 + qJDD(5) * t111 + t113 * t78;
t31 = qJD(6) * t74 + qJDD(5) * t113 - t111 * t78;
t13 = -mrSges(7,1) * t31 + mrSges(7,2) * t30;
t223 = -qJDD(5) * mrSges(6,1) + mrSges(6,3) * t78 + t13;
t270 = t223 * t114;
t108 = sin(pkin(9));
t109 = cos(pkin(9));
t196 = qJ(2) * qJD(1);
t115 = -pkin(1) - pkin(2);
t87 = qJD(1) * t115 + qJD(2);
t54 = -t108 * t196 + t109 * t87;
t49 = qJD(1) * pkin(3) + qJD(4) - t54;
t44 = qJD(1) * pkin(7) + t49;
t34 = qJD(3) * t114 + t112 * t44;
t25 = qJD(5) * pkin(8) + t34;
t106 = qJD(1) * qJ(4);
t230 = pkin(8) * t114;
t156 = pkin(5) * t112 - t230;
t205 = qJ(2) * t109;
t128 = -t156 + t205;
t216 = t108 * t87;
t35 = qJD(1) * t128 - t106 + t216;
t10 = t111 * t35 + t113 * t25;
t86 = qJDD(1) * t115 + qJDD(2);
t183 = qJD(1) * qJD(2);
t88 = qJDD(1) * qJ(2) + t183;
t46 = t108 * t86 + t109 * t88;
t37 = -qJDD(1) * qJ(4) - qJD(1) * qJD(4) + t46;
t79 = qJDD(1) * t112 + t114 * t182;
t14 = -pkin(5) * t79 - pkin(8) * t78 + t37;
t190 = qJD(5) * t114;
t192 = qJD(5) * t112;
t45 = -t108 * t88 + t109 * t86;
t40 = qJDD(1) * pkin(3) + qJDD(4) - t45;
t36 = qJDD(1) * pkin(7) + t40;
t11 = -qJD(3) * t192 + t114 * qJDD(3) + t112 * t36 + t44 * t190;
t7 = qJDD(5) * pkin(8) + t11;
t9 = -t111 * t25 + t113 * t35;
t1 = qJD(6) * t9 + t111 * t14 + t113 * t7;
t2 = -qJD(6) * t10 - t111 * t7 + t113 * t14;
t155 = t1 * t113 - t111 * t2;
t188 = qJD(6) * t113;
t189 = qJD(6) * t111;
t269 = -t10 * t189 - t9 * t188 + t155;
t232 = sin(qJ(1));
t233 = cos(qJ(1));
t70 = -t108 * t232 - t109 * t233;
t71 = t108 * t233 - t109 * t232;
t256 = -g(1) * t71 + g(2) * t70;
t243 = t30 / 0.2e1;
t242 = t31 / 0.2e1;
t69 = qJDD(6) - t79;
t241 = t69 / 0.2e1;
t268 = t78 / 0.2e1;
t267 = t79 / 0.2e1;
t266 = -m(6) - m(7);
t33 = -qJD(3) * t112 + t114 * t44;
t24 = -qJD(5) * pkin(5) - t33;
t265 = m(7) * t24;
t264 = mrSges(3,1) + mrSges(2,1);
t225 = -mrSges(4,1) + mrSges(5,2);
t224 = mrSges(4,2) - mrSges(5,3);
t263 = -mrSges(3,3) + mrSges(2,2);
t169 = t109 * t190;
t199 = t112 * t113;
t201 = t111 * t112;
t59 = t108 * t113 + t109 * t201;
t262 = qJD(6) * t59 - t113 * t169 - (t108 * t199 + t109 * t111) * qJD(1);
t134 = -t108 * t111 + t109 * t199;
t261 = qJD(6) * t134 + t111 * t169 - (-t108 * t201 + t109 * t113) * qJD(1);
t67 = Ifges(7,4) * t73;
t23 = -Ifges(7,1) * t74 + Ifges(7,5) * t90 + t67;
t95 = Ifges(6,4) * t195;
t260 = -Ifges(6,1) * t194 + Ifges(6,5) * qJD(5) + t113 * t23 + t95;
t204 = qJD(5) * t34;
t12 = -qJDD(3) * t112 + t114 * t36 - t204;
t193 = qJD(2) * t108;
t80 = -t108 * qJ(2) + t109 * t115;
t75 = pkin(3) - t80;
t68 = pkin(7) + t75;
t259 = t112 * t193 + t68 * t190;
t140 = -t11 * t112 - t114 * t12;
t17 = mrSges(7,1) * t69 - mrSges(7,3) * t30;
t18 = -mrSges(7,2) * t69 + mrSges(7,3) * t31;
t257 = -t111 * t17 + t113 * t18;
t254 = mrSges(6,3) - t225;
t220 = Ifges(6,4) * t114;
t147 = t112 * Ifges(6,2) - t220;
t253 = -Ifges(6,6) * qJD(5) / 0.2e1 + t147 * t272 + Ifges(7,5) * t239 + Ifges(7,6) * t274 + Ifges(7,3) * t273;
t251 = qJD(5) * (t112 * t33 - t114 * t34) + t140;
t250 = t222 + t265;
t249 = -t2 * mrSges(7,1) + t1 * mrSges(7,2);
t85 = mrSges(6,1) * t112 + mrSges(6,2) * t114;
t247 = t114 * mrSges(7,3) - t237 * qJ(4) + t224 - t85;
t47 = -mrSges(7,2) * t90 + mrSges(7,3) * t73;
t48 = mrSges(7,1) * t90 + mrSges(7,3) * t74;
t8 = -qJDD(5) * pkin(5) - t12;
t177 = mrSges(6,3) * t195;
t83 = -qJD(5) * mrSges(6,2) + t177;
t246 = m(6) * (t12 + t204) + m(7) * (t10 * t191 - t184 * t9 - t8) - qJD(5) * (t111 * t48 - t113 * t47 - t83) - t223;
t116 = qJD(1) ^ 2;
t245 = Ifges(7,4) * t243 + Ifges(7,2) * t242 + Ifges(7,6) * t241;
t231 = Ifges(7,4) * t74;
t22 = Ifges(7,2) * t73 + Ifges(7,6) * t90 - t231;
t244 = -t22 / 0.2e1;
t238 = t74 / 0.2e1;
t221 = Ifges(6,4) * t112;
t219 = Ifges(7,4) * t111;
t218 = Ifges(7,4) * t113;
t217 = t108 * t37;
t55 = t109 * t196 + t216;
t50 = -t106 + t55;
t215 = t109 * t50;
t212 = t112 * t34;
t211 = t112 * t68;
t203 = t108 * t114;
t202 = t108 * t115;
t200 = t111 * t114;
t198 = t113 * t114;
t197 = t233 * pkin(1) + t232 * qJ(2);
t187 = qJD(6) * t114;
t186 = qJDD(1) * mrSges(3,1);
t185 = qJDD(1) * mrSges(5,2);
t180 = Ifges(7,5) * t30 + Ifges(7,6) * t31 + Ifges(7,3) * t69;
t178 = m(7) * pkin(8) + mrSges(7,3);
t172 = t233 * pkin(2) + t197;
t168 = t112 * t184;
t167 = t108 * t183;
t163 = -t182 / 0.2e1;
t161 = -t70 * pkin(3) + t172;
t89 = qJD(2) * t109 - qJD(4);
t159 = -pkin(1) * t232 + t233 * qJ(2);
t157 = -pkin(5) * t114 - pkin(8) * t112;
t158 = qJD(5) * t157 - qJD(6) * t211 + t89;
t81 = t202 + t205;
t72 = -qJ(4) + t81;
t154 = t37 * t72 + t50 * t89;
t153 = -t10 * t111 - t113 * t9;
t152 = mrSges(6,1) * t114 - mrSges(6,2) * t112;
t151 = -mrSges(7,1) * t113 + mrSges(7,2) * t111;
t150 = t111 * mrSges(7,1) + t113 * mrSges(7,2);
t149 = Ifges(7,1) * t113 - t219;
t148 = Ifges(7,1) * t111 + t218;
t146 = -Ifges(7,2) * t111 + t218;
t145 = Ifges(7,2) * t113 + t219;
t144 = Ifges(6,5) * t112 + Ifges(6,6) * t114;
t143 = Ifges(7,5) * t113 - Ifges(7,6) * t111;
t142 = Ifges(7,5) * t111 + Ifges(7,6) * t113;
t141 = -t108 * t54 + t109 * t55;
t139 = -t111 * t47 - t113 * t48;
t133 = t24 * t150;
t132 = t50 * t152;
t131 = t114 * (Ifges(6,1) * t112 + t220);
t130 = m(7) * pkin(5) - t151;
t127 = t111 * t187 + t112 * t191;
t126 = -t113 * t187 + t168;
t124 = -pkin(2) * t232 + t159;
t123 = t130 * t112;
t121 = -Ifges(7,5) * t114 + t112 * t149;
t120 = -Ifges(7,6) * t114 + t112 * t146;
t119 = -Ifges(7,3) * t114 + t112 * t143;
t51 = -qJ(4) + t128 + t202;
t118 = qJD(6) * t51 + t259;
t58 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t79;
t117 = t58 + t139 * qJD(6) + t222 * qJD(5) + m(6) * (-qJD(5) * t33 + t11) + m(7) * (qJD(5) * t24 + t269) + t257;
t96 = -qJDD(1) * pkin(1) + qJDD(2);
t77 = t157 * qJD(1);
t76 = t85 * qJD(1);
t66 = t150 * t114;
t43 = -mrSges(6,1) * t79 + mrSges(6,2) * t78;
t29 = -t111 * t70 + t199 * t71;
t28 = -t113 * t70 - t201 * t71;
t20 = t111 * t51 + t199 * t68;
t19 = t113 * t51 - t201 * t68;
t16 = t111 * t77 + t113 * t33;
t15 = -t111 * t33 + t113 * t77;
t6 = t30 * Ifges(7,1) + t31 * Ifges(7,4) + t69 * Ifges(7,5);
t4 = -t111 * t118 + t113 * t158;
t3 = t111 * t158 + t113 * t118;
t5 = [(t167 - t45) * mrSges(4,1) + (-mrSges(4,1) * t80 + mrSges(4,2) * t81 - mrSges(5,3) * t72 + Ifges(5,1) + Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1) + (qJD(1) * (Ifges(6,2) * t114 + t221) + t260) * t192 / 0.2e1 - t68 * t270 - t193 * t271 + qJDD(5) * (-Ifges(6,5) * t114 + Ifges(6,6) * t112) + 0.2e1 * t88 * mrSges(3,3) - t89 * t76 - t96 * mrSges(3,1) - t37 * t85 + t72 * t43 - t8 * t66 + t3 * t47 + t4 * t48 + t19 * t17 + t20 * t18 + (m(7) * (-t193 * t24 - t68 * t8) - t143 * t241 - t146 * t242 - t149 * t243 - Ifges(6,1) * t78 + m(6) * t33 * t193 - Ifges(6,4) * t79 / 0.2e1) * t114 + (t23 * t234 + t113 * t22 / 0.2e1 + t148 * t239 + t142 * t273 + t145 * t274) * t187 + (-m(3) * t159 - m(4) * t124 + t263 * t233 + t264 * t232 - t255 * (t71 * pkin(3) + t124) + (pkin(7) * t266 - t150 - t254) * t71 + (-m(7) * (qJ(4) - t230) - t123 + t247) * t70) * g(1) + (-m(3) * t197 - m(4) * t172 - m(5) * t161 - t29 * mrSges(7,1) - t28 * mrSges(7,2) - t264 * t233 + t263 * t232 + t266 * (-pkin(7) * t70 + t161) + t254 * t70 + (-m(7) * (qJ(4) + t156) + t247) * t71) * g(2) + (Ifges(6,4) * t268 + Ifges(6,2) * t267 - Ifges(7,3) * t241 - Ifges(7,6) * t242 - Ifges(7,5) * t243 - t180 / 0.2e1 + t249) * t112 + t24 * (mrSges(7,1) * t126 + mrSges(7,2) * t127) + (t1 * t200 - t10 * t126 - t127 * t9 + t198 * t2) * mrSges(7,3) + (t119 * t273 + t120 * t274 + t121 * t239 - t132 + t144 * qJD(5) / 0.2e1) * qJD(5) - t6 * t198 / 0.2e1 + m(5) * (t193 * t49 + t40 * t75 + t154) + t147 * t267 + t221 * t268 + t168 * t244 + t200 * t245 + (-t9 * mrSges(7,1) + t10 * mrSges(7,2) - t253) * t190 + (-t167 - t40) * mrSges(5,2) + m(7) * (t1 * t20 + t10 * t3 + t19 * t2 + t4 * t9) + (-qJD(1) * t89 - t37) * mrSges(5,3) + m(6) * (t212 * t193 - t251 * t68 + t154) + m(4) * (qJD(2) * t141 + t45 * t80 + t46 * t81) + (t190 * t34 - t192 * t33 - t140) * mrSges(6,3) + t259 * t83 + t250 * t68 * t192 + t131 * t163 + pkin(1) * t186 + t58 * t211 + m(3) * (-pkin(1) * t96 + (t88 + t183) * qJ(2)) - t75 * t185 + (t109 * t183 + t46) * mrSges(4,2); -t186 + t59 * t17 - t134 * t18 + t261 * t48 + t262 * t47 + (-m(3) * qJ(2) - mrSges(3,3)) * t116 + (t224 * qJDD(1) + t225 * t116 + t43) * t108 + (-t112 * t58 - t224 * t116 + t270 + t225 * qJDD(1) + (-t112 * t222 - t114 * t83) * qJD(5)) * t109 + m(6) * (t109 * t251 + t217) + m(5) * (-t109 * t40 + t217) + m(4) * (t108 * t46 + t109 * t45) + m(3) * t96 + (-t1 * t134 + t2 * t59 + (t114 * t8 - t192 * t24) * t109 + t261 * t9 + t262 * t10) * m(7) + ((-t112 * t83 + t271) * t108 + t76 * t109 - m(6) * (t108 * t212 + t203 * t33 + t215) + t203 * t265 - m(5) * (t108 * t49 + t215) - m(4) * t141) * qJD(1) + (g(1) * t232 - g(2) * t233) * (-m(3) - t275); (m(4) + m(5)) * qJDD(3) + t275 * g(3) - t246 * t112 + t117 * t114; -t185 - t116 * mrSges(5,3) + t40 * m(5) + t246 * t114 + t117 * t112 + (-m(7) * t153 + t237 * t50 - t139 - t76) * qJD(1) + t255 * t256; (t22 * t234 - t133) * t195 + (t112 * t178 + t114 * t130 + t152) * t256 + (-pkin(5) * t8 - t10 * t16 - t15 * t9) * m(7) - (Ifges(6,2) * t194 + t260 + t95) * t195 / 0.2e1 + (t149 * t239 + t133) * qJD(6) + t269 * mrSges(7,3) + (t143 * t90 + t146 * t73) * qJD(6) / 0.2e1 + (-t10 * (mrSges(7,2) * t114 - mrSges(7,3) * t201) - t9 * (-mrSges(7,1) * t114 - mrSges(7,3) * t199) + t121 * t238 + t132) * qJD(1) + Ifges(6,5) * t78 + Ifges(6,6) * t79 - t16 * t47 - t15 * t48 - t11 * mrSges(6,2) + t12 * mrSges(6,1) - pkin(5) * t13 + t116 * t131 / 0.2e1 + t6 * t234 + t142 * t241 + t145 * t242 + t148 * t243 + t189 * t244 + t113 * t245 + t253 * t194 + (-t176 - t250) * t34 + t8 * t151 + (-t48 * t188 - t47 * t189 + m(7) * (qJD(6) * t153 + t155) + t257) * pkin(8) + (t114 * t178 - t123 - t85) * g(3) + t144 * t163 + (t90 * t119 + t73 * t120) * t272 + (t177 - t83) * t33 + t23 * t188 / 0.2e1 + Ifges(6,3) * qJDD(5); -t24 * (-mrSges(7,1) * t74 + mrSges(7,2) * t73) + (Ifges(7,1) * t73 + t231) * t238 + t22 * t239 - t90 * (Ifges(7,5) * t73 + Ifges(7,6) * t74) / 0.2e1 - t9 * t47 + t10 * t48 - g(1) * (mrSges(7,1) * t28 - mrSges(7,2) * t29) - g(2) * ((-t113 * t71 + t201 * t70) * mrSges(7,1) + (t111 * t71 + t199 * t70) * mrSges(7,2)) - g(3) * t66 + (-t10 * t74 + t73 * t9) * mrSges(7,3) + t180 - (Ifges(7,2) * t74 + t23 + t67) * t73 / 0.2e1 - t249;];
tau  = t5;
