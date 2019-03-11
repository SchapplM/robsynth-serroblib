% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRRR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:17:57
% EndTime: 2019-03-09 02:18:06
% DurationCPUTime: 4.42s
% Computational Cost: add. (8951->401), mult. (22079->541), div. (0->0), fcn. (16678->10), ass. (0->196)
t148 = qJD(4) + qJD(5);
t153 = sin(qJ(6));
t156 = cos(qJ(6));
t149 = sin(pkin(11));
t155 = sin(qJ(4));
t151 = cos(pkin(11));
t157 = cos(qJ(4));
t209 = t151 * t157;
t136 = -t149 * t155 + t209;
t129 = t136 * qJD(1);
t137 = t149 * t157 + t151 * t155;
t130 = t137 * qJD(1);
t154 = sin(qJ(5));
t251 = cos(qJ(5));
t171 = t154 * t129 + t130 * t251;
t89 = t148 * t156 - t153 * t171;
t284 = t89 / 0.2e1;
t90 = t148 * t153 + t156 * t171;
t258 = t90 / 0.2e1;
t191 = t251 * t129 - t154 * t130;
t271 = qJD(6) - t191;
t283 = t271 / 0.2e1;
t282 = -t148 * Ifges(6,6) / 0.2e1 - Ifges(6,4) * t171 / 0.2e1;
t281 = t171 * Ifges(6,1) / 0.2e1;
t235 = t191 * Ifges(6,2);
t280 = -t235 / 0.2e1;
t190 = qJD(1) * (t149 ^ 2 + t151 ^ 2);
t279 = mrSges(4,3) * t190;
t131 = t136 * qJD(4);
t121 = qJD(1) * t131;
t132 = t137 * qJD(4);
t122 = qJD(1) * t132;
t64 = qJD(5) * t191 + t251 * t121 - t154 * t122;
t45 = qJD(6) * t89 + t156 * t64;
t65 = qJD(5) * t171 + t154 * t121 + t122 * t251;
t23 = mrSges(7,1) * t65 - mrSges(7,3) * t45;
t46 = -qJD(6) * t90 - t153 * t64;
t24 = -mrSges(7,2) * t65 + mrSges(7,3) * t46;
t176 = -t153 * t23 + t156 * t24;
t202 = qJD(6) * t156;
t203 = qJD(6) * t153;
t56 = -mrSges(7,2) * t271 + mrSges(7,3) * t89;
t57 = mrSges(7,1) * t271 - mrSges(7,3) * t90;
t278 = -t57 * t202 - t56 * t203 + t176;
t277 = t148 * Ifges(6,5) / 0.2e1 + Ifges(6,4) * t191 / 0.2e1;
t276 = t281 + t277;
t180 = Ifges(7,5) * t156 - Ifges(7,6) * t153;
t228 = Ifges(7,4) * t156;
t182 = -Ifges(7,2) * t153 + t228;
t229 = Ifges(7,4) * t153;
t184 = Ifges(7,1) * t156 - t229;
t252 = t156 / 0.2e1;
t253 = -t153 / 0.2e1;
t250 = Ifges(7,4) * t90;
t41 = Ifges(7,2) * t89 + Ifges(7,6) * t271 + t250;
t88 = Ifges(7,4) * t89;
t42 = Ifges(7,1) * t90 + Ifges(7,5) * t271 + t88;
t275 = t180 * t283 + t182 * t284 + t184 * t258 + t42 * t252 + t41 * t253;
t274 = Ifges(7,5) * t258 + Ifges(7,6) * t284 + Ifges(7,3) * t283 + t282;
t233 = -mrSges(6,1) * t148 - mrSges(7,1) * t89 + mrSges(7,2) * t90 + mrSges(6,3) * t171;
t141 = sin(pkin(10)) * pkin(1) + qJ(3);
t234 = pkin(7) + t141;
t133 = t234 * t149;
t134 = t234 * t151;
t102 = -t155 * t133 + t157 * t134;
t138 = t141 * qJD(1);
t145 = t151 * qJD(2);
t226 = pkin(7) * qJD(1);
t109 = t145 + (-t138 - t226) * t149;
t117 = t149 * qJD(2) + t151 * t138;
t110 = t151 * t226 + t117;
t81 = t109 * t155 + t110 * t157;
t74 = pkin(8) * t129 + t81;
t199 = t251 * t74;
t80 = t157 * t109 - t110 * t155;
t73 = -pkin(8) * t130 + t80;
t71 = qJD(4) * pkin(4) + t73;
t33 = t154 * t71 + t199;
t31 = t148 * pkin(9) + t33;
t174 = -cos(pkin(10)) * pkin(1) - pkin(3) * t151 - pkin(2);
t127 = qJD(1) * t174 + qJD(3);
t103 = -pkin(4) * t129 + t127;
t47 = -pkin(5) * t191 - pkin(9) * t171 + t103;
t14 = -t153 * t31 + t156 * t47;
t15 = t153 * t47 + t156 * t31;
t272 = -t14 * t153 + t15 * t156;
t169 = t137 * qJD(3);
t69 = -qJD(1) * t169 - qJD(4) * t81;
t161 = -pkin(8) * t121 + t69;
t215 = t154 * t74;
t32 = t251 * t71 - t215;
t140 = qJD(3) * t209;
t204 = qJD(4) * t157;
t205 = qJD(3) * t149;
t68 = t109 * t204 + qJD(1) * t140 + (-qJD(1) * t205 - qJD(4) * t110) * t155;
t53 = -pkin(8) * t122 + t68;
t12 = qJD(5) * t32 + t154 * t161 + t251 * t53;
t249 = pkin(4) * t122;
t27 = pkin(5) * t65 - pkin(9) * t64 + t249;
t2 = qJD(6) * t14 + t12 * t156 + t153 * t27;
t3 = -qJD(6) * t15 - t12 * t153 + t156 * t27;
t270 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t45 + Ifges(7,6) * t46;
t91 = -mrSges(6,2) * t148 + mrSges(6,3) * t191;
t269 = -t153 * t57 + t156 * t56 + t91;
t244 = t153 * t3;
t268 = -t14 * t202 - t15 * t203 - t244;
t67 = pkin(5) * t171 - pkin(9) * t191;
t178 = t14 * t156 + t15 * t153;
t185 = mrSges(7,1) * t153 + mrSges(7,2) * t156;
t30 = -t148 * pkin(5) - t32;
t172 = t30 * t185;
t266 = t103 * mrSges(6,2) + t172 + t276;
t267 = -t178 * mrSges(7,3) + t266 + t275 + t277;
t265 = t103 * mrSges(6,1) + t14 * mrSges(7,1) - t15 * mrSges(7,2) + t274 + t280;
t264 = m(6) * pkin(4);
t263 = t45 / 0.2e1;
t262 = t46 / 0.2e1;
t261 = t65 / 0.2e1;
t260 = -t89 / 0.2e1;
t259 = -t90 / 0.2e1;
t257 = -t271 / 0.2e1;
t255 = t131 / 0.2e1;
t254 = -t132 / 0.2e1;
t248 = pkin(4) * t130;
t247 = pkin(4) * t132;
t246 = pkin(4) * t154;
t13 = qJD(5) * t33 + t154 * t53 - t161 * t251;
t101 = -t157 * t133 - t134 * t155;
t86 = -pkin(8) * t137 + t101;
t87 = pkin(8) * t136 + t102;
t173 = -t154 * t87 + t251 * t86;
t245 = t13 * t173;
t243 = t156 * t2;
t242 = t32 * mrSges(6,3);
t241 = t33 * mrSges(6,3);
t240 = t64 * mrSges(6,3);
t239 = t65 * mrSges(6,3);
t232 = mrSges(5,3) * t129;
t231 = mrSges(5,3) * t130;
t227 = pkin(4) * qJD(5);
t170 = t136 * t251 - t154 * t137;
t224 = t170 * t13;
t223 = t130 * Ifges(5,4);
t216 = t153 * t191;
t212 = t156 * t191;
t211 = t136 * t121;
t210 = t137 * t122;
t201 = t251 * pkin(4);
t192 = t65 * mrSges(6,1) + t64 * mrSges(6,2);
t187 = -t153 * t2 - t156 * t3;
t186 = mrSges(7,1) * t156 - mrSges(7,2) * t153;
t183 = Ifges(7,1) * t153 + t228;
t181 = Ifges(7,2) * t156 + t229;
t179 = Ifges(7,5) * t153 + Ifges(7,6) * t156;
t49 = t154 * t86 + t251 * t87;
t105 = t154 * t136 + t137 * t251;
t108 = -pkin(4) * t136 + t174;
t55 = -pkin(5) * t170 - pkin(9) * t105 + t108;
t26 = t153 * t55 + t156 * t49;
t25 = -t153 * t49 + t156 * t55;
t175 = -(-t138 * t149 + t145) * t149 + t117 * t151;
t165 = -qJD(6) * t178 - t244;
t82 = -t133 * t204 + t140 + (-qJD(4) * t134 - t205) * t155;
t164 = t165 + t243;
t83 = -qJD(4) * t102 - t169;
t10 = t45 * Ifges(7,1) + t46 * Ifges(7,4) + t65 * Ifges(7,5);
t9 = t45 * Ifges(7,4) + t46 * Ifges(7,2) + t65 * Ifges(7,6);
t163 = -t12 * mrSges(6,2) + mrSges(7,3) * t243 + t183 * t263 + t181 * t262 + t179 * t261 + t153 * t10 / 0.2e1 - Ifges(6,6) * t65 + Ifges(6,5) * t64 + t9 * t252 + (-mrSges(6,1) - t186) * t13 + (t172 + t275) * qJD(6);
t162 = -pkin(8) * t131 + t83;
t160 = t265 + t274;
t143 = -t201 - pkin(5);
t126 = Ifges(5,4) * t129;
t114 = t121 * mrSges(5,2);
t113 = qJD(4) * mrSges(5,1) - t231;
t112 = -qJD(4) * mrSges(5,2) + t232;
t96 = t130 * Ifges(5,1) + Ifges(5,5) * qJD(4) + t126;
t95 = t129 * Ifges(5,2) + Ifges(5,6) * qJD(4) + t223;
t77 = -pkin(8) * t132 + t82;
t76 = qJD(5) * t105 + t154 * t131 + t132 * t251;
t75 = qJD(5) * t170 + t131 * t251 - t154 * t132;
t66 = -mrSges(6,1) * t191 + mrSges(6,2) * t171;
t61 = Ifges(7,3) * t65;
t52 = t248 + t67;
t35 = t251 * t73 - t215;
t34 = t154 * t73 + t199;
t29 = pkin(5) * t76 - pkin(9) * t75 + t247;
t22 = t153 * t67 + t156 * t32;
t21 = -t153 * t32 + t156 * t67;
t20 = t153 * t52 + t156 * t35;
t19 = -t153 * t35 + t156 * t52;
t18 = qJD(5) * t49 + t154 * t77 - t162 * t251;
t17 = qJD(5) * t173 + t154 * t162 + t251 * t77;
t16 = -mrSges(7,1) * t46 + mrSges(7,2) * t45;
t5 = -qJD(6) * t26 - t153 * t17 + t156 * t29;
t4 = qJD(6) * t25 + t153 * t29 + t156 * t17;
t1 = [(t137 * t121 + t130 * t255) * Ifges(5,1) + t66 * t247 + t233 * t18 + t108 * t192 + (t10 * t252 + t9 * t253 + mrSges(6,2) * t249 + Ifges(6,1) * t64 - Ifges(6,4) * t65 + t184 * t263 + t182 * t262 + t180 * t261 + (mrSges(6,3) + t185) * t13 + t187 * mrSges(7,3) + (t42 * t253 - t156 * t41 / 0.2e1 + t30 * t186 + t181 * t260 + t183 * t259 + t179 * t257 - t272 * mrSges(7,3)) * qJD(6)) * t105 + m(7) * (t14 * t5 + t15 * t4 + t18 * t30 + t2 * t26 + t25 * t3 - t245) + qJD(4) * (Ifges(5,5) * t131 - Ifges(5,6) * t132) / 0.2e1 + t127 * (mrSges(5,1) * t132 + mrSges(5,2) * t131) + t95 * t254 + t96 * t255 + (t280 + t160) * t76 + (t281 + t267) * t75 + (m(4) * (t141 * t190 + t175) + 0.2e1 * t279) * qJD(3) + m(6) * (t12 * t49 - t245 + t17 * t33 - t18 * t32 + (t103 * t132 + t108 * t122) * pkin(4)) + (-t101 * t121 - t102 * t122 - t131 * t80 - t132 * t81 + t136 * t68 - t137 * t69) * mrSges(5,3) + t174 * (t122 * mrSges(5,1) + t114) + (t129 * t255 + t130 * t254 - t210 + t211) * Ifges(5,4) - (mrSges(6,1) * t249 - t12 * mrSges(6,3) + t61 / 0.2e1 - Ifges(6,4) * t64 + (Ifges(6,2) + Ifges(7,3) / 0.2e1) * t65 + t270) * t170 + (-t173 * t64 - t32 * t75 - t33 * t76 - t49 * t65) * mrSges(6,3) - t173 * t16 + (-t136 * t122 + t129 * t254) * Ifges(5,2) + m(5) * (t101 * t69 + t102 * t68 + t80 * t83 + t81 * t82) + t25 * t23 + t26 * t24 + t4 * t56 + t5 * t57 + t17 * t91 + t82 * t112 + t83 * t113; t131 * t112 - t132 * t113 + t233 * t76 - (t16 + t240) * t170 + (-t210 - t211) * mrSges(5,3) + t269 * t75 + (-t239 + (-t153 * t56 - t156 * t57) * qJD(6) + t176) * t105 + m(7) * (t105 * t164 + t272 * t75 + t30 * t76 - t224) + m(5) * (t131 * t81 - t132 * t80 + t136 * t69 + t137 * t68) + m(6) * (t105 * t12 - t32 * t76 + t33 * t75 - t224); -t129 * t112 + t130 * t113 - t191 * t91 + t114 - t233 * t171 - (-mrSges(5,1) - t264) * t122 + (t271 * t56 + t23) * t156 + (-t271 * t57 + t24) * t153 - m(5) * (t129 * t81 - t130 * t80) - m(6) * (-t171 * t32 + t191 * t33) + t192 + (-m(4) * t175 - t279) * qJD(1) + (-t30 * t171 + t271 * t272 - t187) * m(7); (m(6) * t32 - t233) * t34 + (t231 + t113) * t81 - t239 * t246 - t66 * t248 - t130 * (Ifges(5,1) * t129 - t223) / 0.2e1 + (-t14 * t19 - t15 * t20 - t30 * t34 + t13 * t143 + (t154 * t30 + t251 * t272) * t227) * m(7) + (t14 * t212 + t15 * t216 + t268) * mrSges(7,3) + t41 * t216 / 0.2e1 + (t232 - t112) * t80 + (qJD(5) * t269 - t240) * t201 + (t180 * t257 + t182 * t260 + t184 * t259 + t242 - t266 - t276) * t191 - t42 * t212 / 0.2e1 + (Ifges(7,3) * t257 + Ifges(7,5) * t259 + Ifges(7,6) * t260 + t241 + t235 / 0.2e1 - t265 - t282) * t171 - m(6) * (t103 * t248 + t33 * t35) - (-Ifges(5,2) * t130 + t126 + t96) * t129 / 0.2e1 + (m(7) * t164 + t278) * (pkin(9) + t246) + t233 * t154 * t227 + (-t251 * t13 + t12 * t154 + (-t154 * t32 + t251 * t33) * qJD(5)) * t264 + t163 - t20 * t56 - t19 * t57 - t68 * mrSges(5,2) + t69 * mrSges(5,1) - t35 * t91 + Ifges(5,5) * t121 - Ifges(5,6) * t122 + t130 * t95 / 0.2e1 - qJD(4) * (Ifges(5,5) * t129 - Ifges(5,6) * t130) / 0.2e1 - t127 * (mrSges(5,1) * t130 + mrSges(5,2) * t129) + t143 * t16; t165 * mrSges(7,3) - t233 * t33 + t163 + (t242 + (Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t171 - t267) * t191 + (-t160 + t241) * t171 + t278 * pkin(9) - pkin(5) * t16 - t22 * t56 - t21 * t57 - t32 * t91 + (-pkin(5) * t13 - t14 * t21 - t15 * t22 - t30 * t33 + (t243 + t268) * pkin(9)) * m(7); t61 - t30 * (mrSges(7,1) * t90 + mrSges(7,2) * t89) + (Ifges(7,1) * t89 - t250) * t259 + t41 * t258 + (Ifges(7,5) * t89 - Ifges(7,6) * t90) * t257 - t14 * t56 + t15 * t57 + (t14 * t89 + t15 * t90) * mrSges(7,3) + (-Ifges(7,2) * t90 + t42 + t88) * t260 + t270;];
tauc  = t1(:);
