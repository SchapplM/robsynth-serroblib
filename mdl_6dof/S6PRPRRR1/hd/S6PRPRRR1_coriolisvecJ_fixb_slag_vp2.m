% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:03
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRPRRR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:03:20
% EndTime: 2018-11-23 15:03:25
% DurationCPUTime: 5.40s
% Computational Cost: add. (6222->421), mult. (15278->587), div. (0->0), fcn. (11438->12), ass. (0->212)
t167 = sin(pkin(12));
t168 = sin(pkin(6));
t169 = cos(pkin(12));
t174 = sin(qJ(2));
t178 = cos(qJ(2));
t135 = (t167 * t178 + t169 * t174) * t168;
t129 = qJD(1) * t135;
t173 = sin(qJ(4));
t220 = qJD(4) * t173;
t217 = pkin(4) * t220;
t295 = -t129 + t217;
t223 = qJD(1) * t168;
t213 = t174 * t223;
t153 = t167 * t213;
t212 = t178 * t223;
t131 = t169 * t212 - t153;
t160 = pkin(2) * t167 + pkin(8);
t255 = pkin(9) + t160;
t209 = qJD(4) * t255;
t137 = t173 * t209;
t176 = cos(qJ(5));
t177 = cos(qJ(4));
t172 = sin(qJ(5));
t227 = t172 * t173;
t148 = -t176 * t177 + t227;
t146 = t255 * t173;
t147 = t255 * t177;
t192 = -t146 * t176 - t147 * t172;
t206 = t177 * t209;
t286 = qJD(5) * t192 + t131 * t148 - t137 * t176 - t172 * t206;
t166 = qJD(4) + qJD(5);
t113 = t166 * t148;
t149 = t172 * t177 + t173 * t176;
t114 = t166 * t149;
t294 = pkin(5) * t114 + pkin(10) * t113 + t295;
t293 = t166 * Ifges(6,6) / 0.2e1;
t171 = sin(qJ(6));
t175 = cos(qJ(6));
t152 = qJD(2) * pkin(2) + t212;
t122 = t152 * t167 + t169 * t213;
t118 = qJD(2) * pkin(8) + t122;
t208 = pkin(9) * qJD(2) + t118;
t170 = cos(pkin(6));
t158 = qJD(1) * t170 + qJD(3);
t226 = t173 * t158;
t83 = t177 * t208 + t226;
t235 = t176 * t83;
t154 = t177 * t158;
t194 = t208 * t173;
t82 = t154 - t194;
t79 = qJD(4) * pkin(4) + t82;
t30 = t172 * t79 + t235;
t28 = pkin(10) * t166 + t30;
t121 = t152 * t169 - t153;
t214 = -pkin(4) * t177 - pkin(3);
t106 = qJD(2) * t214 - t121;
t142 = t148 * qJD(2);
t143 = t149 * qJD(2);
t65 = pkin(5) * t142 - pkin(10) * t143 + t106;
t14 = -t171 * t28 + t175 * t65;
t104 = t113 * qJD(2);
t105 = t114 * qJD(2);
t130 = qJD(2) * t135;
t125 = qJD(1) * t130;
t111 = qJD(2) * t217 + t125;
t33 = pkin(5) * t105 + pkin(10) * t104 + t111;
t191 = t167 * t174 - t169 * t178;
t275 = qJD(2) * t168;
t132 = t191 * t275;
t126 = qJD(1) * t132;
t182 = t83 * qJD(4);
t237 = t172 * t83;
t29 = t176 * t79 - t237;
t225 = qJD(4) * t154 - t126 * t177;
t48 = -qJD(4) * t194 + t225;
t8 = qJD(5) * t29 + t126 * t227 - t172 * t182 + t176 * t48;
t2 = qJD(6) * t14 + t171 * t33 + t175 * t8;
t15 = t171 * t65 + t175 * t28;
t3 = -qJD(6) * t15 - t171 * t8 + t175 * t33;
t292 = -t171 * t3 + t175 * t2;
t239 = t166 * Ifges(6,5);
t291 = t106 * mrSges(6,2) + t239 / 0.2e1;
t241 = t143 * Ifges(6,4);
t290 = t293 + t241 / 0.2e1 - t142 * Ifges(6,2) / 0.2e1;
t103 = -t146 * t172 + t147 * t176;
t260 = pkin(2) * t169;
t155 = t214 - t260;
t96 = pkin(5) * t148 - pkin(10) * t149 + t155;
t46 = -t103 * t171 + t175 * t96;
t288 = qJD(6) * t46 + t171 * t294 + t175 * t286;
t47 = t103 * t175 + t171 * t96;
t287 = -qJD(6) * t47 - t171 * t286 + t175 * t294;
t285 = qJD(5) * t103 - t131 * t149 - t137 * t172 + t176 * t206;
t123 = -t143 * t171 + t166 * t175;
t66 = qJD(6) * t123 - t104 * t175;
t124 = t143 * t175 + t166 * t171;
t67 = -qJD(6) * t124 + t104 * t171;
t24 = -mrSges(7,1) * t67 + mrSges(7,2) * t66;
t230 = t126 * t173;
t9 = t172 * t48 - t176 * (-t182 + t230) + t30 * qJD(5);
t283 = m(7) * t9 + t24;
t252 = mrSges(6,3) * t143;
t234 = mrSges(6,1) * t166 + mrSges(7,1) * t123 - mrSges(7,2) * t124 - t252;
t197 = -t14 * t175 - t15 * t171;
t282 = qJD(6) * t197 + t292;
t139 = qJD(6) + t142;
t204 = mrSges(7,1) * t171 + mrSges(7,2) * t175;
t27 = -pkin(5) * t166 - t29;
t189 = t27 * t204;
t199 = Ifges(7,5) * t175 - Ifges(7,6) * t171;
t248 = Ifges(7,4) * t175;
t201 = -Ifges(7,2) * t171 + t248;
t249 = Ifges(7,4) * t171;
t203 = Ifges(7,1) * t175 - t249;
t261 = t175 / 0.2e1;
t262 = -t171 / 0.2e1;
t266 = t124 / 0.2e1;
t250 = Ifges(7,4) * t124;
t58 = Ifges(7,2) * t123 + Ifges(7,6) * t139 + t250;
t120 = Ifges(7,4) * t123;
t59 = t124 * Ifges(7,1) + t139 * Ifges(7,5) + t120;
t281 = t59 * t261 + t58 * t262 + t139 * t199 / 0.2e1 + t203 * t266 + t123 * t201 / 0.2e1 + t189;
t280 = -t106 * mrSges(6,1) - t14 * mrSges(7,1) + t15 * mrSges(7,2) + t290;
t279 = -Ifges(5,1) / 0.2e1;
t278 = -t58 / 0.2e1;
t221 = qJD(2) * t177;
t277 = -Ifges(5,4) * t221 / 0.2e1;
t276 = -t14 * t171 + t15 * t175;
t274 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t66 + Ifges(7,6) * t67;
t272 = t66 / 0.2e1;
t271 = t67 / 0.2e1;
t115 = -t135 * t173 + t170 * t177;
t116 = t135 * t177 + t170 * t173;
t193 = t115 * t176 - t116 * t172;
t270 = t193 * t9;
t269 = t105 / 0.2e1;
t268 = -t123 / 0.2e1;
t267 = -t124 / 0.2e1;
t265 = -t139 / 0.2e1;
t264 = t142 / 0.2e1;
t263 = -t143 / 0.2e1;
t259 = t192 * t9;
t258 = t148 * t9;
t253 = mrSges(6,3) * t142;
t251 = Ifges(5,4) * t173;
t138 = Ifges(6,4) * t142;
t247 = t123 * Ifges(7,6);
t246 = t124 * Ifges(7,5);
t245 = t139 * Ifges(7,3);
t242 = t143 * Ifges(6,1);
t87 = -t118 * t173 + t154;
t236 = t173 * t87;
t233 = Ifges(5,5) * qJD(4);
t232 = Ifges(5,6) * qJD(4);
t134 = t191 * t168;
t231 = t125 * t134;
t229 = t142 * t171;
t228 = t142 * t175;
t108 = mrSges(6,1) * t142 + mrSges(6,2) * t143;
t224 = (-mrSges(5,1) * t177 + mrSges(5,2) * t173) * qJD(2) + t108;
t222 = qJD(2) * t173;
t219 = qJD(6) * t171;
t218 = qJD(6) * t175;
t211 = t233 / 0.2e1;
t210 = -t232 / 0.2e1;
t109 = pkin(5) * t143 + pkin(10) * t142;
t205 = mrSges(7,1) * t175 - mrSges(7,2) * t171;
t202 = Ifges(7,1) * t171 + t248;
t200 = Ifges(7,2) * t175 + t249;
t198 = Ifges(7,5) * t171 + Ifges(7,6) * t175;
t31 = mrSges(7,1) * t105 - mrSges(7,3) * t66;
t32 = -mrSges(7,2) * t105 + mrSges(7,3) * t67;
t195 = -t171 * t31 + t175 * t32;
t69 = t115 * t172 + t116 * t176;
t38 = t134 * t175 - t171 * t69;
t39 = t134 * t171 + t175 * t69;
t88 = t118 * t177 + t226;
t127 = -mrSges(6,2) * t166 - t253;
t84 = -mrSges(7,2) * t139 + mrSges(7,3) * t123;
t85 = mrSges(7,1) * t139 - mrSges(7,3) * t124;
t190 = -t171 * t85 + t175 * t84 + t127;
t117 = -qJD(2) * pkin(3) - t121;
t188 = t87 * mrSges(5,3) + t222 * t279 + t277 - t233 / 0.2e1 - t117 * mrSges(5,2);
t187 = t88 * mrSges(5,3) + t232 / 0.2e1 + (t177 * Ifges(5,2) + t251) * qJD(2) / 0.2e1 - t117 * mrSges(5,1);
t181 = m(7) * (-t14 * t218 - t15 * t219 + t292) - t84 * t219 - t85 * t218 + t195;
t18 = Ifges(7,4) * t66 + Ifges(7,2) * t67 + Ifges(7,6) * t105;
t19 = Ifges(7,1) * t66 + Ifges(7,4) * t67 + Ifges(7,5) * t105;
t57 = t245 + t246 + t247;
t94 = -t138 + t239 + t242;
t180 = t59 * t228 / 0.2e1 + t229 * t278 + t171 * t19 / 0.2e1 + t198 * t269 + t200 * t271 + t202 * t272 + t18 * t261 - t29 * t253 - t8 * mrSges(6,2) - Ifges(6,5) * t104 - Ifges(6,6) * t105 + (-t205 - mrSges(6,1)) * t9 + (-t138 + t94) * t264 + (-t241 + t57) * t263 + (-Ifges(6,1) * t263 - t199 * t265 - t201 * t268 - t203 * t267 + t189 + t291) * t142 + (Ifges(7,5) * t267 - Ifges(6,2) * t264 + Ifges(7,6) * t268 + Ifges(7,3) * t265 + t280 + t293) * t143 + (-t14 * t228 - t15 * t229 + t282) * mrSges(7,3) + t281 * qJD(6);
t161 = -pkin(3) - t260;
t157 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t221;
t156 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t222;
t144 = (mrSges(5,1) * t173 + mrSges(5,2) * t177) * qJD(4) * qJD(2);
t99 = Ifges(7,3) * t105;
t92 = pkin(4) * t222 + t109;
t72 = -qJD(4) * t116 + t132 * t173;
t71 = qJD(4) * t115 - t132 * t177;
t61 = -qJD(4) * t88 + t230;
t60 = -t118 * t220 + t225;
t56 = mrSges(6,1) * t105 - mrSges(6,2) * t104;
t35 = t176 * t82 - t237;
t34 = t172 * t82 + t235;
t23 = t109 * t171 + t175 * t29;
t22 = t109 * t175 - t171 * t29;
t21 = t171 * t92 + t175 * t35;
t20 = -t171 * t35 + t175 * t92;
t13 = qJD(5) * t69 + t172 * t71 - t176 * t72;
t12 = qJD(5) * t193 + t172 * t72 + t176 * t71;
t5 = -qJD(6) * t39 - t12 * t171 + t130 * t175;
t4 = qJD(6) * t38 + t12 * t175 + t130 * t171;
t1 = [t12 * t127 + t72 * t156 + t71 * t157 - t193 * t24 + t38 * t31 + t39 * t32 + t4 * t84 + t5 * t85 + (t144 + t56) * t134 + t224 * t130 - t234 * t13 + (t104 * t193 - t105 * t69) * mrSges(6,3) + m(7) * (t13 * t27 + t14 * t5 + t15 * t4 + t2 * t39 + t3 * t38 - t270) + m(6) * (t106 * t130 + t111 * t134 + t12 * t30 - t13 * t29 + t69 * t8 - t270) + m(5) * (t115 * t61 + t116 * t60 + t117 * t130 + t71 * t88 + t72 * t87 + t231) + m(4) * (-t121 * t130 - t122 * t132 - t126 * t135 + t231) + (-t130 * mrSges(4,1) + t132 * mrSges(4,2) + (-t115 * t177 - t116 * t173) * qJD(4) * mrSges(5,3) + (-mrSges(3,1) * t174 - mrSges(3,2) * t178) * t275) * qJD(2); (t246 / 0.2e1 + t245 / 0.2e1 + t247 / 0.2e1 + t57 / 0.2e1 - t280 - t290) * t114 + t286 * t127 + (-t103 * t105 + t104 * t192 + t113 * t29 - t114 * t30) * mrSges(6,3) - (t242 / 0.2e1 - t138 / 0.2e1 + t94 / 0.2e1 + t197 * mrSges(7,3) + t281 + t291) * t113 - m(5) * (t117 * t129 - t131 * t236) + (t125 * mrSges(5,2) - t61 * mrSges(5,3) + t131 * t156 + (pkin(4) * t108 - t160 * t157 - 0.3e1 / 0.2e1 * Ifges(5,4) * t222 + t210 - t187) * qJD(4)) * t173 + (mrSges(4,1) * qJD(2) - t224) * t129 + (-t125 * mrSges(5,1) + (t211 + (-m(5) * t87 - t156) * t160 + (0.3e1 / 0.2e1 * Ifges(5,4) * t177 + (0.3e1 / 0.2e1 * Ifges(5,1) - 0.3e1 / 0.2e1 * Ifges(5,2)) * t173) * qJD(2) - t188) * qJD(4) + (m(5) * t160 + mrSges(5,3)) * t60 + (-m(5) * t88 - t157) * t131) * t177 + t287 * t85 + t288 * t84 + (t203 * t272 + t201 * t271 + t199 * t269 - Ifges(6,1) * t104 - Ifges(6,4) * t105 + t111 * mrSges(6,2) + t18 * t262 + t19 * t261 + (mrSges(6,3) + t204) * t9 + (-t171 * t2 - t175 * t3) * mrSges(7,3) + (-mrSges(7,3) * t276 + t175 * t278 + t198 * t265 + t200 * t268 + t202 * t267 + t27 * t205 + t59 * t262) * qJD(6)) * t149 + (t99 / 0.2e1 + Ifges(6,4) * t104 + t111 * mrSges(6,1) - t8 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2)) * t105 + t274) * t148 + (qJD(2) * t131 + t126) * mrSges(4,2) + t155 * t56 + t161 * t144 + m(5) * (t125 * t161 + (-t173 * t61 - t220 * t88) * t160) + t46 * t31 + t47 * t32 - t192 * t24 - t125 * mrSges(4,1) - t234 * t285 + (t14 * t287 + t15 * t288 + t2 * t47 + t27 * t285 + t3 * t46 - t259) * m(7) + (t103 * t8 + t106 * t295 + t111 * t155 - t285 * t29 + t286 * t30 - t259) * m(6) + ((-t125 * t169 - t126 * t167) * pkin(2) + t121 * t129 - t122 * t131) * m(4); (-t104 * mrSges(6,3) + t24) * t148 - t234 * t114 - t190 * t113 + (-t173 * t156 + t177 * t157 + (-t173 ^ 2 - t177 ^ 2) * qJD(2) * mrSges(5,3)) * qJD(4) + (-t105 * mrSges(6,3) + (-t171 * t84 - t175 * t85) * qJD(6) + t195) * t149 + m(7) * (-t113 * t276 + t114 * t27 + t149 * t282 + t258) + m(5) * (t173 * t60 + t177 * t61 + (t177 * t88 - t236) * qJD(4)) + m(6) * (-t113 * t30 - t114 * t29 + t149 * t8 + t258); -m(7) * (t14 * t20 + t15 * t21 + t27 * t34) + t181 * (pkin(4) * t172 + pkin(10)) + ((t211 + t277 + t188) * t177 + (t210 + (t251 / 0.2e1 + (t279 + Ifges(5,2) / 0.2e1) * t177) * qJD(2) + (-m(6) * t106 - t108) * pkin(4) + t187) * t173) * qJD(2) + (m(6) * (t172 * t8 - t176 * t9) + (t104 * t176 - t105 * t172) * mrSges(6,3) + ((-m(6) * t29 + m(7) * t27 - t234) * t172 + (m(6) * t30 + m(7) * t276 + t190) * t176) * qJD(5)) * pkin(4) + t180 + t234 * t34 - m(6) * (-t29 * t34 + t30 * t35) + t88 * t156 - t87 * t157 + t30 * t252 - t60 * mrSges(5,2) + t61 * mrSges(5,1) - t21 * t84 - t20 * t85 - t35 * t127 + t283 * (-pkin(4) * t176 - pkin(5)); t181 * pkin(10) + t180 + (t234 + t252) * t30 - m(7) * (t14 * t22 + t15 * t23 + t27 * t30) - t23 * t84 - t22 * t85 - t29 * t127 - t283 * pkin(5); t99 - t27 * (mrSges(7,1) * t124 + mrSges(7,2) * t123) + (Ifges(7,1) * t123 - t250) * t267 + t58 * t266 + (Ifges(7,5) * t123 - Ifges(7,6) * t124) * t265 - t14 * t84 + t15 * t85 + (t123 * t14 + t124 * t15) * mrSges(7,3) + (-Ifges(7,2) * t124 + t120 + t59) * t268 + t274;];
tauc  = t1(:);
