% Calculate time derivative of joint inertia matrix for
% S6RRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:57:36
% EndTime: 2019-03-09 21:57:45
% DurationCPUTime: 4.09s
% Computational Cost: add. (13281->450), mult. (28951->674), div. (0->0), fcn. (28750->10), ass. (0->207)
t276 = 2 * mrSges(5,3);
t195 = sin(qJ(3));
t196 = sin(qJ(2));
t199 = cos(qJ(3));
t200 = cos(qJ(2));
t162 = t195 * t200 + t196 * t199;
t262 = -pkin(8) - pkin(7);
t173 = t262 * t196;
t174 = t262 * t200;
t139 = t173 * t199 + t174 * t195;
t112 = -pkin(9) * t162 + t139;
t165 = t195 * t173;
t140 = -t174 * t199 + t165;
t228 = t199 * t200;
t161 = -t195 * t196 + t228;
t113 = pkin(9) * t161 + t140;
t194 = sin(qJ(4));
t198 = cos(qJ(4));
t275 = t112 * t198 - t113 * t194;
t191 = sin(pkin(11));
t192 = cos(pkin(11));
t127 = t161 * t194 + t162 * t198;
t186 = -pkin(2) * t200 - pkin(1);
t142 = -pkin(3) * t161 + t186;
t208 = t161 * t198 - t162 * t194;
t77 = -pkin(4) * t208 - qJ(5) * t127 + t142;
t79 = t112 * t194 + t113 * t198;
t49 = -t191 * t79 + t192 * t77;
t50 = t191 * t77 + t192 * t79;
t209 = -t191 * t49 + t192 * t50;
t193 = sin(qJ(6));
t197 = cos(qJ(6));
t207 = t191 * t193 - t192 * t197;
t150 = t207 * qJD(6);
t274 = qJD(2) + qJD(3);
t170 = (-pkin(10) - qJ(5)) * t191;
t188 = t192 * pkin(10);
t172 = qJ(5) * t192 + t188;
t136 = t170 * t193 + t172 * t197;
t160 = t191 * t197 + t192 * t193;
t151 = t160 * qJD(6);
t254 = mrSges(7,3) * t151;
t107 = t136 * t254;
t119 = t151 * mrSges(7,1) - mrSges(7,2) * t150;
t259 = pkin(5) * t192;
t181 = -pkin(4) - t259;
t108 = t181 * t119;
t189 = t191 ^ 2;
t240 = mrSges(6,3) * qJD(5);
t182 = t189 * t240;
t190 = t192 ^ 2;
t183 = t190 * t240;
t135 = t170 * t197 - t172 * t193;
t105 = -qJD(5) * t207 + qJD(6) * t135;
t253 = mrSges(7,3) * t207;
t93 = t105 * t253;
t273 = -t107 + t108 + t182 + t183 - t93;
t260 = pkin(3) * t198;
t169 = t181 - t260;
t103 = t169 * t119;
t130 = mrSges(7,1) * t207 + mrSges(7,2) * t160;
t223 = qJD(4) * t194;
t218 = pkin(3) * t223;
t122 = t130 * t218;
t171 = -mrSges(6,1) * t192 + mrSges(6,2) * t191;
t153 = t171 * t218;
t222 = qJD(4) * t198;
t217 = pkin(3) * t222;
t177 = qJD(5) + t217;
t234 = t177 * t189;
t163 = mrSges(6,3) * t234;
t233 = t177 * t190;
t164 = mrSges(6,3) * t233;
t180 = pkin(3) * t194 + qJ(5);
t156 = (-pkin(10) - t180) * t191;
t232 = t180 * t192;
t157 = t188 + t232;
t116 = t156 * t197 - t157 * t193;
t88 = qJD(6) * t116 - t177 * t207;
t80 = t88 * t253;
t117 = t156 * t193 + t157 * t197;
t98 = t117 * t254;
t272 = t103 + t122 + t153 + t163 + t164 - t80 - t98;
t271 = 2 * m(5);
t270 = 2 * m(6);
t269 = 2 * m(7);
t133 = t274 * t161;
t257 = t133 * pkin(9);
t134 = t274 * t162;
t224 = qJD(3) * t199;
t225 = qJD(3) * t195;
t96 = qJD(2) * t162 * t262 + t173 * t224 + t174 * t225;
t75 = -pkin(9) * t134 + t96;
t204 = (t228 * t262 - t165) * qJD(2);
t97 = -qJD(3) * t140 + t204;
t39 = qJD(4) * t79 + t194 * t75 - t198 * (t97 - t257);
t268 = 0.2e1 * t39;
t114 = pkin(2) * qJD(2) * t196 + pkin(3) * t134;
t267 = 0.2e1 * t114;
t266 = 0.2e1 * t186;
t265 = m(5) / 0.2e1;
t264 = m(6) / 0.2e1;
t261 = m(7) * t194;
t70 = qJD(4) * t208 + t133 * t198 - t134 * t194;
t71 = qJD(4) * t127 + t133 * t194 + t134 * t198;
t36 = pkin(4) * t71 - qJ(5) * t70 - qJD(5) * t127 + t114;
t38 = t198 * t75 + t275 * qJD(4) + (-t173 * t225 + t174 * t224 + t204 - t257) * t194;
t12 = -t191 * t38 + t192 * t36;
t258 = t12 * mrSges(6,3);
t256 = t39 * t275;
t13 = t191 * t36 + t192 * t38;
t242 = t192 * t70;
t244 = t191 * t70;
t42 = mrSges(6,1) * t244 + mrSges(6,2) * t242;
t82 = t160 * t127;
t83 = t207 * t127;
t52 = mrSges(7,1) * t82 - mrSges(7,2) * t83;
t84 = (mrSges(6,1) * t191 + mrSges(6,2) * t192) * t127;
t255 = t84 + t52;
t252 = Ifges(6,4) * t191;
t251 = Ifges(6,4) * t192;
t250 = Ifges(6,2) * t191;
t249 = t12 * t191;
t178 = t194 * t195 * pkin(2);
t185 = pkin(2) * t199 + pkin(3);
t123 = t185 * t222 + t198 * pkin(2) * t224 + (-qJD(3) - qJD(4)) * t178;
t248 = t123 * mrSges(5,2);
t231 = t195 * t198;
t124 = t185 * t223 + (t195 * t222 + (t194 * t199 + t231) * qJD(3)) * pkin(2);
t247 = t124 * t275;
t246 = t13 * t192;
t241 = t194 * t275;
t118 = qJD(5) + t123;
t238 = t118 * t189;
t237 = t118 * t190;
t236 = t127 * t191;
t235 = t127 * t192;
t227 = -Ifges(7,5) * t150 - Ifges(7,6) * t151;
t149 = pkin(2) * t231 + t185 * t194;
t226 = t189 + t190;
t221 = 0.2e1 * mrSges(7,3);
t220 = 0.2e1 * t200;
t30 = -t127 * t151 - t207 * t70;
t31 = t127 * t150 - t160 * t70;
t219 = Ifges(7,5) * t30 + Ifges(7,6) * t31 + Ifges(7,3) * t71;
t10 = -t31 * mrSges(7,1) + mrSges(7,2) * t30;
t215 = t226 * t118;
t214 = t226 * t177;
t148 = t185 * t198 - t178;
t213 = qJD(5) * t226;
t120 = -Ifges(7,4) * t150 - Ifges(7,2) * t151;
t121 = -Ifges(7,1) * t150 - Ifges(7,4) * t151;
t131 = Ifges(7,4) * t160 - Ifges(7,2) * t207;
t132 = Ifges(7,1) * t160 - Ifges(7,4) * t207;
t212 = -t120 * t207 + t121 * t160 - t131 * t151 - t132 * t150;
t147 = -pkin(4) - t148;
t211 = Ifges(6,5) * t192 - Ifges(6,6) * t191;
t210 = t246 - t249;
t40 = -pkin(5) * t208 - pkin(10) * t235 + t49;
t41 = -pkin(10) * t236 + t50;
t15 = -t193 * t41 + t197 * t40;
t16 = t193 * t40 + t197 * t41;
t143 = qJ(5) + t149;
t137 = (-pkin(10) - t143) * t191;
t138 = t143 * t192 + t188;
t91 = t137 * t197 - t138 * t193;
t92 = t137 * t193 + t138 * t197;
t206 = (-mrSges(4,1) * t195 - mrSges(4,2) * t199) * qJD(3) * pkin(2);
t205 = (-mrSges(5,1) * t194 - mrSges(5,2) * t198) * qJD(4) * pkin(3);
t104 = t124 * t171;
t110 = mrSges(6,3) * t238;
t111 = mrSges(6,3) * t237;
t115 = t124 * mrSges(5,1);
t54 = qJD(6) * t91 - t118 * t207;
t51 = t54 * t253;
t81 = t92 * t254;
t90 = t124 * t130;
t141 = t147 - t259;
t95 = t141 * t119;
t203 = t104 + t110 + t111 - t51 - t81 + t90 + t95 - t115 + t212;
t4 = pkin(5) * t71 - pkin(10) * t242 + t12;
t9 = -pkin(10) * t244 + t13;
t2 = qJD(6) * t15 + t193 * t4 + t197 * t9;
t20 = pkin(5) * t244 + t39;
t25 = t71 * Ifges(6,6) + (-t250 + t251) * t70;
t26 = t71 * Ifges(6,5) + (Ifges(6,1) * t192 - t252) * t70;
t3 = -qJD(6) * t16 - t193 * t9 + t197 * t4;
t47 = -Ifges(7,4) * t83 - Ifges(7,2) * t82 - Ifges(7,6) * t208;
t48 = -Ifges(7,1) * t83 - Ifges(7,4) * t82 - Ifges(7,5) * t208;
t57 = pkin(5) * t236 - t275;
t7 = Ifges(7,4) * t30 + Ifges(7,2) * t31 + Ifges(7,6) * t71;
t8 = Ifges(7,1) * t30 + Ifges(7,4) * t31 + Ifges(7,5) * t71;
t202 = mrSges(6,3) * t246 - t2 * t253 - t16 * t254 + (Ifges(6,1) * t191 + t251) * t242 / 0.2e1 - (Ifges(6,2) * t192 + t252) * t244 / 0.2e1 - t208 * t227 / 0.2e1 + (t15 * t150 - t160 * t3) * mrSges(7,3) + t191 * t26 / 0.2e1 + t192 * t25 / 0.2e1 - t151 * t47 / 0.2e1 - t207 * t7 / 0.2e1 + t160 * t8 / 0.2e1 - t150 * t48 / 0.2e1 + t20 * t130 + t31 * t131 / 0.2e1 + t30 * t132 / 0.2e1 + t57 * t119 - t82 * t120 / 0.2e1 - t83 * t121 / 0.2e1 - Ifges(5,6) * t71 + Ifges(5,5) * t70 - t38 * mrSges(5,2) + (t171 - mrSges(5,1)) * t39 + (Ifges(6,5) * t191 + Ifges(7,5) * t160 + Ifges(6,6) * t192 - Ifges(7,6) * t207) * t71 / 0.2e1;
t201 = mrSges(4,1) * t97 - t96 * mrSges(4,2) + Ifges(4,5) * t133 - Ifges(4,6) * t134 + t202;
t184 = -pkin(4) - t260;
t106 = -qJD(5) * t160 - qJD(6) * t136;
t89 = -qJD(6) * t117 - t160 * t177;
t86 = -mrSges(6,1) * t208 - mrSges(6,3) * t235;
t85 = mrSges(6,2) * t208 - mrSges(6,3) * t236;
t63 = -mrSges(7,1) * t208 + mrSges(7,3) * t83;
t62 = mrSges(7,2) * t208 - mrSges(7,3) * t82;
t55 = -qJD(6) * t92 - t118 * t160;
t46 = mrSges(6,1) * t71 - mrSges(6,3) * t242;
t45 = -mrSges(6,2) * t71 - mrSges(6,3) * t244;
t19 = -mrSges(7,2) * t71 + mrSges(7,3) * t31;
t18 = mrSges(7,1) * t71 - mrSges(7,3) * t30;
t1 = [t84 * t268 + (t15 * t3 + t16 * t2 + t20 * t57) * t269 + 0.2e1 * m(4) * (t139 * t97 + t140 * t96) + (mrSges(5,2) * t267 + mrSges(5,3) * t268 - t191 * t25 + t192 * t26 + (-(2 * Ifges(5,4)) + t211) * t71 + (Ifges(6,1) * t190 + (2 * Ifges(5,1)) + (t250 - 0.2e1 * t251) * t191) * t70) * t127 + ((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t200) * t220 + (-0.2e1 * pkin(1) * mrSges(3,1) + m(4) * pkin(2) * t266 + 0.2e1 * pkin(2) * (-mrSges(4,1) * t161 + mrSges(4,2) * t162) - 0.2e1 * Ifges(3,4) * t196 + (Ifges(3,1) - Ifges(3,2)) * t220) * t196) * qJD(2) - 0.2e1 * t161 * Ifges(4,2) * t134 + 0.2e1 * t133 * t162 * Ifges(4,1) + t71 * (-Ifges(7,5) * t83 - Ifges(7,6) * t82) + (-0.2e1 * (-Ifges(5,4) + t211) * t70 - mrSges(5,1) * t267 + t38 * t276 - ((2 * Ifges(5,2)) + (2 * Ifges(6,3)) + Ifges(7,3)) * t71 - t219) * t208 + 0.2e1 * t142 * (mrSges(5,1) * t71 + mrSges(5,2) * t70) - t82 * t7 - t83 * t8 + 0.2e1 * t13 * t85 + 0.2e1 * t12 * t86 + 0.2e1 * t2 * t62 + 0.2e1 * t3 * t63 + 0.2e1 * t49 * t46 + 0.2e1 * t50 * t45 + 0.2e1 * t20 * t52 + 0.2e1 * t57 * t10 + t31 * t47 + t30 * t48 + 0.2e1 * t15 * t18 + 0.2e1 * t16 * t19 - 0.2e1 * t275 * t42 + (-t275 * t70 - t71 * t79) * t276 + (mrSges(4,1) * t134 + mrSges(4,2) * t133) * t266 + 0.2e1 * (t133 * t161 - t134 * t162) * Ifges(4,4) + 0.2e1 * (-t133 * t139 - t134 * t140 + t161 * t96 - t162 * t97) * mrSges(4,3) + (t12 * t49 + t13 * t50 - t256) * t270 + (t114 * t142 + t38 * t79 - t256) * t271; t201 + (Ifges(3,5) * t200 - Ifges(3,6) * t196 + (-mrSges(3,1) * t200 + mrSges(3,2) * t196) * pkin(7)) * qJD(2) + (t118 * t85 + t143 * t45) * t192 + (-t118 * t86 - t143 * t46 - t258) * t191 + t255 * t124 + m(7) * (t124 * t57 + t141 * t20 + t15 * t55 + t16 * t54 + t2 * t92 + t3 * t91) + m(5) * (t123 * t79 - t148 * t39 + t149 * t38 - t247) + (t123 * t208 + t124 * t127 - t148 * t70 - t149 * t71) * mrSges(5,3) + t141 * t10 + t147 * t42 + t91 * t18 + t92 * t19 + t54 * t62 + t55 * t63 + m(6) * (t118 * t209 + t143 * t210 + t147 * t39 - t247) + (m(4) * (t195 * t96 + t199 * t97 + (-t139 * t195 + t140 * t199) * qJD(3)) + (-t199 * t133 - t195 * t134 + (t161 * t199 + t162 * t195) * qJD(3)) * mrSges(4,3)) * pkin(2); -0.2e1 * t248 + 0.2e1 * t104 + 0.2e1 * t110 + 0.2e1 * t111 - 0.2e1 * t115 - 0.2e1 * t51 - 0.2e1 * t81 + 0.2e1 * t90 + 0.2e1 * t95 + 0.2e1 * t206 + (t91 * t150 - t160 * t55) * t221 + (t124 * t141 + t54 * t92 + t55 * t91) * t269 + (t124 * t147 + t143 * t215) * t270 + (t123 * t149 - t124 * t148) * t271 + t212; t201 + t184 * t42 + t169 * t10 + t116 * t18 + t117 * t19 + t88 * t62 + t89 * t63 + (m(5) * (t194 * t38 - t198 * t39) + (-t194 * t71 - t198 * t70) * mrSges(5,3) + (t198 * t208 * mrSges(5,3) + (mrSges(5,3) * t127 + t255) * t194 + t57 * t261 - m(6) * t241 + m(5) * (t198 * t79 - t241)) * qJD(4)) * pkin(3) + m(7) * (t116 * t3 + t117 * t2 + t15 * t89 + t16 * t88 + t169 * t20) + m(6) * (t13 * t232 + t177 * t209 - t180 * t249 + t184 * t39) + (-t177 * t86 - t180 * t46 - t258) * t191 + (t177 * t85 + t180 * t45) * t192; t203 + ((-t55 - t89) * t160 - (-t116 - t91) * t150) * mrSges(7,3) + 0.2e1 * ((t123 * t194 - t124 * t198) * t265 + (t141 * t261 / 0.2e1 + t147 * t194 * t264 + (-t148 * t194 + t149 * t198) * t265) * qJD(4)) * pkin(3) + m(7) * (t116 * t55 + t117 * t54 + t124 * t169 + t88 * t92 + t89 * t91) + m(6) * (t124 * t184 + (t237 + t238) * t180 + (t233 + t234) * t143) + (-t123 - t217) * mrSges(5,2) + t206 - mrSges(5,1) * t218 + t272; 0.2e1 * t103 + 0.2e1 * t122 + 0.2e1 * t153 + 0.2e1 * t163 + 0.2e1 * t164 - 0.2e1 * t80 - 0.2e1 * t98 + 0.2e1 * t205 + (t116 * t150 - t160 * t89) * t221 + (t116 * t89 + t117 * t88 + t169 * t218) * t269 + (t180 * t214 + t184 * t218) * t270 + t212; t202 + (qJ(5) * t45 + qJD(5) * t85) * t192 + (-qJ(5) * t46 - qJD(5) * t86 - t258) * t191 + m(7) * (t105 * t16 + t106 * t15 + t135 * t3 + t136 * t2 + t181 * t20) + t181 * t10 + t135 * t18 + t136 * t19 + t105 * t62 + t106 * t63 + m(6) * (-pkin(4) * t39 + qJ(5) * t210 + qJD(5) * t209) - pkin(4) * t42; t203 + ((-t106 - t55) * t160 - (-t135 - t91) * t150) * mrSges(7,3) + m(7) * (t105 * t92 + t106 * t91 + t124 * t181 + t135 * t55 + t136 * t54) + m(6) * (-pkin(4) * t124 + qJ(5) * t215 + t143 * t213) - t248 + t273; t212 + m(6) * (-pkin(4) * t218 + qJ(5) * t214 + t180 * t213) + m(7) * (t105 * t117 + t106 * t116 + t135 * t89 + t136 * t88 + t181 * t218) + t205 + ((-t106 - t89) * t160 - (-t116 - t135) * t150) * mrSges(7,3) + t272 + t273; -0.2e1 * t107 + 0.2e1 * t108 + 0.2e1 * t182 + 0.2e1 * t183 - 0.2e1 * t93 + (-t106 * t160 + t135 * t150) * t221 + (t105 * t136 + t106 * t135) * t269 + qJ(5) * t213 * t270 + t212; m(6) * t39 + m(7) * t20 + t10 + t42; 0.2e1 * (m(7) / 0.2e1 + t264) * t124 + t119; (m(6) + m(7)) * t218 + t119; t119; 0; mrSges(7,1) * t3 - mrSges(7,2) * t2 + t219; mrSges(7,1) * t55 - mrSges(7,2) * t54 + t227; mrSges(7,1) * t89 - mrSges(7,2) * t88 + t227; mrSges(7,1) * t106 - mrSges(7,2) * t105 + t227; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
