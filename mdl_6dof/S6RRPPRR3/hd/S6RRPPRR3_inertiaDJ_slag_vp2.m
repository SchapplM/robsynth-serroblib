% Calculate time derivative of joint inertia matrix for
% S6RRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 09:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:54:37
% EndTime: 2019-03-09 08:54:48
% DurationCPUTime: 4.60s
% Computational Cost: add. (10213->515), mult. (27848->772), div. (0->0), fcn. (28573->12), ass. (0->226)
t181 = sin(pkin(6));
t283 = 0.2e1 * t181;
t179 = sin(pkin(12));
t182 = cos(pkin(12));
t186 = sin(qJ(5));
t189 = cos(qJ(5));
t153 = t179 * t186 - t189 * t182;
t148 = t153 * qJD(5);
t154 = t179 * t189 + t182 * t186;
t188 = cos(qJ(6));
t185 = sin(qJ(6));
t220 = qJD(6) * t185;
t193 = t148 * t188 + t154 * t220;
t282 = t185 / 0.2e1;
t261 = t188 / 0.2e1;
t219 = qJD(6) * t188;
t194 = t148 * t185 - t154 * t219;
t180 = sin(pkin(11));
t183 = cos(pkin(11));
t187 = sin(qJ(2));
t190 = cos(qJ(2));
t139 = (t180 * t190 + t183 * t187) * t181;
t184 = cos(pkin(6));
t121 = t139 * t182 + t179 * t184;
t227 = t183 * t190;
t229 = t181 * t187;
t138 = t180 * t229 - t181 * t227;
t260 = pkin(1) * t184;
t170 = t190 * t260;
t254 = -pkin(8) - qJ(3);
t210 = t254 * t187;
t122 = pkin(2) * t184 + t181 * t210 + t170;
t169 = t187 * t260;
t228 = t181 * t190;
t152 = pkin(8) * t228 + t169;
t132 = qJ(3) * t228 + t152;
t90 = t180 * t122 + t183 * t132;
t83 = qJ(4) * t184 + t90;
t155 = (-pkin(2) * t190 - pkin(1)) * t181;
t93 = t138 * pkin(3) - t139 * qJ(4) + t155;
t57 = -t179 * t83 + t182 * t93;
t40 = pkin(4) * t138 - pkin(9) * t121 + t57;
t120 = -t139 * t179 + t182 * t184;
t58 = t179 * t93 + t182 * t83;
t48 = pkin(9) * t120 + t58;
t250 = t186 * t40 + t189 * t48;
t133 = qJD(2) * t139;
t224 = qJD(2) * t181;
t134 = (-t180 * t187 + t227) * t224;
t233 = t134 * t182;
t167 = qJD(2) * t170;
t113 = t167 + (qJD(2) * t210 + qJD(3) * t190) * t181;
t191 = -qJD(3) * t229 + (t228 * t254 - t169) * qJD(2);
t79 = t183 * t113 + t180 * t191;
t71 = qJD(4) * t184 + t79;
t214 = t187 * t224;
t207 = pkin(2) * t214;
t75 = pkin(3) * t133 - qJ(4) * t134 - qJD(4) * t139 + t207;
t41 = -t179 * t71 + t182 * t75;
t32 = pkin(4) * t133 - pkin(9) * t233 + t41;
t234 = t134 * t179;
t42 = t179 * t75 + t182 * t71;
t34 = -pkin(9) * t234 + t42;
t6 = -qJD(5) * t250 - t186 * t34 + t189 * t32;
t281 = 2 * m(5);
t280 = 2 * m(6);
t279 = 2 * m(7);
t278 = -2 * mrSges(3,3);
t277 = -2 * mrSges(4,3);
t276 = -2 * mrSges(6,3);
t259 = pkin(2) * t180;
t171 = qJ(4) + t259;
t253 = pkin(9) + t171;
t150 = t253 * t182;
t209 = qJD(5) * t253;
t217 = t189 * qJD(4);
t218 = t186 * qJD(4);
t222 = qJD(5) * t189;
t92 = t182 * t218 + t150 * t222 + (-t186 * t209 + t217) * t179;
t275 = 0.2e1 * t92;
t274 = 0.2e1 * t155;
t273 = m(4) * pkin(2);
t196 = t189 * t120 - t121 * t186;
t55 = qJD(5) * t196 - t134 * t153;
t81 = t120 * t186 + t121 * t189;
t64 = t138 * t185 + t188 * t81;
t29 = -qJD(6) * t64 + t133 * t188 - t185 * t55;
t272 = t29 / 0.2e1;
t149 = t154 * qJD(5);
t62 = -Ifges(7,1) * t193 + Ifges(7,4) * t194 + Ifges(7,5) * t149;
t271 = t62 / 0.2e1;
t63 = t138 * t188 - t185 * t81;
t270 = t63 / 0.2e1;
t246 = Ifges(7,4) * t185;
t204 = Ifges(7,1) * t188 - t246;
t101 = Ifges(7,5) * t153 + t154 * t204;
t269 = t101 / 0.2e1;
t268 = -t154 / 0.2e1;
t174 = Ifges(7,5) * t219;
t267 = -Ifges(7,6) * t220 / 0.2e1 + t174 / 0.2e1;
t159 = t204 * qJD(6);
t266 = t159 / 0.2e1;
t265 = Ifges(7,5) * t282 + Ifges(7,6) * t261;
t163 = Ifges(7,2) * t188 + t246;
t264 = -t163 / 0.2e1;
t245 = Ifges(7,4) * t188;
t164 = Ifges(7,1) * t185 + t245;
t263 = t164 / 0.2e1;
t262 = -t185 / 0.2e1;
t258 = pkin(2) * t183;
t257 = pkin(5) * t149;
t78 = t113 * t180 - t183 * t191;
t256 = t78 * mrSges(4,1);
t255 = t79 * mrSges(4,2);
t28 = qJD(6) * t63 + t133 * t185 + t188 * t55;
t12 = -mrSges(7,1) * t29 + mrSges(7,2) * t28;
t49 = mrSges(6,1) * t133 - mrSges(6,3) * t55;
t252 = t12 - t49;
t35 = -mrSges(7,1) * t63 + mrSges(7,2) * t64;
t67 = mrSges(6,1) * t138 - mrSges(6,3) * t81;
t251 = t35 - t67;
t249 = mrSges(7,3) * t154;
t248 = Ifges(5,4) * t179;
t247 = Ifges(5,4) * t182;
t244 = Ifges(7,6) * t185;
t211 = t179 * t253;
t103 = t186 * t150 + t189 * t211;
t243 = t103 * t92;
t242 = t133 * Ifges(5,5);
t241 = t133 * Ifges(6,5);
t240 = t133 * Ifges(5,6);
t239 = t133 * Ifges(6,6);
t142 = -pkin(8) * t214 + t167;
t238 = t142 * mrSges(3,2);
t143 = t152 * qJD(2);
t237 = t143 * mrSges(3,1);
t173 = -pkin(3) - t258;
t160 = -pkin(4) * t182 + t173;
t102 = pkin(5) * t153 - pkin(10) * t154 + t160;
t104 = t189 * t150 - t186 * t211;
t68 = t102 * t188 - t104 * t185;
t236 = qJD(6) * t68;
t69 = t102 * t185 + t104 * t188;
t235 = qJD(6) * t69;
t230 = t149 * t153;
t94 = mrSges(5,1) * t234 + mrSges(5,2) * t233;
t226 = -Ifges(6,5) * t148 - Ifges(6,6) * t149;
t223 = qJD(5) * t186;
t221 = qJD(6) * t154;
t56 = qJD(5) * t81 + t134 * t154;
t7 = Ifges(7,5) * t28 + Ifges(7,6) * t29 + Ifges(7,3) * t56;
t216 = Ifges(6,5) * t55 - Ifges(6,6) * t56 + Ifges(6,3) * t133;
t215 = Ifges(3,5) * t190 * t224 + Ifges(4,5) * t134 - Ifges(4,6) * t133;
t25 = t56 * mrSges(6,1) + t55 * mrSges(6,2);
t108 = t149 * mrSges(6,1) - t148 * mrSges(6,2);
t208 = t148 * (t185 ^ 2 + t188 ^ 2);
t89 = t122 * t183 - t180 * t132;
t16 = pkin(10) * t138 + t250;
t84 = -pkin(3) * t184 - t89;
t65 = -pkin(4) * t120 + t84;
t30 = -pkin(5) * t196 - pkin(10) * t81 + t65;
t10 = -t16 * t185 + t188 * t30;
t59 = pkin(4) * t234 + t78;
t17 = pkin(5) * t56 - pkin(10) * t55 + t59;
t5 = t186 * t32 + t189 * t34 + t40 * t222 - t223 * t48;
t3 = pkin(10) * t133 + t5;
t1 = qJD(6) * t10 + t17 * t185 + t188 * t3;
t11 = t16 * t188 + t185 * t30;
t2 = -qJD(6) * t11 + t17 * t188 - t185 * t3;
t206 = t1 * t188 - t185 * t2;
t161 = -mrSges(7,1) * t188 + mrSges(7,2) * t185;
t205 = mrSges(7,1) * t185 + mrSges(7,2) * t188;
t203 = -Ifges(7,2) * t185 + t245;
t202 = -t10 * t185 + t11 * t188;
t13 = mrSges(7,1) * t56 - mrSges(7,3) * t28;
t14 = -mrSges(7,2) * t56 + mrSges(7,3) * t29;
t201 = -t185 * t13 + t188 * t14;
t43 = mrSges(7,2) * t196 + mrSges(7,3) * t63;
t44 = -mrSges(7,1) * t196 - mrSges(7,3) * t64;
t200 = -t185 * t44 + t188 * t43;
t18 = -t186 * t48 + t189 * t40;
t197 = t103 * t149 + t153 * t92;
t23 = Ifges(7,4) * t64 + Ifges(7,2) * t63 - Ifges(7,6) * t196;
t24 = Ifges(7,1) * t64 + Ifges(7,4) * t63 - Ifges(7,5) * t196;
t195 = t23 * t262 + t24 * t261;
t192 = (-t10 * t188 - t11 * t185) * qJD(6) + t206;
t60 = -Ifges(7,5) * t193 + Ifges(7,6) * t194 + Ifges(7,3) * t149;
t158 = t203 * qJD(6);
t156 = t205 * qJD(6);
t151 = -pkin(8) * t229 + t170;
t128 = t134 * mrSges(4,2);
t118 = Ifges(6,1) * t154 - Ifges(6,4) * t153;
t117 = Ifges(6,4) * t154 - Ifges(6,2) * t153;
t116 = mrSges(7,1) * t153 - t188 * t249;
t115 = -mrSges(7,2) * t153 - t185 * t249;
t112 = pkin(10) * t148 + t257;
t111 = t205 * t154;
t110 = -Ifges(6,1) * t148 - Ifges(6,4) * t149;
t109 = -Ifges(6,4) * t148 - Ifges(6,2) * t149;
t100 = Ifges(7,6) * t153 + t154 * t203;
t99 = Ifges(7,3) * t153 + (Ifges(7,5) * t188 - t244) * t154;
t98 = mrSges(5,1) * t133 - mrSges(5,3) * t233;
t97 = -mrSges(5,2) * t133 - mrSges(5,3) * t234;
t96 = mrSges(5,1) * t138 - mrSges(5,3) * t121;
t95 = -mrSges(5,2) * t138 + mrSges(5,3) * t120;
t91 = t182 * t217 - t150 * t223 + (-t189 * t209 - t218) * t179;
t86 = -mrSges(7,2) * t149 + mrSges(7,3) * t194;
t85 = mrSges(7,1) * t149 + mrSges(7,3) * t193;
t77 = t242 + (t182 * Ifges(5,1) - t248) * t134;
t76 = t240 + (-t179 * Ifges(5,2) + t247) * t134;
t74 = t194 * mrSges(7,1) + t193 * mrSges(7,2);
t66 = -mrSges(6,2) * t138 + mrSges(6,3) * t196;
t61 = -Ifges(7,4) * t193 + Ifges(7,2) * t194 + Ifges(7,6) * t149;
t50 = -mrSges(6,2) * t133 - mrSges(6,3) * t56;
t46 = Ifges(6,1) * t81 + Ifges(6,4) * t196 + Ifges(6,5) * t138;
t45 = Ifges(6,4) * t81 + Ifges(6,2) * t196 + Ifges(6,6) * t138;
t37 = t112 * t188 - t185 * t91 - t235;
t36 = t112 * t185 + t188 * t91 + t236;
t22 = Ifges(7,5) * t64 + Ifges(7,6) * t63 - Ifges(7,3) * t196;
t21 = Ifges(6,1) * t55 - Ifges(6,4) * t56 + t241;
t20 = Ifges(6,4) * t55 - Ifges(6,2) * t56 + t239;
t15 = -pkin(5) * t138 - t18;
t9 = Ifges(7,1) * t28 + Ifges(7,4) * t29 + Ifges(7,5) * t56;
t8 = Ifges(7,4) * t28 + Ifges(7,2) * t29 + Ifges(7,6) * t56;
t4 = -pkin(5) * t133 - t6;
t19 = [-t196 * t7 + t196 * t20 + 0.2e1 * t59 * (-mrSges(6,1) * t196 + mrSges(6,2) * t81) + (mrSges(4,1) * t274 + t90 * t277 - 0.2e1 * t139 * Ifges(4,4) + Ifges(5,5) * t121 + Ifges(6,5) * t81 - Ifges(4,6) * t184 + Ifges(5,6) * t120 + Ifges(6,6) * t196 + ((2 * Ifges(4,2)) + (2 * Ifges(5,3)) + Ifges(6,3)) * t138) * t133 + (t41 * t57 + t42 * t58 + t78 * t84) * t281 + (t1 * t11 + t10 * t2 + t15 * t4) * t279 + t138 * t216 + t120 * t76 + 0.2e1 * t78 * (-mrSges(5,1) * t120 + mrSges(5,2) * t121) + t121 * t77 + 0.2e1 * t58 * t97 + 0.2e1 * t57 * t98 + 0.2e1 * t84 * t94 + 0.2e1 * t42 * t95 + 0.2e1 * t41 * t96 + t81 * t21 + t64 * t9 + 0.2e1 * t65 * t25 + 0.2e1 * t5 * t66 + 0.2e1 * t6 * t67 + t63 * t8 + t55 * t46 + 0.2e1 * t18 * t49 + 0.2e1 * t1 * t43 + 0.2e1 * t2 * t44 + 0.2e1 * t4 * t35 + t28 * t24 + t29 * t23 + 0.2e1 * t10 * t13 + 0.2e1 * t11 * t14 + 0.2e1 * t15 * t12 + t128 * t274 + (0.2e1 * (t142 * t190 + t143 * t187) * mrSges(3,3) + ((t151 * t278 + Ifges(3,5) * t184 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t190) * t283) * t190 + (0.2e1 * pkin(2) * (mrSges(4,1) * t138 + mrSges(4,2) * t139) + t152 * t278 + t273 * t274 - 0.2e1 * Ifges(3,6) * t184 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t187 + (Ifges(3,1) - Ifges(3,2)) * t190) * t283) * t187) * qJD(2)) * t181 + (t182 * (Ifges(5,1) * t121 + Ifges(5,4) * t120) - t179 * (Ifges(5,4) * t121 + Ifges(5,2) * t120) + t89 * t277 + 0.2e1 * Ifges(4,1) * t139 + Ifges(4,5) * t184 + 0.2e1 * (Ifges(5,5) * t182 - Ifges(5,6) * t179 - Ifges(4,4)) * t138) * t134 + (t215 - 0.2e1 * t237 - 0.2e1 * t238 - 0.2e1 * t255 - 0.2e1 * t256) * t184 + (t18 * t6 + t250 * t5 + t59 * t65) * t280 + 0.2e1 * t250 * t50 + 0.2e1 * m(3) * (t142 * t152 - t143 * t151) + 0.2e1 * m(4) * (-t78 * t89 + t79 * t90) + 0.2e1 * (-t138 * t79 + t139 * t78) * mrSges(4,3) + (t22 - t45) * t56; -(-t109 / 0.2e1 + t60 / 0.2e1) * t196 + t251 * t92 + t252 * t103 + (-t5 * mrSges(6,3) + t7 / 0.2e1 - t20 / 0.2e1 - t239 / 0.2e1 + t59 * mrSges(6,1)) * t153 + (t171 * t97 + qJD(4) * t95 + t42 * mrSges(5,3) + t240 / 0.2e1 - t78 * mrSges(5,1) + t76 / 0.2e1) * t182 + (-t171 * t98 - qJD(4) * t96 - t41 * mrSges(5,3) + t242 / 0.2e1 + t78 * mrSges(5,2) + t77 / 0.2e1) * t179 + m(5) * (t173 * t78 + (-t41 * t179 + t42 * t182) * t171 + (-t179 * t57 + t182 * t58) * qJD(4)) - Ifges(3,6) * t214 + (-t117 / 0.2e1 + t99 / 0.2e1) * t56 - (-t18 * mrSges(6,3) + t46 / 0.2e1 + t195) * t148 + t173 * t94 + t160 * t25 + t55 * t118 / 0.2e1 + t65 * t108 + t81 * t110 / 0.2e1 + t4 * t111 + t1 * t115 + t2 * t116 + t104 * t50 + t11 * t86 + t91 * t66 + t10 * t85 - t15 * t74 + t68 * t13 + t69 * t14 + t36 * t43 + t37 * t44 + t28 * t269 + t61 * t270 + t64 * t271 + t100 * t272 + (t180 * t79 - t183 * t78) * t273 - t133 * mrSges(4,3) * t259 + t215 + t138 * t226 / 0.2e1 + m(7) * (t1 * t69 + t10 * t37 + t103 * t4 + t11 * t36 + t15 * t92 + t2 * t68) + (t182 * (Ifges(5,1) * t179 + t247) / 0.2e1 - t179 * (Ifges(5,2) * t182 + t248) / 0.2e1 - mrSges(4,3) * t258) * t134 + (t9 * t261 + t8 * t262 - t6 * mrSges(6,3) + t241 / 0.2e1 + t59 * mrSges(6,2) + t21 / 0.2e1 + (-t188 * t23 / 0.2e1 + t24 * t262) * qJD(6)) * t154 + (-t250 * mrSges(6,3) + t22 / 0.2e1 - t45 / 0.2e1) * t149 + m(6) * (-t103 * t6 + t104 * t5 + t160 * t59 - t18 * t92 + t250 * t91) - t237 - t238 - t255 - t256; -0.2e1 * t103 * t74 + 0.2e1 * t160 * t108 + t111 * t275 + 0.2e1 * t36 * t115 + 0.2e1 * t37 * t116 + 0.2e1 * t68 * t85 + 0.2e1 * t69 * t86 + (t36 * t69 + t37 * t68 + t243) * t279 + (t104 * t91 + t243) * t280 + (t276 * t91 - t109 + t60) * t153 + (t104 * t276 - t117 + t99) * t149 - (0.2e1 * mrSges(6,3) * t103 - t185 * t100 + t188 * t101 + t118) * t148 + (mrSges(6,3) * t275 - t185 * t61 + t188 * t62 + t110 + (-t100 * t188 - t101 * t185) * qJD(6)) * t154 + (t171 * t281 + 0.2e1 * mrSges(5,3)) * qJD(4) * (t179 ^ 2 + t182 ^ 2); m(4) * t207 + t133 * mrSges(4,1) + t179 * t97 + t182 * t98 + t128 + t252 * t153 + t251 * t149 - (t200 + t66) * t148 + (t50 + (-t185 * t43 - t188 * t44) * qJD(6) + t201) * t154 + m(7) * (-t148 * t202 + t149 * t15 + t153 * t4 + t154 * t192) + m(6) * (-t148 * t250 - t149 * t18 - t153 * t6 + t154 * t5) + m(5) * (t179 * t42 + t182 * t41); t149 * t111 - t153 * t74 + m(7) * t197 + m(6) * (-t104 * t148 + t154 * t91 + t197) + (m(7) * (-t148 * t69 + t154 * t36 - t221 * t68) - t148 * t115 + t154 * t86 - t116 * t221) * t188 + (m(7) * (t148 * t68 - t154 * t37 - t221 * t69) - t115 * t221 + t148 * t116 - t154 * t85) * t185; 0.2e1 * m(6) * (-t148 * t154 + t230) + 0.2e1 * m(7) * (-t154 * t208 + t230); t188 * t13 + t185 * t14 + t200 * qJD(6) + m(7) * (qJD(6) * t202 + t1 * t185 + t188 * t2) + m(6) * t59 + m(5) * t78 + t25 + t94; m(7) * (t185 * t36 + t188 * t37 + (-t185 * t68 + t188 * t69) * qJD(6)) + t115 * t219 + t185 * t86 - t116 * t220 + t188 * t85 + t108; 0; 0; t8 * t261 + t9 * t282 + t15 * t156 - t196 * t267 + t158 * t270 + t64 * t266 + t4 * t161 + t56 * t265 + t163 * t272 + t28 * t263 - t5 * mrSges(6,2) + t6 * mrSges(6,1) + t195 * qJD(6) + (-m(7) * t4 - t12) * pkin(5) + t192 * mrSges(7,3) + (m(7) * (-t10 * t219 - t11 * t220 + t206) - t44 * t219 - t43 * t220 + t201) * pkin(10) + t216; t103 * t156 + t153 * t267 + t149 * t265 - t91 * mrSges(6,2) + pkin(5) * t74 + (-m(7) * pkin(5) - mrSges(6,1) + t161) * t92 + (t154 * t266 - t148 * t263 + t36 * mrSges(7,3) + t61 / 0.2e1 + (-t68 * mrSges(7,3) + t154 * t264 + t269) * qJD(6) + (m(7) * (t36 - t236) + t86 - qJD(6) * t116) * pkin(10)) * t188 + (t158 * t268 - t148 * t264 - t37 * mrSges(7,3) + t271 + (-t100 / 0.2e1 + t164 * t268 - t69 * mrSges(7,3)) * qJD(6) + (m(7) * (-t37 - t235) - qJD(6) * t115 - t85) * pkin(10)) * t185 + t226; t149 * t161 + t153 * t156 + m(7) * (-pkin(10) * t208 - t257) - mrSges(7,3) * t208 - t108; 0; -0.2e1 * pkin(5) * t156 + t158 * t188 + t159 * t185 + (-t163 * t185 + t164 * t188) * qJD(6); mrSges(7,1) * t2 - mrSges(7,2) * t1 + t7; mrSges(7,1) * t37 - mrSges(7,2) * t36 + t60; t74; -t156; t174 + (pkin(10) * t161 - t244) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t19(1) t19(2) t19(4) t19(7) t19(11) t19(16); t19(2) t19(3) t19(5) t19(8) t19(12) t19(17); t19(4) t19(5) t19(6) t19(9) t19(13) t19(18); t19(7) t19(8) t19(9) t19(10) t19(14) t19(19); t19(11) t19(12) t19(13) t19(14) t19(15) t19(20); t19(16) t19(17) t19(18) t19(19) t19(20) t19(21);];
Mq  = res;
