% Calculate time derivative of joint inertia matrix for
% S6RPRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR12_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR12_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR12_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR12_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR12_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:48:13
% EndTime: 2019-03-09 05:48:27
% DurationCPUTime: 5.98s
% Computational Cost: add. (10459->575), mult. (30459->823), div. (0->0), fcn. (31891->12), ass. (0->228)
t184 = sin(qJ(6));
t187 = cos(qJ(6));
t188 = cos(qJ(4));
t232 = qJD(6) * t188;
t185 = sin(qJ(4));
t234 = qJD(4) * t185;
t197 = t184 * t232 + t187 * t234;
t196 = t184 * t234 - t187 * t232;
t178 = sin(pkin(12));
t180 = sin(pkin(6));
t181 = cos(pkin(12));
t183 = cos(pkin(6));
t186 = sin(qJ(3));
t182 = cos(pkin(7));
t189 = cos(qJ(3));
t246 = t182 * t189;
t179 = sin(pkin(7));
t249 = t179 * t189;
t291 = t180 * (-t178 * t186 + t181 * t246) + t183 * t249;
t290 = -2 * Ifges(4,4);
t273 = -t184 / 0.2e1;
t154 = mrSges(7,1) * t184 + mrSges(7,2) * t187;
t289 = -mrSges(6,3) - t154;
t248 = t180 * t181;
t271 = pkin(1) * t183;
t238 = qJ(2) * t248 + t178 * t271;
t109 = (t179 * t183 + t182 * t248) * pkin(9) + t238;
t171 = t181 * t271;
t252 = t178 * t180;
t112 = pkin(2) * t183 + t171 + (-pkin(9) * t182 - qJ(2)) * t252;
t123 = (-pkin(9) * t178 * t179 - pkin(2) * t181 - pkin(1)) * t180;
t206 = t112 * t182 + t123 * t179;
t66 = -t186 * t109 + t206 * t189;
t250 = t179 * t186;
t227 = t185 * t250;
t131 = -t188 * t182 + t227;
t118 = t131 * t187 + t184 * t249;
t200 = -t131 * t184 + t187 * t249;
t132 = t182 * t185 + t188 * t250;
t235 = qJD(3) * t189;
t222 = t179 * t235;
t117 = qJD(4) * t132 + t185 * t222;
t236 = qJD(3) * t186;
t223 = t179 * t236;
t86 = qJD(6) * t118 + t117 * t184 + t187 * t223;
t87 = qJD(6) * t200 + t117 * t187 - t184 * t223;
t286 = qJD(6) * (t118 * t184 + t187 * t200) - t86 * t184 - t87 * t187;
t253 = qJ(5) * t185;
t278 = pkin(4) + pkin(11);
t134 = -t188 * t278 - pkin(3) - t253;
t277 = pkin(5) + pkin(10);
t163 = t277 * t185;
t113 = -t134 * t184 + t163 * t187;
t114 = t134 * t187 + t163 * t184;
t265 = mrSges(7,3) * t188;
t146 = mrSges(7,1) * t185 + t184 * t265;
t147 = -mrSges(7,2) * t185 - t187 * t265;
t285 = m(7) * (-t113 * t184 + t114 * t187) + t187 * t147 - t184 * t146;
t284 = 0.2e1 * m(7);
t283 = m(6) / 0.2e1;
t282 = m(5) * pkin(3);
t247 = t182 * t186;
t111 = t183 * t250 + (t178 * t189 + t181 * t247) * t180;
t107 = t111 * qJD(3);
t106 = t291 * qJD(3);
t130 = -t179 * t248 + t182 * t183;
t90 = t111 * t188 + t130 * t185;
t73 = qJD(4) * t90 + t106 * t185;
t89 = t111 * t185 - t130 * t188;
t75 = t184 * t291 + t187 * t89;
t33 = qJD(6) * t75 + t107 * t187 + t184 * t73;
t281 = t33 / 0.2e1;
t280 = t75 / 0.2e1;
t262 = Ifges(7,4) * t184;
t160 = Ifges(7,1) * t187 - t262;
t261 = Ifges(7,4) * t187;
t213 = Ifges(7,1) * t184 + t261;
t96 = -t160 * t232 + (Ifges(7,5) * t188 + t185 * t213) * qJD(4);
t279 = t96 / 0.2e1;
t211 = -Ifges(7,5) * t184 - Ifges(7,6) * t187;
t276 = t211 * qJD(6) / 0.2e1;
t275 = Ifges(7,5) * t187 / 0.2e1 + Ifges(7,6) * t273;
t158 = -Ifges(7,2) * t184 + t261;
t274 = t158 / 0.2e1;
t272 = -t187 / 0.2e1;
t47 = mrSges(6,1) * t73 - mrSges(6,3) * t107;
t49 = -mrSges(5,2) * t107 - mrSges(5,3) * t73;
t270 = -t47 + t49;
t74 = -t111 * t234 + (qJD(4) * t130 + t106) * t188;
t48 = t74 * mrSges(6,1) + t107 * mrSges(6,2);
t50 = mrSges(5,1) * t107 - mrSges(5,3) * t74;
t269 = t48 - t50;
t88 = -t112 * t179 + t182 * t123;
t58 = -pkin(3) * t291 - pkin(10) * t111 + t88;
t99 = t189 * t109;
t67 = t112 * t247 + t123 * t250 + t99;
t63 = pkin(10) * t130 + t67;
t25 = t185 * t58 + t188 * t63;
t77 = mrSges(6,1) * t89 + mrSges(6,3) * t291;
t79 = mrSges(5,2) * t291 - mrSges(5,3) * t89;
t268 = t77 - t79;
t78 = mrSges(6,1) * t90 - mrSges(6,2) * t291;
t80 = -mrSges(5,1) * t291 - mrSges(5,3) * t90;
t267 = t78 - t80;
t266 = mrSges(4,3) * t107;
t264 = Ifges(5,4) * t185;
t263 = Ifges(5,4) * t188;
t260 = Ifges(6,6) * t185;
t259 = Ifges(6,6) * t188;
t258 = t106 * mrSges(4,3);
t257 = t291 * Ifges(6,4);
t256 = t291 * Ifges(5,6);
t237 = qJD(2) * t180;
t54 = (t178 * t246 + t181 * t186) * t237 + (t186 * t206 + t99) * qJD(3);
t255 = t189 * t54;
t153 = -t188 * mrSges(5,1) + t185 * mrSges(5,2);
t254 = t153 - mrSges(4,1);
t233 = qJD(4) * t188;
t116 = qJD(4) * t227 - t182 * t233 - t188 * t222;
t93 = t132 * t116;
t244 = t184 * t160;
t243 = t184 * t278;
t241 = t187 * t158;
t240 = t187 * t278;
t239 = Ifges(4,5) * t106 - Ifges(4,6) * t107;
t76 = t184 * t89 - t187 * t291;
t32 = -qJD(6) * t76 - t107 * t184 + t187 * t73;
t7 = Ifges(7,5) * t33 + Ifges(7,6) * t32 + Ifges(7,3) * t74;
t230 = Ifges(6,1) * t107 - Ifges(6,4) * t74 + Ifges(6,5) * t73;
t229 = Ifges(5,5) * t74 - Ifges(5,6) * t73 + Ifges(5,3) * t107;
t228 = m(6) * pkin(10) + mrSges(6,1);
t164 = t277 * t188;
t224 = t178 * t237;
t24 = -t185 * t63 + t188 * t58;
t217 = pkin(4) * t234 - qJD(5) * t185;
t216 = t179 * t224;
t22 = qJ(5) * t291 - t25;
t17 = pkin(5) * t90 + t278 * t291 - t24;
t62 = -pkin(3) * t130 - t66;
t193 = -qJ(5) * t90 + t62;
t21 = t278 * t89 + t193;
t3 = t17 * t187 - t184 * t21;
t4 = t17 * t184 + t187 * t21;
t215 = t184 * t3 - t187 * t4;
t214 = mrSges(7,1) * t187 - mrSges(7,2) * t184;
t152 = t188 * mrSges(6,2) - t185 * mrSges(6,3);
t212 = Ifges(7,2) * t187 + t262;
t210 = -pkin(4) * t188 - t253;
t41 = -mrSges(7,2) * t90 + mrSges(7,3) * t75;
t42 = mrSges(7,1) * t90 - mrSges(7,3) * t76;
t209 = -t184 * t42 + t187 * t41;
t124 = (pkin(11) * t185 - qJ(5) * t188) * qJD(4) + t217;
t149 = qJD(4) * t164;
t84 = qJD(6) * t113 + t124 * t187 + t149 * t184;
t85 = -qJD(6) * t114 - t124 * t184 + t149 * t187;
t208 = -t184 * t84 - t187 * t85;
t53 = (-t178 * t247 + t181 * t189) * t237 + t66 * qJD(3);
t82 = pkin(3) * t107 - pkin(10) * t106 + t216;
t14 = -t185 * t53 + t188 * t82 - t63 * t233 - t58 * t234;
t204 = -t116 * t188 + t117 * t185;
t202 = -qJ(5) * t116 + qJD(5) * t132;
t27 = Ifges(7,4) * t76 + Ifges(7,2) * t75 + Ifges(7,6) * t90;
t28 = Ifges(7,1) * t76 + Ifges(7,4) * t75 + Ifges(7,5) * t90;
t201 = t27 * t272 + t273 * t28;
t13 = t185 * t82 + t188 * t53 + t58 * t233 - t234 * t63;
t94 = t196 * Ifges(7,5) + Ifges(7,6) * t197 + Ifges(7,3) * t233;
t10 = -qJ(5) * t107 + qJD(5) * t291 - t13;
t192 = m(7) * t286;
t191 = -qJ(5) * t74 - qJD(5) * t90 + t54;
t175 = Ifges(5,5) * t233;
t174 = Ifges(6,5) * t234;
t161 = Ifges(5,1) * t185 + t263;
t159 = Ifges(5,2) * t188 + t264;
t156 = -Ifges(6,2) * t185 - t259;
t155 = -Ifges(6,3) * t188 - t260;
t150 = -pkin(3) + t210;
t148 = t277 * t234;
t144 = (Ifges(5,1) * t188 - t264) * qJD(4);
t143 = t213 * qJD(6);
t142 = (-Ifges(5,2) * t185 + t263) * qJD(4);
t141 = t212 * qJD(6);
t139 = (-Ifges(6,2) * t188 + t260) * qJD(4);
t138 = (Ifges(6,3) * t185 - t259) * qJD(4);
t137 = (mrSges(5,1) * t185 + mrSges(5,2) * t188) * qJD(4);
t136 = (-mrSges(6,2) * t185 - mrSges(6,3) * t188) * qJD(4);
t135 = t214 * qJD(6);
t133 = t214 * t188;
t129 = -qJ(5) * t233 + t217;
t128 = Ifges(7,5) * t185 - t188 * t213;
t127 = Ifges(7,6) * t185 - t188 * t212;
t126 = Ifges(7,3) * t185 + t188 * t211;
t122 = mrSges(7,1) * t233 - mrSges(7,3) * t196;
t121 = -mrSges(7,2) * t233 + mrSges(7,3) * t197;
t108 = -mrSges(7,1) * t197 + mrSges(7,2) * t196;
t95 = -t158 * t232 + (Ifges(7,6) * t188 + t185 * t212) * qJD(4);
t92 = mrSges(4,1) * t130 - mrSges(4,3) * t111;
t91 = -mrSges(4,2) * t130 + mrSges(4,3) * t291;
t83 = mrSges(4,1) * t107 + mrSges(4,2) * t106;
t65 = -mrSges(6,2) * t89 - mrSges(6,3) * t90;
t64 = mrSges(5,1) * t89 + mrSges(5,2) * t90;
t46 = Ifges(5,1) * t90 - Ifges(5,4) * t89 - Ifges(5,5) * t291;
t45 = Ifges(5,4) * t90 - Ifges(5,2) * t89 - t256;
t44 = -Ifges(6,2) * t90 + Ifges(6,6) * t89 - t257;
t43 = -Ifges(6,5) * t291 - Ifges(6,6) * t90 + Ifges(6,3) * t89;
t40 = -mrSges(7,1) * t75 + mrSges(7,2) * t76;
t39 = -mrSges(6,2) * t73 - mrSges(6,3) * t74;
t38 = mrSges(5,1) * t73 + mrSges(5,2) * t74;
t37 = Ifges(5,1) * t74 - Ifges(5,4) * t73 + t107 * Ifges(5,5);
t36 = Ifges(5,4) * t74 - Ifges(5,2) * t73 + t107 * Ifges(5,6);
t35 = t107 * Ifges(6,4) - Ifges(6,2) * t74 + Ifges(6,6) * t73;
t34 = t107 * Ifges(6,5) - Ifges(6,6) * t74 + Ifges(6,3) * t73;
t31 = pkin(4) * t89 + t193;
t26 = Ifges(7,5) * t76 + Ifges(7,6) * t75 + Ifges(7,3) * t90;
t23 = pkin(4) * t291 - t24;
t20 = mrSges(7,1) * t74 - mrSges(7,3) * t33;
t19 = -mrSges(7,2) * t74 + mrSges(7,3) * t32;
t18 = -pkin(5) * t89 - t22;
t16 = pkin(4) * t73 + t191;
t15 = -mrSges(7,1) * t32 + mrSges(7,2) * t33;
t12 = t278 * t73 + t191;
t11 = -pkin(4) * t107 - t14;
t9 = Ifges(7,1) * t33 + Ifges(7,4) * t32 + Ifges(7,5) * t74;
t8 = Ifges(7,4) * t33 + Ifges(7,2) * t32 + Ifges(7,6) * t74;
t6 = -pkin(5) * t73 - t10;
t5 = pkin(5) * t74 - t107 * t278 - t14;
t2 = -qJD(6) * t4 - t12 * t184 + t187 * t5;
t1 = qJD(6) * t3 + t12 * t187 + t184 * t5;
t29 = [0.2e1 * t53 * t91 - 0.2e1 * t54 * t92 + 0.2e1 * t88 * t83 + 0.2e1 * t13 * t79 + 0.2e1 * t14 * t80 + t75 * t8 + t76 * t9 + 0.2e1 * t10 * t77 + 0.2e1 * t11 * t78 + 0.2e1 * t54 * t64 + 0.2e1 * t16 * t65 + 0.2e1 * t62 * t38 + 0.2e1 * t23 * t48 + 0.2e1 * t25 * t49 + 0.2e1 * t24 * t50 + 0.2e1 * t1 * t41 + 0.2e1 * t2 * t42 + 0.2e1 * t22 * t47 + 0.2e1 * t31 * t39 + 0.2e1 * t6 * t40 + t32 * t27 + t33 * t28 + 0.2e1 * t18 * t15 + 0.2e1 * t4 * t19 + 0.2e1 * t3 * t20 + (t37 + t7 - t35) * t90 + (-t36 + t34) * t89 + (t111 * t290 - Ifges(4,6) * t130 + (-Ifges(6,4) + Ifges(5,5)) * t90 + (Ifges(6,5) - Ifges(5,6)) * t89 - (Ifges(6,1) + (2 * Ifges(4,2)) + Ifges(5,3)) * t291) * t107 + (0.2e1 * Ifges(4,1) * t111 + Ifges(4,5) * t130 - t290 * t291) * t106 + 0.2e1 * (-mrSges(4,1) * t291 + mrSges(4,2) * t111) * t216 - t291 * t229 - t291 * t230 + 0.2e1 * m(5) * (t13 * t25 + t14 * t24 + t54 * t62) + 0.2e1 * m(6) * (t10 * t22 + t11 * t23 + t16 * t31) + (t46 + t26 - t44) * t74 + (t43 - t45) * t73 + 0.2e1 * (t181 * (-mrSges(3,2) * t183 + mrSges(3,3) * t248) + m(3) * (t238 * t181 + (qJ(2) * t252 - t171) * t178)) * t237 + 0.2e1 * m(4) * (t216 * t88 + t53 * t67 - t54 * t66) + (t1 * t4 + t18 * t6 + t2 * t3) * t284 + t130 * t239 - 0.2e1 * (mrSges(3,1) * t183 - mrSges(3,3) * t252) * t224 - 0.2e1 * t66 * t258 - 0.2e1 * t67 * t266; t118 * t20 - t200 * t19 + t182 * t83 + t86 * t41 + t87 * t42 + t269 * t131 + t267 * t117 + (t15 + t270) * t132 + (-t40 + t268) * t116 + m(7) * (-t1 * t200 - t116 * t18 + t118 * t2 + t132 * t6 + t3 * t87 + t4 * t86) + m(5) * (-t116 * t25 - t117 * t24 + t13 * t132 - t131 * t14) + m(6) * (-t10 * t132 + t11 * t131 + t116 * t22 + t117 * t23) + (-t186 * t266 + (-t38 - t39 - t258) * t189 + (t189 * t91 + (t64 + t65 - t92) * t186) * qJD(3) + m(5) * (t236 * t62 - t255) + m(6) * (-t16 * t189 + t236 * t31) + m(4) * (t182 * t224 + t186 * t53 + t235 * t67 - t236 * t66 - t255)) * t179; 0.2e1 * m(7) * (t118 * t87 - t200 * t86 - t93) + 0.2e1 * (m(6) + m(5)) * (-t179 ^ 2 * t186 * t235 + t117 * t131 - t93); (-t139 / 0.2e1 + t144 / 0.2e1 + t94 / 0.2e1) * t90 + (t138 / 0.2e1 - t142 / 0.2e1) * t89 + (t254 - t282) * t54 + t150 * t39 + t16 * t152 + t164 * t15 + t31 * t136 + t62 * t137 + t2 * t146 + t1 * t147 - t148 * t40 + t6 * t133 + t4 * t121 + t3 * t122 + t32 * t127 / 0.2e1 + t129 * t65 + t114 * t19 + t113 * t20 + t18 * t108 + t84 * t41 + t85 * t42 - t53 * mrSges(4,2) - pkin(3) * t38 + t239 + (t155 / 0.2e1 - t159 / 0.2e1) * t73 + (-t156 / 0.2e1 + t161 / 0.2e1 + t126 / 0.2e1) * t74 - (t175 / 0.2e1 + t174 / 0.2e1) * t291 + (-t35 / 0.2e1 + t37 / 0.2e1 + t7 / 0.2e1 + t11 * mrSges(6,1) - t14 * mrSges(5,3) + (-Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * t107 + (t22 * mrSges(6,1) - t25 * mrSges(5,3) + t256 / 0.2e1 + t43 / 0.2e1 - t45 / 0.2e1 - t201) * qJD(4) + (t268 * qJD(4) + m(5) * (-t25 * qJD(4) - t14) + m(6) * (t22 * qJD(4) + t11) + t269) * pkin(10)) * t185 + (-t34 / 0.2e1 + t36 / 0.2e1 + t13 * mrSges(5,3) - t10 * mrSges(6,1) + t9 * t273 + t8 * t272 + (-Ifges(6,5) / 0.2e1 + Ifges(5,6) / 0.2e1) * t107 + (t184 * t27 / 0.2e1 + t28 * t272) * qJD(6) + (t257 / 0.2e1 - t44 / 0.2e1 + t46 / 0.2e1 + t26 / 0.2e1 + t23 * mrSges(6,1) - t24 * mrSges(5,3)) * qJD(4) + (t267 * qJD(4) + m(5) * (-t24 * qJD(4) + t13) + m(6) * (t23 * qJD(4) - t10) + t270) * pkin(10)) * t188 + m(7) * (t1 * t114 + t113 * t2 - t148 * t18 + t164 * t6 + t3 * t85 + t4 * t84) + t76 * t279 + t95 * t280 + t128 * t281 + m(6) * (t129 * t31 + t150 * t16); t132 * t108 - t116 * t133 + t118 * t122 - t200 * t121 + t87 * t146 + t86 * t147 + m(7) * (t113 * t87 + t114 * t86 - t116 * t164 + t118 * t85 - t132 * t148 - t200 * t84) + 0.2e1 * (t283 + m(5) / 0.2e1) * (t131 * t233 - t132 * t234 + t204) * pkin(10) + (mrSges(5,3) + mrSges(6,1)) * ((t131 * t188 - t132 * t185) * qJD(4) + t204) + ((-t136 - t137) * t189 + (-mrSges(4,2) * t189 + (t152 + t254) * t186) * qJD(3) + 0.2e1 * (-t129 * t189 + t150 * t236) * t283 - t236 * t282) * t179; 0.2e1 * t150 * t136 + 0.2e1 * t164 * t108 - 0.2e1 * pkin(3) * t137 + 0.2e1 * t85 * t146 + 0.2e1 * t84 * t147 - 0.2e1 * t148 * t133 + 0.2e1 * t114 * t121 + 0.2e1 * t113 * t122 + (t113 * t85 + t114 * t84 - t148 * t164) * t284 + 0.2e1 * (m(6) * t150 + t152) * t129 + (-t139 + t144 + t94 + (t127 * t187 + t128 * t184 + t155 - t159) * qJD(4)) * t185 + (-t184 * t96 - t187 * t95 - t138 + t142 + (t127 * t184 - t128 * t187) * qJD(6) + (t126 - t156 + t161) * qJD(4)) * t188; (-t77 + t40) * qJD(5) + (-t47 + t15) * qJ(5) + m(7) * (qJ(5) * t6 + qJD(5) * t18 - t1 * t243 - t2 * t240) + t229 + t230 + (t9 / 0.2e1 - t2 * mrSges(7,3) - t278 * t20) * t187 + (-t8 / 0.2e1 - t1 * mrSges(7,3) - t278 * t19) * t184 + t6 * t154 + t74 * t275 + t32 * t274 + t160 * t281 + t18 * t135 + t90 * t276 - t141 * t280 - t76 * t143 / 0.2e1 - pkin(4) * t48 - t13 * mrSges(5,2) + t14 * mrSges(5,1) - t10 * mrSges(6,3) + t11 * mrSges(6,2) + m(6) * (-pkin(4) * t11 - qJ(5) * t10 - qJD(5) * t22) + (t215 * mrSges(7,3) - (-m(7) * t215 + t209) * t278 + t201) * qJD(6); t132 * t135 + (-mrSges(5,1) + mrSges(6,2)) * t117 + (mrSges(5,2) + t289) * t116 + m(6) * (-pkin(4) * t117 + t202) + m(7) * t202 + t278 * t192 + t286 * mrSges(7,3); t185 * t276 + t187 * t279 + t95 * t273 - t148 * t154 + t164 * t135 + qJD(5) * t133 + qJ(5) * t108 + t175 + t174 - t121 * t243 - t122 * t240 + m(7) * (-qJ(5) * t148 + qJD(5) * t164 - t240 * t85 - t243 * t84) + t208 * mrSges(7,3) + (qJD(5) * t228 - t141 * t272 - t143 * t273) * t188 + ((-t127 / 0.2e1 - t114 * mrSges(7,3) - t188 * t160 / 0.2e1) * t187 + (t113 * mrSges(7,3) - t128 / 0.2e1 + t188 * t274) * t184 - t285 * t278) * qJD(6) + ((-pkin(4) * mrSges(6,1) - Ifges(6,4) + t275) * t188 + (t244 / 0.2e1 - qJ(5) * mrSges(6,1) - Ifges(5,6) + t241 / 0.2e1) * t185 + (m(6) * t210 + t152 + t153) * pkin(10)) * qJD(4); 0.2e1 * qJ(5) * t135 + t141 * t184 - t143 * t187 + (-t241 - t244) * qJD(6) + 0.2e1 * ((m(6) + m(7)) * qJ(5) - t289) * qJD(5); t184 * t19 + t187 * t20 + t209 * qJD(6) + m(7) * (-qJD(6) * t215 + t1 * t184 + t187 * t2) + m(6) * t11 + t48; m(6) * t117 - t192; -m(7) * t208 + qJD(6) * t285 + t184 * t121 + t187 * t122 + t228 * t233; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t7; mrSges(7,1) * t87 - mrSges(7,2) * t86; mrSges(7,1) * t85 - mrSges(7,2) * t84 + t94; ((mrSges(7,2) * t278 - Ifges(7,6)) * t187 + (mrSges(7,1) * t278 - Ifges(7,5)) * t184) * qJD(6); -t154 * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t29(1) t29(2) t29(4) t29(7) t29(11) t29(16); t29(2) t29(3) t29(5) t29(8) t29(12) t29(17); t29(4) t29(5) t29(6) t29(9) t29(13) t29(18); t29(7) t29(8) t29(9) t29(10) t29(14) t29(19); t29(11) t29(12) t29(13) t29(14) t29(15) t29(20); t29(16) t29(17) t29(18) t29(19) t29(20) t29(21);];
Mq  = res;
