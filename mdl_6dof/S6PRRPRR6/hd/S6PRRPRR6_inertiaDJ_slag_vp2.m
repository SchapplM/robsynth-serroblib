% Calculate time derivative of joint inertia matrix for
% S6PRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR6_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_inertiaDJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR6_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR6_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR6_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:22:59
% EndTime: 2019-03-08 22:23:10
% DurationCPUTime: 4.73s
% Computational Cost: add. (7085->546), mult. (19907->832), div. (0->0), fcn. (20414->14), ass. (0->227)
t179 = sin(pkin(13));
t182 = cos(pkin(13));
t186 = sin(qJ(5));
t190 = cos(qJ(5));
t149 = t179 * t186 - t190 * t182;
t144 = t149 * qJD(5);
t150 = t179 * t190 + t182 * t186;
t189 = cos(qJ(6));
t185 = sin(qJ(6));
t222 = qJD(6) * t185;
t194 = t189 * t144 + t150 * t222;
t279 = t185 / 0.2e1;
t257 = t189 / 0.2e1;
t183 = cos(pkin(7));
t191 = cos(qJ(3));
t192 = cos(qJ(2));
t230 = t191 * t192;
t187 = sin(qJ(3));
t188 = sin(qJ(2));
t234 = t187 * t188;
t278 = t183 * t230 - t234;
t180 = sin(pkin(7));
t237 = t180 * t187;
t143 = t179 * t183 + t182 * t237;
t236 = t180 * t191;
t147 = t183 * t187 * pkin(2) + pkin(9) * t236;
t130 = qJ(4) * t183 + t147;
t131 = (-pkin(3) * t191 - qJ(4) * t187 - pkin(2)) * t180;
t90 = -t130 * t179 + t182 * t131;
t61 = -pkin(4) * t236 - pkin(10) * t143 + t90;
t141 = -t179 * t237 + t182 * t183;
t91 = t182 * t130 + t179 * t131;
t73 = pkin(10) * t141 + t91;
t251 = t186 * t61 + t190 * t73;
t226 = qJD(3) * t180;
t235 = t182 * t191;
t121 = (-qJD(4) * t187 + (pkin(3) * t187 - qJ(4) * t191) * qJD(3)) * t180;
t225 = qJD(3) * t187;
t214 = t180 * t225;
t256 = pkin(2) * t191;
t218 = t183 * t256;
t135 = -pkin(9) * t214 + qJD(3) * t218;
t126 = qJD(4) * t183 + t135;
t79 = t182 * t121 - t126 * t179;
t57 = (pkin(4) * t187 - pkin(10) * t235) * t226 + t79;
t213 = t191 * t226;
t208 = t179 * t213;
t80 = t179 * t121 + t182 * t126;
t72 = -pkin(10) * t208 + t80;
t11 = -qJD(5) * t251 - t186 * t72 + t190 * t57;
t163 = -mrSges(7,1) * t189 + mrSges(7,2) * t185;
t277 = -m(7) * pkin(5) - mrSges(6,1) + t163;
t276 = 2 * m(5);
t275 = 2 * m(6);
t274 = 0.2e1 * m(7);
t273 = -2 * mrSges(4,3);
t272 = -2 * mrSges(6,3);
t254 = pkin(10) + qJ(4);
t162 = t254 * t182;
t210 = qJD(5) * t254;
t219 = t190 * qJD(4);
t220 = t186 * qJD(4);
t223 = qJD(5) * t190;
t93 = t182 * t220 + t162 * t223 + (-t186 * t210 + t219) * t179;
t271 = 0.2e1 * t93;
t211 = t179 * t254;
t116 = t186 * t162 + t190 * t211;
t270 = 0.2e1 * t116;
t95 = t141 * t186 + t143 * t190;
t197 = t185 * t236 - t189 * t95;
t199 = t190 * t141 - t143 * t186;
t69 = qJD(5) * t199 - t149 * t213;
t37 = qJD(6) * t197 - t185 * t69 + t189 * t214;
t269 = t37 / 0.2e1;
t145 = t150 * qJD(5);
t221 = qJD(6) * t189;
t195 = -t185 * t144 + t150 * t221;
t46 = -Ifges(7,1) * t194 - Ifges(7,4) * t195 + Ifges(7,5) * t145;
t268 = t46 / 0.2e1;
t83 = -t185 * t95 - t189 * t236;
t267 = t83 / 0.2e1;
t247 = Ifges(7,4) * t185;
t204 = Ifges(7,1) * t189 - t247;
t87 = Ifges(7,5) * t149 + t150 * t204;
t266 = t87 / 0.2e1;
t265 = -t150 / 0.2e1;
t176 = Ifges(7,5) * t221;
t264 = -Ifges(7,6) * t222 / 0.2e1 + t176 / 0.2e1;
t156 = t204 * qJD(6);
t263 = t156 / 0.2e1;
t262 = Ifges(7,5) * t279 + Ifges(7,6) * t257;
t165 = Ifges(7,2) * t189 + t247;
t261 = -t165 / 0.2e1;
t246 = Ifges(7,4) * t189;
t166 = Ifges(7,1) * t185 + t246;
t260 = t166 / 0.2e1;
t259 = t182 / 0.2e1;
t258 = -t185 / 0.2e1;
t181 = sin(pkin(6));
t184 = cos(pkin(6));
t232 = t188 * t191;
t233 = t187 * t192;
t196 = t183 * t233 + t232;
t109 = t181 * t196 + t184 * t237;
t142 = -t180 * t181 * t192 + t184 * t183;
t81 = -t109 * t179 + t142 * t182;
t82 = t109 * t182 + t142 * t179;
t40 = t186 * t81 + t190 * t82;
t227 = qJD(2) * t181;
t215 = t188 * t227;
t209 = t180 * t215;
t76 = t184 * t213 + (t278 * qJD(3) + (-t183 * t234 + t230) * qJD(2)) * t181;
t58 = -t179 * t76 + t182 * t209;
t59 = t179 * t209 + t182 * t76;
t16 = qJD(5) * t40 + t186 * t59 - t190 * t58;
t39 = t186 * t82 - t190 * t81;
t255 = t16 * t39;
t36 = qJD(6) * t83 + t185 * t214 + t189 * t69;
t12 = -mrSges(7,1) * t37 + mrSges(7,2) * t36;
t51 = mrSges(6,1) * t214 - mrSges(6,3) * t69;
t253 = t12 - t51;
t41 = -mrSges(7,1) * t83 - mrSges(7,2) * t197;
t89 = -mrSges(6,1) * t236 - mrSges(6,3) * t95;
t252 = t41 - t89;
t250 = mrSges(7,3) * t150;
t249 = Ifges(5,4) * t179;
t248 = Ifges(5,4) * t182;
t245 = Ifges(7,6) * t185;
t108 = -t181 * t278 - t184 * t236;
t75 = t184 * t214 + (t196 * qJD(3) + (t183 * t232 + t233) * qJD(2)) * t181;
t49 = t108 * t75;
t244 = t116 * t93;
t243 = -mrSges(5,1) * t182 + mrSges(5,2) * t179 - mrSges(4,1);
t122 = t182 * mrSges(5,2) * t213 + mrSges(5,1) * t208;
t70 = qJD(5) * t95 + t150 * t213;
t33 = t70 * mrSges(6,1) + t69 * mrSges(6,2);
t242 = t33 + t122;
t241 = -mrSges(4,1) * t183 - mrSges(5,1) * t141 + mrSges(5,2) * t143 + mrSges(4,3) * t237;
t175 = -pkin(4) * t182 - pkin(3);
t103 = pkin(5) * t149 - pkin(11) * t150 + t175;
t117 = t190 * t162 - t186 * t211;
t67 = t103 * t189 - t117 * t185;
t240 = qJD(6) * t67;
t68 = t103 * t185 + t117 * t189;
t239 = qJD(6) * t68;
t136 = t147 * qJD(3);
t238 = t108 * t136;
t229 = -Ifges(6,5) * t144 - Ifges(6,6) * t145;
t224 = qJD(5) * t186;
t5 = Ifges(7,5) * t36 + Ifges(7,6) * t37 + Ifges(7,3) * t70;
t217 = -Ifges(6,5) * t69 + Ifges(6,6) * t70 - Ifges(6,3) * t214;
t98 = t145 * mrSges(6,1) - t144 * mrSges(6,2);
t25 = -pkin(11) * t236 + t251;
t171 = pkin(9) * t237;
t133 = t171 + (-pkin(3) - t256) * t183;
t96 = -pkin(4) * t141 + t133;
t38 = -pkin(5) * t199 - pkin(11) * t95 + t96;
t13 = -t185 * t25 + t189 * t38;
t115 = pkin(4) * t208 + t136;
t19 = pkin(5) * t70 - pkin(11) * t69 + t115;
t10 = t186 * t57 + t190 * t72 + t61 * t223 - t224 * t73;
t8 = pkin(11) * t214 + t10;
t1 = qJD(6) * t13 + t185 * t19 + t189 * t8;
t14 = t185 * t38 + t189 * t25;
t2 = -qJD(6) * t14 - t185 * t8 + t189 * t19;
t207 = t1 * t189 - t2 * t185;
t206 = t116 * t16 + t39 * t93;
t205 = mrSges(7,1) * t185 + mrSges(7,2) * t189;
t203 = -Ifges(7,2) * t185 + t246;
t202 = -t179 * t58 + t182 * t59;
t31 = -t186 * t73 + t190 * t61;
t22 = t108 * t189 - t185 * t40;
t23 = t108 * t185 + t189 * t40;
t27 = -Ifges(7,4) * t197 + Ifges(7,2) * t83 - Ifges(7,6) * t199;
t28 = -Ifges(7,1) * t197 + Ifges(7,4) * t83 - Ifges(7,5) * t199;
t198 = t257 * t28 + t258 * t27;
t44 = -Ifges(7,5) * t194 - Ifges(7,6) * t195 + Ifges(7,3) * t145;
t169 = Ifges(4,5) * t213;
t155 = t203 * qJD(6);
t153 = t205 * qJD(6);
t152 = -mrSges(4,2) * t183 + mrSges(4,3) * t236;
t146 = -t171 + t218;
t134 = (mrSges(4,1) * t187 + mrSges(4,2) * t191) * t226;
t128 = (mrSges(5,1) * t187 - mrSges(5,3) * t235) * t226;
t127 = (-mrSges(5,3) * t179 * t191 - mrSges(5,2) * t187) * t226;
t120 = -mrSges(5,1) * t236 - mrSges(5,3) * t143;
t119 = mrSges(5,2) * t236 + mrSges(5,3) * t141;
t112 = Ifges(6,1) * t150 - Ifges(6,4) * t149;
t111 = Ifges(6,4) * t150 - Ifges(6,2) * t149;
t110 = mrSges(6,1) * t149 + mrSges(6,2) * t150;
t107 = (t187 * Ifges(5,5) + (t182 * Ifges(5,1) - t249) * t191) * t226;
t106 = (t187 * Ifges(5,6) + (-t179 * Ifges(5,2) + t248) * t191) * t226;
t105 = mrSges(7,1) * t149 - t189 * t250;
t104 = -mrSges(7,2) * t149 - t185 * t250;
t102 = pkin(5) * t145 + pkin(11) * t144;
t101 = t205 * t150;
t100 = -Ifges(6,1) * t144 - Ifges(6,4) * t145;
t99 = -Ifges(6,4) * t144 - Ifges(6,2) * t145;
t92 = t182 * t219 - t162 * t224 + (-t190 * t210 - t220) * t179;
t88 = mrSges(6,2) * t236 + mrSges(6,3) * t199;
t86 = Ifges(7,6) * t149 + t150 * t203;
t85 = Ifges(7,3) * t149 + (Ifges(7,5) * t189 - t245) * t150;
t78 = -mrSges(7,2) * t145 - mrSges(7,3) * t195;
t77 = mrSges(7,1) * t145 + mrSges(7,3) * t194;
t55 = mrSges(7,1) * t195 - mrSges(7,2) * t194;
t52 = -mrSges(6,2) * t214 - mrSges(6,3) * t70;
t50 = -mrSges(6,1) * t199 + mrSges(6,2) * t95;
t48 = Ifges(6,1) * t95 + Ifges(6,4) * t199 - Ifges(6,5) * t236;
t47 = Ifges(6,4) * t95 + Ifges(6,2) * t199 - Ifges(6,6) * t236;
t45 = -Ifges(7,4) * t194 - Ifges(7,2) * t195 + Ifges(7,6) * t145;
t43 = -mrSges(7,1) * t199 + mrSges(7,3) * t197;
t42 = mrSges(7,2) * t199 + mrSges(7,3) * t83;
t30 = Ifges(6,1) * t69 - Ifges(6,4) * t70 + Ifges(6,5) * t214;
t29 = Ifges(6,4) * t69 - Ifges(6,2) * t70 + Ifges(6,6) * t214;
t26 = -Ifges(7,5) * t197 + Ifges(7,6) * t83 - Ifges(7,3) * t199;
t24 = pkin(5) * t236 - t31;
t21 = t102 * t189 - t185 * t92 - t239;
t20 = t102 * t185 + t189 * t92 + t240;
t18 = -mrSges(7,2) * t70 + mrSges(7,3) * t37;
t17 = mrSges(7,1) * t70 - mrSges(7,3) * t36;
t15 = -qJD(5) * t39 + t186 * t58 + t190 * t59;
t9 = -pkin(5) * t214 - t11;
t7 = Ifges(7,1) * t36 + Ifges(7,4) * t37 + Ifges(7,5) * t70;
t6 = Ifges(7,4) * t36 + Ifges(7,2) * t37 + Ifges(7,6) * t70;
t4 = -qJD(6) * t23 - t15 * t185 + t189 * t75;
t3 = qJD(6) * t22 + t15 * t189 + t185 * t75;
t32 = [0.2e1 * m(7) * (t22 * t4 + t23 * t3 + t255) + 0.2e1 * m(6) * (t15 * t40 + t255 + t49) + 0.2e1 * m(5) * (t58 * t81 + t59 * t82 + t49) + 0.2e1 * m(4) * (t109 * t76 + t142 * t209 + t49); t59 * t119 + t58 * t120 + t82 * t127 + t81 * t128 + t142 * t134 + t15 * t88 + t76 * t152 + t22 * t17 + t23 * t18 + t3 * t42 + t4 * t43 + t40 * t52 + t253 * t39 + t252 * t16 + t242 * t108 + (-mrSges(3,1) * t188 - mrSges(3,2) * t192) * t227 + (t50 + t241) * t75 + ((-mrSges(4,1) * t191 + mrSges(4,2) * t187) * t209 + (t108 * t191 - t109 * t187) * qJD(3) * mrSges(4,3)) * t180 + m(4) * (-pkin(2) * t180 ^ 2 * t215 + t109 * t135 - t146 * t75 + t147 * t76 + t238) + m(5) * (t133 * t75 + t58 * t90 + t59 * t91 + t79 * t81 + t80 * t82 + t238) + m(6) * (t10 * t40 + t108 * t115 - t11 * t39 + t15 * t251 - t16 * t31 + t75 * t96) + m(7) * (t1 * t23 + t13 * t4 + t14 * t3 + t16 * t24 + t2 * t22 + t39 * t9); (t26 - t47) * t70 - t197 * t7 + (-0.2e1 * pkin(2) * t134 + (-0.2e1 * Ifges(4,4) * t237 + Ifges(5,5) * t143 + Ifges(6,5) * t95 - 0.2e1 * Ifges(4,6) * t183 + Ifges(5,6) * t141 + Ifges(6,6) * t199 + t147 * t273) * t225 + ((t182 * (Ifges(5,1) * t143 + Ifges(5,4) * t141) - t179 * (Ifges(5,4) * t143 + Ifges(5,2) * t141) + t146 * t273 + Ifges(4,5) * t183 + 0.2e1 * (-Ifges(5,5) * t182 + Ifges(5,6) * t179 + Ifges(4,4)) * t236 + ((2 * Ifges(4,1)) - (2 * Ifges(4,2)) - (2 * Ifges(5,3)) - Ifges(6,3)) * t237) * qJD(3) + t217) * t191) * t180 - (t5 - t29) * t199 + 0.2e1 * t241 * t136 + 0.2e1 * m(4) * (t135 * t147 - t136 * t146) + t183 * t169 + (t1 * t14 + t13 * t2 + t24 * t9) * t274 + (t133 * t136 + t79 * t90 + t80 * t91) * t276 + 0.2e1 * t13 * t17 + 0.2e1 * t14 * t18 + 0.2e1 * t24 * t12 + t36 * t28 + t37 * t27 + 0.2e1 * t9 * t41 + 0.2e1 * t1 * t42 + 0.2e1 * t2 * t43 + 0.2e1 * t31 * t51 + t69 * t48 + t83 * t6 + 0.2e1 * t10 * t88 + 0.2e1 * t11 * t89 + t95 * t30 + 0.2e1 * t96 * t33 + 0.2e1 * t115 * t50 + 0.2e1 * t80 * t119 + 0.2e1 * t79 * t120 + 0.2e1 * t91 * t127 + 0.2e1 * t90 * t128 + 0.2e1 * t133 * t122 + (t10 * t251 + t11 * t31 + t115 * t96) * t275 + 0.2e1 * t251 * t52 + t141 * t106 + t143 * t107 + 0.2e1 * t135 * t152; -t76 * mrSges(4,2) + t16 * t101 + t3 * t104 + t4 * t105 + t108 * t98 + t22 * t77 + t23 * t78 + t39 * t55 + t202 * mrSges(5,3) + (t110 + t243) * t75 + m(7) * (t20 * t23 + t21 * t22 + t3 * t68 + t4 * t67 + t206) + m(6) * (t117 * t15 + t175 * t75 + t40 * t92 + t206) + m(5) * (-pkin(3) * t75 + (-t179 * t81 + t182 * t82) * qJD(4) + t202 * qJ(4)) + (-t144 * t39 - t145 * t40 - t149 * t15 + t150 * t16) * mrSges(6,3); t169 - t197 * t268 - (t44 / 0.2e1 - t99 / 0.2e1) * t199 + (-t10 * mrSges(6,3) + t5 / 0.2e1 - t29 / 0.2e1) * t149 + m(7) * (t1 * t68 + t116 * t9 + t13 * t21 + t14 * t20 + t2 * t67 + t24 * t93) + m(5) * (-pkin(3) * t136 + (-t179 * t90 + t182 * t91) * qJD(4) + (-t79 * t179 + t80 * t182) * qJ(4)) + (t85 / 0.2e1 - t111 / 0.2e1) * t70 + t175 * t33 - (-t31 * mrSges(6,3) + t48 / 0.2e1 + t198) * t144 + t36 * t266 + t45 * t267 + t86 * t269 + (-qJD(4) * t120 - qJ(4) * t128 - t79 * mrSges(5,3) + t107 / 0.2e1) * t179 + (qJD(4) * t119 + qJ(4) * t127 + t80 * mrSges(5,3) + t106 / 0.2e1) * t182 + t243 * t136 + t252 * t93 + t253 * t116 + (t7 * t257 + t6 * t258 - t11 * mrSges(6,3) + t30 / 0.2e1 + (t28 * t258 - t189 * t27 / 0.2e1) * qJD(6)) * t150 + (-t191 * t229 / 0.2e1 + ((Ifges(6,5) * t150 / 0.2e1 - Ifges(6,6) * t149 / 0.2e1 - Ifges(4,6) + Ifges(5,5) * t179 / 0.2e1 + Ifges(5,6) * t259) * t187 + ((Ifges(5,1) * t179 + t248) * t259 - t179 * (Ifges(5,2) * t182 + t249) / 0.2e1) * t191) * qJD(3)) * t180 + t20 * t42 + t21 * t43 + t24 * t55 + t67 * t17 + t68 * t18 + t13 * t77 + t14 * t78 + t92 * t88 + t96 * t98 + t95 * t100 / 0.2e1 + t9 * t101 + t1 * t104 + t2 * t105 + t69 * t112 / 0.2e1 + t115 * t110 + t117 * t52 - pkin(3) * t122 + m(6) * (t10 * t117 - t11 * t116 + t115 * t175 + t251 * t92 - t31 * t93) + (-t251 * mrSges(6,3) + t26 / 0.2e1 - t47 / 0.2e1) * t145 - t135 * mrSges(4,2); t101 * t271 + 0.2e1 * t20 * t104 + 0.2e1 * t21 * t105 + t55 * t270 + 0.2e1 * t175 * t98 + 0.2e1 * t67 * t77 + 0.2e1 * t68 * t78 + (t20 * t68 + t21 * t67 + t244) * t274 + (t117 * t92 + t244) * t275 + (t272 * t92 + t44 - t99) * t149 + (t117 * t272 - t111 + t85) * t145 - (mrSges(6,3) * t270 - t185 * t86 + t189 * t87 + t112) * t144 + (mrSges(6,3) * t271 - t185 * t45 + t189 * t46 + t100 + (-t185 * t87 - t189 * t86) * qJD(6)) * t150 + (qJ(4) * t276 + 0.2e1 * mrSges(5,3)) * qJD(4) * (t179 ^ 2 + t182 ^ 2); m(7) * (t185 * t3 + t189 * t4 + (-t185 * t22 + t189 * t23) * qJD(6)) + 0.2e1 * (m(6) / 0.2e1 + m(5) / 0.2e1) * t75; t189 * t17 + t185 * t18 + (-t185 * t43 + t189 * t42) * qJD(6) + m(7) * (t1 * t185 + t189 * t2 + (-t13 * t185 + t14 * t189) * qJD(6)) + m(6) * t115 + m(5) * t136 + t242; m(7) * (t185 * t20 + t189 * t21 + (-t185 * t67 + t189 * t68) * qJD(6)) + t104 * t221 + t185 * t78 - t105 * t222 + t189 * t77 + t98; 0; -t15 * mrSges(6,2) + t39 * t153 + (m(7) * pkin(11) + mrSges(7,3)) * (-t185 * t4 + t189 * t3 + (-t185 * t23 - t189 * t22) * qJD(6)) + t277 * t16; t6 * t257 + t7 * t279 + t9 * t163 + t70 * t262 + t165 * t269 + t36 * t260 - t10 * mrSges(6,2) + t11 * mrSges(6,1) + t24 * t153 - t199 * t264 + t155 * t267 - t197 * t263 + t198 * qJD(6) + (-m(7) * t9 - t12) * pkin(5) + ((-t13 * t189 - t14 * t185) * qJD(6) + t207) * mrSges(7,3) + (m(7) * (-t13 * t221 - t14 * t222 + t207) + t189 * t18 - t185 * t17 - t42 * t222 - t43 * t221) * pkin(11) - t217; t145 * t262 - pkin(5) * t55 - t92 * mrSges(6,2) + t116 * t153 + t149 * t264 + t277 * t93 + (t150 * t263 - t144 * t260 + t20 * mrSges(7,3) + t45 / 0.2e1 + (-t67 * mrSges(7,3) + t150 * t261 + t266) * qJD(6) + (m(7) * (t20 - t240) + t78 - qJD(6) * t105) * pkin(11)) * t189 + (t155 * t265 - t144 * t261 - t21 * mrSges(7,3) + t268 + (-t86 / 0.2e1 + t166 * t265 - t68 * mrSges(7,3)) * qJD(6) + (m(7) * (-t21 - t239) - t77 - qJD(6) * t104) * pkin(11)) * t185 + t229; 0; -0.2e1 * pkin(5) * t153 + t155 * t189 + t156 * t185 + (-t185 * t165 + t189 * t166) * qJD(6); mrSges(7,1) * t4 - mrSges(7,2) * t3; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t5; mrSges(7,1) * t21 - mrSges(7,2) * t20 + t44; -t153; t176 + (pkin(11) * t163 - t245) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t32(1) t32(2) t32(4) t32(7) t32(11) t32(16); t32(2) t32(3) t32(5) t32(8) t32(12) t32(17); t32(4) t32(5) t32(6) t32(9) t32(13) t32(18); t32(7) t32(8) t32(9) t32(10) t32(14) t32(19); t32(11) t32(12) t32(13) t32(14) t32(15) t32(20); t32(16) t32(17) t32(18) t32(19) t32(20) t32(21);];
Mq  = res;
