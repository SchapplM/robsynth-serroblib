% Calculate time derivative of joint inertia matrix for
% S6RRPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:01:16
% EndTime: 2019-03-09 09:01:25
% DurationCPUTime: 4.14s
% Computational Cost: add. (6067->461), mult. (16549->685), div. (0->0), fcn. (15982->10), ass. (0->217)
t153 = sin(pkin(6));
t265 = 0.2e1 * t153;
t264 = -Ifges(5,4) + Ifges(4,5);
t263 = Ifges(5,5) - Ifges(4,6);
t161 = cos(qJ(2));
t202 = qJD(2) * t153;
t192 = t161 * t202;
t156 = sin(qJ(6));
t159 = cos(qJ(6));
t160 = cos(qJ(5));
t196 = qJD(6) * t160;
t157 = sin(qJ(5));
t201 = qJD(5) * t157;
t259 = t156 * t201 - t159 * t196;
t262 = t156 / 0.2e1;
t238 = t159 / 0.2e1;
t204 = t156 ^ 2 + t159 ^ 2;
t261 = -0.1e1 + t204;
t237 = m(7) * t157;
t260 = t157 ^ 2 - t160 ^ 2;
t152 = sin(pkin(11));
t154 = cos(pkin(11));
t158 = sin(qJ(2));
t103 = (t152 * t161 + t154 * t158) * t153;
t155 = cos(pkin(6));
t246 = pkin(3) + pkin(9);
t236 = pkin(1) * t155;
t140 = t161 * t236;
t229 = -pkin(8) - qJ(3);
t187 = t229 * t158;
t84 = pkin(2) * t155 + t153 * t187 + t140;
t139 = t158 * t236;
t208 = t153 * t161;
t110 = pkin(8) * t208 + t139;
t98 = qJ(3) * t208 + t110;
t63 = -t152 * t98 + t154 * t84;
t41 = pkin(4) * t103 - t155 * t246 - t63;
t209 = t153 * t158;
t102 = t152 * t209 - t154 * t208;
t114 = (-pkin(2) * t161 - pkin(1)) * t153;
t168 = -t103 * qJ(4) + t114;
t52 = t102 * t246 + t168;
t226 = t157 * t41 + t160 * t52;
t258 = qJD(5) * t226;
t100 = qJD(2) * t103;
t193 = t158 * t202;
t101 = -t152 * t193 + t154 * t192;
t134 = pkin(2) * t193;
t170 = -qJ(4) * t101 - qJD(4) * t103 + t134;
t35 = t100 * t246 + t170;
t164 = -qJD(3) * t209 + (t208 * t229 - t139) * qJD(2);
t135 = qJD(2) * t140;
t76 = t135 + (qJD(2) * t187 + qJD(3) * t161) * t153;
t50 = t152 * t76 - t154 * t164;
t36 = pkin(4) * t101 + t50;
t6 = -t157 * t35 + t160 * t36 - t258;
t17 = pkin(10) * t103 + t226;
t173 = t102 * t160 - t155 * t157;
t64 = t152 * t84 + t154 * t98;
t56 = -qJ(4) * t155 - t64;
t46 = -pkin(4) * t102 - t56;
t82 = t102 * t157 + t155 * t160;
t25 = -pkin(5) * t173 - pkin(10) * t82 + t46;
t10 = -t156 * t17 + t159 * t25;
t11 = t156 * t25 + t159 * t17;
t200 = qJD(5) * t159;
t61 = qJD(5) * t173 + t100 * t157;
t66 = t103 * t156 + t159 * t82;
t26 = -qJD(6) * t66 + t101 * t159 - t156 * t61;
t65 = t103 * t159 - t156 * t82;
t27 = qJD(6) * t65 + t101 * t156 + t159 * t61;
t12 = -mrSges(7,1) * t26 + mrSges(7,2) * t27;
t44 = mrSges(6,1) * t101 - mrSges(6,3) * t61;
t227 = t44 - t12;
t39 = mrSges(7,2) * t173 + mrSges(7,3) * t65;
t4 = -pkin(5) * t101 - t6;
t40 = -mrSges(7,1) * t173 - mrSges(7,3) * t66;
t68 = -mrSges(6,2) * t103 + mrSges(6,3) * t173;
t257 = m(6) * (t6 + t258) - m(7) * (qJD(5) * t10 * t156 - t11 * t200 + t4) - qJD(5) * (t156 * t40 - t159 * t39 - t68) + t227;
t256 = 0.2e1 * m(7);
t255 = 2 * mrSges(5,1);
t254 = -2 * mrSges(3,3);
t253 = -2 * mrSges(4,3);
t67 = pkin(3) * t102 + t168;
t252 = -0.2e1 * t67;
t251 = 0.2e1 * t114;
t143 = -pkin(2) * t154 - pkin(3);
t141 = -pkin(9) + t143;
t250 = 0.2e1 * t141;
t249 = m(4) * pkin(2);
t222 = Ifges(7,4) * t156;
t128 = Ifges(7,2) * t159 + t222;
t221 = Ifges(7,4) * t159;
t177 = -Ifges(7,2) * t156 + t221;
t71 = -t128 * t196 + (Ifges(7,6) * t160 - t157 * t177) * qJD(5);
t248 = t71 / 0.2e1;
t130 = Ifges(7,1) * t156 + t221;
t178 = Ifges(7,1) * t159 - t222;
t72 = -t130 * t196 + (Ifges(7,5) * t160 - t157 * t178) * qJD(5);
t247 = t72 / 0.2e1;
t199 = qJD(5) * t160;
t245 = t261 * t199 * t237;
t106 = Ifges(7,5) * t157 + t160 * t178;
t244 = t106 / 0.2e1;
t197 = qJD(6) * t159;
t145 = Ifges(7,5) * t197;
t198 = qJD(6) * t156;
t243 = -Ifges(7,6) * t198 / 0.2e1 + t145 / 0.2e1;
t242 = Ifges(7,5) * t262 + Ifges(7,6) * t238;
t241 = t128 / 0.2e1;
t240 = -t156 / 0.2e1;
t239 = -t159 / 0.2e1;
t235 = pkin(5) * t157;
t234 = pkin(5) * t160;
t233 = pkin(10) * t157;
t232 = pkin(10) * t160;
t51 = t152 * t164 + t154 * t76;
t231 = t51 * mrSges(4,2);
t230 = mrSges(5,2) - mrSges(4,1);
t31 = -mrSges(7,1) * t65 + mrSges(7,2) * t66;
t69 = mrSges(6,1) * t103 - mrSges(6,3) * t82;
t228 = t31 - t69;
t225 = mrSges(7,3) * t160;
t224 = Ifges(6,4) * t157;
t223 = Ifges(6,4) * t160;
t220 = Ifges(7,6) * t156;
t219 = t101 * Ifges(6,5);
t218 = t101 * Ifges(6,6);
t217 = t103 * Ifges(6,5);
t216 = t103 * Ifges(6,6);
t215 = t103 * t50;
t107 = -pkin(8) * t193 + t135;
t214 = t107 * mrSges(3,2);
t108 = t110 * qJD(2);
t213 = t108 * mrSges(3,1);
t142 = pkin(2) * t152 + qJ(4);
t111 = t142 - t232 + t235;
t207 = t156 * t157;
t78 = t111 * t159 - t141 * t207;
t212 = qJD(6) * t78;
t206 = t157 * t159;
t79 = t111 * t156 + t141 * t206;
t211 = qJD(6) * t79;
t210 = t141 * t160;
t122 = -mrSges(7,2) * t157 - t156 * t225;
t205 = t159 * t122;
t62 = qJD(5) * t82 - t100 * t160;
t7 = Ifges(7,5) * t27 + Ifges(7,6) * t26 + Ifges(7,3) * t62;
t195 = Ifges(6,5) * t61 - Ifges(6,6) * t62 + Ifges(6,3) * t101;
t47 = -qJD(4) * t155 - t51;
t194 = t230 * t50;
t191 = t156 * t199;
t190 = t159 * t199;
t113 = qJD(4) + (t233 + t234) * qJD(5);
t53 = t113 * t156 + t141 * t190 + t212;
t185 = t53 - t212;
t54 = t113 * t159 - t141 * t191 - t211;
t184 = -t54 - t211;
t183 = -2 * Ifges(4,4) - 2 * Ifges(5,6);
t182 = Ifges(3,5) * t192 + t263 * t100 + t264 * t101;
t32 = -pkin(4) * t100 - t47;
t13 = pkin(5) * t62 - pkin(10) * t61 + t32;
t5 = t157 * t36 + t160 * t35 + t199 * t41 - t201 * t52;
t3 = pkin(10) * t101 + t5;
t1 = qJD(6) * t10 + t13 * t156 + t159 * t3;
t2 = -qJD(6) * t11 + t13 * t159 - t156 * t3;
t181 = t1 * t159 - t156 * t2;
t180 = t160 * mrSges(6,1) - t157 * mrSges(6,2);
t126 = t157 * mrSges(6,1) + t160 * mrSges(6,2);
t125 = -mrSges(7,1) * t159 + mrSges(7,2) * t156;
t179 = mrSges(7,1) * t156 + mrSges(7,2) * t159;
t14 = -mrSges(7,2) * t62 + mrSges(7,3) * t26;
t15 = mrSges(7,1) * t62 - mrSges(7,3) * t27;
t176 = t159 * t14 - t156 * t15;
t18 = -t157 * t52 + t160 * t41;
t21 = Ifges(7,4) * t66 + Ifges(7,2) * t65 - Ifges(7,6) * t173;
t22 = Ifges(7,1) * t66 + Ifges(7,4) * t65 - Ifges(7,5) * t173;
t171 = t21 * t240 + t22 * t238;
t169 = t156 * t196 + t157 * t200;
t166 = -t10 * t197 - t11 * t198 + t181;
t165 = -t156 * t54 + t159 * t53 - t197 * t78 - t198 * t79;
t112 = t179 * t160;
t123 = mrSges(7,1) * t157 - t159 * t225;
t96 = mrSges(7,1) * t199 + mrSges(7,3) * t169;
t97 = -mrSges(7,2) * t199 + mrSges(7,3) * t259;
t163 = qJD(5) * t112 - t122 * t198 - t123 * t197 - t156 * t96 + t159 * t97;
t70 = -Ifges(7,5) * t169 + Ifges(7,6) * t259 + Ifges(7,3) * t199;
t16 = -pkin(5) * t103 - t18;
t45 = -mrSges(6,2) * t101 - mrSges(6,3) * t62;
t162 = t45 + (-t156 * t39 - t159 * t40) * qJD(6) + t228 * qJD(5) + m(7) * (qJD(5) * t16 + t166) + m(6) * (-qJD(5) * t18 + t5) + t176;
t131 = Ifges(6,1) * t160 - t224;
t129 = -Ifges(6,2) * t157 + t223;
t121 = (-Ifges(6,1) * t157 - t223) * qJD(5);
t120 = t178 * qJD(6);
t119 = (-Ifges(6,2) * t160 - t224) * qJD(5);
t118 = t177 * qJD(6);
t116 = t180 * qJD(5);
t115 = t179 * qJD(6);
t109 = -pkin(8) * t209 + t140;
t105 = Ifges(7,6) * t157 + t160 * t177;
t104 = Ifges(7,3) * t157 + (Ifges(7,5) * t159 - t220) * t160;
t90 = t101 * mrSges(5,3);
t89 = t101 * mrSges(4,2);
t85 = mrSges(5,1) * t102 - mrSges(5,3) * t155;
t77 = mrSges(7,1) * t259 + mrSges(7,2) * t169;
t57 = -pkin(3) * t155 - t63;
t55 = -mrSges(6,1) * t173 + mrSges(6,2) * t82;
t49 = pkin(3) * t100 + t170;
t43 = Ifges(6,1) * t82 + Ifges(6,4) * t173 + t217;
t42 = Ifges(6,4) * t82 + Ifges(6,2) * t173 + t216;
t30 = mrSges(6,1) * t62 + mrSges(6,2) * t61;
t29 = Ifges(6,1) * t61 - Ifges(6,4) * t62 + t219;
t28 = Ifges(6,4) * t61 - Ifges(6,2) * t62 + t218;
t20 = Ifges(7,5) * t66 + Ifges(7,6) * t65 - Ifges(7,3) * t173;
t9 = Ifges(7,1) * t27 + Ifges(7,4) * t26 + Ifges(7,5) * t62;
t8 = Ifges(7,4) * t27 + Ifges(7,2) * t26 + Ifges(7,6) * t62;
t19 = [-(t7 - t28) * t173 + 0.2e1 * (-t102 * t51 + t215) * mrSges(4,3) + 0.2e1 * t194 * t155 + (t57 * t255 + t63 * t253 + Ifges(6,5) * t82 + Ifges(6,6) * t173 + t264 * t155 + t102 * t183 + ((2 * Ifges(4,1)) + (2 * Ifges(5,2)) + Ifges(6,3)) * t103) * t101 + 0.2e1 * m(5) * (t47 * t56 + t49 * t67 + t50 * t57) + (t20 - t42) * t62 + ((t107 * t161 + t108 * t158) * mrSges(3,3) + (-mrSges(3,2) * pkin(1) + (Ifges(3,1) - Ifges(3,2)) * t158 + Ifges(3,4) * t161) * t192) * t265 + (mrSges(4,1) * t251 + mrSges(5,2) * t252 + t103 * t183 + t263 * t155 + t64 * t253 + t56 * t255 + 0.2e1 * (Ifges(4,2) + Ifges(5,3)) * t102) * t100 + 0.2e1 * m(6) * (t18 * t6 + t226 * t5 + t32 * t46) + 0.2e1 * t226 * t45 + 0.2e1 * m(3) * (t107 * t110 - t108 * t109) + (t182 - 0.2e1 * t213 - 0.2e1 * t214 - 0.2e1 * t231) * t155 + 0.2e1 * m(4) * (-t50 * t63 + t51 * t64) + ((Ifges(3,5) * t155 + t109 * t254) * t161 + (-0.2e1 * Ifges(3,6) * t155 + t249 * t251 + t110 * t254 + 0.2e1 * pkin(2) * (mrSges(4,1) * t102 + mrSges(4,2) * t103) + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t158) * t265) * t158) * t202 + t89 * t251 + t90 * t252 + t215 * t255 + (t1 * t11 + t10 * t2 + t16 * t4) * t256 + 0.2e1 * t11 * t14 + 0.2e1 * t10 * t15 + 0.2e1 * t16 * t12 + t26 * t21 + t27 * t22 + 0.2e1 * t4 * t31 + 0.2e1 * t1 * t39 + 0.2e1 * t2 * t40 + 0.2e1 * t18 * t44 + 0.2e1 * t46 * t30 + 0.2e1 * t32 * t55 + t61 * t43 + t65 * t8 + t66 * t9 + 0.2e1 * t5 * t68 + 0.2e1 * t6 * t69 + t82 * t29 + 0.2e1 * t47 * t85 + 0.2e1 * t49 * (-mrSges(5,2) * t102 - mrSges(5,3) * t103) + t103 * t195; (t55 - t85) * qJD(4) - (t70 / 0.2e1 - t119 / 0.2e1) * t173 + t182 + m(5) * (-qJD(4) * t56 - t142 * t47 + t143 * t50) + (t104 / 0.2e1 - t129 / 0.2e1) * t62 - t213 - t214 + t194 + ((-t226 * mrSges(6,3) - t42 / 0.2e1 + t20 / 0.2e1 - t216 / 0.2e1) * t160 + (t18 * mrSges(6,3) - t43 / 0.2e1 - t217 / 0.2e1 - t171) * t157 + (t160 * t68 + t228 * t157 + t16 * t237 + m(6) * (-t157 * t18 + t160 * t226)) * t141) * qJD(5) + (-t100 * t142 + t101 * t143) * mrSges(5,1) + (-t100 * t152 - t101 * t154) * pkin(2) * mrSges(4,3) - t231 + t27 * t244 + t66 * t247 + t65 * t248 + (t152 * t51 - t154 * t50) * t249 - t47 * mrSges(5,3) + t53 * t39 + t54 * t40 - t16 * t77 + t78 * t15 + t79 * t14 + t10 * t96 + t11 * t97 + t26 * t105 / 0.2e1 + t4 * t112 + t46 * t116 + t82 * t121 / 0.2e1 + t1 * t122 + t2 * t123 + t32 * t126 + t61 * t131 / 0.2e1 + t142 * t30 - Ifges(3,6) * t193 + m(7) * (t1 * t79 + t10 * t54 + t11 * t53 + t2 * t78 - t210 * t4) + m(6) * (t141 * t157 * t5 + qJD(4) * t46 + t142 * t32 + t210 * t6) + (t141 * t45 - t5 * mrSges(6,3) + t7 / 0.2e1 - t28 / 0.2e1 - t218 / 0.2e1) * t157 + (t9 * t238 + t8 * t240 - t6 * mrSges(6,3) + t29 / 0.2e1 + t219 / 0.2e1 + t227 * t141 + (t21 * t239 + t22 * t240) * qJD(6)) * t160; 0.2e1 * t53 * t122 + 0.2e1 * t79 * t97 + (t53 * t79 + t54 * t78) * t256 + 0.2e1 * t54 * t123 + 0.2e1 * t78 * t96 + 0.2e1 * t142 * t116 + 0.2e1 * (mrSges(5,3) + t126 + (m(5) + m(6)) * t142) * qJD(4) + (-t119 + t70 + (t105 * t156 - t106 * t159 + t112 * t250 - t131) * qJD(5)) * t157 + (t77 * t250 - t156 * t71 + t159 * t72 + t121 + (-t105 * t159 - t106 * t156) * qJD(6) + (-0.2e1 * t141 ^ 2 * t237 + t104 - t129) * qJD(5)) * t160; m(4) * t134 + m(5) * t49 - t100 * t230 - t157 * t257 + t160 * t162 + t89 - t90; -t157 * t77 + (m(7) * (t141 * t260 - t206 * t79 + t207 * t78) - t157 * t205 + t123 * t207) * qJD(5) + (m(7) * t165 + t163) * t160; -0.2e1 * t245; m(5) * t50 + t101 * mrSges(5,1) + t157 * t162 + t160 * t257; (t77 + (m(7) * (-t156 * t78 + t159 * t79) + t205 - t156 * t123) * qJD(5)) * t160 + (m(7) * (-0.2e1 * t141 * t199 + t165) + t163) * t157; -m(7) * t261 * t260 * qJD(5); 0.2e1 * t245; -t5 * mrSges(6,2) + t6 * mrSges(6,1) + t16 * t115 - t173 * t243 + t65 * t118 / 0.2e1 + t66 * t120 / 0.2e1 + t4 * t125 + t62 * t242 + t26 * t241 + t27 * t130 / 0.2e1 + t9 * t262 + t8 * t238 + t171 * qJD(6) + (-m(7) * t4 - t12) * pkin(5) + ((-t10 * t159 - t11 * t156) * qJD(6) + t181) * mrSges(7,3) + (m(7) * t166 - t197 * t40 - t198 * t39 + t176) * pkin(10) + t195; pkin(5) * t77 + (t243 + (-Ifges(6,5) + (-m(7) * pkin(5) - mrSges(6,1) + t125) * t141) * qJD(5)) * t157 + (-t130 * t201 / 0.2e1 + qJD(6) * t244 + t248 + t185 * mrSges(7,3) + (m(7) * t185 - qJD(6) * t123 + t97) * pkin(10)) * t159 + (t201 * t241 - qJD(6) * t105 / 0.2e1 + t247 + t184 * mrSges(7,3) + (m(7) * t184 - qJD(6) * t122 - t96) * pkin(10)) * t156 + (t120 * t238 + t118 * t240 - t141 * t115 + (t128 * t239 + t130 * t240) * qJD(6) + (-mrSges(6,2) * t141 - Ifges(6,6) + t242) * qJD(5)) * t160; t157 * t115 + (t160 * t125 + m(7) * (-t204 * t233 - t234) - t204 * t157 * mrSges(7,3) - t180) * qJD(5); -t160 * t115 + (t157 * t125 + m(7) * (t204 * t232 - t235) + t204 * t225 - t126) * qJD(5); -0.2e1 * pkin(5) * t115 + t118 * t159 + t120 * t156 + (-t128 * t156 + t130 * t159) * qJD(6); mrSges(7,1) * t2 - mrSges(7,2) * t1 + t7; mrSges(7,1) * t54 - mrSges(7,2) * t53 + t70; t77; (t157 * t198 - t190) * mrSges(7,2) + (-t157 * t197 - t191) * mrSges(7,1); t145 + (pkin(10) * t125 - t220) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t19(1) t19(2) t19(4) t19(7) t19(11) t19(16); t19(2) t19(3) t19(5) t19(8) t19(12) t19(17); t19(4) t19(5) t19(6) t19(9) t19(13) t19(18); t19(7) t19(8) t19(9) t19(10) t19(14) t19(19); t19(11) t19(12) t19(13) t19(14) t19(15) t19(20); t19(16) t19(17) t19(18) t19(19) t19(20) t19(21);];
Mq  = res;
