% Calculate time derivative of joint inertia matrix for
% S6RRPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR11_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR11_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR11_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR11_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR11_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:39:06
% EndTime: 2019-03-09 09:39:14
% DurationCPUTime: 4.60s
% Computational Cost: add. (6730->531), mult. (16678->780), div. (0->0), fcn. (16145->10), ass. (0->221)
t172 = sin(pkin(11));
t174 = cos(pkin(11));
t248 = cos(qJ(5));
t208 = qJD(5) * t248;
t178 = sin(qJ(5));
t218 = qJD(5) * t178;
t129 = -t172 * t208 - t174 * t218;
t134 = t172 * t178 - t174 * t248;
t180 = cos(qJ(6));
t177 = sin(qJ(6));
t217 = qJD(6) * t177;
t183 = t180 * t129 + t134 * t217;
t250 = t177 / 0.2e1;
t249 = t180 / 0.2e1;
t202 = (t172 ^ 2 + t174 ^ 2) * qJD(4);
t175 = cos(pkin(6));
t173 = sin(pkin(6));
t181 = cos(qJ(2));
t227 = t173 * t181;
t128 = -t172 * t227 + t175 * t174;
t179 = sin(qJ(2));
t228 = t173 * t179;
t176 = -pkin(2) - qJ(4);
t204 = -qJ(3) * t179 - pkin(1);
t107 = (t176 * t181 + t204) * t173;
t158 = pkin(8) * t228;
t247 = pkin(1) * t181;
t212 = -pkin(2) - t247;
t98 = pkin(3) * t228 + t158 + (-qJ(4) + t212) * t175;
t64 = -t107 * t172 + t174 * t98;
t47 = pkin(4) * t228 - pkin(9) * t128 + t64;
t127 = -t175 * t172 - t174 * t227;
t65 = t174 * t107 + t172 * t98;
t49 = pkin(9) * t127 + t65;
t242 = t178 * t47 + t248 * t49;
t221 = qJD(2) * t173;
t229 = t172 * t179;
t161 = t175 * t179 * pkin(1);
t258 = pkin(3) + pkin(8);
t100 = -t175 * qJD(4) + (t227 * t258 + t161) * qJD(2);
t211 = t179 * t221;
t152 = pkin(2) * t211;
t219 = qJD(3) * t179;
t80 = t152 + (-t219 - qJD(4) * t181 + (-qJ(3) * t181 + qJ(4) * t179) * qJD(2)) * t173;
t51 = t174 * t100 - t172 * t80;
t35 = (pkin(4) * t181 - pkin(9) * t229) * t221 + t51;
t200 = t174 * t211;
t52 = t172 * t100 + t174 * t80;
t44 = pkin(9) * t200 + t52;
t6 = -qJD(5) * t242 - t178 * t44 + t248 * t35;
t145 = -mrSges(7,1) * t180 + mrSges(7,2) * t177;
t271 = m(7) * pkin(5) + mrSges(6,1) - t145;
t270 = 2 * m(5);
t269 = 2 * m(6);
t268 = 0.2e1 * m(7);
t267 = -0.2e1 * pkin(1);
t266 = 2 * mrSges(4,1);
t265 = -2 * mrSges(3,3);
t264 = -2 * mrSges(6,3);
t245 = -pkin(9) + t176;
t142 = t245 * t172;
t199 = t245 * t248;
t102 = t142 * t178 - t174 * t199;
t263 = 0.2e1 * t102;
t118 = (-pkin(2) * t181 + t204) * t173;
t262 = -0.2e1 * t118;
t130 = -t172 * t218 + t174 * t208;
t216 = qJD(6) * t180;
t184 = -t177 * t129 + t134 * t216;
t39 = Ifges(7,4) * t183 + Ifges(7,2) * t184 + Ifges(7,6) * t130;
t261 = t39 / 0.2e1;
t82 = t178 * t127 + t128 * t248;
t69 = t177 * t228 + t180 * t82;
t260 = t69 / 0.2e1;
t135 = t172 * t248 + t178 * t174;
t238 = Ifges(7,4) * t177;
t196 = Ifges(7,1) * t180 - t238;
t72 = Ifges(7,5) * t135 - t134 * t196;
t259 = t72 / 0.2e1;
t257 = -t134 / 0.2e1;
t164 = Ifges(7,5) * t216;
t256 = -Ifges(7,6) * t217 / 0.2e1 + t164 / 0.2e1;
t237 = Ifges(7,4) * t180;
t195 = -Ifges(7,2) * t177 + t237;
t140 = t195 * qJD(6);
t255 = t140 / 0.2e1;
t254 = Ifges(7,5) * t250 + Ifges(7,6) * t249;
t147 = Ifges(7,2) * t180 + t238;
t253 = t147 / 0.2e1;
t148 = Ifges(7,1) * t177 + t237;
t252 = t148 / 0.2e1;
t251 = t174 / 0.2e1;
t246 = Ifges(4,6) + Ifges(3,4);
t220 = qJD(2) * t181;
t210 = t173 * t220;
t185 = t127 * t248 - t178 * t128;
t62 = qJD(5) * t185 + t135 * t211;
t68 = -t177 * t82 + t180 * t228;
t30 = qJD(6) * t68 + t177 * t210 + t180 * t62;
t31 = -qJD(6) * t69 - t177 * t62 + t180 * t210;
t12 = -mrSges(7,1) * t31 + mrSges(7,2) * t30;
t53 = mrSges(6,1) * t210 - t62 * mrSges(6,3);
t244 = t12 - t53;
t33 = -mrSges(7,1) * t68 + mrSges(7,2) * t69;
t74 = mrSges(6,1) * t228 - mrSges(6,3) * t82;
t243 = t33 - t74;
t241 = mrSges(7,3) * t134;
t240 = Ifges(5,4) * t172;
t239 = Ifges(5,4) * t174;
t236 = Ifges(7,6) * t177;
t205 = t248 * qJD(4);
t207 = t178 * t245;
t215 = t178 * qJD(4);
t76 = t142 * t208 - t172 * t215 + (qJD(5) * t207 + t205) * t174;
t235 = t102 * t76;
t214 = t175 * t247;
t153 = qJD(2) * t214;
t121 = -pkin(8) * t211 + t153;
t234 = t121 * mrSges(3,2);
t103 = t142 * t248 + t174 * t207;
t162 = t172 * pkin(4) + qJ(3);
t88 = pkin(5) * t135 + pkin(10) * t134 + t162;
t56 = -t103 * t177 + t180 * t88;
t233 = qJD(6) * t56;
t57 = t103 * t180 + t177 * t88;
t232 = qJD(6) * t57;
t132 = pkin(8) * t227 + t161;
t122 = t132 * qJD(2);
t231 = t122 * t179;
t230 = t129 * t134;
t225 = Ifges(6,5) * t129 - Ifges(6,6) * t130;
t224 = Ifges(3,5) * t210 + Ifges(4,5) * t211;
t165 = t175 * qJD(3);
t223 = t153 + t165;
t201 = t172 * t211;
t63 = qJD(5) * t82 + t178 * t201 - t200 * t248;
t7 = Ifges(7,5) * t30 + Ifges(7,6) * t31 + Ifges(7,3) * t63;
t213 = Ifges(6,5) * t62 - Ifges(6,6) * t63 + Ifges(6,3) * t210;
t117 = -t175 * qJ(3) - t132;
t27 = t63 * mrSges(6,1) + t62 * mrSges(6,2);
t206 = (mrSges(4,2) - mrSges(3,1)) * t122;
t85 = t130 * mrSges(6,1) + t129 * mrSges(6,2);
t203 = (t177 ^ 2 + t180 ^ 2) * t130;
t108 = pkin(3) * t227 - t117;
t16 = pkin(10) * t228 + t242;
t78 = -pkin(4) * t127 + t108;
t32 = -pkin(5) * t185 - pkin(10) * t82 + t78;
t10 = -t16 * t177 + t180 * t32;
t79 = (-pkin(4) * t174 - t258) * t211 + t223;
t17 = pkin(5) * t63 - pkin(10) * t62 + t79;
t5 = t178 * t35 + t47 * t208 - t218 * t49 + t248 * t44;
t3 = pkin(10) * t210 + t5;
t1 = qJD(6) * t10 + t17 * t177 + t180 * t3;
t11 = t16 * t180 + t177 * t32;
t2 = -qJD(6) * t11 + t17 * t180 - t177 * t3;
t198 = t1 * t180 - t177 * t2;
t197 = mrSges(7,1) * t177 + mrSges(7,2) * t180;
t194 = -t10 * t177 + t11 * t180;
t13 = mrSges(7,1) * t63 - mrSges(7,3) * t30;
t14 = -mrSges(7,2) * t63 + mrSges(7,3) * t31;
t193 = -t177 * t13 + t180 * t14;
t192 = t52 * t172 + t51 * t174;
t36 = mrSges(7,2) * t185 + mrSges(7,3) * t68;
t37 = -mrSges(7,1) * t185 - mrSges(7,3) * t69;
t191 = -t177 * t37 + t180 * t36;
t190 = -t177 * t56 + t180 * t57;
t189 = -t102 * t129 + t134 * t76;
t18 = -t178 * t49 + t248 * t47;
t23 = Ifges(7,4) * t69 + Ifges(7,2) * t68 - Ifges(7,6) * t185;
t24 = Ifges(7,1) * t69 + Ifges(7,4) * t68 - Ifges(7,5) * t185;
t188 = t24 * t249 - t177 * t23 / 0.2e1;
t112 = -mrSges(5,1) * t200 + mrSges(5,2) * t201;
t182 = (-t10 * t180 - t11 * t177) * qJD(6) + t198;
t38 = Ifges(7,5) * t183 + Ifges(7,6) * t184 + Ifges(7,3) * t130;
t144 = mrSges(5,1) * t172 + mrSges(5,2) * t174;
t141 = t196 * qJD(6);
t138 = t197 * qJD(6);
t137 = -mrSges(4,1) * t227 - t175 * mrSges(4,3);
t131 = -t158 + t214;
t120 = t175 * t212 + t158;
t115 = (mrSges(5,3) * t174 * t179 - mrSges(5,2) * t181) * t221;
t114 = (mrSges(5,1) * t181 - mrSges(5,3) * t229) * t221;
t113 = -t121 - t165;
t111 = t152 + (-qJ(3) * t220 - t219) * t173;
t110 = mrSges(5,1) * t228 - mrSges(5,3) * t128;
t109 = -mrSges(5,2) * t228 + mrSges(5,3) * t127;
t106 = -Ifges(6,1) * t134 - Ifges(6,4) * t135;
t105 = -Ifges(6,4) * t134 - Ifges(6,2) * t135;
t104 = mrSges(6,1) * t135 - mrSges(6,2) * t134;
t99 = -t211 * t258 + t223;
t97 = (t181 * Ifges(5,5) + (t172 * Ifges(5,1) + t239) * t179) * t221;
t96 = (t181 * Ifges(5,6) + (t174 * Ifges(5,2) + t240) * t179) * t221;
t95 = mrSges(7,1) * t135 + t180 * t241;
t94 = -mrSges(7,2) * t135 + t177 * t241;
t89 = t197 * t134;
t87 = Ifges(6,1) * t129 - Ifges(6,4) * t130;
t86 = Ifges(6,4) * t129 - Ifges(6,2) * t130;
t84 = -mrSges(5,1) * t127 + mrSges(5,2) * t128;
t83 = pkin(5) * t130 - pkin(10) * t129 + qJD(3);
t75 = -t172 * t205 - t142 * t218 + (qJD(5) * t199 - t215) * t174;
t73 = -mrSges(6,2) * t228 + mrSges(6,3) * t185;
t71 = Ifges(7,6) * t135 - t134 * t195;
t70 = Ifges(7,3) * t135 + (-Ifges(7,5) * t180 + t236) * t134;
t67 = -mrSges(7,2) * t130 + mrSges(7,3) * t184;
t66 = mrSges(7,1) * t130 - mrSges(7,3) * t183;
t55 = -mrSges(7,1) * t184 + mrSges(7,2) * t183;
t54 = -mrSges(6,2) * t210 - t63 * mrSges(6,3);
t50 = -mrSges(6,1) * t185 + mrSges(6,2) * t82;
t46 = Ifges(6,1) * t82 + Ifges(6,4) * t185 + Ifges(6,5) * t228;
t45 = Ifges(6,4) * t82 + Ifges(6,2) * t185 + Ifges(6,6) * t228;
t40 = Ifges(7,1) * t183 + Ifges(7,4) * t184 + Ifges(7,5) * t130;
t26 = Ifges(6,1) * t62 - Ifges(6,4) * t63 + Ifges(6,5) * t210;
t25 = Ifges(6,4) * t62 - Ifges(6,2) * t63 + Ifges(6,6) * t210;
t22 = Ifges(7,5) * t69 + Ifges(7,6) * t68 - Ifges(7,3) * t185;
t21 = -t177 * t75 + t180 * t83 - t232;
t20 = t177 * t83 + t180 * t75 + t233;
t15 = -pkin(5) * t228 - t18;
t9 = Ifges(7,1) * t30 + Ifges(7,4) * t31 + Ifges(7,5) * t63;
t8 = Ifges(7,4) * t30 + Ifges(7,2) * t31 + Ifges(7,6) * t63;
t4 = -pkin(5) * t210 - t6;
t19 = [(t1 * t11 + t10 * t2 + t15 * t4) * t268 + (t108 * t99 + t51 * t64 + t52 * t65) * t270 + 0.2e1 * m(3) * (t121 * t132 - t122 * t131) + 0.2e1 * m(4) * (t111 * t118 + t113 * t117 + t120 * t122) + (0.2e1 * t206 + t224 - 0.2e1 * t234) * t175 + (t18 * t6 + t242 * t5 + t78 * t79) * t269 + 0.2e1 * t242 * t54 + (t22 - t45) * t63 + (t179 * t213 + t231 * t266 + 0.2e1 * t111 * (mrSges(4,2) * t181 - mrSges(4,3) * t179) + 0.2e1 * (t121 * t181 + t231) * mrSges(3,3) + ((t120 * t266 + t131 * t265 + mrSges(4,3) * t262 + Ifges(5,5) * t128 + Ifges(6,5) * t82 + Ifges(5,6) * t127 + Ifges(6,6) * t185 + (-(2 * Ifges(4,4)) + Ifges(3,5)) * t175 + (mrSges(3,2) * t267 + 0.2e1 * t181 * t246) * t173) * t181 + (t172 * (Ifges(5,1) * t128 + Ifges(5,4) * t127) + t132 * t265 + t117 * t266 + t174 * (Ifges(5,4) * t128 + Ifges(5,2) * t127) + mrSges(4,2) * t262 + (Ifges(4,5) - (2 * Ifges(3,6))) * t175 + (mrSges(3,1) * t267 + 0.2e1 * (Ifges(5,5) * t172 + Ifges(5,6) * t174 - t246) * t179) * t173 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) + (2 * Ifges(4,2)) - (2 * Ifges(4,3)) + (2 * Ifges(5,3)) + Ifges(6,3)) * t227) * t179) * qJD(2)) * t173 - (t7 - t25) * t185 + 0.2e1 * t10 * t13 + 0.2e1 * t11 * t14 + 0.2e1 * t15 * t12 + t30 * t24 + t31 * t23 + 0.2e1 * t4 * t33 + 0.2e1 * t1 * t36 + 0.2e1 * t2 * t37 + 0.2e1 * t18 * t53 + t62 * t46 + t68 * t8 + t69 * t9 + 0.2e1 * t5 * t73 + 0.2e1 * t6 * t74 + 0.2e1 * t78 * t27 + 0.2e1 * t79 * t50 + t82 * t26 + 0.2e1 * t99 * t84 + 0.2e1 * t52 * t109 + 0.2e1 * t51 * t110 + 0.2e1 * t108 * t112 + 0.2e1 * t64 * t114 + 0.2e1 * t65 * t115 + t127 * t96 + t128 * t97 + 0.2e1 * t113 * t137; (t50 + t84 - t137) * qJD(3) + t224 + t206 + (-t5 * mrSges(6,3) + t7 / 0.2e1 - t25 / 0.2e1) * t135 + m(4) * (-pkin(2) * t122 - qJ(3) * t113 - qJD(3) * t117) + m(7) * (t1 * t57 + t10 * t21 + t102 * t4 + t11 * t20 + t15 * t76 + t2 * t56) + t30 * t259 + t40 * t260 + t68 * t261 + (t70 / 0.2e1 - t105 / 0.2e1) * t63 + m(6) * (qJD(3) * t78 - t102 * t6 + t103 * t5 + t162 * t79 - t18 * t76 + t242 * t75) + (-t242 * mrSges(6,3) + t22 / 0.2e1 - t45 / 0.2e1) * t130 - t234 + (-t18 * mrSges(6,3) + t46 / 0.2e1 + t188) * t129 + (t176 * t115 - qJD(4) * t109 - t52 * mrSges(5,3) - t96 / 0.2e1) * t172 + (-qJD(4) * t110 - t51 * mrSges(5,3) + t176 * t114 + t97 / 0.2e1) * t174 + m(5) * (qJ(3) * t99 + qJD(3) * t108 + t192 * t176 + (-t172 * t65 - t174 * t64) * qJD(4)) - (t38 / 0.2e1 - t86 / 0.2e1) * t185 + t20 * t36 + t21 * t37 + t15 * t55 + t56 * t13 + t57 * t14 + t10 * t66 + t11 * t67 + t31 * t71 / 0.2e1 + t75 * t73 + t78 * t85 + t82 * t87 / 0.2e1 - t4 * t89 + t1 * t94 + t2 * t95 + t103 * t54 + t79 * t104 + t62 * t106 / 0.2e1 + qJ(3) * t112 - t113 * mrSges(4,3) + t99 * t144 + t162 * t27 + t243 * t76 + t244 * t102 + (-t180 * t9 / 0.2e1 + t8 * t250 + t6 * mrSges(6,3) - t26 / 0.2e1 + (t23 * t249 + t24 * t250) * qJD(6)) * t134 + (t179 * t225 / 0.2e1 + ((Ifges(5,5) * t251 - Ifges(5,6) * t172 / 0.2e1 + Ifges(6,5) * t257 - Ifges(6,6) * t135 / 0.2e1 - Ifges(4,4) - pkin(2) * mrSges(4,1)) * t181 + (-Ifges(3,6) - qJ(3) * mrSges(4,1) + t172 * (Ifges(5,1) * t174 - t240) / 0.2e1 + (-Ifges(5,2) * t172 + t239) * t251) * t179) * qJD(2)) * t173; t55 * t263 + 0.2e1 * t162 * t85 + 0.2e1 * t20 * t94 + 0.2e1 * t21 * t95 + 0.2e1 * t56 * t66 + 0.2e1 * t57 * t67 - 0.2e1 * t76 * t89 + (t264 * t75 + t38 - t86) * t135 + (t103 * t264 - t105 + t70) * t130 + (qJ(3) * qJD(3) - t176 * t202) * t270 + (qJD(3) * t162 + t103 * t75 + t235) * t269 + (t20 * t57 + t21 * t56 + t235) * t268 + (mrSges(6,3) * t263 - t177 * t71 + t180 * t72 + t106) * t129 + (t76 * t264 + t177 * t39 - t180 * t40 - t87 + (t177 * t72 + t180 * t71) * qJD(6)) * t134 + 0.2e1 * mrSges(5,3) * t202 + 0.2e1 * (m(4) * qJ(3) + mrSges(4,3) + t104 + t144) * qJD(3); mrSges(4,1) * t210 + t174 * t114 + t172 * t115 + t244 * t134 - t243 * t129 + (t191 + t73) * t130 + (t54 + (-t177 * t36 - t180 * t37) * qJD(6) + t193) * t135 + m(7) * (-t129 * t15 + t130 * t194 + t134 * t4 + t135 * t182) + m(6) * (t18 * t129 + t130 * t242 - t134 * t6 + t5 * t135) + m(5) * t192 + m(4) * t122; t134 * t55 + (-t177 * t95 + t180 * t94) * t130 + (0.2e1 * t134 * mrSges(6,3) + t89) * t129 + (t130 * t264 - t177 * t66 + t180 * t67 + (-t177 * t94 - t180 * t95) * qJD(6)) * t135 + m(7) * (t190 * t130 + (-t177 * t21 + t180 * t20 + (-t177 * t57 - t180 * t56) * qJD(6)) * t135 + t189) + m(6) * (t103 * t130 + t135 * t75 + t189) - m(5) * t202; 0.2e1 * m(6) * (t130 * t135 - t230) + 0.2e1 * m(7) * (t135 * t203 - t230); t180 * t13 + t177 * t14 + t191 * qJD(6) + m(7) * (qJD(6) * t194 + t1 * t177 + t180 * t2) + m(6) * t79 + m(5) * t99 + t112 + t27; m(7) * (qJD(6) * t190 + t177 * t20 + t180 * t21) + t94 * t216 + t177 * t67 - t95 * t217 + t180 * t66 + (m(6) + m(5)) * qJD(3) + t85; 0; 0; -t5 * mrSges(6,2) + t6 * mrSges(6,1) + t15 * t138 - t185 * t256 + t68 * t255 + t141 * t260 + t4 * t145 + t63 * t254 + t31 * t253 + t30 * t252 + t9 * t250 + t8 * t249 + t188 * qJD(6) + (-m(7) * t4 - t12) * pkin(5) + t182 * mrSges(7,3) + (-t37 * t216 - t36 * t217 + m(7) * (-t10 * t216 - t11 * t217 + t198) + t193) * pkin(10) + t213; -pkin(5) * t55 - t75 * mrSges(6,2) + t102 * t138 + t135 * t256 + t130 * t254 - t271 * t76 + (t141 * t257 + t129 * t252 + t20 * mrSges(7,3) + t261 + (-t56 * mrSges(7,3) + t134 * t253 + t259) * qJD(6) + (m(7) * (t20 - t233) + t67 - qJD(6) * t95) * pkin(10)) * t180 + (t134 * t255 - t129 * t147 / 0.2e1 - t21 * mrSges(7,3) + t40 / 0.2e1 + (t134 * t252 - t57 * mrSges(7,3) - t71 / 0.2e1) * qJD(6) + (-qJD(6) * t94 + m(7) * (-t21 - t232) - t66) * pkin(10)) * t177 + t225; -t130 * mrSges(6,2) + t134 * t138 + (m(7) * pkin(10) + mrSges(7,3)) * t203 + t271 * t129; 0; -0.2e1 * pkin(5) * t138 + t140 * t180 + t141 * t177 + (-t147 * t177 + t148 * t180) * qJD(6); mrSges(7,1) * t2 - mrSges(7,2) * t1 + t7; mrSges(7,1) * t21 - mrSges(7,2) * t20 + t38; (-t130 * t180 + t135 * t217) * mrSges(7,2) + (-t130 * t177 - t135 * t216) * mrSges(7,1); -t138; t164 + (pkin(10) * t145 - t236) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t19(1) t19(2) t19(4) t19(7) t19(11) t19(16); t19(2) t19(3) t19(5) t19(8) t19(12) t19(17); t19(4) t19(5) t19(6) t19(9) t19(13) t19(18); t19(7) t19(8) t19(9) t19(10) t19(14) t19(19); t19(11) t19(12) t19(13) t19(14) t19(15) t19(20); t19(16) t19(17) t19(18) t19(19) t19(20) t19(21);];
Mq  = res;
