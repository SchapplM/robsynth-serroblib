% Calculate time derivative of joint inertia matrix for
% S5RRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP11_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP11_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP11_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP11_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP11_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:14:28
% EndTime: 2019-12-31 22:14:40
% DurationCPUTime: 4.57s
% Computational Cost: add. (3784->522), mult. (10440->736), div. (0->0), fcn. (9161->8), ass. (0->211)
t178 = sin(qJ(4));
t266 = (Ifges(6,4) + Ifges(5,5)) * t178;
t179 = sin(qJ(3));
t182 = cos(qJ(3));
t230 = qJD(3) * t182;
t210 = t178 * t230;
t181 = cos(qJ(4));
t227 = qJD(4) * t181;
t185 = t179 * t227 + t210;
t177 = cos(pkin(5));
t176 = sin(pkin(5));
t180 = sin(qJ(2));
t240 = t176 * t180;
t119 = t177 * t179 + t182 * t240;
t232 = qJD(2) * t180;
t214 = t176 * t232;
t183 = cos(qJ(2));
t239 = t176 * t183;
t215 = t178 * t239;
t118 = -t177 * t182 + t179 * t240;
t233 = qJD(2) * t176;
t213 = t183 * t233;
t82 = -qJD(3) * t118 + t182 * t213;
t41 = -qJD(4) * t215 + t119 * t227 + t178 * t82 - t181 * t214;
t83 = t119 * t178 + t181 * t239;
t42 = -qJD(4) * t83 + t178 * t214 + t82 * t181;
t81 = qJD(3) * t119 + t179 * t213;
t7 = Ifges(5,5) * t42 - Ifges(5,6) * t41 + Ifges(5,3) * t81;
t8 = Ifges(6,4) * t42 + Ifges(6,2) * t81 + Ifges(6,6) * t41;
t264 = t7 + t8;
t13 = mrSges(5,1) * t41 + mrSges(5,2) * t42;
t122 = t177 * t180 * pkin(1) + pkin(7) * t239;
t105 = pkin(8) * t177 + t122;
t106 = (-pkin(2) * t183 - pkin(8) * t180 - pkin(1)) * t176;
t114 = (pkin(2) * t180 - pkin(8) * t183) * t233;
t164 = pkin(7) * t240;
t255 = pkin(1) * t183;
t121 = t177 * t255 - t164;
t115 = t121 * qJD(2);
t231 = qJD(3) * t179;
t27 = -t105 * t230 - t106 * t231 + t114 * t182 - t179 * t115;
t25 = -pkin(3) * t214 - t27;
t262 = -m(5) * t25 - t13;
t142 = (pkin(3) * t179 - pkin(9) * t182) * qJD(3);
t146 = -pkin(3) * t182 - pkin(9) * t179 - pkin(2);
t226 = qJD(4) * t182;
t229 = qJD(4) * t178;
t62 = pkin(8) * (t178 * t231 - t181 * t226) + t142 * t181 - t146 * t229;
t26 = -t105 * t231 + t106 * t230 + t179 * t114 + t182 * t115;
t24 = pkin(9) * t214 + t26;
t104 = t164 + (-pkin(2) - t255) * t177;
t52 = t118 * pkin(3) - t119 * pkin(9) + t104;
t64 = t182 * t105 + t179 * t106;
t54 = -pkin(9) * t239 + t64;
t253 = t178 * t52 + t181 * t54;
t116 = t122 * qJD(2);
t35 = t81 * pkin(3) - t82 * pkin(9) + t116;
t4 = -qJD(4) * t253 - t178 * t24 + t181 * t35;
t192 = pkin(4) * t181 + qJ(5) * t178;
t225 = qJD(5) * t181;
t261 = qJD(4) * t192 - t225;
t260 = 0.2e1 * m(5);
t259 = 2 * m(6);
t258 = 0.2e1 * pkin(8);
t257 = -2 * mrSges(3,3);
t254 = pkin(8) * t182;
t252 = Ifges(4,4) * t179;
t251 = Ifges(4,4) * t182;
t250 = Ifges(5,4) * t178;
t249 = Ifges(5,4) * t181;
t248 = Ifges(6,5) * t178;
t247 = Ifges(6,5) * t181;
t246 = Ifges(5,6) * t181;
t245 = t115 * mrSges(3,2);
t244 = t116 * mrSges(3,1);
t243 = t116 * mrSges(4,1);
t242 = t116 * mrSges(4,2);
t241 = t146 * t181;
t238 = t178 * t179;
t237 = t179 * t181;
t193 = Ifges(6,3) * t178 + t247;
t107 = -Ifges(6,6) * t182 + t179 * t193;
t194 = -Ifges(5,2) * t178 + t249;
t110 = -Ifges(5,6) * t182 + t179 * t194;
t236 = t107 - t110;
t195 = Ifges(6,1) * t181 + t248;
t111 = -Ifges(6,4) * t182 + t179 * t195;
t196 = Ifges(5,1) * t181 - t250;
t112 = -Ifges(5,5) * t182 + t179 * t196;
t235 = t111 + t112;
t209 = t181 * t230;
t234 = Ifges(5,5) * t209 + Ifges(5,3) * t231;
t103 = t178 * t146 + t181 * t254;
t132 = Ifges(6,4) * t227 + Ifges(6,6) * t229;
t228 = qJD(4) * t179;
t6 = Ifges(6,5) * t42 + Ifges(6,6) * t81 + Ifges(6,3) * t41;
t9 = Ifges(5,4) * t42 - Ifges(5,2) * t41 + Ifges(5,6) * t81;
t224 = t6 / 0.2e1 - t9 / 0.2e1;
t222 = Ifges(4,5) * t82 - Ifges(4,6) * t81 + Ifges(4,3) * t214;
t221 = Ifges(4,6) * t239;
t10 = Ifges(6,1) * t42 + Ifges(6,4) * t81 + Ifges(6,5) * t41;
t11 = Ifges(5,1) * t42 - Ifges(5,4) * t41 + Ifges(5,5) * t81;
t220 = t10 / 0.2e1 + t11 / 0.2e1;
t84 = t119 * t181 - t215;
t29 = Ifges(6,5) * t84 + Ifges(6,6) * t118 + Ifges(6,3) * t83;
t32 = Ifges(5,4) * t84 - Ifges(5,2) * t83 + Ifges(5,6) * t118;
t219 = t29 / 0.2e1 - t32 / 0.2e1;
t33 = Ifges(6,1) * t84 + Ifges(6,4) * t118 + Ifges(6,5) * t83;
t34 = Ifges(5,1) * t84 - Ifges(5,4) * t83 + Ifges(5,5) * t118;
t218 = t33 / 0.2e1 + t34 / 0.2e1;
t149 = -Ifges(6,3) * t181 + t248;
t69 = -t149 * t228 + (Ifges(6,6) * t179 + t182 * t193) * qJD(3);
t152 = Ifges(5,2) * t181 + t250;
t72 = -t152 * t228 + (Ifges(5,6) * t179 + t182 * t194) * qJD(3);
t217 = t69 / 0.2e1 - t72 / 0.2e1;
t154 = Ifges(6,1) * t178 - t247;
t73 = -t154 * t228 + (Ifges(6,4) * t179 + t182 * t195) * qJD(3);
t155 = Ifges(5,1) * t178 + t249;
t74 = -t155 * t228 + (Ifges(5,5) * t179 + t182 * t196) * qJD(3);
t216 = t73 / 0.2e1 + t74 / 0.2e1;
t212 = t178 * t228;
t208 = t107 / 0.2e1 - t110 / 0.2e1;
t207 = t111 / 0.2e1 + t112 / 0.2e1;
t130 = t193 * qJD(4);
t133 = t194 * qJD(4);
t206 = t130 / 0.2e1 - t133 / 0.2e1;
t131 = Ifges(5,5) * t227 - Ifges(5,6) * t229;
t205 = t131 / 0.2e1 + t132 / 0.2e1;
t135 = t195 * qJD(4);
t136 = t196 * qJD(4);
t204 = t135 / 0.2e1 + t136 / 0.2e1;
t203 = t149 / 0.2e1 - t152 / 0.2e1;
t202 = t246 / 0.2e1 - Ifges(6,6) * t181 / 0.2e1 + t266 / 0.2e1;
t201 = t154 / 0.2e1 + t155 / 0.2e1;
t21 = -t81 * mrSges(6,1) + t42 * mrSges(6,2);
t63 = -t179 * t105 + t106 * t182;
t200 = Ifges(6,4) * t209 + Ifges(6,2) * t231 + t185 * Ifges(6,6);
t199 = t214 / 0.2e1;
t53 = pkin(3) * t239 - t63;
t148 = -t181 * mrSges(5,1) + t178 * mrSges(5,2);
t198 = mrSges(5,1) * t178 + mrSges(5,2) * t181;
t147 = -t181 * mrSges(6,1) - t178 * mrSges(6,3);
t197 = mrSges(6,1) * t178 - mrSges(6,3) * t181;
t191 = pkin(4) * t178 - qJ(5) * t181;
t16 = -t178 * t54 + t181 * t52;
t188 = pkin(8) + t191;
t3 = t178 * t35 + t181 * t24 + t52 * t227 - t229 * t54;
t186 = t209 - t212;
t93 = mrSges(6,2) * t209 + (-mrSges(6,1) * qJD(3) - mrSges(6,2) * t229) * t179;
t61 = t178 * t142 + t146 * t227 + (-t178 * t226 - t181 * t231) * pkin(8);
t174 = Ifges(4,5) * t230;
t158 = Ifges(3,5) * t213;
t156 = Ifges(4,1) * t179 + t251;
t153 = Ifges(4,2) * t182 + t252;
t143 = -pkin(3) - t192;
t141 = -mrSges(6,2) * t238 - mrSges(6,3) * t182;
t140 = mrSges(6,1) * t182 + mrSges(6,2) * t237;
t139 = -mrSges(5,1) * t182 - mrSges(5,3) * t237;
t138 = mrSges(5,2) * t182 - mrSges(5,3) * t238;
t137 = (Ifges(4,1) * t182 - t252) * qJD(3);
t134 = (-Ifges(4,2) * t179 + t251) * qJD(3);
t129 = (mrSges(4,1) * t179 + mrSges(4,2) * t182) * qJD(3);
t128 = t198 * qJD(4);
t127 = t197 * qJD(4);
t124 = t198 * t179;
t123 = t197 * t179;
t117 = qJD(4) * t191 - qJD(5) * t178;
t113 = t188 * t179;
t109 = -Ifges(6,2) * t182 + (Ifges(6,4) * t181 + Ifges(6,6) * t178) * t179;
t108 = -Ifges(5,3) * t182 + (Ifges(5,5) * t181 - Ifges(5,6) * t178) * t179;
t102 = -t178 * t254 + t241;
t95 = -mrSges(6,2) * t185 + mrSges(6,3) * t231;
t94 = -mrSges(5,2) * t231 - mrSges(5,3) * t185;
t92 = mrSges(5,1) * t231 - mrSges(5,3) * t186;
t88 = -t241 + (pkin(8) * t178 + pkin(4)) * t182;
t87 = -qJ(5) * t182 + t103;
t86 = -mrSges(4,1) * t239 - t119 * mrSges(4,3);
t85 = mrSges(4,2) * t239 - t118 * mrSges(4,3);
t76 = mrSges(5,1) * t185 + mrSges(5,2) * t186;
t75 = mrSges(6,1) * t185 - mrSges(6,3) * t186;
t71 = -Ifges(6,4) * t212 + t200;
t70 = -Ifges(5,5) * t212 - Ifges(5,6) * t185 + t234;
t68 = mrSges(4,1) * t214 - mrSges(4,3) * t82;
t67 = -mrSges(4,2) * t214 - mrSges(4,3) * t81;
t66 = Ifges(4,1) * t119 - Ifges(4,4) * t118 - Ifges(4,5) * t239;
t65 = Ifges(4,4) * t119 - Ifges(4,2) * t118 - t221;
t60 = t261 * t179 + t188 * t230;
t59 = -mrSges(6,1) * t118 + mrSges(6,2) * t84;
t58 = mrSges(5,1) * t118 - mrSges(5,3) * t84;
t57 = -mrSges(5,2) * t118 - mrSges(5,3) * t83;
t56 = -mrSges(6,2) * t83 + mrSges(6,3) * t118;
t55 = -pkin(4) * t231 - t62;
t50 = qJ(5) * t231 - qJD(5) * t182 + t61;
t47 = mrSges(5,1) * t83 + mrSges(5,2) * t84;
t46 = mrSges(6,1) * t83 - mrSges(6,3) * t84;
t45 = mrSges(4,1) * t81 + mrSges(4,2) * t82;
t44 = Ifges(4,1) * t82 - Ifges(4,4) * t81 + Ifges(4,5) * t214;
t43 = Ifges(4,4) * t82 - Ifges(4,2) * t81 + Ifges(4,6) * t214;
t31 = Ifges(6,4) * t84 + Ifges(6,2) * t118 + Ifges(6,6) * t83;
t30 = Ifges(5,5) * t84 - Ifges(5,6) * t83 + Ifges(5,3) * t118;
t22 = -mrSges(6,2) * t41 + mrSges(6,3) * t81;
t20 = mrSges(5,1) * t81 - mrSges(5,3) * t42;
t19 = -mrSges(5,2) * t81 - mrSges(5,3) * t41;
t18 = pkin(4) * t83 - qJ(5) * t84 + t53;
t15 = -pkin(4) * t118 - t16;
t14 = qJ(5) * t118 + t253;
t12 = mrSges(6,1) * t41 - mrSges(6,3) * t42;
t5 = pkin(4) * t41 - qJ(5) * t42 - qJD(5) * t84 + t25;
t2 = -pkin(4) * t81 - t4;
t1 = qJ(5) * t81 + qJD(5) * t118 + t3;
t17 = [(t16 * t4 + t25 * t53 + t253 * t3) * t260 + 0.2e1 * t253 * t19 + (t1 * t14 + t15 * t2 + t18 * t5) * t259 + (-t43 + 0.2e1 * t243 + t264) * t118 + (-t183 * t222 + 0.2e1 * (t115 * t183 + t116 * t180) * mrSges(3,3) + ((t121 * t257 + Ifges(3,5) * t177 + 0.2e1 * (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t183) * t176) * t183 + (t122 * t257 + Ifges(4,5) * t119 - 0.2e1 * Ifges(3,6) * t177 - Ifges(4,6) * t118 + (-0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t180 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3)) * t183) * t176) * t180) * qJD(2)) * t176 + 0.2e1 * m(4) * (t104 * t116 + t26 * t64 + t27 * t63) + 0.2e1 * m(3) * (t115 * t122 - t116 * t121) + (t10 + t11) * t84 + (t6 - t9) * t83 + (t30 + t31 - t65) * t81 + (t33 + t34) * t42 + (t29 - t32) * t41 + 0.2e1 * t18 * t12 + 0.2e1 * t16 * t20 + 0.2e1 * t15 * t21 + 0.2e1 * t14 * t22 + 0.2e1 * t5 * t46 + 0.2e1 * t25 * t47 + 0.2e1 * t53 * t13 + 0.2e1 * t1 * t56 + 0.2e1 * t3 * t57 + 0.2e1 * t4 * t58 + 0.2e1 * t2 * t59 + 0.2e1 * t64 * t67 + 0.2e1 * t63 * t68 + t82 * t66 + 0.2e1 * t26 * t85 + 0.2e1 * t27 * t86 + 0.2e1 * t104 * t45 + (t44 + 0.2e1 * t242) * t119 + (t158 - 0.2e1 * t244 - 0.2e1 * t245) * t177; (Ifges(4,5) * t199 - t27 * mrSges(4,3) + t44 / 0.2e1 + t242 + t220 * t181 + t224 * t178 + (-t178 * t218 + t181 * t219) * qJD(4) + (-t64 * mrSges(4,3) - t65 / 0.2e1 + t30 / 0.2e1 + t31 / 0.2e1 + t221 / 0.2e1) * qJD(3) + (-qJD(3) * t85 - t68 + m(4) * (-qJD(3) * t64 - t27) - t262) * pkin(8)) * t179 + (-m(4) * t116 - t45) * pkin(2) + t158 + m(5) * (t102 * t4 + t103 * t3 + t16 * t62 + t253 * t61) + t253 * t94 - t244 - t245 + m(6) * (t1 * t87 + t113 * t5 + t14 * t50 + t15 * t55 + t18 * t60 + t2 * t88) + (t108 / 0.2e1 + t109 / 0.2e1 - t153 / 0.2e1) * t81 + (t70 / 0.2e1 + t71 / 0.2e1 - t134 / 0.2e1) * t118 + t207 * t42 + t208 * t41 + t216 * t84 + t217 * t83 + t50 * t56 + t55 * t59 + t60 * t46 + t61 * t57 + t62 * t58 + t18 * t75 + t53 * t76 + t87 * t22 + t88 * t21 + t16 * t92 + t15 * t93 + t14 * t95 + t102 * t20 + t103 * t19 + t113 * t12 + (-t183 * t174 / 0.2e1 - Ifges(3,6) * t232) * t176 + t5 * t123 + t25 * t124 + t104 * t129 + t119 * t137 / 0.2e1 + t3 * t138 + t4 * t139 + t2 * t140 + t1 * t141 + t82 * t156 / 0.2e1 + (Ifges(4,6) * t199 + t26 * mrSges(4,3) - t7 / 0.2e1 - t8 / 0.2e1 + t43 / 0.2e1 - t243 + (m(4) * t26 + t67) * pkin(8) + (-t63 * mrSges(4,3) + t66 / 0.2e1 + t218 * t181 + t219 * t178 + (-m(4) * t63 + m(5) * t53 + t47 - t86) * pkin(8)) * qJD(3)) * t182; -0.2e1 * pkin(2) * t129 + 0.2e1 * t102 * t92 + 0.2e1 * t103 * t94 + 0.2e1 * t113 * t75 + 0.2e1 * t60 * t123 + 0.2e1 * t61 * t138 + 0.2e1 * t62 * t139 + 0.2e1 * t55 * t140 + 0.2e1 * t50 * t141 + 0.2e1 * t87 * t95 + 0.2e1 * t88 * t93 + (t113 * t60 + t50 * t87 + t55 * t88) * t259 + (t102 * t62 + t103 * t61) * t260 + (t134 - t70 - t71 + (t124 * t258 + t178 * t236 + t181 * t235 + t156) * qJD(3)) * t182 + (t76 * t258 + t137 + (t73 + t74) * t181 + (t69 - t72) * t178 + (-t178 * t235 + t181 * t236) * qJD(4) + (pkin(8) ^ 2 * t182 * t260 + t108 + t109 - t153) * qJD(3)) * t179; t222 + (t1 * mrSges(6,2) + t3 * mrSges(5,3) - t224) * t181 + m(6) * (t117 * t18 + t143 * t5) + t201 * t42 + t202 * t81 + t203 * t41 + t204 * t84 + t205 * t118 + t206 * t83 - t26 * mrSges(4,2) + t27 * mrSges(4,1) + ((t15 * mrSges(6,2) - t16 * mrSges(5,3) + t218) * t181 + (-t14 * mrSges(6,2) - mrSges(5,3) * t253 + t219) * t178) * qJD(4) + (t2 * mrSges(6,2) - t4 * mrSges(5,3) + t220) * t178 + ((t19 + t22) * t181 + (-t20 + t21) * t178 + ((-t58 + t59) * t181 + (-t56 - t57) * t178) * qJD(4) + m(5) * (-t16 * t227 - t178 * t4 + t181 * t3 - t229 * t253) + m(6) * (t1 * t181 - t14 * t229 + t15 * t227 + t178 * t2)) * pkin(9) + t117 * t46 + t18 * t127 + t53 * t128 + t143 * t12 + t5 * t147 + t25 * t148 + t262 * pkin(3); t174 + m(6) * (t113 * t117 + t143 * t60) + t179 * pkin(8) * t128 - pkin(3) * t76 + t117 * t123 + t113 * t127 + t143 * t75 + t60 * t147 - t205 * t182 + ((-Ifges(4,6) + t202) * t179 + (t179 * mrSges(4,2) + (-m(5) * pkin(3) - mrSges(4,1) + t148) * t182) * pkin(8)) * qJD(3) + (t61 * mrSges(5,3) + t50 * mrSges(6,2) + t204 * t179 + t201 * t230 + (t88 * mrSges(6,2) - t102 * mrSges(5,3) + t179 * t203 + t207) * qJD(4) + (t94 + t95 + (-t139 + t140) * qJD(4) + m(6) * (qJD(4) * t88 + t50) + m(5) * (-qJD(4) * t102 + t61)) * pkin(9) - t217) * t181 + (-t62 * mrSges(5,3) + t55 * mrSges(6,2) + t206 * t179 + t203 * t230 + (-t87 * mrSges(6,2) - t103 * mrSges(5,3) - t179 * t201 + t208) * qJD(4) + (-t92 + t93 + (-t138 - t141) * qJD(4) + m(6) * (-qJD(4) * t87 + t55) + m(5) * (-qJD(4) * t103 - t62)) * pkin(9) + t216) * t178; -0.2e1 * pkin(3) * t128 + 0.2e1 * t127 * t143 + (-t130 + t133) * t181 + (t135 + t136) * t178 + 0.2e1 * (m(6) * t143 + t147) * t117 + ((t154 + t155) * t181 + (t149 - t152) * t178) * qJD(4); m(6) * (-pkin(4) * t2 + qJ(5) * t1 + qJD(5) * t14) + t1 * mrSges(6,3) + qJD(5) * t56 + qJ(5) * t22 - t3 * mrSges(5,2) + t4 * mrSges(5,1) - t2 * mrSges(6,1) - pkin(4) * t21 + t264; -Ifges(5,6) * t210 - t61 * mrSges(5,2) - pkin(4) * t93 + m(6) * (-pkin(4) * t55 + qJ(5) * t50 + qJD(5) * t87) + qJD(5) * t141 + qJ(5) * t95 + t50 * mrSges(6,3) - t55 * mrSges(6,1) + t62 * mrSges(5,1) + (-t246 - t266) * t228 + t200 + t234; -t261 * mrSges(6,2) + (m(6) * t225 + (-m(6) * t192 + t147 + t148) * qJD(4)) * pkin(9) + t131 + t132; 0.2e1 * (m(6) * qJ(5) + mrSges(6,3)) * qJD(5); m(6) * t2 + t21; m(6) * t55 + t93; (m(6) * pkin(9) + mrSges(6,2)) * t227; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t17(1), t17(2), t17(4), t17(7), t17(11); t17(2), t17(3), t17(5), t17(8), t17(12); t17(4), t17(5), t17(6), t17(9), t17(13); t17(7), t17(8), t17(9), t17(10), t17(14); t17(11), t17(12), t17(13), t17(14), t17(15);];
Mq = res;
