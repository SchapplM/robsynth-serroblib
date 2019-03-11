% Calculate time derivative of joint inertia matrix for
% S6RRRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:05:18
% EndTime: 2019-03-09 21:05:29
% DurationCPUTime: 4.79s
% Computational Cost: add. (4548->517), mult. (11130->700), div. (0->0), fcn. (9126->6), ass. (0->212)
t265 = Ifges(6,4) + Ifges(5,5);
t271 = Ifges(7,5) - t265;
t264 = -Ifges(5,6) + Ifges(6,6);
t270 = Ifges(7,6) - t264;
t269 = mrSges(6,3) + mrSges(7,2);
t183 = sin(qJ(2));
t182 = sin(qJ(3));
t186 = cos(qJ(2));
t223 = qJD(2) * t186;
t211 = t182 * t223;
t185 = cos(qJ(3));
t220 = qJD(3) * t185;
t193 = t183 * t220 + t211;
t268 = -Ifges(6,2) - Ifges(5,3) - Ifges(7,3);
t181 = sin(qJ(4));
t184 = cos(qJ(4));
t244 = -mrSges(6,1) - mrSges(7,1);
t267 = ((-mrSges(5,1) + t244) * t181 - mrSges(5,2) * t184) * pkin(3) * qJD(4);
t266 = 2 * mrSges(6,2) + 2 * mrSges(5,3);
t137 = t181 * t185 + t182 * t184;
t124 = t137 * t183;
t252 = -pkin(9) - pkin(8);
t155 = t252 * t182;
t156 = t252 * t185;
t110 = t181 * t155 - t184 * t156;
t214 = t185 * t223;
t224 = qJD(2) * t183;
t263 = -Ifges(4,5) * t214 - Ifges(4,3) * t224;
t149 = -pkin(2) * t186 - pkin(8) * t183 - pkin(1);
t230 = t185 * t186;
t164 = pkin(7) * t230;
t121 = t182 * t149 + t164;
t221 = qJD(3) * t183;
t213 = t182 * t221;
t194 = -t213 + t214;
t262 = qJD(3) + qJD(4);
t261 = 2 * m(4);
t260 = 2 * m(5);
t259 = 2 * m(6);
t258 = 2 * m(7);
t257 = -0.2e1 * pkin(1);
t256 = 0.2e1 * pkin(7);
t255 = -2 * mrSges(7,3);
t254 = 2 * mrSges(7,3);
t253 = m(6) + m(7);
t187 = -pkin(4) - pkin(5);
t251 = -t182 / 0.2e1;
t250 = m(7) * t181;
t136 = t181 * t182 - t184 * t185;
t97 = t262 * t136;
t249 = Ifges(7,5) * t97;
t98 = t262 * t137;
t248 = Ifges(7,6) * t98;
t247 = pkin(7) * t182;
t60 = -t124 * t262 - t136 * t223;
t246 = t60 * mrSges(7,3);
t242 = mrSges(6,2) - mrSges(7,3);
t219 = qJD(4) * t181;
t232 = t183 * t185;
t233 = t182 * t183;
t61 = -t219 * t233 + (t232 * t262 + t211) * t184 + t194 * t181;
t34 = -mrSges(6,2) * t61 + mrSges(6,3) * t224;
t38 = mrSges(7,2) * t224 + mrSges(7,3) * t61;
t241 = t34 + t38;
t108 = -pkin(9) * t233 + t121;
t135 = t185 * t149;
t89 = -pkin(9) * t232 + t135 + (-pkin(3) - t247) * t186;
t41 = t184 * t108 + t181 * t89;
t239 = Ifges(4,4) * t182;
t238 = Ifges(4,4) * t185;
t237 = Ifges(4,6) * t182;
t235 = t186 * Ifges(4,6);
t109 = -t184 * t155 - t156 * t181;
t234 = t109 * t181;
t154 = Ifges(4,1) * t182 + t238;
t231 = t185 * t154;
t111 = -mrSges(6,2) * t124 - mrSges(6,3) * t186;
t112 = -mrSges(7,2) * t186 + mrSges(7,3) * t124;
t229 = t111 + t112;
t125 = t136 * t183;
t115 = -mrSges(5,1) * t186 + mrSges(5,3) * t125;
t116 = mrSges(6,1) * t186 - mrSges(6,2) * t125;
t228 = -t115 + t116;
t147 = (pkin(2) * t183 - pkin(8) * t186) * qJD(2);
t227 = t185 * t147 + t224 * t247;
t218 = qJD(4) * t184;
t162 = pkin(3) * t218 + qJD(5);
t166 = pkin(3) * t181 + qJ(5);
t226 = qJ(5) * t162 + qJD(5) * t166;
t148 = pkin(3) * t233 + t183 * pkin(7);
t225 = t269 * qJD(5);
t222 = qJD(3) * t182;
t217 = pkin(3) * t222;
t216 = pkin(3) * t219;
t176 = pkin(7) * t223;
t119 = pkin(3) * t193 + t176;
t170 = -pkin(3) * t185 - pkin(2);
t169 = -pkin(3) * t184 - pkin(4);
t215 = qJD(3) * t252;
t19 = -t61 * mrSges(7,1) + t60 * mrSges(7,2);
t45 = -t98 * mrSges(7,1) - t97 * mrSges(7,2);
t210 = t269 * t162;
t209 = (2 * Ifges(3,4)) + t237;
t146 = t182 * t215;
t204 = t185 * t215;
t64 = t184 * t146 + t155 * t218 + t156 * t219 + t181 * t204;
t65 = qJD(4) * t110 + t146 * t181 - t184 * t204;
t208 = t109 * t65 + t110 * t64;
t40 = -t181 * t108 + t184 * t89;
t30 = -qJ(5) * t186 + t41;
t31 = t186 * pkin(4) - t40;
t203 = -qJ(5) * t125 - t148;
t202 = -mrSges(4,1) * t185 + mrSges(4,2) * t182;
t201 = mrSges(4,1) * t182 + mrSges(4,2) * t185;
t200 = Ifges(4,1) * t185 - t239;
t199 = -Ifges(4,2) * t182 + t238;
t153 = Ifges(4,2) * t185 + t239;
t198 = Ifges(4,5) * t182 + Ifges(4,6) * t185;
t197 = t136 * t162 + t166 * t98;
t33 = (pkin(3) * t183 - pkin(9) * t230) * qJD(2) + (-t164 + (pkin(9) * t183 - t149) * t182) * qJD(3) + t227;
t68 = t182 * t147 + t149 * t220 + (-t185 * t224 - t186 * t222) * pkin(7);
t53 = -pkin(9) * t193 + t68;
t10 = -t108 * t218 - t181 * t53 + t184 * t33 - t89 * t219;
t196 = -qJ(5) * t98 - qJD(5) * t136;
t195 = qJ(5) * t137 - t170;
t9 = -t108 * t219 + t181 * t33 + t184 * t53 + t89 * t218;
t192 = qJ(5) * t60 - qJD(5) * t125 - t119;
t191 = -qJ(5) * t97 + qJD(5) * t137 - t217;
t6 = qJ(5) * t224 - qJD(5) * t186 + t9;
t190 = t268 * t224 + t270 * t61 + t271 * t60;
t23 = qJ(6) * t98 + qJD(6) * t136 + t64;
t24 = qJ(6) * t97 - qJD(6) * t137 + t65;
t91 = Ifges(6,6) * t98;
t92 = Ifges(5,6) * t98;
t93 = Ifges(5,5) * t97;
t94 = Ifges(6,4) * t97;
t189 = -t24 * mrSges(7,1) + t23 * mrSges(7,2) - t248 + t249 + t91 - t92 - t93 - t94 + (-mrSges(5,1) - mrSges(6,1)) * t65 + (-mrSges(5,2) + mrSges(6,3)) * t64;
t2 = -qJ(6) * t60 + qJD(6) * t125 + t187 * t224 - t10;
t3 = qJ(6) * t61 + qJD(6) * t124 + t6;
t7 = -pkin(4) * t224 - t10;
t188 = t10 * mrSges(5,1) - t7 * mrSges(6,1) - t2 * mrSges(7,1) - t9 * mrSges(5,2) + t3 * mrSges(7,2) + t6 * mrSges(6,3) - t190;
t175 = Ifges(4,5) * t220;
t165 = -pkin(5) + t169;
t145 = -mrSges(4,1) * t186 - mrSges(4,3) * t232;
t144 = mrSges(4,2) * t186 - mrSges(4,3) * t233;
t143 = t200 * qJD(3);
t142 = t199 * qJD(3);
t141 = t201 * qJD(3);
t132 = t166 * t162;
t123 = -Ifges(4,5) * t186 + t183 * t200;
t122 = t183 * t199 - t235;
t120 = -t186 * t247 + t135;
t118 = -mrSges(4,2) * t224 - mrSges(4,3) * t193;
t117 = mrSges(4,1) * t224 - mrSges(4,3) * t194;
t114 = mrSges(7,1) * t186 + mrSges(7,3) * t125;
t113 = mrSges(5,2) * t186 - mrSges(5,3) * t124;
t107 = Ifges(5,1) * t137 - Ifges(5,4) * t136;
t106 = Ifges(6,1) * t137 + Ifges(6,5) * t136;
t105 = Ifges(7,1) * t137 + Ifges(7,4) * t136;
t104 = Ifges(5,4) * t137 - Ifges(5,2) * t136;
t103 = Ifges(7,4) * t137 + Ifges(7,2) * t136;
t102 = Ifges(6,5) * t137 + Ifges(6,3) * t136;
t101 = mrSges(5,1) * t136 + mrSges(5,2) * t137;
t100 = -mrSges(7,1) * t136 + mrSges(7,2) * t137;
t99 = mrSges(6,1) * t136 - mrSges(6,3) * t137;
t87 = mrSges(4,1) * t193 + mrSges(4,2) * t194;
t86 = pkin(4) * t136 - t195;
t82 = mrSges(5,1) * t124 - mrSges(5,2) * t125;
t81 = -mrSges(7,1) * t124 - mrSges(7,2) * t125;
t80 = mrSges(6,1) * t124 + mrSges(6,3) * t125;
t79 = qJ(6) * t136 + t110;
t78 = -qJ(6) * t137 + t109;
t77 = -t154 * t221 + (Ifges(4,5) * t183 + t186 * t200) * qJD(2);
t76 = -t153 * t221 + (Ifges(4,6) * t183 + t186 * t199) * qJD(2);
t75 = -Ifges(5,1) * t125 - Ifges(5,4) * t124 - Ifges(5,5) * t186;
t74 = -Ifges(6,1) * t125 - Ifges(6,4) * t186 + Ifges(6,5) * t124;
t73 = -Ifges(7,1) * t125 + Ifges(7,4) * t124 + t186 * Ifges(7,5);
t72 = -Ifges(5,4) * t125 - Ifges(5,2) * t124 - Ifges(5,6) * t186;
t71 = -Ifges(7,4) * t125 + Ifges(7,2) * t124 + t186 * Ifges(7,6);
t70 = -Ifges(6,5) * t125 - Ifges(6,6) * t186 + Ifges(6,3) * t124;
t69 = -qJD(3) * t121 + t227;
t67 = t136 * t187 + t195;
t66 = pkin(4) * t124 - t203;
t55 = t60 * mrSges(6,2);
t52 = -Ifges(5,1) * t97 - Ifges(5,4) * t98;
t51 = -Ifges(6,1) * t97 + Ifges(6,5) * t98;
t50 = -Ifges(7,1) * t97 + Ifges(7,4) * t98;
t49 = -Ifges(5,4) * t97 - Ifges(5,2) * t98;
t48 = -Ifges(7,4) * t97 + Ifges(7,2) * t98;
t47 = -Ifges(6,5) * t97 + Ifges(6,3) * t98;
t46 = mrSges(5,1) * t98 - mrSges(5,2) * t97;
t44 = mrSges(6,1) * t98 + mrSges(6,3) * t97;
t39 = -mrSges(5,2) * t224 - mrSges(5,3) * t61;
t37 = -mrSges(6,1) * t224 + t55;
t36 = mrSges(5,1) * t224 - mrSges(5,3) * t60;
t35 = -mrSges(7,1) * t224 - t246;
t32 = t124 * t187 + t203;
t28 = pkin(4) * t98 - t191;
t26 = qJ(6) * t124 + t30;
t25 = pkin(5) * t186 + qJ(6) * t125 + t31;
t21 = t187 * t98 + t191;
t20 = mrSges(5,1) * t61 + mrSges(5,2) * t60;
t18 = mrSges(6,1) * t61 - mrSges(6,3) * t60;
t17 = Ifges(5,1) * t60 - Ifges(5,4) * t61 + Ifges(5,5) * t224;
t16 = Ifges(6,1) * t60 + Ifges(6,4) * t224 + Ifges(6,5) * t61;
t15 = Ifges(7,1) * t60 + Ifges(7,4) * t61 - Ifges(7,5) * t224;
t14 = Ifges(5,4) * t60 - Ifges(5,2) * t61 + Ifges(5,6) * t224;
t13 = Ifges(7,4) * t60 + Ifges(7,2) * t61 - Ifges(7,6) * t224;
t12 = Ifges(6,5) * t60 + Ifges(6,6) * t224 + Ifges(6,3) * t61;
t11 = pkin(4) * t61 - t192;
t4 = t187 * t61 + t192;
t1 = [(t87 * t256 - t182 * t76 + t185 * t77 + (-t185 * t122 - t182 * t123 + t186 * t198) * qJD(3) + (mrSges(3,1) * t257 + (Ifges(4,5) * t185 - t209) * t183 + t271 * t125 - t270 * t124 + (pkin(7) ^ 2 * t261 + t201 * t256 + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3) + t268) * t186) * qJD(2)) * t183 + 0.2e1 * t68 * t144 + 0.2e1 * t69 * t145 + 0.2e1 * t148 * t20 + ((mrSges(3,2) * t257 - t182 * t122 + t185 * t123 + t186 * t209) * qJD(2) + t190 + t263) * t186 + 0.2e1 * t6 * t111 + 0.2e1 * t3 * t112 + 0.2e1 * t9 * t113 + 0.2e1 * t2 * t114 + 0.2e1 * t10 * t115 + 0.2e1 * t7 * t116 + 0.2e1 * t119 * t82 + 0.2e1 * t120 * t117 + 0.2e1 * t121 * t118 + 0.2e1 * t11 * t80 + 0.2e1 * t4 * t81 + 0.2e1 * t66 * t18 + 0.2e1 * t31 * t37 + 0.2e1 * t26 * t38 + 0.2e1 * t40 * t36 + 0.2e1 * t41 * t39 + 0.2e1 * t32 * t19 + 0.2e1 * t30 * t34 + 0.2e1 * t25 * t35 + (t70 + t71 - t72) * t61 + (t73 + t74 + t75) * t60 + (t12 + t13 - t14) * t124 - (t15 + t16 + t17) * t125 + (t2 * t25 + t26 * t3 + t32 * t4) * t258 + (t11 * t66 + t30 * t6 + t31 * t7) * t259 + (t10 * t40 + t119 * t148 + t41 * t9) * t260 + (t120 * t69 + t121 * t68) * t261; m(7) * (t2 * t78 + t21 * t32 + t23 * t26 + t24 * t25 + t3 * t79 + t4 * t67) + m(6) * (t109 * t7 + t11 * t86 + t110 * t6 + t28 * t66 + t30 * t64 + t31 * t65) - (t50 / 0.2e1 + t51 / 0.2e1 + t52 / 0.2e1) * t125 + (t47 / 0.2e1 + t48 / 0.2e1 - t49 / 0.2e1) * t124 + (-pkin(2) * t176 + (-t69 * t182 + t68 * t185 + (-t120 * t185 - t121 * t182) * qJD(3)) * pkin(8)) * m(4) + (t70 / 0.2e1 + t71 / 0.2e1 - t72 / 0.2e1 + t26 * mrSges(7,3) - t30 * mrSges(6,2) - t41 * mrSges(5,3)) * t98 - (t74 / 0.2e1 + t75 / 0.2e1 + t73 / 0.2e1 - t25 * mrSges(7,3) + t31 * mrSges(6,2) - t40 * mrSges(5,3)) * t97 + (t76 / 0.2e1 + pkin(8) * t118 + t68 * mrSges(4,3)) * t185 + (t77 / 0.2e1 - pkin(8) * t117 - t69 * mrSges(4,3)) * t182 + (t185 * t143 / 0.2e1 + t142 * t251 - qJD(2) * (Ifges(7,5) * t137 + Ifges(7,6) * t136) / 0.2e1 - Ifges(3,6) * qJD(2) + (-t185 * t153 / 0.2e1 + t154 * t251) * qJD(3) + (mrSges(3,2) * qJD(2) + t141) * pkin(7) + (t264 * t136 + t265 * t137 + t198) * qJD(2) / 0.2e1) * t183 + t170 * t20 + t148 * t46 + t23 * t112 + t24 * t114 + t119 * t101 + t11 * t99 + t4 * t100 + t78 * t35 + t79 * t38 + t28 * t80 + t21 * t81 + t86 * t18 - pkin(2) * t87 + t66 * t44 + t67 * t19 + t32 * t45 + (t106 / 0.2e1 + t107 / 0.2e1 + t105 / 0.2e1) * t60 + (t111 + t113) * t64 + (t39 + t34) * t110 + (t37 - t36) * t109 + m(5) * (-t10 * t109 + t110 * t9 + t119 * t170 + t148 * t217 - t40 * t65 + t41 * t64) + (-t175 / 0.2e1 - t249 / 0.2e1 + t248 / 0.2e1 + t93 / 0.2e1 + t92 / 0.2e1 + t94 / 0.2e1 - t91 / 0.2e1 + (t231 / 0.2e1 + t153 * t251 + Ifges(3,5) + (-mrSges(3,1) + t202) * pkin(7)) * qJD(2)) * t186 + t228 * t65 + ((-pkin(8) * t145 - t120 * mrSges(4,3) + t123 / 0.2e1) * t185 + (-pkin(8) * t144 - t121 * mrSges(4,3) + pkin(3) * t82 + t235 / 0.2e1 - t122 / 0.2e1) * t182) * qJD(3) + (t102 / 0.2e1 + t103 / 0.2e1 - t104 / 0.2e1) * t61 + (t12 / 0.2e1 + t13 / 0.2e1 - t14 / 0.2e1 + t3 * mrSges(7,3) - t6 * mrSges(6,2) - t9 * mrSges(5,3)) * t136 + (t15 / 0.2e1 + t16 / 0.2e1 + t17 / 0.2e1 - t2 * mrSges(7,3) + t7 * mrSges(6,2) - t10 * mrSges(5,3)) * t137; -0.2e1 * pkin(2) * t141 + 0.2e1 * t21 * t100 + t185 * t142 + t182 * t143 + 0.2e1 * t170 * t46 + 0.2e1 * t28 * t99 + 0.2e1 * t86 * t44 + 0.2e1 * t67 * t45 + (t231 + (0.2e1 * pkin(3) * t101 - t153) * t182) * qJD(3) + (t170 * t217 + t208) * t260 + (t21 * t67 + t23 * t79 + t24 * t78) * t258 + (t28 * t86 + t208) * t259 + (-t110 * t266 + t254 * t79 + t102 + t103 - t104) * t98 - (t109 * t266 + t255 * t78 + t105 + t106 + t107) * t97 + (t24 * t255 + t266 * t65 + t50 + t51 + t52) * t137 + (t23 * t254 - t266 * t64 + t47 + t48 - t49) * t136; t188 + (t25 * t250 * qJD(4) + (m(5) * t10 + t36 + (m(5) * t41 + t113) * qJD(4)) * t184 + (m(5) * t9 + t39 + (-m(5) * t40 + m(6) * t31 + t114 + t228) * qJD(4)) * t181) * pkin(3) - Ifges(4,5) * t213 + t165 * t35 + t169 * t37 - t68 * mrSges(4,2) + t69 * mrSges(4,1) + t241 * t166 + t229 * t162 + m(7) * (t162 * t26 + t165 * t2 + t166 * t3) + m(6) * (t162 * t30 + t166 * t6 + t169 * t7) - t193 * Ifges(4,6) - t263; t189 + (m(5) * (t181 * t64 - t184 * t65) + (-t181 * t98 + t184 * t97) * mrSges(5,3) + (-t184 * t136 * mrSges(5,3) + (mrSges(5,3) + t242) * t181 * t137 + m(5) * (t110 * t184 + t234) + t78 * t250 + m(6) * t234) * qJD(4)) * pkin(3) + t175 + (pkin(8) * t202 - t237) * qJD(3) + (-t169 * t97 - t197) * mrSges(6,2) + (t165 * t97 + t197) * mrSges(7,3) + m(7) * (t162 * t79 + t165 * t24 + t166 * t23) + m(6) * (t110 * t162 + t166 * t64 + t169 * t65); (t169 * t216 + t132) * t259 + (t165 * t216 + t132) * t258 + 0.2e1 * t210 + 0.2e1 * t267; t188 + t229 * qJD(5) + t241 * qJ(5) + m(6) * (-pkin(4) * t7 + qJ(5) * t6 + qJD(5) * t30) + m(7) * (qJ(5) * t3 + qJD(5) * t26 + t187 * t2) + t187 * t35 - pkin(4) * t37; t189 + m(6) * (-pkin(4) * t65 + qJ(5) * t64 + qJD(5) * t110) + m(7) * (qJ(5) * t23 + qJD(5) * t79 + t187 * t24) + (t187 * t97 - t196) * mrSges(7,3) + (pkin(4) * t97 + t196) * mrSges(6,2); t210 + t267 + m(6) * (-pkin(4) * t216 + t226) + m(7) * (t187 * t216 + t226) + t225; 0.2e1 * qJ(5) * qJD(5) * t253 + 0.2e1 * t225; m(6) * t7 + m(7) * t2 + t224 * t244 - t246 + t55; m(6) * t65 + m(7) * t24 - t242 * t97; t253 * t216; 0; 0; m(7) * t4 + t19; m(7) * t21 + t45; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
