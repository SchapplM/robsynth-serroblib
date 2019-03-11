% Calculate time derivative of joint inertia matrix for
% S6RRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:15:58
% EndTime: 2019-03-09 10:16:09
% DurationCPUTime: 4.38s
% Computational Cost: add. (8772->457), mult. (19483->683), div. (0->0), fcn. (19334->10), ass. (0->188)
t239 = Ifges(5,3) + Ifges(6,3);
t172 = sin(pkin(11));
t174 = cos(pkin(11));
t177 = sin(qJ(4));
t180 = cos(qJ(4));
t155 = t172 * t180 + t174 * t177;
t147 = t155 * qJD(4);
t188 = t172 * t177 - t174 * t180;
t149 = t188 * qJD(4);
t206 = qJD(4) * t180;
t238 = Ifges(5,5) * t206 - Ifges(6,5) * t149 - Ifges(6,6) * t147;
t173 = sin(pkin(10));
t175 = cos(pkin(10));
t178 = sin(qJ(2));
t181 = cos(qJ(2));
t154 = t173 * t178 - t175 * t181;
t150 = t154 * qJD(2);
t156 = t173 * t181 + t175 * t178;
t199 = t156 * t206;
t186 = -t177 * t150 + t199;
t176 = sin(qJ(6));
t179 = cos(qJ(6));
t170 = -pkin(2) * t181 - pkin(1);
t109 = t154 * pkin(3) - t156 * pkin(8) + t170;
t223 = -qJ(3) - pkin(7);
t162 = t223 * t178;
t163 = t223 * t181;
t119 = t162 * t173 - t163 * t175;
t117 = t180 * t119;
t148 = t156 * qJD(2);
t187 = qJ(5) * t150 - qJD(5) * t156;
t196 = qJD(2) * t223;
t146 = qJD(3) * t181 + t178 * t196;
t182 = -t178 * qJD(3) + t181 * t196;
t97 = t175 * t146 + t173 * t182;
t202 = pkin(2) * qJD(2) * t178;
t98 = pkin(3) * t148 + pkin(8) * t150 + t202;
t197 = -t177 * t97 + t180 * t98;
t27 = pkin(4) * t148 + t187 * t180 + (-t117 + (qJ(5) * t156 - t109) * t177) * qJD(4) + t197;
t205 = t109 * t206 + t177 * t98 + t180 * t97;
t29 = -qJ(5) * t199 + (-qJD(4) * t119 + t187) * t177 + t205;
t11 = -t172 * t29 + t174 * t27;
t63 = -t147 * t156 + t150 * t188;
t4 = pkin(5) * t148 - pkin(9) * t63 + t11;
t12 = t172 * t27 + t174 * t29;
t62 = t149 * t156 + t150 * t155;
t7 = pkin(9) * t62 + t12;
t101 = t188 * t156;
t213 = t156 * t180;
t80 = t180 * t109 - t119 * t177;
t52 = pkin(4) * t154 - qJ(5) * t213 + t80;
t214 = t156 * t177;
t81 = t177 * t109 + t117;
t61 = -qJ(5) * t214 + t81;
t30 = -t172 * t61 + t174 * t52;
t23 = pkin(5) * t154 + pkin(9) * t101 + t30;
t100 = t155 * t156;
t31 = t172 * t52 + t174 * t61;
t26 = -pkin(9) * t100 + t31;
t9 = -t176 * t26 + t179 * t23;
t2 = qJD(6) * t9 + t176 * t4 + t179 * t7;
t10 = t176 * t23 + t179 * t26;
t3 = -qJD(6) * t10 - t176 * t7 + t179 * t4;
t237 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t193 = mrSges(5,1) * t177 + mrSges(5,2) * t180;
t158 = t193 * qJD(4);
t236 = 2 * m(6);
t235 = 2 * m(7);
t234 = -2 * mrSges(4,3);
t118 = -t175 * t162 - t163 * t173;
t232 = 0.2e1 * t118;
t231 = 0.2e1 * t170;
t230 = m(6) * pkin(4);
t226 = -t156 / 0.2e1;
t221 = Ifges(5,4) * t177;
t164 = Ifges(5,2) * t180 + t221;
t225 = -t164 / 0.2e1;
t224 = pkin(4) * t172;
t110 = -t155 * t176 - t179 * t188;
t73 = qJD(6) * t110 - t147 * t176 - t149 * t179;
t111 = t155 * t179 - t176 * t188;
t74 = -qJD(6) * t111 - t147 * t179 + t149 * t176;
t222 = Ifges(7,5) * t73 + Ifges(7,6) * t74;
t220 = Ifges(5,4) * t180;
t219 = Ifges(5,6) * t177;
t96 = t146 * t173 - t175 * t182;
t218 = t118 * t96;
t217 = t148 * Ifges(5,5);
t216 = t148 * Ifges(5,6);
t215 = t154 * Ifges(5,6);
t210 = t180 * t150;
t167 = pkin(2) * t173 + pkin(8);
t209 = qJ(5) + t167;
t194 = qJD(4) * t209;
t125 = qJD(5) * t180 - t177 * t194;
t126 = -qJD(5) * t177 - t180 * t194;
t87 = t174 * t125 + t172 * t126;
t151 = t209 * t177;
t152 = t209 * t180;
t105 = -t172 * t151 + t174 * t152;
t207 = qJD(4) * t177;
t64 = -t100 * t179 + t101 * t176;
t18 = qJD(6) * t64 + t176 * t62 + t179 * t63;
t65 = -t100 * t176 - t101 * t179;
t19 = -qJD(6) * t65 - t176 * t63 + t179 * t62;
t203 = Ifges(7,5) * t18 + Ifges(7,6) * t19 + Ifges(7,3) * t148;
t201 = pkin(4) * t207;
t169 = -pkin(2) * t175 - pkin(3);
t200 = t156 * t207;
t36 = -t62 * mrSges(6,1) + t63 * mrSges(6,2);
t8 = -t19 * mrSges(7,1) + t18 * mrSges(7,2);
t41 = -t74 * mrSges(7,1) + t73 * mrSges(7,2);
t198 = -(2 * Ifges(4,4)) - t219;
t195 = 0.2e1 * t202;
t106 = t147 * mrSges(6,1) - t149 * mrSges(6,2);
t86 = -t125 * t172 + t174 * t126;
t104 = -t174 * t151 - t152 * t172;
t95 = pkin(4) * t214 + t118;
t168 = pkin(4) * t174 + pkin(5);
t144 = t168 * t179 - t176 * t224;
t130 = t144 * qJD(6);
t145 = t168 * t176 + t179 * t224;
t131 = t145 * qJD(6);
t192 = -t131 * mrSges(7,1) - t130 * mrSges(7,2);
t191 = Ifges(5,1) * t180 - t221;
t190 = -Ifges(5,2) * t177 + t220;
t88 = -pkin(9) * t155 + t104;
t89 = -pkin(9) * t188 + t105;
t44 = -t176 * t89 + t179 * t88;
t45 = t176 * t88 + t179 * t89;
t78 = pkin(9) * t149 + t86;
t79 = -pkin(9) * t147 + t87;
t21 = qJD(6) * t44 + t176 * t78 + t179 * t79;
t22 = -qJD(6) * t45 - t176 * t79 + t179 * t78;
t189 = t22 * mrSges(7,1) - t21 * mrSges(7,2) + t222;
t161 = -pkin(4) * t180 + t169;
t68 = pkin(4) * t186 + t96;
t185 = t200 + t210;
t184 = -t106 - t41;
t183 = -Ifges(5,5) * t210 + Ifges(6,5) * t63 + Ifges(6,6) * t62 + t148 * t239 + t203;
t165 = Ifges(5,1) * t177 + t220;
t160 = t191 * qJD(4);
t159 = t190 * qJD(4);
t137 = t150 * mrSges(4,2);
t124 = pkin(5) * t147 + t201;
t123 = pkin(5) * t188 + t161;
t116 = Ifges(6,1) * t155 - Ifges(6,4) * t188;
t115 = Ifges(6,4) * t155 - Ifges(6,2) * t188;
t114 = mrSges(6,1) * t188 + mrSges(6,2) * t155;
t113 = mrSges(5,1) * t154 - mrSges(5,3) * t213;
t112 = -mrSges(5,2) * t154 - mrSges(5,3) * t214;
t108 = -Ifges(6,1) * t149 - Ifges(6,4) * t147;
t107 = -Ifges(6,4) * t149 - Ifges(6,2) * t147;
t91 = Ifges(5,5) * t154 + t156 * t191;
t90 = t156 * t190 + t215;
t85 = mrSges(6,1) * t154 + mrSges(6,3) * t101;
t84 = -mrSges(6,2) * t154 - mrSges(6,3) * t100;
t83 = -mrSges(5,2) * t148 - mrSges(5,3) * t186;
t82 = mrSges(5,1) * t148 + mrSges(5,3) * t185;
t77 = Ifges(7,1) * t111 + Ifges(7,4) * t110;
t76 = Ifges(7,4) * t111 + Ifges(7,2) * t110;
t75 = -mrSges(7,1) * t110 + mrSges(7,2) * t111;
t69 = mrSges(5,1) * t186 - mrSges(5,2) * t185;
t67 = mrSges(6,1) * t100 - mrSges(6,2) * t101;
t66 = pkin(5) * t100 + t95;
t56 = -Ifges(5,1) * t185 - Ifges(5,4) * t186 + t217;
t55 = -Ifges(5,4) * t185 - Ifges(5,2) * t186 + t216;
t54 = -Ifges(6,1) * t101 - Ifges(6,4) * t100 + Ifges(6,5) * t154;
t53 = -Ifges(6,4) * t101 - Ifges(6,2) * t100 + Ifges(6,6) * t154;
t51 = mrSges(7,1) * t154 - mrSges(7,3) * t65;
t50 = -mrSges(7,2) * t154 + mrSges(7,3) * t64;
t47 = mrSges(6,1) * t148 - mrSges(6,3) * t63;
t46 = -mrSges(6,2) * t148 + mrSges(6,3) * t62;
t43 = Ifges(7,1) * t73 + Ifges(7,4) * t74;
t42 = Ifges(7,4) * t73 + Ifges(7,2) * t74;
t40 = -qJD(4) * t81 + t197;
t39 = -t119 * t207 + t205;
t38 = -pkin(5) * t62 + t68;
t37 = -mrSges(7,1) * t64 + mrSges(7,2) * t65;
t35 = Ifges(7,1) * t65 + Ifges(7,4) * t64 + Ifges(7,5) * t154;
t34 = Ifges(7,4) * t65 + Ifges(7,2) * t64 + Ifges(7,6) * t154;
t33 = Ifges(6,1) * t63 + Ifges(6,4) * t62 + t148 * Ifges(6,5);
t32 = Ifges(6,4) * t63 + Ifges(6,2) * t62 + t148 * Ifges(6,6);
t14 = -mrSges(7,2) * t148 + mrSges(7,3) * t19;
t13 = mrSges(7,1) * t148 - mrSges(7,3) * t18;
t6 = Ifges(7,1) * t18 + Ifges(7,4) * t19 + t148 * Ifges(7,5);
t5 = Ifges(7,4) * t18 + Ifges(7,2) * t19 + t148 * Ifges(7,6);
t1 = [(t180 * t56 - t177 * t55 - 0.2e1 * Ifges(4,1) * t150 + mrSges(4,2) * t195 + (Ifges(5,5) * t180 + t198) * t148 + (t154 * (-Ifges(5,5) * t177 - Ifges(5,6) * t180) - t177 * t91 - t180 * t90) * qJD(4) + 0.2e1 * (mrSges(4,3) + t193) * t96) * t156 + (mrSges(4,1) * t195 + t97 * t234 - t198 * t150 + ((2 * Ifges(4,2)) + Ifges(7,3) + t239) * t148 + t183) * t154 + (mrSges(4,1) * t231 - Ifges(6,5) * t101 + Ifges(7,5) * t65 - Ifges(6,6) * t100 + Ifges(7,6) * t64 + t119 * t234) * t148 - (mrSges(4,3) * t232 - t177 * t90 + t180 * t91) * t150 + 0.2e1 * m(5) * (t39 * t81 + t40 * t80 + t218) + 0.2e1 * m(4) * (t119 * t97 + t170 * t202 + t218) + 0.2e1 * (-pkin(1) * (mrSges(3,1) * t178 + mrSges(3,2) * t181) + (-Ifges(3,2) + Ifges(3,1)) * t178 * t181 + (-t178 ^ 2 + t181 ^ 2) * Ifges(3,4)) * qJD(2) + 0.2e1 * t39 * t112 + 0.2e1 * t40 * t113 - t100 * t32 - t101 * t33 + 0.2e1 * t95 * t36 + 0.2e1 * t80 * t82 + 0.2e1 * t81 * t83 + 0.2e1 * t12 * t84 + 0.2e1 * t11 * t85 + 0.2e1 * t68 * t67 + t62 * t53 + t63 * t54 + t64 * t5 + t65 * t6 + 0.2e1 * t66 * t8 + 0.2e1 * t2 * t50 + 0.2e1 * t3 * t51 + 0.2e1 * t31 * t46 + 0.2e1 * t30 * t47 + 0.2e1 * t38 * t37 + t19 * t34 + t18 * t35 - t137 * t231 + t69 * t232 + (t10 * t2 + t3 * t9 + t38 * t66) * t235 + (t11 * t30 + t12 * t31 + t68 * t95) * t236 + 0.2e1 * t9 * t13 + 0.2e1 * t10 * t14; (m(5) * t96 + t69) * t169 + m(6) * (t104 * t11 + t105 * t12 + t161 * t68 + t30 * t86 + t31 * t87) + (t10 * t74 + t110 * t2 - t111 * t3 - t73 * t9) * mrSges(7,3) + (Ifges(3,5) * t181 - Ifges(3,6) * t178 + (-mrSges(3,1) * t181 + mrSges(3,2) * t178) * pkin(7)) * qJD(2) + (t159 * t226 - t150 * t225 - t40 * mrSges(5,3) + t56 / 0.2e1 + t217 / 0.2e1 + t96 * mrSges(5,2) + (-m(5) * t40 - t82) * t167 + (-t90 / 0.2e1 - t215 / 0.2e1 + pkin(4) * t67 - t81 * mrSges(5,3) + t165 * t226 + t95 * t230 + (-m(5) * t81 - t112) * t167) * qJD(4)) * t177 + (t156 * t160 / 0.2e1 - t150 * t165 / 0.2e1 + t39 * mrSges(5,3) + t216 / 0.2e1 - t96 * mrSges(5,1) + t55 / 0.2e1 + (t91 / 0.2e1 - t80 * mrSges(5,3) + t156 * t225) * qJD(4) + (m(5) * (-qJD(4) * t80 + t39) + t83 - qJD(4) * t113) * t167) * t180 + (t222 + t238) * t154 / 0.2e1 + ((-t148 * t173 + t150 * t175) * mrSges(4,3) + m(4) * (t173 * t97 - t175 * t96)) * pkin(2) + (-t11 * t155 - t12 * t188 - t147 * t31 + t149 * t30) * mrSges(6,3) + (Ifges(6,5) * t155 + Ifges(7,5) * t111 - Ifges(6,6) * t188 + Ifges(7,6) * t110) * t148 / 0.2e1 - t188 * t32 / 0.2e1 + t118 * t158 + t161 * t36 + t155 * t33 / 0.2e1 - t149 * t54 / 0.2e1 - Ifges(4,5) * t150 - t147 * t53 / 0.2e1 - Ifges(4,6) * t148 + t123 * t8 + t124 * t37 + t110 * t5 / 0.2e1 + t111 * t6 / 0.2e1 + t68 * t114 + t62 * t115 / 0.2e1 + t63 * t116 / 0.2e1 + t104 * t47 + t105 * t46 + t95 * t106 - t100 * t107 / 0.2e1 - t101 * t108 / 0.2e1 - t96 * mrSges(4,1) - t97 * mrSges(4,2) + t86 * t85 + t87 * t84 + t74 * t34 / 0.2e1 + t38 * t75 + t19 * t76 / 0.2e1 + t18 * t77 / 0.2e1 + t73 * t35 / 0.2e1 + t64 * t42 / 0.2e1 + t65 * t43 / 0.2e1 + t66 * t41 + t22 * t51 + t44 * t13 + t45 * t14 + t21 * t50 + m(7) * (t10 * t21 + t123 * t38 + t124 * t66 + t2 * t45 + t22 * t9 + t3 * t44); 0.2e1 * t161 * t106 - t188 * t107 + t155 * t108 + t110 * t42 + t111 * t43 - t147 * t115 - t149 * t116 + 0.2e1 * t123 * t41 + 0.2e1 * t124 * t75 + 0.2e1 * t169 * t158 + t180 * t159 + t177 * t160 + t73 * t77 + t74 * t76 + (t180 * t165 + (0.2e1 * pkin(4) * t114 - t164) * t177) * qJD(4) + (t104 * t86 + t105 * t87 + t161 * t201) * t236 + (t123 * t124 + t21 * t45 + t22 * t44) * t235 + 0.2e1 * (t110 * t21 - t111 * t22 - t44 * t73 + t45 * t74) * mrSges(7,3) + 0.2e1 * (t104 * t149 - t105 * t147 - t155 * t86 - t188 * t87) * mrSges(6,3); m(4) * t202 + t148 * mrSges(4,1) + t110 * t13 + t111 * t14 - t147 * t85 - t149 * t84 - t188 * t47 + t155 * t46 + t177 * t83 + t180 * t82 + t73 * t50 + t74 * t51 - t137 + (t112 * t180 - t113 * t177) * qJD(4) + m(7) * (t10 * t73 + t110 * t3 + t111 * t2 + t74 * t9) + m(6) * (-t11 * t188 + t12 * t155 - t147 * t30 - t149 * t31) + m(5) * (t177 * t39 + t180 * t40 + (-t177 * t80 + t180 * t81) * qJD(4)); m(7) * (t110 * t22 + t111 * t21 + t44 * t74 + t45 * t73) + m(6) * (-t104 * t147 - t105 * t149 + t155 * t87 - t188 * t86); 0.2e1 * m(6) * (t147 * t188 - t149 * t155) + 0.2e1 * m(7) * (t110 * t74 + t111 * t73); (m(6) * (t11 * t174 + t12 * t172) + t174 * t47 + t172 * t46) * pkin(4) + t183 + m(7) * (t10 * t130 - t131 * t9 + t144 * t3 + t145 * t2) - t186 * Ifges(5,6) + t144 * t13 + t145 * t14 + t130 * t50 - t131 * t51 - t39 * mrSges(5,2) + t40 * mrSges(5,1) - Ifges(5,5) * t200 + t11 * mrSges(6,1) - t12 * mrSges(6,2) + t237; m(7) * (t130 * t45 - t131 * t44 + t144 * t22 + t145 * t21) - t87 * mrSges(6,2) + t86 * mrSges(6,1) + (-t219 + (-mrSges(5,1) * t180 + mrSges(5,2) * t177) * t167) * qJD(4) + (m(6) * (t172 * t87 + t174 * t86) + (-t147 * t172 + t149 * t174) * mrSges(6,3)) * pkin(4) + (t110 * t130 + t111 * t131 - t144 * t73 + t145 * t74) * mrSges(7,3) + t189 + t238; -t158 + m(7) * (-t110 * t131 + t111 * t130 + t144 * t74 + t145 * t73) + (-t147 * t174 - t149 * t172) * t230 + t184; 0.2e1 * m(7) * (t130 * t145 - t131 * t144) + 0.2e1 * t192; m(6) * t68 + m(7) * t38 + t36 + t8; m(6) * t201 + m(7) * t124 - t184; 0; 0; 0; t203 + t237; t189; -t41; t192; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
