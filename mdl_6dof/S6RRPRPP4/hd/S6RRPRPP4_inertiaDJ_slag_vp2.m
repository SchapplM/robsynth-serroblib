% Calculate time derivative of joint inertia matrix for
% S6RRPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
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
% Datum: 2019-03-09 10:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPP4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:58:53
% EndTime: 2019-03-09 09:59:00
% DurationCPUTime: 3.62s
% Computational Cost: add. (2783->392), mult. (6247->546), div. (0->0), fcn. (4845->6), ass. (0->174)
t237 = Ifges(7,4) + Ifges(6,5);
t236 = Ifges(6,6) - Ifges(7,6);
t146 = cos(qJ(4));
t144 = sin(qJ(4));
t147 = cos(qJ(2));
t183 = qJD(4) * t147;
t177 = t144 * t183;
t145 = sin(qJ(2));
t187 = qJD(2) * t145;
t155 = t146 * t187 + t177;
t143 = sin(pkin(9));
t194 = cos(pkin(9));
t170 = t194 * t146;
t185 = qJD(4) * t144;
t96 = -qJD(4) * t170 + t143 * t185;
t157 = -t143 * t146 - t144 * t194;
t97 = t157 * qJD(4);
t235 = t236 * t96 + t237 * t97;
t234 = Ifges(7,2) + Ifges(5,3) + Ifges(6,3);
t233 = m(4) * pkin(7);
t232 = m(6) + m(7);
t231 = -mrSges(6,1) - mrSges(7,1);
t230 = mrSges(6,3) + mrSges(7,2);
t105 = t143 * t144 - t170;
t229 = t105 * t97 - t157 * t96;
t184 = qJD(4) * t146;
t228 = -mrSges(5,1) * t185 - mrSges(5,2) * t184;
t148 = -pkin(2) - pkin(8);
t193 = qJ(3) * t145;
t103 = t147 * t148 - pkin(1) - t193;
t186 = qJD(2) * t147;
t216 = pkin(3) + pkin(7);
t115 = t216 * t186;
t122 = t216 * t145;
t167 = pkin(2) * t187 - qJD(3) * t145;
t81 = (pkin(8) * t145 - qJ(3) * t147) * qJD(2) + t167;
t21 = -t103 * t185 + t144 * t115 + t122 * t184 + t146 * t81;
t100 = t146 * t115;
t63 = t146 * t103 + t144 * t122;
t22 = -qJD(4) * t63 - t144 * t81 + t100;
t227 = t21 * t144 + t22 * t146;
t226 = t143 * t96 - t194 * t97;
t181 = t146 * qJD(5);
t188 = qJ(5) - t148;
t153 = t185 * t188 - t181;
t169 = t188 * t146;
t85 = -qJD(4) * t169 - t144 * qJD(5);
t36 = t143 * t85 - t153 * t194;
t37 = t143 * t153 + t194 * t85;
t116 = t188 * t144;
t70 = -t116 * t143 + t169 * t194;
t71 = -t116 * t194 - t143 * t169;
t225 = t105 * t36 - t157 * t37 - t70 * t97 - t96 * t71;
t212 = pkin(4) * t143;
t131 = qJ(6) + t212;
t178 = t194 * pkin(4);
t133 = -t178 - pkin(5);
t224 = qJD(6) * t157 + t131 * t96 + t133 * t97;
t223 = m(7) * t131 + mrSges(7,3);
t222 = 2 * m(6);
t221 = 0.2e1 * m(7);
t220 = -0.2e1 * pkin(1);
t95 = -qJ(3) * t186 + t167;
t219 = 0.2e1 * t95;
t161 = -pkin(2) * t147 - t193;
t117 = -pkin(1) + t161;
t218 = -0.2e1 * t117;
t217 = m(6) * pkin(4);
t168 = qJ(5) * t147 - t103;
t7 = pkin(4) * t186 + t100 + t168 * t184 + (-qJ(5) * t187 - qJD(4) * t122 + qJD(5) * t147 - t81) * t144;
t9 = qJ(5) * t155 - t147 * t181 + t21;
t4 = t143 * t7 + t194 * t9;
t213 = -t147 / 0.2e1;
t141 = t147 * pkin(7);
t175 = t144 * t187;
t46 = -t143 * t175 - t147 * t97 + t170 * t187;
t26 = mrSges(7,2) * t46 + mrSges(7,3) * t186;
t27 = -mrSges(6,2) * t186 + mrSges(6,3) * t46;
t211 = t26 + t27;
t47 = t105 * t183 - t157 * t187;
t28 = mrSges(6,1) * t186 - mrSges(6,3) * t47;
t29 = -mrSges(7,1) * t186 + t47 * mrSges(7,2);
t210 = -t28 + t29;
t108 = t146 * t122;
t45 = pkin(4) * t145 + t144 * t168 + t108;
t189 = t146 * t147;
t53 = -qJ(5) * t189 + t63;
t17 = t143 * t45 + t194 * t53;
t191 = t144 * t147;
t86 = t143 * t191 - t147 * t170;
t75 = mrSges(7,2) * t86 + mrSges(7,3) * t145;
t76 = -mrSges(6,2) * t145 + mrSges(6,3) * t86;
t209 = t75 + t76;
t87 = t157 * t147;
t77 = mrSges(6,1) * t145 - mrSges(6,3) * t87;
t78 = -mrSges(7,1) * t145 + mrSges(7,2) * t87;
t208 = -t77 + t78;
t205 = Ifges(5,4) * t144;
t204 = Ifges(5,4) * t146;
t164 = Ifges(5,1) * t144 + t204;
t198 = t145 * Ifges(5,5);
t84 = -t147 * t164 + t198;
t199 = t144 * t84;
t163 = Ifges(5,2) * t146 + t205;
t83 = t145 * Ifges(5,6) - t147 * t163;
t197 = t146 * t83;
t121 = Ifges(5,1) * t146 - t205;
t192 = t144 * t121;
t120 = -Ifges(5,2) * t144 + t204;
t190 = t146 * t120;
t134 = t144 * pkin(4) + qJ(3);
t123 = t147 * pkin(3) + t141;
t128 = pkin(4) * t184 + qJD(3);
t98 = pkin(4) * t189 + t123;
t176 = t146 * t183;
t19 = -t46 * mrSges(6,1) + t47 * mrSges(6,2);
t55 = -t96 * mrSges(6,1) + t97 * mrSges(6,2);
t18 = -t46 * mrSges(7,1) - t47 * mrSges(7,3);
t54 = -t96 * mrSges(7,1) - t97 * mrSges(7,3);
t173 = t36 * t70 + t71 * t37;
t165 = mrSges(5,1) * t146 - mrSges(5,2) * t144;
t162 = -Ifges(5,5) * t144 - Ifges(5,6) * t146;
t62 = -t103 * t144 + t108;
t160 = t144 * t62 - t146 * t63;
t112 = mrSges(5,1) * t145 + mrSges(5,3) * t191;
t113 = -mrSges(5,2) * t145 - mrSges(5,3) * t189;
t159 = -t144 * t112 + t146 * t113;
t3 = -t143 * t9 + t194 * t7;
t16 = -t143 * t53 + t194 * t45;
t154 = t175 - t176;
t1 = qJ(6) * t186 + qJD(6) * t145 + t4;
t14 = qJ(6) * t145 + t17;
t15 = -t145 * pkin(5) - t16;
t2 = -pkin(5) * t186 - t3;
t152 = -t1 * t157 + t105 * t2 - t14 * t96 - t15 * t97;
t151 = t105 * t3 + t157 * t4 - t16 * t97 + t17 * t96;
t150 = Ifges(5,5) * t175 + t155 * Ifges(5,6) + t234 * t186 + t236 * t46 + t237 * t47;
t73 = -pkin(4) * t177 + (-pkin(4) * t146 - t216) * t187;
t119 = mrSges(5,1) * t144 + mrSges(5,2) * t146;
t114 = t216 * t187;
t111 = t164 * qJD(4);
t110 = t163 * qJD(4);
t109 = t165 * qJD(4);
t101 = t165 * t147;
t80 = mrSges(5,1) * t186 - mrSges(5,3) * t154;
t79 = -mrSges(5,2) * t186 + mrSges(5,3) * t155;
t69 = -Ifges(6,1) * t105 + Ifges(6,4) * t157;
t68 = -Ifges(7,1) * t105 - Ifges(7,5) * t157;
t67 = -Ifges(6,4) * t105 + Ifges(6,2) * t157;
t66 = -Ifges(7,5) * t105 - Ifges(7,3) * t157;
t65 = -mrSges(6,1) * t157 - mrSges(6,2) * t105;
t64 = -mrSges(7,1) * t157 + mrSges(7,3) * t105;
t61 = -mrSges(5,1) * t155 + mrSges(5,2) * t154;
t60 = -pkin(5) * t157 + qJ(6) * t105 + t134;
t59 = Ifges(6,1) * t97 + Ifges(6,4) * t96;
t58 = Ifges(7,1) * t97 - Ifges(7,5) * t96;
t57 = Ifges(6,4) * t97 + Ifges(6,2) * t96;
t56 = Ifges(7,5) * t97 - Ifges(7,3) * t96;
t52 = -t121 * t183 + (t147 * Ifges(5,5) + t145 * t164) * qJD(2);
t51 = -t120 * t183 + (t147 * Ifges(5,6) + t145 * t163) * qJD(2);
t50 = -mrSges(6,1) * t86 + mrSges(6,2) * t87;
t49 = -mrSges(7,1) * t86 - mrSges(7,3) * t87;
t35 = Ifges(6,1) * t87 + Ifges(6,4) * t86 + Ifges(6,5) * t145;
t34 = Ifges(7,1) * t87 + Ifges(7,4) * t145 - Ifges(7,5) * t86;
t33 = Ifges(6,4) * t87 + Ifges(6,2) * t86 + Ifges(6,6) * t145;
t32 = Ifges(7,5) * t87 + Ifges(7,6) * t145 - Ifges(7,3) * t86;
t25 = -pkin(5) * t86 - qJ(6) * t87 + t98;
t24 = -pkin(5) * t96 - qJ(6) * t97 + qJD(6) * t105 + t128;
t13 = Ifges(6,1) * t47 + Ifges(6,4) * t46 + Ifges(6,5) * t186;
t12 = Ifges(7,1) * t47 + Ifges(7,4) * t186 - Ifges(7,5) * t46;
t11 = Ifges(6,4) * t47 + Ifges(6,2) * t46 + Ifges(6,6) * t186;
t10 = Ifges(7,5) * t47 + Ifges(7,6) * t186 - Ifges(7,3) * t46;
t5 = -pkin(5) * t46 - qJ(6) * t47 - qJD(6) * t87 + t73;
t6 = [(t34 + t35) * t47 + (t33 - t32) * t46 + 0.2e1 * m(5) * (-t114 * t123 + t21 * t63 + t22 * t62) + (mrSges(4,2) * t219 - t144 * t52 - t146 * t51 + (t144 * t83 + (-t84 - t198) * t146) * qJD(4)) * t147 + (-0.2e1 * mrSges(4,3) * t95 + t150) * t145 + m(4) * t117 * t219 + 0.2e1 * t25 * t18 + 0.2e1 * t14 * t26 + 0.2e1 * t17 * t27 + 0.2e1 * t16 * t28 + 0.2e1 * t15 * t29 + 0.2e1 * t5 * t49 + 0.2e1 * t73 * t50 + 0.2e1 * t1 * t75 + 0.2e1 * t4 * t76 + 0.2e1 * t3 * t77 + 0.2e1 * t2 * t78 + 0.2e1 * t63 * t79 + 0.2e1 * t62 * t80 + 0.2e1 * t98 * t19 + (t1 * t14 + t15 * t2 + t25 * t5) * t221 + (t16 * t3 + t17 * t4 + t73 * t98) * t222 + ((mrSges(3,1) * t220 + mrSges(4,2) * t218 + t199 + t197 + 0.2e1 * (-Ifges(4,6) - Ifges(3,4)) * t145) * t145 + (mrSges(3,2) * t220 + mrSges(4,3) * t218 + (0.2e1 * Ifges(3,4) + 0.2e1 * Ifges(4,6) + t162) * t147 + t237 * t87 + t236 * t86 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) + (2 * Ifges(4,2)) - (2 * Ifges(4,3)) + t234) * t145) * t147) * qJD(2) + 0.2e1 * t22 * t112 + 0.2e1 * t21 * t113 - 0.2e1 * t114 * t101 + 0.2e1 * t123 * t61 + (t11 - t10) * t86 + (t12 + t13) * t87; (t68 / 0.2e1 + t69 / 0.2e1) * t47 + (-t66 / 0.2e1 + t67 / 0.2e1) * t46 + m(7) * (t1 * t71 + t14 * t37 + t15 * t36 + t2 * t70 + t24 * t25 + t5 * t60) + m(6) * (t128 * t98 + t134 * t73 - t16 * t36 + t17 * t37 - t3 * t70 + t4 * t71) + (-t22 * mrSges(5,3) - t110 * t213 + t148 * t80 + t52 / 0.2e1) * t146 + (-t21 * mrSges(5,3) - t111 * t213 + t148 * t79 - t51 / 0.2e1) * t144 - (t10 / 0.2e1 - t11 / 0.2e1) * t157 + (t161 * t233 + (-pkin(2) * mrSges(4,1) + Ifges(5,5) * t146 / 0.2e1 - Ifges(5,6) * t144 / 0.2e1 - Ifges(4,4) + Ifges(3,5) - (Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1) * t157 + (-Ifges(7,4) / 0.2e1 - Ifges(6,5) / 0.2e1) * t105 + (-mrSges(3,1) + mrSges(4,2)) * pkin(7)) * t147 + (t192 / 0.2e1 - qJ(3) * mrSges(4,1) + t190 / 0.2e1 + Ifges(4,5) - Ifges(3,6) + (mrSges(3,2) - mrSges(4,3)) * pkin(7)) * t145) * qJD(2) + (t58 / 0.2e1 + t59 / 0.2e1) * t87 + t208 * t36 + t209 * t37 + t210 * t70 + t211 * t71 + (t57 / 0.2e1 - t56 / 0.2e1) * t86 + t151 * mrSges(6,3) - t152 * mrSges(7,2) - (-t34 / 0.2e1 - t35 / 0.2e1) * t97 + (t33 / 0.2e1 - t32 / 0.2e1) * t96 + (-t12 / 0.2e1 - t13 / 0.2e1) * t105 + (-t199 / 0.2e1 - t197 / 0.2e1 + (-t146 * t121 / 0.2e1 + t144 * t120 / 0.2e1) * t147 + t160 * mrSges(5,3) + (-m(5) * t160 + t159) * t148) * qJD(4) + (m(4) * t141 + m(5) * t123 + t147 * mrSges(4,1) + t101) * qJD(3) + t24 * t49 + t25 * t54 + t60 * t18 + qJ(3) * t61 + t5 * t64 + t73 * t65 + t98 * t55 + (qJD(4) * t162 + t235) * t145 / 0.2e1 - t114 * t119 + t123 * t109 + t128 * t50 + t134 * t19 + m(5) * (-qJ(3) * t114 + t148 * t227); 0.2e1 * qJ(3) * t109 + t144 * t110 - t146 * t111 + 0.2e1 * t128 * t65 + 0.2e1 * t134 * t55 + 0.2e1 * t24 * t64 + 0.2e1 * t60 * t54 - (-t68 - t69) * t97 + (t67 - t66) * t96 - (t56 - t57) * t157 + (-t58 - t59) * t105 + (-t190 - t192) * qJD(4) + (t24 * t60 + t173) * t221 + (t128 * t134 + t173) * t222 + 0.2e1 * (mrSges(4,3) + t119 + (m(4) + m(5)) * qJ(3)) * qJD(3) - 0.2e1 * t230 * t225; t144 * t79 + t146 * t80 - t208 * t97 - t209 * t96 - t211 * t157 + t210 * t105 + t159 * qJD(4) + (mrSges(4,1) + t233) * t186 + m(7) * t152 - m(6) * t151 + m(5) * (-qJD(4) * t160 + t227); t225 * t232 + 0.2e1 * t229 * t230; -0.2e1 * t232 * t229; t150 + t1 * mrSges(7,3) - t2 * mrSges(7,1) + t3 * mrSges(6,1) - t4 * mrSges(6,2) + t28 * t178 - Ifges(5,5) * t176 + (t143 * t4 + t194 * t3) * t217 + m(7) * (qJD(6) * t14 + t1 * t131 + t133 * t2) + t27 * t212 - t21 * mrSges(5,2) + t22 * mrSges(5,1) + qJD(6) * t75 + t131 * t26 + t133 * t29; m(7) * qJD(6) * t71 - Ifges(5,5) * t185 - Ifges(5,6) * t184 + (t143 * t217 - mrSges(6,2) + t223) * t37 + (m(7) * t133 - t194 * t217 + t231) * t36 + t226 * mrSges(6,3) * pkin(4) + t228 * t148 + t224 * mrSges(7,2) + t235; -t226 * t217 - m(7) * t224 - t231 * t97 + (mrSges(6,2) - mrSges(7,3)) * t96 + t228; 0.2e1 * t223 * qJD(6); m(6) * t73 + m(7) * t5 + t18 + t19; m(6) * t128 + m(7) * t24 + t54 + t55; 0; 0; 0; m(7) * t2 + t29; m(7) * t36 + t97 * mrSges(7,2); -m(7) * t97; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
