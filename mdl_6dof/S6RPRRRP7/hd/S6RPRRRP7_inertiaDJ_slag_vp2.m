% Calculate time derivative of joint inertia matrix for
% S6RPRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP7_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP7_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP7_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP7_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:18:51
% EndTime: 2019-03-09 06:19:02
% DurationCPUTime: 5.40s
% Computational Cost: add. (6298->419), mult. (14134->594), div. (0->0), fcn. (13693->8), ass. (0->172)
t149 = sin(pkin(10));
t150 = cos(pkin(10));
t153 = sin(qJ(3));
t156 = cos(qJ(3));
t126 = t149 * t156 + t153 * t150;
t122 = t126 * qJD(3);
t234 = Ifges(7,4) + Ifges(6,5);
t238 = t122 * t234;
t232 = -Ifges(6,6) + Ifges(7,6);
t151 = sin(qJ(5));
t152 = sin(qJ(4));
t154 = cos(qJ(5));
t155 = cos(qJ(4));
t128 = t151 * t155 + t152 * t154;
t226 = qJD(4) + qJD(5);
t100 = t226 * t128;
t125 = t149 * t153 - t156 * t150;
t121 = t125 * qJD(3);
t127 = t151 * t152 - t154 * t155;
t33 = -t100 * t126 + t127 * t121;
t190 = qJD(5) * t151;
t195 = t152 * t121;
t197 = t126 * t155;
t198 = t126 * t152;
t193 = qJD(4) * t152;
t182 = t126 * t193;
t194 = t155 * t121;
t227 = -t194 - t182;
t34 = -t190 * t198 + (t197 * t226 - t195) * t154 + t227 * t151;
t237 = (-Ifges(6,4) + Ifges(7,5)) * t34 + (Ifges(6,1) + Ifges(7,1)) * t33 + t238;
t99 = t226 * t127;
t236 = t232 * t100 - t234 * t99;
t235 = -mrSges(6,1) - mrSges(7,1);
t233 = Ifges(7,2) + Ifges(6,3);
t192 = qJD(4) * t155;
t165 = t126 * t192 - t195;
t230 = m(6) + m(7);
t229 = -Ifges(5,5) * t194 + Ifges(5,3) * t122;
t213 = pkin(7) + qJ(2);
t134 = t213 * t149;
t135 = t213 * t150;
t228 = -t156 * t134 - t135 * t153;
t171 = mrSges(5,1) * t152 + mrSges(5,2) * t155;
t130 = t171 * qJD(4);
t76 = -t125 * qJD(2) + qJD(3) * t228;
t87 = pkin(3) * t122 + pkin(8) * t121;
t176 = -t152 * t76 + t155 * t87;
t184 = -pkin(2) * t150 - pkin(1);
t88 = pkin(3) * t125 - pkin(8) * t126 + t184;
t108 = -t153 * t134 + t135 * t156;
t98 = t155 * t108;
t15 = pkin(9) * t194 + pkin(4) * t122 + (-t98 + (pkin(9) * t126 - t88) * t152) * qJD(4) + t176;
t22 = -t108 * t193 + t152 * t87 + t155 * t76 + t88 * t192;
t19 = -pkin(9) * t165 + t22;
t54 = -t108 * t152 + t155 * t88;
t38 = pkin(4) * t125 - pkin(9) * t197 + t54;
t55 = t152 * t88 + t98;
t46 = -pkin(9) * t198 + t55;
t210 = t151 * t38 + t154 * t46;
t6 = -qJD(5) * t210 + t15 * t154 - t151 * t19;
t225 = 2 * m(6);
t224 = 2 * m(7);
t223 = 0.2e1 * pkin(4);
t222 = -2 * mrSges(4,3);
t220 = -0.2e1 * t228;
t219 = -pkin(9) - pkin(8);
t218 = t122 / 0.2e1;
t216 = -t126 / 0.2e1;
t205 = Ifges(5,4) * t152;
t136 = Ifges(5,2) * t155 + t205;
t214 = -t136 / 0.2e1;
t24 = mrSges(6,1) * t122 - mrSges(6,3) * t33;
t25 = -t122 * mrSges(7,1) + t33 * mrSges(7,2);
t212 = -t24 + t25;
t26 = -mrSges(6,2) * t122 - mrSges(6,3) * t34;
t27 = -mrSges(7,2) * t34 + mrSges(7,3) * t122;
t211 = t26 + t27;
t83 = t128 * t126;
t68 = -mrSges(6,2) * t125 - mrSges(6,3) * t83;
t71 = -mrSges(7,2) * t83 + mrSges(7,3) * t125;
t209 = t68 + t71;
t84 = t127 * t126;
t69 = mrSges(6,1) * t125 + mrSges(6,3) * t84;
t70 = -mrSges(7,1) * t125 - mrSges(7,2) * t84;
t208 = -t69 + t70;
t204 = Ifges(5,4) * t155;
t203 = Ifges(5,6) * t152;
t77 = qJD(2) * t126 + qJD(3) * t108;
t202 = t228 * t77;
t201 = t122 * Ifges(5,5);
t200 = t122 * Ifges(5,6);
t199 = t125 * Ifges(5,6);
t189 = qJD(5) * t154;
t188 = pkin(4) * t193;
t187 = pkin(4) * t190;
t186 = pkin(4) * t189;
t185 = t219 * t152;
t144 = -pkin(4) * t155 - pkin(3);
t183 = qJD(4) * t219;
t138 = t219 * t155;
t109 = -t151 * t138 - t154 * t185;
t180 = t109 * t190;
t179 = t127 * t190;
t57 = t100 * mrSges(6,1) - t99 * mrSges(6,2);
t56 = t100 * mrSges(7,1) + t99 * mrSges(7,3);
t178 = -(2 * Ifges(4,4)) - t203;
t110 = -t154 * t138 + t151 * t185;
t133 = t152 * t183;
t172 = t155 * t183;
t64 = -qJD(5) * t109 + t154 * t133 + t151 * t172;
t65 = qJD(5) * t110 + t151 * t133 - t154 * t172;
t177 = t109 * t65 + t110 * t64;
t175 = t122 * mrSges(4,1) - t121 * mrSges(4,2);
t78 = pkin(4) * t198 - t228;
t170 = Ifges(5,1) * t155 - t205;
t169 = -Ifges(5,2) * t152 + t204;
t20 = -t151 * t46 + t154 * t38;
t166 = t122 * t233 + t232 * t34 + t234 * t33;
t5 = t151 * t15 + t154 * t19 + t38 * t189 - t190 * t46;
t163 = -t56 - t57;
t161 = -pkin(5) * t100 - qJ(6) * t99 + qJD(6) * t128;
t160 = t235 * t65 + (-mrSges(6,2) + mrSges(7,3)) * t64 + t236;
t2 = qJ(6) * t122 + qJD(6) * t125 + t5;
t3 = -pkin(5) * t122 - t6;
t159 = t6 * mrSges(6,1) - t3 * mrSges(7,1) - t5 * mrSges(6,2) + t2 * mrSges(7,3) + t166;
t139 = qJD(6) + t186;
t157 = -mrSges(6,2) * t186 + t139 * mrSges(7,3) + t187 * t235;
t47 = pkin(4) * t165 + t77;
t148 = qJD(6) * mrSges(7,3);
t145 = Ifges(5,5) * t192;
t143 = -pkin(4) * t154 - pkin(5);
t141 = pkin(4) * t151 + qJ(6);
t137 = Ifges(5,1) * t152 + t204;
t132 = t170 * qJD(4);
t131 = t169 * qJD(4);
t106 = Ifges(6,1) * t128 - Ifges(6,4) * t127;
t105 = Ifges(7,1) * t128 + Ifges(7,5) * t127;
t104 = Ifges(6,4) * t128 - Ifges(6,2) * t127;
t103 = Ifges(7,5) * t128 + Ifges(7,3) * t127;
t102 = mrSges(6,1) * t127 + mrSges(6,2) * t128;
t101 = mrSges(7,1) * t127 - mrSges(7,3) * t128;
t91 = mrSges(5,1) * t125 - mrSges(5,3) * t197;
t90 = -mrSges(5,2) * t125 - mrSges(5,3) * t198;
t89 = pkin(5) * t127 - qJ(6) * t128 + t144;
t74 = Ifges(5,5) * t125 + t126 * t170;
t73 = t126 * t169 + t199;
t67 = -mrSges(5,2) * t122 - mrSges(5,3) * t165;
t66 = mrSges(5,1) * t122 - mrSges(5,3) * t227;
t61 = -Ifges(6,1) * t99 - Ifges(6,4) * t100;
t60 = -Ifges(7,1) * t99 + Ifges(7,5) * t100;
t59 = -Ifges(6,4) * t99 - Ifges(6,2) * t100;
t58 = -Ifges(7,5) * t99 + Ifges(7,3) * t100;
t52 = mrSges(5,1) * t165 + mrSges(5,2) * t227;
t51 = mrSges(6,1) * t83 - mrSges(6,2) * t84;
t50 = mrSges(7,1) * t83 + mrSges(7,3) * t84;
t49 = -t161 + t188;
t44 = -Ifges(6,1) * t84 - Ifges(6,4) * t83 + Ifges(6,5) * t125;
t43 = -Ifges(7,1) * t84 + Ifges(7,4) * t125 + Ifges(7,5) * t83;
t42 = -Ifges(6,4) * t84 - Ifges(6,2) * t83 + Ifges(6,6) * t125;
t41 = -Ifges(7,5) * t84 + Ifges(7,6) * t125 + Ifges(7,3) * t83;
t40 = Ifges(5,1) * t227 - Ifges(5,4) * t165 + t201;
t39 = Ifges(5,4) * t227 - Ifges(5,2) * t165 + t200;
t35 = pkin(5) * t83 + qJ(6) * t84 + t78;
t23 = -t55 * qJD(4) + t176;
t17 = -pkin(5) * t125 - t20;
t16 = qJ(6) * t125 + t210;
t13 = mrSges(6,1) * t34 + mrSges(6,2) * t33;
t12 = mrSges(7,1) * t34 - mrSges(7,3) * t33;
t9 = Ifges(6,4) * t33 - Ifges(6,2) * t34 + t122 * Ifges(6,6);
t8 = Ifges(7,5) * t33 + t122 * Ifges(7,6) + Ifges(7,3) * t34;
t7 = t34 * pkin(5) - t33 * qJ(6) + t84 * qJD(6) + t47;
t1 = [0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * qJD(2) * (t149 ^ 2 + t150 ^ 2) + (t122 * t232 + t8 - t9) * t83 + (t76 * t222 - t178 * t121 + ((2 * Ifges(4,2)) + Ifges(5,3) + t233) * t122 + t166 + t229) * t125 + (t41 - t42) * t34 + (t43 + t44) * t33 + t108 * t122 * t222 + t52 * t220 + (t16 * t2 + t17 * t3 + t35 * t7) * t224 - (mrSges(4,3) * t220 - t152 * t73 + t155 * t74) * t121 + 0.2e1 * m(5) * (t22 * t55 + t23 * t54 - t202) + 0.2e1 * m(4) * (t108 * t76 - t202) + 0.2e1 * t184 * t175 - (t237 + t238) * t84 + (t20 * t6 + t210 * t5 + t47 * t78) * t225 + 0.2e1 * t210 * t26 + (t155 * t40 - t152 * t39 - 0.2e1 * Ifges(4,1) * t121 + (Ifges(5,5) * t155 + t178) * t122 + (-t155 * t73 - t152 * t74 + t125 * (-Ifges(5,5) * t152 - Ifges(5,6) * t155)) * qJD(4) + 0.2e1 * (mrSges(4,3) + t171) * t77) * t126 + 0.2e1 * t20 * t24 + 0.2e1 * t17 * t25 + 0.2e1 * t16 * t27 + 0.2e1 * t35 * t12 + 0.2e1 * t7 * t50 + 0.2e1 * t47 * t51 + 0.2e1 * t54 * t66 + 0.2e1 * t55 * t67 + 0.2e1 * t5 * t68 + 0.2e1 * t6 * t69 + 0.2e1 * t3 * t70 + 0.2e1 * t2 * t71 + 0.2e1 * t78 * t13 + 0.2e1 * t22 * t90 + 0.2e1 * t23 * t91; t152 * t67 + t155 * t66 - t209 * t99 + t211 * t128 + t212 * t127 + t208 * t100 + (-t152 * t91 + t155 * t90) * qJD(4) + m(7) * (t100 * t17 + t127 * t3 + t128 * t2 - t16 * t99) + m(6) * (-t100 * t20 - t127 * t6 + t128 * t5 - t210 * t99) + m(5) * (t152 * t22 + t155 * t23 + (-t152 * t54 + t155 * t55) * qJD(4)) + t175; 0.2e1 * t230 * (t100 * t127 - t128 * t99); (-t100 * t16 - t17 * t99) * mrSges(7,2) + (-t100 * t210 + t20 * t99) * mrSges(6,3) + (t237 / 0.2e1 + t3 * mrSges(7,2) - t6 * mrSges(6,3) + t234 * t218) * t128 - t228 * t130 + m(7) * (t109 * t3 + t110 * t2 + t16 * t64 + t17 * t65 + t35 * t49 + t7 * t89) + (t145 + t236) * t125 / 0.2e1 + (-t5 * mrSges(6,3) - t2 * mrSges(7,2) + t8 / 0.2e1 - t9 / 0.2e1 + t232 * t218) * t127 - (t43 / 0.2e1 + t44 / 0.2e1) * t99 + (t58 / 0.2e1 - t59 / 0.2e1) * t83 + (t41 / 0.2e1 - t42 / 0.2e1) * t100 + (t103 / 0.2e1 - t104 / 0.2e1) * t34 + (t105 / 0.2e1 + t106 / 0.2e1) * t33 + t211 * t110 + t212 * t109 + t208 * t65 + t209 * t64 + m(6) * (-t109 * t6 + t110 * t5 + t144 * t47 - t20 * t65 + t210 * t64) - (t60 / 0.2e1 + t61 / 0.2e1) * t84 + (t126 * t132 / 0.2e1 - t121 * t137 / 0.2e1 + t22 * mrSges(5,3) + t39 / 0.2e1 - t77 * mrSges(5,1) + t200 / 0.2e1 + (t126 * t214 - t54 * mrSges(5,3) + t74 / 0.2e1) * qJD(4) + (t67 - qJD(4) * t91 + m(5) * (-qJD(4) * t54 + t22)) * pkin(8)) * t155 + (t131 * t216 - t121 * t214 - t23 * mrSges(5,3) + t40 / 0.2e1 + t77 * mrSges(5,2) + t201 / 0.2e1 + (t137 * t216 - t55 * mrSges(5,3) - t73 / 0.2e1 - t199 / 0.2e1 + (m(6) * t78 + t51) * pkin(4)) * qJD(4) + (-m(5) * t23 - t66 + (-m(5) * t55 - t90) * qJD(4)) * pkin(8)) * t152 + t49 * t50 + t35 * t56 + (-m(5) * t77 - t52) * pkin(3) - t76 * mrSges(4,2) - t77 * mrSges(4,1) + t78 * t57 + t89 * t12 + t7 * t101 + t47 * t102 - Ifges(4,5) * t121 - Ifges(4,6) * t122 + t144 * t13; t230 * (t100 * t109 - t110 * t99 + t127 * t65 + t64 * t128); -0.2e1 * pkin(3) * t130 + 0.2e1 * t49 * t101 + t155 * t131 + t152 * t132 + 0.2e1 * t144 * t57 + 0.2e1 * t89 * t56 - (t105 + t106) * t99 + (t60 + t61) * t128 + (t58 - t59) * t127 + (t103 - t104) * t100 + (t155 * t137 + (t102 * t223 - t136) * t152) * qJD(4) + (t144 * t188 + t177) * t225 + (t49 * t89 + t177) * t224 + 0.2e1 * (mrSges(6,3) + mrSges(7,2)) * (-t100 * t110 - t109 * t99 - t127 * t64 + t128 * t65); t159 - Ifges(5,5) * t182 - t165 * Ifges(5,6) - t22 * mrSges(5,2) + t23 * mrSges(5,1) + t139 * t71 + t141 * t27 + t143 * t25 + m(7) * (t139 * t16 + t141 * t2 + t143 * t3) + ((m(6) * t6 + t24 + (m(6) * t210 + t68) * qJD(5)) * t154 + (m(6) * t5 + t26 + (-m(6) * t20 + m(7) * t17 + t208) * qJD(5)) * t151) * pkin(4) + t229; -t130 + m(7) * (t100 * t143 + t128 * t139 - t141 * t99) + (m(6) * (-t100 * t154 + t128 * t189 - t151 * t99 + t179) / 0.2e1 + m(7) * t179 / 0.2e1) * t223 + t163; m(7) * (t110 * t139 + t141 * t64 + t143 * t65) + t145 + (-t203 + (-mrSges(5,1) * t155 + mrSges(5,2) * t152) * pkin(8)) * qJD(4) + (-t100 * t141 - t127 * t139 - t143 * t99) * mrSges(7,2) + (t128 * mrSges(7,2) * t190 + m(6) * (t110 * t189 + t151 * t64 - t154 * t65 + t180) + m(7) * t180 + (-t151 * t100 + t154 * t99 + (-t127 * t154 + t128 * t151) * qJD(5)) * mrSges(6,3)) * pkin(4) + t160; 0.2e1 * m(7) * (t139 * t141 + t143 * t187) + 0.2e1 * t157; t159 - pkin(5) * t25 + qJ(6) * t27 + qJD(6) * t71 + m(7) * (-pkin(5) * t3 + qJ(6) * t2 + qJD(6) * t16); m(7) * t161 + t163; m(7) * (-pkin(5) * t65 + qJ(6) * t64 + qJD(6) * t110) + (pkin(5) * t99 - qJ(6) * t100 - qJD(6) * t127) * mrSges(7,2) + t160; t148 + m(7) * (-pkin(5) * t187 + qJ(6) * t139 + qJD(6) * t141) + t157; 0.2e1 * m(7) * qJ(6) * qJD(6) + 0.2e1 * t148; m(7) * t3 + t25; m(7) * t100; m(7) * t65 - t99 * mrSges(7,2); m(7) * t187; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
