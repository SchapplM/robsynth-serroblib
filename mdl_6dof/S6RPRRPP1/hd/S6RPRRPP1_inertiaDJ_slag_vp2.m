% Calculate time derivative of joint inertia matrix for
% S6RPRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-03-09 04:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:27:55
% EndTime: 2019-03-09 04:28:03
% DurationCPUTime: 3.70s
% Computational Cost: add. (2676->414), mult. (6505->583), div. (0->0), fcn. (5248->8), ass. (0->175)
t137 = sin(pkin(10));
t141 = cos(qJ(4));
t186 = cos(pkin(10));
t157 = t186 * t141;
t139 = sin(qJ(4));
t158 = t186 * t139;
t140 = sin(qJ(3));
t175 = qJD(4) * t140;
t142 = cos(qJ(3));
t177 = qJD(3) * t142;
t163 = t141 * t177;
t165 = t139 * t175;
t213 = t165 - t163;
t44 = t213 * t137 - t157 * t175 - t158 * t177;
t104 = t137 * t141 + t158;
t183 = t137 * t139;
t150 = t157 - t183;
t45 = -t104 * t175 + t150 * t177;
t15 = -t44 * mrSges(7,1) - t45 * mrSges(7,3);
t16 = -t44 * mrSges(6,1) + t45 * mrSges(6,2);
t152 = -t15 - t16;
t174 = qJD(4) * t141;
t146 = t139 * t177 + t140 * t174;
t62 = mrSges(5,1) * t146 - t213 * mrSges(5,2);
t222 = -t62 + t152;
t221 = mrSges(6,3) + mrSges(7,2);
t220 = Ifges(7,4) + Ifges(6,5);
t219 = Ifges(7,6) - Ifges(6,6);
t218 = -Ifges(7,2) - Ifges(5,3) - Ifges(6,3);
t170 = -cos(pkin(9)) * pkin(1) - pkin(2);
t217 = 0.2e1 * t170;
t216 = m(6) + m(7);
t96 = t104 * qJD(4);
t97 = t150 * qJD(4);
t55 = t96 * mrSges(7,1) - t97 * mrSges(7,3);
t56 = t96 * mrSges(6,1) + t97 * mrSges(6,2);
t215 = -t55 - t56;
t179 = t139 ^ 2 + t141 ^ 2;
t203 = pkin(4) * t137;
t126 = qJ(6) + t203;
t212 = m(7) * t126 + mrSges(7,3);
t211 = 2 * m(5);
t210 = 2 * m(6);
t209 = 0.2e1 * m(7);
t128 = sin(pkin(9)) * pkin(1) + pkin(7);
t208 = 0.2e1 * t128;
t207 = m(6) * pkin(4);
t194 = Ifges(5,4) * t139;
t119 = Ifges(5,2) * t141 + t194;
t206 = -t119 / 0.2e1;
t205 = -t139 / 0.2e1;
t101 = -pkin(3) * t142 - t140 * pkin(8) + t170;
t180 = t141 * t142;
t111 = t128 * t180;
t173 = qJD(5) * t141;
t202 = pkin(8) * t142;
t204 = pkin(3) * t140;
t112 = (-t202 + t204) * qJD(3);
t178 = qJD(3) * t140;
t168 = t128 * t178;
t187 = t141 * t112 + t139 * t168;
t12 = -t140 * t173 + (pkin(4) * t140 - qJ(5) * t180) * qJD(3) + (-t111 + (qJ(5) * t140 - t101) * t139) * qJD(4) + t187;
t181 = t140 * t141;
t195 = t101 * t174 + t139 * t112;
t19 = (-qJ(5) * qJD(4) - qJD(3) * t128) * t181 + (-qJD(5) * t140 + (-qJ(5) * qJD(3) - qJD(4) * t128) * t142) * t139 + t195;
t4 = t137 * t12 + t186 * t19;
t201 = t97 * mrSges(7,2);
t200 = -qJ(5) - pkin(8);
t28 = mrSges(7,2) * t44 + mrSges(7,3) * t178;
t29 = -mrSges(6,2) * t178 + mrSges(6,3) * t44;
t199 = t28 + t29;
t30 = mrSges(6,1) * t178 - mrSges(6,3) * t45;
t31 = -mrSges(7,1) * t178 + t45 * mrSges(7,2);
t198 = -t30 + t31;
t87 = t141 * t101;
t46 = -qJ(5) * t181 + t87 + (-t128 * t139 - pkin(4)) * t142;
t182 = t139 * t140;
t70 = t139 * t101 + t111;
t54 = -qJ(5) * t182 + t70;
t18 = t137 * t46 + t186 * t54;
t84 = t104 * t140;
t74 = -t84 * mrSges(7,2) - mrSges(7,3) * t142;
t75 = mrSges(6,2) * t142 - t84 * mrSges(6,3);
t197 = t74 + t75;
t85 = t150 * t140;
t76 = -mrSges(6,1) * t142 - t85 * mrSges(6,3);
t77 = mrSges(7,1) * t142 + t85 * mrSges(7,2);
t196 = t76 - t77;
t193 = Ifges(5,4) * t141;
t192 = Ifges(5,5) * t139;
t191 = Ifges(5,6) * t141;
t190 = t142 * Ifges(5,6);
t23 = -qJD(4) * t70 + t187;
t189 = t23 * t139;
t188 = -mrSges(5,1) * t141 + mrSges(5,2) * t139 - mrSges(4,1);
t185 = qJD(6) * t85;
t184 = t128 * t142;
t95 = pkin(4) * t182 + t140 * t128;
t176 = qJD(4) * t139;
t172 = pkin(4) * t176;
t171 = pkin(8) * t176;
t113 = t128 * t177;
t71 = pkin(4) * t146 + t113;
t130 = -pkin(4) * t141 - pkin(3);
t169 = t186 * pkin(4);
t166 = t140 * t177;
t162 = t142 * t176;
t159 = qJD(4) * t200;
t145 = -qJD(5) * t139 + t141 * t159;
t94 = t139 * t159 + t173;
t47 = t137 * t94 - t145 * t186;
t48 = t137 * t145 + t186 * t94;
t118 = t200 * t141;
t72 = -t118 * t137 - t158 * t200;
t73 = -t118 * t186 + t183 * t200;
t161 = t47 * t72 + t73 * t48;
t160 = Ifges(5,6) * t139 + (2 * Ifges(4,4));
t22 = (-t141 * t178 - t162) * t128 + t195;
t69 = -t139 * t184 + t87;
t156 = -qJD(4) * t69 + t22;
t155 = mrSges(5,1) * t139 + mrSges(5,2) * t141;
t154 = Ifges(5,1) * t141 - t194;
t120 = Ifges(5,1) * t139 + t193;
t153 = -Ifges(5,2) * t139 + t193;
t3 = t12 * t186 - t137 * t19;
t17 = -t137 * t54 + t186 * t46;
t151 = -t72 * t44 + t73 * t45 + t47 * t84 + t48 * t85;
t148 = -Ifges(5,5) * t163 + t178 * t218 + t219 * t44 - t220 * t45;
t134 = Ifges(5,5) * t174;
t129 = -t169 - pkin(5);
t110 = -mrSges(5,1) * t142 - mrSges(5,3) * t181;
t109 = mrSges(5,2) * t142 - mrSges(5,3) * t182;
t108 = t154 * qJD(4);
t107 = t153 * qJD(4);
t106 = t155 * qJD(4);
t100 = t155 * t140;
t93 = Ifges(7,4) * t97;
t92 = Ifges(6,5) * t97;
t91 = Ifges(6,6) * t96;
t90 = Ifges(7,6) * t96;
t83 = -Ifges(5,5) * t142 + t140 * t154;
t82 = t140 * t153 - t190;
t79 = -mrSges(5,2) * t178 - mrSges(5,3) * t146;
t78 = mrSges(5,1) * t178 + mrSges(5,3) * t213;
t68 = Ifges(6,1) * t104 + Ifges(6,4) * t150;
t67 = Ifges(7,1) * t104 - Ifges(7,5) * t150;
t66 = Ifges(6,4) * t104 + Ifges(6,2) * t150;
t65 = Ifges(7,5) * t104 - Ifges(7,3) * t150;
t64 = -mrSges(6,1) * t150 + mrSges(6,2) * t104;
t63 = -mrSges(7,1) * t150 - mrSges(7,3) * t104;
t61 = -pkin(5) * t150 - qJ(6) * t104 + t130;
t60 = Ifges(6,1) * t97 - Ifges(6,4) * t96;
t59 = Ifges(7,1) * t97 + Ifges(7,5) * t96;
t58 = Ifges(6,4) * t97 - Ifges(6,2) * t96;
t57 = Ifges(7,5) * t97 + Ifges(7,3) * t96;
t53 = -t120 * t175 + (Ifges(5,5) * t140 + t142 * t154) * qJD(3);
t52 = -t119 * t175 + (Ifges(5,6) * t140 + t142 * t153) * qJD(3);
t51 = mrSges(6,1) * t84 + mrSges(6,2) * t85;
t50 = mrSges(7,1) * t84 - mrSges(7,3) * t85;
t35 = Ifges(6,1) * t85 - Ifges(6,4) * t84 - Ifges(6,5) * t142;
t34 = Ifges(7,1) * t85 - Ifges(7,4) * t142 + Ifges(7,5) * t84;
t33 = Ifges(6,4) * t85 - Ifges(6,2) * t84 - Ifges(6,6) * t142;
t32 = Ifges(7,5) * t85 - Ifges(7,6) * t142 + Ifges(7,3) * t84;
t27 = pkin(5) * t96 - qJ(6) * t97 - qJD(6) * t104 + t172;
t26 = pkin(5) * t84 - qJ(6) * t85 + t95;
t13 = t142 * pkin(5) - t17;
t11 = -qJ(6) * t142 + t18;
t10 = Ifges(6,1) * t45 + Ifges(6,4) * t44 + Ifges(6,5) * t178;
t9 = Ifges(7,1) * t45 + Ifges(7,4) * t178 - Ifges(7,5) * t44;
t8 = Ifges(6,4) * t45 + Ifges(6,2) * t44 + Ifges(6,6) * t178;
t7 = Ifges(7,5) * t45 + Ifges(7,6) * t178 - Ifges(7,3) * t44;
t5 = -pkin(5) * t44 - qJ(6) * t45 - t185 + t71;
t2 = -pkin(5) * t178 - t3;
t1 = qJ(6) * t178 - qJD(6) * t142 + t4;
t6 = [(t34 + t35) * t45 + (t33 - t32) * t44 + 0.2e1 * t17 * t30 + 0.2e1 * t13 * t31 + 0.2e1 * t26 * t15 + 0.2e1 * t11 * t28 + 0.2e1 * t18 * t29 + (t1 * t11 + t13 * t2 + t26 * t5) * t209 + (t17 * t3 + t18 * t4 + t71 * t95) * t210 + (t70 * t22 + t69 * t23) * t211 + ((mrSges(4,2) * t217 + t100 * t208 - t139 * t82 + t141 * t83 + t142 * t160) * qJD(3) + t148) * t142 + (t62 * t208 - t139 * t52 + t141 * t53 + (-t139 * t83 - t141 * t82 - t142 * (-t191 - t192)) * qJD(4) + (mrSges(4,1) * t217 + (Ifges(5,5) * t141 - t160) * t140 + t220 * t85 + t219 * t84 + (t128 ^ 2 * t211 + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) + t218) * t142) * qJD(3)) * t140 + (t9 + t10) * t85 + (t7 - t8) * t84 + 0.2e1 * t5 * t50 + 0.2e1 * t71 * t51 + 0.2e1 * t1 * t74 + 0.2e1 * t4 * t75 + 0.2e1 * t3 * t76 + 0.2e1 * t2 * t77 + 0.2e1 * t69 * t78 + 0.2e1 * t70 * t79 + 0.2e1 * t95 * t16 + 0.2e1 * t22 * t109 + 0.2e1 * t23 * t110; t199 * t85 + t198 * t84 + t197 * t45 + t196 * t44 + ((t109 * t141 - t110 * t139) * qJD(3) + t222) * t142 + (-t139 * t78 + t141 * t79 + (-t109 * t139 - t110 * t141) * qJD(4) + (t100 + t50 + t51) * qJD(3)) * t140 + m(6) * (-t142 * t71 + t17 * t44 + t178 * t95 + t18 * t45 - t3 * t84 + t4 * t85) + m(7) * (t1 * t85 + t11 * t45 - t13 * t44 - t142 * t5 + t178 * t26 + t2 * t84) + m(5) * ((-t139 * t69 + t141 * t70 - t184) * t177 + (t168 - t189 + t141 * t22 + (-t139 * t70 - t141 * t69) * qJD(4)) * t140); 0.2e1 * m(5) * (-0.1e1 + t179) * t166 + 0.2e1 * t216 * (-t84 * t44 + t85 * t45 - t166); (t59 / 0.2e1 + t60 / 0.2e1) * t85 + (t34 / 0.2e1 + t35 / 0.2e1) * t97 + (t57 / 0.2e1 - t58 / 0.2e1) * t84 + (t67 / 0.2e1 + t68 / 0.2e1) * t45 + (t32 / 0.2e1 - t33 / 0.2e1) * t96 + (-t65 / 0.2e1 + t66 / 0.2e1) * t44 + m(7) * (t1 * t73 + t11 * t48 + t13 * t47 + t2 * t72 + t26 * t27 + t5 * t61) - t196 * t47 + t197 * t48 + t198 * t72 + t199 * t73 + m(5) * (-pkin(3) * t113 - pkin(8) * t189 - t171 * t70) + (-t134 / 0.2e1 - t92 / 0.2e1 + t91 / 0.2e1 - t93 / 0.2e1 - t90 / 0.2e1 + (t128 * t188 + Ifges(4,5)) * qJD(3)) * t142 + (qJD(4) * t83 / 0.2e1 + t120 * t177 / 0.2e1 + t52 / 0.2e1 + t156 * mrSges(5,3) + (m(5) * t156 - qJD(4) * t110 + t79) * pkin(8)) * t141 + m(6) * (t130 * t71 - t17 * t47 + t172 * t95 + t18 * t48 - t3 * t72 + t4 * t73) + (-t104 * t3 + t150 * t4 - t17 * t97 - t18 * t96) * mrSges(6,3) + (t1 * t150 + t104 * t2 - t11 * t96 + t13 * t97) * mrSges(7,2) + (t128 * t106 + t107 * t205 + t141 * t108 / 0.2e1 + (t120 * t205 + t141 * t206) * qJD(4) + (-Ifges(4,6) + t192 / 0.2e1 + t191 / 0.2e1 + t128 * mrSges(4,2) + (Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1) * t104 - (Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1) * t150) * qJD(3)) * t140 - (t7 / 0.2e1 - t8 / 0.2e1) * t150 + (t9 / 0.2e1 + t10 / 0.2e1) * t104 + t27 * t50 + t26 * t55 + t61 * t15 - pkin(3) * t62 + t5 * t63 + (-pkin(8) * t78 - t23 * mrSges(5,3) + t177 * t206 + t53 / 0.2e1 + (-t82 / 0.2e1 - pkin(8) * t109 + pkin(4) * t51 - t70 * mrSges(5,3) + t190 / 0.2e1) * qJD(4)) * t139 + t71 * t64 + t95 * t56 + t130 * t16; (-t106 + t215) * t142 + m(6) * (-pkin(4) * t162 + t151) + m(7) * (-t142 * t27 + t151) + ((mrSges(5,3) * t179 - mrSges(4,2)) * t142 + m(5) * (t179 * t202 - t204) + (m(6) * t130 + m(7) * t61 + t188 + t63 + t64) * t140) * qJD(3) + t221 * (-t104 * t44 + t150 * t45 + t84 * t97 - t85 * t96); -0.2e1 * pkin(3) * t106 + t141 * t107 + t139 * t108 + 0.2e1 * t130 * t56 + 0.2e1 * t27 * t63 + 0.2e1 * t61 * t55 + (t67 + t68) * t97 + (t65 - t66) * t96 + (t59 + t60) * t104 - (t57 - t58) * t150 + (t141 * t120 + (0.2e1 * pkin(4) * t64 - t119) * t139) * qJD(4) + (t27 * t61 + t161) * t209 + (t130 * t172 + t161) * t210 + 0.2e1 * t221 * (t104 * t47 + t150 * t48 + t72 * t97 - t73 * t96); -t22 * mrSges(5,2) + t23 * mrSges(5,1) + t3 * mrSges(6,1) - t4 * mrSges(6,2) + t1 * mrSges(7,3) - t2 * mrSges(7,1) + (t137 * t4 + t186 * t3) * t207 + m(7) * (qJD(6) * t11 + t1 * t126 + t129 * t2) + t30 * t169 - t148 + t29 * t203 + qJD(6) * t74 - Ifges(5,5) * t165 + t126 * t28 + t129 * t31 - t146 * Ifges(5,6); (t137 * t45 + t186 * t44) * t207 + m(7) * (t126 * t45 - t129 * t44 + t185) + t222; m(7) * qJD(6) * t73 - pkin(8) * mrSges(5,1) * t174 + mrSges(5,2) * t171 - Ifges(5,6) * t176 + t129 * t201 + t134 + t90 - t91 + t92 + t93 + (t137 * t207 - mrSges(6,2) + t212) * t48 + (m(7) * t129 - t186 * t207 - mrSges(6,1) - mrSges(7,1)) * t47 + (-t169 * t97 - t203 * t96) * mrSges(6,3) + (qJD(6) * t150 - t126 * t96) * mrSges(7,2); 0.2e1 * t212 * qJD(6); m(6) * t71 + m(7) * t5 - t152; t216 * t178; m(6) * t172 + m(7) * t27 - t215; 0; 0; m(7) * t2 + t31; -m(7) * t44; m(7) * t47 + t201; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
