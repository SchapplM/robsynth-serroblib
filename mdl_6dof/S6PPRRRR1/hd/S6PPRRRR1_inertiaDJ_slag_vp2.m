% Calculate time derivative of joint inertia matrix for
% S6PPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRRRR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_inertiaDJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:59:47
% EndTime: 2019-03-08 18:59:52
% DurationCPUTime: 3.04s
% Computational Cost: add. (4183->341), mult. (11764->535), div. (0->0), fcn. (12391->14), ass. (0->162)
t116 = sin(pkin(7));
t128 = cos(qJ(3));
t175 = t116 * t128;
t115 = sin(pkin(13));
t117 = sin(pkin(6));
t120 = cos(pkin(6));
t124 = sin(qJ(3));
t118 = cos(pkin(13));
t119 = cos(pkin(7));
t174 = t118 * t119;
t218 = (-t115 * t124 + t128 * t174) * t117 + t120 * t175;
t59 = t218 * qJD(3);
t121 = sin(qJ(6));
t125 = cos(qJ(6));
t126 = cos(qJ(5));
t215 = (t121 ^ 2 + t125 ^ 2) * t126;
t101 = -mrSges(7,1) * t125 + mrSges(7,2) * t121;
t220 = t101 - mrSges(6,1);
t170 = qJD(3) * t124;
t156 = t116 * t170;
t122 = sin(qJ(5));
t123 = sin(qJ(4));
t127 = cos(qJ(4));
t176 = t116 * t124;
t86 = t119 * t127 - t123 * t176;
t87 = t119 * t123 + t127 * t176;
t55 = t122 * t87 - t126 * t86;
t169 = qJD(3) * t128;
t155 = t116 * t169;
t77 = qJD(4) * t86 + t127 * t155;
t78 = -qJD(4) * t87 - t123 * t155;
t25 = -qJD(5) * t55 + t122 * t78 + t126 * t77;
t56 = t122 * t86 + t126 * t87;
t47 = -t121 * t56 - t125 * t175;
t12 = qJD(6) * t47 + t121 * t156 + t125 * t25;
t138 = t121 * t175 - t125 * t56;
t13 = qJD(6) * t138 - t121 * t25 + t125 * t156;
t219 = -t13 * t121 + (t121 * t138 - t125 * t47) * qJD(6) + t12 * t125;
t167 = qJD(6) * t125;
t214 = qJD(4) + qJD(5);
t92 = t122 * t123 - t126 * t127;
t72 = t214 * t92;
t93 = t122 * t127 + t123 * t126;
t140 = -t121 * t72 + t93 * t167;
t217 = m(5) * pkin(9) + mrSges(5,3);
t213 = 2 * m(7);
t212 = 2 * pkin(4);
t211 = -2 * mrSges(6,3);
t205 = -pkin(10) - pkin(9);
t157 = qJD(4) * t205;
t100 = t123 * t157;
t150 = t127 * t157;
t106 = t205 * t127;
t163 = t205 * t123;
t83 = -t126 * t106 + t122 * t163;
t42 = qJD(5) * t83 + t122 * t100 - t126 * t150;
t210 = 0.2e1 * t42;
t82 = -t122 * t106 - t126 * t163;
t209 = 0.2e1 * t82;
t208 = m(6) / 0.2e1;
t207 = m(7) / 0.2e1;
t71 = t120 * t176 + (t115 * t128 + t124 * t174) * t117;
t60 = qJD(3) * t71;
t204 = m(5) * t60;
t85 = -t116 * t117 * t118 + t119 * t120;
t43 = -t123 * t71 + t127 * t85;
t44 = t123 * t85 + t127 * t71;
t21 = t122 * t44 - t126 * t43;
t22 = t122 * t43 + t126 * t44;
t29 = -qJD(4) * t44 - t123 * t59;
t30 = qJD(4) * t43 + t127 * t59;
t7 = qJD(5) * t22 + t122 * t30 - t126 * t29;
t203 = t21 * t7;
t201 = t42 * t82;
t26 = qJD(5) * t56 + t122 * t77 - t126 * t78;
t200 = t55 * t26;
t35 = t218 * t60;
t73 = t214 * t93;
t38 = mrSges(6,1) * t73 - mrSges(6,2) * t72;
t96 = (mrSges(5,1) * t123 + mrSges(5,2) * t127) * qJD(4);
t199 = t38 + t96;
t184 = t125 * t72;
t198 = -Ifges(7,5) * t184 + Ifges(7,3) * t73;
t197 = mrSges(7,3) * t125;
t196 = Ifges(7,4) * t121;
t195 = Ifges(7,4) * t125;
t194 = Ifges(7,6) * t121;
t193 = pkin(4) * qJD(5);
t164 = pkin(4) * qJD(4) * t123;
t33 = pkin(5) * t73 + pkin(11) * t72 + t164;
t110 = -pkin(4) * t127 - pkin(3);
t63 = pkin(5) * t92 - pkin(11) * t93 + t110;
t37 = t121 * t63 + t125 * t83;
t41 = -qJD(5) * t82 + t126 * t100 + t122 * t150;
t10 = -qJD(6) * t37 - t121 * t41 + t125 * t33;
t192 = t10 * t121;
t189 = t121 * t93;
t188 = t122 * mrSges(6,1);
t187 = t122 * t21;
t186 = t122 * t55;
t185 = t122 * t82;
t183 = t125 * t93;
t182 = t126 * mrSges(6,2);
t108 = pkin(4) * t122 + pkin(11);
t181 = t108 * t121;
t180 = t108 * t125;
t173 = t121 * t126;
t172 = t122 * t101;
t171 = t125 * t126;
t168 = qJD(6) * t121;
t166 = 0.2e1 * t123;
t102 = -mrSges(5,1) * t127 + mrSges(5,2) * t123;
t76 = mrSges(6,1) * t92 + mrSges(6,2) * t93;
t165 = -mrSges(4,1) + t102 + t76;
t162 = t93 * t168;
t139 = t162 + t184;
t27 = mrSges(7,1) * t140 - mrSges(7,2) * t139;
t158 = m(7) * t42 + t27;
t154 = t108 * t168;
t153 = t108 * t167;
t152 = -t168 / 0.2e1;
t151 = -(2 * Ifges(6,4)) - t194;
t149 = t116 ^ 2 * t124 * t169;
t148 = mrSges(7,3) * t215;
t147 = t26 * t21 + t55 * t7;
t146 = t21 * t42 + t7 * t82;
t145 = t82 * t26 + t42 * t55;
t144 = mrSges(7,1) * t121 + mrSges(7,2) * t125;
t143 = Ifges(7,1) * t125 - t196;
t142 = -Ifges(7,2) * t121 + t195;
t141 = Ifges(7,5) * t121 + Ifges(7,6) * t125;
t15 = -t121 * t218 + t125 * t22;
t14 = -t121 * t22 - t125 * t218;
t36 = -t121 * t83 + t125 * t63;
t103 = Ifges(7,2) * t125 + t196;
t104 = Ifges(7,1) * t121 + t195;
t97 = t142 * qJD(6);
t98 = t143 * qJD(6);
t137 = -t103 * t168 + t104 * t167 + t121 * t98 + t125 * t97;
t6 = -qJD(5) * t21 + t122 * t29 + t126 * t30;
t3 = -qJD(6) * t15 - t121 * t6 + t125 * t60;
t136 = -t3 * t121 + (-t121 * t15 - t125 * t14) * qJD(6);
t2 = qJD(6) * t14 + t121 * t60 + t125 * t6;
t95 = t144 * qJD(6);
t132 = -t6 * mrSges(6,2) + mrSges(7,3) * t136 + t2 * t197 + t21 * t95 + t220 * t7;
t131 = -t25 * mrSges(6,2) + t219 * mrSges(7,3) + t220 * t26 + t55 * t95;
t31 = mrSges(7,1) * t73 + mrSges(7,3) * t139;
t32 = -mrSges(7,2) * t73 - mrSges(7,3) * t140;
t65 = -mrSges(7,2) * t92 - mrSges(7,3) * t189;
t66 = mrSges(7,1) * t92 - mrSges(7,3) * t183;
t9 = qJD(6) * t36 + t121 * t33 + t125 * t41;
t130 = -t66 * t167 - t65 * t168 + t125 * t32 + m(7) * (t125 * t9 - t167 * t36 - t168 * t37 - t192) - t121 * t31;
t111 = Ifges(7,5) * t167;
t18 = -Ifges(7,4) * t139 - Ifges(7,2) * t140 + t73 * Ifges(7,6);
t19 = -Ifges(7,1) * t139 - Ifges(7,4) * t140 + t73 * Ifges(7,5);
t50 = t92 * Ifges(7,6) + t142 * t93;
t51 = t92 * Ifges(7,5) + t143 * t93;
t129 = t121 * t19 / 0.2e1 + t125 * t18 / 0.2e1 + t50 * t152 + t51 * t167 / 0.2e1 + t82 * t95 - Ifges(6,5) * t72 - t97 * t189 / 0.2e1 + t98 * t183 / 0.2e1 + t9 * t197 + t92 * (-Ifges(7,6) * t168 + t111) / 0.2e1 + (-t192 + (-t121 * t37 - t125 * t36) * qJD(6)) * mrSges(7,3) - t41 * mrSges(6,2) + (t141 / 0.2e1 - Ifges(6,6)) * t73 + t220 * t42 - t140 * t103 / 0.2e1 + (-t184 / 0.2e1 + t93 * t152) * t104;
t109 = -pkin(4) * t126 - pkin(5);
t61 = t144 * t93;
t57 = t218 * t156;
t1 = [0.2e1 * m(7) * (t14 * t3 + t15 * t2 + t203) + 0.2e1 * m(6) * (t22 * t6 + t203 - t35) + 0.2e1 * m(5) * (t29 * t43 + t30 * t44 - t35) + 0.2e1 * m(4) * (t59 * t71 - t35); m(7) * (t12 * t15 + t13 * t14 - t138 * t2 + t3 * t47 + t147) + m(5) * (t86 * t29 + t87 * t30 + t78 * t43 + t77 * t44 - t57) - t204 * t175 + (-t60 * t175 + t25 * t22 + t56 * t6 + t147 - t57) * m(6); 0.2e1 * m(7) * (-t12 * t138 + t13 * t47 + t200) + 0.2e1 * m(6) * (t56 * t25 - t149 + t200) + 0.2e1 * m(5) * (t87 * t77 + t86 * t78 - t149); -t59 * mrSges(4,2) + t14 * t31 + t15 * t32 + t2 * t65 + t21 * t27 + t3 * t66 + t7 * t61 - t199 * t218 + t165 * t60 + m(7) * (t10 * t14 + t15 * t9 + t2 * t37 + t3 * t36 + t146) + m(6) * (t110 * t60 - t164 * t218 + t22 * t41 + t6 * t83 + t146) - pkin(3) * t204 + (-t21 * t72 - t22 * t73 - t6 * t92 + t7 * t93) * mrSges(6,3) + t217 * (-t29 * t123 + t30 * t127 + (-t123 * t44 - t127 * t43) * qJD(4)); t12 * t65 + t13 * t66 + t26 * t61 + t55 * t27 + t47 * t31 - t138 * t32 + (-t199 * t128 + (-mrSges(4,2) * t128 + t124 * t165) * qJD(3)) * t116 + m(7) * (t10 * t47 + t12 * t37 + t13 * t36 - t138 * t9 + t145) + m(6) * (t83 * t25 + t41 * t56 + (t110 * t170 - t128 * t164) * t116 + t145) - m(5) * pkin(3) * t156 + (-t25 * t92 + t26 * t93 - t55 * t72 - t56 * t73) * mrSges(6,3) + t217 * (-t78 * t123 + t77 * t127 + (-t123 * t87 - t127 * t86) * qJD(4)); t83 * t73 * t211 - 0.2e1 * pkin(3) * t96 + 0.2e1 * t10 * t66 + 0.2e1 * t110 * t38 + t27 * t209 + 0.2e1 * t36 * t31 + 0.2e1 * t37 * t32 + t61 * t210 + 0.2e1 * t9 * t65 + 0.2e1 * m(6) * (t110 * t164 + t41 * t83 + t201) + (t10 * t36 + t37 * t9 + t201) * t213 - (mrSges(6,3) * t209 - t121 * t50 + t125 * t51) * t72 + (t41 * t211 + ((2 * Ifges(6,2)) + Ifges(7,3)) * t73 - t151 * t72 + t198) * t92 + (mrSges(6,3) * t210 - 0.2e1 * Ifges(6,1) * t72 - t121 * t18 + t125 * t19 + (Ifges(7,5) * t125 + t151) * t73 + (-t121 * t51 - t125 * t50 - t141 * t92) * qJD(6)) * t93 + ((-Ifges(5,4) * t123 + pkin(4) * t76) * t166 + (0.2e1 * Ifges(5,4) * t127 + (Ifges(5,1) - Ifges(5,2)) * t166) * t127) * qJD(4); m(7) * (t109 * t7 - t14 * t153 - t15 * t154 + t180 * t2 - t181 * t3) - t30 * mrSges(5,2) + t29 * mrSges(5,1) + ((t122 * t6 - t126 * t7) * t208 + ((-t14 * t173 + t15 * t171 + t187) * t207 + (t126 * t22 + t187) * t208) * qJD(5)) * t212 + t132; m(7) * (t109 * t26 + t12 * t180 - t13 * t181 + t138 * t154 - t153 * t47) - t77 * mrSges(5,2) + t78 * mrSges(5,1) + ((t122 * t25 - t126 * t26) * t208 + ((-t138 * t171 - t173 * t47 + t186) * t207 + (t126 * t56 + t186) * t208) * qJD(5)) * t212 + t131; t158 * t109 + t130 * t108 + (m(6) * (t122 * t41 - t126 * t42) + (-t122 * t73 + t126 * t72) * mrSges(6,3) + ((t93 * mrSges(6,3) + t61) * t122 + (-t92 * mrSges(6,3) - t121 * t66 + t125 * t65) * t126 + m(7) * (t171 * t37 - t173 * t36 + t185) + m(6) * (t126 * t83 + t185)) * qJD(5)) * pkin(4) + (Ifges(5,5) * t127 - Ifges(5,6) * t123 + pkin(9) * t102) * qJD(4) + t129; 0.2e1 * t109 * t95 + (0.2e1 * t172 + (t215 * t108 + t109 * t122) * t213 - 0.2e1 * t182 - 0.2e1 * t188 + 0.2e1 * t148) * t193 + t137; m(7) * (-pkin(5) * t7 + (t125 * t2 + t136) * pkin(11)) + t132; m(7) * (-pkin(5) * t26 + t219 * pkin(11)) + t131; -pkin(5) * t158 + pkin(11) * t130 + t129; (-pkin(5) + t109) * t95 + (t172 + m(7) * (-pkin(5) * t122 + t215 * pkin(11)) - t182 - t188 + t148) * t193 + t137; -0.2e1 * pkin(5) * t95 + t137; mrSges(7,1) * t3 - mrSges(7,2) * t2; mrSges(7,1) * t13 - mrSges(7,2) * t12; mrSges(7,1) * t10 - mrSges(7,2) * t9 - Ifges(7,5) * t162 - Ifges(7,6) * t140 + t198; t111 - t144 * t126 * t193 + (t101 * t108 - t194) * qJD(6); t111 + (pkin(11) * t101 - t194) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
