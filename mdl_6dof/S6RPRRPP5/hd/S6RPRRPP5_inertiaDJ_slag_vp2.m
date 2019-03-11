% Calculate time derivative of joint inertia matrix for
% S6RPRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-03-09 04:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:42:44
% EndTime: 2019-03-09 04:42:51
% DurationCPUTime: 2.99s
% Computational Cost: add. (2644->380), mult. (6113->510), div. (0->0), fcn. (5437->6), ass. (0->162)
t199 = Ifges(6,4) + Ifges(5,5);
t197 = Ifges(6,2) + Ifges(5,3);
t117 = sin(qJ(4));
t119 = cos(qJ(4));
t152 = qJD(4) * t119;
t115 = sin(pkin(9));
t116 = cos(pkin(9));
t118 = sin(qJ(3));
t120 = cos(qJ(3));
t84 = t115 * t120 + t116 * t118;
t143 = t84 * t152;
t83 = t115 * t118 - t116 * t120;
t77 = t83 * qJD(3);
t124 = -t117 * t77 + t143;
t153 = qJD(4) * t117;
t144 = t84 * t153;
t159 = t119 * t77;
t123 = t144 + t159;
t156 = qJ(5) * t119;
t185 = pkin(4) + pkin(5);
t196 = -t117 * t185 + t156;
t78 = t84 * qJD(3);
t195 = qJ(5) * t78 + qJD(5) * t83;
t176 = pkin(7) + qJ(2);
t95 = t176 * t115;
t96 = t176 * t116;
t194 = -t118 * t96 - t120 * t95;
t193 = Ifges(6,6) * t153 + t152 * t199;
t60 = -t118 * t95 + t120 * t96;
t41 = t84 * qJD(2) + t60 * qJD(3);
t192 = 2 * m(7);
t191 = -2 * mrSges(4,3);
t190 = -2 * mrSges(7,3);
t189 = -2 * Ifges(4,4);
t188 = -0.2e1 * t194;
t186 = m(6) + m(7);
t184 = pkin(8) * t117;
t183 = pkin(8) * t119;
t182 = t41 * mrSges(5,1);
t181 = t41 * mrSges(5,2);
t180 = t41 * t194;
t179 = t83 * Ifges(7,5);
t178 = -mrSges(6,2) + mrSges(7,3);
t177 = -Ifges(5,6) - Ifges(7,6);
t175 = pkin(8) - qJ(6);
t27 = mrSges(5,1) * t78 + mrSges(5,3) * t123;
t28 = -t78 * mrSges(6,1) - mrSges(6,2) * t123;
t174 = t28 - t27;
t30 = -mrSges(5,2) * t78 - mrSges(5,3) * t124;
t31 = -mrSges(6,2) * t124 + mrSges(6,3) * t78;
t173 = t30 + t31;
t161 = t117 * t84;
t51 = -mrSges(5,2) * t83 - mrSges(5,3) * t161;
t55 = -mrSges(6,2) * t161 + mrSges(6,3) * t83;
t172 = -t51 - t55;
t158 = t119 * t84;
t53 = mrSges(5,1) * t83 - mrSges(5,3) * t158;
t54 = -mrSges(6,1) * t83 + mrSges(6,2) * t158;
t171 = -t53 + t54;
t142 = -pkin(2) * t116 - pkin(1);
t49 = pkin(3) * t83 - pkin(8) * t84 + t142;
t24 = t117 * t49 + t119 * t60;
t170 = Ifges(5,4) * t117;
t169 = Ifges(5,4) * t119;
t168 = Ifges(7,4) * t117;
t167 = Ifges(7,4) * t119;
t166 = Ifges(6,5) * t117;
t165 = Ifges(6,5) * t119;
t164 = qJ(5) * t77;
t163 = t117 * mrSges(7,1);
t157 = qJ(5) * t117;
t155 = qJD(4) * t84;
t151 = t117 * qJD(5);
t150 = -Ifges(7,5) + t199;
t132 = Ifges(7,1) * t119 + t168;
t35 = t132 * t84 - t179;
t133 = Ifges(6,1) * t119 + t166;
t36 = Ifges(6,4) * t83 + t133 * t84;
t134 = Ifges(5,1) * t119 - t170;
t37 = Ifges(5,5) * t83 + t134 * t84;
t149 = -t35 - t36 - t37;
t40 = -qJD(2) * t83 + qJD(3) * t194;
t148 = t117 * t40 + t152 * t60 + t153 * t49;
t48 = pkin(3) * t78 + pkin(8) * t77;
t147 = t117 * t48 + t119 * t40 + t152 * t49;
t17 = qJ(5) * t83 + t24;
t146 = -Ifges(7,6) / 0.2e1 - Ifges(5,6) / 0.2e1;
t145 = qJ(5) * t155;
t141 = mrSges(4,1) * t78 - mrSges(4,2) * t77;
t100 = t175 * t119;
t57 = t117 * t60;
t23 = t119 * t49 - t57;
t139 = -pkin(4) * t153 + t151;
t108 = mrSges(7,2) * t152;
t86 = -mrSges(7,1) * t153 + t108;
t101 = -Ifges(6,3) * t119 + t166;
t102 = -Ifges(7,2) * t119 + t168;
t103 = Ifges(5,2) * t119 + t170;
t138 = t101 / 0.2e1 + t102 / 0.2e1 - t103 / 0.2e1;
t104 = Ifges(7,1) * t117 - t167;
t105 = Ifges(6,1) * t117 - t165;
t106 = Ifges(5,1) * t117 + t169;
t137 = t104 / 0.2e1 + t105 / 0.2e1 + t106 / 0.2e1;
t136 = mrSges(5,1) * t117 + mrSges(5,2) * t119;
t98 = -t119 * mrSges(6,1) - t117 * mrSges(6,3);
t135 = mrSges(6,1) * t117 - mrSges(6,3) * t119;
t131 = -Ifges(5,2) * t117 + t169;
t130 = Ifges(7,2) * t117 + t167;
t129 = Ifges(6,3) * t117 + t165;
t128 = -pkin(4) * t119 - t157;
t127 = qJ(6) * t77 - qJD(6) * t84;
t26 = -t78 * mrSges(7,1) + mrSges(7,3) * t123;
t7 = t119 * t48 - t148;
t126 = Ifges(6,6) * t124 - t159 * t199 + t197 * t78;
t125 = (m(6) * pkin(8) - t178) * t119;
t6 = -t153 * t60 + t147;
t32 = Ifges(6,6) * t83 + t129 * t84;
t33 = -t83 * Ifges(7,6) + t130 * t84;
t34 = t83 * Ifges(5,6) + t131 * t84;
t122 = t177 * t83 + t32 + t33 - t34;
t21 = -mrSges(7,1) * t124 - mrSges(7,2) * t123;
t99 = mrSges(7,1) * t119 + mrSges(7,2) * t117;
t97 = t175 * t117;
t94 = -pkin(3) + t128;
t93 = t134 * qJD(4);
t92 = t133 * qJD(4);
t91 = t132 * qJD(4);
t90 = t131 * qJD(4);
t89 = t130 * qJD(4);
t88 = t129 * qJD(4);
t87 = t136 * qJD(4);
t85 = t135 * qJD(4);
t81 = t119 * t185 + pkin(3) + t157;
t76 = qJD(4) * t100 - qJD(6) * t117;
t75 = qJ(5) * t152 + t139;
t74 = -qJD(6) * t119 - t153 * t175;
t61 = (-pkin(5) * t117 + t156) * qJD(4) + t139;
t52 = -mrSges(7,1) * t83 - mrSges(7,3) * t158;
t50 = mrSges(7,2) * t83 + mrSges(7,3) * t161;
t47 = (mrSges(7,2) * t119 - t163) * t84;
t46 = t135 * t84;
t29 = mrSges(7,2) * t78 + mrSges(7,3) * t124;
t25 = (pkin(4) * t117 - t156) * t84 - t194;
t22 = mrSges(5,1) * t124 - mrSges(5,2) * t123;
t20 = mrSges(6,1) * t124 + mrSges(6,3) * t123;
t19 = t196 * t84 + t194;
t18 = -pkin(4) * t83 - t23;
t16 = -Ifges(5,1) * t123 - Ifges(5,4) * t124 + t78 * Ifges(5,5);
t15 = -Ifges(6,1) * t123 + t78 * Ifges(6,4) + Ifges(6,5) * t124;
t14 = -Ifges(7,1) * t123 + Ifges(7,4) * t124 - t78 * Ifges(7,5);
t13 = -Ifges(5,4) * t123 - Ifges(5,2) * t124 + t78 * Ifges(5,6);
t12 = -Ifges(7,4) * t123 + Ifges(7,2) * t124 - t78 * Ifges(7,6);
t11 = -Ifges(6,5) * t123 + t78 * Ifges(6,6) + Ifges(6,3) * t124;
t10 = qJ(6) * t161 + t17;
t9 = t57 - t185 * t83 + (-qJ(6) * t84 - t49) * t119;
t8 = (-pkin(4) * t77 + t145) * t117 + (t164 + (pkin(4) * qJD(4) - qJD(5)) * t84) * t119 + t41;
t5 = -pkin(4) * t78 - t7;
t4 = (t185 * t77 - t145) * t117 + (-t164 + (-qJD(4) * t185 + qJD(5)) * t84) * t119 - t41;
t3 = t6 + t195;
t2 = qJ(6) * t143 + (-qJD(4) * t60 - t127) * t117 + t147 + t195;
t1 = qJ(6) * t144 - t185 * t78 + (t127 - t48) * t119 + t148;
t38 = [0.2e1 * t19 * t21 + 0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * (t115 ^ 2 + t116 ^ 2) * qJD(2) + 0.2e1 * m(5) * (t23 * t7 + t24 * t6 - t180) + 0.2e1 * m(4) * (t40 * t60 - t180) + 0.2e1 * m(6) * (t17 * t3 + t18 * t5 + t25 * t8) + (t1 * t9 + t10 * t2 + t19 * t4) * t192 + t22 * t188 - (mrSges(4,3) * t188 + t83 * t189 + (-t149 - t179) * t119 + t122 * t117) * t77 + 0.2e1 * t142 * t141 + (t40 * t191 + ((2 * Ifges(4,2)) + (2 * Ifges(7,3)) + t197) * t78 + t126) * t83 + (0.2e1 * t41 * mrSges(4,3) - 0.2e1 * Ifges(4,1) * t77 + (t14 + t15 + t16 + 0.2e1 * t181) * t119 + (t11 + t12 - t13 + 0.2e1 * t182) * t117 + (t189 + t150 * t119 + (Ifges(6,6) + t177) * t117) * t78 + (t122 * t119 + (-t150 * t83 + t149) * t117) * qJD(4)) * t84 + 0.2e1 * t25 * t20 + 0.2e1 * t9 * t26 + 0.2e1 * t23 * t27 + 0.2e1 * t18 * t28 + 0.2e1 * t10 * t29 + 0.2e1 * t24 * t30 + 0.2e1 * t17 * t31 + 0.2e1 * t8 * t46 + 0.2e1 * t4 * t47 + 0.2e1 * t2 * t50 + 0.2e1 * t6 * t51 + 0.2e1 * t1 * t52 + 0.2e1 * t7 * t53 + 0.2e1 * t5 * t54 + 0.2e1 * t3 * t55 + t60 * t78 * t191; (-t26 - t174) * t119 + (t29 + t173) * t117 + ((t50 - t172) * t119 + (t52 + t171) * t117) * qJD(4) + m(7) * (-t1 * t119 + t117 * t2 + (t10 * t119 + t117 * t9) * qJD(4)) + m(6) * (t117 * t3 - t119 * t5 + (t117 * t18 + t119 * t17) * qJD(4)) + m(5) * (t117 * t6 + t119 * t7 + (-t117 * t23 + t119 * t24) * qJD(4)) + t141; 0; -pkin(3) * t22 + m(7) * (t1 * t97 + t10 * t74 + t100 * t2 + t19 * t61 + t4 * t81 + t76 * t9) + m(5) * (-pkin(3) * t41 + t6 * t183 - t7 * t184) + m(6) * (t3 * t183 + t5 * t184 - t25 * t75 + t8 * t94) + ((t35 / 0.2e1 + t36 / 0.2e1 + t37 / 0.2e1 - t179 / 0.2e1 - t9 * mrSges(7,3) + t18 * mrSges(6,2) - t23 * mrSges(5,3) + t138 * t84) * t119 + (t33 / 0.2e1 - t34 / 0.2e1 + t32 / 0.2e1 + t10 * mrSges(7,3) - t17 * mrSges(6,2) - t24 * mrSges(5,3) + t146 * t83 - t137 * t84) * t117 + (t171 * t119 + t172 * t117 + m(5) * (-t117 * t24 - t119 * t23) + m(6) * (-t117 * t17 + t119 * t18)) * pkin(8)) * qJD(4) + (t5 * mrSges(6,2) - t7 * mrSges(5,3) - t1 * mrSges(7,3) + t14 / 0.2e1 + t15 / 0.2e1 + t16 / 0.2e1 + t181 + t174 * pkin(8) + (Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1 - Ifges(7,5) / 0.2e1) * t78) * t117 + (t3 * mrSges(6,2) + t6 * mrSges(5,3) - t2 * mrSges(7,3) - t11 / 0.2e1 - t12 / 0.2e1 + t13 / 0.2e1 - t182 + t173 * pkin(8) + (-Ifges(6,6) / 0.2e1 - t146) * t78) * t119 - (t117 * t138 + t119 * t137 + Ifges(4,5)) * t77 + t193 * t83 / 0.2e1 - t194 * t87 + ((t91 / 0.2e1 + t92 / 0.2e1 + t93 / 0.2e1) * t119 + (t88 / 0.2e1 + t89 / 0.2e1 - t90 / 0.2e1) * t117) * t84 - t40 * mrSges(4,2) - t41 * mrSges(4,1) + t61 * t47 + t74 * t50 - t75 * t46 + t76 * t52 - Ifges(4,6) * t78 + t81 * t21 + t25 * t85 + t19 * t86 + t94 * t20 + t97 * t26 + t8 * t98 + t4 * t99 + t100 * t29; m(7) * (t117 * t74 - t119 * t76 + (t100 * t119 + t117 * t97) * qJD(4)); 0.2e1 * t94 * t85 + 0.2e1 * t61 * t99 + 0.2e1 * t81 * t86 + (t100 * t74 + t61 * t81 + t76 * t97) * t192 - 0.2e1 * pkin(3) * t87 + 0.2e1 * (-m(6) * t94 - t98) * t75 + (t74 * t190 - t88 - t89 + t90) * t119 + (t76 * t190 + t91 + t92 + t93) * t117 + ((t190 * t97 + t104 + t105 + t106) * t119 + (0.2e1 * mrSges(7,3) * t100 + t101 + t102 - t103) * t117) * qJD(4); -t5 * mrSges(6,1) - t6 * mrSges(5,2) + t7 * mrSges(5,1) - t1 * mrSges(7,1) + t2 * mrSges(7,2) + t3 * mrSges(6,3) + t126 - (-Ifges(7,5) * t119 + t177 * t117) * t77 + (t50 + t55) * qJD(5) + (t29 + t31) * qJ(5) + m(6) * (-pkin(4) * t5 + qJ(5) * t3 + qJD(5) * t17) + m(7) * (qJ(5) * t2 + qJD(5) * t10 - t1 * t185) + (-t117 * t150 + t119 * t177) * t155 - pkin(4) * t28 + Ifges(7,3) * t78 - t185 * t26; t108 + m(7) * t151 + m(6) * t75 + (m(7) * t196 - t135 - t136 - t163) * qJD(4); m(7) * (qJ(5) * t74 + qJD(5) * t100 - t185 * t76) - t76 * mrSges(7,1) + t74 * mrSges(7,2) + qJD(5) * t125 + ((-mrSges(6,2) * pkin(4) + mrSges(7,3) * t185 - Ifges(7,5)) * t119 + (qJ(5) * t178 + t177) * t117 + (m(6) * t128 - t119 * mrSges(5,1) + t117 * mrSges(5,2) + t98) * pkin(8)) * qJD(4) + t193; 0.2e1 * (qJ(5) * t186 + mrSges(7,2) + mrSges(6,3)) * qJD(5); m(6) * t5 + m(7) * t1 + t26 + t28; t186 * t153; m(7) * t76 + qJD(4) * t125; 0; 0; m(7) * t4 + t21; 0; m(7) * t61 + t86; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t38(1) t38(2) t38(4) t38(7) t38(11) t38(16); t38(2) t38(3) t38(5) t38(8) t38(12) t38(17); t38(4) t38(5) t38(6) t38(9) t38(13) t38(18); t38(7) t38(8) t38(9) t38(10) t38(14) t38(19); t38(11) t38(12) t38(13) t38(14) t38(15) t38(20); t38(16) t38(17) t38(18) t38(19) t38(20) t38(21);];
Mq  = res;
