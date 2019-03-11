% Calculate time derivative of joint inertia matrix for
% S6PRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:26:48
% EndTime: 2019-03-08 20:26:55
% DurationCPUTime: 2.84s
% Computational Cost: add. (3324->397), mult. (8748->603), div. (0->0), fcn. (8246->12), ass. (0->173)
t113 = sin(qJ(6));
t114 = sin(qJ(5));
t117 = cos(qJ(6));
t118 = cos(qJ(5));
t127 = t113 * t114 - t117 * t118;
t189 = qJD(5) + qJD(6);
t51 = t189 * t127;
t84 = t113 * t118 + t114 * t117;
t52 = t189 * t84;
t21 = mrSges(7,1) * t52 - mrSges(7,2) * t51;
t132 = mrSges(6,1) * t114 + mrSges(6,2) * t118;
t85 = t132 * qJD(5);
t193 = t21 + t85;
t115 = sin(qJ(4));
t119 = cos(qJ(4));
t159 = cos(pkin(6));
t110 = sin(pkin(12));
t111 = sin(pkin(6));
t112 = cos(pkin(12));
t116 = sin(qJ(2));
t120 = cos(qJ(2));
t68 = (t110 * t120 + t112 * t116) * t111;
t124 = -t68 * t115 + t159 * t119;
t192 = t124 * qJD(4);
t149 = qJD(4) * t119;
t29 = -t52 * t115 - t127 * t149;
t73 = t127 * t115;
t30 = -t84 * t149 + t189 * t73;
t11 = -t30 * mrSges(7,1) + t29 * mrSges(7,2);
t146 = qJD(5) * t118;
t122 = t114 * t149 + t115 * t146;
t137 = t118 * t149;
t147 = qJD(5) * t115;
t138 = t114 * t147;
t123 = t137 - t138;
t45 = t122 * mrSges(6,1) + t123 * mrSges(6,2);
t171 = t11 + t45;
t103 = -pkin(2) * t112 - pkin(3);
t80 = -pkin(4) * t119 - pkin(9) * t115 + t103;
t102 = pkin(2) * t110 + pkin(8);
t154 = t118 * t119;
t91 = t102 * t154;
t50 = t114 * t80 + t91;
t191 = qJD(5) * t50;
t151 = qJD(4) * t115;
t190 = -Ifges(6,5) * t137 - Ifges(6,3) * t151;
t153 = t114 ^ 2 + t118 ^ 2;
t104 = -pkin(5) * t118 - pkin(4);
t95 = -mrSges(6,1) * t118 + mrSges(6,2) * t114;
t168 = t95 - mrSges(5,1);
t55 = mrSges(7,1) * t127 + mrSges(7,2) * t84;
t188 = m(7) * t104 + t168 + t55;
t187 = 2 * m(6);
t186 = 0.2e1 * m(7);
t185 = 0.2e1 * t102;
t184 = m(5) / 0.2e1;
t183 = m(6) / 0.2e1;
t182 = -m(7) / 0.2e1;
t181 = m(6) * pkin(4);
t180 = m(7) * pkin(5);
t179 = -t127 / 0.2e1;
t178 = t84 / 0.2e1;
t167 = Ifges(6,4) * t114;
t96 = Ifges(6,2) * t118 + t167;
t177 = -t96 / 0.2e1;
t176 = -pkin(10) - pkin(9);
t175 = -t114 / 0.2e1;
t174 = pkin(4) * t115;
t173 = pkin(9) * t119;
t54 = t159 * t115 + t68 * t119;
t158 = qJD(4) * t54;
t128 = t110 * t116 - t112 * t120;
t152 = qJD(2) * t111;
t66 = t128 * t152;
t34 = -t66 * t115 + t158;
t14 = t124 * t34;
t65 = qJD(2) * t68;
t67 = t128 * t111;
t172 = t65 * t67;
t72 = t84 * t115;
t43 = mrSges(7,1) * t72 - mrSges(7,2) * t73;
t79 = t132 * t115;
t170 = t43 + t79;
t157 = t102 * t114;
t94 = (-t173 + t174) * qJD(4);
t169 = t118 * t94 + t151 * t157;
t166 = Ifges(6,4) * t118;
t165 = Ifges(6,5) * t114;
t164 = Ifges(6,6) * t114;
t163 = Ifges(6,6) * t118;
t162 = Ifges(6,6) * t119;
t161 = t34 * t115;
t35 = -t66 * t119 + t192;
t160 = t35 * t119;
t156 = t114 * t115;
t155 = t115 * t118;
t150 = qJD(4) * t118;
t148 = qJD(5) * t114;
t145 = qJD(6) * t113;
t144 = qJD(6) * t117;
t142 = -Ifges(7,5) * t29 - Ifges(7,6) * t30 - Ifges(7,3) * t151;
t141 = pkin(5) * t148;
t36 = -t114 * t54 + t118 * t67;
t37 = t114 * t67 + t118 * t54;
t12 = -t113 * t37 + t117 * t36;
t7 = -t37 * qJD(5) - t114 * t35 + t118 * t65;
t8 = t36 * qJD(5) + t114 * t65 + t118 * t35;
t2 = t12 * qJD(6) + t113 * t7 + t117 * t8;
t13 = t113 * t36 + t117 * t37;
t3 = -t13 * qJD(6) - t113 * t8 + t117 * t7;
t140 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t139 = qJD(5) * t176;
t136 = t119 * t148;
t135 = (2 * Ifges(5,4)) + t164;
t27 = t80 * t146 + t114 * t94 + (-t115 * t150 - t136) * t102;
t75 = t118 * t80;
t49 = -t119 * t157 + t75;
t134 = -qJD(5) * t49 + t27;
t133 = -t114 * t7 + t118 * t8;
t131 = Ifges(6,1) * t118 - t167;
t97 = Ifges(6,1) * t114 + t166;
t130 = -Ifges(6,2) * t114 + t166;
t40 = -pkin(10) * t155 + t75 + (-pkin(5) - t157) * t119;
t44 = -pkin(10) * t156 + t50;
t15 = -t113 * t44 + t117 * t40;
t16 = t113 * t40 + t117 * t44;
t98 = t176 * t114;
t99 = t176 * t118;
t59 = t113 * t99 + t117 * t98;
t60 = t113 * t98 - t117 * t99;
t92 = t114 * t139;
t93 = t118 * t139;
t32 = t59 * qJD(6) + t113 * t93 + t117 * t92;
t33 = -t60 * qJD(6) - t113 * t92 + t117 * t93;
t47 = Ifges(7,6) * t52;
t48 = Ifges(7,5) * t51;
t129 = t33 * mrSges(7,1) - t32 * mrSges(7,2) - t47 - t48;
t17 = (pkin(5) * t115 - pkin(10) * t154) * qJD(4) + (-t91 + (pkin(10) * t115 - t80) * t114) * qJD(5) + t169;
t18 = -t122 * pkin(10) + t27;
t5 = t15 * qJD(6) + t113 * t17 + t117 * t18;
t6 = -t16 * qJD(6) - t113 * t18 + t117 * t17;
t126 = t6 * mrSges(7,1) - t5 * mrSges(7,2) - t142;
t125 = -t124 * t149 + t161;
t107 = Ifges(6,5) * t146;
t90 = -mrSges(6,1) * t119 - mrSges(6,3) * t155;
t89 = mrSges(6,2) * t119 - mrSges(6,3) * t156;
t88 = t131 * qJD(5);
t87 = t130 * qJD(5);
t86 = (t115 * mrSges(5,1) + t119 * mrSges(5,2)) * qJD(4);
t81 = (-mrSges(7,1) * t113 - mrSges(7,2) * t117) * qJD(6) * pkin(5);
t76 = (pkin(5) * t114 + t102) * t115;
t71 = -Ifges(6,5) * t119 + t131 * t115;
t70 = t130 * t115 - t162;
t64 = -mrSges(6,2) * t151 - t122 * mrSges(6,3);
t63 = mrSges(6,1) * t151 - t123 * mrSges(6,3);
t62 = -mrSges(7,1) * t119 + mrSges(7,3) * t73;
t61 = mrSges(7,2) * t119 - mrSges(7,3) * t72;
t58 = t122 * pkin(5) + t102 * t149;
t57 = Ifges(7,1) * t84 - Ifges(7,4) * t127;
t56 = Ifges(7,4) * t84 - Ifges(7,2) * t127;
t46 = t124 * t151;
t42 = -t97 * t147 + (Ifges(6,5) * t115 + t131 * t119) * qJD(4);
t41 = -t96 * t147 + (Ifges(6,6) * t115 + t130 * t119) * qJD(4);
t39 = -Ifges(7,1) * t73 - Ifges(7,4) * t72 - Ifges(7,5) * t119;
t38 = -Ifges(7,4) * t73 - Ifges(7,2) * t72 - Ifges(7,6) * t119;
t28 = t169 - t191;
t23 = -Ifges(7,1) * t51 - Ifges(7,4) * t52;
t22 = -Ifges(7,4) * t51 - Ifges(7,2) * t52;
t20 = -mrSges(7,2) * t151 + mrSges(7,3) * t30;
t19 = mrSges(7,1) * t151 - mrSges(7,3) * t29;
t10 = Ifges(7,1) * t29 + Ifges(7,4) * t30 + Ifges(7,5) * t151;
t9 = Ifges(7,4) * t29 + Ifges(7,2) * t30 + Ifges(7,6) * t151;
t1 = [0.2e1 * m(7) * (t12 * t3 + t13 * t2 - t14) + 0.2e1 * m(6) * (t36 * t7 + t37 * t8 - t14) + 0.2e1 * m(5) * (t35 * t54 - t14 + t172) + 0.2e1 * m(4) * (-t66 * t68 + t172); t66 * mrSges(4,2) + t12 * t19 + t13 * t20 + t2 * t61 + t3 * t62 + t36 * t63 + t37 * t64 + t67 * t86 + t7 * t90 + t8 * t89 + (-mrSges(5,1) * t119 + mrSges(5,2) * t115 - mrSges(4,1)) * t65 - t171 * t124 + t170 * t34 + (-mrSges(3,1) * t116 - mrSges(3,2) * t120) * t152 + (t161 + t160 + (-t115 * t54 - t119 * t124) * qJD(4)) * mrSges(5,3) + m(7) * (t12 * t6 - t124 * t58 + t13 * t5 + t15 * t3 + t16 * t2 + t34 * t76) + m(6) * (t27 * t37 + t28 * t36 + t49 * t7 + t50 * t8) + m(5) * t103 * t65 + (t125 * t183 + (-t54 * t151 + t125 + t160) * t184) * t185 + m(4) * (-t110 * t66 - t112 * t65) * pkin(2); -t73 * t10 + 0.2e1 * t103 * t86 + 0.2e1 * t76 * t11 + 0.2e1 * t15 * t19 + 0.2e1 * t16 * t20 + 0.2e1 * t27 * t89 + 0.2e1 * t28 * t90 + t29 * t39 + t30 * t38 + 0.2e1 * t58 * t43 + 0.2e1 * t49 * t63 + 0.2e1 * t5 * t61 + 0.2e1 * t50 * t64 + 0.2e1 * t6 * t62 - t72 * t9 + (t15 * t6 + t16 * t5 + t58 * t76) * t186 + (t27 * t50 + t28 * t49) * t187 + ((-t114 * t70 + t118 * t71 + t135 * t119 + t79 * t185) * qJD(4) + t142 + t190) * t119 + (t45 * t185 - t114 * t41 + t118 * t42 + (-t118 * t70 - t114 * t71 - t119 * (-t163 - t165)) * qJD(5) + (-Ifges(7,5) * t73 - Ifges(7,6) * t72 + (Ifges(6,5) * t118 - t135) * t115 + (t102 ^ 2 * t187 + (2 * Ifges(5,1)) - (2 * Ifges(5,2)) - Ifges(6,3) - Ifges(7,3)) * t119) * qJD(4)) * t115; m(7) * (t12 * t30 + t13 * t29 - t2 * t73 - t3 * t72 - t46) - m(6) * t46 + 0.2e1 * (t34 * t182 + (-qJD(4) * t114 * t36 + t37 * t150 - t34) * t183 + (-t34 + t158) * t184) * t119 + 0.2e1 * ((-t36 * t146 - t37 * t148 + t133) * t183 + (t35 - t192) * t184) * t115; -t72 * t19 - t73 * t20 + t29 * t61 + t30 * t62 - t171 * t119 + t170 * t151 + m(7) * (-t119 * t58 + t15 * t30 + t76 * t151 + t16 * t29 - t5 * t73 - t6 * t72) + m(6) * (t115 ^ 2 - t119 ^ 2) * qJD(4) * t102 + (m(6) * (t115 * t27 - t49 * t147 + t50 * t149) + t89 * t149 + t115 * t64 - t90 * t147) * t118 + (m(6) * (-t115 * t28 - t50 * t147 - t49 * t149) - t89 * t147 - t90 * t149 - t115 * t63) * t114; (-t29 * t73 - t30 * t72) * t186 + 0.4e1 * ((-0.1e1 + t153) * t183 + t182) * t115 * t149; -t35 * mrSges(5,2) - t193 * t124 + m(7) * (t12 * t33 - t124 * t141 + t13 * t32 + t2 * t60 + t3 * t59) + (t12 * t51 - t127 * t2 - t13 * t52 - t3 * t84) * mrSges(7,3) + (-t181 + t188) * t34 + (m(6) * pkin(9) + mrSges(6,3)) * ((-t114 * t37 - t118 * t36) * qJD(5) + t133); m(7) * (t104 * t58 + t15 * t33 + t16 * t32 + t5 * t60 + t59 * t6) - pkin(4) * t45 - t51 * t39 / 0.2e1 - t52 * t38 / 0.2e1 + t30 * t56 / 0.2e1 + t29 * t57 / 0.2e1 + t58 * t55 + t59 * t19 + t60 * t20 + t32 * t61 + t33 * t62 - t72 * t22 / 0.2e1 - t73 * t23 / 0.2e1 + t76 * t21 + t9 * t179 + t10 * t178 + t104 * t11 + (t48 / 0.2e1 + t47 / 0.2e1 - t107 / 0.2e1 + (Ifges(5,5) + (t168 - t181) * t102) * qJD(4)) * t119 + (t149 * t177 - t28 * mrSges(6,3) + t42 / 0.2e1 + (-t50 * mrSges(6,3) - t70 / 0.2e1 + t162 / 0.2e1 + (m(7) * t76 + t43) * pkin(5)) * qJD(5) + (m(6) * (-t28 - t191) - qJD(5) * t89 - t63) * pkin(9)) * t114 + (-t127 * t5 + t15 * t51 - t16 * t52 - t6 * t84) * mrSges(7,3) + (t97 * t149 / 0.2e1 + qJD(5) * t71 / 0.2e1 + t41 / 0.2e1 + t134 * mrSges(6,3) + (m(6) * t134 - qJD(5) * t90 + t64) * pkin(9)) * t118 + (t102 * t85 + t87 * t175 + t118 * t88 / 0.2e1 + (t118 * t177 + t175 * t97) * qJD(5) + (t102 * mrSges(5,2) + Ifges(7,5) * t178 + Ifges(7,6) * t179 - Ifges(5,6) + t165 / 0.2e1 + t163 / 0.2e1) * qJD(4)) * t115; m(7) * (-pkin(5) * t136 + t29 * t60 + t30 * t59 - t32 * t73 - t33 * t72) + (-t127 * t29 - t30 * t84 - t51 * t72 + t52 * t73) * mrSges(7,3) + (m(6) * (t153 * t173 - t174) + t188 * t115) * qJD(4) + ((mrSges(6,3) * t153 - mrSges(5,2)) * qJD(4) - t193) * t119; -t51 * t57 + t84 * t23 - t52 * t56 - t127 * t22 + (t104 * t141 + t32 * t60 + t33 * t59) * t186 + 0.2e1 * t55 * t141 + 0.2e1 * t104 * t21 - 0.2e1 * pkin(4) * t85 - t96 * t148 + t114 * t88 + (qJD(5) * t97 + t87) * t118 + 0.2e1 * (-t127 * t32 - t33 * t84 + t51 * t59 - t52 * t60) * mrSges(7,3); t7 * mrSges(6,1) - t8 * mrSges(6,2) + (t113 * t2 + t117 * t3 + (-t113 * t12 + t117 * t13) * qJD(6)) * t180 + t140; -Ifges(6,5) * t138 + t28 * mrSges(6,1) - t27 * mrSges(6,2) - t122 * Ifges(6,6) + (m(7) * (t113 * t5 + t117 * t6 + t144 * t16 - t145 * t15) + t61 * t144 + t113 * t20 - t62 * t145 + t117 * t19) * pkin(5) + t126 - t190; (t113 * t29 + t117 * t30 + (t113 * t72 - t117 * t73) * qJD(6)) * t180 - t171; t107 + (pkin(9) * t95 - t164) * qJD(5) + (m(7) * (t113 * t32 + t117 * t33 + (-t113 * t59 + t117 * t60) * qJD(6)) + (-t113 * t52 + t117 * t51 + (t113 * t84 - t117 * t127) * qJD(6)) * mrSges(7,3)) * pkin(5) + t129; 0.2e1 * t81; t140; t126; -t11; t129; t81; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
