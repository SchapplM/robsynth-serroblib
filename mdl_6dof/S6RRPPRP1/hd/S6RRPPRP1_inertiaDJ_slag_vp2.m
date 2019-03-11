% Calculate time derivative of joint inertia matrix for
% S6RRPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-03-09 08:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:25:28
% EndTime: 2019-03-09 08:25:34
% DurationCPUTime: 2.78s
% Computational Cost: add. (4219->353), mult. (9568->513), div. (0->0), fcn. (9094->8), ass. (0->143)
t197 = Ifges(6,5) + Ifges(7,4);
t196 = -Ifges(6,6) + Ifges(7,6);
t195 = Ifges(7,2) + Ifges(6,3);
t194 = m(6) + m(7);
t133 = sin(pkin(10));
t135 = cos(pkin(10));
t137 = sin(qJ(5));
t139 = cos(qJ(5));
t119 = t133 * t137 - t139 * t135;
t115 = t119 * qJD(5);
t121 = t133 * t139 + t135 * t137;
t116 = t121 * qJD(5);
t63 = pkin(5) * t116 + qJ(6) * t115 - qJD(6) * t121;
t78 = t116 * mrSges(7,1) + t115 * mrSges(7,3);
t79 = t116 * mrSges(6,1) - t115 * mrSges(6,2);
t193 = -m(7) * t63 - t78 - t79;
t192 = -t197 * t115 + t196 * t116;
t191 = m(7) * qJ(6) + mrSges(7,3);
t134 = sin(pkin(9));
t136 = cos(pkin(9));
t138 = sin(qJ(2));
t140 = cos(qJ(2));
t120 = t134 * t140 + t136 * t138;
t113 = t120 * qJD(2);
t118 = t134 * t138 - t136 * t140;
t114 = t118 * qJD(2);
t166 = t114 * t135;
t155 = pkin(2) * qJD(2) * t138;
t55 = pkin(3) * t113 + qJ(4) * t114 - qJD(4) * t120 + t155;
t178 = -qJ(3) - pkin(7);
t152 = qJD(2) * t178;
t112 = qJD(3) * t140 + t138 * t152;
t142 = -t138 * qJD(3) + t140 * t152;
t68 = t136 * t112 + t134 * t142;
t26 = -t133 * t68 + t135 * t55;
t17 = pkin(4) * t113 + pkin(8) * t166 + t26;
t164 = t120 * t135;
t130 = -pkin(2) * t140 - pkin(1);
t84 = t118 * pkin(3) - t120 * qJ(4) + t130;
t124 = t178 * t138;
t125 = t178 * t140;
t95 = t124 * t134 - t125 * t136;
t47 = -t133 * t95 + t135 * t84;
t29 = pkin(4) * t118 - pkin(8) * t164 + t47;
t165 = t120 * t133;
t48 = t133 * t84 + t135 * t95;
t43 = -pkin(8) * t165 + t48;
t175 = t137 * t29 + t139 * t43;
t167 = t114 * t133;
t27 = t133 * t55 + t135 * t68;
t20 = pkin(8) * t167 + t27;
t4 = -qJD(5) * t175 - t137 * t20 + t139 * t17;
t190 = 2 * m(5);
t189 = 2 * m(6);
t188 = 0.2e1 * m(7);
t132 = t135 ^ 2;
t186 = 0.2e1 * t130;
t185 = m(4) * pkin(2);
t182 = pkin(2) * t134;
t181 = pkin(2) * t136;
t67 = t112 * t134 - t136 * t142;
t94 = -t136 * t124 - t125 * t134;
t180 = t67 * t94;
t126 = qJ(4) + t182;
t179 = pkin(8) + t126;
t44 = t114 * t119 - t116 * t120;
t22 = mrSges(6,1) * t113 - mrSges(6,3) * t44;
t23 = -t113 * mrSges(7,1) + t44 * mrSges(7,2);
t177 = -t22 + t23;
t159 = qJD(5) * t139;
t160 = qJD(5) * t137;
t45 = -t114 * t121 + t159 * t164 - t160 * t165;
t24 = -mrSges(6,2) * t113 - mrSges(6,3) * t45;
t25 = -mrSges(7,2) * t45 + mrSges(7,3) * t113;
t176 = t24 + t25;
t70 = t121 * t120;
t56 = -mrSges(6,2) * t118 - mrSges(6,3) * t70;
t59 = -mrSges(7,2) * t70 + mrSges(7,3) * t118;
t174 = t56 + t59;
t71 = t119 * t120;
t57 = mrSges(6,1) * t118 + mrSges(6,3) * t71;
t58 = -mrSges(7,1) * t118 - mrSges(7,2) * t71;
t173 = -t57 + t58;
t66 = -mrSges(5,1) * t167 - mrSges(5,2) * t166;
t172 = Ifges(5,4) * t133;
t171 = Ifges(5,4) * t135;
t170 = t113 * Ifges(5,5);
t169 = t113 * Ifges(5,6);
t168 = t133 * Ifges(5,2);
t158 = t137 * qJD(4);
t157 = t139 * qJD(4);
t156 = 0.2e1 * t140;
t129 = -pkin(3) - t181;
t15 = t45 * mrSges(6,1) + t44 * mrSges(6,2);
t14 = t45 * mrSges(7,1) - t44 * mrSges(7,3);
t117 = t179 * t135;
t151 = qJD(5) * t179;
t60 = t135 * t157 - t117 * t160 + (-t139 * t151 - t158) * t133;
t61 = t135 * t158 + t117 * t159 + (-t137 * t151 + t157) * t133;
t153 = t179 * t133;
t76 = t137 * t117 + t139 * t153;
t77 = t139 * t117 - t137 * t153;
t154 = t77 * t60 + t61 * t76;
t150 = t113 * mrSges(4,1) - t114 * mrSges(4,2);
t50 = -pkin(4) * t167 + t67;
t65 = pkin(4) * t165 + t94;
t148 = Ifges(5,5) * t135 - Ifges(5,6) * t133;
t8 = -t137 * t43 + t139 * t29;
t123 = -pkin(4) * t135 + t129;
t145 = t195 * t113 + t196 * t45 + t197 * t44;
t3 = t137 * t17 + t139 * t20 + t29 * t159 - t160 * t43;
t91 = Ifges(6,1) * t121 - Ifges(6,4) * t119;
t90 = Ifges(7,1) * t121 + Ifges(7,5) * t119;
t89 = Ifges(6,4) * t121 - Ifges(6,2) * t119;
t88 = Ifges(7,5) * t121 + Ifges(7,3) * t119;
t87 = mrSges(7,1) * t119 - mrSges(7,3) * t121;
t86 = mrSges(5,1) * t118 - mrSges(5,3) * t164;
t85 = -mrSges(5,2) * t118 - mrSges(5,3) * t165;
t83 = -Ifges(6,1) * t115 - Ifges(6,4) * t116;
t82 = -Ifges(7,1) * t115 + Ifges(7,5) * t116;
t81 = -Ifges(6,4) * t115 - Ifges(6,2) * t116;
t80 = -Ifges(7,5) * t115 + Ifges(7,3) * t116;
t75 = mrSges(5,1) * t113 + mrSges(5,3) * t166;
t74 = -mrSges(5,2) * t113 + mrSges(5,3) * t167;
t69 = pkin(5) * t119 - qJ(6) * t121 + t123;
t54 = t170 - (t135 * Ifges(5,1) - t172) * t114;
t53 = t169 - (-t168 + t171) * t114;
t46 = mrSges(7,1) * t70 + mrSges(7,3) * t71;
t33 = -Ifges(6,1) * t71 - Ifges(6,4) * t70 + Ifges(6,5) * t118;
t32 = -Ifges(7,1) * t71 + Ifges(7,4) * t118 + Ifges(7,5) * t70;
t31 = -Ifges(6,4) * t71 - Ifges(6,2) * t70 + Ifges(6,6) * t118;
t30 = -Ifges(7,5) * t71 + Ifges(7,6) * t118 + Ifges(7,3) * t70;
t19 = pkin(5) * t70 + qJ(6) * t71 + t65;
t13 = Ifges(6,1) * t44 - Ifges(6,4) * t45 + t113 * Ifges(6,5);
t12 = Ifges(7,1) * t44 + t113 * Ifges(7,4) + Ifges(7,5) * t45;
t11 = Ifges(6,4) * t44 - Ifges(6,2) * t45 + t113 * Ifges(6,6);
t10 = Ifges(7,5) * t44 + t113 * Ifges(7,6) + Ifges(7,3) * t45;
t7 = -pkin(5) * t118 - t8;
t6 = qJ(6) * t118 + t175;
t5 = pkin(5) * t45 - qJ(6) * t44 + qJD(6) * t71 + t50;
t2 = -pkin(5) * t113 - t4;
t1 = qJ(6) * t113 + qJD(6) * t118 + t3;
t9 = [0.2e1 * (-t113 * t95 - t114 * t94) * mrSges(4,3) + t113 * (-Ifges(6,5) * t71 - Ifges(6,6) * t70) + t113 * (-Ifges(7,4) * t71 + Ifges(7,6) * t70) + (-t133 * t53 + t135 * t54 + (-(2 * Ifges(4,4)) + t148) * t113 - (Ifges(5,1) * t132 + (2 * Ifges(4,1)) + (t168 - 0.2e1 * t171) * t133) * t114 + 0.2e1 * (mrSges(5,1) * t133 + mrSges(5,2) * t135 + mrSges(4,3)) * t67) * t120 + (t175 * t3 + t4 * t8 + t50 * t65) * t189 + 0.2e1 * t175 * t24 + 0.2e1 * t50 * (mrSges(6,1) * t70 - mrSges(6,2) * t71) + t150 * t186 + (t1 * t6 + t19 * t5 + t2 * t7) * t188 + (t26 * t47 + t27 * t48 + t180) * t190 + (-0.2e1 * t68 * mrSges(4,3) - 0.2e1 * (-Ifges(4,4) + t148) * t114 + ((2 * Ifges(4,2)) + (2 * Ifges(5,3)) + t195) * t113 + t145) * t118 + (t30 - t31) * t45 + (t32 + t33) * t44 + 0.2e1 * t94 * t66 + 0.2e1 * t27 * t85 + 0.2e1 * t26 * t86 - t70 * t11 - t71 * t12 - t71 * t13 + 0.2e1 * t48 * t74 + 0.2e1 * t47 * t75 + t70 * t10 + 0.2e1 * t65 * t15 + 0.2e1 * t3 * t56 + 0.2e1 * t4 * t57 + 0.2e1 * t2 * t58 + 0.2e1 * t1 * t59 + 0.2e1 * t5 * t46 + 0.2e1 * t8 * t22 + 0.2e1 * t7 * t23 + 0.2e1 * t6 * t25 + 0.2e1 * t19 * t14 + ((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t140) * t156 + (-0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * pkin(2) * (mrSges(4,1) * t118 + mrSges(4,2) * t120) + t185 * t186 - 0.2e1 * Ifges(3,4) * t138 + (-Ifges(3,2) + Ifges(3,1)) * t156) * t138) * qJD(2) + 0.2e1 * m(4) * (t68 * t95 + t180); (-t31 / 0.2e1 + t30 / 0.2e1) * t116 - (t32 / 0.2e1 + t33 / 0.2e1) * t115 + (t115 * t8 - t116 * t175 - t119 * t3 - t121 * t4) * mrSges(6,3) + (-t1 * t119 - t115 * t7 - t116 * t6 + t121 * t2) * mrSges(7,2) + m(6) * (t123 * t50 + t175 * t60 + t3 * t77 - t4 * t76 - t61 * t8) + (t134 * t68 - t136 * t67) * t185 + t192 * t118 / 0.2e1 + (t50 * mrSges(6,2) + t12 / 0.2e1 + t13 / 0.2e1) * t121 + m(5) * (t129 * t67 + (-t133 * t26 + t135 * t27) * t126 + (-t133 * t47 + t135 * t48) * qJD(4)) + (-t11 / 0.2e1 + t50 * mrSges(6,1) + t10 / 0.2e1) * t119 + m(7) * (t1 * t77 + t19 * t63 + t2 * t76 + t5 * t69 + t6 * t60 + t61 * t7) + (Ifges(3,5) * t140 - Ifges(3,6) * t138 + (-mrSges(3,1) * t140 + mrSges(3,2) * t138) * pkin(7)) * qJD(2) - (t82 / 0.2e1 + t83 / 0.2e1) * t71 + (t80 / 0.2e1 - t81 / 0.2e1) * t70 + t129 * t66 + t123 * t15 + t19 * t78 + t65 * t79 + t5 * t87 - t67 * mrSges(4,1) - t68 * mrSges(4,2) + t69 * t14 + t63 * t46 + (t90 / 0.2e1 + t91 / 0.2e1) * t44 + (t88 / 0.2e1 - t89 / 0.2e1) * t45 + (t169 / 0.2e1 - t67 * mrSges(5,1) + t53 / 0.2e1 + t27 * mrSges(5,3) + qJD(4) * t85 + t126 * t74) * t135 + (t54 / 0.2e1 + t170 / 0.2e1 + t67 * mrSges(5,2) - t26 * mrSges(5,3) - qJD(4) * t86 - t126 * t75) * t133 + t173 * t61 + t174 * t60 + t176 * t77 + t177 * t76 - (-mrSges(4,3) * t181 + Ifges(4,5) + t135 * (Ifges(5,1) * t133 + t171) / 0.2e1 - t133 * (Ifges(5,2) * t135 + t172) / 0.2e1) * t114 + (-mrSges(4,3) * t182 - Ifges(4,6) + (Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1) * t121 + (-Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1) * t119) * t113; 0.2e1 * t123 * t79 + 0.2e1 * t63 * t87 + 0.2e1 * t69 * t78 + (t82 + t83) * t121 + (-t81 + t80) * t119 + (-t89 + t88) * t116 - (t90 + t91) * t115 + (t63 * t69 + t154) * t188 + t154 * t189 + (t126 * t190 + 0.2e1 * mrSges(5,3)) * qJD(4) * (t133 ^ 2 + t132) + 0.2e1 * (mrSges(6,3) + mrSges(7,2)) * (-t115 * t76 - t116 * t77 - t119 * t60 + t121 * t61); m(4) * t155 + t133 * t74 + t135 * t75 + t176 * t121 + t177 * t119 + t173 * t116 - t174 * t115 + m(7) * (t1 * t121 - t115 * t6 + t116 * t7 + t119 * t2) + m(6) * (-t115 * t175 - t116 * t8 - t119 * t4 + t121 * t3) + m(5) * (t133 * t27 + t135 * t26) + t150; t194 * (-t115 * t77 + t116 * t76 + t119 * t61 + t121 * t60); 0.2e1 * t194 * (-t121 * t115 + t116 * t119); m(5) * t67 + m(6) * t50 + m(7) * t5 + t14 + t15 + t66; -t193; 0; 0; qJD(6) * t59 + qJ(6) * t25 + m(7) * (-pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t6) + t1 * mrSges(7,3) - t2 * mrSges(7,1) + t4 * mrSges(6,1) - t3 * mrSges(6,2) - pkin(5) * t23 + t145; m(7) * qJD(6) * t77 + (pkin(5) * t115 - qJ(6) * t116 - qJD(6) * t119) * mrSges(7,2) + (-m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1)) * t61 + (-mrSges(6,2) + t191) * t60 + t192; t193; 0; 0.2e1 * t191 * qJD(6); m(7) * t2 + t23; m(7) * t61 - t115 * mrSges(7,2); m(7) * t116; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t9(1) t9(2) t9(4) t9(7) t9(11) t9(16); t9(2) t9(3) t9(5) t9(8) t9(12) t9(17); t9(4) t9(5) t9(6) t9(9) t9(13) t9(18); t9(7) t9(8) t9(9) t9(10) t9(14) t9(19); t9(11) t9(12) t9(13) t9(14) t9(15) t9(20); t9(16) t9(17) t9(18) t9(19) t9(20) t9(21);];
Mq  = res;
