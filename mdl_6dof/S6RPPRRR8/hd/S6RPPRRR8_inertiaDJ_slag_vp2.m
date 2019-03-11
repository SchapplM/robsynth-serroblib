% Calculate time derivative of joint inertia matrix for
% S6RPPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 02:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR8_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR8_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR8_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR8_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:35:25
% EndTime: 2019-03-09 02:35:30
% DurationCPUTime: 2.34s
% Computational Cost: add. (4475->359), mult. (9244->542), div. (0->0), fcn. (8959->8), ass. (0->146)
t114 = sin(pkin(10));
t115 = cos(pkin(10));
t167 = sin(qJ(4));
t168 = cos(qJ(4));
t85 = t167 * t114 - t168 * t115;
t181 = 2 * qJD(2);
t117 = sin(qJ(6));
t118 = sin(qJ(5));
t119 = cos(qJ(6));
t120 = cos(qJ(5));
t88 = t117 * t120 + t119 * t118;
t55 = t88 * t85;
t124 = t117 * t118 - t119 * t120;
t177 = qJD(5) + qJD(6);
t69 = t177 * t124;
t70 = t177 * t88;
t33 = t70 * mrSges(7,1) - t69 * mrSges(7,2);
t130 = mrSges(6,1) * t118 + mrSges(6,2) * t120;
t89 = t130 * qJD(5);
t180 = t33 + t89;
t136 = qJD(4) * t167;
t137 = qJD(4) * t168;
t81 = -t114 * t137 - t115 * t136;
t155 = t120 * t81;
t82 = -t114 * t136 + t115 * t137;
t179 = Ifges(6,5) * t155 + Ifges(6,3) * t82;
t86 = t114 * t168 + t115 * t167;
t56 = t124 * t86;
t116 = -pkin(1) - qJ(3);
t161 = -pkin(7) + t116;
t94 = t161 * t114;
t95 = t161 * t115;
t178 = -t167 * t94 + t168 * t95;
t132 = (t114 ^ 2 + t115 ^ 2) * qJD(3);
t176 = 2 * m(7);
t175 = -2 * mrSges(5,3);
t174 = -0.2e1 * t178;
t173 = m(7) * pkin(5);
t172 = t85 / 0.2e1;
t158 = Ifges(6,4) * t120;
t98 = Ifges(6,1) * t118 + t158;
t170 = t98 / 0.2e1;
t169 = -pkin(9) - pkin(8);
t166 = Ifges(6,5) * t82;
t165 = Ifges(6,6) * t82;
t164 = Ifges(6,6) * t86;
t68 = t167 * t95 + t168 * t94;
t47 = -t85 * qJD(3) + qJD(4) * t68;
t163 = t47 * t178;
t71 = t85 * t81;
t96 = -mrSges(6,1) * t120 + t118 * mrSges(6,2);
t162 = mrSges(5,1) - t96;
t105 = t114 * pkin(3) + qJ(2);
t59 = pkin(4) * t86 + pkin(8) * t85 + t105;
t61 = t120 * t68;
t32 = t118 * t59 + t61;
t160 = -Ifges(7,5) * t69 - Ifges(7,6) * t70;
t159 = Ifges(6,4) * t118;
t157 = Ifges(6,6) * t118;
t156 = t118 * t85;
t154 = t120 * t85;
t152 = t118 ^ 2 + t120 ^ 2;
t151 = qJD(5) * t118;
t150 = qJD(5) * t120;
t149 = qJD(6) * t117;
t148 = qJD(6) * t119;
t147 = t82 * t175;
t19 = -t124 * t81 + t177 * t55;
t21 = -t69 * t85 - t81 * t88;
t146 = Ifges(7,5) * t19 + Ifges(7,6) * t21 + Ifges(7,3) * t82;
t145 = pkin(5) * t151;
t143 = t85 * t151;
t142 = qJD(5) * t169;
t139 = t82 * mrSges(5,1) + t81 * mrSges(5,2);
t18 = -t124 * t82 - t70 * t86;
t20 = t177 * t56 - t88 * t82;
t138 = t20 * mrSges(7,1) - t18 * mrSges(7,2);
t135 = -(2 * Ifges(5,4)) - t157;
t46 = -t86 * qJD(3) + t178 * qJD(4);
t58 = pkin(4) * t82 - pkin(8) * t81 + qJD(2);
t134 = -t118 * t46 + t120 * t58;
t31 = -t118 * t68 + t120 * t59;
t133 = t152 * t82;
t131 = t178 * t81 + t85 * t47;
t129 = Ifges(6,1) * t120 - t159;
t128 = -Ifges(6,2) * t118 + t158;
t22 = t86 * pkin(5) + pkin(9) * t154 + t31;
t27 = pkin(9) * t156 + t32;
t9 = -t117 * t27 + t119 * t22;
t10 = t117 * t22 + t119 * t27;
t127 = -t118 * t31 + t120 * t32;
t62 = -mrSges(6,2) * t86 + mrSges(6,3) * t156;
t63 = t86 * mrSges(6,1) + mrSges(6,3) * t154;
t126 = -t118 * t63 + t120 * t62;
t100 = t169 * t120;
t99 = t169 * t118;
t75 = t100 * t117 + t119 * t99;
t92 = t118 * t142;
t93 = t120 * t142;
t37 = qJD(6) * t75 + t117 * t93 + t119 * t92;
t76 = -t100 * t119 + t117 * t99;
t38 = -qJD(6) * t76 - t117 * t92 + t119 * t93;
t125 = t38 * mrSges(7,1) - t37 * mrSges(7,2) + t160;
t7 = -pkin(9) * t155 + t82 * pkin(5) + (-t61 + (-pkin(9) * t85 - t59) * t118) * qJD(5) + t134;
t11 = t118 * t58 + t120 * t46 + t59 * t150 - t151 * t68;
t122 = -t118 * t81 + t150 * t85;
t8 = pkin(9) * t122 + t11;
t2 = qJD(6) * t9 + t117 * t7 + t119 * t8;
t3 = -qJD(6) * t10 - t117 * t8 + t119 * t7;
t123 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + t146;
t121 = t143 + t155;
t108 = Ifges(6,5) * t150;
t107 = -pkin(5) * t120 - pkin(4);
t97 = Ifges(6,2) * t120 + t159;
t91 = t129 * qJD(5);
t90 = t128 * qJD(5);
t84 = (-mrSges(7,1) * t117 - mrSges(7,2) * t119) * qJD(6) * pkin(5);
t74 = Ifges(7,1) * t88 - Ifges(7,4) * t124;
t73 = Ifges(7,4) * t88 - Ifges(7,2) * t124;
t72 = mrSges(7,1) * t124 + mrSges(7,2) * t88;
t60 = t130 * t85;
t57 = t124 * t85;
t54 = t88 * t86;
t48 = -pkin(5) * t156 - t178;
t45 = Ifges(6,5) * t86 - t129 * t85;
t44 = -t128 * t85 + t164;
t42 = mrSges(7,1) * t86 - mrSges(7,3) * t57;
t41 = -mrSges(7,2) * t86 + mrSges(7,3) * t55;
t40 = -t82 * mrSges(6,2) + mrSges(6,3) * t122;
t39 = t82 * mrSges(6,1) - mrSges(6,3) * t121;
t35 = -Ifges(7,1) * t69 - Ifges(7,4) * t70;
t34 = -Ifges(7,4) * t69 - Ifges(7,2) * t70;
t30 = -mrSges(6,1) * t122 + mrSges(6,2) * t121;
t29 = -mrSges(7,1) * t55 + mrSges(7,2) * t57;
t28 = -pkin(5) * t122 + t47;
t26 = Ifges(7,1) * t57 + Ifges(7,4) * t55 + Ifges(7,5) * t86;
t25 = Ifges(7,4) * t57 + Ifges(7,2) * t55 + Ifges(7,6) * t86;
t24 = Ifges(6,1) * t121 + Ifges(6,4) * t122 + t166;
t23 = Ifges(6,4) * t121 + Ifges(6,2) * t122 + t165;
t14 = -mrSges(7,2) * t82 + mrSges(7,3) * t21;
t13 = mrSges(7,1) * t82 - mrSges(7,3) * t19;
t12 = -qJD(5) * t32 + t134;
t6 = -mrSges(7,1) * t21 + mrSges(7,2) * t19;
t5 = Ifges(7,1) * t19 + Ifges(7,4) * t21 + t82 * Ifges(7,5);
t4 = Ifges(7,4) * t19 + Ifges(7,2) * t21 + t82 * Ifges(7,6);
t1 = [t30 * t174 + (t10 * t2 + t28 * t48 + t3 * t9) * t176 + (mrSges(5,1) * t181 + t46 * t175 + t135 * t81 + ((2 * Ifges(5,2)) + Ifges(6,3) + Ifges(7,3)) * t82 + t146 + t179) * t86 + t68 * t147 + 0.2e1 * mrSges(4,3) * t132 + (m(3) * qJ(2) + mrSges(4,1) * t114 + mrSges(4,2) * t115 + mrSges(3,3)) * t181 + 0.2e1 * m(6) * (t11 * t32 + t12 * t31 - t163) + 0.2e1 * m(5) * (qJD(2) * t105 + t46 * t68 - t163) + t82 * (Ifges(7,5) * t57 + Ifges(7,6) * t55) + 0.2e1 * m(4) * (qJ(2) * qJD(2) - t116 * t132) + t55 * t4 + t57 * t5 - 0.2e1 * t47 * t60 + 0.2e1 * t11 * t62 + 0.2e1 * t12 * t63 + 0.2e1 * t48 * t6 + 0.2e1 * t31 * t39 + 0.2e1 * t32 * t40 + 0.2e1 * t2 * t41 + 0.2e1 * t3 * t42 + t21 * t25 + t19 * t26 + 0.2e1 * t28 * t29 + 0.2e1 * t9 * t13 + 0.2e1 * t10 * t14 + 0.2e1 * t105 * t139 + (mrSges(5,3) * t174 - t118 * t44 + t120 * t45) * t81 + (-0.2e1 * qJD(2) * mrSges(5,2) + t47 * t175 - 0.2e1 * Ifges(5,1) * t81 + t118 * t23 - t120 * t24 + (-Ifges(6,5) * t120 - t135) * t82 + (t120 * t44 + t118 * t45 + t86 * (Ifges(6,5) * t118 + Ifges(6,6) * t120)) * qJD(5)) * t85; -t54 * t13 - t56 * t14 + t18 * t41 + t20 * t42 + (t30 + t6) * t85 + t126 * t82 + (0.2e1 * t85 * mrSges(5,3) - t29 + t60) * t81 + (t147 - t118 * t39 + t120 * t40 + (-t118 * t62 - t120 * t63) * qJD(5)) * t86 + m(7) * (t10 * t18 - t2 * t56 + t20 * t9 + t28 * t85 - t3 * t54 - t48 * t81) + m(6) * (t127 * t82 + (t11 * t120 - t118 * t12 + (-t118 * t32 - t120 * t31) * qJD(5)) * t86 + t131) + m(5) * (t46 * t86 + t68 * t82 + t131) - m(4) * t132; 0.2e1 * m(7) * (-t18 * t56 - t20 * t54 - t71) + 0.2e1 * m(6) * (t133 * t86 - t71) + 0.2e1 * m(5) * (t82 * t86 - t71); t118 * t40 + t120 * t39 - t124 * t13 + t88 * t14 - t69 * t41 - t70 * t42 + t126 * qJD(5) + (m(5) + m(4)) * qJD(2) + m(7) * (-t10 * t69 - t124 * t3 + t2 * t88 - t70 * t9) + m(6) * (qJD(5) * t127 + t118 * t11 + t12 * t120) + t139; m(7) * (-t124 * t20 + t18 * t88 + t54 * t70 + t56 * t69); (t124 * t70 - t69 * t88) * t176; (-t10 * t70 - t124 * t2 - t3 * t88 + t69 * t9) * mrSges(7,3) + (t160 + t108) * t86 / 0.2e1 + t82 * (Ifges(7,5) * t88 - Ifges(7,6) * t124) / 0.2e1 - t124 * t4 / 0.2e1 - t178 * t89 + t107 * t6 + t88 * t5 / 0.2e1 + Ifges(5,5) * t81 - Ifges(5,6) * t82 + t75 * t13 + t76 * t14 - t69 * t26 / 0.2e1 - t70 * t25 / 0.2e1 + t28 * t72 + t21 * t73 / 0.2e1 + t19 * t74 / 0.2e1 + t57 * t35 / 0.2e1 - t46 * mrSges(5,2) + t48 * t33 + t55 * t34 / 0.2e1 + t37 * t41 + t38 * t42 - pkin(4) * t30 + (t24 / 0.2e1 + t166 / 0.2e1 - t12 * mrSges(6,3) + t90 * t172 - t81 * t97 / 0.2e1 + (t85 * t170 + pkin(5) * t29 - t32 * mrSges(6,3) - t164 / 0.2e1 - t44 / 0.2e1 + t48 * t173) * qJD(5) + (-m(6) * t12 - t39 + (-m(6) * t32 - t62) * qJD(5)) * pkin(8)) * t118 + m(7) * (t10 * t37 + t107 * t28 + t2 * t76 + t3 * t75 + t38 * t9) + (-m(6) * pkin(4) - t162) * t47 + (t165 / 0.2e1 + t23 / 0.2e1 + t11 * mrSges(6,3) - t85 * t91 / 0.2e1 + t81 * t170 + (t97 * t172 - t31 * mrSges(6,3) + t45 / 0.2e1) * qJD(5) + (-qJD(5) * t63 + t40 + m(6) * (-t31 * qJD(5) + t11)) * pkin(8)) * t120; t180 * t85 + (mrSges(6,3) * t152 - mrSges(5,2)) * t82 + (-t72 + t162) * t81 + m(7) * (pkin(5) * t143 - t107 * t81 + t18 * t76 + t20 * t75 - t37 * t56 - t38 * t54) + m(6) * (pkin(4) * t81 + pkin(8) * t133) + (-t124 * t18 - t20 * t88 - t54 * t69 + t56 * t70) * mrSges(7,3); m(7) * (-t124 * t38 + t37 * t88 - t69 * t76 - t70 * t75); (t107 * t145 + t37 * t76 + t38 * t75) * t176 - t69 * t74 + t88 * t35 - t70 * t73 - t124 * t34 + 0.2e1 * t72 * t145 + 0.2e1 * t107 * t33 - t97 * t151 + t118 * t91 - 0.2e1 * pkin(4) * t89 + (qJD(5) * t98 + t90) * t120 + 0.2e1 * (-t124 * t37 - t38 * t88 + t69 * t75 - t70 * t76) * mrSges(7,3); Ifges(6,5) * t143 + t12 * mrSges(6,1) - t11 * mrSges(6,2) + t122 * Ifges(6,6) + (m(7) * (t10 * t148 + t117 * t2 + t119 * t3 - t149 * t9) + t41 * t148 + t117 * t14 - t42 * t149 + t119 * t13) * pkin(5) + t123 + t179; (-t120 * t82 + t151 * t86) * mrSges(6,2) + (-t118 * t82 - t150 * t86) * mrSges(6,1) + (t117 * t18 + t119 * t20 + (t117 * t54 - t119 * t56) * qJD(6)) * t173 + t138; (-t117 * t69 - t119 * t70 + (t117 * t124 + t119 * t88) * qJD(6)) * t173 - t180; t108 + (pkin(8) * t96 - t157) * qJD(5) + (m(7) * (t117 * t37 + t119 * t38 + (-t117 * t75 + t119 * t76) * qJD(6)) + (-t117 * t70 + t119 * t69 + (t117 * t88 - t119 * t124) * qJD(6)) * mrSges(7,3)) * pkin(5) + t125; 0.2e1 * t84; t123; t138; -t33; t125; t84; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
