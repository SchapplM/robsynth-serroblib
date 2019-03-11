% Calculate time derivative of joint inertia matrix for
% S6RPPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 02:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:22:56
% EndTime: 2019-03-09 02:23:00
% DurationCPUTime: 2.41s
% Computational Cost: add. (2596->346), mult. (5686->524), div. (0->0), fcn. (4626->8), ass. (0->151)
t162 = 2 * qJD(3);
t101 = sin(qJ(5));
t102 = sin(qJ(4));
t138 = qJD(4) * t102;
t122 = t101 * t138;
t104 = cos(qJ(5));
t105 = cos(qJ(4));
t134 = qJD(5) * t105;
t123 = t104 * t134;
t108 = t122 - t123;
t124 = t101 * t134;
t109 = t104 * t138 + t124;
t36 = -t108 * mrSges(6,1) - t109 * mrSges(6,2);
t100 = sin(qJ(6));
t103 = cos(qJ(6));
t111 = t100 * t101 - t103 * t104;
t166 = qJD(5) + qJD(6);
t73 = t100 * t104 + t101 * t103;
t42 = t166 * t73;
t23 = -t105 * t42 + t111 * t138;
t173 = t166 * t111;
t25 = t105 * t173 + t138 * t73;
t6 = -mrSges(7,1) * t25 + t23 * mrSges(7,2);
t174 = -t36 - t6;
t149 = t101 ^ 2 + t104 ^ 2;
t136 = qJD(5) * t101;
t125 = t102 * t136;
t137 = qJD(4) * t105;
t172 = t104 * t137 - t125;
t15 = mrSges(7,1) * t42 - mrSges(7,2) * t173;
t115 = mrSges(6,1) * t101 + mrSges(6,2) * t104;
t74 = t115 * qJD(5);
t171 = -t74 - t15 + (mrSges(6,3) * t149 - mrSges(5,2)) * qJD(4);
t170 = Ifges(6,6) * t122 + Ifges(6,3) * t137;
t96 = t102 ^ 2;
t98 = t105 ^ 2;
t169 = t96 - t98;
t151 = t102 * pkin(4);
t152 = pkin(8) * t105;
t89 = sin(pkin(10)) * pkin(1) + qJ(3);
t66 = t151 + t89 - t152;
t140 = t102 * t104;
t88 = -cos(pkin(10)) * pkin(1) - pkin(2) - pkin(7);
t71 = t88 * t140;
t38 = t101 * t66 + t71;
t168 = t38 * qJD(5);
t81 = -mrSges(6,1) * t104 + mrSges(6,2) * t101;
t150 = t81 - mrSges(5,1);
t43 = mrSges(7,1) * t111 + mrSges(7,2) * t73;
t90 = -pkin(5) * t104 - pkin(4);
t167 = m(7) * t90 + t150 + t43;
t164 = 0.2e1 * m(7);
t163 = -0.2e1 * t88;
t161 = m(7) * pkin(5);
t160 = -t111 / 0.2e1;
t159 = t73 / 0.2e1;
t158 = -pkin(9) - pkin(8);
t87 = t102 * t137;
t157 = m(6) * (t149 - 0.1e1) * t87;
t155 = -t101 / 0.2e1;
t154 = pkin(4) * t105;
t153 = pkin(8) * t102;
t148 = Ifges(6,4) * t101;
t147 = Ifges(6,4) * t104;
t146 = Ifges(6,5) * t101;
t145 = Ifges(6,6) * t101;
t144 = Ifges(6,6) * t102;
t143 = Ifges(6,6) * t104;
t142 = t101 * t102;
t141 = t101 * t105;
t139 = t104 * t105;
t135 = qJD(5) * t104;
t133 = qJD(6) * t100;
t132 = qJD(6) * t103;
t130 = Ifges(7,5) * t23 + Ifges(7,6) * t25 + Ifges(7,3) * t137;
t129 = pkin(5) * t136;
t60 = t73 * t105;
t62 = t111 * t105;
t34 = mrSges(7,1) * t60 - mrSges(7,2) * t62;
t63 = (pkin(5) * t101 - t88) * t105;
t128 = m(7) * t63 + t34;
t127 = qJD(5) * t158;
t126 = t101 * t137;
t22 = -qJD(4) * t62 - t102 * t42;
t24 = -qJD(4) * t60 + t102 * t173;
t120 = mrSges(7,1) * t24 - t22 * mrSges(7,2);
t119 = -t101 * t88 + pkin(5);
t118 = -Ifges(6,5) * t104 + (2 * Ifges(5,4));
t70 = qJD(3) + (t153 + t154) * qJD(4);
t13 = t101 * t70 + t66 * t135 + t172 * t88;
t58 = t104 * t66;
t37 = -t142 * t88 + t58;
t117 = -t37 * qJD(5) + t13;
t67 = t115 * t105;
t116 = t128 + t67;
t114 = Ifges(6,1) * t104 - t148;
t83 = Ifges(6,1) * t101 + t147;
t113 = -Ifges(6,2) * t101 + t147;
t82 = Ifges(6,2) * t104 + t148;
t31 = -pkin(9) * t139 + t102 * t119 + t58;
t35 = -pkin(9) * t141 + t38;
t7 = -t100 * t35 + t103 * t31;
t8 = t100 * t31 + t103 * t35;
t84 = t158 * t101;
t85 = t158 * t104;
t48 = t100 * t85 + t103 * t84;
t49 = t100 * t84 - t103 * t85;
t79 = t101 * t127;
t80 = t104 * t127;
t27 = qJD(6) * t48 + t100 * t80 + t103 * t79;
t28 = -qJD(6) * t49 - t100 * t79 + t103 * t80;
t39 = Ifges(7,6) * t42;
t40 = Ifges(7,5) * t173;
t112 = mrSges(7,1) * t28 - t27 * mrSges(7,2) - t39 - t40;
t10 = pkin(9) * t108 + t13;
t65 = t104 * t70;
t9 = t65 + (-t71 + (pkin(9) * t105 - t66) * t101) * qJD(5) + (pkin(9) * t140 + t105 * t119) * qJD(4);
t2 = qJD(6) * t7 + t10 * t103 + t100 * t9;
t3 = -qJD(6) * t8 - t10 * t100 + t103 * t9;
t110 = mrSges(7,1) * t3 - mrSges(7,2) * t2 + t130;
t14 = -t126 * t88 - t168 + t65;
t52 = mrSges(6,1) * t137 + mrSges(6,3) * t109;
t53 = -mrSges(6,2) * t137 + mrSges(6,3) * t108;
t77 = -mrSges(6,2) * t102 - mrSges(6,3) * t141;
t78 = mrSges(6,1) * t102 - mrSges(6,3) * t139;
t106 = m(6) * (-t101 * t14 + t104 * t13 - t135 * t37 - t136 * t38) + t104 * t53 - t101 * t52 - t77 * t136 - t78 * t135;
t94 = qJD(4) * t96;
t93 = Ifges(6,5) * t135;
t76 = t114 * qJD(5);
t75 = t113 * qJD(5);
t69 = (-mrSges(7,1) * t100 - mrSges(7,2) * t103) * qJD(6) * pkin(5);
t61 = t111 * t102;
t59 = t73 * t102;
t56 = Ifges(6,5) * t102 + t105 * t114;
t55 = t105 * t113 + t144;
t51 = mrSges(7,1) * t102 + mrSges(7,3) * t62;
t50 = -mrSges(7,2) * t102 - mrSges(7,3) * t60;
t46 = -pkin(5) * t108 + t138 * t88;
t45 = Ifges(7,1) * t73 - Ifges(7,4) * t111;
t44 = Ifges(7,4) * t73 - Ifges(7,2) * t111;
t33 = -t83 * t134 + (Ifges(6,5) * t105 - t102 * t114) * qJD(4);
t32 = -t82 * t134 + (Ifges(6,6) * t105 - t102 * t113) * qJD(4);
t30 = -Ifges(7,1) * t62 - Ifges(7,4) * t60 + Ifges(7,5) * t102;
t29 = -Ifges(7,4) * t62 - Ifges(7,2) * t60 + Ifges(7,6) * t102;
t17 = -Ifges(7,1) * t173 - Ifges(7,4) * t42;
t16 = -Ifges(7,4) * t173 - Ifges(7,2) * t42;
t12 = -mrSges(7,2) * t137 + mrSges(7,3) * t25;
t11 = mrSges(7,1) * t137 - mrSges(7,3) * t23;
t5 = Ifges(7,1) * t23 + Ifges(7,4) * t25 + Ifges(7,5) * t137;
t4 = Ifges(7,4) * t23 + Ifges(7,2) * t25 + Ifges(7,6) * t137;
t1 = [0.2e1 * t7 * t11 + 0.2e1 * t8 * t12 + 0.2e1 * t13 * t77 + 0.2e1 * t14 * t78 + 0.2e1 * t2 * t50 + t23 * t30 + t25 * t29 + 0.2e1 * t3 * t51 + 0.2e1 * t46 * t34 + 0.2e1 * t37 * t52 + 0.2e1 * t38 * t53 - t60 * t4 - t62 * t5 + 0.2e1 * t63 * t6 + (t2 * t8 + t3 * t7 + t46 * t63) * t164 + 0.2e1 * m(6) * (t13 * t38 + t14 * t37) + (mrSges(4,3) + (m(4) + m(5)) * t89) * t162 + (mrSges(5,1) * t162 + (-0.2e1 * t89 * mrSges(5,2) + t101 * t55 + t102 * t118 - t104 * t56 + 0.2e1 * t88 * t67) * qJD(4) + t130 + t170) * t102 + (mrSges(5,2) * t162 - t101 * t32 + t104 * t33 + t36 * t163 + (-t104 * t55 - t101 * t56 + t102 * (-t143 - t146)) * qJD(5) + (-Ifges(7,5) * t62 - Ifges(7,6) * t60 + 0.2e1 * t89 * mrSges(5,1) + (-t118 - t145) * t105 + (-0.2e1 * m(6) * t88 ^ 2 - (2 * Ifges(5,1)) + (2 * Ifges(5,2)) + Ifges(6,3) + Ifges(7,3)) * t102) * qJD(4)) * t105; t23 * t50 - t62 * t12 + t25 * t51 - t60 * t11 + m(7) * (-t2 * t62 + t23 * t8 + t25 * t7 - t3 * t60) + (-t77 * t140 + t78 * t142 + m(6) * (-t38 * t140 + t37 * t142 + t169 * t88)) * qJD(4) + (qJD(4) * t116 + t106) * t105 + (m(7) * t46 - t174) * t102; -0.2e1 * t157 + 0.2e1 * m(7) * (-t23 * t62 - t25 * t60 + t87); m(7) * (-t105 * t46 - t2 * t61 + t22 * t8 + t24 * t7 - t3 * t59) + t22 * t50 - t61 * t12 + t24 * t51 - t59 * t11 - t105 * t6 - t105 * t36 + (m(6) * (-t101 * t37 + t104 * t38) + t104 * t77 - t101 * t78) * t137 + ((m(6) * t105 * t163 + t116) * qJD(4) + t106) * t102; (t94 + (-t149 * t169 - t98) * qJD(4)) * m(6) + m(7) * (-qJD(4) * t98 - t22 * t62 - t23 * t61 - t24 * t60 - t25 * t59 + t94); 0.2e1 * t157 + 0.2e1 * m(7) * (-t22 * t61 - t24 * t59 - t87); t23 * t45 / 0.2e1 + t46 * t43 + t48 * t11 + t49 * t12 + t27 * t50 + t28 * t51 - pkin(4) * t36 - t173 * t30 / 0.2e1 - t42 * t29 / 0.2e1 + t25 * t44 / 0.2e1 + m(7) * (t2 * t49 + t27 * t8 + t28 * t7 + t3 * t48 + t46 * t90) - t60 * t16 / 0.2e1 - t62 * t17 / 0.2e1 + t63 * t15 + t4 * t160 + t5 * t159 + t90 * t6 + (-t40 / 0.2e1 - t39 / 0.2e1 + t93 / 0.2e1 + (-Ifges(5,5) + (-m(6) * pkin(4) + t150) * t88) * qJD(4)) * t102 + (-t14 * mrSges(6,3) + t33 / 0.2e1 + t82 * t138 / 0.2e1 + (-t55 / 0.2e1 - t144 / 0.2e1 - t38 * mrSges(6,3) + t128 * pkin(5)) * qJD(5) + (m(6) * (-t14 - t168) - qJD(5) * t77 - t52) * pkin(8)) * t101 + (-t111 * t2 + t173 * t7 - t3 * t73 - t8 * t42) * mrSges(7,3) + (qJD(5) * t56 / 0.2e1 - t83 * t138 / 0.2e1 + t32 / 0.2e1 + t117 * mrSges(6,3) + (m(6) * t117 - qJD(5) * t78 + t53) * pkin(8)) * t104 + (-t88 * t74 + t75 * t155 + t104 * t76 / 0.2e1 + (-t104 * t82 / 0.2e1 + t83 * t155) * qJD(5) + (Ifges(7,5) * t159 + Ifges(7,6) * t160 + t146 / 0.2e1 + t143 / 0.2e1 - Ifges(5,6) - t88 * mrSges(5,2)) * qJD(4)) * t105; m(7) * (pkin(5) * t125 + t23 * t49 + t25 * t48 - t27 * t62 - t28 * t60) + (-t111 * t23 - t173 * t60 - t25 * t73 + t42 * t62) * mrSges(7,3) + (m(6) * (-t149 * t153 - t154) + t167 * t105) * qJD(4) - t171 * t102; m(7) * (-pkin(5) * t124 + t22 * t49 + t24 * t48 - t27 * t61 - t28 * t59) + (-t111 * t22 - t173 * t59 - t24 * t73 + t42 * t61) * mrSges(7,3) + (m(6) * (t149 * t152 - t151) + t167 * t102) * qJD(4) + t171 * t105; -t173 * t45 + t73 * t17 - t42 * t44 - t111 * t16 + 0.2e1 * t43 * t129 + 0.2e1 * t90 * t15 + (t129 * t90 + t27 * t49 + t28 * t48) * t164 - t82 * t136 - 0.2e1 * pkin(4) * t74 + t101 * t76 + (qJD(5) * t83 + t75) * t104 + 0.2e1 * (-t111 * t27 + t173 * t48 - t28 * t73 - t42 * t49) * mrSges(7,3); -Ifges(6,6) * t123 + t14 * mrSges(6,1) - t13 * mrSges(6,2) - t109 * Ifges(6,5) + (m(7) * (t100 * t2 + t103 * t3 + t132 * t8 - t133 * t7) + t50 * t132 + t100 * t12 - t51 * t133 + t103 * t11) * pkin(5) + t110 + t170; (t100 * t23 + t103 * t25 + (t100 * t60 - t103 * t62) * qJD(6)) * t161 + t174; -t172 * mrSges(6,2) + (-t102 * t135 - t126) * mrSges(6,1) + (t100 * t22 + t103 * t24 + (t100 * t59 - t103 * t61) * qJD(6)) * t161 + t120; t93 + (pkin(8) * t81 - t145) * qJD(5) + (m(7) * (t100 * t27 + t103 * t28 + (-t100 * t48 + t103 * t49) * qJD(6)) + (-t100 * t42 + t103 * t173 + (t100 * t73 - t103 * t111) * qJD(6)) * mrSges(7,3)) * pkin(5) + t112; 0.2e1 * t69; t110; -t6; t120; t112; t69; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
