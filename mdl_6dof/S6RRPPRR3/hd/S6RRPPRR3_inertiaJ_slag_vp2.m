% Calculate joint inertia matrix for
% S6RRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
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
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:50
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:50:05
% EndTime: 2018-11-23 16:50:06
% DurationCPUTime: 1.47s
% Computational Cost: add. (3516->321), mult. (8308->466), div. (0->0), fcn. (9365->12), ass. (0->129)
t178 = Ifges(3,3) + Ifges(4,3);
t123 = sin(pkin(11));
t126 = cos(pkin(11));
t124 = sin(pkin(6));
t132 = cos(qJ(2));
t145 = t124 * t132;
t130 = sin(qJ(2));
t146 = t124 * t130;
t79 = t123 * t146 - t126 * t145;
t80 = (t123 * t132 + t126 * t130) * t124;
t177 = Ifges(4,5) * t80 - Ifges(4,6) * t79;
t128 = sin(qJ(6));
t131 = cos(qJ(6));
t122 = sin(pkin(12));
t125 = cos(pkin(12));
t129 = sin(qJ(5));
t162 = cos(qJ(5));
t93 = t122 * t162 + t129 * t125;
t150 = t128 * t93;
t91 = t129 * t122 - t125 * t162;
t60 = -mrSges(7,2) * t91 - mrSges(7,3) * t150;
t147 = t131 * t93;
t61 = mrSges(7,1) * t91 - mrSges(7,3) * t147;
t176 = -t128 * t61 + t131 * t60;
t127 = cos(pkin(6));
t67 = -t122 * t80 + t125 * t127;
t68 = t122 * t127 + t125 * t80;
t39 = t129 * t67 + t162 * t68;
t25 = -t128 * t39 + t131 * t79;
t38 = t129 * t68 - t162 * t67;
t14 = -mrSges(7,2) * t38 + mrSges(7,3) * t25;
t26 = t128 * t79 + t131 * t39;
t15 = mrSges(7,1) * t38 - mrSges(7,3) * t26;
t175 = -t128 * t15 + t131 * t14;
t110 = pkin(2) * t123 + qJ(4);
t157 = pkin(9) + t110;
t141 = t157 * t122;
t83 = t157 * t125;
t55 = t129 * t83 + t141 * t162;
t174 = t55 ^ 2;
t173 = t91 ^ 2;
t172 = 0.2e1 * t55;
t94 = (-pkin(2) * t132 - pkin(1)) * t124;
t171 = 0.2e1 * t94;
t170 = t25 / 0.2e1;
t169 = t26 / 0.2e1;
t167 = pkin(5) * t91;
t100 = Ifges(7,5) * t128 + Ifges(7,6) * t131;
t166 = t100 / 0.2e1;
t165 = -t128 / 0.2e1;
t164 = t128 / 0.2e1;
t163 = t131 / 0.2e1;
t161 = pkin(1) * t127;
t108 = t132 * t161;
t84 = -pkin(8) * t146 + t108;
t160 = t84 * mrSges(3,1);
t85 = pkin(8) * t145 + t130 * t161;
t159 = t85 * mrSges(3,2);
t158 = t91 * t55;
t11 = -mrSges(7,1) * t25 + mrSges(7,2) * t26;
t29 = mrSges(6,1) * t79 - mrSges(6,3) * t39;
t156 = t11 - t29;
t69 = t127 * pkin(2) + t108 + (-pkin(8) - qJ(3)) * t146;
t74 = qJ(3) * t145 + t85;
t47 = t123 * t69 + t126 * t74;
t42 = qJ(4) * t127 + t47;
t48 = t79 * pkin(3) - t80 * qJ(4) + t94;
t21 = -t122 * t42 + t125 * t48;
t13 = pkin(4) * t79 - pkin(9) * t68 + t21;
t22 = t122 * t48 + t125 * t42;
t19 = pkin(9) * t67 + t22;
t6 = t129 * t13 + t162 * t19;
t155 = Ifges(6,5) * t93 - Ifges(6,6) * t91;
t154 = Ifges(7,4) * t128;
t153 = Ifges(7,4) * t131;
t144 = t122 ^ 2 + t125 ^ 2;
t143 = t128 ^ 2 + t131 ^ 2;
t7 = Ifges(7,5) * t26 + Ifges(7,6) * t25 + Ifges(7,3) * t38;
t142 = Ifges(6,5) * t39 - Ifges(6,6) * t38 + Ifges(6,3) * t79;
t112 = -pkin(2) * t126 - pkin(3);
t41 = -t67 * mrSges(5,1) + t68 * mrSges(5,2);
t20 = t38 * mrSges(6,1) + t39 * mrSges(6,2);
t63 = t91 * mrSges(6,1) + t93 * mrSges(6,2);
t46 = -t123 * t74 + t126 * t69;
t140 = t93 * t143;
t96 = -t125 * mrSges(5,1) + t122 * mrSges(5,2);
t43 = -pkin(3) * t127 - t46;
t27 = -pkin(4) * t67 + t43;
t10 = pkin(5) * t38 - pkin(10) * t39 + t27;
t4 = pkin(10) * t79 + t6;
t1 = t10 * t131 - t128 * t4;
t2 = t10 * t128 + t131 * t4;
t139 = -t1 * t128 + t131 * t2;
t138 = mrSges(7,1) * t128 + mrSges(7,2) * t131;
t137 = -t122 * t21 + t125 * t22;
t95 = -pkin(4) * t125 + t112;
t54 = -pkin(10) * t93 + t167 + t95;
t57 = -t129 * t141 + t162 * t83;
t30 = -t128 * t57 + t131 * t54;
t31 = t128 * t54 + t131 * t57;
t136 = -t128 * t30 + t131 * t31;
t51 = Ifges(7,5) * t147 - Ifges(7,6) * t150 + Ifges(7,3) * t91;
t135 = Ifges(3,5) * t146 + Ifges(3,6) * t145 + t127 * t178 + t177;
t5 = -t129 * t19 + t13 * t162;
t102 = Ifges(7,1) * t128 + t153;
t101 = Ifges(7,2) * t131 + t154;
t99 = -mrSges(7,1) * t131 + mrSges(7,2) * t128;
t98 = Ifges(5,1) * t122 + Ifges(5,4) * t125;
t97 = Ifges(5,4) * t122 + Ifges(5,2) * t125;
t90 = t93 ^ 2;
t75 = t80 * mrSges(4,2);
t71 = mrSges(4,1) * t127 - mrSges(4,3) * t80;
t70 = -mrSges(4,2) * t127 - mrSges(4,3) * t79;
t65 = Ifges(6,1) * t93 - Ifges(6,4) * t91;
t64 = Ifges(6,4) * t93 - Ifges(6,2) * t91;
t58 = t138 * t93;
t53 = Ifges(7,5) * t91 + (Ifges(7,1) * t131 - t154) * t93;
t52 = Ifges(7,6) * t91 + (-Ifges(7,2) * t128 + t153) * t93;
t50 = mrSges(5,1) * t79 - mrSges(5,3) * t68;
t49 = -mrSges(5,2) * t79 + mrSges(5,3) * t67;
t33 = Ifges(5,1) * t68 + Ifges(5,4) * t67 + Ifges(5,5) * t79;
t32 = Ifges(5,4) * t68 + Ifges(5,2) * t67 + Ifges(5,6) * t79;
t28 = -mrSges(6,2) * t79 - mrSges(6,3) * t38;
t17 = Ifges(6,1) * t39 - Ifges(6,4) * t38 + Ifges(6,5) * t79;
t16 = Ifges(6,4) * t39 - Ifges(6,2) * t38 + Ifges(6,6) * t79;
t9 = Ifges(7,1) * t26 + Ifges(7,4) * t25 + Ifges(7,5) * t38;
t8 = Ifges(7,4) * t26 + Ifges(7,2) * t25 + Ifges(7,6) * t38;
t3 = -t79 * pkin(5) - t5;
t12 = [(mrSges(4,1) * t171 - 0.2e1 * Ifges(4,4) * t80 + Ifges(5,5) * t68 + Ifges(5,6) * t67 + (Ifges(5,3) + Ifges(4,2)) * t79 + t142) * t79 + ((Ifges(3,5) * t130 + Ifges(3,6) * t132) * t127 + 0.2e1 * (-t130 * t84 + t132 * t85) * mrSges(3,3) + (m(3) * pkin(1) ^ 2 + t130 * (Ifges(3,1) * t130 + Ifges(3,4) * t132) + t132 * (Ifges(3,4) * t130 + Ifges(3,2) * t132) - 0.2e1 * pkin(1) * (-mrSges(3,1) * t132 + mrSges(3,2) * t130)) * t124) * t124 + Ifges(4,1) * t80 ^ 2 + 0.2e1 * t43 * t41 + 0.2e1 * t22 * t49 + 0.2e1 * t21 * t50 + t39 * t17 + t25 * t8 + t26 * t9 + 0.2e1 * t27 * t20 + 0.2e1 * t6 * t28 + 0.2e1 * t5 * t29 + 0.2e1 * t2 * t14 + 0.2e1 * t1 * t15 + 0.2e1 * t3 * t11 + m(3) * (t84 ^ 2 + t85 ^ 2) + m(4) * (t46 ^ 2 + t47 ^ 2 + t94 ^ 2) + m(5) * (t21 ^ 2 + t22 ^ 2 + t43 ^ 2) + m(6) * (t27 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + (t7 - t16) * t38 + t75 * t171 + (t135 - 0.2e1 * t159 + 0.2e1 * t160 + t177) * t127 + t67 * t32 + t68 * t33 + 0.2e1 * t47 * t70 + 0.2e1 * t46 * t71 + Ifges(2,3); t57 * t28 + t3 * t58 + t2 * t60 + t1 * t61 + t27 * t63 + t46 * mrSges(4,1) - t47 * mrSges(4,2) + (m(4) * (t123 * t47 + t126 * t46) + t123 * t70 + t126 * t71) * pkin(2) + t30 * t15 + t31 * t14 + (-t122 * t50 + t125 * t49) * t110 + m(6) * (t27 * t95 - t5 * t55 + t57 * t6) + m(7) * (t1 * t30 + t2 * t31 + t3 * t55) + (-t6 * mrSges(6,3) + t7 / 0.2e1 - t16 / 0.2e1) * t91 + t53 * t169 + t52 * t170 + t160 - t159 + (t51 / 0.2e1 - t64 / 0.2e1) * t38 + t135 + (Ifges(5,5) * t122 + Ifges(5,6) * t125 + t155) * t79 / 0.2e1 + (t9 * t163 + t8 * t165 - t5 * mrSges(6,3) + t17 / 0.2e1) * t93 + t156 * t55 + m(5) * (t110 * t137 + t112 * t43) + t137 * mrSges(5,3) + t39 * t65 / 0.2e1 + t95 * t20 + t43 * t96 + t67 * t97 / 0.2e1 + t68 * t98 / 0.2e1 + t112 * t41 + t122 * t33 / 0.2e1 + t125 * t32 / 0.2e1; 0.2e1 * t112 * t96 + t122 * t98 + t125 * t97 + 0.2e1 * t30 * t61 + 0.2e1 * t31 * t60 + t58 * t172 + 0.2e1 * t95 * t63 + (-0.2e1 * mrSges(6,3) * t57 + t51 - t64) * t91 + (mrSges(6,3) * t172 - t128 * t52 + t131 * t53 + t65) * t93 + m(7) * (t30 ^ 2 + t31 ^ 2 + t174) + m(6) * (t57 ^ 2 + t95 ^ 2 + t174) + m(5) * (t110 ^ 2 * t144 + t112 ^ 2) + m(4) * (t123 ^ 2 + t126 ^ 2) * pkin(2) ^ 2 + 0.2e1 * (mrSges(4,1) * t126 - mrSges(4,2) * t123) * pkin(2) + 0.2e1 * t144 * t110 * mrSges(5,3) + t178; t79 * mrSges(4,1) + t122 * t49 + t125 * t50 + t75 + t156 * t91 + (t28 + t175) * t93 + m(7) * (t139 * t93 + t91 * t3) + m(6) * (-t5 * t91 + t6 * t93) + m(5) * (t122 * t22 + t125 * t21) + m(4) * t94; t91 * t58 + t176 * t93 + m(7) * (t136 * t93 + t158) + m(6) * (t57 * t93 + t158); m(4) + m(5) * t144 + m(6) * (t90 + t173) + m(7) * (t143 * t90 + t173); t128 * t14 + t131 * t15 + m(7) * (t1 * t131 + t128 * t2) + m(6) * t27 + m(5) * t43 + t20 + t41; t128 * t60 + t131 * t61 + m(7) * (t128 * t31 + t131 * t30) + m(6) * t95 + m(5) * t112 + t96 + t63; 0; m(7) * t143 + m(5) + m(6); t5 * mrSges(6,1) - t6 * mrSges(6,2) + t139 * mrSges(7,3) + t101 * t170 + t102 * t169 + t8 * t163 + t9 * t164 + t38 * t166 + t3 * t99 + t142 + (-m(7) * t3 - t11) * pkin(5) + (m(7) * t139 + t175) * pkin(10); t91 * t166 + t52 * t163 + t53 * t164 - pkin(5) * t58 - t57 * mrSges(6,2) + (t101 * t165 + t102 * t163) * t93 + t136 * mrSges(7,3) + t155 + (-m(7) * pkin(5) - mrSges(6,1) + t99) * t55 + (m(7) * t136 + t176) * pkin(10); t91 * t99 + m(7) * (pkin(10) * t140 - t167) + mrSges(7,3) * t140 - t63; 0; Ifges(6,3) - 0.2e1 * pkin(5) * t99 + t131 * t101 + t128 * t102 + m(7) * (pkin(10) ^ 2 * t143 + pkin(5) ^ 2) + 0.2e1 * t143 * pkin(10) * mrSges(7,3); mrSges(7,1) * t1 - mrSges(7,2) * t2 + t7; mrSges(7,1) * t30 - mrSges(7,2) * t31 + t51; -t58; -t99; -pkin(10) * t138 + t100; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t12(1) t12(2) t12(4) t12(7) t12(11) t12(16); t12(2) t12(3) t12(5) t12(8) t12(12) t12(17); t12(4) t12(5) t12(6) t12(9) t12(13) t12(18); t12(7) t12(8) t12(9) t12(10) t12(14) t12(19); t12(11) t12(12) t12(13) t12(14) t12(15) t12(20); t12(16) t12(17) t12(18) t12(19) t12(20) t12(21);];
Mq  = res;
