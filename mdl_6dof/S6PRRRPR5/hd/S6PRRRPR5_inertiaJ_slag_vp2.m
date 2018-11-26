% Calculate joint inertia matrix for
% S6PRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2018-11-23 15:25
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_inertiaJ_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:24:46
% EndTime: 2018-11-23 15:24:47
% DurationCPUTime: 1.40s
% Computational Cost: add. (2599->352), mult. (6220->524), div. (0->0), fcn. (7054->14), ass. (0->129)
t124 = sin(pkin(7));
t135 = cos(qJ(3));
t151 = t124 * t135;
t130 = sin(qJ(4));
t134 = cos(qJ(4));
t123 = sin(pkin(13));
t126 = cos(pkin(13));
t88 = t123 * t130 - t126 * t134;
t89 = t123 * t134 + t126 * t130;
t175 = Ifges(5,5) * t130 + Ifges(6,5) * t89 + Ifges(5,6) * t134 - Ifges(6,6) * t88;
t125 = sin(pkin(6));
t128 = cos(pkin(6));
t131 = sin(qJ(3));
t132 = sin(qJ(2));
t127 = cos(pkin(7));
t136 = cos(qJ(2));
t150 = t127 * t136;
t152 = t124 * t131;
t59 = t128 * t152 + (t131 * t150 + t132 * t135) * t125;
t78 = -t124 * t125 * t136 + t127 * t128;
t29 = -t130 * t59 + t134 * t78;
t30 = t130 * t78 + t134 * t59;
t13 = t123 * t30 - t126 * t29;
t174 = t13 ^ 2;
t57 = -t128 * t151 + (t131 * t132 - t135 * t150) * t125;
t56 = t57 ^ 2;
t161 = -qJ(5) - pkin(10);
t145 = t161 * t130;
t96 = t161 * t134;
t63 = -t123 * t96 - t126 * t145;
t173 = t63 ^ 2;
t172 = 0.2e1 * t63;
t129 = sin(qJ(6));
t133 = cos(qJ(6));
t79 = t127 * t134 - t130 * t152;
t80 = t127 * t130 + t134 * t152;
t49 = t123 * t79 + t126 * t80;
t33 = -t129 * t49 - t133 * t151;
t171 = t33 / 0.2e1;
t34 = -t129 * t151 + t133 * t49;
t170 = t34 / 0.2e1;
t97 = Ifges(7,5) * t129 + Ifges(7,6) * t133;
t169 = t97 / 0.2e1;
t168 = -t129 / 0.2e1;
t167 = t129 / 0.2e1;
t166 = t133 / 0.2e1;
t164 = pkin(2) * t135;
t163 = t13 * t63;
t162 = Ifges(5,3) + Ifges(6,3);
t16 = -mrSges(7,1) * t33 + mrSges(7,2) * t34;
t39 = -mrSges(6,1) * t151 - mrSges(6,3) * t49;
t160 = t16 - t39;
t83 = t127 * t131 * pkin(2) + pkin(9) * t151;
t72 = pkin(10) * t127 + t83;
t73 = (-pkin(3) * t135 - pkin(10) * t131 - pkin(2)) * t124;
t40 = -t130 * t72 + t134 * t73;
t23 = -pkin(4) * t151 - qJ(5) * t80 + t40;
t41 = t130 * t73 + t134 * t72;
t27 = qJ(5) * t79 + t41;
t11 = t123 * t23 + t126 * t27;
t158 = Ifges(7,4) * t129;
t157 = Ifges(7,4) * t133;
t156 = t129 * t89;
t155 = t133 * t89;
t110 = pkin(4) * t123 + pkin(11);
t154 = t110 * t129;
t153 = t110 * t133;
t148 = t129 ^ 2 + t133 ^ 2;
t147 = t130 ^ 2 + t134 ^ 2;
t48 = t123 * t80 - t126 * t79;
t5 = Ifges(7,5) * t34 + Ifges(7,6) * t33 + Ifges(7,3) * t48;
t146 = Ifges(4,5) * t152 + Ifges(4,6) * t151 + Ifges(4,3) * t127;
t112 = -pkin(4) * t134 - pkin(3);
t21 = t48 * mrSges(6,1) + t49 * mrSges(6,2);
t60 = t88 * mrSges(6,1) + t89 * mrSges(6,2);
t144 = -Ifges(5,5) * t80 - Ifges(6,5) * t49 - Ifges(5,6) * t79 + Ifges(6,6) * t48;
t105 = pkin(9) * t152;
t71 = t105 + (-pkin(3) - t164) * t127;
t50 = -pkin(4) * t79 + t71;
t12 = pkin(5) * t48 - pkin(11) * t49 + t50;
t9 = -pkin(11) * t151 + t11;
t1 = t12 * t133 - t129 * t9;
t2 = t12 * t129 + t133 * t9;
t143 = -t1 * t129 + t2 * t133;
t15 = t123 * t29 + t126 * t30;
t3 = -t129 * t15 + t133 * t57;
t4 = t129 * t57 + t133 * t15;
t142 = -t129 * t3 + t133 * t4;
t141 = mrSges(7,1) * t129 + mrSges(7,2) * t133;
t10 = -t123 * t27 + t126 * t23;
t53 = pkin(5) * t88 - pkin(11) * t89 + t112;
t65 = t123 * t145 - t126 * t96;
t24 = -t129 * t65 + t133 * t53;
t25 = t129 * t53 + t133 * t65;
t140 = -t24 * t129 + t25 * t133;
t139 = -t130 * t29 + t134 * t30;
t35 = Ifges(7,5) * t155 - Ifges(7,6) * t156 + Ifges(7,3) * t88;
t111 = -pkin(4) * t126 - pkin(5);
t101 = Ifges(5,1) * t130 + Ifges(5,4) * t134;
t100 = Ifges(7,1) * t129 + t157;
t99 = Ifges(5,4) * t130 + Ifges(5,2) * t134;
t98 = Ifges(7,2) * t133 + t158;
t95 = -mrSges(5,1) * t134 + mrSges(5,2) * t130;
t94 = -mrSges(7,1) * t133 + mrSges(7,2) * t129;
t92 = -mrSges(4,2) * t127 + mrSges(4,3) * t151;
t91 = mrSges(4,1) * t127 - mrSges(4,3) * t152;
t82 = t127 * t164 - t105;
t81 = (-mrSges(4,1) * t135 + mrSges(4,2) * t131) * t124;
t67 = -mrSges(5,1) * t151 - mrSges(5,3) * t80;
t66 = mrSges(5,2) * t151 + mrSges(5,3) * t79;
t62 = Ifges(6,1) * t89 - Ifges(6,4) * t88;
t61 = Ifges(6,4) * t89 - Ifges(6,2) * t88;
t55 = mrSges(7,1) * t88 - mrSges(7,3) * t155;
t54 = -mrSges(7,2) * t88 - mrSges(7,3) * t156;
t52 = -mrSges(5,1) * t79 + mrSges(5,2) * t80;
t51 = t141 * t89;
t43 = Ifges(5,1) * t80 + Ifges(5,4) * t79 - Ifges(5,5) * t151;
t42 = Ifges(5,4) * t80 + Ifges(5,2) * t79 - Ifges(5,6) * t151;
t38 = mrSges(6,2) * t151 - mrSges(6,3) * t48;
t37 = Ifges(7,5) * t88 + (Ifges(7,1) * t133 - t158) * t89;
t36 = Ifges(7,6) * t88 + (-Ifges(7,2) * t129 + t157) * t89;
t20 = Ifges(6,1) * t49 - Ifges(6,4) * t48 - Ifges(6,5) * t151;
t19 = Ifges(6,4) * t49 - Ifges(6,2) * t48 - Ifges(6,6) * t151;
t18 = mrSges(7,1) * t48 - mrSges(7,3) * t34;
t17 = -mrSges(7,2) * t48 + mrSges(7,3) * t33;
t8 = pkin(5) * t151 - t10;
t7 = Ifges(7,1) * t34 + Ifges(7,4) * t33 + Ifges(7,5) * t48;
t6 = Ifges(7,4) * t34 + Ifges(7,2) * t33 + Ifges(7,6) * t48;
t14 = [m(2) + m(6) * (t15 ^ 2 + t174 + t56) + m(7) * (t3 ^ 2 + t4 ^ 2 + t174) + m(5) * (t29 ^ 2 + t30 ^ 2 + t56) + m(4) * (t59 ^ 2 + t78 ^ 2 + t56) + m(3) * (t128 ^ 2 + (t132 ^ 2 + t136 ^ 2) * t125 ^ 2); t15 * t38 + t4 * t17 + t3 * t18 + t29 * t67 + t30 * t66 + t59 * t92 + t78 * t81 + t160 * t13 + (mrSges(3,1) * t136 - mrSges(3,2) * t132) * t125 + (t21 + t52 - t91) * t57 + m(6) * (-t10 * t13 + t11 * t15 + t50 * t57) + m(7) * (t1 * t3 + t13 * t8 + t2 * t4) + m(5) * (t29 * t40 + t30 * t41 + t57 * t71) + m(4) * (-pkin(2) * t124 * t78 - t57 * t82 + t59 * t83); t127 * t146 + 0.2e1 * t82 * t91 + 0.2e1 * t83 * t92 + 0.2e1 * t71 * t52 + t79 * t42 + t80 * t43 + 0.2e1 * t41 * t66 + 0.2e1 * t40 * t67 + t49 * t20 + 0.2e1 * t50 * t21 + 0.2e1 * t10 * t39 + t33 * t6 + t34 * t7 + 0.2e1 * t11 * t38 + 0.2e1 * t8 * t16 + 0.2e1 * t2 * t17 + 0.2e1 * t1 * t18 + Ifges(3,3) + (-t19 + t5) * t48 + (-0.2e1 * pkin(2) * t81 + (Ifges(4,1) * t152 + Ifges(4,5) * t127) * t131 + (Ifges(4,6) * t127 + (0.2e1 * Ifges(4,4) * t131 + (Ifges(4,2) + t162) * t135) * t124 + t144) * t135) * t124 + m(4) * (pkin(2) ^ 2 * t124 ^ 2 + t82 ^ 2 + t83 ^ 2) + m(5) * (t40 ^ 2 + t41 ^ 2 + t71 ^ 2) + m(6) * (t10 ^ 2 + t11 ^ 2 + t50 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t8 ^ 2); -t59 * mrSges(4,2) + t13 * t51 + t3 * t55 + t4 * t54 + (t13 * t89 - t15 * t88) * mrSges(6,3) + t139 * mrSges(5,3) + (-mrSges(4,1) + t60 + t95) * t57 + m(6) * (t112 * t57 + t15 * t65 + t163) + m(7) * (t24 * t3 + t25 * t4 + t163) + m(5) * (-pkin(3) * t57 + t139 * pkin(10)); -t175 * t151 / 0.2e1 + (pkin(10) * t66 + t41 * mrSges(5,3) + t42 / 0.2e1) * t134 + (-pkin(10) * t67 - t40 * mrSges(5,3) + t43 / 0.2e1) * t130 + t71 * t95 + t79 * t99 / 0.2e1 + t80 * t101 / 0.2e1 + t112 * t21 + t82 * mrSges(4,1) - t83 * mrSges(4,2) + t50 * t60 + t49 * t62 / 0.2e1 + t65 * t38 + t8 * t51 - pkin(3) * t52 + t2 * t54 + t1 * t55 + t24 * t18 + t25 * t17 + m(6) * (-t10 * t63 + t11 * t65 + t112 * t50) + m(7) * (t1 * t24 + t2 * t25 + t63 * t8) + (t5 / 0.2e1 - t19 / 0.2e1 - t11 * mrSges(6,3)) * t88 + (-t61 / 0.2e1 + t35 / 0.2e1) * t48 + m(5) * (-pkin(3) * t71 + (-t40 * t130 + t41 * t134) * pkin(10)) + (t7 * t166 + t6 * t168 + t20 / 0.2e1 - t10 * mrSges(6,3)) * t89 + t160 * t63 + t146 + t37 * t170 + t36 * t171; -0.2e1 * pkin(3) * t95 + t130 * t101 + 0.2e1 * t112 * t60 + t134 * t99 + 0.2e1 * t24 * t55 + 0.2e1 * t25 * t54 + t51 * t172 + Ifges(4,3) + 0.2e1 * t147 * pkin(10) * mrSges(5,3) + (-0.2e1 * t65 * mrSges(6,3) + t35 - t61) * t88 + m(7) * (t24 ^ 2 + t25 ^ 2 + t173) + m(6) * (t112 ^ 2 + t65 ^ 2 + t173) + m(5) * (t147 * pkin(10) ^ 2 + pkin(3) ^ 2) + (mrSges(6,3) * t172 - t129 * t36 + t133 * t37 + t62) * t89; t29 * mrSges(5,1) - t30 * mrSges(5,2) - t15 * mrSges(6,2) + (-mrSges(6,1) + t94) * t13 + t142 * mrSges(7,3) + m(7) * (t142 * t110 + t111 * t13) + m(6) * (t123 * t15 - t126 * t13) * pkin(4); t17 * t153 - t18 * t154 + m(7) * (t143 * t110 + t111 * t8) + t143 * mrSges(7,3) + t6 * t166 + t7 * t167 + t8 * t94 + t48 * t169 + t98 * t171 + t100 * t170 + t111 * t16 + t40 * mrSges(5,1) - t41 * mrSges(5,2) + t10 * mrSges(6,1) - t11 * mrSges(6,2) + (m(6) * (t10 * t126 + t11 * t123) + t126 * t39 + t123 * t38) * pkin(4) - t144 - t162 * t151; -t55 * t154 + m(7) * (t140 * t110 + t111 * t63) + t37 * t167 + t36 * t166 + t111 * t51 - t63 * mrSges(6,1) - t65 * mrSges(6,2) + t63 * t94 + t88 * t169 + t54 * t153 + (t100 * t166 + t98 * t168) * t89 + (-t130 * mrSges(5,1) - t134 * mrSges(5,2)) * pkin(10) + t140 * mrSges(7,3) + (m(6) * (t123 * t65 - t126 * t63) + (-t123 * t88 - t126 * t89) * mrSges(6,3)) * pkin(4) + t175; t129 * t100 + 0.2e1 * t111 * t94 + t133 * t98 + m(7) * (t148 * t110 ^ 2 + t111 ^ 2) + m(6) * (t123 ^ 2 + t126 ^ 2) * pkin(4) ^ 2 + t162 + 0.2e1 * (mrSges(6,1) * t126 - mrSges(6,2) * t123) * pkin(4) + 0.2e1 * t148 * t110 * mrSges(7,3); m(6) * t57 + m(7) * (t129 * t4 + t133 * t3); t129 * t17 + t133 * t18 + m(7) * (t1 * t133 + t129 * t2) + m(6) * t50 + t21; t129 * t54 + t133 * t55 + m(7) * (t129 * t25 + t133 * t24) + m(6) * t112 + t60; 0; m(7) * t148 + m(6); mrSges(7,1) * t3 - mrSges(7,2) * t4; mrSges(7,1) * t1 - mrSges(7,2) * t2 + t5; mrSges(7,1) * t24 - mrSges(7,2) * t25 + t35; -t141 * t110 + t97; -t94; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t14(1) t14(2) t14(4) t14(7) t14(11) t14(16); t14(2) t14(3) t14(5) t14(8) t14(12) t14(17); t14(4) t14(5) t14(6) t14(9) t14(13) t14(18); t14(7) t14(8) t14(9) t14(10) t14(14) t14(19); t14(11) t14(12) t14(13) t14(14) t14(15) t14(20); t14(16) t14(17) t14(18) t14(19) t14(20) t14(21);];
Mq  = res;
