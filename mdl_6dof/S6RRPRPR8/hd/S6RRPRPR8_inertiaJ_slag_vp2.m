% Calculate joint inertia matrix for
% S6RRPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2018-11-23 17:05
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR8_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR8_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR8_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:05:15
% EndTime: 2018-11-23 17:05:16
% DurationCPUTime: 1.19s
% Computational Cost: add. (1764->306), mult. (3509->413), div. (0->0), fcn. (3580->8), ass. (0->112)
t149 = 2 * pkin(7);
t116 = sin(qJ(2));
t112 = sin(pkin(10));
t113 = cos(pkin(10));
t115 = sin(qJ(4));
t145 = cos(qJ(4));
t87 = t112 * t145 + t113 * t115;
t72 = t87 * t116;
t123 = -t112 * t115 + t113 * t145;
t73 = t123 * t116;
t90 = (pkin(3) * t112 + pkin(7)) * t116;
t19 = pkin(4) * t72 - qJ(5) * t73 + t90;
t148 = (-Ifges(6,4) - Ifges(5,5)) * t73 + (Ifges(5,6) - Ifges(6,6)) * t72;
t119 = -pkin(4) - pkin(5);
t147 = -t112 / 0.2e1;
t146 = t113 / 0.2e1;
t144 = pkin(7) * t116;
t118 = cos(qJ(2));
t143 = pkin(7) * t118;
t114 = sin(qJ(6));
t117 = cos(qJ(6));
t91 = -qJ(5) * t114 + t117 * t119;
t142 = t91 * mrSges(7,1);
t92 = qJ(5) * t117 + t114 * t119;
t141 = t92 * mrSges(7,2);
t139 = Ifges(6,2) + Ifges(5,3);
t138 = pkin(8) + qJ(3);
t133 = t113 * t116;
t93 = -pkin(2) * t118 - qJ(3) * t116 - pkin(1);
t82 = t113 * t93;
t41 = -pkin(8) * t133 + t82 + (-pkin(7) * t112 - pkin(3)) * t118;
t134 = t112 * t116;
t60 = t112 * t93 + t113 * t143;
t50 = -pkin(8) * t134 + t60;
t18 = t115 * t41 + t145 * t50;
t128 = t138 * t112;
t95 = t138 * t113;
t54 = -t115 * t128 + t145 * t95;
t74 = mrSges(4,1) * t134 + mrSges(4,2) * t133;
t136 = Ifges(4,4) * t112;
t135 = Ifges(4,4) * t113;
t58 = mrSges(6,1) * t118 + mrSges(6,2) * t73;
t132 = t112 ^ 2 + t113 ^ 2;
t52 = t115 * t95 + t128 * t145;
t131 = t52 ^ 2 + t54 ^ 2;
t29 = -t114 * t73 + t117 * t72;
t30 = t114 * t72 + t117 * t73;
t130 = Ifges(7,5) * t30 + Ifges(7,6) * t29 + Ifges(7,3) * t118;
t102 = -pkin(3) * t113 - pkin(2);
t34 = t72 * mrSges(5,1) + mrSges(5,2) * t73;
t45 = -mrSges(5,1) * t123 + mrSges(5,2) * t87;
t33 = mrSges(6,1) * t72 - t73 * mrSges(6,3);
t44 = -mrSges(6,1) * t123 - t87 * mrSges(6,3);
t94 = -t113 * mrSges(4,1) + mrSges(4,2) * t112;
t10 = -qJ(5) * t118 + t18;
t17 = -t115 * t50 + t145 * t41;
t7 = -t29 * mrSges(7,1) + t30 * mrSges(7,2);
t39 = -t114 * t87 - t117 * t123;
t40 = -t114 * t123 + t117 * t87;
t14 = -t39 * mrSges(7,1) + t40 * mrSges(7,2);
t127 = t117 * mrSges(7,1) - t114 * mrSges(7,2);
t59 = -t112 * t143 + t82;
t126 = -t112 * t59 + t113 * t60;
t11 = pkin(4) * t118 - t17;
t125 = qJ(5) * t87 - t102;
t37 = Ifges(7,6) * t39;
t38 = Ifges(7,5) * t40;
t31 = -pkin(9) * t87 + t52;
t32 = -pkin(9) * t123 + t54;
t8 = -t114 * t32 + t117 * t31;
t9 = t114 * t31 + t117 * t32;
t124 = t8 * mrSges(7,1) - t9 * mrSges(7,2) + t37 + t38;
t3 = pkin(5) * t118 - pkin(9) * t73 + t11;
t6 = pkin(9) * t72 + t10;
t1 = -t114 * t6 + t117 * t3;
t2 = t114 * t3 + t117 * t6;
t122 = t1 * mrSges(7,1) - t2 * mrSges(7,2) + t130;
t121 = pkin(7) ^ 2;
t111 = t118 ^ 2;
t110 = t116 ^ 2;
t106 = t110 * t121;
t97 = Ifges(4,1) * t112 + t135;
t96 = Ifges(4,2) * t113 + t136;
t89 = -mrSges(4,1) * t118 - mrSges(4,3) * t133;
t88 = mrSges(4,2) * t118 - mrSges(4,3) * t134;
t80 = Ifges(6,4) * t87;
t79 = Ifges(5,5) * t87;
t78 = Ifges(5,6) * t123;
t77 = Ifges(6,6) * t123;
t71 = -t118 * Ifges(4,5) + (Ifges(4,1) * t113 - t136) * t116;
t70 = -t118 * Ifges(4,6) + (-Ifges(4,2) * t112 + t135) * t116;
t57 = -mrSges(5,1) * t118 - mrSges(5,3) * t73;
t56 = mrSges(5,2) * t118 - mrSges(5,3) * t72;
t55 = -mrSges(6,2) * t72 - mrSges(6,3) * t118;
t49 = Ifges(5,1) * t87 + Ifges(5,4) * t123;
t48 = Ifges(6,1) * t87 - Ifges(6,5) * t123;
t47 = Ifges(5,4) * t87 + Ifges(5,2) * t123;
t46 = Ifges(6,5) * t87 - Ifges(6,3) * t123;
t36 = -pkin(4) * t123 - t125;
t28 = Ifges(5,1) * t73 - Ifges(5,4) * t72 - Ifges(5,5) * t118;
t27 = Ifges(6,1) * t73 - Ifges(6,4) * t118 + Ifges(6,5) * t72;
t26 = Ifges(5,4) * t73 - Ifges(5,2) * t72 - Ifges(5,6) * t118;
t25 = Ifges(6,5) * t73 - Ifges(6,6) * t118 + Ifges(6,3) * t72;
t22 = -t119 * t123 + t125;
t21 = mrSges(7,1) * t118 - mrSges(7,3) * t30;
t20 = -mrSges(7,2) * t118 + mrSges(7,3) * t29;
t16 = Ifges(7,1) * t40 + Ifges(7,4) * t39;
t15 = Ifges(7,4) * t40 + Ifges(7,2) * t39;
t12 = pkin(5) * t72 + t19;
t5 = Ifges(7,1) * t30 + Ifges(7,4) * t29 + Ifges(7,5) * t118;
t4 = Ifges(7,4) * t30 + Ifges(7,2) * t29 + Ifges(7,6) * t118;
t13 = [0.2e1 * t1 * t21 + 0.2e1 * t10 * t55 + 0.2e1 * t11 * t58 - 0.2e1 * t12 * t7 + 0.2e1 * t17 * t57 + 0.2e1 * t18 * t56 + 0.2e1 * t19 * t33 + 0.2e1 * t2 * t20 + t29 * t4 + t30 * t5 + 0.2e1 * t90 * t34 + 0.2e1 * t59 * t89 + 0.2e1 * t60 * t88 + Ifges(2,3) + (t27 + t28) * t73 + (t25 - t26) * t72 + (t110 + t111) * mrSges(3,3) * t149 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t116 - t112 * t70 + t113 * t71 + t74 * t149) * t116 + m(3) * (pkin(1) ^ 2 + t111 * t121 + t106) + m(4) * (t59 ^ 2 + t60 ^ 2 + t106) + m(5) * (t17 ^ 2 + t18 ^ 2 + t90 ^ 2) + m(7) * (t1 ^ 2 + t12 ^ 2 + t2 ^ 2) + m(6) * (t10 ^ 2 + t11 ^ 2 + t19 ^ 2) + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(4,3) + Ifges(3,2) + t139) * t118 + (-Ifges(4,5) * t113 + Ifges(4,6) * t112 + (2 * Ifges(3,4))) * t116 + t130 + t148) * t118; t70 * t146 + (t55 + t56) * t54 + (-t57 + t58) * t52 + (-t112 * t89 + t113 * t88) * qJ(3) + m(5) * (t102 * t90 - t17 * t52 + t18 * t54) + m(6) * (t10 * t54 + t11 * t52 + t19 * t36) + m(7) * (t1 * t8 - t12 * t22 + t2 * t9) + (t27 / 0.2e1 + t28 / 0.2e1 + t11 * mrSges(6,2) - t17 * mrSges(5,3)) * t87 + (t48 / 0.2e1 + t49 / 0.2e1) * t73 + (t46 / 0.2e1 - t47 / 0.2e1) * t72 + t112 * t71 / 0.2e1 + t90 * t45 + t102 * t34 - pkin(2) * t74 + t39 * t4 / 0.2e1 + t40 * t5 / 0.2e1 + t19 * t44 + t29 * t15 / 0.2e1 + t30 * t16 / 0.2e1 + t36 * t33 + t9 * t20 + t8 * t21 + t22 * t7 - t12 * t14 + (Ifges(3,5) + t96 * t147 + t97 * t146 + (t94 - mrSges(3,1)) * pkin(7)) * t116 + m(4) * (-pkin(2) * t144 + qJ(3) * t126) + t126 * mrSges(4,3) + (-t1 * t40 + t2 * t39) * mrSges(7,3) - (t25 / 0.2e1 - t26 / 0.2e1 - t18 * mrSges(5,3) - t10 * mrSges(6,2)) * t123 + (Ifges(4,5) * t147 - Ifges(4,6) * t113 / 0.2e1 - t80 / 0.2e1 + t77 / 0.2e1 - t79 / 0.2e1 - t78 / 0.2e1 + t38 / 0.2e1 + t37 / 0.2e1 + Ifges(3,6) - pkin(7) * mrSges(3,2)) * t118; -0.2e1 * pkin(2) * t94 + 0.2e1 * t102 * t45 + t112 * t97 + t113 * t96 + 0.2e1 * t22 * t14 + t39 * t15 + t40 * t16 + 0.2e1 * t36 * t44 + Ifges(3,3) + (t48 + t49) * t87 - (t46 - t47) * t123 + m(4) * (qJ(3) ^ 2 * t132 + pkin(2) ^ 2) + m(5) * (t102 ^ 2 + t131) + m(6) * (t36 ^ 2 + t131) + m(7) * (t22 ^ 2 + t8 ^ 2 + t9 ^ 2) + 0.2e1 * (t39 * t9 - t40 * t8) * mrSges(7,3) + 0.2e1 * t132 * qJ(3) * mrSges(4,3) + 0.2e1 * (t123 * t54 + t52 * t87) * (mrSges(6,2) + mrSges(5,3)); m(4) * t144 + m(5) * t90 + m(6) * t19 + m(7) * t12 + t33 + t34 - t7 + t74; -m(4) * pkin(2) + m(5) * t102 + m(6) * t36 - m(7) * t22 - t14 + t44 + t45 + t94; m(4) + m(5) + m(6) + m(7); m(7) * (t1 * t91 + t2 * t92) + m(6) * (-pkin(4) * t11 + qJ(5) * t10) - t139 * t118 + t91 * t21 + t92 * t20 + qJ(5) * t55 - pkin(4) * t58 + t10 * mrSges(6,3) - t11 * mrSges(6,1) + t17 * mrSges(5,1) - t18 * mrSges(5,2) - t122 - t148; -t77 + t78 + t79 + t80 + (-mrSges(5,2) + mrSges(6,3)) * t54 + (-mrSges(5,1) - mrSges(6,1)) * t52 + (t39 * t92 - t40 * t91) * mrSges(7,3) + (-pkin(4) * t87 + qJ(5) * t123) * mrSges(6,2) + m(7) * (t8 * t91 + t9 * t92) + m(6) * (-pkin(4) * t52 + qJ(5) * t54) - t124; 0; 0.2e1 * pkin(4) * mrSges(6,1) - 0.2e1 * t142 + 0.2e1 * t141 + 0.2e1 * qJ(5) * mrSges(6,3) + Ifges(7,3) + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2) + m(7) * (t91 ^ 2 + t92 ^ 2) + t139; t114 * t20 + t117 * t21 + m(7) * (t1 * t117 + t114 * t2) + m(6) * t11 + t58; t87 * mrSges(6,2) + (t114 * t39 - t117 * t40) * mrSges(7,3) + m(7) * (t114 * t9 + t117 * t8) + m(6) * t52; 0; -mrSges(6,1) - m(6) * pkin(4) + m(7) * (t114 * t92 + t117 * t91) - t127; m(6) + m(7) * (t114 ^ 2 + t117 ^ 2); t122; t124; 0; -Ifges(7,3) - t141 + t142; t127; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t13(1) t13(2) t13(4) t13(7) t13(11) t13(16); t13(2) t13(3) t13(5) t13(8) t13(12) t13(17); t13(4) t13(5) t13(6) t13(9) t13(13) t13(18); t13(7) t13(8) t13(9) t13(10) t13(14) t13(19); t13(11) t13(12) t13(13) t13(14) t13(15) t13(20); t13(16) t13(17) t13(18) t13(19) t13(20) t13(21);];
Mq  = res;
