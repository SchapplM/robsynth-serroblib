% Calculate joint inertia matrix for
% S6RRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2018-11-23 17:24
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR6_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR6_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR6_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:24:26
% EndTime: 2018-11-23 17:24:28
% DurationCPUTime: 1.18s
% Computational Cost: add. (1900->252), mult. (3338->350), div. (0->0), fcn. (3477->8), ass. (0->101)
t90 = sin(qJ(2));
t94 = cos(qJ(2));
t156 = t90 ^ 2 + t94 ^ 2;
t87 = sin(qJ(6));
t136 = Ifges(7,4) * t87;
t91 = cos(qJ(6));
t68 = Ifges(7,2) * t91 + t136;
t135 = Ifges(7,4) * t91;
t69 = Ifges(7,1) * t87 + t135;
t117 = t91 * t68 + t87 * t69 + Ifges(6,3);
t155 = Ifges(5,3) + t117;
t85 = t91 ^ 2;
t125 = t85 * mrSges(7,3);
t83 = t87 ^ 2;
t126 = t83 * mrSges(7,3);
t88 = sin(qJ(5));
t89 = sin(qJ(4));
t92 = cos(qJ(5));
t93 = cos(qJ(4));
t55 = t88 * t89 - t92 * t93;
t59 = t88 * t93 + t89 * t92;
t66 = -t91 * mrSges(7,1) + mrSges(7,2) * t87;
t154 = (-mrSges(6,2) + t125 + t126) * t59 + (t66 - mrSges(6,1)) * t55;
t139 = pkin(7) - pkin(8);
t115 = t139 * t94;
t116 = t139 * t90;
t38 = -t89 * t115 + t93 * t116;
t58 = -t89 * t94 + t90 * t93;
t100 = -t58 * pkin(9) + t38;
t39 = t93 * t115 + t89 * t116;
t57 = -t89 * t90 - t93 * t94;
t22 = pkin(9) * t57 + t39;
t11 = -t92 * t100 + t22 * t88;
t29 = -t92 * t57 + t58 * t88;
t67 = Ifges(7,5) * t87 + Ifges(7,6) * t91;
t153 = t11 * t66 + t29 * t67 / 0.2e1;
t119 = t83 + t85;
t113 = mrSges(7,3) * t119;
t152 = t93 * mrSges(5,1) - t89 * mrSges(5,2) + t154;
t124 = t87 * mrSges(7,3);
t30 = t57 * t88 + t58 * t92;
t15 = -mrSges(7,2) * t29 - t30 * t124;
t131 = t30 * t91;
t16 = mrSges(7,1) * t29 - mrSges(7,3) * t131;
t151 = t91 * t15 - t87 * t16;
t150 = -m(4) * pkin(2) - mrSges(4,1);
t107 = mrSges(7,1) * t87 + mrSges(7,2) * t91;
t14 = t107 * t30;
t149 = m(7) * t11 + t14;
t13 = t88 * t100 + t92 * t22;
t65 = -t94 * pkin(2) - t90 * qJ(3) - pkin(1);
t45 = t94 * pkin(3) - t65;
t37 = -pkin(4) * t57 + t45;
t8 = pkin(5) * t29 - pkin(10) * t30 + t37;
t3 = t13 * t91 + t8 * t87;
t137 = t3 * t91;
t2 = -t13 * t87 + t8 * t91;
t109 = -t2 * t87 + t137;
t148 = m(7) * t109 + t151;
t132 = t30 * t87;
t6 = Ifges(7,6) * t29 + (-Ifges(7,2) * t87 + t135) * t30;
t7 = Ifges(7,5) * t29 + (Ifges(7,1) * t91 - t136) * t30;
t147 = -t68 * t132 / 0.2e1 + t69 * t131 / 0.2e1 + t87 * t7 / 0.2e1 + t91 * t6 / 0.2e1 - t2 * t124 + t153;
t146 = t11 ^ 2;
t145 = t55 ^ 2;
t144 = 0.2e1 * t11;
t143 = 0.2e1 * t37;
t142 = -0.2e1 * t65;
t138 = pkin(5) * t66;
t134 = t11 * t55;
t95 = -pkin(2) - pkin(3);
t63 = -qJ(3) * t89 + t93 * t95;
t62 = -pkin(4) + t63;
t64 = qJ(3) * t93 + t89 * t95;
t35 = t62 * t92 - t64 * t88;
t33 = pkin(5) - t35;
t21 = t33 * t66;
t36 = t88 * t62 + t92 * t64;
t130 = t36 * mrSges(6,2);
t129 = t63 * mrSges(5,1);
t128 = t64 * mrSges(5,2);
t73 = -pkin(4) * t92 - pkin(5);
t127 = t73 * t66;
t121 = Ifges(7,5) * t131 + Ifges(7,3) * t29;
t120 = t156 * pkin(7) ^ 2;
t114 = t119 * pkin(10);
t112 = t119 * t59;
t72 = pkin(4) * t88 + pkin(10);
t111 = t119 * t72;
t108 = t92 * mrSges(6,1) - t88 * mrSges(6,2);
t106 = 0.2e1 * t113;
t104 = -t11 * mrSges(6,1) - t13 * mrSges(6,2) + mrSges(7,3) * t137 + Ifges(6,5) * t30 - Ifges(6,6) * t29;
t103 = t108 * pkin(4);
t34 = -pkin(10) + t36;
t24 = t34 * t126;
t25 = t34 * t125;
t31 = t35 * mrSges(6,1);
t101 = t21 + t24 + t25 + t31 - t130 - t117;
t99 = t38 * mrSges(5,1) - t39 * mrSges(5,2) + Ifges(5,5) * t58 + Ifges(5,6) * t57 + t104;
t52 = t59 ^ 2;
t1 = [0.2e1 * t45 * (-t57 * mrSges(5,1) + t58 * mrSges(5,2)) + t58 * (Ifges(5,1) * t58 + Ifges(5,4) * t57) + t57 * (Ifges(5,4) * t58 + Ifges(5,2) * t57) + t14 * t144 + 0.2e1 * t3 * t15 + 0.2e1 * t2 * t16 + Ifges(2,3) + (0.2e1 * pkin(1) * mrSges(3,1) + mrSges(4,1) * t142 + (Ifges(3,2) + Ifges(4,3)) * t94) * t94 + (mrSges(6,1) * t143 - 0.2e1 * t13 * mrSges(6,3) + Ifges(6,2) * t29 + t121) * t29 + (mrSges(6,2) * t143 + mrSges(6,3) * t144 + Ifges(6,1) * t30 - t87 * t6 + t91 * t7 + (-Ifges(7,6) * t87 - (2 * Ifges(6,4))) * t29) * t30 + m(4) * (t65 ^ 2 + t120) + m(3) * (pkin(1) ^ 2 + t120) + m(6) * (t13 ^ 2 + t37 ^ 2 + t146) + m(5) * (t38 ^ 2 + t39 ^ 2 + t45 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t146) + 0.2e1 * (-t38 * t58 + t39 * t57) * mrSges(5,3) + 0.2e1 * (mrSges(4,2) + mrSges(3,3)) * pkin(7) * t156 + (-0.2e1 * pkin(1) * mrSges(3,2) + mrSges(4,3) * t142 + 0.2e1 * (Ifges(3,4) - Ifges(4,5)) * t94 + (Ifges(3,1) + Ifges(4,1)) * t90) * t90; (-pkin(2) * mrSges(4,2) + Ifges(4,4) + Ifges(3,5)) * t90 + m(5) * (t38 * t63 + t39 * t64) + m(6) * (-t11 * t35 + t13 * t36) + (-t36 * t29 - t35 * t30) * mrSges(6,3) + (t64 * t57 - t63 * t58) * mrSges(5,3) + ((m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3)) * t94 + (-mrSges(3,1) + t150) * t90) * pkin(7) + m(7) * (t109 * t34 + t11 * t33) - t99 + t33 * t14 + (-t30 * t69 / 0.2e1 + t34 * t15 - t6 / 0.2e1) * t91 + (t30 * t68 / 0.2e1 + t2 * mrSges(7,3) - t34 * t16 - t7 / 0.2e1) * t87 + (qJ(3) * mrSges(4,2) + Ifges(3,6) - Ifges(4,6)) * t94 - t153; 0.2e1 * pkin(2) * mrSges(4,1) - 0.2e1 * t129 + 0.2e1 * t128 + 0.2e1 * t130 + 0.2e1 * qJ(3) * mrSges(4,3) - 0.2e1 * t21 + Ifges(4,2) + Ifges(3,3) - 0.2e1 * t24 - 0.2e1 * t25 - 0.2e1 * t31 + m(7) * (t119 * t34 ^ 2 + t33 ^ 2) + m(6) * (t35 ^ 2 + t36 ^ 2) + m(5) * (t63 ^ 2 + t64 ^ 2) + m(4) * (pkin(2) ^ 2 + qJ(3) ^ 2) + t155; (m(4) * pkin(7) + mrSges(4,2)) * t90 + (t30 * mrSges(6,3) + t14) * t55 + (t57 * t89 - t58 * t93) * mrSges(5,3) + (-t29 * mrSges(6,3) + t151) * t59 + m(7) * (t109 * t59 + t134) + m(6) * (t13 * t59 + t134) + m(5) * (t38 * t93 + t39 * t89); m(7) * (t34 * t112 + t33 * t55) + m(6) * (-t35 * t55 + t36 * t59) + m(5) * (t63 * t93 + t64 * t89) + t150 - t152; m(4) + m(5) * (t89 ^ 2 + t93 ^ 2) + m(6) * (t52 + t145) + m(7) * (t119 * t52 + t145); (m(6) * (-t11 * t92 + t13 * t88) + (-t29 * t88 - t30 * t92) * mrSges(6,3)) * pkin(4) + t99 + t149 * t73 + t148 * t72 + t147; -t72 * t113 + m(7) * (t34 * t111 + t33 * t73) + (m(6) * (t35 * t92 + t36 * t88) - t108) * pkin(4) + t101 - t127 + t129 - t128 - Ifges(5,3); m(7) * (t59 * t111 + t73 * t55) + m(6) * (-t55 * t92 + t59 * t88) * pkin(4) + t152; 0.2e1 * t127 + 0.2e1 * t103 + t72 * t106 + m(7) * (t119 * t72 ^ 2 + t73 ^ 2) + m(6) * (t88 ^ 2 + t92 ^ 2) * pkin(4) ^ 2 + t155; -t149 * pkin(5) + t148 * pkin(10) + t104 + t147; m(7) * (-pkin(5) * t33 + t34 * t114) + t138 - pkin(10) * t113 + t101; m(7) * (-pkin(5) * t55 + pkin(10) * t112) + t154; m(7) * (-pkin(5) * t73 + pkin(10) * t111) + (-pkin(5) + t73) * t66 + t103 + (t111 + t114) * mrSges(7,3) + t117; -0.2e1 * t138 + m(7) * (t119 * pkin(10) ^ 2 + pkin(5) ^ 2) + pkin(10) * t106 + t117; mrSges(7,1) * t2 - mrSges(7,2) * t3 - Ifges(7,6) * t132 + t121; -t107 * t34 - t67; -t107 * t59; -t107 * t72 + t67; -t107 * pkin(10) + t67; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
