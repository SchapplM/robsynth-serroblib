% Calculate joint inertia matrix for
% S6RRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:16:43
% EndTime: 2019-03-09 13:16:45
% DurationCPUTime: 1.03s
% Computational Cost: add. (3228->278), mult. (6055->403), div. (0->0), fcn. (7082->10), ass. (0->110)
t118 = cos(qJ(5));
t115 = sin(qJ(4));
t158 = cos(qJ(4));
t111 = sin(pkin(11));
t112 = cos(pkin(11));
t116 = sin(qJ(2));
t119 = cos(qJ(2));
t83 = -t111 * t116 + t112 * t119;
t84 = t111 * t119 + t112 * t116;
t59 = t115 * t83 + t158 * t84;
t143 = t118 * t59;
t58 = t115 * t84 - t158 * t83;
t166 = Ifges(6,5) * t143 + Ifges(6,3) * t58;
t113 = sin(qJ(6));
t114 = sin(qJ(5));
t141 = Ifges(6,5) * t114 + Ifges(6,6) * t118;
t117 = cos(qJ(6));
t89 = -t113 * t114 + t117 * t118;
t152 = t89 * mrSges(7,3);
t165 = t113 * pkin(5) * t152 + t141;
t150 = -qJ(3) - pkin(7);
t92 = t150 * t116;
t93 = t150 * t119;
t64 = t111 * t93 + t112 * t92;
t125 = -pkin(8) * t84 + t64;
t65 = t111 * t92 - t112 * t93;
t51 = pkin(8) * t83 + t65;
t31 = t115 * t51 - t158 * t125;
t164 = t31 ^ 2;
t163 = 0.2e1 * t31;
t90 = t113 * t118 + t114 * t117;
t60 = -t89 * mrSges(7,1) + t90 * mrSges(7,2);
t162 = 0.2e1 * t60;
t103 = -pkin(2) * t119 - pkin(1);
t71 = -pkin(3) * t83 + t103;
t161 = 0.2e1 * t71;
t160 = 0.2e1 * t83;
t157 = Ifges(6,6) * t58;
t156 = pkin(2) * t111;
t155 = t118 * pkin(5);
t101 = pkin(2) * t112 + pkin(3);
t76 = t158 * t101 - t115 * t156;
t154 = t76 * mrSges(5,1);
t77 = t115 * t101 + t158 * t156;
t153 = t77 * mrSges(5,2);
t151 = t90 * mrSges(7,3);
t30 = pkin(4) * t58 - pkin(9) * t59 + t71;
t33 = t115 * t125 + t158 * t51;
t13 = t114 * t30 + t118 * t33;
t149 = Ifges(7,5) * t90 + Ifges(7,6) * t89;
t148 = Ifges(6,4) * t114;
t147 = Ifges(6,4) * t118;
t12 = -t114 * t33 + t118 * t30;
t146 = t114 * t12;
t145 = t114 * t59;
t144 = t118 * t13;
t74 = pkin(9) + t77;
t142 = t118 * t74;
t140 = t114 ^ 2 + t118 ^ 2;
t139 = t116 ^ 2 + t119 ^ 2;
t138 = 0.2e1 * mrSges(7,3);
t36 = t90 * t59;
t37 = t89 * t59;
t137 = Ifges(7,5) * t37 - Ifges(7,6) * t36 + Ifges(7,3) * t58;
t136 = t117 * t151;
t135 = -t83 * mrSges(4,1) + t84 * mrSges(4,2);
t134 = t140 * t74;
t61 = Ifges(7,4) * t90 + Ifges(7,2) * t89;
t62 = Ifges(7,1) * t90 + Ifges(7,4) * t89;
t94 = Ifges(6,2) * t118 + t148;
t95 = Ifges(6,1) * t114 + t147;
t133 = t114 * t95 + t118 * t94 + t89 * t61 + t90 * t62 + Ifges(5,3);
t91 = -t118 * mrSges(6,1) + t114 * mrSges(6,2);
t132 = mrSges(6,1) * t114 + mrSges(6,2) * t118;
t131 = t144 - t146;
t66 = (-pkin(10) - t74) * t114;
t106 = t118 * pkin(10);
t67 = t106 + t142;
t43 = -t113 * t67 + t117 * t66;
t44 = t113 * t66 + t117 * t67;
t130 = t43 * mrSges(7,1) - t44 * mrSges(7,2) + t149;
t97 = (-pkin(10) - pkin(9)) * t114;
t98 = pkin(9) * t118 + t106;
t68 = -t113 * t98 + t117 * t97;
t69 = t113 * t97 + t117 * t98;
t129 = t68 * mrSges(7,1) - t69 * mrSges(7,2) + t149;
t128 = 0.2e1 * t140 * mrSges(6,3);
t73 = -pkin(4) - t76;
t5 = pkin(5) * t58 - pkin(10) * t143 + t12;
t6 = -pkin(10) * t145 + t13;
t3 = -t113 * t6 + t117 * t5;
t4 = t113 * t5 + t117 * t6;
t127 = t3 * mrSges(7,1) - t4 * mrSges(7,2) + t137;
t126 = (mrSges(7,1) * t117 - mrSges(7,2) * t113) * pkin(5);
t10 = Ifges(7,1) * t37 - Ifges(7,4) * t36 + Ifges(7,5) * t58;
t16 = pkin(5) * t145 + t31;
t21 = t157 + (-Ifges(6,2) * t114 + t147) * t59;
t22 = Ifges(6,5) * t58 + (Ifges(6,1) * t118 - t148) * t59;
t9 = Ifges(7,4) * t37 - Ifges(7,2) * t36 + Ifges(7,6) * t58;
t124 = -t33 * mrSges(5,2) - t3 * t151 + t4 * t152 + mrSges(6,3) * t144 + t16 * t60 + t114 * t22 / 0.2e1 + t118 * t21 / 0.2e1 - t36 * t61 / 0.2e1 + t37 * t62 / 0.2e1 - t94 * t145 / 0.2e1 + t95 * t143 / 0.2e1 - Ifges(5,6) * t58 + Ifges(5,5) * t59 + t89 * t9 / 0.2e1 + t90 * t10 / 0.2e1 + (t91 - mrSges(5,1)) * t31 + (t149 + t141) * t58 / 0.2e1;
t102 = -pkin(4) - t155;
t70 = t73 - t155;
t53 = t59 * mrSges(5,2);
t40 = mrSges(6,1) * t58 - mrSges(6,3) * t143;
t39 = -mrSges(6,2) * t58 - mrSges(6,3) * t145;
t38 = t132 * t59;
t18 = mrSges(7,1) * t58 - mrSges(7,3) * t37;
t17 = -mrSges(7,2) * t58 - mrSges(7,3) * t36;
t14 = mrSges(7,1) * t36 + mrSges(7,2) * t37;
t1 = [-0.2e1 * pkin(1) * (-mrSges(3,1) * t119 + mrSges(3,2) * t116) + t116 * (Ifges(3,1) * t116 + Ifges(3,4) * t119) + t119 * (Ifges(3,4) * t116 + Ifges(3,2) * t119) + 0.2e1 * t103 * t135 + 0.2e1 * t139 * pkin(7) * mrSges(3,3) + 0.2e1 * t12 * t40 - t36 * t9 + t37 * t10 + 0.2e1 * t13 * t39 + 0.2e1 * t16 * t14 + 0.2e1 * t4 * t17 + 0.2e1 * t3 * t18 + m(4) * (t103 ^ 2 + t64 ^ 2 + t65 ^ 2) + m(7) * (t16 ^ 2 + t3 ^ 2 + t4 ^ 2) + Ifges(4,2) * t83 ^ 2 + (mrSges(5,1) * t161 - 0.2e1 * t33 * mrSges(5,3) + Ifges(5,2) * t58 + t137 + t166) * t58 + m(3) * (t139 * pkin(7) ^ 2 + pkin(1) ^ 2) + Ifges(2,3) + (-0.2e1 * t64 * mrSges(4,3) + Ifges(4,1) * t84 + Ifges(4,4) * t160) * t84 + (mrSges(5,3) * t163 + Ifges(5,1) * t59 - 0.2e1 * Ifges(5,4) * t58 + t118 * t22 + (-t21 - t157) * t114) * t59 + m(5) * (t33 ^ 2 + t71 ^ 2 + t164) + m(6) * (t12 ^ 2 + t13 ^ 2 + t164) + t65 * mrSges(4,3) * t160 + t53 * t161 + t38 * t163; Ifges(3,6) * t119 + t124 + (-t12 * mrSges(6,3) - t74 * t40) * t114 + m(5) * (-t31 * t76 + t33 * t77) + m(7) * (t16 * t70 + t3 * t43 + t4 * t44) + t39 * t142 + Ifges(3,5) * t116 + Ifges(4,6) * t83 + Ifges(4,5) * t84 + t64 * mrSges(4,1) - t65 * mrSges(4,2) + t70 * t14 + t73 * t38 + t43 * t18 + t44 * t17 + (-t116 * mrSges(3,1) - t119 * mrSges(3,2)) * pkin(7) + (-t77 * t58 - t76 * t59) * mrSges(5,3) + m(6) * (t131 * t74 + t31 * t73) + (m(4) * (t111 * t65 + t112 * t64) + (t111 * t83 - t112 * t84) * mrSges(4,3)) * pkin(2); 0.2e1 * t154 - 0.2e1 * t153 + t70 * t162 + 0.2e1 * t73 * t91 + Ifges(3,3) + Ifges(4,3) + (-t43 * t90 + t44 * t89) * t138 + t74 * t128 + m(7) * (t43 ^ 2 + t44 ^ 2 + t70 ^ 2) + m(6) * (t140 * t74 ^ 2 + t73 ^ 2) + m(5) * (t76 ^ 2 + t77 ^ 2) + t133 + (0.2e1 * t112 * mrSges(4,1) - 0.2e1 * t111 * mrSges(4,2) + m(4) * (t111 ^ 2 + t112 ^ 2) * pkin(2)) * pkin(2); t58 * mrSges(5,1) + t114 * t39 + t118 * t40 + t90 * t17 + t89 * t18 + t53 + m(7) * (t3 * t89 + t4 * t90) + m(6) * (t114 * t13 + t118 * t12) + m(5) * t71 + m(4) * t103 + t135; m(7) * (t43 * t89 + t44 * t90); m(4) + m(5) + m(6) * t140 + m(7) * (t89 ^ 2 + t90 ^ 2); (-t114 * t40 + t118 * t39) * pkin(9) + m(7) * (t102 * t16 + t3 * t68 + t4 * t69) + t124 - mrSges(6,3) * t146 + t102 * t14 + t68 * t18 + t69 * t17 - pkin(4) * t38 + m(6) * (-pkin(4) * t31 + t131 * pkin(9)); t154 - t153 + (t73 - pkin(4)) * t91 + (t70 + t102) * t60 + m(7) * (t102 * t70 + t43 * t68 + t44 * t69) + m(6) * (-pkin(4) * t73 + pkin(9) * t134) + ((-t43 - t68) * t90 + (t44 + t69) * t89) * mrSges(7,3) + (t140 * pkin(9) + t134) * mrSges(6,3) + t133; m(7) * (t68 * t89 + t69 * t90); -0.2e1 * pkin(4) * t91 + t102 * t162 + (-t68 * t90 + t69 * t89) * t138 + pkin(9) * t128 + m(7) * (t102 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(6) * (t140 * pkin(9) ^ 2 + pkin(4) ^ 2) + t133; -Ifges(6,6) * t145 + t12 * mrSges(6,1) - t13 * mrSges(6,2) + (m(7) * (t113 * t4 + t117 * t3) + t113 * t17 + t117 * t18) * pkin(5) + t127 + t166; -t132 * t74 + (m(7) * (t113 * t44 + t117 * t43) - t136) * pkin(5) + t130 + t165; m(7) * (t113 * t90 + t117 * t89) * pkin(5) - t91 - t60; -t132 * pkin(9) + (m(7) * (t113 * t69 + t117 * t68) - t136) * pkin(5) + t129 + t165; Ifges(6,3) + Ifges(7,3) + m(7) * (t113 ^ 2 + t117 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t126; t127; t130; -t60; t129; Ifges(7,3) + t126; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
