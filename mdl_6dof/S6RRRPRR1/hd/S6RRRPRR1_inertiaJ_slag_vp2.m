% Calculate joint inertia matrix for
% S6RRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:02:23
% EndTime: 2019-03-09 18:02:25
% DurationCPUTime: 0.97s
% Computational Cost: add. (3429->239), mult. (6379->351), div. (0->0), fcn. (7510->10), ass. (0->107)
t95 = sin(qJ(6));
t127 = t95 * mrSges(7,3);
t138 = cos(qJ(5));
t100 = cos(qJ(3));
t101 = cos(qJ(2));
t97 = sin(qJ(3));
t98 = sin(qJ(2));
t72 = t100 * t101 - t97 * t98;
t73 = t100 * t98 + t101 * t97;
t93 = sin(pkin(11));
t94 = cos(pkin(11));
t48 = t72 * t94 - t73 * t93;
t49 = t72 * t93 + t73 * t94;
t96 = sin(qJ(5));
t30 = -t138 * t48 + t49 * t96;
t31 = t138 * t49 + t96 * t48;
t15 = -mrSges(7,2) * t30 - t31 * t127;
t99 = cos(qJ(6));
t132 = t31 * t99;
t16 = mrSges(7,1) * t30 - mrSges(7,3) * t132;
t113 = t99 * t15 - t95 * t16;
t142 = pkin(2) * t97;
t83 = pkin(2) * t100 + pkin(3);
t64 = -t93 * t142 + t94 * t83;
t61 = pkin(4) + t64;
t66 = t94 * t142 + t83 * t93;
t42 = t138 * t61 - t96 * t66;
t38 = -pkin(5) - t42;
t74 = -mrSges(7,1) * t99 + mrSges(7,2) * t95;
t33 = t38 * t74;
t89 = t95 ^ 2;
t137 = mrSges(7,3) * t89;
t43 = t138 * t66 + t96 * t61;
t39 = pkin(10) + t43;
t34 = t39 * t137;
t91 = t99 ^ 2;
t136 = mrSges(7,3) * t91;
t35 = t39 * t136;
t36 = t42 * mrSges(6,1);
t149 = t33 + t34 + t35 + t36;
t139 = t93 * pkin(3);
t82 = t94 * pkin(3) + pkin(4);
t65 = t138 * t82 - t96 * t139;
t62 = -pkin(5) - t65;
t50 = t62 * t74;
t67 = t138 * t139 + t96 * t82;
t63 = pkin(10) + t67;
t54 = t63 * t137;
t55 = t63 * t136;
t59 = t65 * mrSges(6,1);
t148 = t50 + t54 + t55 + t59;
t143 = -pkin(8) - pkin(7);
t79 = t143 * t98;
t80 = t143 * t101;
t52 = t100 * t79 + t80 * t97;
t109 = -qJ(4) * t73 + t52;
t53 = -t100 * t80 + t79 * t97;
t110 = qJ(4) * t72 + t53;
t20 = t94 * t109 - t93 * t110;
t107 = -t49 * pkin(9) + t20;
t21 = t93 * t109 + t94 * t110;
t18 = pkin(9) * t48 + t21;
t6 = -t138 * t107 + t18 * t96;
t147 = t6 ^ 2;
t146 = 0.2e1 * t6;
t84 = -pkin(2) * t101 - pkin(1);
t58 = -pkin(3) * t72 + t84;
t32 = -pkin(4) * t48 + t58;
t145 = 0.2e1 * t32;
t144 = 0.2e1 * t48;
t141 = pkin(5) * t74;
t13 = pkin(5) * t30 - pkin(10) * t31 + t32;
t8 = t96 * t107 + t138 * t18;
t3 = t13 * t95 + t8 * t99;
t140 = t3 * t99;
t135 = Ifges(7,4) * t95;
t134 = Ifges(7,4) * t99;
t133 = t31 * t95;
t131 = t43 * mrSges(6,2);
t130 = t64 * mrSges(5,1);
t129 = t66 * mrSges(5,2);
t128 = t67 * mrSges(6,2);
t124 = Ifges(7,5) * t132 + Ifges(7,3) * t30;
t123 = Ifges(7,5) * t95 + Ifges(7,6) * t99;
t122 = t89 + t91;
t121 = t101 ^ 2 + t98 ^ 2;
t75 = Ifges(7,2) * t99 + t135;
t76 = Ifges(7,1) * t95 + t134;
t120 = t99 * t75 + t95 * t76 + Ifges(6,3);
t119 = t122 * t63;
t118 = -t48 * mrSges(5,1) + t49 * mrSges(5,2);
t117 = Ifges(4,3) + Ifges(5,3) + t120;
t2 = t13 * t99 - t8 * t95;
t116 = -t2 * t95 + t140;
t115 = t94 * mrSges(5,1) - t93 * mrSges(5,2);
t114 = mrSges(7,1) * t95 + mrSges(7,2) * t99;
t112 = (t100 * mrSges(4,1) - t97 * mrSges(4,2)) * pkin(2);
t85 = pkin(10) * t137;
t86 = pkin(10) * t136;
t111 = t120 + t85 + t86 - t141;
t11 = Ifges(7,6) * t30 + (-Ifges(7,2) * t95 + t134) * t31;
t12 = Ifges(7,5) * t30 + (Ifges(7,1) * t99 - t135) * t31;
t108 = -t8 * mrSges(6,2) + mrSges(7,3) * t140 - t2 * t127 + t99 * t11 / 0.2e1 - t75 * t133 / 0.2e1 + t76 * t132 / 0.2e1 + Ifges(6,5) * t31 + t95 * t12 / 0.2e1 + (-mrSges(6,1) + t74) * t6 + (t123 / 0.2e1 - Ifges(6,6)) * t30;
t106 = t52 * mrSges(4,1) + t20 * mrSges(5,1) - t53 * mrSges(4,2) - t21 * mrSges(5,2) + Ifges(4,5) * t73 + Ifges(5,5) * t49 + Ifges(4,6) * t72 + Ifges(5,6) * t48 + t108;
t26 = t31 * mrSges(6,2);
t14 = t114 * t31;
t1 = [-0.2e1 * pkin(1) * (-mrSges(3,1) * t101 + mrSges(3,2) * t98) + t98 * (Ifges(3,1) * t98 + Ifges(3,4) * t101) + t101 * (Ifges(3,4) * t98 + Ifges(3,2) * t101) + 0.2e1 * t84 * (-mrSges(4,1) * t72 + mrSges(4,2) * t73) + t72 * (Ifges(4,4) * t73 + Ifges(4,2) * t72) + t73 * (Ifges(4,1) * t73 + Ifges(4,4) * t72) + Ifges(5,2) * t48 ^ 2 + 0.2e1 * t58 * t118 + t21 * mrSges(5,3) * t144 + t26 * t145 + t14 * t146 + 0.2e1 * t3 * t15 + 0.2e1 * t2 * t16 + Ifges(2,3) + (mrSges(6,1) * t145 - 0.2e1 * t8 * mrSges(6,3) + Ifges(6,2) * t30 + t124) * t30 + (-0.2e1 * t20 * mrSges(5,3) + Ifges(5,1) * t49 + Ifges(5,4) * t144) * t49 + (mrSges(6,3) * t146 + Ifges(6,1) * t31 - t95 * t11 + t99 * t12 + (-Ifges(7,6) * t95 - (2 * Ifges(6,4))) * t30) * t31 + m(3) * (t121 * pkin(7) ^ 2 + pkin(1) ^ 2) + m(4) * (t52 ^ 2 + t53 ^ 2 + t84 ^ 2) + m(5) * (t20 ^ 2 + t21 ^ 2 + t58 ^ 2) + m(6) * (t32 ^ 2 + t8 ^ 2 + t147) + m(7) * (t2 ^ 2 + t3 ^ 2 + t147) + 0.2e1 * (-t52 * t73 + t53 * t72) * mrSges(4,3) + 0.2e1 * t121 * pkin(7) * mrSges(3,3); Ifges(3,6) * t101 + Ifges(3,5) * t98 + (-t98 * mrSges(3,1) - t101 * mrSges(3,2)) * pkin(7) + (-t43 * t30 - t42 * t31) * mrSges(6,3) + (t66 * t48 - t64 * t49) * mrSges(5,3) + m(7) * (t116 * t39 + t38 * t6) + t38 * t14 + t113 * t39 + t106 + m(5) * (t20 * t64 + t21 * t66) + m(6) * (-t42 * t6 + t43 * t8) + (m(4) * (t100 * t52 + t53 * t97) + (-t100 * t73 + t72 * t97) * mrSges(4,3)) * pkin(2); 0.2e1 * t130 - 0.2e1 * t129 - 0.2e1 * t131 + Ifges(3,3) + 0.2e1 * t33 + 0.2e1 * t34 + 0.2e1 * t35 + 0.2e1 * t36 + 0.2e1 * t112 + m(7) * (t122 * t39 ^ 2 + t38 ^ 2) + m(6) * (t42 ^ 2 + t43 ^ 2) + m(5) * (t64 ^ 2 + t66 ^ 2) + m(4) * (t100 ^ 2 + t97 ^ 2) * pkin(2) ^ 2 + t117; t113 * t63 + t62 * t14 + (m(5) * (t20 * t94 + t21 * t93) + (t93 * t48 - t94 * t49) * mrSges(5,3)) * pkin(3) + m(7) * (t116 * t63 + t6 * t62) + m(6) * (-t6 * t65 + t67 * t8) + t106 + (-t67 * t30 - t65 * t31) * mrSges(6,3); t130 - t129 + m(7) * (t39 * t119 + t38 * t62) + (m(5) * (t64 * t94 + t66 * t93) + t115) * pkin(3) + t112 + t117 + (-t43 - t67) * mrSges(6,2) + m(6) * (t42 * t65 + t43 * t67) + t148 + t149; -0.2e1 * t128 + 0.2e1 * t50 + 0.2e1 * t54 + 0.2e1 * t55 + 0.2e1 * t59 + m(7) * (t122 * t63 ^ 2 + t62 ^ 2) + m(6) * (t65 ^ 2 + t67 ^ 2) + t117 + (0.2e1 * t115 + m(5) * (t93 ^ 2 + t94 ^ 2) * pkin(3)) * pkin(3); t30 * mrSges(6,1) + t95 * t15 + t99 * t16 + t26 + m(7) * (t2 * t99 + t3 * t95) + m(6) * t32 + m(5) * t58 + t118; 0; 0; m(7) * t122 + m(5) + m(6); t108 + (-m(7) * t6 - t14) * pkin(5) + (m(7) * t116 + t113) * pkin(10); m(7) * (t122 * t39 * pkin(10) - pkin(5) * t38) - t131 + t111 + t149; m(7) * (-pkin(5) * t62 + pkin(10) * t119) - t128 + t111 + t148; 0; -0.2e1 * t141 + m(7) * (t122 * pkin(10) ^ 2 + pkin(5) ^ 2) + 0.2e1 * t85 + 0.2e1 * t86 + t120; mrSges(7,1) * t2 - mrSges(7,2) * t3 - Ifges(7,6) * t133 + t124; -t114 * t39 + t123; -t114 * t63 + t123; -t74; -t114 * pkin(10) + t123; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
