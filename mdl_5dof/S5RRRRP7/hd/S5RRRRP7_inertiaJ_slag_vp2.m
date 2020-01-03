% Calculate joint inertia matrix for
% S5RRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP7_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP7_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP7_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:56:08
% EndTime: 2019-12-31 21:56:11
% DurationCPUTime: 0.86s
% Computational Cost: add. (896->183), mult. (1736->250), div. (0->0), fcn. (1621->6), ass. (0->80)
t142 = Ifges(5,1) + Ifges(6,1);
t141 = Ifges(6,4) + Ifges(5,5);
t140 = Ifges(6,2) + Ifges(5,3);
t81 = sin(qJ(4));
t120 = Ifges(6,5) * t81;
t122 = Ifges(5,4) * t81;
t82 = sin(qJ(3));
t83 = sin(qJ(2));
t85 = cos(qJ(3));
t86 = cos(qJ(2));
t50 = t82 * t83 - t85 * t86;
t51 = t82 * t86 + t85 * t83;
t84 = cos(qJ(4));
t139 = (t142 * t84 + t120 - t122) * t51 + t141 * t50;
t119 = Ifges(6,5) * t84;
t121 = Ifges(5,4) * t84;
t138 = t142 * t81 - t119 + t121;
t137 = t81 ^ 2 + t84 ^ 2;
t136 = (Ifges(5,6) - Ifges(6,6)) * t84 + t141 * t81;
t135 = (mrSges(5,3) + mrSges(6,2)) * t137;
t128 = -pkin(7) - pkin(6);
t107 = t128 * t83;
t62 = t128 * t86;
t33 = -t85 * t107 - t82 * t62;
t134 = t33 ^ 2;
t133 = 0.2e1 * t33;
t56 = -t84 * mrSges(6,1) - t81 * mrSges(6,3);
t132 = 0.2e1 * t56;
t69 = -t86 * pkin(2) - pkin(1);
t131 = 0.2e1 * t69;
t127 = m(6) * t81;
t22 = t50 * pkin(3) - t51 * pkin(8) + t69;
t35 = t82 * t107 - t85 * t62;
t8 = t81 * t22 + t84 * t35;
t3 = t50 * qJ(5) + t8;
t125 = t3 * t84;
t124 = t8 * t84;
t123 = t85 * pkin(2);
t118 = Ifges(5,6) * t50;
t117 = t51 * t81;
t116 = t51 * t84;
t71 = t81 * mrSges(6,2);
t115 = t81 * mrSges(5,3);
t67 = t82 * pkin(2) + pkin(8);
t114 = t137 * pkin(8) * t67;
t113 = t137 * t67 ^ 2;
t111 = t137 * pkin(8) ^ 2;
t109 = t83 ^ 2 + t86 ^ 2;
t108 = qJ(5) * t84;
t25 = -t50 * mrSges(6,1) + mrSges(6,2) * t116;
t103 = Ifges(6,6) * t117 + t141 * t116 + t140 * t50;
t7 = t84 * t22 - t81 * t35;
t4 = -t50 * pkin(4) - t7;
t102 = t4 * t81 + t125;
t101 = -t7 * t81 + t124;
t100 = t81 * mrSges(5,1) + t84 * mrSges(5,2);
t99 = t81 * mrSges(6,1) - t84 * mrSges(6,3);
t98 = pkin(4) * t81 - t108;
t53 = -t84 * pkin(4) - t81 * qJ(5) - pkin(3);
t58 = -Ifges(6,3) * t84 + t120;
t59 = Ifges(5,2) * t84 + t122;
t97 = Ifges(4,3) + (-t58 + t59) * t84 + t138 * t81;
t96 = (t85 * mrSges(4,1) - t82 * mrSges(4,2)) * pkin(2);
t23 = -t50 * mrSges(5,2) - t51 * t115;
t24 = t50 * mrSges(5,1) - mrSges(5,3) * t116;
t26 = t50 * mrSges(6,3) - t51 * t71;
t94 = (t23 + t26) * t84 + (-t24 + t25) * t81;
t93 = mrSges(6,2) * t108 - pkin(4) * t71 + t136;
t92 = 0.2e1 * t135;
t91 = -m(6) * t98 - t100 - t99;
t14 = Ifges(6,6) * t50 + (Ifges(6,3) * t81 + t119) * t51;
t15 = t118 + (-Ifges(5,2) * t81 + t121) * t51;
t57 = -t84 * mrSges(5,1) + t81 * mrSges(5,2);
t9 = t98 * t51 + t33;
t90 = -t35 * mrSges(4,2) + mrSges(6,2) * t125 + mrSges(5,3) * t124 + Ifges(4,5) * t51 - t7 * t115 + t4 * t71 + t9 * t56 + (-t14 / 0.2e1 + t15 / 0.2e1) * t84 + (t57 - mrSges(4,1)) * t33 + t139 * t81 / 0.2e1 + (t58 / 0.2e1 - t59 / 0.2e1) * t117 + t138 * t116 / 0.2e1 + (-Ifges(4,6) + t136 / 0.2e1) * t50;
t68 = -pkin(3) - t123;
t42 = t53 - t123;
t21 = t100 * t51;
t20 = t99 * t51;
t1 = [-0.2e1 * pkin(1) * (-t86 * mrSges(3,1) + t83 * mrSges(3,2)) + t83 * (Ifges(3,1) * t83 + Ifges(3,4) * t86) + t86 * (Ifges(3,4) * t83 + Ifges(3,2) * t86) + 0.2e1 * t9 * t20 + 0.2e1 * t8 * t23 + 0.2e1 * t7 * t24 + 0.2e1 * t4 * t25 + 0.2e1 * t3 * t26 + t21 * t133 + Ifges(2,3) + 0.2e1 * t109 * pkin(6) * mrSges(3,3) + (mrSges(4,1) * t131 - 0.2e1 * t35 * mrSges(4,3) + Ifges(4,2) * t50 + t103) * t50 + m(6) * (t3 ^ 2 + t4 ^ 2 + t9 ^ 2) + m(5) * (t7 ^ 2 + t8 ^ 2 + t134) + m(4) * (t35 ^ 2 + t69 ^ 2 + t134) + m(3) * (t109 * pkin(6) ^ 2 + pkin(1) ^ 2) + (mrSges(4,2) * t131 + mrSges(4,3) * t133 + Ifges(4,1) * t51 - 0.2e1 * Ifges(4,4) * t50 + t139 * t84 + (t14 - t15 - t118) * t81) * t51; (m(4) * (-t33 * t85 + t35 * t82) + (-t82 * t50 - t85 * t51) * mrSges(4,3)) * pkin(2) + m(5) * (t101 * t67 + t68 * t33) + m(6) * (t102 * t67 + t42 * t9) + t94 * t67 + (-t83 * mrSges(3,1) - t86 * mrSges(3,2)) * pkin(6) + Ifges(3,6) * t86 + Ifges(3,5) * t83 + t68 * t21 + t42 * t20 + t90; t42 * t132 + 0.2e1 * t68 * t57 + Ifges(3,3) + 0.2e1 * t96 + m(6) * (t42 ^ 2 + t113) + m(5) * (t68 ^ 2 + t113) + m(4) * (t82 ^ 2 + t85 ^ 2) * pkin(2) ^ 2 + t92 * t67 + t97; m(6) * (t102 * pkin(8) + t53 * t9) + m(5) * (-pkin(3) * t33 + t101 * pkin(8)) + t94 * pkin(8) + t53 * t20 - pkin(3) * t21 + t90; (t68 - pkin(3)) * t57 + (t42 + t53) * t56 + t96 + m(6) * (t53 * t42 + t114) + m(5) * (-pkin(3) * t68 + t114) + t97 + (pkin(8) + t67) * t135; -0.2e1 * pkin(3) * t57 + t53 * t132 + m(6) * (t53 ^ 2 + t111) + m(5) * (pkin(3) ^ 2 + t111) + t92 * pkin(8) + t97; -Ifges(5,6) * t117 - pkin(4) * t25 + t3 * mrSges(6,3) + m(6) * (-pkin(4) * t4 + qJ(5) * t3) + qJ(5) * t26 - t4 * mrSges(6,1) + t7 * mrSges(5,1) - t8 * mrSges(5,2) + t103; t91 * t67 + t93; t91 * pkin(8) + t93; 0.2e1 * pkin(4) * mrSges(6,1) + 0.2e1 * qJ(5) * mrSges(6,3) + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2) + t140; m(6) * t4 + t25; t67 * t127 + t71; pkin(8) * t127 + t71; -m(6) * pkin(4) - mrSges(6,1); m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
