% Calculate joint inertia matrix for
% S5RRRRP6
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
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP6_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP6_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP6_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:52:50
% EndTime: 2019-12-31 21:52:53
% DurationCPUTime: 0.98s
% Computational Cost: add. (892->192), mult. (1715->264), div. (0->0), fcn. (1624->6), ass. (0->76)
t141 = Ifges(5,4) + Ifges(6,4);
t136 = Ifges(5,6) + Ifges(6,6);
t82 = sin(qJ(3));
t83 = sin(qJ(2));
t85 = cos(qJ(3));
t86 = cos(qJ(2));
t52 = t82 * t83 - t85 * t86;
t140 = t136 * t52;
t139 = Ifges(5,1) + Ifges(6,1);
t138 = Ifges(5,5) + Ifges(6,5);
t137 = Ifges(5,2) + Ifges(6,2);
t84 = cos(qJ(4));
t135 = t141 * t84;
t81 = sin(qJ(4));
t134 = t141 * t81;
t133 = Ifges(5,3) + Ifges(6,3);
t53 = t82 * t86 + t85 * t83;
t132 = (-t137 * t81 + t135) * t53 + t140;
t131 = (t139 * t84 - t134) * t53 + t138 * t52;
t130 = t137 * t84 + t134;
t129 = t139 * t81 + t135;
t97 = t136 * t84 + t138 * t81;
t121 = -pkin(7) - pkin(6);
t102 = t121 * t83;
t64 = t121 * t86;
t33 = -t85 * t102 - t82 * t64;
t128 = t33 ^ 2;
t127 = 0.2e1 * t33;
t57 = -t84 * mrSges(6,1) + t81 * mrSges(6,2);
t126 = 0.2e1 * t57;
t70 = -t86 * pkin(2) - pkin(1);
t125 = 0.2e1 * t70;
t22 = t52 * pkin(3) - t53 * pkin(8) + t70;
t35 = t102 * t82 - t85 * t64;
t6 = t81 * t22 + t84 * t35;
t120 = t6 * t84;
t119 = t84 * pkin(8);
t118 = t85 * pkin(2);
t113 = t53 * t81;
t112 = t53 * t84;
t111 = t81 * mrSges(6,3);
t110 = t84 * mrSges(6,3);
t67 = t82 * pkin(2) + pkin(8);
t109 = t84 * t67;
t20 = mrSges(6,1) * t113 + mrSges(6,2) * t112;
t105 = t81 ^ 2 + t84 ^ 2;
t104 = t83 ^ 2 + t86 ^ 2;
t71 = t84 * qJ(5);
t103 = 0.2e1 * mrSges(6,3);
t69 = -t84 * pkin(4) - pkin(3);
t99 = t105 * t67;
t5 = t84 * t22 - t81 * t35;
t98 = t138 * t112 + t133 * t52;
t96 = t129 * t81 + t130 * t84 + Ifges(4,3);
t1 = t52 * pkin(4) - t53 * t71 + t5;
t95 = -t5 * mrSges(5,3) - t1 * mrSges(6,3);
t94 = -t5 * t81 + t120;
t93 = mrSges(5,1) * t81 + mrSges(5,2) * t84;
t92 = 0.2e1 * mrSges(5,3) * t105;
t91 = (t85 * mrSges(4,1) - t82 * mrSges(4,2)) * pkin(2);
t16 = pkin(4) * t113 + t33;
t3 = -qJ(5) * t113 + t6;
t58 = -t84 * mrSges(5,1) + t81 * mrSges(5,2);
t90 = -t35 * mrSges(4,2) + mrSges(5,3) * t120 + Ifges(4,5) * t53 + t3 * t110 + t16 * t57 + (-mrSges(4,1) + t58) * t33 + t131 * t81 / 0.2e1 + t132 * t84 / 0.2e1 - t130 * t113 / 0.2e1 + t129 * t112 / 0.2e1 + (-Ifges(4,6) + t97 / 0.2e1) * t52;
t68 = -pkin(3) - t118;
t59 = t71 + t119;
t56 = (-qJ(5) - pkin(8)) * t81;
t55 = t69 - t118;
t43 = t71 + t109;
t42 = (-qJ(5) - t67) * t81;
t26 = t52 * mrSges(5,1) - mrSges(5,3) * t112;
t25 = t52 * mrSges(6,1) - t110 * t53;
t24 = -t52 * mrSges(5,2) - mrSges(5,3) * t113;
t23 = -t52 * mrSges(6,2) - t111 * t53;
t21 = t93 * t53;
t2 = [t83 * (Ifges(3,1) * t83 + Ifges(3,4) * t86) + t86 * (Ifges(3,4) * t83 + Ifges(3,2) * t86) - 0.2e1 * pkin(1) * (-t86 * mrSges(3,1) + t83 * mrSges(3,2)) + 0.2e1 * t3 * t23 + 0.2e1 * t6 * t24 + 0.2e1 * t1 * t25 + 0.2e1 * t5 * t26 + t21 * t127 + 0.2e1 * t16 * t20 + Ifges(2,3) + 0.2e1 * t104 * pkin(6) * mrSges(3,3) + (mrSges(4,1) * t125 - 0.2e1 * t35 * mrSges(4,3) + Ifges(4,2) * t52 + t98) * t52 + m(3) * (pkin(6) ^ 2 * t104 + pkin(1) ^ 2) + m(4) * (t35 ^ 2 + t70 ^ 2 + t128) + m(6) * (t1 ^ 2 + t16 ^ 2 + t3 ^ 2) + m(5) * (t5 ^ 2 + t6 ^ 2 + t128) + (mrSges(4,2) * t125 + mrSges(4,3) * t127 + Ifges(4,1) * t53 - 0.2e1 * Ifges(4,4) * t52 + t131 * t84 + (-t132 - t140) * t81) * t53; (m(4) * (-t33 * t85 + t35 * t82) + (-t82 * t52 - t85 * t53) * mrSges(4,3)) * pkin(2) + m(5) * (t68 * t33 + t67 * t94) + (-t67 * t26 + t95) * t81 + m(6) * (t42 * t1 + t55 * t16 + t43 * t3) + (-t83 * mrSges(3,1) - t86 * mrSges(3,2)) * pkin(6) + t24 * t109 + Ifges(3,6) * t86 + Ifges(3,5) * t83 + t68 * t21 + t55 * t20 + t42 * t25 + t43 * t23 + t90; t55 * t126 + 0.2e1 * t68 * t58 + Ifges(3,3) + 0.2e1 * t91 + (-t42 * t81 + t43 * t84) * t103 + t67 * t92 + m(5) * (t105 * t67 ^ 2 + t68 ^ 2) + m(6) * (t42 ^ 2 + t43 ^ 2 + t55 ^ 2) + m(4) * (t82 ^ 2 + t85 ^ 2) * pkin(2) ^ 2 + t96; m(5) * (-pkin(3) * t33 + pkin(8) * t94) + (-pkin(8) * t26 + t95) * t81 + m(6) * (t56 * t1 + t69 * t16 + t59 * t3) + t24 * t119 + t69 * t20 + t56 * t25 + t59 * t23 - pkin(3) * t21 + t90; (t68 - pkin(3)) * t58 + (t55 + t69) * t57 + t91 + m(5) * (-pkin(3) * t68 + pkin(8) * t99) + m(6) * (t56 * t42 + t59 * t43 + t69 * t55) + ((t43 + t59) * t84 + (-t42 - t56) * t81) * mrSges(6,3) + (pkin(8) * t105 + t99) * mrSges(5,3) + t96; -0.2e1 * pkin(3) * t58 + t69 * t126 + (-t56 * t81 + t59 * t84) * t103 + pkin(8) * t92 + m(6) * (t56 ^ 2 + t59 ^ 2 + t69 ^ 2) + m(5) * (pkin(8) ^ 2 * t105 + pkin(3) ^ 2) + t96; t5 * mrSges(5,1) + t1 * mrSges(6,1) - t6 * mrSges(5,2) - t3 * mrSges(6,2) - t136 * t113 + (m(6) * t1 + t25) * pkin(4) + t98; t42 * mrSges(6,1) - t43 * mrSges(6,2) - t93 * t67 + (m(6) * t42 - t111) * pkin(4) + t97; t56 * mrSges(6,1) - t59 * mrSges(6,2) - t93 * pkin(8) + (m(6) * t56 - t111) * pkin(4) + t97; (m(6) * pkin(4) + 0.2e1 * mrSges(6,1)) * pkin(4) + t133; m(6) * t16 + t20; m(6) * t55 + t57; m(6) * t69 + t57; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
