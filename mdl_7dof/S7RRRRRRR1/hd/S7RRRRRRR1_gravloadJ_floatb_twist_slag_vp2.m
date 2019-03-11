% Calculate Gravitation load on the joints for
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% m_mdh [8x1]
%   mass of all robot links (including the base)
% mrSges [8x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [7x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 08:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S7RRRRRRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(3,1),zeros(4,1),zeros(8,1),zeros(8,3)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [7x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S7RRRRRRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [8 1]), ...
  'S7RRRRRRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [8x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [8,3]), ...
  'S7RRRRRRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [8x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 06:29:16
% EndTime: 2019-03-10 06:29:21
% DurationCPUTime: 1.77s
% Computational Cost: add. (1119->187), mult. (3035->291), div. (0->0), fcn. (3749->14), ass. (0->92)
t87 = sin(qJ(2));
t88 = sin(qJ(1));
t134 = t87 * t88;
t93 = cos(qJ(3));
t94 = cos(qJ(2));
t127 = t93 * t94;
t86 = sin(qJ(3));
t95 = cos(qJ(1));
t69 = t88 * t127 + t86 * t95;
t85 = sin(qJ(4));
t92 = cos(qJ(4));
t50 = t85 * t134 + t69 * t92;
t124 = t95 * t93;
t130 = t88 * t94;
t68 = t86 * t130 - t124;
t84 = sin(qJ(5));
t91 = cos(qJ(5));
t24 = t50 * t91 - t68 * t84;
t133 = t87 * t92;
t49 = -t88 * t133 + t69 * t85;
t83 = sin(qJ(6));
t90 = cos(qJ(6));
t153 = t24 * t83 - t49 * t90;
t4 = t24 * t90 + t49 * t83;
t82 = sin(qJ(7));
t89 = cos(qJ(7));
t152 = t89 * mrSges(8,1) - t82 * mrSges(8,2) + mrSges(7,1);
t108 = m(8) * pkin(4) + mrSges(7,2) + mrSges(8,3);
t119 = mrSges(6,3) - mrSges(5,2);
t143 = -m(6) - m(7);
t146 = (m(8) - t143) * pkin(3) + t119;
t99 = mrSges(5,1) * t92 + t146 * t85 + mrSges(4,1);
t120 = mrSges(6,2) - mrSges(7,3);
t149 = t82 * mrSges(8,1) + t89 * mrSges(8,2) + t120;
t23 = t50 * t84 + t68 * t91;
t148 = -t108 * t83 + t152 * t90 + mrSges(6,1);
t121 = mrSges(4,2) + mrSges(5,3);
t145 = -mrSges(4,1) * t93 + t121 * t86 - mrSges(3,1);
t144 = -m(4) - m(5);
t140 = t83 * t85;
t139 = t83 * t91;
t138 = t84 * t92;
t137 = t85 * t90;
t136 = t86 * t87;
t135 = t86 * t94;
t132 = t87 * t93;
t131 = t87 * t95;
t129 = t90 * t91;
t128 = t91 * t92;
t126 = t94 * t85;
t125 = t94 * t95;
t123 = mrSges(2,2) - mrSges(3,3);
t122 = -mrSges(3,2) - mrSges(4,3);
t117 = t84 * t136;
t116 = t85 * t136;
t115 = t91 * t136;
t114 = t86 * t131;
t112 = -mrSges(3,1) * t94 - mrSges(2,1);
t109 = m(8) * pkin(3) + t119;
t67 = t92 * t132 - t126;
t66 = t85 * t132 + t92 * t94;
t103 = (m(8) - t144) * pkin(2) - t122;
t100 = t103 * t87;
t97 = mrSges(6,1) * t91 - t149 * t84 + mrSges(5,1);
t96 = t103 * t94 - t145 * t87;
t77 = pkin(2) * t134;
t73 = t94 * t124 - t88 * t86;
t72 = -t86 * t125 - t88 * t93;
t71 = t92 * t127 + t85 * t87;
t70 = t93 * t126 - t133;
t63 = t67 * t95;
t62 = t66 * t95;
t61 = t67 * t88;
t60 = t66 * t88;
t59 = (-t86 * t128 - t84 * t93) * t87;
t56 = t85 * t131 + t73 * t92;
t55 = -t92 * t131 + t73 * t85;
t54 = -t84 * t135 + t71 * t91;
t48 = t67 * t91 - t117;
t47 = t67 * t84 + t115;
t44 = t84 * t114 - t63 * t91;
t42 = t88 * t117 - t61 * t91;
t36 = t72 * t128 - t73 * t84;
t34 = -t68 * t128 - t69 * t84;
t30 = t56 * t91 + t72 * t84;
t29 = t56 * t84 - t72 * t91;
t22 = t48 * t90 + t66 * t83;
t8 = t30 * t90 + t55 * t83;
t7 = -t30 * t83 + t55 * t90;
t2 = -t29 * t82 + t8 * t89;
t1 = -t29 * t89 - t8 * t82;
t3 = [(-t73 * mrSges(4,1) - t56 * mrSges(5,1) - mrSges(6,1) * t30 - mrSges(7,1) * t8 - t2 * mrSges(8,1) - t1 * mrSges(8,2) + t123 * t88 - t121 * t72 - t108 * t7 - t109 * t55 + t143 * (-pkin(2) * t131 + t55 * pkin(3)) + t120 * t29 + (t100 + t112) * t95) * g(2) + (-m(8) * t77 + mrSges(4,1) * t69 + mrSges(5,1) * t50 + mrSges(6,1) * t24 + t123 * t95 - t121 * t68 + t109 * t49 - t108 * t153 + t143 * (-pkin(3) * t49 + t77) + t152 * t4 - t149 * t23 + ((t144 * pkin(2) + t122) * t87 - t112) * t88) * g(1) (-mrSges(5,1) * t71 - mrSges(6,1) * t54 + t143 * (-pkin(2) * t87 + pkin(3) * t70) - t152 * (t54 * t90 + t70 * t83) + t149 * (t91 * t135 + t71 * t84) - t109 * t70 - t108 * (-t54 * t83 + t70 * t90) + t145 * t94 + t100) * g(3) + (mrSges(5,1) * t61 - mrSges(6,1) * t42 + t143 * (-pkin(2) * t130 - pkin(3) * t60) - t152 * (t42 * t90 - t60 * t83) + t149 * (-t88 * t115 - t61 * t84) + t109 * t60 - t108 * (-t42 * t83 - t60 * t90) + t96 * t88) * g(2) + (t63 * mrSges(5,1) - mrSges(6,1) * t44 + t143 * (-pkin(2) * t125 - t62 * pkin(3)) - t152 * (t44 * t90 - t62 * t83) + t149 * (-t91 * t114 - t63 * t84) + t109 * t62 - t108 * (-t44 * t83 - t62 * t90) + t96 * t95) * g(1) (-mrSges(6,1) * t59 - t152 * (-t83 * t116 + t59 * t90) - t149 * (t92 * t117 - t91 * t132) - t108 * (-t90 * t116 - t59 * t83) + (t121 * t93 + t99 * t86) * t87) * g(3) + (-mrSges(6,1) * t34 + t121 * t69 - t152 * (-t68 * t140 + t34 * t90) + t149 * (-t68 * t138 + t69 * t91) - t108 * (-t68 * t137 - t34 * t83) + t99 * t68) * g(2) + (-mrSges(6,1) * t36 + t121 * t73 - t152 * (t72 * t140 + t36 * t90) + t149 * (t72 * t138 + t73 * t91) - t108 * (t72 * t137 - t36 * t83) - t99 * t72) * g(1) (-t152 * (-t66 * t129 + t67 * t83) - t108 * (t66 * t139 + t67 * t90) - t146 * t67 + t97 * t66) * g(3) + (-t108 * (t49 * t139 + t50 * t90) - t152 * (-t49 * t129 + t50 * t83) - t146 * t50 + t97 * t49) * g(2) + (-t152 * (-t55 * t129 + t56 * t83) - t108 * (t55 * t139 + t56 * t90) - t146 * t56 + t97 * t55) * g(1) (t148 * t47 + t149 * t48) * g(3) + (t148 * t23 + t149 * t24) * g(2) + (t148 * t29 + t149 * t30) * g(1) (t108 * t22 - t152 * (-t48 * t83 + t66 * t90)) * g(3) + (t108 * t4 + t152 * t153) * g(2) + (t108 * t8 - t152 * t7) * g(1), -g(1) * (mrSges(8,1) * t1 - mrSges(8,2) * t2) - g(2) * ((-t23 * t89 - t4 * t82) * mrSges(8,1) + (t23 * t82 - t4 * t89) * mrSges(8,2)) - g(3) * ((-t22 * t82 - t47 * t89) * mrSges(8,1) + (-t22 * t89 + t47 * t82) * mrSges(8,2))];
taug  = t3(:);
