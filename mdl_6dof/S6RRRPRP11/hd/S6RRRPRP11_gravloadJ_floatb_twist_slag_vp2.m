% Calculate Gravitation load on the joints for
% S6RRRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP11_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP11_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP11_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP11_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:41:35
% EndTime: 2019-03-09 17:41:39
% DurationCPUTime: 1.32s
% Computational Cost: add. (653->128), mult. (1583->173), div. (0->0), fcn. (1860->10), ass. (0->70)
t128 = -m(5) - m(7);
t123 = m(6) - t128;
t127 = mrSges(4,2) - mrSges(5,3);
t137 = m(7) * pkin(5);
t55 = sin(qJ(5));
t64 = -qJ(4) * t123 - t55 * t137 + t127;
t136 = m(7) * (-qJ(6) - pkin(10)) - mrSges(7,3) - mrSges(4,1) + mrSges(5,2);
t92 = -mrSges(6,1) - mrSges(7,1);
t126 = mrSges(6,2) + mrSges(7,2);
t117 = pkin(3) * t123;
t56 = sin(qJ(3));
t60 = cos(qJ(3));
t99 = t55 * t56;
t133 = t127 * t56 + t136 * t60 - t99 * t137 - mrSges(3,1);
t119 = -m(6) * pkin(10) - mrSges(6,3);
t113 = -t119 - t136;
t132 = t113 + t117;
t59 = cos(qJ(5));
t131 = mrSges(3,2) - mrSges(5,1) - mrSges(4,3) - m(7) * (pkin(5) * t59 + pkin(4));
t57 = sin(qJ(2));
t58 = sin(qJ(1));
t61 = cos(qJ(2));
t108 = cos(qJ(1));
t87 = cos(pkin(6));
t75 = t87 * t108;
t37 = t57 * t75 + t58 * t61;
t53 = sin(pkin(6));
t83 = t53 * t108;
t17 = t37 * t56 + t60 * t83;
t36 = t57 * t58 - t61 * t75;
t130 = -t17 * t55 - t36 * t59;
t129 = t17 * t59 - t36 * t55;
t105 = t36 * t60;
t88 = qJ(4) * t56;
t125 = -pkin(3) * t105 - t36 * t88;
t79 = t58 * t87;
t38 = t108 * t57 + t61 * t79;
t104 = t38 * t60;
t124 = -pkin(3) * t104 - t38 * t88;
t120 = t92 - t137;
t67 = -m(6) * pkin(4) + t131;
t118 = pkin(9) * (m(4) + t123) - t67;
t116 = t60 * mrSges(6,3) - t133;
t115 = -t126 * t59 + t92 * t55 + t64;
t114 = -m(6) * (pkin(4) + pkin(9)) + t131;
t103 = t53 * t57;
t102 = t53 * t58;
t101 = t53 * t60;
t100 = t53 * t61;
t98 = t55 * t61;
t97 = t56 * t59;
t96 = t59 * t61;
t90 = pkin(2) * t100 + pkin(9) * t103;
t89 = t108 * pkin(1) + pkin(8) * t102;
t39 = t108 * t61 - t57 * t79;
t85 = t39 * pkin(2) + t89;
t82 = -pkin(1) * t58 + pkin(8) * t83;
t30 = t36 * pkin(2);
t81 = pkin(9) * t37 - t30;
t32 = t38 * pkin(2);
t80 = pkin(9) * t39 - t32;
t18 = t37 * t60 - t56 * t83;
t76 = -t37 * pkin(2) + t82;
t21 = -t101 * t58 + t39 * t56;
t5 = t21 * t59 - t38 * t55;
t35 = t101 * t57 + t56 * t87;
t34 = t103 * t56 - t60 * t87;
t22 = t102 * t56 + t39 * t60;
t6 = t21 * t55 + t38 * t59;
t1 = [(-t108 * mrSges(2,1) - m(3) * t89 - t39 * mrSges(3,1) - m(4) * t85 + t92 * t6 + (-mrSges(3,3) * t53 + mrSges(2,2)) * t58 - t126 * t5 + t64 * t21 - t113 * t22 - t118 * t38 - t123 * (t22 * pkin(3) + t85)) * g(2) + (t58 * mrSges(2,1) + t108 * mrSges(2,2) - m(3) * t82 + t37 * mrSges(3,1) - mrSges(3,3) * t83 - m(4) * t76 + t92 * t130 + t126 * t129 - t64 * t17 + t113 * t18 + t118 * t36 + t123 * (pkin(3) * t18 - t76)) * g(1) (-m(4) * t81 - m(6) * (-pkin(10) * t105 + t125 - t30) + t92 * (-t36 * t99 + t37 * t59) + t128 * (t81 + t125) - t126 * (-t36 * t97 - t37 * t55) + t114 * t37 + t116 * t36) * g(2) + (-m(4) * t80 - m(6) * (-pkin(10) * t104 + t124 - t32) - t126 * (-t38 * t97 - t39 * t55) + t128 * (t80 + t124) + t92 * (-t38 * t99 + t39 * t59) + t114 * t39 + t116 * t38) * g(1) + (-m(4) * t90 - t123 * (t88 * t100 + t90) + (t92 * (t56 * t98 + t57 * t59) - t126 * (-t55 * t57 + t56 * t96) + t67 * t57 + (t133 + (t119 - t117) * t60) * t61) * t53) * g(3) (t115 * t35 + t132 * t34) * g(3) + (t115 * t18 + t132 * t17) * g(2) + (t115 * t22 + t132 * t21) * g(1), t123 * (-g(1) * t21 - g(2) * t17 - g(3) * t34) (-t126 * (-t34 * t55 + t53 * t96) + t120 * (t34 * t59 + t53 * t98)) * g(3) + (t120 * t129 - t126 * t130) * g(2) + (t120 * t5 + t126 * t6) * g(1) (-g(1) * t22 - g(2) * t18 - g(3) * t35) * m(7)];
taug  = t1(:);
