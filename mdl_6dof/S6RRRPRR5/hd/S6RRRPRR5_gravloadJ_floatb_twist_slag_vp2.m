% Calculate Gravitation load on the joints for
% S6RRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:20:39
% EndTime: 2019-03-09 18:20:41
% DurationCPUTime: 0.90s
% Computational Cost: add. (502->153), mult. (574->171), div. (0->0), fcn. (503->10), ass. (0->81)
t143 = m(5) + m(7);
t54 = qJ(5) + qJ(6);
t51 = cos(t54);
t55 = qJ(2) + qJ(3);
t52 = cos(t55);
t111 = t51 * t52;
t49 = sin(t54);
t117 = t49 * t52;
t59 = cos(qJ(5));
t118 = mrSges(6,2) * t59;
t146 = -mrSges(7,1) * t117 - mrSges(7,2) * t111 - t52 * t118;
t58 = sin(qJ(1));
t61 = cos(qJ(1));
t134 = g(1) * t61 + g(2) * t58;
t50 = sin(t55);
t144 = (-mrSges(4,1) + mrSges(5,2)) * t52 + (mrSges(4,2) - mrSges(5,3)) * t50;
t62 = -pkin(10) - pkin(9);
t112 = t50 * t62;
t56 = sin(qJ(5));
t105 = t56 * t58;
t91 = t52 * t105;
t142 = pkin(5) * t91 + t58 * t112;
t104 = t56 * t61;
t90 = t52 * t104;
t141 = pkin(5) * t90 + t61 * t112;
t107 = t52 * t61;
t40 = t50 * qJ(4);
t140 = pkin(3) * t107 + t61 * t40;
t138 = -m(7) * pkin(5) - mrSges(6,1);
t115 = t50 * t56;
t137 = pkin(5) * t115 - t52 * t62;
t113 = t50 * t61;
t96 = qJ(4) * t52;
t30 = t61 * t96;
t133 = -m(6) * t30 - mrSges(6,1) * t90 - mrSges(5,2) * t113 - mrSges(5,3) * t107 + t146 * t61;
t114 = t50 * t58;
t28 = t58 * t96;
t132 = -m(6) * t28 - mrSges(6,1) * t91 - mrSges(5,2) * t114 + (-mrSges(5,3) * t52 + t146) * t58;
t57 = sin(qJ(2));
t60 = cos(qJ(2));
t75 = t60 * mrSges(3,1) - t57 * mrSges(3,2);
t131 = -m(3) * pkin(1) - mrSges(2,1) + t144 - t75;
t42 = t52 * mrSges(7,3);
t130 = -mrSges(6,1) * t115 - t52 * mrSges(6,3) + t144 - t42 + (-mrSges(7,1) * t49 - mrSges(7,2) * t51 - t118) * t50;
t124 = pkin(2) * t57;
t77 = m(6) * (-pkin(3) - pkin(9)) - mrSges(6,3);
t68 = t77 * t50;
t129 = m(6) * t124 + t50 * mrSges(7,3) - t68 - t143 * (-pkin(3) * t50 - t124);
t123 = pkin(5) * t59;
t63 = -pkin(8) - pkin(7);
t128 = -m(6) * (pkin(4) - t63) - m(7) * (pkin(4) + t123) - mrSges(5,1) + mrSges(2,2) - mrSges(4,3) - m(3) * pkin(7) - mrSges(3,3);
t109 = t51 * t61;
t5 = t109 * t50 - t49 * t58;
t110 = t51 * t58;
t6 = t113 * t49 + t110;
t126 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t7 = t110 * t50 + t49 * t61;
t8 = -t114 * t49 + t109;
t125 = t7 * mrSges(7,1) + t8 * mrSges(7,2);
t120 = g(3) * t52;
t46 = t52 * pkin(3);
t53 = t60 * pkin(2);
t103 = t58 * t59;
t102 = t59 * t61;
t97 = t46 + t40;
t87 = t53 + t97;
t48 = t53 + pkin(1);
t31 = t61 * t48;
t83 = -t58 * t63 + t31;
t82 = -t48 - t40;
t81 = -pkin(3) * t114 + t28;
t80 = -pkin(3) * t113 + t30;
t72 = mrSges(4,1) * t50 + mrSges(4,2) * t52;
t69 = t97 + t137;
t9 = t102 * t50 - t105;
t11 = t103 * t50 + t104;
t45 = t52 * pkin(9);
t25 = mrSges(7,2) * t117;
t12 = -t105 * t50 + t102;
t10 = t104 * t50 + t103;
t1 = [(-m(4) * t83 - m(6) * (pkin(9) * t107 + t140 + t31) - t10 * mrSges(6,1) - t9 * mrSges(6,2) - mrSges(6,3) * t107 - t6 * mrSges(7,1) - t5 * mrSges(7,2) - t143 * (t83 + t140) + t128 * t58 + (-m(7) * t137 + t131 - t42) * t61) * g(2) + (-t12 * mrSges(6,1) - t8 * mrSges(7,1) + t11 * mrSges(6,2) + t7 * mrSges(7,2) + ((m(4) + t143) * t63 + t128) * t61 + (m(4) * t48 - m(5) * (t82 - t46) - m(6) * t82 - m(7) * (-t48 + (-pkin(5) * t56 - qJ(4)) * t50) + (-t77 - m(7) * (-pkin(3) + t62) + mrSges(7,3)) * t52 - t131) * t58) * g(1) (-m(5) * t28 - m(7) * (t28 + t142) + t129 * t58 + t132) * g(2) + (-m(5) * t30 - m(7) * (t30 + t141) + t129 * t61 + t133) * g(1) + (-t75 - m(4) * t53 - m(5) * t87 - m(6) * (t45 + t87) - m(7) * (t53 + t69) + t130) * g(3) + t134 * (m(4) * t124 + mrSges(3,1) * t57 + mrSges(3,2) * t60 + t72) t134 * t72 + (-m(5) * t81 - t58 * t68 - m(7) * (t81 + t142) + mrSges(7,3) * t114 + t132) * g(2) + (-m(5) * t80 - t61 * t68 - m(7) * (t80 + t141) + mrSges(7,3) * t113 + t133) * g(1) + (-m(5) * t97 - m(6) * (t45 + t97) - m(7) * t69 + t130) * g(3) (-t134 * t50 + t120) * (m(6) + t143) -(-mrSges(6,1) * t59 + mrSges(6,2) * t56) * t120 - g(3) * (t25 + (-m(7) * t123 - mrSges(7,1) * t51) * t52) + (-t12 * mrSges(6,2) + t138 * t11 - t125) * g(2) + (t10 * mrSges(6,2) + t138 * t9 - t126) * g(1), -g(1) * t126 - g(2) * t125 - g(3) * (-mrSges(7,1) * t111 + t25)];
taug  = t1(:);
