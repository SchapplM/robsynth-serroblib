% Calculate Gravitation load on the joints for
% S6RRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 18:16
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:16:03
% EndTime: 2018-11-23 18:16:04
% DurationCPUTime: 1.24s
% Computational Cost: add. (1619->172), mult. (1550->217), div. (0->0), fcn. (1503->18), ass. (0->90)
t177 = -mrSges(7,3) + mrSges(6,2);
t91 = sin(qJ(6));
t95 = cos(qJ(6));
t175 = -mrSges(7,1) * t95 + mrSges(7,2) * t91 - mrSges(6,1);
t176 = m(7) * pkin(5) - t175;
t126 = -m(7) * pkin(11) + t177;
t172 = m(6) + m(7);
t89 = sin(pkin(6));
t94 = sin(qJ(1));
t154 = t89 * t94;
t139 = pkin(6) + qJ(2);
t128 = sin(t139);
t119 = t128 / 0.2e1;
t140 = pkin(6) - qJ(2);
t129 = sin(t140);
t111 = t119 - t129 / 0.2e1;
t97 = cos(qJ(2));
t98 = cos(qJ(1));
t151 = t98 * t97;
t52 = -t111 * t94 + t151;
t88 = qJ(3) + qJ(4);
t83 = sin(t88);
t84 = cos(t88);
t23 = t84 * t154 - t52 * t83;
t121 = cos(t139) / 0.2e1;
t130 = cos(t140);
t66 = t121 - t130 / 0.2e1;
t90 = cos(pkin(6));
t174 = t66 * t83 + t84 * t90;
t171 = -m(5) * pkin(3) - mrSges(4,1);
t99 = -pkin(10) - pkin(9);
t102 = -m(4) * pkin(9) + m(5) * t99 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t170 = mrSges(7,1) * t91 + mrSges(7,2) * t95 - t102;
t96 = cos(qJ(3));
t85 = t96 * pkin(3);
t81 = t85 + pkin(2);
t92 = sin(qJ(3));
t169 = m(4) * pkin(2) + m(5) * t81 + t96 * mrSges(4,1) + t84 * mrSges(5,1) - t92 * mrSges(4,2) - t83 * mrSges(5,2) + mrSges(3,1);
t82 = pkin(12) + t88;
t78 = sin(t82);
t79 = cos(t82);
t38 = t66 * t78 + t79 * t90;
t39 = -t66 * t79 + t78 * t90;
t168 = -t174 * mrSges(5,1) - (t66 * t84 - t83 * t90) * mrSges(5,2) + t177 * t39 + t175 * t38;
t17 = -t154 * t79 + t52 * t78;
t18 = t154 * t78 + t52 * t79;
t24 = t154 * t83 + t52 * t84;
t167 = -t23 * mrSges(5,1) + t24 * mrSges(5,2) - t175 * t17 + t177 * t18;
t153 = t89 * t98;
t152 = t94 * t97;
t49 = t111 * t98 + t152;
t110 = -t153 * t84 - t49 * t83;
t13 = -t79 * t153 - t49 * t78;
t14 = -t78 * t153 + t49 * t79;
t70 = t83 * t153;
t166 = -t110 * mrSges(5,1) - (-t49 * t84 + t70) * mrSges(5,2) + t177 * t14 + t175 * t13;
t164 = -t126 * t78 + t176 * t79 + t169;
t161 = pkin(3) * t92;
t68 = pkin(4) * t83 + t161;
t69 = pkin(4) * t84 + t85;
t146 = t69 * t154 - t52 * t68;
t142 = t66 * t68 + t90 * t69;
t141 = t98 * pkin(1) + pkin(8) * t154;
t134 = t94 * pkin(1) - pkin(8) * t153;
t133 = t13 * pkin(5) + pkin(11) * t14;
t132 = -t17 * pkin(5) + pkin(11) * t18;
t131 = t38 * pkin(5) + pkin(11) * t39;
t127 = -m(5) * t161 - mrSges(3,3);
t125 = t23 * pkin(4);
t124 = t174 * pkin(4);
t123 = -t153 * t69 - t49 * t68;
t101 = t130 / 0.2e1 + t121;
t93 = sin(qJ(2));
t51 = t101 * t94 + t98 * t93;
t64 = pkin(2) + t69;
t87 = -qJ(5) + t99;
t122 = t68 * t154 - t51 * t87 + t52 * t64 + t141;
t120 = t129 / 0.2e1;
t112 = t120 - t128 / 0.2e1;
t25 = t154 * t96 - t52 * t92;
t104 = t110 * pkin(4);
t72 = t92 * t153;
t65 = t119 + t120;
t53 = t112 * t94 + t151;
t50 = -t112 * t98 + t152;
t48 = -t101 * t98 + t93 * t94;
t26 = t154 * t92 + t52 * t96;
t2 = t18 * t95 + t51 * t91;
t1 = -t18 * t91 + t51 * t95;
t3 = [(-t98 * mrSges(2,1) - m(3) * t141 - t52 * mrSges(3,1) - m(4) * (pkin(2) * t52 + t141) - t26 * mrSges(4,1) - t25 * mrSges(4,2) - m(5) * (t52 * t81 + t141) - t24 * mrSges(5,1) - t23 * mrSges(5,2) - m(6) * t122 - t18 * mrSges(6,1) - m(7) * (pkin(5) * t18 + t122) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + (t127 * t89 + mrSges(2,2)) * t94 + t126 * t17 + t102 * t51) * g(2) + (t94 * mrSges(2,1) - t72 * mrSges(4,1) - t70 * mrSges(5,1) + t169 * t49 + (mrSges(2,2) + (-mrSges(4,2) * t96 - mrSges(5,2) * t84 + t127) * t89) * t98 + t126 * t13 + t176 * t14 + t170 * t48 + (m(3) + m(4) + m(5)) * t134 + t172 * (-t68 * t153 - t48 * t87 + t49 * t64 + t134)) * g(1) (-t172 * (t65 * t64 + t66 * t87) + t170 * t66 - t164 * t65) * g(3) + (-t172 * (-t48 * t64 - t50 * t87) - t170 * t50 + t164 * t48) * g(2) + (-t172 * (-t51 * t64 - t53 * t87) - t170 * t53 + t164 * t51) * g(1) (-(t66 * t96 - t90 * t92) * mrSges(4,2) - m(6) * t142 - m(7) * (t131 + t142) + t171 * (t66 * t92 + t90 * t96) + t168) * g(3) + (-(-t49 * t96 + t72) * mrSges(4,2) - m(6) * t123 - m(7) * (t123 + t133) + t171 * (-t153 * t96 - t49 * t92) + t166) * g(2) + (mrSges(4,2) * t26 - m(6) * t146 - m(7) * (t132 + t146) + t171 * t25 + t167) * g(1) (-m(6) * t124 - m(7) * (t124 + t131) + t168) * g(3) + (-m(6) * t104 - m(7) * (t104 + t133) + t166) * g(2) + (-m(6) * t125 - m(7) * (t125 + t132) + t167) * g(1), t172 * (-g(1) * t51 - g(2) * t48 + g(3) * t65) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t14 * t91 + t48 * t95) * mrSges(7,1) + (-t14 * t95 - t48 * t91) * mrSges(7,2)) - g(3) * ((-t39 * t91 - t65 * t95) * mrSges(7,1) + (-t39 * t95 + t65 * t91) * mrSges(7,2))];
taug  = t3(:);
