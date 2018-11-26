% Calculate Gravitation load on the joints for
% S6RRRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
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
% Datum: 2018-11-23 18:09
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRRPP7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:08:52
% EndTime: 2018-11-23 18:08:54
% DurationCPUTime: 1.45s
% Computational Cost: add. (1868->144), mult. (2167->185), div. (0->0), fcn. (2191->16), ass. (0->77)
t114 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t170 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t81 = sin(qJ(4));
t85 = cos(qJ(4));
t172 = m(5) * pkin(3) + t85 * mrSges(5,1) - t81 * mrSges(5,2) + mrSges(4,1);
t79 = qJ(4) + pkin(11);
t76 = sin(t79);
t77 = cos(t79);
t156 = t114 * t77 + t170 * t76 + t172;
t160 = m(6) + m(7);
t173 = pkin(4) * t160;
t99 = -m(5) * pkin(10) + mrSges(4,2) - mrSges(7,2) - mrSges(5,3) - mrSges(6,3);
t75 = pkin(4) * t85 + pkin(3);
t80 = -qJ(5) - pkin(10);
t82 = sin(qJ(3));
t86 = cos(qJ(3));
t171 = t160 * (-t75 * t86 + t80 * t82) + t99 * t82 - mrSges(3,1) - t156 * t86;
t169 = -m(4) - t160;
t130 = cos(pkin(6));
t127 = pkin(6) + qJ(2);
t107 = cos(t127) / 0.2e1;
t128 = pkin(6) - qJ(2);
t117 = cos(t128);
t67 = t107 - t117 / 0.2e1;
t51 = t130 * t82 - t67 * t86;
t115 = sin(t127);
t105 = t115 / 0.2e1;
t116 = sin(t128);
t106 = t116 / 0.2e1;
t66 = t105 + t106;
t167 = t51 * t81 + t66 * t85;
t129 = sin(pkin(6));
t84 = sin(qJ(1));
t118 = t84 * t129;
t100 = t105 - t116 / 0.2e1;
t153 = cos(qJ(1));
t87 = cos(qJ(2));
t123 = t153 * t87;
t56 = -t100 * t84 + t123;
t30 = t118 * t82 + t56 * t86;
t83 = sin(qJ(2));
t91 = t117 / 0.2e1 + t107;
t55 = t153 * t83 + t84 * t91;
t7 = -t30 * t81 + t55 * t85;
t108 = t153 * t129;
t136 = t84 * t87;
t53 = t100 * t153 + t136;
t26 = -t82 * t108 + t53 * t86;
t52 = -t153 * t91 + t83 * t84;
t166 = -t26 * t81 + t52 * t85;
t1 = t26 * t76 - t52 * t77;
t165 = t26 * t77 + t52 * t76;
t103 = t81 * mrSges(5,1) + t85 * mrSges(5,2);
t135 = mrSges(3,2) - mrSges(4,3);
t163 = m(5) * pkin(9) + t114 * t76 - t170 * t77 + t81 * t173 + t103 - t135;
t161 = m(4) + m(5);
t157 = -pkin(9) * (t160 + t161) + t135;
t146 = t52 * t81;
t144 = t55 * t81;
t131 = t153 * pkin(1) + pkin(8) * t118;
t126 = t56 * pkin(2) + t131;
t122 = -pkin(1) * t84 + pkin(8) * t108;
t110 = t53 * pkin(2) - t122;
t25 = t108 * t86 + t53 * t82;
t90 = t106 - t115 / 0.2e1;
t65 = t66 * pkin(2);
t57 = t84 * t90 + t123;
t54 = -t153 * t90 + t136;
t50 = -t130 * t86 - t67 * t82;
t48 = t55 * pkin(2);
t46 = t52 * pkin(2);
t29 = -t118 * t86 + t56 * t82;
t13 = t51 * t76 + t66 * t77;
t8 = t30 * t85 + t144;
t6 = t30 * t77 + t55 * t76;
t5 = t30 * t76 - t55 * t77;
t2 = [(-t153 * mrSges(2,1) - m(3) * t131 - t56 * mrSges(3,1) - m(4) * t126 - t30 * mrSges(4,1) - m(5) * (pkin(3) * t30 + t126) - t8 * mrSges(5,1) - t7 * mrSges(5,2) + (-mrSges(3,3) * t129 + mrSges(2,2)) * t84 - t114 * t6 + t157 * t55 - t170 * t5 + t99 * t29 - t160 * (pkin(4) * t144 - t29 * t80 + t30 * t75 + t126)) * g(2) + (t84 * mrSges(2,1) + t153 * mrSges(2,2) - m(3) * t122 + t53 * mrSges(3,1) - mrSges(3,3) * t108 + t172 * t26 + (t103 - t157) * t52 + t114 * t165 + t170 * t1 - t99 * t25 + t160 * (pkin(4) * t146 - t25 * t80 + t26 * t75 + t110) + t161 * t110) * g(1) (-m(5) * t65 + t169 * (-pkin(9) * t67 + t65) + t163 * t67 + t171 * t66) * g(3) + (m(5) * t46 + t169 * (pkin(9) * t54 - t46) - t163 * t54 - t171 * t52) * g(2) + (m(5) * t48 + t169 * (pkin(9) * t57 - t48) - t163 * t57 - t171 * t55) * g(1) (-t160 * (-t50 * t75 - t51 * t80) + t99 * t51 + t156 * t50) * g(3) + (-t160 * (-t25 * t75 - t26 * t80) + t99 * t26 + t156 * t25) * g(2) + (-t160 * (-t29 * t75 - t30 * t80) + t99 * t30 + t156 * t29) * g(1) (t167 * mrSges(5,1) - (-t51 * t85 + t66 * t81) * mrSges(5,2) - t170 * (t51 * t77 - t66 * t76) + t114 * t13) * g(3) + (-t166 * mrSges(5,1) - (-t26 * t85 - t146) * mrSges(5,2) - t170 * t165 + t114 * t1) * g(2) + (-t7 * mrSges(5,1) + t8 * mrSges(5,2) + t114 * t5 - t170 * t6) * g(1) + (-g(1) * t7 - g(2) * t166 + g(3) * t167) * t173, t160 * (-g(1) * t29 - g(2) * t25 - g(3) * t50) (-g(1) * t5 - g(2) * t1 - g(3) * t13) * m(7)];
taug  = t2(:);
