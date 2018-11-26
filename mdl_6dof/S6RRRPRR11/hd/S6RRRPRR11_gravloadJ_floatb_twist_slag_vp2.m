% Calculate Gravitation load on the joints for
% S6RRRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2018-11-23 17:59
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR11_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR11_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR11_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR11_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:59:01
% EndTime: 2018-11-23 17:59:02
% DurationCPUTime: 1.77s
% Computational Cost: add. (1931->161), mult. (2476->203), div. (0->0), fcn. (2585->16), ass. (0->88)
t175 = m(6) + m(7);
t171 = -mrSges(3,2) + mrSges(5,2) + mrSges(4,3) - mrSges(6,3) + t175 * (pkin(9) - pkin(10));
t141 = cos(qJ(1));
t82 = sin(qJ(2));
t83 = sin(qJ(1));
t126 = pkin(6) + qJ(2);
t101 = cos(t126) / 0.2e1;
t127 = pkin(6) - qJ(2);
t111 = cos(t127);
t89 = t111 / 0.2e1 + t101;
t55 = -t141 * t89 + t82 * t83;
t128 = sin(pkin(6));
t102 = t141 * t128;
t110 = sin(t127);
t109 = sin(t126);
t74 = t109 / 0.2e1;
t132 = t74 - t110 / 0.2e1;
t87 = cos(qJ(2));
t133 = t83 * t87;
t56 = t132 * t141 + t133;
t81 = sin(qJ(3));
t86 = cos(qJ(3));
t33 = t102 * t86 + t56 * t81;
t34 = -t81 * t102 + t56 * t86;
t80 = sin(qJ(5));
t85 = cos(qJ(5));
t6 = t33 * t80 + t34 * t85;
t79 = sin(qJ(6));
t84 = cos(qJ(6));
t174 = t55 * t84 + t6 * t79;
t173 = t55 * t79 - t6 * t84;
t161 = mrSges(4,1) + mrSges(5,1);
t159 = mrSges(4,2) - mrSges(5,3);
t172 = -t79 * mrSges(7,1) - t84 * mrSges(7,2) + t171;
t155 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t160 = m(7) * pkin(5) + t84 * mrSges(7,1) - t79 * mrSges(7,2) + mrSges(6,1);
t164 = (-t80 * t86 + t81 * t85) * t155 + (t80 * t81 + t85 * t86) * t160 + mrSges(3,1) + t161 * t86 - t159 * t81;
t163 = t33 * t85 - t34 * t80;
t170 = t155 * t6 - t160 * t163;
t118 = t83 * t128;
t123 = t141 * t87;
t59 = -t132 * t83 + t123;
t37 = -t118 * t86 + t59 * t81;
t38 = t118 * t81 + t59 * t86;
t12 = t37 * t80 + t38 * t85;
t97 = -t37 * t85 + t38 * t80;
t169 = t12 * t155 + t160 * t97;
t129 = cos(pkin(6));
t69 = t101 - t111 / 0.2e1;
t53 = -t129 * t86 - t69 * t81;
t54 = t129 * t81 - t69 * t86;
t20 = t53 * t80 + t54 * t85;
t168 = t155 * t20 - t160 * (t53 * t85 - t54 * t80);
t130 = qJ(4) * t81;
t138 = t55 * t86;
t158 = -pkin(3) * t138 - t55 * t130;
t58 = t141 * t82 + t83 * t89;
t137 = t58 * t86;
t157 = -pkin(3) * t137 - t58 * t130;
t100 = t110 / 0.2e1;
t68 = t74 + t100;
t136 = t68 * t86;
t156 = pkin(3) * t136 + t68 * t130;
t143 = pkin(9) * t55;
t142 = pkin(9) * t58;
t131 = t141 * pkin(1) + pkin(8) * t118;
t125 = t59 * pkin(2) + t131;
t122 = -pkin(1) * t83 + pkin(8) * t102;
t49 = t55 * pkin(2);
t88 = t100 - t109 / 0.2e1;
t57 = -t141 * t88 + t133;
t121 = pkin(9) * t57 - t49;
t51 = t58 * pkin(2);
t60 = t83 * t88 + t123;
t120 = pkin(9) * t60 - t51;
t67 = t68 * pkin(2);
t119 = -pkin(9) * t69 + t67;
t117 = -t33 * pkin(3) + qJ(4) * t34;
t116 = -t37 * pkin(3) + qJ(4) * t38;
t115 = -t53 * pkin(3) + qJ(4) * t54;
t108 = -t56 * pkin(2) + t122;
t95 = t38 * pkin(3) + qJ(4) * t37 + t125;
t92 = t38 * pkin(4) + t95;
t91 = -pkin(3) * t34 - qJ(4) * t33 + t108;
t90 = -pkin(4) * t34 + t91;
t2 = t12 * t84 - t58 * t79;
t1 = -t12 * t79 - t58 * t84;
t3 = [(-t141 * mrSges(2,1) + t83 * mrSges(2,2) - m(3) * t131 - t59 * mrSges(3,1) - mrSges(3,3) * t118 - m(4) * (t125 + t142) - m(5) * (t95 + t142) - m(6) * t92 - t12 * mrSges(6,1) - m(7) * (pkin(5) * t12 + t92) - t2 * mrSges(7,1) - t1 * mrSges(7,2) - t161 * t38 + t159 * t37 + t155 * t97 - t171 * t58) * g(2) + (t83 * mrSges(2,1) + t141 * mrSges(2,2) - m(3) * t122 + t56 * mrSges(3,1) - mrSges(3,3) * t102 - m(4) * (t108 - t143) - m(5) * (t91 - t143) - m(6) * t90 + t6 * mrSges(6,1) - m(7) * (-pkin(5) * t6 + t90) - t173 * mrSges(7,1) - t174 * mrSges(7,2) + t155 * t163 + t161 * t34 - t159 * t33 + t171 * t55) * g(1) (-m(4) * t119 - m(5) * (t119 + t156) - t175 * (pkin(4) * t136 + t156 + t67) + t172 * t69 - t164 * t68) * g(3) + (-m(4) * t121 - m(5) * (t121 + t158) - t175 * (-pkin(4) * t138 + t158 - t49) - t172 * t57 + t164 * t55) * g(2) + (-m(4) * t120 - m(5) * (t120 + t157) - t175 * (-pkin(4) * t137 + t157 - t51) - t172 * t60 + t164 * t58) * g(1) (-m(5) * t115 + t159 * t54 + t161 * t53 - t175 * (-t53 * pkin(4) + t115) - t168) * g(3) + (-m(5) * t117 + t159 * t34 + t161 * t33 - t175 * (-t33 * pkin(4) + t117) - t170) * g(2) + (-m(5) * t116 - t175 * (-t37 * pkin(4) + t116) + t159 * t38 + t161 * t37 - t169) * g(1) (m(5) + t175) * (-g(1) * t37 - g(2) * t33 - g(3) * t53) g(1) * t169 + g(2) * t170 + g(3) * t168, -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * (-t174 * mrSges(7,1) + t173 * mrSges(7,2)) - g(3) * ((-t20 * t79 + t68 * t84) * mrSges(7,1) + (-t20 * t84 - t68 * t79) * mrSges(7,2))];
taug  = t3(:);
