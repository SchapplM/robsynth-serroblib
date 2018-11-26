% Calculate Gravitation load on the joints for
% S6PRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
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
% Datum: 2018-11-23 15:27
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRRPR8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR8_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:26:48
% EndTime: 2018-11-23 15:26:50
% DurationCPUTime: 1.47s
% Computational Cost: add. (3626->156), mult. (3761->215), div. (0->0), fcn. (3658->22), ass. (0->104)
t168 = m(6) + m(7);
t79 = sin(qJ(6));
t82 = cos(qJ(6));
t173 = -t79 * mrSges(7,1) - t82 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t177 = -qJ(5) * t168 + t173;
t176 = mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t161 = m(7) * pkin(11) + t176;
t172 = pkin(4) * t168 + t161;
t129 = pkin(7) + qJ(3);
t112 = sin(t129) / 0.2e1;
t130 = pkin(7) - qJ(3);
t117 = sin(t130);
t163 = t112 - t117 / 0.2e1;
t131 = pkin(6) + qJ(2);
t116 = cos(t131) / 0.2e1;
t132 = pkin(6) - qJ(2);
t121 = cos(t132);
t73 = t116 - t121 / 0.2e1;
t84 = cos(qJ(3));
t113 = sin(t131) / 0.2e1;
t118 = sin(t132);
t97 = t113 + t118 / 0.2e1;
t171 = t163 * t97 - t73 * t84;
t134 = sin(pkin(12));
t136 = cos(pkin(12));
t72 = t113 - t118 / 0.2e1;
t85 = cos(qJ(2));
t64 = -t134 * t72 + t136 * t85;
t100 = t121 / 0.2e1 + t116;
t154 = sin(qJ(2));
t89 = t134 * t100 + t136 * t154;
t170 = -t163 * t89 + t64 * t84;
t62 = t134 * t85 + t136 * t72;
t88 = -t136 * t100 + t134 * t154;
t169 = -t163 * t88 + t62 * t84;
t80 = sin(qJ(4));
t139 = qJ(5) * t80;
t81 = sin(qJ(3));
t135 = sin(pkin(6));
t95 = t112 + t117 / 0.2e1;
t92 = t95 * t135;
t119 = cos(t129);
t114 = t119 / 0.2e1;
t120 = cos(t130);
t115 = t120 / 0.2e1;
t99 = t115 + t114;
t22 = t136 * t92 + t62 * t81 + t88 * t99;
t83 = cos(qJ(4));
t153 = t22 * t83;
t166 = -pkin(4) * t153 - t22 * t139;
t25 = -t134 * t92 + t64 * t81 + t89 * t99;
t152 = t25 * t83;
t165 = -pkin(4) * t152 - t25 * t139;
t138 = cos(pkin(6));
t40 = -t138 * t95 - t73 * t81 - t97 * t99;
t151 = t40 * t83;
t164 = -pkin(4) * t151 - t40 * t139;
t159 = -t82 * mrSges(7,1) + t79 * mrSges(7,2) - mrSges(6,1) + mrSges(4,2) - mrSges(5,3);
t158 = -t173 * t80 + t176 * t83 + mrSges(4,1);
t157 = -m(7) * (pkin(5) + pkin(10)) + t159;
t78 = sin(pkin(7));
t147 = t78 * t80;
t146 = t78 * t83;
t34 = -t163 * t62 - t88 * t84;
t60 = t88 * pkin(2);
t142 = t34 * pkin(3) - t60;
t36 = -t163 * t64 - t89 * t84;
t61 = t89 * pkin(2);
t141 = t36 * pkin(3) - t61;
t45 = t163 * t73 + t97 * t84;
t71 = t97 * pkin(2);
t140 = t45 * pkin(3) + t71;
t137 = cos(pkin(7));
t133 = -m(5) - t168;
t125 = m(4) - t133;
t20 = t22 * pkin(3);
t98 = t114 - t120 / 0.2e1;
t93 = t98 * t135;
t24 = t136 * t93 + t169;
t124 = pkin(10) * t24 - t20;
t21 = t25 * pkin(3);
t27 = -t134 * t93 + t170;
t123 = pkin(10) * t27 - t21;
t39 = t40 * pkin(3);
t42 = -t138 * t98 + t171;
t122 = pkin(10) * t42 - t39;
t107 = t137 * t135;
t106 = t115 - t119 / 0.2e1;
t102 = t106 * t135;
t101 = -mrSges(3,2) + (t125 * pkin(9) + mrSges(4,3)) * t78;
t91 = -m(7) * pkin(5) + t133 * pkin(10) + t159;
t90 = t138 * t137 - t97 * t78;
t87 = t134 * t107 + t89 * t78;
t86 = -t136 * t107 + t88 * t78;
t41 = t138 * t106 + t171;
t32 = -t73 * t147 + t45 * t83;
t26 = t134 * t102 + t170;
t23 = -t136 * t102 + t169;
t14 = t41 * t80 - t90 * t83;
t12 = t147 * t64 + t36 * t83;
t10 = t147 * t62 + t34 * t83;
t5 = t26 * t80 - t87 * t83;
t3 = t23 * t80 - t86 * t83;
t1 = [(-m(2) - m(3) - t125) * g(3) (-t97 * mrSges(3,1) - m(4) * t71 - t45 * mrSges(4,1) - m(5) * t140 + t101 * t73 - t161 * t32 + t177 * (t73 * t146 + t45 * t80) + t91 * (-t73 * t99 + t97 * t81) - t168 * (t32 * pkin(4) + t140)) * g(3) + (t88 * mrSges(3,1) + m(4) * t60 - t34 * mrSges(4,1) - m(5) * t142 - t101 * t62 - t161 * t10 + t177 * (-t146 * t62 + t34 * t80) + t91 * (t62 * t99 - t88 * t81) - t168 * (t10 * pkin(4) + t142)) * g(2) + (t89 * mrSges(3,1) + m(4) * t61 - t36 * mrSges(4,1) - m(5) * t141 - t101 * t64 - t161 * t12 + t177 * (-t146 * t64 + t36 * t80) + t91 * (t64 * t99 - t89 * t81) - t168 * (t12 * pkin(4) + t141)) * g(1) (-m(5) * t122 - m(6) * (t122 + t164) - m(7) * (-pkin(11) * t151 + t164 - t39) + t157 * t42 + t158 * t40) * g(3) + (-m(5) * t124 - m(6) * (t124 + t166) - m(7) * (-pkin(11) * t153 + t166 - t20) + t157 * t24 + t158 * t22) * g(2) + (-m(5) * t123 - m(6) * (t123 + t165) - m(7) * (-pkin(11) * t152 + t165 - t21) + t157 * t27 + t158 * t25) * g(1) (t177 * (t41 * t83 + t90 * t80) + t172 * t14) * g(3) + (t177 * (t23 * t83 + t86 * t80) + t172 * t3) * g(2) + (t177 * (t26 * t83 + t87 * t80) + t172 * t5) * g(1), t168 * (-g(1) * t5 - g(2) * t3 - g(3) * t14) -g(1) * ((-t25 * t79 + t5 * t82) * mrSges(7,1) + (-t25 * t82 - t5 * t79) * mrSges(7,2)) - g(2) * ((-t22 * t79 + t3 * t82) * mrSges(7,1) + (-t22 * t82 - t3 * t79) * mrSges(7,2)) - g(3) * ((t14 * t82 - t40 * t79) * mrSges(7,1) + (-t14 * t79 - t40 * t82) * mrSges(7,2))];
taug  = t1(:);
