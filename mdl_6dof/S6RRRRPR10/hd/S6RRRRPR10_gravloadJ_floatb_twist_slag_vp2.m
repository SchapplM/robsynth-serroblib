% Calculate Gravitation load on the joints for
% S6RRRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
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
% Datum: 2018-11-23 18:20
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR10_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:19:18
% EndTime: 2018-11-23 18:19:19
% DurationCPUTime: 1.31s
% Computational Cost: add. (1645->163), mult. (1737->205), div. (0->0), fcn. (1712->16), ass. (0->86)
t86 = sin(qJ(6));
t90 = cos(qJ(6));
t170 = mrSges(7,1) * t86 + t90 * mrSges(7,2);
t173 = mrSges(6,3) - mrSges(5,2);
t174 = t170 + t173;
t172 = mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t157 = cos(qJ(1));
t133 = pkin(6) + qJ(2);
t112 = sin(t133) / 0.2e1;
t134 = pkin(6) - qJ(2);
t117 = sin(t134);
t64 = t112 - t117 / 0.2e1;
t89 = sin(qJ(1));
t92 = cos(qJ(2));
t104 = t157 * t92 - t89 * t64;
t84 = sin(pkin(6));
t147 = t84 * t89;
t87 = sin(qJ(3));
t91 = cos(qJ(3));
t27 = -t104 * t87 + t91 * t147;
t113 = cos(t133) / 0.2e1;
t118 = cos(t134);
t65 = t113 - t118 / 0.2e1;
t85 = cos(pkin(6));
t169 = t65 * t87 + t85 * t91;
t168 = m(6) + m(7);
t83 = qJ(3) + qJ(4);
t80 = sin(t83);
t81 = cos(t83);
t45 = -t65 * t80 - t85 * t81;
t46 = -t65 * t81 + t80 * t85;
t167 = t172 * t45 - t174 * t46;
t25 = t104 * t80 - t147 * t81;
t26 = t104 * t81 + t147 * t80;
t166 = t172 * t25 - t174 * t26;
t105 = t157 * t64 + t89 * t92;
t127 = t84 * t157;
t21 = t105 * t80 + t127 * t81;
t22 = t105 * t81 - t80 * t127;
t165 = t172 * t21 - t174 * t22;
t97 = -m(4) * pkin(9) - m(7) * pkin(5) - mrSges(6,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t163 = mrSges(7,1) * t90 - mrSges(7,2) * t86 - t97;
t164 = m(4) * pkin(2) + t91 * mrSges(4,1) - t87 * mrSges(4,2) + t172 * t81 + t174 * t80 + mrSges(3,1);
t162 = pkin(11) * t25;
t161 = pkin(11) * t45;
t160 = t21 * pkin(11);
t88 = sin(qJ(2));
t95 = t118 / 0.2e1 + t113;
t51 = -t157 * t95 + t88 * t89;
t155 = t51 * t81;
t54 = t157 * t88 + t89 * t95;
t154 = t54 * t81;
t63 = t112 + t117 / 0.2e1;
t152 = t63 * t81;
t79 = pkin(3) * t91 + pkin(2);
t93 = -pkin(10) - pkin(9);
t142 = -t105 * t93 - t51 * t79;
t141 = -t104 * t93 - t54 * t79;
t138 = t63 * t79 + t65 * t93;
t137 = t157 * pkin(1) + pkin(8) * t147;
t136 = qJ(5) * t80;
t131 = t87 * t147;
t126 = -pkin(1) * t89 + pkin(8) * t127;
t71 = t87 * t127;
t125 = -t105 * t91 + t71;
t124 = -t21 * pkin(4) + t22 * qJ(5);
t123 = -t25 * pkin(4) + qJ(5) * t26;
t122 = -t45 * pkin(4) + qJ(5) * t46;
t121 = -pkin(4) * t155 - t51 * t136 + t142;
t120 = -pkin(4) * t154 - t54 * t136 + t141;
t119 = pkin(4) * t152 + t63 * t136 + t138;
t116 = t27 * pkin(3);
t115 = t169 * pkin(3);
t114 = pkin(3) * t131 + t104 * t79 - t54 * t93 + t137;
t109 = -m(7) * pkin(11) - t172;
t107 = -t168 * qJ(5) - t173;
t106 = pkin(3) * t71 - t105 * t79 + t51 * t93 + t126;
t100 = t105 * t87 + t127 * t91;
t99 = t116 + t123;
t98 = t115 + t122;
t96 = t100 * pkin(3);
t94 = t124 - t96;
t28 = t104 * t91 + t131;
t2 = t25 * t86 + t54 * t90;
t1 = t25 * t90 - t54 * t86;
t3 = [(-t157 * mrSges(2,1) - m(3) * t137 - t104 * mrSges(3,1) - m(4) * (pkin(2) * t104 + t137) - t28 * mrSges(4,1) - t27 * mrSges(4,2) - m(5) * t114 - t2 * mrSges(7,1) - t1 * mrSges(7,2) + (-mrSges(3,3) * t84 + mrSges(2,2)) * t89 + t107 * t25 + t109 * t26 + t97 * t54 - t168 * (t26 * pkin(4) + t114)) * g(2) + (t89 * mrSges(2,1) + t157 * mrSges(2,2) - m(3) * t126 + t105 * mrSges(3,1) - mrSges(3,3) * t127 - m(4) * (-pkin(2) * t105 + t126) - t125 * mrSges(4,1) - t100 * mrSges(4,2) - m(5) * t106 - t109 * t22 - (t107 - t170) * t21 + t163 * t51 + t168 * (pkin(4) * t22 - t106)) * g(1) (-m(5) * t138 - m(6) * t119 - m(7) * (pkin(11) * t152 + t119) + t163 * t65 - t164 * t63) * g(3) + (-m(5) * t142 - m(6) * t121 - m(7) * (-pkin(11) * t155 + t121) - t163 * t105 + t164 * t51) * g(2) + (-m(5) * t141 - m(6) * t120 - m(7) * (-pkin(11) * t154 + t120) - t163 * t104 + t164 * t54) * g(1) (-t169 * mrSges(4,1) - (t65 * t91 - t85 * t87) * mrSges(4,2) - m(5) * t115 - m(6) * t98 - m(7) * (t98 - t161) + t167) * g(3) + (t100 * mrSges(4,1) - t125 * mrSges(4,2) + m(5) * t96 - m(6) * t94 - m(7) * (t94 - t160) + t165) * g(2) + (-mrSges(4,1) * t27 + mrSges(4,2) * t28 - m(5) * t116 - m(6) * t99 - m(7) * (t99 - t162) + t166) * g(1) (-m(6) * t122 - m(7) * (t122 - t161) + t167) * g(3) + (-m(6) * t124 - m(7) * (t124 - t160) + t165) * g(2) + (-m(6) * t123 - m(7) * (t123 - t162) + t166) * g(1), t168 * (-g(1) * t25 - g(2) * t21 - g(3) * t45) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t21 * t90 - t51 * t86) * mrSges(7,1) + (-t21 * t86 - t51 * t90) * mrSges(7,2)) - g(3) * ((t45 * t90 + t63 * t86) * mrSges(7,1) + (-t45 * t86 + t63 * t90) * mrSges(7,2))];
taug  = t3(:);
