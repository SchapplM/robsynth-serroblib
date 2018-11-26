% Calculate Gravitation load on the joints for
% S6RRRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2018-11-23 18:34
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP10_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:33:52
% EndTime: 2018-11-23 18:33:54
% DurationCPUTime: 1.65s
% Computational Cost: add. (2056->174), mult. (2350->226), div. (0->0), fcn. (2391->16), ass. (0->95)
t189 = mrSges(6,2) - mrSges(7,3);
t126 = -m(7) * qJ(6) + t189;
t190 = mrSges(6,1) + mrSges(7,1);
t194 = m(7) * pkin(5) + t190;
t100 = cos(qJ(4));
t96 = sin(qJ(4));
t191 = m(5) * pkin(3) + t100 * mrSges(5,1) - t96 * mrSges(5,2) + mrSges(4,1);
t114 = -m(5) * pkin(10) + mrSges(4,2) - mrSges(7,2) - mrSges(5,3) - mrSges(6,3);
t94 = qJ(4) + qJ(5);
t91 = sin(t94);
t92 = cos(t94);
t188 = -t126 * t91 + t194 * t92 + t191;
t101 = cos(qJ(3));
t151 = pkin(6) + qJ(2);
t124 = cos(t151) / 0.2e1;
t152 = pkin(6) - qJ(2);
t134 = cos(t152);
t82 = t124 - t134 / 0.2e1;
t95 = cos(pkin(6));
t97 = sin(qJ(3));
t66 = -t101 * t82 + t95 * t97;
t132 = sin(t151);
t122 = t132 / 0.2e1;
t133 = sin(t152);
t123 = t133 / 0.2e1;
t81 = t122 + t123;
t187 = -t100 * t81 - t66 * t96;
t153 = sin(pkin(6));
t99 = sin(qJ(1));
t138 = t99 * t153;
t115 = t122 - t133 / 0.2e1;
t102 = cos(qJ(2));
t176 = cos(qJ(1));
t141 = t176 * t102;
t71 = -t115 * t99 + t141;
t45 = t71 * t101 + t138 * t97;
t107 = t134 / 0.2e1 + t124;
t98 = sin(qJ(2));
t70 = t107 * t99 + t176 * t98;
t17 = t100 * t70 - t45 * t96;
t125 = t176 * t153;
t154 = t99 * t102;
t68 = t115 * t176 + t154;
t41 = t101 * t68 - t97 * t125;
t67 = -t107 * t176 + t98 * t99;
t186 = t100 * t67 - t41 * t96;
t11 = t41 * t91 - t67 * t92;
t12 = t41 * t92 + t67 * t91;
t185 = m(4) + m(5);
t184 = m(6) + m(7);
t166 = mrSges(3,2) - mrSges(4,3);
t117 = t96 * mrSges(5,1) + t100 * mrSges(5,2);
t182 = -m(5) * pkin(9) - t117 + t166;
t180 = -pkin(9) * (t184 + t185) + t166;
t178 = t191 * t101 - t114 * t97 + mrSges(3,1);
t177 = pkin(4) * t96;
t170 = t67 * t96;
t169 = t70 * t96;
t162 = t176 * pkin(1) + pkin(8) * t138;
t158 = t101 * t70;
t103 = -pkin(11) - pkin(10);
t157 = t103 * t97;
t156 = t67 * t101;
t155 = t81 * t101;
t150 = t71 * pkin(2) + t162;
t148 = t11 * t190 + t189 * t12;
t15 = t45 * t91 - t70 * t92;
t16 = t45 * t92 + t70 * t91;
t147 = t15 * t190 + t189 * t16;
t145 = -pkin(1) * t99 + pkin(8) * t125;
t61 = t67 * pkin(2);
t106 = t123 - t132 / 0.2e1;
t69 = -t106 * t176 + t154;
t144 = t69 * pkin(9) - t61;
t63 = t70 * pkin(2);
t72 = t106 * t99 + t141;
t143 = t72 * pkin(9) - t63;
t80 = t81 * pkin(2);
t142 = -t82 * pkin(9) + t80;
t140 = -t11 * pkin(5) + qJ(6) * t12;
t139 = -t101 * t125 - t68 * t97;
t137 = -t15 * pkin(5) + qJ(6) * t16;
t28 = t66 * t91 + t81 * t92;
t29 = t66 * t92 - t81 * t91;
t136 = -t28 * pkin(5) + qJ(6) * t29;
t135 = t189 * t29 + t190 * t28;
t130 = t186 * pkin(4);
t129 = t17 * pkin(4);
t128 = t187 * pkin(4);
t127 = t68 * pkin(2) - t145;
t90 = pkin(4) * t100 + pkin(3);
t65 = t101 * t95 + t82 * t97;
t44 = -t101 * t138 + t71 * t97;
t18 = t100 * t45 + t169;
t1 = [(-t176 * mrSges(2,1) - m(3) * t162 - t71 * mrSges(3,1) - m(4) * t150 - t45 * mrSges(4,1) - m(5) * (pkin(3) * t45 + t150) - t18 * mrSges(5,1) - t17 * mrSges(5,2) + (-mrSges(3,3) * t153 + mrSges(2,2)) * t99 + t180 * t70 - t194 * t16 + t126 * t15 + t114 * t44 - t184 * (pkin(4) * t169 - t44 * t103 + t45 * t90 + t150)) * g(2) + (t99 * mrSges(2,1) + t176 * mrSges(2,2) - m(3) * t145 + t68 * mrSges(3,1) - mrSges(3,3) * t125 + t191 * t41 + (t117 - t180) * t67 + t194 * t12 - t126 * t11 + t114 * t139 + t185 * t127 + t184 * (pkin(4) * t170 + t103 * t139 + t41 * t90 + t127)) * g(1) (-m(4) * t142 - m(5) * t80 - t184 * (t90 * t155 - t157 * t81 - t177 * t82 + t142) - t194 * (t155 * t92 - t82 * t91) + t126 * (t155 * t91 + t82 * t92) - t182 * t82 - t178 * t81) * g(3) + (-m(4) * t144 + m(5) * t61 - t184 * (-t90 * t156 + t157 * t67 + t177 * t69 + t144) - t194 * (-t156 * t92 + t69 * t91) + t126 * (-t156 * t91 - t69 * t92) + t182 * t69 + t178 * t67) * g(2) + (-m(4) * t143 + m(5) * t63 - t184 * (t157 * t70 - t90 * t158 + t177 * t72 + t143) - t194 * (-t158 * t92 + t72 * t91) + t126 * (-t158 * t91 - t72 * t92) + t182 * t72 + t178 * t70) * g(1) (-t184 * (-t66 * t103 + t65 * t90) + t114 * t66 - t188 * t65) * g(3) + (-t184 * (-t41 * t103 + t139 * t90) + t114 * t41 - t188 * t139) * g(2) + (-t184 * (-t45 * t103 - t44 * t90) + t114 * t45 + t188 * t44) * g(1) (-t187 * mrSges(5,1) - (-t100 * t66 + t81 * t96) * mrSges(5,2) - m(6) * t128 - m(7) * (t128 + t136) + t135) * g(3) + (-t186 * mrSges(5,1) - (-t100 * t41 - t170) * mrSges(5,2) - m(6) * t130 - m(7) * (t130 + t140) + t148) * g(2) + (-t17 * mrSges(5,1) + t18 * mrSges(5,2) - m(6) * t129 - m(7) * (t129 + t137) + t147) * g(1) (-m(7) * t136 + t135) * g(3) + (-m(7) * t140 + t148) * g(2) + (-m(7) * t137 + t147) * g(1) (-g(1) * t15 - g(2) * t11 - g(3) * t28) * m(7)];
taug  = t1(:);
