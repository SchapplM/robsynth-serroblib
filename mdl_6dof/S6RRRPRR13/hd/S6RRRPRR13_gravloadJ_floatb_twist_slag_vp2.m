% Calculate Gravitation load on the joints for
% S6RRRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2018-11-23 18:01
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR13_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR13_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR13_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR13_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:00:48
% EndTime: 2018-11-23 18:00:50
% DurationCPUTime: 1.85s
% Computational Cost: add. (4116->181), mult. (4188->248), div. (0->0), fcn. (4099->24), ass. (0->89)
t106 = cos(qJ(3));
t108 = cos(qJ(1));
t160 = pkin(6) + qJ(2);
t144 = cos(t160) / 0.2e1;
t161 = pkin(6) - qJ(2);
t153 = cos(t161);
t130 = t153 / 0.2e1 + t144;
t174 = sin(qJ(2));
t175 = sin(qJ(1));
t124 = -t108 * t130 + t175 * t174;
t162 = sin(pkin(6));
t154 = t108 * t162;
t107 = cos(qJ(2));
t142 = sin(t160) / 0.2e1;
t151 = sin(t161);
t81 = t142 - t151 / 0.2e1;
t71 = t175 * t107 + t108 * t81;
t158 = pkin(7) + qJ(3);
t141 = sin(t158) / 0.2e1;
t159 = pkin(7) - qJ(3);
t150 = sin(t159);
t80 = t141 - t150 / 0.2e1;
t143 = cos(t158) / 0.2e1;
t152 = cos(t159);
t82 = t143 - t152 / 0.2e1;
t27 = -t71 * t106 + t124 * t80 - t82 * t154;
t101 = cos(pkin(7));
t99 = sin(pkin(7));
t55 = -t101 * t154 + t124 * t99;
t97 = pkin(13) + qJ(5);
t94 = sin(t97);
t95 = cos(t97);
t4 = -t27 * t95 + t55 * t94;
t186 = t27 * t94 + t55 * t95;
t103 = sin(qJ(6));
t105 = cos(qJ(6));
t180 = m(7) * pkin(5) + mrSges(7,1) * t105 - mrSges(7,2) * t103 + mrSges(6,1);
t149 = -m(7) * pkin(12) + mrSges(6,2) - mrSges(7,3);
t176 = m(6) + m(7);
t183 = m(5) + t176;
t121 = t108 * t174 + t175 * t130;
t145 = t162 * t175;
t182 = t101 * t145 + t121 * t99;
t157 = m(4) + t183;
t181 = pkin(2) * t157 + mrSges(3,1);
t100 = cos(pkin(13));
t98 = sin(pkin(13));
t132 = -m(5) * pkin(3) - t100 * mrSges(5,1) + t98 * mrSges(5,2) - mrSges(4,1);
t138 = -m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t179 = -mrSges(7,1) * t103 - mrSges(7,2) * t105 + t138;
t104 = sin(qJ(3));
t127 = t141 + t150 / 0.2e1;
t125 = t127 * t162;
t129 = t152 / 0.2e1 + t143;
t23 = t71 * t104 + t108 * t125 + t124 * t129;
t177 = -t149 * t94 + t180 * t95 - t132;
t173 = t55 * t98;
t172 = t94 * t99;
t171 = t95 * t99;
t164 = t108 * pkin(1) + pkin(9) * t145;
t163 = cos(pkin(6));
t147 = -t175 * pkin(1) + pkin(9) * t154;
t73 = t108 * t107 - t175 * t81;
t128 = t142 + t151 / 0.2e1;
t119 = -mrSges(3,2) + (t100 * mrSges(5,2) + mrSges(4,3) + (t176 * pkin(4) + mrSges(5,1)) * t98 + t157 * pkin(10)) * t99;
t83 = t144 - t153 / 0.2e1;
t118 = -t83 * t106 + t128 * t80 - t163 * t82;
t116 = -t71 * pkin(2) - t55 * pkin(10) + t147;
t114 = t73 * pkin(2) + t182 * pkin(10) + t164;
t102 = -pkin(11) - qJ(4);
t109 = t73 * t106 - t121 * t80 - t82 * t145;
t110 = t182 * t98;
t28 = t104 * t73 + t121 * t129 - t175 * t125;
t92 = pkin(4) * t100 + pkin(3);
t111 = pkin(4) * t110 - t28 * t102 + t109 * t92 + t114;
t70 = t163 * t101 - t128 * t99;
t51 = t128 * t106 + t83 * t80;
t50 = t128 * t104 - t83 * t129;
t43 = -t104 * t83 - t163 * t127 - t128 * t129;
t41 = -t121 * t106 - t73 * t80;
t40 = -t121 * t104 + t129 * t73;
t39 = -t124 * t106 - t71 * t80;
t38 = -t124 * t104 + t129 * t71;
t14 = t118 * t95 + t70 * t94;
t8 = t109 * t95 + t182 * t94;
t7 = t109 * t94 - t182 * t95;
t2 = t103 * t28 + t105 * t8;
t1 = -t103 * t8 + t105 * t28;
t3 = [(-t108 * mrSges(2,1) + t175 * mrSges(2,2) - m(3) * t164 - t73 * mrSges(3,1) + t121 * mrSges(3,2) - mrSges(3,3) * t145 - m(4) * t114 - t109 * mrSges(4,1) - t182 * mrSges(4,3) - m(5) * (pkin(3) * t109 + t114) - (t100 * t109 + t110) * mrSges(5,1) - (t100 * t182 - t109 * t98) * mrSges(5,2) - m(6) * t111 - t8 * mrSges(6,1) - m(7) * (t8 * pkin(5) + t111) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t149 * t7 + t138 * t28) * g(2) + (t175 * mrSges(2,1) + t108 * mrSges(2,2) - m(3) * t147 + t71 * mrSges(3,1) - t124 * mrSges(3,2) - mrSges(3,3) * t154 - m(4) * t116 - t27 * mrSges(4,1) + t55 * mrSges(4,3) - m(5) * (t27 * pkin(3) + t116) - (t100 * t27 - t173) * mrSges(5,1) - (-t100 * t55 - t27 * t98) * mrSges(5,2) + t149 * t186 + t180 * t4 - t179 * t23 - t176 * (-pkin(4) * t173 + t102 * t23 + t27 * t92 + t116)) * g(1) (t149 * (t83 * t171 + t51 * t94) + t132 * t51 - t180 * (-t83 * t172 + t51 * t95) + t179 * t50 + t119 * t83 - t176 * (-t50 * t102 + t51 * t92) - t181 * t128) * g(3) + (t149 * (-t171 * t71 + t39 * t94) + t132 * t39 - t180 * (t172 * t71 + t39 * t95) + t179 * t38 - t119 * t71 - t176 * (-t38 * t102 + t39 * t92) + t181 * t124) * g(2) + (t149 * (-t171 * t73 + t41 * t94) + t132 * t41 - t180 * (t172 * t73 + t41 * t95) + t179 * t40 - t119 * t73 - t176 * (-t40 * t102 + t41 * t92) + t181 * t121) * g(1) (-t176 * (-t102 * t118 - t43 * t92) + t179 * t118 + t177 * t43) * g(3) + (-t176 * (t102 * t27 - t23 * t92) - t179 * t27 + t177 * t23) * g(2) + (-t176 * (-t102 * t109 - t28 * t92) + t179 * t109 + t177 * t28) * g(1), t183 * (-g(1) * t28 - g(2) * t23 - g(3) * t43) (t149 * t14 - t180 * (-t118 * t94 + t70 * t95)) * g(3) + (t149 * t4 - t180 * t186) * g(2) + (t149 * t8 + t180 * t7) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t103 * t4 + t105 * t23) * mrSges(7,1) + (-t103 * t23 - t105 * t4) * mrSges(7,2)) - g(3) * ((-t103 * t14 + t105 * t43) * mrSges(7,1) + (-t103 * t43 - t105 * t14) * mrSges(7,2))];
taug  = t3(:);
