% Calculate Gravitation load on the joints for
% S6RRPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2018-11-23 17:29
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR12_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR12_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR12_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR12_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:29:10
% EndTime: 2018-11-23 17:29:11
% DurationCPUTime: 1.12s
% Computational Cost: add. (1349->146), mult. (1459->185), div. (0->0), fcn. (1403->16), ass. (0->85)
t143 = m(6) + m(7);
t144 = m(4) + m(5);
t152 = -t144 - t143;
t142 = t152 * qJ(3) + mrSges(3,2) - mrSges(4,3);
t77 = sin(qJ(4));
t81 = cos(qJ(4));
t153 = t77 * mrSges(5,1) + t81 * mrSges(5,2) - t142;
t76 = sin(qJ(6));
t80 = cos(qJ(6));
t145 = mrSges(7,1) * t80 - mrSges(7,2) * t76 + mrSges(6,1);
t149 = mrSges(6,2) - mrSges(7,3);
t106 = -m(7) * pkin(11) + t149;
t150 = -m(7) * pkin(5) - t145;
t74 = sin(pkin(6));
t83 = cos(qJ(1));
t130 = t74 * t83;
t78 = sin(qJ(2));
t79 = sin(qJ(1));
t123 = pkin(6) - qJ(2);
t103 = cos(t123) / 0.2e1;
t122 = pkin(6) + qJ(2);
t110 = cos(t122);
t86 = t103 + t110 / 0.2e1;
t45 = t78 * t79 - t83 * t86;
t148 = t77 * t130 + t45 * t81;
t131 = t74 * t79;
t48 = t83 * t78 + t79 * t86;
t19 = -t77 * t131 + t48 * t81;
t108 = sin(t122);
t101 = t108 / 0.2e1;
t109 = sin(t123);
t102 = t109 / 0.2e1;
t60 = t101 + t102;
t75 = cos(pkin(6));
t147 = -t60 * t81 - t75 * t77;
t94 = -m(5) * pkin(9) - mrSges(3,1) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t141 = mrSges(7,1) * t76 + mrSges(7,2) * t80 - t94;
t73 = qJ(4) + qJ(5);
t70 = sin(t73);
t71 = cos(t73);
t140 = -t106 * t71 + t150 * t70 - t153;
t139 = pkin(4) * t77;
t136 = t45 * t77;
t134 = t48 * t77;
t82 = cos(qJ(2));
t128 = t79 * t82;
t127 = t83 * t82;
t125 = t148 * pkin(4);
t124 = t83 * pkin(1) + pkin(8) * t131;
t95 = t101 - t109 / 0.2e1;
t49 = -t79 * t95 + t127;
t117 = t49 * pkin(2) + t124;
t114 = -pkin(1) * t79 + pkin(8) * t130;
t13 = t70 * t131 - t48 * t71;
t14 = t71 * t131 + t48 * t70;
t113 = -t13 * pkin(5) + pkin(11) * t14;
t15 = t70 * t130 + t45 * t71;
t17 = t71 * t130 - t45 * t70;
t112 = t15 * pkin(5) - pkin(11) * t17;
t29 = -t60 * t71 - t70 * t75;
t30 = -t60 * t70 + t71 * t75;
t111 = t29 * pkin(5) + pkin(11) * t30;
t107 = -m(5) * pkin(3) - mrSges(4,1) - mrSges(3,3);
t105 = t147 * pkin(4);
t46 = t83 * t95 + t128;
t104 = t46 * pkin(2) - t114;
t98 = t19 * pkin(4);
t69 = pkin(4) * t81 + pkin(3);
t84 = -pkin(10) - pkin(9);
t97 = pkin(4) * t134 + t69 * t131 - t49 * t84 + t117;
t96 = t102 - t108 / 0.2e1;
t93 = t145 * t13 + t149 * t14;
t92 = -t145 * t15 - t149 * t17;
t89 = -t145 * t29 + t149 * t30;
t63 = t81 * t130;
t61 = t103 - t110 / 0.2e1;
t57 = t60 * pkin(2);
t50 = t79 * t96 + t127;
t47 = -t83 * t96 + t128;
t43 = t48 * pkin(2);
t41 = t45 * pkin(2);
t20 = t81 * t131 + t134;
t2 = t14 * t80 + t49 * t76;
t1 = -t14 * t76 + t49 * t80;
t3 = [(-t83 * mrSges(2,1) - m(3) * t124 - t20 * mrSges(5,1) - t19 * mrSges(5,2) - m(6) * t97 - t14 * mrSges(6,1) - m(7) * (pkin(5) * t14 + t97) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t142 * t48 + t106 * t13 + (t107 * t74 + mrSges(2,2)) * t79 + t94 * t49 - t144 * t117) * g(2) + (t79 * mrSges(2,1) - m(3) * t114 - t63 * mrSges(5,1) + t106 * t15 + t153 * t45 + (mrSges(2,2) + (mrSges(5,2) * t77 + t107) * t74) * t83 + t150 * t17 + t141 * t46 + t143 * (pkin(4) * t136 - t69 * t130 - t46 * t84 + t104) + t144 * t104) * g(1) (-t144 * t57 - t143 * (t61 * t139 - t60 * t84 + t57) + t140 * t61 - t141 * t60) * g(3) + (t144 * t41 - t143 * (t47 * t139 + t45 * t84 - t41) + t140 * t47 + t141 * t45) * g(2) + (t144 * t43 - t143 * (t50 * t139 + t48 * t84 - t43) + t140 * t50 + t141 * t48) * g(1) -(-g(1) * t48 - g(2) * t45 + g(3) * t60) * t152 (-t147 * mrSges(5,1) - (t60 * t77 - t75 * t81) * mrSges(5,2) - m(6) * t105 - m(7) * (t105 + t111) + t89) * g(3) + (-t148 * mrSges(5,1) - (t63 - t136) * mrSges(5,2) - m(6) * t125 - m(7) * (t112 + t125) + t92) * g(2) + (-t19 * mrSges(5,1) + t20 * mrSges(5,2) - m(6) * t98 - m(7) * (t113 + t98) + t93) * g(1) (-m(7) * t111 + t89) * g(3) + (-m(7) * t112 + t92) * g(2) + (-m(7) * t113 + t93) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t17 * t76 + t46 * t80) * mrSges(7,1) + (t17 * t80 - t46 * t76) * mrSges(7,2)) - g(3) * ((-t30 * t76 + t61 * t80) * mrSges(7,1) + (-t30 * t80 - t61 * t76) * mrSges(7,2))];
taug  = t3(:);
