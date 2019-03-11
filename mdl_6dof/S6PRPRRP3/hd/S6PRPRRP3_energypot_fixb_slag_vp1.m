% Calculate potential energy for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRP3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:04:36
% EndTime: 2019-03-08 20:04:36
% DurationCPUTime: 0.38s
% Computational Cost: add. (343->123), mult. (545->160), div. (0->0), fcn. (633->12), ass. (0->53)
t136 = rSges(6,3) + pkin(9);
t106 = sin(pkin(11));
t135 = pkin(3) * t106;
t114 = sin(qJ(5));
t107 = sin(pkin(10));
t110 = cos(pkin(10));
t115 = sin(qJ(2));
t111 = cos(pkin(6));
t117 = cos(qJ(2));
t124 = t111 * t117;
t86 = t107 * t115 - t110 * t124;
t134 = t86 * t114;
t88 = t107 * t124 + t110 * t115;
t133 = t88 * t114;
t132 = rSges(7,3) + qJ(6) + pkin(9);
t131 = qJ(3) + rSges(4,3);
t108 = sin(pkin(6));
t129 = t107 * t108;
t130 = t110 * pkin(1) + pkin(7) * t129;
t128 = t108 * t110;
t127 = t108 * t115;
t126 = t108 * t117;
t125 = t111 * t115;
t123 = t111 * pkin(7) + qJ(1);
t122 = t106 * t129;
t121 = t114 * t126;
t113 = -pkin(8) - qJ(3);
t89 = -t107 * t125 + t110 * t117;
t109 = cos(pkin(11));
t98 = t109 * pkin(3) + pkin(2);
t120 = pkin(3) * t122 - t88 * t113 + t89 * t98 + t130;
t119 = t111 * t135 + t113 * t126 + t98 * t127 + t123;
t102 = t107 * pkin(1);
t87 = t107 * t117 + t110 * t125;
t118 = t102 + t87 * t98 + (-pkin(7) - t135) * t128 - t86 * t113;
t116 = cos(qJ(5));
t105 = pkin(11) + qJ(4);
t101 = cos(t105);
t100 = sin(t105);
t99 = t116 * pkin(5) + pkin(4);
t83 = t111 * t100 + t101 * t127;
t82 = t100 * t127 - t111 * t101;
t79 = t83 * t116 - t121;
t78 = -t83 * t114 - t116 * t126;
t77 = t100 * t129 + t89 * t101;
t76 = t89 * t100 - t101 * t129;
t75 = -t100 * t128 + t87 * t101;
t74 = t87 * t100 + t101 * t128;
t73 = t77 * t116 + t133;
t72 = -t77 * t114 + t88 * t116;
t71 = t75 * t116 + t134;
t70 = -t75 * t114 + t86 * t116;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t110 * rSges(2,1) - t107 * rSges(2,2)) + g(2) * (t107 * rSges(2,1) + t110 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t89 * rSges(3,1) - t88 * rSges(3,2) + t130) + g(2) * (t87 * rSges(3,1) - t86 * rSges(3,2) + t102) + g(3) * (t111 * rSges(3,3) + t123) + (g(1) * rSges(3,3) * t107 + g(3) * (rSges(3,1) * t115 + rSges(3,2) * t117) + g(2) * (-rSges(3,3) - pkin(7)) * t110) * t108) - m(4) * (g(1) * (t89 * pkin(2) + (t89 * t109 + t122) * rSges(4,1) + (-t89 * t106 + t109 * t129) * rSges(4,2) + t131 * t88 + t130) + g(2) * (t87 * pkin(2) + t102 - pkin(7) * t128 + (-t106 * t128 + t87 * t109) * rSges(4,1) + (-t87 * t106 - t109 * t128) * rSges(4,2) + t131 * t86) + g(3) * ((t106 * rSges(4,1) + t109 * rSges(4,2)) * t111 + (-t131 * t117 + (t109 * rSges(4,1) - t106 * rSges(4,2) + pkin(2)) * t115) * t108 + t123)) - m(5) * (g(1) * (t77 * rSges(5,1) - t76 * rSges(5,2) + t88 * rSges(5,3) + t120) + g(2) * (t75 * rSges(5,1) - t74 * rSges(5,2) + t86 * rSges(5,3) + t118) + g(3) * (t83 * rSges(5,1) - t82 * rSges(5,2) - rSges(5,3) * t126 + t119)) - m(6) * (g(1) * (t73 * rSges(6,1) + t72 * rSges(6,2) + t77 * pkin(4) + t136 * t76 + t120) + g(2) * (t71 * rSges(6,1) + t70 * rSges(6,2) + t75 * pkin(4) + t136 * t74 + t118) + g(3) * (t79 * rSges(6,1) + t78 * rSges(6,2) + t83 * pkin(4) + t136 * t82 + t119)) - m(7) * (g(1) * (t73 * rSges(7,1) + t72 * rSges(7,2) + pkin(5) * t133 + t132 * t76 + t77 * t99 + t120) + g(2) * (t71 * rSges(7,1) + t70 * rSges(7,2) + pkin(5) * t134 + t132 * t74 + t75 * t99 + t118) + g(3) * (t79 * rSges(7,1) + t78 * rSges(7,2) - pkin(5) * t121 + t132 * t82 + t83 * t99 + t119));
U  = t1;
