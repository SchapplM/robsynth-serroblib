% Calculate potential energy for
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPPRR2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR2_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPPRR2_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:17:35
% EndTime: 2019-03-08 19:17:36
% DurationCPUTime: 0.56s
% Computational Cost: add. (336->119), mult. (724->162), div. (0->0), fcn. (881->12), ass. (0->48)
t101 = sin(pkin(11));
t104 = cos(pkin(11));
t109 = sin(qJ(2));
t112 = cos(qJ(2));
t91 = -t109 * t101 + t112 * t104;
t133 = rSges(6,3) + pkin(8);
t132 = pkin(9) + rSges(7,3);
t102 = sin(pkin(10));
t105 = cos(pkin(10));
t106 = cos(pkin(6));
t118 = t112 * t101 + t109 * t104;
t117 = t118 * t106;
t79 = -t102 * t117 + t105 * t91;
t97 = t112 * pkin(2) + pkin(1);
t93 = t105 * t97;
t131 = t79 * pkin(3) + t93;
t103 = sin(pkin(6));
t125 = t106 * t109;
t89 = pkin(2) * t125 + (-pkin(7) - qJ(3)) * t103;
t130 = t102 * t97 + t105 * t89;
t129 = rSges(5,3) + qJ(4);
t128 = t106 * pkin(7) + qJ(1);
t127 = t102 * t103;
t126 = t105 * t103;
t124 = t106 * t112;
t77 = t102 * t91 + t105 * t117;
t121 = t77 * pkin(3) + t130;
t120 = t103 * t109 * pkin(2) + t106 * qJ(3) + t128;
t88 = t118 * t103;
t119 = t88 * pkin(3) + t120;
t116 = t91 * t106;
t78 = -t102 * t116 - t105 * t118;
t115 = pkin(4) * t127 - t78 * qJ(4) - t102 * t89 + t131;
t87 = t91 * t103;
t114 = t106 * pkin(4) - t87 * qJ(4) + t119;
t76 = -t102 * t118 + t105 * t116;
t113 = -pkin(4) * t126 - t76 * qJ(4) + t121;
t111 = cos(qJ(5));
t110 = cos(qJ(6));
t108 = sin(qJ(5));
t107 = sin(qJ(6));
t81 = t106 * t111 - t87 * t108;
t80 = t106 * t108 + t87 * t111;
t72 = -t76 * t108 - t111 * t126;
t71 = t108 * t126 - t76 * t111;
t70 = -t78 * t108 + t111 * t127;
t69 = t108 * t127 + t78 * t111;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t105 * rSges(2,1) - t102 * rSges(2,2)) + g(2) * (t102 * rSges(2,1) + t105 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t105 * pkin(1) + (-t102 * t125 + t105 * t112) * rSges(3,1) + (-t102 * t124 - t105 * t109) * rSges(3,2)) + g(2) * (t102 * pkin(1) + (t102 * t112 + t105 * t125) * rSges(3,1) + (-t102 * t109 + t105 * t124) * rSges(3,2)) + g(3) * (t106 * rSges(3,3) + t128) + (g(3) * (rSges(3,1) * t109 + rSges(3,2) * t112) + (g(1) * t102 - g(2) * t105) * (rSges(3,3) + pkin(7))) * t103) - m(4) * (g(1) * (t79 * rSges(4,1) + t78 * rSges(4,2) + t93 + (rSges(4,3) * t103 - t89) * t102) + g(2) * (t77 * rSges(4,1) + t76 * rSges(4,2) - rSges(4,3) * t126 + t130) + g(3) * (t88 * rSges(4,1) + t87 * rSges(4,2) + t106 * rSges(4,3) + t120)) - m(5) * (g(1) * (-t79 * rSges(5,2) - t129 * t78 + (rSges(5,1) * t103 - t89) * t102 + t131) + g(2) * (-rSges(5,1) * t126 - t77 * rSges(5,2) - t129 * t76 + t121) + g(3) * (t106 * rSges(5,1) - t88 * rSges(5,2) - t129 * t87 + t119)) - m(6) * (g(1) * (t70 * rSges(6,1) - t69 * rSges(6,2) + t133 * t79 + t115) + g(2) * (t72 * rSges(6,1) + t71 * rSges(6,2) + t133 * t77 + t113) + g(3) * (t81 * rSges(6,1) - t80 * rSges(6,2) + t133 * t88 + t114)) - m(7) * (g(1) * (t70 * pkin(5) + t79 * pkin(8) + (t79 * t107 + t70 * t110) * rSges(7,1) + (-t70 * t107 + t79 * t110) * rSges(7,2) + t132 * t69 + t115) + g(2) * (t72 * pkin(5) + t77 * pkin(8) + (t77 * t107 + t72 * t110) * rSges(7,1) + (-t72 * t107 + t77 * t110) * rSges(7,2) - t132 * t71 + t113) + g(3) * (t81 * pkin(5) + t88 * pkin(8) + (t88 * t107 + t81 * t110) * rSges(7,1) + (-t81 * t107 + t88 * t110) * rSges(7,2) + t132 * t80 + t114));
U  = t1;
