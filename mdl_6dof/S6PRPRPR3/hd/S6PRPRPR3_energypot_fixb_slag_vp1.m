% Calculate potential energy for
% S6PRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRPR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:34:20
% EndTime: 2019-03-08 19:34:21
% DurationCPUTime: 0.50s
% Computational Cost: add. (372->109), mult. (818->144), div. (0->0), fcn. (1011->12), ass. (0->52)
t106 = sin(pkin(11));
t113 = sin(qJ(2));
t116 = cos(qJ(2));
t134 = cos(pkin(11));
t95 = -t113 * t106 + t116 * t134;
t141 = rSges(6,1) + pkin(8);
t139 = rSges(5,3) + pkin(8);
t138 = pkin(9) + rSges(7,3);
t103 = pkin(2) * t116 + pkin(1);
t107 = sin(pkin(10));
t109 = cos(pkin(10));
t108 = sin(pkin(6));
t110 = cos(pkin(6));
t131 = t110 * t113;
t93 = pkin(2) * t131 + (-pkin(7) - qJ(3)) * t108;
t137 = t107 * t103 + t109 * t93;
t136 = rSges(4,3) * t108;
t135 = rSges(6,3) + qJ(5);
t112 = sin(qJ(4));
t133 = t108 * t112;
t115 = cos(qJ(4));
t132 = t108 * t115;
t130 = t110 * t116;
t128 = t110 * pkin(7) + qJ(1);
t94 = -t116 * t106 - t113 * t134;
t92 = t94 * t110;
t81 = t107 * t95 - t109 * t92;
t127 = t81 * pkin(3) + t137;
t75 = -t109 * t133 + t115 * t81;
t125 = t75 * pkin(4) + t127;
t124 = t108 * t113 * pkin(2) + t110 * qJ(3) + t128;
t83 = t107 * t92 + t109 * t95;
t97 = t109 * t103;
t123 = t83 * pkin(3) - t107 * t93 + t97;
t91 = t94 * t108;
t122 = -t91 * pkin(3) + t124;
t77 = t107 * t133 + t115 * t83;
t121 = t77 * pkin(4) + t123;
t86 = t110 * t112 - t115 * t91;
t120 = t86 * pkin(4) + t122;
t111 = sin(qJ(6));
t114 = cos(qJ(6));
t119 = t111 * rSges(7,1) + t114 * rSges(7,2) + qJ(5);
t118 = t114 * rSges(7,1) - t111 * rSges(7,2) + pkin(5) + pkin(8);
t117 = t110 * t95;
t90 = t95 * t108;
t85 = -t110 * t115 - t112 * t91;
t82 = -t107 * t117 + t109 * t94;
t80 = t107 * t94 + t109 * t117;
t76 = -t107 * t132 + t112 * t83;
t74 = t109 * t132 + t112 * t81;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t109 - rSges(2,2) * t107) + g(2) * (rSges(2,1) * t107 + rSges(2,2) * t109) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t109 * pkin(1) + (-t107 * t131 + t109 * t116) * rSges(3,1) + (-t107 * t130 - t109 * t113) * rSges(3,2)) + g(2) * (t107 * pkin(1) + (t107 * t116 + t109 * t131) * rSges(3,1) + (-t107 * t113 + t109 * t130) * rSges(3,2)) + g(3) * (t110 * rSges(3,3) + t128) + (g(3) * (rSges(3,1) * t113 + rSges(3,2) * t116) + (g(1) * t107 - g(2) * t109) * (rSges(3,3) + pkin(7))) * t108) - m(4) * (g(1) * (rSges(4,1) * t83 + rSges(4,2) * t82 + t97 + (-t93 + t136) * t107) + g(2) * (rSges(4,1) * t81 + rSges(4,2) * t80 - t109 * t136 + t137) + g(3) * (-rSges(4,1) * t91 + rSges(4,2) * t90 + rSges(4,3) * t110 + t124)) - m(5) * (g(1) * (rSges(5,1) * t77 - rSges(5,2) * t76 - t139 * t82 + t123) + g(2) * (rSges(5,1) * t75 - rSges(5,2) * t74 - t139 * t80 + t127) + g(3) * (rSges(5,1) * t86 - rSges(5,2) * t85 - t139 * t90 + t122)) - m(6) * (g(1) * (-rSges(6,2) * t77 + t135 * t76 - t141 * t82 + t121) + g(2) * (-rSges(6,2) * t75 + t135 * t74 - t141 * t80 + t125) + g(3) * (-t86 * rSges(6,2) + t135 * t85 - t141 * t90 + t120)) - m(7) * (g(1) * (-t118 * t82 + t119 * t76 + t138 * t77 + t121) + g(2) * (-t118 * t80 + t119 * t74 + t138 * t75 + t125) + g(3) * (-t118 * t90 + t119 * t85 + t138 * t86 + t120));
U  = t1;
