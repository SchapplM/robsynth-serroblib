% Calculate potential energy for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRP5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP5_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP5_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:13:57
% EndTime: 2019-03-08 20:13:57
% DurationCPUTime: 0.51s
% Computational Cost: add. (260->120), mult. (536->157), div. (0->0), fcn. (622->10), ass. (0->53)
t132 = rSges(7,3) + qJ(6) + pkin(9);
t98 = sin(pkin(10));
t131 = g(1) * t98;
t130 = rSges(6,3) + pkin(9);
t100 = cos(pkin(10));
t129 = g(2) * t100;
t99 = sin(pkin(6));
t128 = t98 * t99;
t127 = t100 * pkin(1) + pkin(7) * t128;
t103 = sin(qJ(5));
t108 = cos(qJ(2));
t101 = cos(pkin(6));
t105 = sin(qJ(2));
t119 = t101 * t105;
t82 = t100 * t119 + t108 * t98;
t126 = t103 * t82;
t104 = sin(qJ(4));
t125 = t104 * t99;
t124 = t105 * t99;
t107 = cos(qJ(4));
t123 = t107 * t99;
t122 = t108 * t99;
t121 = t101 * pkin(7) + qJ(1);
t120 = qJ(3) * t108;
t118 = t101 * t108;
t117 = t103 * t105;
t116 = pkin(2) * t124 + t121;
t115 = t100 * (-pkin(3) - pkin(7));
t114 = g(2) * t115;
t81 = -t100 * t118 + t105 * t98;
t94 = t98 * pkin(1);
t113 = t82 * pkin(2) + qJ(3) * t81 + t94;
t112 = t101 * pkin(3) + pkin(8) * t124 + t116;
t83 = t100 * t105 + t98 * t118;
t84 = t100 * t108 - t98 * t119;
t111 = t84 * pkin(2) + t83 * qJ(3) + t127;
t110 = pkin(3) * t128 + t111;
t109 = pkin(8) * t82 + t113;
t106 = cos(qJ(5));
t93 = pkin(5) * t106 + pkin(4);
t86 = t101 * t107 - t104 * t122;
t85 = t101 * t104 + t107 * t122;
t77 = t106 * t86 + t99 * t117;
t76 = -t103 * t86 + t106 * t124;
t75 = -t100 * t123 + t104 * t81;
t74 = t100 * t125 + t107 * t81;
t73 = t104 * t83 + t98 * t123;
t72 = -t83 * t107 + t98 * t125;
t71 = t106 * t75 + t126;
t70 = -t103 * t75 + t106 * t82;
t69 = t103 * t84 + t106 * t73;
t68 = -t103 * t73 + t106 * t84;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t100 - rSges(2,2) * t98) + g(2) * (rSges(2,1) * t98 + rSges(2,2) * t100) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t84 - rSges(3,2) * t83 + t127) + g(2) * (rSges(3,1) * t82 - rSges(3,2) * t81 + t94) + g(3) * (t101 * rSges(3,3) + t121) + (rSges(3,3) * t131 + g(3) * (rSges(3,1) * t105 + rSges(3,2) * t108) + (-rSges(3,3) - pkin(7)) * t129) * t99) - m(4) * (g(1) * (-rSges(4,2) * t84 + rSges(4,3) * t83 + t111) + g(2) * (-rSges(4,2) * t82 + rSges(4,3) * t81 + t113) + g(3) * (t101 * rSges(4,1) + t116) + (rSges(4,1) * t131 + g(3) * (-rSges(4,2) * t105 - rSges(4,3) * t108 - t120) + (-rSges(4,1) - pkin(7)) * t129) * t99) - m(5) * (g(1) * (rSges(5,1) * t73 - rSges(5,2) * t72 + (rSges(5,3) + pkin(8)) * t84 + t110) + g(2) * (rSges(5,1) * t75 + rSges(5,2) * t74 + rSges(5,3) * t82 + t109) + g(3) * (t86 * rSges(5,1) - t85 * rSges(5,2) + t112) + (g(3) * (rSges(5,3) * t105 - t120) + t114) * t99) - m(6) * (g(1) * (rSges(6,1) * t69 + rSges(6,2) * t68 + pkin(4) * t73 + pkin(8) * t84 + t130 * t72 + t110) + g(2) * (rSges(6,1) * t71 + rSges(6,2) * t70 + pkin(4) * t75 + t99 * t115 - t130 * t74 + t109) + g(3) * (t77 * rSges(6,1) + t76 * rSges(6,2) + t86 * pkin(4) - t99 * t120 + t130 * t85 + t112)) - m(7) * (g(1) * (rSges(7,1) * t69 + rSges(7,2) * t68 + t73 * t93 + (pkin(5) * t103 + pkin(8)) * t84 + t132 * t72 + t110) + g(2) * (rSges(7,1) * t71 + rSges(7,2) * t70 + pkin(5) * t126 - t132 * t74 + t75 * t93 + t109) + g(3) * (t77 * rSges(7,1) + t76 * rSges(7,2) + t132 * t85 + t86 * t93 + t112) + (g(3) * (pkin(5) * t117 - t120) + t114) * t99);
U  = t1;
