% Calculate potential energy for
% S6PRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-03-08 22:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRPP2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP2_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPP2_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:48:39
% EndTime: 2019-03-08 22:48:39
% DurationCPUTime: 0.43s
% Computational Cost: add. (335->108), mult. (739->138), div. (0->0), fcn. (903->10), ass. (0->52)
t140 = rSges(7,1) + pkin(5);
t139 = rSges(7,2) + qJ(5);
t138 = rSges(6,2) + pkin(9);
t137 = rSges(4,3) + pkin(8);
t136 = rSges(5,3) + pkin(9);
t135 = cos(qJ(3));
t107 = sin(pkin(6));
t134 = pkin(7) * t107;
t133 = rSges(6,3) + qJ(5);
t132 = cos(pkin(6));
t110 = sin(qJ(3));
t131 = t107 * t110;
t111 = sin(qJ(2));
t130 = t107 * t111;
t113 = cos(qJ(2));
t129 = t107 * t113;
t128 = t132 * pkin(7) + qJ(1);
t106 = sin(pkin(10));
t108 = cos(pkin(10));
t127 = t108 * pkin(1) + t106 * t134;
t122 = t111 * t132;
t95 = -t106 * t122 + t108 * t113;
t125 = t95 * pkin(2) + t127;
t124 = pkin(2) * t130 + t128;
t123 = t107 * t135;
t121 = t113 * t132;
t103 = t106 * pkin(1);
t93 = t106 * t113 + t108 * t122;
t120 = t93 * pkin(2) - t108 * t134 + t103;
t84 = t106 * t131 + t95 * t135;
t94 = t106 * t121 + t108 * t111;
t119 = t84 * pkin(3) + t94 * pkin(8) + t125;
t109 = sin(qJ(4));
t112 = cos(qJ(4));
t77 = t109 * t94 + t112 * t84;
t118 = t77 * pkin(4) + t119;
t97 = t132 * t110 + t111 * t123;
t117 = t97 * pkin(3) - pkin(8) * t129 + t124;
t86 = -t109 * t129 + t97 * t112;
t116 = t86 * pkin(4) + t117;
t82 = -t108 * t131 + t93 * t135;
t92 = t106 * t111 - t108 * t121;
t115 = t82 * pkin(3) + t92 * pkin(8) + t120;
t75 = t109 * t92 + t112 * t82;
t114 = t75 * pkin(4) + t115;
t96 = t110 * t130 - t132 * t135;
t85 = t97 * t109 + t112 * t129;
t83 = -t106 * t123 + t110 * t95;
t81 = t108 * t123 + t93 * t110;
t76 = t109 * t84 - t94 * t112;
t74 = t109 * t82 - t92 * t112;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t108 - rSges(2,2) * t106) + g(2) * (rSges(2,1) * t106 + rSges(2,2) * t108) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t95 - rSges(3,2) * t94 + t127) + g(2) * (t93 * rSges(3,1) - t92 * rSges(3,2) + t103) + g(3) * (t132 * rSges(3,3) + t128) + (g(1) * rSges(3,3) * t106 + g(3) * (rSges(3,1) * t111 + rSges(3,2) * t113) + g(2) * (-rSges(3,3) - pkin(7)) * t108) * t107) - m(4) * (g(1) * (t84 * rSges(4,1) - t83 * rSges(4,2) + t137 * t94 + t125) + g(2) * (t82 * rSges(4,1) - t81 * rSges(4,2) + t137 * t92 + t120) + g(3) * (t97 * rSges(4,1) - t96 * rSges(4,2) - t137 * t129 + t124)) - m(5) * (g(1) * (t77 * rSges(5,1) - t76 * rSges(5,2) + t136 * t83 + t119) + g(2) * (t75 * rSges(5,1) - t74 * rSges(5,2) + t136 * t81 + t115) + g(3) * (t86 * rSges(5,1) - t85 * rSges(5,2) + t136 * t96 + t117)) - m(6) * (g(1) * (t77 * rSges(6,1) + t133 * t76 + t138 * t83 + t118) + g(2) * (t75 * rSges(6,1) + t133 * t74 + t138 * t81 + t114) + g(3) * (t86 * rSges(6,1) + t133 * t85 + t138 * t96 + t116)) - m(7) * (g(1) * (t139 * t76 + t140 * t77 + t118) + g(2) * (t139 * t74 + t140 * t75 + t114) + g(3) * (t139 * t85 + t140 * t86 + t116) + (g(1) * t83 + g(2) * t81 + g(3) * t96) * (-rSges(7,3) + pkin(9) - qJ(6)));
U  = t1;
