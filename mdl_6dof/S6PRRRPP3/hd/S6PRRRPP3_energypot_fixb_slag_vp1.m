% Calculate potential energy for
% S6PRRRPP3
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
% Datum: 2019-03-08 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRPP3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPP3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:54:21
% EndTime: 2019-03-08 22:54:22
% DurationCPUTime: 0.43s
% Computational Cost: add. (335->108), mult. (739->138), div. (0->0), fcn. (903->10), ass. (0->52)
t142 = rSges(7,2) + qJ(5);
t141 = rSges(7,3) + qJ(6);
t140 = rSges(6,1) + pkin(9);
t139 = rSges(4,3) + pkin(8);
t138 = rSges(5,3) + pkin(9);
t137 = cos(qJ(3));
t109 = sin(pkin(6));
t136 = pkin(7) * t109;
t135 = rSges(6,3) + qJ(5);
t134 = cos(pkin(6));
t112 = sin(qJ(3));
t133 = t109 * t112;
t113 = sin(qJ(2));
t132 = t109 * t113;
t115 = cos(qJ(2));
t131 = t109 * t115;
t130 = t134 * pkin(7) + qJ(1);
t108 = sin(pkin(10));
t110 = cos(pkin(10));
t129 = t110 * pkin(1) + t108 * t136;
t124 = t113 * t134;
t96 = -t108 * t124 + t110 * t115;
t127 = t96 * pkin(2) + t129;
t126 = pkin(2) * t132 + t130;
t125 = t109 * t137;
t123 = t115 * t134;
t105 = t108 * pkin(1);
t94 = t108 * t115 + t110 * t124;
t122 = t94 * pkin(2) - t110 * t136 + t105;
t85 = t108 * t133 + t96 * t137;
t95 = t108 * t123 + t110 * t113;
t121 = t85 * pkin(3) + t95 * pkin(8) + t127;
t111 = sin(qJ(4));
t114 = cos(qJ(4));
t78 = t111 * t95 + t114 * t85;
t120 = t78 * pkin(4) + t121;
t98 = t134 * t112 + t113 * t125;
t119 = t98 * pkin(3) - pkin(8) * t131 + t126;
t87 = -t111 * t131 + t114 * t98;
t118 = t87 * pkin(4) + t119;
t83 = -t110 * t133 + t94 * t137;
t93 = t108 * t113 - t110 * t123;
t117 = t83 * pkin(3) + t93 * pkin(8) + t122;
t76 = t111 * t93 + t114 * t83;
t116 = t76 * pkin(4) + t117;
t97 = t112 * t132 - t134 * t137;
t86 = t98 * t111 + t114 * t131;
t84 = -t108 * t125 + t112 * t96;
t82 = t110 * t125 + t94 * t112;
t77 = t111 * t85 - t95 * t114;
t75 = t111 * t83 - t93 * t114;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t110 - rSges(2,2) * t108) + g(2) * (rSges(2,1) * t108 + rSges(2,2) * t110) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t96 - rSges(3,2) * t95 + t129) + g(2) * (t94 * rSges(3,1) - t93 * rSges(3,2) + t105) + g(3) * (t134 * rSges(3,3) + t130) + (g(1) * rSges(3,3) * t108 + g(3) * (rSges(3,1) * t113 + rSges(3,2) * t115) + g(2) * (-rSges(3,3) - pkin(7)) * t110) * t109) - m(4) * (g(1) * (t85 * rSges(4,1) - t84 * rSges(4,2) + t139 * t95 + t127) + g(2) * (t83 * rSges(4,1) - t82 * rSges(4,2) + t139 * t93 + t122) + g(3) * (t98 * rSges(4,1) - t97 * rSges(4,2) - t139 * t131 + t126)) - m(5) * (g(1) * (t78 * rSges(5,1) - t77 * rSges(5,2) + t138 * t84 + t121) + g(2) * (t76 * rSges(5,1) - t75 * rSges(5,2) + t138 * t82 + t117) + g(3) * (t87 * rSges(5,1) - t86 * rSges(5,2) + t138 * t97 + t119)) - m(6) * (g(1) * (-t78 * rSges(6,2) + t135 * t77 + t140 * t84 + t120) + g(2) * (-t76 * rSges(6,2) + t135 * t75 + t140 * t82 + t116) + g(3) * (-t87 * rSges(6,2) + t135 * t86 + t140 * t97 + t118)) - m(7) * (g(1) * (t141 * t78 + t142 * t77 + t120) + g(2) * (t141 * t76 + t142 * t75 + t116) + g(3) * (t141 * t87 + t142 * t86 + t118) + (g(1) * t84 + g(2) * t82 + g(3) * t97) * (rSges(7,1) + pkin(5) + pkin(9)));
U  = t1;
