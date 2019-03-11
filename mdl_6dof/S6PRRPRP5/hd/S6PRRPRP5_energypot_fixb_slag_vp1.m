% Calculate potential energy for
% S6PRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2019-03-08 21:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRP5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP5_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRP5_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:45:22
% EndTime: 2019-03-08 21:45:22
% DurationCPUTime: 0.36s
% Computational Cost: add. (301->108), mult. (650->137), div. (0->0), fcn. (780->10), ass. (0->54)
t142 = pkin(4) + pkin(8);
t141 = rSges(5,1) + pkin(8);
t140 = rSges(7,1) + pkin(5);
t139 = rSges(7,2) + pkin(9);
t138 = rSges(4,3) + pkin(8);
t137 = rSges(6,3) + pkin(9);
t136 = cos(qJ(3));
t108 = sin(pkin(6));
t135 = pkin(7) * t108;
t134 = rSges(5,3) + qJ(4);
t133 = rSges(7,3) + qJ(6);
t132 = cos(pkin(6));
t111 = sin(qJ(3));
t131 = t108 * t111;
t112 = sin(qJ(2));
t130 = t108 * t112;
t114 = cos(qJ(2));
t129 = t108 * t114;
t128 = t132 * pkin(7) + qJ(1);
t107 = sin(pkin(10));
t109 = cos(pkin(10));
t127 = t109 * pkin(1) + t107 * t135;
t123 = t112 * t132;
t95 = -t107 * t123 + t109 * t114;
t126 = t95 * pkin(2) + t127;
t125 = pkin(2) * t130 + t128;
t124 = t108 * t136;
t122 = t114 * t132;
t84 = t107 * t131 + t95 * t136;
t121 = t84 * pkin(3) + t126;
t97 = t132 * t111 + t112 * t124;
t120 = t97 * pkin(3) + t125;
t104 = t107 * pkin(1);
t93 = t107 * t114 + t109 * t123;
t119 = t93 * pkin(2) - t109 * t135 + t104;
t82 = -t109 * t131 + t93 * t136;
t118 = t82 * pkin(3) + t119;
t83 = -t107 * t124 + t95 * t111;
t94 = t107 * t122 + t109 * t112;
t117 = t83 * qJ(4) + t142 * t94 + t121;
t96 = t111 * t130 - t132 * t136;
t116 = t96 * qJ(4) - t142 * t129 + t120;
t81 = t109 * t124 + t93 * t111;
t92 = t107 * t112 - t109 * t122;
t115 = t81 * qJ(4) + t142 * t92 + t118;
t113 = cos(qJ(5));
t110 = sin(qJ(5));
t86 = t96 * t110 - t113 * t129;
t85 = t110 * t129 + t96 * t113;
t76 = t83 * t110 + t94 * t113;
t75 = t94 * t110 - t83 * t113;
t74 = t81 * t110 + t92 * t113;
t73 = t92 * t110 - t81 * t113;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t109 * rSges(2,1) - t107 * rSges(2,2)) + g(2) * (t107 * rSges(2,1) + t109 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t95 * rSges(3,1) - t94 * rSges(3,2) + t127) + g(2) * (t93 * rSges(3,1) - t92 * rSges(3,2) + t104) + g(3) * (t132 * rSges(3,3) + t128) + (g(1) * rSges(3,3) * t107 + g(3) * (rSges(3,1) * t112 + rSges(3,2) * t114) + g(2) * (-rSges(3,3) - pkin(7)) * t109) * t108) - m(4) * (g(1) * (t84 * rSges(4,1) - t83 * rSges(4,2) + t138 * t94 + t126) + g(2) * (t82 * rSges(4,1) - t81 * rSges(4,2) + t138 * t92 + t119) + g(3) * (t97 * rSges(4,1) - t96 * rSges(4,2) - t138 * t129 + t125)) - m(5) * (g(1) * (-t84 * rSges(5,2) + t134 * t83 + t141 * t94 + t121) + g(2) * (-t82 * rSges(5,2) + t134 * t81 + t141 * t92 + t118) + g(3) * (-t97 * rSges(5,2) - t141 * t129 + t134 * t96 + t120)) - m(6) * (g(1) * (t76 * rSges(6,1) - t75 * rSges(6,2) + t137 * t84 + t117) + g(2) * (t74 * rSges(6,1) - t73 * rSges(6,2) + t137 * t82 + t115) + g(3) * (t86 * rSges(6,1) + t85 * rSges(6,2) + t137 * t97 + t116)) - m(7) * (g(1) * (t133 * t75 + t139 * t84 + t140 * t76 + t117) + g(2) * (t133 * t73 + t139 * t82 + t140 * t74 + t115) + g(3) * (-t133 * t85 + t139 * t97 + t140 * t86 + t116));
U  = t1;
