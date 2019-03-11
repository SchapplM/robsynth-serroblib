% Calculate potential energy for
% S6RRPRRP14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 13:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP14_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP14_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP14_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP14_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:03:39
% EndTime: 2019-03-09 13:03:40
% DurationCPUTime: 0.40s
% Computational Cost: add. (275->112), mult. (582->143), div. (0->0), fcn. (686->10), ass. (0->55)
t146 = rSges(7,1) + pkin(5);
t145 = rSges(7,2) + pkin(10);
t144 = rSges(6,3) + pkin(10);
t115 = sin(qJ(1));
t143 = g(1) * t115;
t119 = cos(qJ(1));
t142 = g(2) * t119;
t111 = cos(pkin(6));
t141 = t111 * pkin(8) + pkin(7);
t140 = rSges(7,3) + qJ(6);
t118 = cos(qJ(2));
t139 = qJ(3) * t118;
t110 = sin(pkin(6));
t114 = sin(qJ(2));
t138 = t110 * t114;
t137 = t110 * t115;
t136 = t110 * t118;
t135 = t110 * t119;
t134 = t114 * t115;
t133 = t114 * t119;
t132 = t115 * t118;
t131 = t118 * t119;
t130 = t119 * pkin(1) + pkin(8) * t137;
t129 = pkin(2) * t138 + t141;
t128 = (-pkin(3) - pkin(8)) * t119;
t108 = t115 * pkin(1);
t95 = -t111 * t131 + t134;
t96 = t111 * t133 + t132;
t127 = t96 * pkin(2) + t95 * qJ(3) + t108;
t126 = t111 * pkin(3) + pkin(9) * t138 + t129;
t97 = t111 * t132 + t133;
t98 = -t111 * t134 + t131;
t125 = t98 * pkin(2) + qJ(3) * t97 + t130;
t124 = pkin(3) * t137 + t125;
t123 = t96 * pkin(9) + t127;
t113 = sin(qJ(4));
t117 = cos(qJ(4));
t94 = t111 * t117 - t113 * t136;
t122 = t94 * pkin(4) - qJ(3) * t136 + t126;
t84 = t113 * t97 + t117 * t137;
t121 = t84 * pkin(4) + pkin(9) * t98 + t124;
t86 = t95 * t113 - t117 * t135;
t120 = t86 * pkin(4) + t110 * t128 + t123;
t116 = cos(qJ(5));
t112 = sin(qJ(5));
t93 = t111 * t113 + t117 * t136;
t85 = t113 * t135 + t95 * t117;
t83 = t113 * t137 - t97 * t117;
t82 = t112 * t138 + t116 * t94;
t81 = t112 * t94 - t116 * t138;
t78 = t112 * t96 + t116 * t86;
t77 = t112 * t86 - t96 * t116;
t76 = t112 * t98 + t116 * t84;
t75 = t112 * t84 - t98 * t116;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t119 - t115 * rSges(2,2)) + g(2) * (t115 * rSges(2,1) + rSges(2,2) * t119) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t98 - rSges(3,2) * t97 + t130) + g(2) * (t96 * rSges(3,1) - t95 * rSges(3,2) + t108) + g(3) * (rSges(3,3) * t111 + t141) + (rSges(3,3) * t143 + g(3) * (rSges(3,1) * t114 + rSges(3,2) * t118) + (-rSges(3,3) - pkin(8)) * t142) * t110) - m(4) * (g(1) * (-rSges(4,2) * t98 + rSges(4,3) * t97 + t125) + g(2) * (-t96 * rSges(4,2) + t95 * rSges(4,3) + t127) + g(3) * (rSges(4,1) * t111 + t129) + (rSges(4,1) * t143 + g(3) * (-rSges(4,2) * t114 - rSges(4,3) * t118 - t139) + (-rSges(4,1) - pkin(8)) * t142) * t110) - m(5) * (g(1) * (rSges(5,1) * t84 - rSges(5,2) * t83 + (rSges(5,3) + pkin(9)) * t98 + t124) + g(2) * (t86 * rSges(5,1) + t85 * rSges(5,2) + t96 * rSges(5,3) + t123) + g(3) * (rSges(5,1) * t94 - rSges(5,2) * t93 + t126) + (g(3) * (rSges(5,3) * t114 - t139) + g(2) * t128) * t110) - m(6) * (g(1) * (rSges(6,1) * t76 - rSges(6,2) * t75 + t144 * t83 + t121) + g(2) * (t78 * rSges(6,1) - t77 * rSges(6,2) - t144 * t85 + t120) + g(3) * (rSges(6,1) * t82 - t81 * rSges(6,2) + t144 * t93 + t122)) - m(7) * (g(1) * (t140 * t75 + t145 * t83 + t146 * t76 + t121) + g(2) * (t140 * t77 - t145 * t85 + t146 * t78 + t120) + g(3) * (t140 * t81 + t145 * t93 + t146 * t82 + t122));
U  = t1;
