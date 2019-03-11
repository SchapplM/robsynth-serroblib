% Calculate potential energy for
% S6PRPRRP6
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
% Datum: 2019-03-08 20:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRP6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP6_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP6_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:18:07
% EndTime: 2019-03-08 20:18:07
% DurationCPUTime: 0.38s
% Computational Cost: add. (275->112), mult. (582->146), div. (0->0), fcn. (686->10), ass. (0->54)
t144 = rSges(7,1) + pkin(5);
t143 = rSges(7,2) + pkin(9);
t142 = rSges(6,3) + pkin(9);
t109 = sin(pkin(10));
t141 = g(1) * t109;
t111 = cos(pkin(10));
t140 = g(2) * t111;
t139 = rSges(7,3) + qJ(6);
t118 = cos(qJ(2));
t138 = qJ(3) * t118;
t110 = sin(pkin(6));
t137 = t109 * t110;
t114 = sin(qJ(4));
t136 = t110 * t114;
t115 = sin(qJ(2));
t135 = t110 * t115;
t117 = cos(qJ(4));
t134 = t110 * t117;
t133 = t110 * t118;
t112 = cos(pkin(6));
t132 = t112 * t115;
t131 = t112 * t118;
t130 = t112 * pkin(7) + qJ(1);
t129 = t111 * pkin(1) + pkin(7) * t137;
t128 = (-pkin(3) - pkin(7)) * t111;
t127 = pkin(2) * t135 + t130;
t105 = t109 * pkin(1);
t92 = t109 * t115 - t111 * t131;
t93 = t109 * t118 + t111 * t132;
t126 = t93 * pkin(2) + qJ(3) * t92 + t105;
t125 = t112 * pkin(3) + pkin(8) * t135 + t127;
t94 = t109 * t131 + t111 * t115;
t95 = -t109 * t132 + t111 * t118;
t124 = t95 * pkin(2) + qJ(3) * t94 + t129;
t123 = pkin(3) * t137 + t124;
t122 = pkin(8) * t93 + t126;
t81 = t109 * t134 + t114 * t94;
t121 = t81 * pkin(4) + pkin(8) * t95 + t123;
t97 = t112 * t117 - t114 * t133;
t120 = t97 * pkin(4) - qJ(3) * t133 + t125;
t83 = -t111 * t134 + t114 * t92;
t119 = t83 * pkin(4) + t110 * t128 + t122;
t116 = cos(qJ(5));
t113 = sin(qJ(5));
t96 = t112 * t114 + t117 * t133;
t85 = t113 * t135 + t116 * t97;
t84 = t113 * t97 - t116 * t135;
t82 = t111 * t136 + t117 * t92;
t80 = t109 * t136 - t94 * t117;
t77 = t113 * t93 + t116 * t83;
t76 = t113 * t83 - t93 * t116;
t75 = t113 * t95 + t116 * t81;
t74 = t113 * t81 - t95 * t116;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t111 - rSges(2,2) * t109) + g(2) * (rSges(2,1) * t109 + rSges(2,2) * t111) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t95 - rSges(3,2) * t94 + t129) + g(2) * (rSges(3,1) * t93 - rSges(3,2) * t92 + t105) + g(3) * (t112 * rSges(3,3) + t130) + (rSges(3,3) * t141 + g(3) * (rSges(3,1) * t115 + rSges(3,2) * t118) + (-rSges(3,3) - pkin(7)) * t140) * t110) - m(4) * (g(1) * (-rSges(4,2) * t95 + rSges(4,3) * t94 + t124) + g(2) * (-rSges(4,2) * t93 + rSges(4,3) * t92 + t126) + g(3) * (t112 * rSges(4,1) + t127) + (rSges(4,1) * t141 + g(3) * (-rSges(4,2) * t115 - rSges(4,3) * t118 - t138) + (-rSges(4,1) - pkin(7)) * t140) * t110) - m(5) * (g(1) * (rSges(5,1) * t81 - rSges(5,2) * t80 + (rSges(5,3) + pkin(8)) * t95 + t123) + g(2) * (rSges(5,1) * t83 + rSges(5,2) * t82 + rSges(5,3) * t93 + t122) + g(3) * (t97 * rSges(5,1) - t96 * rSges(5,2) + t125) + (g(3) * (rSges(5,3) * t115 - t138) + g(2) * t128) * t110) - m(6) * (g(1) * (rSges(6,1) * t75 - rSges(6,2) * t74 + t142 * t80 + t121) + g(2) * (rSges(6,1) * t77 - rSges(6,2) * t76 - t142 * t82 + t119) + g(3) * (t85 * rSges(6,1) - t84 * rSges(6,2) + t142 * t96 + t120)) - m(7) * (g(1) * (t139 * t74 + t143 * t80 + t144 * t75 + t121) + g(2) * (t139 * t76 - t143 * t82 + t144 * t77 + t119) + g(3) * (t139 * t84 + t143 * t96 + t144 * t85 + t120));
U  = t1;
