% Calculate potential energy for
% S6RRRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
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
% Datum: 2019-03-09 15:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPPP1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPP1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPP1_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPP1_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:14:23
% EndTime: 2019-03-09 15:14:23
% DurationCPUTime: 0.33s
% Computational Cost: add. (292->112), mult. (668->147), div. (0->0), fcn. (786->10), ass. (0->52)
t149 = rSges(7,1) + pkin(5);
t117 = sin(qJ(2));
t148 = t117 * pkin(2) + pkin(7);
t147 = rSges(7,2) + qJ(5);
t146 = rSges(6,3) + qJ(5);
t145 = rSges(7,3) + qJ(6);
t144 = cos(pkin(10));
t114 = sin(pkin(6));
t143 = qJ(4) * t114;
t120 = cos(qJ(2));
t142 = t114 * t120;
t116 = sin(qJ(3));
t141 = t116 * t117;
t118 = sin(qJ(1));
t140 = t117 * t118;
t121 = cos(qJ(1));
t139 = t117 * t121;
t138 = t118 * t120;
t137 = t120 * t121;
t136 = t121 * pkin(1) + t118 * pkin(8);
t135 = t114 * t141;
t115 = cos(pkin(6));
t134 = t115 * t140;
t133 = t115 * t139;
t132 = t115 * t144;
t131 = t117 * t144;
t130 = pkin(2) * t137 + pkin(9) * t139 + t136;
t129 = t114 * t131;
t111 = t118 * pkin(1);
t128 = pkin(2) * t138 - pkin(8) * t121 + pkin(9) * t140 + t111;
t119 = cos(qJ(3));
t127 = qJ(4) * t135 + t117 * t119 * pkin(3) + (-qJ(4) * t115 - pkin(9)) * t120 + t148;
t95 = -t116 * t137 + t118 * t119;
t96 = t118 * t116 + t119 * t137;
t126 = t96 * pkin(3) + qJ(4) * t133 - t95 * t143 + t130;
t113 = sin(pkin(10));
t85 = t119 * t131 + (-t115 * t141 - t142) * t113;
t125 = t85 * pkin(4) + t127;
t82 = t96 * t144 + (t114 * t139 + t115 * t95) * t113;
t124 = t82 * pkin(4) + t126;
t93 = -t116 * t138 - t119 * t121;
t94 = -t116 * t121 + t119 * t138;
t123 = t94 * pkin(3) + qJ(4) * t134 - t93 * t143 + t128;
t80 = t94 * t144 + (t114 * t140 + t115 * t93) * t113;
t122 = t80 * pkin(4) + t123;
t92 = -t115 * t120 + t135;
t87 = -t95 * t114 + t133;
t86 = -t114 * t93 + t134;
t84 = t144 * t142 + (t113 * t119 + t116 * t132) * t117;
t81 = t113 * t96 - t121 * t129 - t95 * t132;
t79 = t113 * t94 - t118 * t129 - t93 * t132;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t121 - t118 * rSges(2,2)) + g(2) * (t118 * rSges(2,1) + rSges(2,2) * t121) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (t118 * rSges(3,3) + t136) + g(2) * (rSges(3,1) * t138 - rSges(3,2) * t140 + t111) + g(3) * (rSges(3,1) * t117 + rSges(3,2) * t120 + pkin(7)) + (g(1) * (rSges(3,1) * t120 - rSges(3,2) * t117) + g(2) * (-rSges(3,3) - pkin(8))) * t121) - m(4) * (g(1) * (t96 * rSges(4,1) + t95 * rSges(4,2) + rSges(4,3) * t139 + t130) + g(2) * (t94 * rSges(4,1) + t93 * rSges(4,2) + rSges(4,3) * t140 + t128) + g(3) * ((-rSges(4,3) - pkin(9)) * t120 + (rSges(4,1) * t119 - rSges(4,2) * t116) * t117 + t148)) - m(5) * (g(1) * (rSges(5,1) * t82 - rSges(5,2) * t81 + rSges(5,3) * t87 + t126) + g(2) * (t80 * rSges(5,1) - t79 * rSges(5,2) + t86 * rSges(5,3) + t123) + g(3) * (rSges(5,1) * t85 - rSges(5,2) * t84 + rSges(5,3) * t92 + t127)) - m(6) * (g(1) * (rSges(6,1) * t87 - rSges(6,2) * t82 + t146 * t81 + t124) + g(2) * (t86 * rSges(6,1) - t80 * rSges(6,2) + t146 * t79 + t122) + g(3) * (rSges(6,1) * t92 - rSges(6,2) * t85 + t146 * t84 + t125)) - m(7) * (g(1) * (t145 * t82 + t147 * t81 + t149 * t87 + t124) + g(2) * (t145 * t80 + t147 * t79 + t149 * t86 + t122) + g(3) * (t145 * t85 + t147 * t84 + t149 * t92 + t125));
U  = t1;
