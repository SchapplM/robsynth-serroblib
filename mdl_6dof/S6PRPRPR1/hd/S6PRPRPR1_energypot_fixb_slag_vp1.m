% Calculate potential energy for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRPR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR1_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR1_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:24:33
% EndTime: 2019-03-08 19:24:34
% DurationCPUTime: 0.64s
% Computational Cost: add. (383->126), mult. (724->174), div. (0->0), fcn. (880->14), ass. (0->53)
t115 = sin(pkin(11));
t118 = cos(pkin(11));
t124 = sin(qJ(2));
t127 = cos(qJ(2));
t99 = -t124 * t115 + t127 * t118;
t147 = pkin(8) + rSges(5,3);
t146 = pkin(9) + rSges(7,3);
t109 = t127 * pkin(2) + pkin(1);
t116 = sin(pkin(10));
t119 = cos(pkin(10));
t117 = sin(pkin(6));
t120 = cos(pkin(6));
t141 = t120 * t124;
t97 = pkin(2) * t141 + (-pkin(7) - qJ(3)) * t117;
t145 = t116 * t109 + t119 * t97;
t144 = t116 * t117;
t143 = t119 * t117;
t123 = sin(qJ(4));
t142 = t120 * t123;
t140 = t120 * t127;
t137 = t120 * pkin(7) + qJ(1);
t136 = t123 * t144;
t135 = t123 * t143;
t103 = t119 * t109;
t134 = -t116 * t97 + t103;
t133 = t117 * t124 * pkin(2) + t120 * qJ(3) + t137;
t132 = t127 * t115 + t124 * t118;
t126 = cos(qJ(4));
t108 = t126 * pkin(4) + pkin(3);
t121 = -qJ(5) - pkin(8);
t130 = t99 * t120;
t86 = -t116 * t130 - t119 * t132;
t96 = t132 * t120;
t87 = -t116 * t96 + t119 * t99;
t131 = pkin(4) * t136 + t87 * t108 + t86 * t121 + t134;
t94 = t99 * t117;
t95 = t132 * t117;
t129 = pkin(4) * t142 + t95 * t108 + t94 * t121 + t133;
t84 = -t116 * t132 + t119 * t130;
t85 = t116 * t99 + t119 * t96;
t128 = -pkin(4) * t135 + t85 * t108 + t84 * t121 + t145;
t125 = cos(qJ(6));
t122 = sin(qJ(6));
t114 = qJ(4) + pkin(12);
t111 = cos(t114);
t110 = sin(t114);
t89 = t120 * t110 + t95 * t111;
t88 = t95 * t110 - t120 * t111;
t79 = t110 * t144 + t87 * t111;
t78 = t87 * t110 - t111 * t144;
t77 = -t110 * t143 + t85 * t111;
t76 = t85 * t110 + t111 * t143;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t119 * rSges(2,1) - t116 * rSges(2,2)) + g(2) * (t116 * rSges(2,1) + t119 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t119 * pkin(1) + (-t116 * t141 + t119 * t127) * rSges(3,1) + (-t116 * t140 - t119 * t124) * rSges(3,2)) + g(2) * (t116 * pkin(1) + (t116 * t127 + t119 * t141) * rSges(3,1) + (-t116 * t124 + t119 * t140) * rSges(3,2)) + g(3) * (t120 * rSges(3,3) + t137) + (g(3) * (rSges(3,1) * t124 + rSges(3,2) * t127) + (g(1) * t116 - g(2) * t119) * (rSges(3,3) + pkin(7))) * t117) - m(4) * (g(1) * (t87 * rSges(4,1) + t86 * rSges(4,2) + t103 + (rSges(4,3) * t117 - t97) * t116) + g(2) * (t85 * rSges(4,1) + t84 * rSges(4,2) - rSges(4,3) * t143 + t145) + g(3) * (t95 * rSges(4,1) + t94 * rSges(4,2) + t120 * rSges(4,3) + t133)) - m(5) * (g(1) * (t87 * pkin(3) + (t87 * t126 + t136) * rSges(5,1) + (-t87 * t123 + t126 * t144) * rSges(5,2) - t147 * t86 + t134) + g(2) * (t85 * pkin(3) + (t85 * t126 - t135) * rSges(5,1) + (-t85 * t123 - t126 * t143) * rSges(5,2) - t147 * t84 + t145) + g(3) * (t95 * pkin(3) + (t95 * t126 + t142) * rSges(5,1) + (t120 * t126 - t95 * t123) * rSges(5,2) - t147 * t94 + t133)) - m(6) * (g(1) * (t79 * rSges(6,1) - t78 * rSges(6,2) - t86 * rSges(6,3) + t131) + g(2) * (t77 * rSges(6,1) - t76 * rSges(6,2) - t84 * rSges(6,3) + t128) + g(3) * (t89 * rSges(6,1) - t88 * rSges(6,2) - t94 * rSges(6,3) + t129)) - m(7) * (g(1) * (t79 * pkin(5) + (-t86 * t122 + t79 * t125) * rSges(7,1) + (-t79 * t122 - t86 * t125) * rSges(7,2) + t146 * t78 + t131) + g(2) * (t77 * pkin(5) + (-t84 * t122 + t77 * t125) * rSges(7,1) + (-t77 * t122 - t84 * t125) * rSges(7,2) + t146 * t76 + t128) + g(3) * (t89 * pkin(5) + (-t94 * t122 + t89 * t125) * rSges(7,1) + (-t89 * t122 - t94 * t125) * rSges(7,2) + t146 * t88 + t129));
U  = t1;
