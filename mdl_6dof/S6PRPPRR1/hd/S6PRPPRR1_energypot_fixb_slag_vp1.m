% Calculate potential energy for
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
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
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPPRR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR1_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPPRR1_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:13:49
% EndTime: 2019-03-08 19:13:50
% DurationCPUTime: 0.59s
% Computational Cost: add. (383->126), mult. (724->174), div. (0->0), fcn. (880->14), ass. (0->53)
t116 = sin(pkin(11));
t120 = cos(pkin(11));
t125 = sin(qJ(2));
t127 = cos(qJ(2));
t149 = t125 * t116 - t120 * t127;
t147 = pkin(9) + rSges(7,3);
t146 = qJ(4) + rSges(5,3);
t109 = pkin(2) * t127 + pkin(1);
t117 = sin(pkin(10));
t121 = cos(pkin(10));
t118 = sin(pkin(6));
t122 = cos(pkin(6));
t140 = t122 * t125;
t97 = pkin(2) * t140 + (-pkin(7) - qJ(3)) * t118;
t145 = t117 * t109 + t121 * t97;
t115 = sin(pkin(12));
t144 = t115 * t122;
t143 = t117 * t118;
t142 = t118 * t121;
t139 = t122 * t127;
t137 = t122 * pkin(7) + qJ(1);
t136 = t115 * t143;
t135 = t115 * t142;
t103 = t121 * t109;
t134 = -t117 * t97 + t103;
t133 = t118 * t125 * pkin(2) + t122 * qJ(3) + t137;
t132 = t116 * t127 + t125 * t120;
t119 = cos(pkin(12));
t108 = pkin(4) * t119 + pkin(3);
t123 = -pkin(8) - qJ(4);
t130 = t149 * t122;
t86 = t117 * t130 - t121 * t132;
t96 = t132 * t122;
t87 = -t117 * t96 - t121 * t149;
t131 = pkin(4) * t136 + t87 * t108 + t86 * t123 + t134;
t94 = t149 * t118;
t95 = t132 * t118;
t129 = pkin(4) * t144 + t95 * t108 - t94 * t123 + t133;
t84 = -t117 * t132 - t121 * t130;
t85 = -t117 * t149 + t121 * t96;
t128 = -pkin(4) * t135 + t85 * t108 + t84 * t123 + t145;
t126 = cos(qJ(6));
t124 = sin(qJ(6));
t114 = pkin(12) + qJ(5);
t111 = cos(t114);
t110 = sin(t114);
t89 = t110 * t122 + t111 * t95;
t88 = t110 * t95 - t122 * t111;
t79 = t110 * t143 + t111 * t87;
t78 = t110 * t87 - t111 * t143;
t77 = -t110 * t142 + t111 * t85;
t76 = t110 * t85 + t111 * t142;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t121 - rSges(2,2) * t117) + g(2) * (rSges(2,1) * t117 + rSges(2,2) * t121) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t121 * pkin(1) + (-t117 * t140 + t121 * t127) * rSges(3,1) + (-t117 * t139 - t121 * t125) * rSges(3,2)) + g(2) * (t117 * pkin(1) + (t117 * t127 + t121 * t140) * rSges(3,1) + (-t117 * t125 + t121 * t139) * rSges(3,2)) + g(3) * (t122 * rSges(3,3) + t137) + (g(3) * (rSges(3,1) * t125 + rSges(3,2) * t127) + (g(1) * t117 - g(2) * t121) * (rSges(3,3) + pkin(7))) * t118) - m(4) * (g(1) * (rSges(4,1) * t87 + rSges(4,2) * t86 + t103 + (rSges(4,3) * t118 - t97) * t117) + g(2) * (rSges(4,1) * t85 + rSges(4,2) * t84 - rSges(4,3) * t142 + t145) + g(3) * (rSges(4,1) * t95 - rSges(4,2) * t94 + rSges(4,3) * t122 + t133)) - m(5) * (g(1) * (t87 * pkin(3) + (t119 * t87 + t136) * rSges(5,1) + (-t115 * t87 + t119 * t143) * rSges(5,2) - t146 * t86 + t134) + g(2) * (t85 * pkin(3) + (t119 * t85 - t135) * rSges(5,1) + (-t115 * t85 - t119 * t142) * rSges(5,2) - t146 * t84 + t145) + g(3) * (t95 * pkin(3) + (t119 * t95 + t144) * rSges(5,1) + (-t115 * t95 + t119 * t122) * rSges(5,2) + t146 * t94 + t133)) - m(6) * (g(1) * (rSges(6,1) * t79 - rSges(6,2) * t78 - rSges(6,3) * t86 + t131) + g(2) * (rSges(6,1) * t77 - rSges(6,2) * t76 - rSges(6,3) * t84 + t128) + g(3) * (rSges(6,1) * t89 - rSges(6,2) * t88 + rSges(6,3) * t94 + t129)) - m(7) * (g(1) * (t79 * pkin(5) + (-t124 * t86 + t126 * t79) * rSges(7,1) + (-t124 * t79 - t126 * t86) * rSges(7,2) + t147 * t78 + t131) + g(2) * (t77 * pkin(5) + (-t124 * t84 + t126 * t77) * rSges(7,1) + (-t124 * t77 - t126 * t84) * rSges(7,2) + t147 * t76 + t128) + g(3) * (t89 * pkin(5) + (t124 * t94 + t126 * t89) * rSges(7,1) + (-t124 * t89 + t126 * t94) * rSges(7,2) + t147 * t88 + t129));
U  = t1;
