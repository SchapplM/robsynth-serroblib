% Calculate potential energy for
% S6PRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRP4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP4_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP4_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:09:08
% EndTime: 2019-03-08 20:09:09
% DurationCPUTime: 0.37s
% Computational Cost: add. (370->118), mult. (591->153), div. (0->0), fcn. (697->12), ass. (0->54)
t146 = rSges(7,1) + pkin(5);
t145 = rSges(7,2) + pkin(9);
t144 = rSges(6,3) + pkin(9);
t115 = sin(pkin(11));
t143 = pkin(3) * t115;
t142 = rSges(7,3) + qJ(6);
t141 = qJ(3) + rSges(4,3);
t116 = sin(pkin(10));
t117 = sin(pkin(6));
t140 = t116 * t117;
t119 = cos(pkin(10));
t139 = t117 * t119;
t123 = sin(qJ(2));
t138 = t117 * t123;
t125 = cos(qJ(2));
t137 = t117 * t125;
t120 = cos(pkin(6));
t136 = t120 * t123;
t135 = t120 * t125;
t134 = t120 * pkin(7) + qJ(1);
t133 = t119 * pkin(1) + pkin(7) * t140;
t132 = t115 * t140;
t118 = cos(pkin(11));
t108 = pkin(3) * t118 + pkin(2);
t121 = -pkin(8) - qJ(3);
t98 = t116 * t135 + t119 * t123;
t99 = -t116 * t136 + t119 * t125;
t131 = pkin(3) * t132 + t99 * t108 - t98 * t121 + t133;
t130 = t108 * t138 + t120 * t143 + t121 * t137 + t134;
t114 = pkin(11) + qJ(4);
t109 = sin(t114);
t110 = cos(t114);
t84 = t109 * t140 + t99 * t110;
t129 = t84 * pkin(4) + t131;
t91 = t109 * t120 + t110 * t138;
t128 = t91 * pkin(4) + t130;
t111 = t116 * pkin(1);
t96 = t116 * t123 - t119 * t135;
t97 = t116 * t125 + t119 * t136;
t127 = t111 + t97 * t108 + (-pkin(7) - t143) * t139 - t96 * t121;
t82 = -t109 * t139 + t97 * t110;
t126 = t82 * pkin(4) + t127;
t124 = cos(qJ(5));
t122 = sin(qJ(5));
t90 = t109 * t138 - t120 * t110;
t86 = -t122 * t137 + t91 * t124;
t85 = t91 * t122 + t124 * t137;
t83 = t109 * t99 - t110 * t140;
t81 = t97 * t109 + t110 * t139;
t78 = t122 * t98 + t124 * t84;
t77 = t122 * t84 - t98 * t124;
t76 = t122 * t96 + t124 * t82;
t75 = t122 * t82 - t96 * t124;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t119 - rSges(2,2) * t116) + g(2) * (rSges(2,1) * t116 + rSges(2,2) * t119) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t99 * rSges(3,1) - t98 * rSges(3,2) + t133) + g(2) * (t97 * rSges(3,1) - t96 * rSges(3,2) + t111) + g(3) * (t120 * rSges(3,3) + t134) + (g(1) * rSges(3,3) * t116 + g(3) * (rSges(3,1) * t123 + rSges(3,2) * t125) + g(2) * (-rSges(3,3) - pkin(7)) * t119) * t117) - m(4) * (g(1) * (t99 * pkin(2) + (t118 * t99 + t132) * rSges(4,1) + (-t99 * t115 + t118 * t140) * rSges(4,2) + t141 * t98 + t133) + g(2) * (t97 * pkin(2) + t111 - pkin(7) * t139 + (-t115 * t139 + t118 * t97) * rSges(4,1) + (-t97 * t115 - t118 * t139) * rSges(4,2) + t141 * t96) + g(3) * ((t115 * rSges(4,1) + t118 * rSges(4,2)) * t120 + (-t141 * t125 + (rSges(4,1) * t118 - rSges(4,2) * t115 + pkin(2)) * t123) * t117 + t134)) - m(5) * (g(1) * (rSges(5,1) * t84 - rSges(5,2) * t83 + rSges(5,3) * t98 + t131) + g(2) * (t82 * rSges(5,1) - t81 * rSges(5,2) + t96 * rSges(5,3) + t127) + g(3) * (t91 * rSges(5,1) - t90 * rSges(5,2) - rSges(5,3) * t137 + t130)) - m(6) * (g(1) * (rSges(6,1) * t78 - rSges(6,2) * t77 + t144 * t83 + t129) + g(2) * (t76 * rSges(6,1) - t75 * rSges(6,2) + t144 * t81 + t126) + g(3) * (rSges(6,1) * t86 - rSges(6,2) * t85 + t144 * t90 + t128)) - m(7) * (g(1) * (t142 * t77 + t145 * t83 + t146 * t78 + t129) + g(2) * (t142 * t75 + t145 * t81 + t146 * t76 + t126) + g(3) * (t142 * t85 + t145 * t90 + t146 * t86 + t128));
U  = t1;
