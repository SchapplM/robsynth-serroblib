% Calculate potential energy for
% S6PRPRRP1
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
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRP1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP1_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP1_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:55:07
% EndTime: 2019-03-08 19:55:08
% DurationCPUTime: 0.54s
% Computational Cost: add. (388->120), mult. (853->163), div. (0->0), fcn. (1059->12), ass. (0->53)
t108 = sin(pkin(11));
t116 = sin(qJ(2));
t119 = cos(qJ(2));
t133 = cos(pkin(11));
t97 = -t116 * t108 + t119 * t133;
t137 = rSges(5,3) + pkin(8);
t136 = rSges(6,3) + pkin(9);
t105 = pkin(2) * t119 + pkin(1);
t109 = sin(pkin(10));
t111 = cos(pkin(10));
t110 = sin(pkin(6));
t112 = cos(pkin(6));
t130 = t112 * t116;
t95 = pkin(2) * t130 + (-pkin(7) - qJ(3)) * t110;
t135 = t109 * t105 + t111 * t95;
t134 = rSges(7,3) + qJ(6) + pkin(9);
t132 = t109 * t110;
t131 = t111 * t110;
t129 = t112 * t119;
t127 = t112 * pkin(7) + qJ(1);
t96 = -t119 * t108 - t116 * t133;
t94 = t96 * t112;
t84 = t109 * t97 - t111 * t94;
t126 = t84 * pkin(3) + t135;
t114 = sin(qJ(5));
t125 = pkin(5) * t114 + pkin(8);
t123 = t110 * t116 * pkin(2) + t112 * qJ(3) + t127;
t86 = t109 * t94 + t111 * t97;
t99 = t111 * t105;
t122 = t86 * pkin(3) - t109 * t95 + t99;
t93 = t96 * t110;
t121 = -t93 * pkin(3) + t123;
t120 = t112 * t97;
t118 = cos(qJ(4));
t117 = cos(qJ(5));
t115 = sin(qJ(4));
t104 = pkin(5) * t117 + pkin(4);
t92 = t97 * t110;
t88 = t112 * t115 - t118 * t93;
t87 = -t112 * t118 - t115 * t93;
t85 = -t109 * t120 + t111 * t96;
t83 = t109 * t96 + t111 * t120;
t80 = t115 * t132 + t118 * t86;
t79 = t115 * t86 - t118 * t132;
t78 = -t115 * t131 + t118 * t84;
t77 = t115 * t84 + t118 * t131;
t76 = -t114 * t92 + t117 * t88;
t75 = -t114 * t88 - t117 * t92;
t74 = -t114 * t85 + t117 * t80;
t73 = -t114 * t80 - t117 * t85;
t72 = -t114 * t83 + t117 * t78;
t71 = -t114 * t78 - t117 * t83;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t111 - rSges(2,2) * t109) + g(2) * (rSges(2,1) * t109 + rSges(2,2) * t111) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t111 * pkin(1) + (-t109 * t130 + t111 * t119) * rSges(3,1) + (-t109 * t129 - t111 * t116) * rSges(3,2)) + g(2) * (t109 * pkin(1) + (t109 * t119 + t111 * t130) * rSges(3,1) + (-t109 * t116 + t111 * t129) * rSges(3,2)) + g(3) * (t112 * rSges(3,3) + t127) + (g(3) * (rSges(3,1) * t116 + rSges(3,2) * t119) + (g(1) * t109 - g(2) * t111) * (rSges(3,3) + pkin(7))) * t110) - m(4) * (g(1) * (rSges(4,1) * t86 + rSges(4,2) * t85 + t99 + (rSges(4,3) * t110 - t95) * t109) + g(2) * (rSges(4,1) * t84 + rSges(4,2) * t83 - rSges(4,3) * t131 + t135) + g(3) * (-rSges(4,1) * t93 + rSges(4,2) * t92 + rSges(4,3) * t112 + t123)) - m(5) * (g(1) * (rSges(5,1) * t80 - rSges(5,2) * t79 - t137 * t85 + t122) + g(2) * (rSges(5,1) * t78 - rSges(5,2) * t77 - t137 * t83 + t126) + g(3) * (rSges(5,1) * t88 - rSges(5,2) * t87 - t137 * t92 + t121)) - m(6) * (g(1) * (rSges(6,1) * t74 + rSges(6,2) * t73 + pkin(4) * t80 - pkin(8) * t85 + t136 * t79 + t122) + g(2) * (rSges(6,1) * t72 + rSges(6,2) * t71 + pkin(4) * t78 - pkin(8) * t83 + t136 * t77 + t126) + g(3) * (rSges(6,1) * t76 + rSges(6,2) * t75 + pkin(4) * t88 - pkin(8) * t92 + t136 * t87 + t121)) - m(7) * (g(1) * (rSges(7,1) * t74 + rSges(7,2) * t73 + t104 * t80 - t125 * t85 + t134 * t79 + t122) + g(2) * (rSges(7,1) * t72 + rSges(7,2) * t71 + t104 * t78 - t125 * t83 + t134 * t77 + t126) + g(3) * (rSges(7,1) * t76 + rSges(7,2) * t75 + t104 * t88 - t125 * t92 + t134 * t87 + t121));
U  = t1;
