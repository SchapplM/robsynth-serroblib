% Calculate potential energy for
% S5RRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRP11_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP11_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP11_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP11_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:14:20
% EndTime: 2019-12-31 22:14:21
% DurationCPUTime: 0.28s
% Computational Cost: add. (222->90), mult. (481->119), div. (0->0), fcn. (575->10), ass. (0->50)
t128 = rSges(6,1) + pkin(4);
t127 = rSges(6,2) + pkin(9);
t126 = rSges(4,3) + pkin(8);
t125 = rSges(5,3) + pkin(9);
t98 = cos(pkin(5));
t124 = t98 * pkin(7) + pkin(6);
t106 = cos(qJ(1));
t102 = sin(qJ(1));
t97 = sin(pkin(5));
t121 = t102 * t97;
t123 = t106 * pkin(1) + pkin(7) * t121;
t101 = sin(qJ(2));
t122 = t101 * t97;
t104 = cos(qJ(3));
t120 = t104 * t97;
t105 = cos(qJ(2));
t119 = t105 * t97;
t118 = t106 * t97;
t117 = rSges(6,3) + qJ(5);
t116 = t102 * t101;
t115 = t102 * t105;
t114 = t106 * t101;
t113 = t106 * t105;
t112 = pkin(2) * t122 + t124;
t88 = -t116 * t98 + t113;
t111 = t88 * pkin(2) + t123;
t86 = t114 * t98 + t115;
t95 = t102 * pkin(1);
t110 = t86 * pkin(2) - pkin(7) * t118 + t95;
t100 = sin(qJ(3));
t77 = t100 * t121 + t104 * t88;
t87 = t115 * t98 + t114;
t109 = t77 * pkin(3) + t87 * pkin(8) + t111;
t84 = t100 * t98 + t101 * t120;
t108 = t84 * pkin(3) - pkin(8) * t119 + t112;
t75 = -t100 * t118 + t104 * t86;
t85 = -t113 * t98 + t116;
t107 = t75 * pkin(3) + t85 * pkin(8) + t110;
t103 = cos(qJ(4));
t99 = sin(qJ(4));
t83 = t100 * t122 - t98 * t104;
t76 = t100 * t88 - t102 * t120;
t74 = t100 * t86 + t104 * t118;
t73 = t103 * t84 - t119 * t99;
t72 = t103 * t119 + t84 * t99;
t69 = t103 * t77 + t87 * t99;
t68 = -t87 * t103 + t77 * t99;
t67 = t103 * t75 + t85 * t99;
t66 = -t85 * t103 + t75 * t99;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + rSges(1,3) * g(3)) - m(2) * (g(1) * (rSges(2,1) * t106 - rSges(2,2) * t102) + g(2) * (rSges(2,1) * t102 + rSges(2,2) * t106) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t88 - rSges(3,2) * t87 + t123) + g(2) * (t86 * rSges(3,1) - t85 * rSges(3,2) + t95) + g(3) * (t98 * rSges(3,3) + t124) + (g(1) * rSges(3,3) * t102 + g(3) * (rSges(3,1) * t101 + rSges(3,2) * t105) + g(2) * (-rSges(3,3) - pkin(7)) * t106) * t97) - m(4) * (g(1) * (t77 * rSges(4,1) - t76 * rSges(4,2) + t126 * t87 + t111) + g(2) * (t75 * rSges(4,1) - t74 * rSges(4,2) + t126 * t85 + t110) + g(3) * (t84 * rSges(4,1) - t83 * rSges(4,2) - t119 * t126 + t112)) - m(5) * (g(1) * (t69 * rSges(5,1) - t68 * rSges(5,2) + t125 * t76 + t109) + g(2) * (t67 * rSges(5,1) - t66 * rSges(5,2) + t125 * t74 + t107) + g(3) * (t73 * rSges(5,1) - t72 * rSges(5,2) + t125 * t83 + t108)) - m(6) * (g(1) * (t117 * t68 + t127 * t76 + t128 * t69 + t109) + g(2) * (t117 * t66 + t127 * t74 + t128 * t67 + t107) + g(3) * (t117 * t72 + t127 * t83 + t128 * t73 + t108));
U = t1;
