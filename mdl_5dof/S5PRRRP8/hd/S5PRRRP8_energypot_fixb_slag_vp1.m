% Calculate potential energy for
% S5PRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRP8_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP8_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP8_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP8_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:58:08
% EndTime: 2019-12-05 16:58:08
% DurationCPUTime: 0.30s
% Computational Cost: add. (222->90), mult. (481->119), div. (0->0), fcn. (575->10), ass. (0->50)
t96 = sin(pkin(5));
t126 = pkin(6) * t96;
t125 = rSges(6,1) + pkin(4);
t124 = rSges(6,2) + pkin(8);
t123 = rSges(4,3) + pkin(7);
t122 = rSges(5,3) + pkin(8);
t95 = sin(pkin(9));
t97 = cos(pkin(9));
t121 = t97 * pkin(1) + t95 * t126;
t100 = sin(qJ(3));
t120 = t100 * t96;
t101 = sin(qJ(2));
t119 = t101 * t96;
t103 = cos(qJ(3));
t118 = t103 * t96;
t104 = cos(qJ(2));
t117 = t104 * t96;
t116 = t95 * t101;
t115 = t95 * t104;
t114 = t97 * t101;
t113 = t97 * t104;
t112 = rSges(6,3) + qJ(5);
t98 = cos(pkin(5));
t111 = t98 * pkin(6) + qJ(1);
t84 = -t98 * t116 + t113;
t110 = t84 * pkin(2) + t121;
t109 = pkin(2) * t119 + t111;
t82 = t98 * t114 + t115;
t92 = t95 * pkin(1);
t108 = t82 * pkin(2) - t97 * t126 + t92;
t73 = t84 * t103 + t95 * t120;
t83 = t98 * t115 + t114;
t107 = t73 * pkin(3) + t83 * pkin(7) + t110;
t86 = t98 * t100 + t101 * t118;
t106 = t86 * pkin(3) - pkin(7) * t117 + t109;
t71 = t82 * t103 - t97 * t120;
t81 = -t98 * t113 + t116;
t105 = t71 * pkin(3) + t81 * pkin(7) + t108;
t102 = cos(qJ(4));
t99 = sin(qJ(4));
t85 = t100 * t119 - t98 * t103;
t75 = t86 * t102 - t99 * t117;
t74 = t102 * t117 + t86 * t99;
t72 = t84 * t100 - t95 * t118;
t70 = t82 * t100 + t97 * t118;
t67 = t73 * t102 + t83 * t99;
t66 = -t83 * t102 + t73 * t99;
t65 = t71 * t102 + t81 * t99;
t64 = -t81 * t102 + t71 * t99;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t97 * rSges(2,1) - t95 * rSges(2,2)) + g(2) * (t95 * rSges(2,1) + t97 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t84 * rSges(3,1) - t83 * rSges(3,2) + t121) + g(2) * (t82 * rSges(3,1) - t81 * rSges(3,2) + t92) + g(3) * (t98 * rSges(3,3) + t111) + (g(1) * rSges(3,3) * t95 + g(3) * (rSges(3,1) * t101 + rSges(3,2) * t104) + g(2) * (-rSges(3,3) - pkin(6)) * t97) * t96) - m(4) * (g(1) * (t73 * rSges(4,1) - t72 * rSges(4,2) + t123 * t83 + t110) + g(2) * (t71 * rSges(4,1) - t70 * rSges(4,2) + t123 * t81 + t108) + g(3) * (t86 * rSges(4,1) - t85 * rSges(4,2) - t123 * t117 + t109)) - m(5) * (g(1) * (t67 * rSges(5,1) - t66 * rSges(5,2) + t122 * t72 + t107) + g(2) * (t65 * rSges(5,1) - t64 * rSges(5,2) + t122 * t70 + t105) + g(3) * (t75 * rSges(5,1) - t74 * rSges(5,2) + t122 * t85 + t106)) - m(6) * (g(1) * (t112 * t66 + t124 * t72 + t125 * t67 + t107) + g(2) * (t112 * t64 + t124 * t70 + t125 * t65 + t105) + g(3) * (t112 * t74 + t124 * t85 + t125 * t75 + t106));
U = t1;
