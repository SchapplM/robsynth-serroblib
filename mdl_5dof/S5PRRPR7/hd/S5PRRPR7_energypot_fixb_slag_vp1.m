% Calculate potential energy for
% S5PRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRPR7_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR7_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR7_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR7_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:34:54
% EndTime: 2019-12-05 16:34:54
% DurationCPUTime: 0.30s
% Computational Cost: add. (243->100), mult. (536->139), div. (0->0), fcn. (651->12), ass. (0->48)
t97 = sin(pkin(5));
t124 = pkin(6) * t97;
t123 = rSges(4,3) + pkin(7);
t122 = pkin(8) + rSges(6,3);
t96 = sin(pkin(9));
t99 = cos(pkin(9));
t121 = t99 * pkin(1) + t96 * t124;
t102 = sin(qJ(3));
t120 = t102 * t97;
t103 = sin(qJ(2));
t119 = t103 * t97;
t105 = cos(qJ(3));
t118 = t105 * t97;
t106 = cos(qJ(2));
t117 = t97 * t106;
t116 = rSges(5,3) + qJ(4);
t100 = cos(pkin(5));
t115 = t100 * pkin(6) + qJ(1);
t114 = t100 * t103;
t113 = t100 * t106;
t84 = t99 * t106 - t96 * t114;
t112 = t84 * pkin(2) + t121;
t111 = pkin(2) * t119 + t115;
t82 = t96 * t106 + t99 * t114;
t92 = t96 * pkin(1);
t110 = t82 * pkin(2) - t99 * t124 + t92;
t75 = t84 * t105 + t96 * t120;
t83 = t99 * t103 + t96 * t113;
t109 = t75 * pkin(3) + t83 * pkin(7) + t112;
t86 = t100 * t102 + t103 * t118;
t108 = t86 * pkin(3) - pkin(7) * t117 + t111;
t73 = t82 * t105 - t99 * t120;
t81 = t96 * t103 - t99 * t113;
t107 = t73 * pkin(3) + t81 * pkin(7) + t110;
t104 = cos(qJ(5));
t101 = sin(qJ(5));
t98 = cos(pkin(10));
t95 = sin(pkin(10));
t85 = -t100 * t105 + t102 * t119;
t74 = t84 * t102 - t96 * t118;
t72 = t82 * t102 + t99 * t118;
t71 = -t95 * t117 + t86 * t98;
t70 = t98 * t117 + t86 * t95;
t67 = t75 * t98 + t83 * t95;
t66 = t75 * t95 - t83 * t98;
t65 = t73 * t98 + t81 * t95;
t64 = t73 * t95 - t81 * t98;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t99 * rSges(2,1) - t96 * rSges(2,2)) + g(2) * (t96 * rSges(2,1) + t99 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t84 * rSges(3,1) - t83 * rSges(3,2) + t121) + g(2) * (t82 * rSges(3,1) - t81 * rSges(3,2) + t92) + g(3) * (t100 * rSges(3,3) + t115) + (g(1) * rSges(3,3) * t96 + g(3) * (rSges(3,1) * t103 + rSges(3,2) * t106) + g(2) * (-rSges(3,3) - pkin(6)) * t99) * t97) - m(4) * (g(1) * (t75 * rSges(4,1) - t74 * rSges(4,2) + t123 * t83 + t112) + g(2) * (t73 * rSges(4,1) - t72 * rSges(4,2) + t123 * t81 + t110) + g(3) * (t86 * rSges(4,1) - t85 * rSges(4,2) - t123 * t117 + t111)) - m(5) * (g(1) * (t67 * rSges(5,1) - t66 * rSges(5,2) + t116 * t74 + t109) + g(2) * (t65 * rSges(5,1) - t64 * rSges(5,2) + t116 * t72 + t107) + g(3) * (t71 * rSges(5,1) - t70 * rSges(5,2) + t116 * t85 + t108)) - m(6) * (g(1) * (t67 * pkin(4) + t74 * qJ(4) + (t74 * t101 + t67 * t104) * rSges(6,1) + (-t67 * t101 + t74 * t104) * rSges(6,2) + t122 * t66 + t109) + g(2) * (t65 * pkin(4) + t72 * qJ(4) + (t72 * t101 + t65 * t104) * rSges(6,1) + (-t65 * t101 + t72 * t104) * rSges(6,2) + t122 * t64 + t107) + g(3) * (t71 * pkin(4) + t85 * qJ(4) + (t85 * t101 + t71 * t104) * rSges(6,1) + (-t71 * t101 + t85 * t104) * rSges(6,2) + t122 * t70 + t108));
U = t1;
