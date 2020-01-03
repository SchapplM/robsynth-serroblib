% Calculate potential energy for
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR10_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR10_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR10_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR10_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:23:42
% EndTime: 2019-12-31 20:23:42
% DurationCPUTime: 0.42s
% Computational Cost: add. (253->100), mult. (549->142), div. (0->0), fcn. (668->12), ass. (0->48)
t124 = rSges(5,3) + pkin(8);
t98 = cos(pkin(5));
t123 = t98 * pkin(7) + pkin(6);
t122 = pkin(9) + rSges(6,3);
t102 = sin(qJ(1));
t106 = cos(qJ(1));
t101 = sin(qJ(2));
t96 = sin(pkin(5));
t83 = t98 * t101 * pkin(2) + (-pkin(7) - qJ(3)) * t96;
t105 = cos(qJ(2));
t92 = t105 * pkin(2) + pkin(1);
t121 = t102 * t92 + t106 * t83;
t120 = t101 * t96;
t119 = t102 * t96;
t97 = cos(pkin(10));
t118 = t105 * t97;
t117 = t106 * t96;
t116 = t102 * t101;
t115 = t102 * t105;
t114 = t106 * t101;
t113 = t106 * t105;
t95 = sin(pkin(10));
t108 = t101 * t97 + t105 * t95;
t82 = t108 * t98;
t85 = -t101 * t95 + t118;
t72 = t102 * t85 + t106 * t82;
t112 = t72 * pkin(3) + t121;
t111 = pkin(2) * t120 + t98 * qJ(3) + t123;
t81 = t108 * t96;
t110 = t81 * pkin(3) + t111;
t74 = -t102 * t82 + t106 * t85;
t88 = t106 * t92;
t109 = t74 * pkin(3) - t102 * t83 + t88;
t107 = t85 * t98;
t104 = cos(qJ(4));
t103 = cos(qJ(5));
t100 = sin(qJ(4));
t99 = sin(qJ(5));
t80 = -t96 * t118 + t95 * t120;
t76 = t98 * t100 + t81 * t104;
t75 = t81 * t100 - t98 * t104;
t73 = -t102 * t107 - t106 * t108;
t71 = -t102 * t108 + t106 * t107;
t68 = t100 * t119 + t74 * t104;
t67 = t74 * t100 - t104 * t119;
t66 = -t100 * t117 + t72 * t104;
t65 = t72 * t100 + t104 * t117;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t106 * rSges(2,1) - t102 * rSges(2,2)) + g(2) * (t102 * rSges(2,1) + t106 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t106 * pkin(1) + (-t98 * t116 + t113) * rSges(3,1) + (-t98 * t115 - t114) * rSges(3,2)) + g(2) * (t102 * pkin(1) + (t98 * t114 + t115) * rSges(3,1) + (t98 * t113 - t116) * rSges(3,2)) + g(3) * (t98 * rSges(3,3) + t123) + (g(3) * (rSges(3,1) * t101 + rSges(3,2) * t105) + (g(1) * t102 - g(2) * t106) * (rSges(3,3) + pkin(7))) * t96) - m(4) * (g(1) * (t74 * rSges(4,1) + t73 * rSges(4,2) + t88 + (rSges(4,3) * t96 - t83) * t102) + g(2) * (t72 * rSges(4,1) + t71 * rSges(4,2) - rSges(4,3) * t117 + t121) + g(3) * (t81 * rSges(4,1) - t80 * rSges(4,2) + t98 * rSges(4,3) + t111)) - m(5) * (g(1) * (t68 * rSges(5,1) - t67 * rSges(5,2) - t124 * t73 + t109) + g(2) * (t66 * rSges(5,1) - t65 * rSges(5,2) - t124 * t71 + t112) + g(3) * (t76 * rSges(5,1) - t75 * rSges(5,2) + t124 * t80 + t110)) - m(6) * (g(1) * (t68 * pkin(4) - t73 * pkin(8) + (t68 * t103 - t73 * t99) * rSges(6,1) + (-t73 * t103 - t68 * t99) * rSges(6,2) + t122 * t67 + t109) + g(2) * (t66 * pkin(4) - t71 * pkin(8) + (t66 * t103 - t71 * t99) * rSges(6,1) + (-t71 * t103 - t66 * t99) * rSges(6,2) + t122 * t65 + t112) + g(3) * (t76 * pkin(4) + t80 * pkin(8) + (t76 * t103 + t80 * t99) * rSges(6,1) + (t80 * t103 - t76 * t99) * rSges(6,2) + t122 * t75 + t110));
U = t1;
