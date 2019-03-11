% Calculate potential energy for
% S6RRPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 14:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR10_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR10_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR10_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR10_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:18:14
% EndTime: 2019-03-09 14:18:14
% DurationCPUTime: 0.42s
% Computational Cost: add. (355->130), mult. (545->170), div. (0->0), fcn. (633->14), ass. (0->52)
t105 = cos(pkin(6));
t134 = t105 * pkin(8) + pkin(7);
t133 = pkin(10) + rSges(6,3);
t102 = sin(pkin(12));
t132 = pkin(3) * t102;
t112 = cos(qJ(1));
t103 = sin(pkin(6));
t109 = sin(qJ(1));
t125 = t103 * t109;
t131 = t112 * pkin(1) + pkin(8) * t125;
t107 = sin(qJ(5));
t111 = cos(qJ(2));
t119 = t111 * t112;
t108 = sin(qJ(2));
t122 = t108 * t109;
t79 = -t105 * t119 + t122;
t130 = t107 * t79;
t120 = t109 * t111;
t121 = t108 * t112;
t81 = t105 * t120 + t121;
t129 = t107 * t81;
t128 = qJ(3) + rSges(4,3);
t127 = pkin(11) + pkin(10) + rSges(7,3);
t126 = t103 * t108;
t124 = t103 * t111;
t123 = t103 * t112;
t118 = t102 * t125;
t117 = t107 * t124;
t106 = -pkin(9) - qJ(3);
t104 = cos(pkin(12));
t91 = pkin(3) * t104 + pkin(2);
t116 = t105 * t132 + t106 * t124 + t91 * t126 + t134;
t82 = -t105 * t122 + t119;
t115 = pkin(3) * t118 - t81 * t106 + t82 * t91 + t131;
t80 = t105 * t121 + t120;
t98 = t109 * pkin(1);
t114 = t80 * t91 + t98 + (-pkin(8) - t132) * t123 - t79 * t106;
t110 = cos(qJ(5));
t101 = qJ(5) + qJ(6);
t100 = pkin(12) + qJ(4);
t96 = cos(t101);
t95 = sin(t101);
t94 = cos(t100);
t93 = sin(t100);
t92 = pkin(5) * t110 + pkin(4);
t76 = t105 * t93 + t126 * t94;
t75 = -t105 * t94 + t126 * t93;
t72 = t125 * t93 + t82 * t94;
t71 = -t125 * t94 + t82 * t93;
t70 = -t123 * t93 + t80 * t94;
t69 = t123 * t94 + t80 * t93;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t112 - rSges(2,2) * t109) + g(2) * (rSges(2,1) * t109 + rSges(2,2) * t112) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t82 - rSges(3,2) * t81 + t131) + g(2) * (rSges(3,1) * t80 - rSges(3,2) * t79 + t98) + g(3) * (rSges(3,3) * t105 + t134) + (g(1) * rSges(3,3) * t109 + g(3) * (rSges(3,1) * t108 + rSges(3,2) * t111) + g(2) * (-rSges(3,3) - pkin(8)) * t112) * t103) - m(4) * (g(1) * (t82 * pkin(2) + (t104 * t82 + t118) * rSges(4,1) + (-t102 * t82 + t104 * t125) * rSges(4,2) + t128 * t81 + t131) + g(2) * (t80 * pkin(2) + t98 - pkin(8) * t123 + (-t102 * t123 + t104 * t80) * rSges(4,1) + (-t102 * t80 - t104 * t123) * rSges(4,2) + t128 * t79) + g(3) * ((t102 * rSges(4,1) + t104 * rSges(4,2)) * t105 + (-t128 * t111 + (t104 * rSges(4,1) - t102 * rSges(4,2) + pkin(2)) * t108) * t103 + t134)) - m(5) * (g(1) * (rSges(5,1) * t72 - rSges(5,2) * t71 + rSges(5,3) * t81 + t115) + g(2) * (rSges(5,1) * t70 - rSges(5,2) * t69 + rSges(5,3) * t79 + t114) + g(3) * (rSges(5,1) * t76 - rSges(5,2) * t75 - rSges(5,3) * t124 + t116)) - m(6) * (g(1) * (t72 * pkin(4) + (t110 * t72 + t129) * rSges(6,1) + (-t107 * t72 + t110 * t81) * rSges(6,2) + t133 * t71 + t115) + g(2) * (t70 * pkin(4) + (t110 * t70 + t130) * rSges(6,1) + (-t107 * t70 + t110 * t79) * rSges(6,2) + t133 * t69 + t114) + g(3) * (t76 * pkin(4) + (t110 * t76 - t117) * rSges(6,1) + (-t107 * t76 - t110 * t124) * rSges(6,2) + t133 * t75 + t116)) - m(7) * (g(1) * (t72 * t92 + pkin(5) * t129 + (t72 * t96 + t81 * t95) * rSges(7,1) + (-t72 * t95 + t81 * t96) * rSges(7,2) + t127 * t71 + t115) + g(2) * (t70 * t92 + pkin(5) * t130 + (t70 * t96 + t79 * t95) * rSges(7,1) + (-t70 * t95 + t79 * t96) * rSges(7,2) + t127 * t69 + t114) + g(3) * (t76 * t92 - pkin(5) * t117 + (-t124 * t95 + t76 * t96) * rSges(7,1) + (-t124 * t96 - t76 * t95) * rSges(7,2) + t127 * t75 + t116));
U  = t1;
