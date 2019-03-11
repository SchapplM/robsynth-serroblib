% Calculate potential energy for
% S6PRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRPR5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR5_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR5_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:42:52
% EndTime: 2019-03-08 19:42:52
% DurationCPUTime: 0.39s
% Computational Cost: add. (332->122), mult. (523->159), div. (0->0), fcn. (603->12), ass. (0->46)
t99 = sin(pkin(11));
t128 = pkin(3) * t99;
t127 = pkin(9) + rSges(7,3);
t103 = cos(pkin(10));
t100 = sin(pkin(10));
t101 = sin(pkin(6));
t122 = t100 * t101;
t126 = pkin(1) * t103 + pkin(7) * t122;
t125 = rSges(6,3) + qJ(5);
t104 = cos(pkin(6));
t124 = pkin(7) * t104 + qJ(1);
t123 = qJ(3) + rSges(4,3);
t121 = t101 * t103;
t107 = sin(qJ(2));
t120 = t101 * t107;
t109 = cos(qJ(2));
t119 = t101 * t109;
t118 = t104 * t107;
t117 = t104 * t109;
t116 = t99 * t122;
t105 = -pkin(8) - qJ(3);
t81 = t100 * t117 + t103 * t107;
t82 = -t100 * t118 + t103 * t109;
t102 = cos(pkin(11));
t92 = pkin(3) * t102 + pkin(2);
t115 = pkin(3) * t116 - t105 * t81 + t82 * t92 + t126;
t114 = t104 * t128 + t105 * t119 + t120 * t92 + t124;
t98 = pkin(11) + qJ(4);
t93 = sin(t98);
t94 = cos(t98);
t71 = t122 * t93 + t82 * t94;
t113 = pkin(4) * t71 + t115;
t76 = t104 * t93 + t120 * t94;
t112 = pkin(4) * t76 + t114;
t79 = t100 * t107 - t103 * t117;
t80 = t100 * t109 + t103 * t118;
t95 = t100 * pkin(1);
t111 = t80 * t92 + t95 + (-pkin(7) - t128) * t121 - t79 * t105;
t69 = -t121 * t93 + t80 * t94;
t110 = pkin(4) * t69 + t111;
t108 = cos(qJ(6));
t106 = sin(qJ(6));
t75 = -t104 * t94 + t120 * t93;
t70 = -t122 * t94 + t82 * t93;
t68 = t121 * t94 + t80 * t93;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t103 - rSges(2,2) * t100) + g(2) * (rSges(2,1) * t100 + rSges(2,2) * t103) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t82 - rSges(3,2) * t81 + t126) + g(2) * (t80 * rSges(3,1) - t79 * rSges(3,2) + t95) + g(3) * (t104 * rSges(3,3) + t124) + (g(1) * rSges(3,3) * t100 + g(3) * (rSges(3,1) * t107 + rSges(3,2) * t109) + g(2) * (-rSges(3,3) - pkin(7)) * t103) * t101) - m(4) * (g(1) * (t82 * pkin(2) + (t102 * t82 + t116) * rSges(4,1) + (t102 * t122 - t82 * t99) * rSges(4,2) + t123 * t81 + t126) + g(2) * (t80 * pkin(2) + t95 - pkin(7) * t121 + (t102 * t80 - t121 * t99) * rSges(4,1) + (-t102 * t121 - t80 * t99) * rSges(4,2) + t123 * t79) + g(3) * ((rSges(4,1) * t99 + rSges(4,2) * t102) * t104 + (-t123 * t109 + (rSges(4,1) * t102 - rSges(4,2) * t99 + pkin(2)) * t107) * t101 + t124)) - m(5) * (g(1) * (rSges(5,1) * t71 - rSges(5,2) * t70 + rSges(5,3) * t81 + t115) + g(2) * (t69 * rSges(5,1) - t68 * rSges(5,2) + t79 * rSges(5,3) + t111) + g(3) * (rSges(5,1) * t76 - rSges(5,2) * t75 - rSges(5,3) * t119 + t114)) - m(6) * (g(1) * (t81 * rSges(6,1) - t71 * rSges(6,2) + t125 * t70 + t113) + g(2) * (t79 * rSges(6,1) - t69 * rSges(6,2) + t125 * t68 + t110) + g(3) * (-rSges(6,1) * t119 - rSges(6,2) * t76 + t125 * t75 + t112)) - m(7) * (g(1) * (t81 * pkin(5) + t70 * qJ(5) + (t106 * t70 + t108 * t81) * rSges(7,1) + (-t106 * t81 + t108 * t70) * rSges(7,2) + t127 * t71 + t113) + g(2) * (t79 * pkin(5) + t68 * qJ(5) + (t106 * t68 + t108 * t79) * rSges(7,1) + (-t106 * t79 + t108 * t68) * rSges(7,2) + t127 * t69 + t110) + g(3) * (-pkin(5) * t119 + t75 * qJ(5) + (t106 * t75 - t108 * t119) * rSges(7,1) + (t106 * t119 + t108 * t75) * rSges(7,2) + t127 * t76 + t112));
U  = t1;
