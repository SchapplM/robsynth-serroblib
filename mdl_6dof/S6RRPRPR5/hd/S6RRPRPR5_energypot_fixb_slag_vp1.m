% Calculate potential energy for
% S6RRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPR5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR5_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR5_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:28:52
% EndTime: 2019-03-09 10:28:52
% DurationCPUTime: 0.55s
% Computational Cost: add. (400->119), mult. (853->159), div. (0->0), fcn. (1059->14), ass. (0->53)
t107 = sin(pkin(11));
t113 = sin(qJ(2));
t116 = cos(qJ(2));
t133 = cos(pkin(11));
t92 = -t113 * t107 + t116 * t133;
t139 = rSges(5,3) + pkin(9);
t138 = pkin(2) * t113;
t110 = cos(pkin(6));
t137 = t110 * pkin(8) + pkin(7);
t100 = t116 * pkin(2) + pkin(1);
t114 = sin(qJ(1));
t117 = cos(qJ(1));
t108 = sin(pkin(6));
t90 = t110 * t138 + (-pkin(8) - qJ(3)) * t108;
t136 = t114 * t100 + t117 * t90;
t135 = qJ(5) + rSges(6,3);
t134 = pkin(10) + qJ(5) + rSges(7,3);
t131 = t114 * t108;
t130 = t114 * t113;
t129 = t114 * t116;
t128 = t117 * t108;
t127 = t117 * t113;
t126 = t117 * t116;
t91 = -t116 * t107 - t113 * t133;
t89 = t91 * t110;
t79 = t114 * t92 - t117 * t89;
t125 = t79 * pkin(3) + t136;
t123 = t110 * qJ(3) + t108 * t138 + t137;
t81 = t114 * t89 + t117 * t92;
t95 = t117 * t100;
t122 = t81 * pkin(3) - t114 * t90 + t95;
t88 = t91 * t108;
t121 = -t88 * pkin(3) + t123;
t105 = pkin(12) + qJ(6);
t101 = sin(t105);
t102 = cos(t105);
t109 = cos(pkin(12));
t120 = t102 * rSges(7,1) - t101 * rSges(7,2) + t109 * pkin(5) + pkin(4);
t106 = sin(pkin(12));
t119 = t101 * rSges(7,1) + t102 * rSges(7,2) + t106 * pkin(5) + pkin(9);
t118 = t92 * t110;
t115 = cos(qJ(4));
t112 = sin(qJ(4));
t87 = t92 * t108;
t83 = t110 * t112 - t88 * t115;
t82 = -t110 * t115 - t88 * t112;
t80 = -t114 * t118 + t117 * t91;
t78 = t114 * t91 + t117 * t118;
t75 = t112 * t131 + t81 * t115;
t74 = t81 * t112 - t115 * t131;
t73 = -t112 * t128 + t79 * t115;
t72 = t79 * t112 + t115 * t128;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t117 * rSges(2,1) - t114 * rSges(2,2)) + g(2) * (t114 * rSges(2,1) + t117 * rSges(2,2)) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (t117 * pkin(1) + (-t110 * t130 + t126) * rSges(3,1) + (-t110 * t129 - t127) * rSges(3,2)) + g(2) * (t114 * pkin(1) + (t110 * t127 + t129) * rSges(3,1) + (t110 * t126 - t130) * rSges(3,2)) + g(3) * (t110 * rSges(3,3) + t137) + (g(3) * (rSges(3,1) * t113 + rSges(3,2) * t116) + (g(1) * t114 - g(2) * t117) * (rSges(3,3) + pkin(8))) * t108) - m(4) * (g(1) * (t81 * rSges(4,1) + t80 * rSges(4,2) + t95 + (rSges(4,3) * t108 - t90) * t114) + g(2) * (t79 * rSges(4,1) + t78 * rSges(4,2) - rSges(4,3) * t128 + t136) + g(3) * (-t88 * rSges(4,1) + t87 * rSges(4,2) + t110 * rSges(4,3) + t123)) - m(5) * (g(1) * (t75 * rSges(5,1) - t74 * rSges(5,2) - t139 * t80 + t122) + g(2) * (t73 * rSges(5,1) - t72 * rSges(5,2) - t139 * t78 + t125) + g(3) * (t83 * rSges(5,1) - t82 * rSges(5,2) - t139 * t87 + t121)) - m(6) * (g(1) * (t75 * pkin(4) - t80 * pkin(9) + (-t80 * t106 + t75 * t109) * rSges(6,1) + (-t75 * t106 - t80 * t109) * rSges(6,2) + t135 * t74 + t122) + g(2) * (t73 * pkin(4) - t78 * pkin(9) + (-t78 * t106 + t73 * t109) * rSges(6,1) + (-t73 * t106 - t78 * t109) * rSges(6,2) + t135 * t72 + t125) + g(3) * (t83 * pkin(4) - t87 * pkin(9) + (-t87 * t106 + t83 * t109) * rSges(6,1) + (-t83 * t106 - t87 * t109) * rSges(6,2) + t135 * t82 + t121)) - m(7) * (g(1) * (-t119 * t80 + t120 * t75 + t134 * t74 + t122) + g(2) * (-t119 * t78 + t120 * t73 + t134 * t72 + t125) + g(3) * (-t119 * t87 + t120 * t83 + t134 * t82 + t121));
U  = t1;
