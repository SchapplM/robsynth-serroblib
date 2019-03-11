% Calculate potential energy for
% S6RRPRRR5
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
% Datum: 2019-03-09 13:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR5_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR5_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:37:56
% EndTime: 2019-03-09 13:37:56
% DurationCPUTime: 0.55s
% Computational Cost: add. (400->119), mult. (853->159), div. (0->0), fcn. (1059->14), ass. (0->53)
t106 = sin(pkin(12));
t111 = sin(qJ(2));
t115 = cos(qJ(2));
t133 = cos(pkin(12));
t92 = -t111 * t106 + t115 * t133;
t139 = rSges(5,3) + pkin(9);
t138 = pkin(10) + rSges(6,3);
t137 = pkin(2) * t111;
t108 = cos(pkin(6));
t136 = t108 * pkin(8) + pkin(7);
t100 = pkin(2) * t115 + pkin(1);
t112 = sin(qJ(1));
t116 = cos(qJ(1));
t107 = sin(pkin(6));
t90 = t108 * t137 + (-pkin(8) - qJ(3)) * t107;
t135 = t112 * t100 + t116 * t90;
t134 = pkin(11) + pkin(10) + rSges(7,3);
t132 = t107 * t112;
t131 = t107 * t116;
t129 = t111 * t112;
t128 = t111 * t116;
t127 = t112 * t115;
t126 = t115 * t116;
t91 = -t115 * t106 - t111 * t133;
t89 = t91 * t108;
t79 = t112 * t92 - t116 * t89;
t125 = t79 * pkin(3) + t135;
t123 = t108 * qJ(3) + t107 * t137 + t136;
t81 = t112 * t89 + t116 * t92;
t95 = t116 * t100;
t122 = t81 * pkin(3) - t112 * t90 + t95;
t88 = t91 * t107;
t121 = -t88 * pkin(3) + t123;
t105 = qJ(5) + qJ(6);
t102 = sin(t105);
t103 = cos(t105);
t113 = cos(qJ(5));
t120 = rSges(7,1) * t103 - rSges(7,2) * t102 + pkin(5) * t113 + pkin(4);
t109 = sin(qJ(5));
t119 = rSges(7,1) * t102 + rSges(7,2) * t103 + pkin(5) * t109 + pkin(9);
t118 = t108 * t92;
t114 = cos(qJ(4));
t110 = sin(qJ(4));
t87 = t92 * t107;
t83 = t108 * t110 - t114 * t88;
t82 = -t108 * t114 - t110 * t88;
t80 = -t112 * t118 + t116 * t91;
t78 = t112 * t91 + t116 * t118;
t75 = t110 * t132 + t114 * t81;
t74 = t110 * t81 - t114 * t132;
t73 = -t110 * t131 + t114 * t79;
t72 = t110 * t79 + t114 * t131;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t116 - rSges(2,2) * t112) + g(2) * (rSges(2,1) * t112 + rSges(2,2) * t116) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (t116 * pkin(1) + (-t108 * t129 + t126) * rSges(3,1) + (-t108 * t127 - t128) * rSges(3,2)) + g(2) * (t112 * pkin(1) + (t108 * t128 + t127) * rSges(3,1) + (t108 * t126 - t129) * rSges(3,2)) + g(3) * (t108 * rSges(3,3) + t136) + (g(3) * (rSges(3,1) * t111 + rSges(3,2) * t115) + (g(1) * t112 - g(2) * t116) * (rSges(3,3) + pkin(8))) * t107) - m(4) * (g(1) * (t81 * rSges(4,1) + t80 * rSges(4,2) + t95 + (rSges(4,3) * t107 - t90) * t112) + g(2) * (rSges(4,1) * t79 + rSges(4,2) * t78 - rSges(4,3) * t131 + t135) + g(3) * (-rSges(4,1) * t88 + rSges(4,2) * t87 + rSges(4,3) * t108 + t123)) - m(5) * (g(1) * (t75 * rSges(5,1) - t74 * rSges(5,2) - t139 * t80 + t122) + g(2) * (t73 * rSges(5,1) - t72 * rSges(5,2) - t139 * t78 + t125) + g(3) * (t83 * rSges(5,1) - t82 * rSges(5,2) - t139 * t87 + t121)) - m(6) * (g(1) * (t75 * pkin(4) - t80 * pkin(9) + (-t109 * t80 + t113 * t75) * rSges(6,1) + (-t109 * t75 - t113 * t80) * rSges(6,2) + t138 * t74 + t122) + g(2) * (t73 * pkin(4) - t78 * pkin(9) + (-t109 * t78 + t113 * t73) * rSges(6,1) + (-t109 * t73 - t113 * t78) * rSges(6,2) + t138 * t72 + t125) + g(3) * (t83 * pkin(4) - t87 * pkin(9) + (-t109 * t87 + t113 * t83) * rSges(6,1) + (-t109 * t83 - t113 * t87) * rSges(6,2) + t138 * t82 + t121)) - m(7) * (g(1) * (-t119 * t80 + t120 * t75 + t134 * t74 + t122) + g(2) * (-t119 * t78 + t120 * t73 + t134 * t72 + t125) + g(3) * (-t119 * t87 + t120 * t83 + t134 * t82 + t121));
U  = t1;
