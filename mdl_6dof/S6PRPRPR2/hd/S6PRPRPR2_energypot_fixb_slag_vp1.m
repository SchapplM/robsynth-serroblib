% Calculate potential energy for
% S6PRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRPR2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR2_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR2_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:29:26
% EndTime: 2019-03-08 19:29:26
% DurationCPUTime: 0.58s
% Computational Cost: add. (400->119), mult. (853->161), div. (0->0), fcn. (1059->14), ass. (0->51)
t104 = sin(pkin(11));
t112 = sin(qJ(2));
t114 = cos(qJ(2));
t129 = cos(pkin(11));
t89 = -t112 * t104 + t114 * t129;
t134 = rSges(5,3) + pkin(8);
t105 = sin(pkin(10));
t108 = cos(pkin(10));
t106 = sin(pkin(6));
t109 = cos(pkin(6));
t126 = t109 * t112;
t87 = pkin(2) * t126 + (-pkin(7) - qJ(3)) * t106;
t97 = pkin(2) * t114 + pkin(1);
t133 = t105 * t97 + t108 * t87;
t132 = rSges(4,3) * t106;
t131 = qJ(5) + rSges(6,3);
t130 = pkin(9) + qJ(5) + rSges(7,3);
t111 = sin(qJ(4));
t128 = t106 * t111;
t113 = cos(qJ(4));
t127 = t106 * t113;
t125 = t109 * t114;
t123 = t109 * pkin(7) + qJ(1);
t88 = -t114 * t104 - t112 * t129;
t86 = t88 * t109;
t76 = t105 * t89 - t108 * t86;
t122 = t76 * pkin(3) + t133;
t120 = t106 * t112 * pkin(2) + t109 * qJ(3) + t123;
t78 = t105 * t86 + t108 * t89;
t91 = t108 * t97;
t119 = t78 * pkin(3) - t105 * t87 + t91;
t85 = t88 * t106;
t118 = -t85 * pkin(3) + t120;
t107 = cos(pkin(12));
t102 = pkin(12) + qJ(6);
t98 = sin(t102);
t99 = cos(t102);
t117 = rSges(7,1) * t99 - rSges(7,2) * t98 + pkin(5) * t107 + pkin(4);
t103 = sin(pkin(12));
t116 = rSges(7,1) * t98 + rSges(7,2) * t99 + pkin(5) * t103 + pkin(8);
t115 = t109 * t89;
t84 = t89 * t106;
t80 = t109 * t111 - t113 * t85;
t79 = -t109 * t113 - t111 * t85;
t77 = -t105 * t115 + t108 * t88;
t75 = t105 * t88 + t108 * t115;
t72 = t105 * t128 + t113 * t78;
t71 = -t105 * t127 + t111 * t78;
t70 = -t108 * t128 + t113 * t76;
t69 = t108 * t127 + t111 * t76;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t108 - rSges(2,2) * t105) + g(2) * (rSges(2,1) * t105 + rSges(2,2) * t108) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t108 * pkin(1) + (-t105 * t126 + t108 * t114) * rSges(3,1) + (-t105 * t125 - t108 * t112) * rSges(3,2)) + g(2) * (t105 * pkin(1) + (t105 * t114 + t108 * t126) * rSges(3,1) + (-t105 * t112 + t108 * t125) * rSges(3,2)) + g(3) * (t109 * rSges(3,3) + t123) + (g(3) * (rSges(3,1) * t112 + rSges(3,2) * t114) + (g(1) * t105 - g(2) * t108) * (rSges(3,3) + pkin(7))) * t106) - m(4) * (g(1) * (t78 * rSges(4,1) + t77 * rSges(4,2) + t91 + (-t87 + t132) * t105) + g(2) * (rSges(4,1) * t76 + rSges(4,2) * t75 - t108 * t132 + t133) + g(3) * (-rSges(4,1) * t85 + rSges(4,2) * t84 + rSges(4,3) * t109 + t120)) - m(5) * (g(1) * (t72 * rSges(5,1) - t71 * rSges(5,2) - t134 * t77 + t119) + g(2) * (t70 * rSges(5,1) - t69 * rSges(5,2) - t134 * t75 + t122) + g(3) * (t80 * rSges(5,1) - t79 * rSges(5,2) - t134 * t84 + t118)) - m(6) * (g(1) * (t72 * pkin(4) - t77 * pkin(8) + (-t103 * t77 + t107 * t72) * rSges(6,1) + (-t103 * t72 - t107 * t77) * rSges(6,2) + t131 * t71 + t119) + g(2) * (t70 * pkin(4) - t75 * pkin(8) + (-t103 * t75 + t107 * t70) * rSges(6,1) + (-t103 * t70 - t107 * t75) * rSges(6,2) + t131 * t69 + t122) + g(3) * (t80 * pkin(4) - t84 * pkin(8) + (-t103 * t84 + t107 * t80) * rSges(6,1) + (-t103 * t80 - t107 * t84) * rSges(6,2) + t131 * t79 + t118)) - m(7) * (g(1) * (-t116 * t77 + t117 * t72 + t130 * t71 + t119) + g(2) * (-t116 * t75 + t117 * t70 + t130 * t69 + t122) + g(3) * (-t116 * t84 + t117 * t80 + t130 * t79 + t118));
U  = t1;
