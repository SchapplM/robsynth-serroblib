% Calculate potential energy for
% S6PRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRR2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR2_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR2_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:26:34
% EndTime: 2019-03-08 20:26:35
% DurationCPUTime: 0.56s
% Computational Cost: add. (400->119), mult. (853->161), div. (0->0), fcn. (1059->14), ass. (0->50)
t103 = sin(pkin(12));
t110 = sin(qJ(2));
t113 = cos(qJ(2));
t129 = cos(pkin(12));
t89 = -t110 * t103 + t113 * t129;
t133 = rSges(5,3) + pkin(8);
t132 = pkin(9) + rSges(6,3);
t104 = sin(pkin(11));
t106 = cos(pkin(11));
t105 = sin(pkin(6));
t107 = cos(pkin(6));
t126 = t107 * t110;
t87 = pkin(2) * t126 + (-pkin(7) - qJ(3)) * t105;
t97 = t113 * pkin(2) + pkin(1);
t131 = t104 * t97 + t106 * t87;
t130 = pkin(10) + pkin(9) + rSges(7,3);
t128 = t104 * t105;
t127 = t106 * t105;
t125 = t107 * t113;
t123 = t107 * pkin(7) + qJ(1);
t88 = -t113 * t103 - t110 * t129;
t86 = t88 * t107;
t76 = t104 * t89 - t106 * t86;
t122 = t76 * pkin(3) + t131;
t120 = t105 * t110 * pkin(2) + t107 * qJ(3) + t123;
t78 = t104 * t86 + t106 * t89;
t91 = t106 * t97;
t119 = t78 * pkin(3) - t104 * t87 + t91;
t85 = t88 * t105;
t118 = -t85 * pkin(3) + t120;
t102 = qJ(5) + qJ(6);
t100 = cos(t102);
t111 = cos(qJ(5));
t99 = sin(t102);
t117 = t100 * rSges(7,1) - t99 * rSges(7,2) + t111 * pkin(5) + pkin(4);
t108 = sin(qJ(5));
t116 = t99 * rSges(7,1) + t100 * rSges(7,2) + t108 * pkin(5) + pkin(8);
t115 = t89 * t107;
t112 = cos(qJ(4));
t109 = sin(qJ(4));
t84 = t89 * t105;
t80 = t107 * t109 - t85 * t112;
t79 = -t107 * t112 - t85 * t109;
t77 = -t104 * t115 + t106 * t88;
t75 = t104 * t88 + t106 * t115;
t72 = t109 * t128 + t78 * t112;
t71 = t78 * t109 - t112 * t128;
t70 = -t109 * t127 + t76 * t112;
t69 = t76 * t109 + t112 * t127;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t106 * rSges(2,1) - t104 * rSges(2,2)) + g(2) * (t104 * rSges(2,1) + t106 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t106 * pkin(1) + (-t104 * t126 + t106 * t113) * rSges(3,1) + (-t104 * t125 - t106 * t110) * rSges(3,2)) + g(2) * (t104 * pkin(1) + (t104 * t113 + t106 * t126) * rSges(3,1) + (-t104 * t110 + t106 * t125) * rSges(3,2)) + g(3) * (t107 * rSges(3,3) + t123) + (g(3) * (rSges(3,1) * t110 + rSges(3,2) * t113) + (g(1) * t104 - g(2) * t106) * (rSges(3,3) + pkin(7))) * t105) - m(4) * (g(1) * (t78 * rSges(4,1) + t77 * rSges(4,2) + t91 + (rSges(4,3) * t105 - t87) * t104) + g(2) * (t76 * rSges(4,1) + t75 * rSges(4,2) - rSges(4,3) * t127 + t131) + g(3) * (-t85 * rSges(4,1) + t84 * rSges(4,2) + t107 * rSges(4,3) + t120)) - m(5) * (g(1) * (t72 * rSges(5,1) - t71 * rSges(5,2) - t133 * t77 + t119) + g(2) * (t70 * rSges(5,1) - t69 * rSges(5,2) - t133 * t75 + t122) + g(3) * (t80 * rSges(5,1) - t79 * rSges(5,2) - t133 * t84 + t118)) - m(6) * (g(1) * (t72 * pkin(4) - t77 * pkin(8) + (-t77 * t108 + t72 * t111) * rSges(6,1) + (-t72 * t108 - t77 * t111) * rSges(6,2) + t132 * t71 + t119) + g(2) * (t70 * pkin(4) - t75 * pkin(8) + (-t75 * t108 + t70 * t111) * rSges(6,1) + (-t70 * t108 - t75 * t111) * rSges(6,2) + t132 * t69 + t122) + g(3) * (t80 * pkin(4) - t84 * pkin(8) + (-t84 * t108 + t80 * t111) * rSges(6,1) + (-t80 * t108 - t84 * t111) * rSges(6,2) + t132 * t79 + t118)) - m(7) * (g(1) * (-t116 * t77 + t117 * t72 + t130 * t71 + t119) + g(2) * (-t116 * t75 + t117 * t70 + t130 * t69 + t122) + g(3) * (-t116 * t84 + t117 * t80 + t130 * t79 + t118));
U  = t1;
