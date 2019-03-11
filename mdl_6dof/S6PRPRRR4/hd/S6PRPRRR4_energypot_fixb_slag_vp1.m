% Calculate potential energy for
% S6PRPRRR4
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
% Datum: 2019-03-08 20:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRR4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR4_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR4_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:35:49
% EndTime: 2019-03-08 20:35:50
% DurationCPUTime: 0.45s
% Computational Cost: add. (355->130), mult. (545->172), div. (0->0), fcn. (633->14), ass. (0->50)
t130 = pkin(9) + rSges(6,3);
t100 = sin(pkin(12));
t129 = pkin(3) * t100;
t104 = cos(pkin(11));
t101 = sin(pkin(11));
t102 = sin(pkin(6));
t122 = t101 * t102;
t128 = t104 * pkin(1) + pkin(7) * t122;
t107 = sin(qJ(5));
t108 = sin(qJ(2));
t105 = cos(pkin(6));
t110 = cos(qJ(2));
t117 = t105 * t110;
t77 = t101 * t108 - t104 * t117;
t127 = t107 * t77;
t79 = t101 * t117 + t104 * t108;
t126 = t107 * t79;
t125 = t105 * pkin(7) + qJ(1);
t124 = qJ(3) + rSges(4,3);
t123 = pkin(10) + pkin(9) + rSges(7,3);
t121 = t102 * t104;
t120 = t102 * t108;
t119 = t102 * t110;
t118 = t105 * t108;
t116 = t100 * t122;
t115 = t107 * t119;
t106 = -pkin(8) - qJ(3);
t80 = -t101 * t118 + t104 * t110;
t103 = cos(pkin(12));
t89 = pkin(3) * t103 + pkin(2);
t114 = pkin(3) * t116 - t79 * t106 + t80 * t89 + t128;
t113 = t105 * t129 + t106 * t119 + t89 * t120 + t125;
t78 = t101 * t110 + t104 * t118;
t95 = t101 * pkin(1);
t112 = t78 * t89 + t95 + (-pkin(7) - t129) * t121 - t77 * t106;
t109 = cos(qJ(5));
t99 = qJ(5) + qJ(6);
t98 = pkin(12) + qJ(4);
t94 = cos(t99);
t93 = sin(t99);
t92 = cos(t98);
t91 = sin(t98);
t90 = pkin(5) * t109 + pkin(4);
t74 = t105 * t91 + t120 * t92;
t73 = -t105 * t92 + t120 * t91;
t70 = t122 * t91 + t80 * t92;
t69 = -t122 * t92 + t80 * t91;
t68 = -t121 * t91 + t78 * t92;
t67 = t121 * t92 + t78 * t91;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t104 - rSges(2,2) * t101) + g(2) * (rSges(2,1) * t101 + rSges(2,2) * t104) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t80 - rSges(3,2) * t79 + t128) + g(2) * (rSges(3,1) * t78 - rSges(3,2) * t77 + t95) + g(3) * (rSges(3,3) * t105 + t125) + (g(1) * rSges(3,3) * t101 + g(3) * (rSges(3,1) * t108 + rSges(3,2) * t110) + g(2) * (-rSges(3,3) - pkin(7)) * t104) * t102) - m(4) * (g(1) * (t80 * pkin(2) + (t103 * t80 + t116) * rSges(4,1) + (-t100 * t80 + t103 * t122) * rSges(4,2) + t124 * t79 + t128) + g(2) * (t78 * pkin(2) + t95 - pkin(7) * t121 + (-t100 * t121 + t103 * t78) * rSges(4,1) + (-t100 * t78 - t103 * t121) * rSges(4,2) + t124 * t77) + g(3) * ((t100 * rSges(4,1) + t103 * rSges(4,2)) * t105 + (-t124 * t110 + (t103 * rSges(4,1) - t100 * rSges(4,2) + pkin(2)) * t108) * t102 + t125)) - m(5) * (g(1) * (rSges(5,1) * t70 - rSges(5,2) * t69 + rSges(5,3) * t79 + t114) + g(2) * (rSges(5,1) * t68 - rSges(5,2) * t67 + rSges(5,3) * t77 + t112) + g(3) * (rSges(5,1) * t74 - rSges(5,2) * t73 - rSges(5,3) * t119 + t113)) - m(6) * (g(1) * (t70 * pkin(4) + (t109 * t70 + t126) * rSges(6,1) + (-t107 * t70 + t109 * t79) * rSges(6,2) + t130 * t69 + t114) + g(2) * (t68 * pkin(4) + (t109 * t68 + t127) * rSges(6,1) + (-t107 * t68 + t109 * t77) * rSges(6,2) + t130 * t67 + t112) + g(3) * (t74 * pkin(4) + (t109 * t74 - t115) * rSges(6,1) + (-t107 * t74 - t109 * t119) * rSges(6,2) + t130 * t73 + t113)) - m(7) * (g(1) * (t70 * t90 + pkin(5) * t126 + (t70 * t94 + t79 * t93) * rSges(7,1) + (-t70 * t93 + t79 * t94) * rSges(7,2) + t123 * t69 + t114) + g(2) * (t68 * t90 + pkin(5) * t127 + (t68 * t94 + t77 * t93) * rSges(7,1) + (-t68 * t93 + t77 * t94) * rSges(7,2) + t123 * t67 + t112) + g(3) * (t74 * t90 - pkin(5) * t115 + (-t119 * t93 + t74 * t94) * rSges(7,1) + (-t119 * t94 - t74 * t93) * rSges(7,2) + t123 * t73 + t113));
U  = t1;
