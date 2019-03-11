% Calculate potential energy for
% S6PRPRPR4
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
% Datum: 2019-03-08 19:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRPR4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR4_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR4_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:38:56
% EndTime: 2019-03-08 19:38:56
% DurationCPUTime: 0.45s
% Computational Cost: add. (355->130), mult. (545->172), div. (0->0), fcn. (633->14), ass. (0->50)
t101 = sin(pkin(11));
t130 = pkin(3) * t101;
t106 = cos(pkin(10));
t102 = sin(pkin(10));
t103 = sin(pkin(6));
t122 = t102 * t103;
t129 = t106 * pkin(1) + pkin(7) * t122;
t100 = sin(pkin(12));
t110 = sin(qJ(2));
t107 = cos(pkin(6));
t111 = cos(qJ(2));
t117 = t107 * t111;
t77 = t102 * t110 - t106 * t117;
t128 = t100 * t77;
t79 = t102 * t117 + t106 * t110;
t127 = t100 * t79;
t126 = t107 * pkin(7) + qJ(1);
t125 = qJ(3) + rSges(4,3);
t124 = qJ(5) + rSges(6,3);
t123 = pkin(9) + qJ(5) + rSges(7,3);
t121 = t103 * t106;
t120 = t103 * t110;
t119 = t103 * t111;
t118 = t107 * t110;
t116 = t100 * t119;
t115 = t101 * t122;
t109 = -pkin(8) - qJ(3);
t80 = -t102 * t118 + t106 * t111;
t105 = cos(pkin(11));
t90 = pkin(3) * t105 + pkin(2);
t114 = pkin(3) * t115 - t79 * t109 + t80 * t90 + t129;
t113 = t107 * t130 + t109 * t119 + t90 * t120 + t126;
t78 = t102 * t111 + t106 * t118;
t95 = t102 * pkin(1);
t112 = t78 * t90 + t95 + (-pkin(7) - t130) * t121 - t77 * t109;
t104 = cos(pkin(12));
t99 = pkin(11) + qJ(4);
t98 = pkin(12) + qJ(6);
t94 = cos(t99);
t93 = cos(t98);
t92 = sin(t99);
t91 = sin(t98);
t89 = pkin(5) * t104 + pkin(4);
t74 = t107 * t92 + t120 * t94;
t73 = -t107 * t94 + t120 * t92;
t70 = t122 * t92 + t80 * t94;
t69 = -t122 * t94 + t80 * t92;
t68 = -t121 * t92 + t78 * t94;
t67 = t121 * t94 + t78 * t92;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t106 - rSges(2,2) * t102) + g(2) * (rSges(2,1) * t102 + rSges(2,2) * t106) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t80 - rSges(3,2) * t79 + t129) + g(2) * (rSges(3,1) * t78 - rSges(3,2) * t77 + t95) + g(3) * (t107 * rSges(3,3) + t126) + (g(1) * rSges(3,3) * t102 + g(3) * (rSges(3,1) * t110 + rSges(3,2) * t111) + g(2) * (-rSges(3,3) - pkin(7)) * t106) * t103) - m(4) * (g(1) * (t80 * pkin(2) + (t105 * t80 + t115) * rSges(4,1) + (-t101 * t80 + t105 * t122) * rSges(4,2) + t125 * t79 + t129) + g(2) * (t78 * pkin(2) + t95 - pkin(7) * t121 + (-t101 * t121 + t105 * t78) * rSges(4,1) + (-t101 * t78 - t105 * t121) * rSges(4,2) + t125 * t77) + g(3) * ((t101 * rSges(4,1) + t105 * rSges(4,2)) * t107 + (-t125 * t111 + (t105 * rSges(4,1) - t101 * rSges(4,2) + pkin(2)) * t110) * t103 + t126)) - m(5) * (g(1) * (rSges(5,1) * t70 - rSges(5,2) * t69 + rSges(5,3) * t79 + t114) + g(2) * (rSges(5,1) * t68 - rSges(5,2) * t67 + rSges(5,3) * t77 + t112) + g(3) * (t74 * rSges(5,1) - t73 * rSges(5,2) - rSges(5,3) * t119 + t113)) - m(6) * (g(1) * (t70 * pkin(4) + (t104 * t70 + t127) * rSges(6,1) + (-t100 * t70 + t104 * t79) * rSges(6,2) + t124 * t69 + t114) + g(2) * (t68 * pkin(4) + (t104 * t68 + t128) * rSges(6,1) + (-t100 * t68 + t104 * t77) * rSges(6,2) + t124 * t67 + t112) + g(3) * (t74 * pkin(4) + (t74 * t104 - t116) * rSges(6,1) + (-t74 * t100 - t104 * t119) * rSges(6,2) + t124 * t73 + t113)) - m(7) * (g(1) * (t70 * t89 + pkin(5) * t127 + (t70 * t93 + t79 * t91) * rSges(7,1) + (-t70 * t91 + t79 * t93) * rSges(7,2) + t123 * t69 + t114) + g(2) * (t68 * t89 + pkin(5) * t128 + (t68 * t93 + t77 * t91) * rSges(7,1) + (-t68 * t91 + t77 * t93) * rSges(7,2) + t123 * t67 + t112) + g(3) * (t74 * t89 - pkin(5) * t116 + (-t119 * t91 + t74 * t93) * rSges(7,1) + (-t119 * t93 - t74 * t91) * rSges(7,2) + t123 * t73 + t113));
U  = t1;
