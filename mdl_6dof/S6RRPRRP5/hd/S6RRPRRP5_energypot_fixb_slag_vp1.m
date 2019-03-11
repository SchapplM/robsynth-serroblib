% Calculate potential energy for
% S6RRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP5_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP5_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:57:41
% EndTime: 2019-03-09 11:57:41
% DurationCPUTime: 0.54s
% Computational Cost: add. (388->120), mult. (853->161), div. (0->0), fcn. (1059->12), ass. (0->56)
t111 = sin(pkin(11));
t117 = sin(qJ(2));
t121 = cos(qJ(2));
t137 = cos(pkin(11));
t100 = -t117 * t111 + t121 * t137;
t143 = rSges(5,3) + pkin(9);
t142 = rSges(6,3) + pkin(10);
t141 = pkin(2) * t117;
t113 = cos(pkin(6));
t140 = t113 * pkin(8) + pkin(7);
t139 = rSges(7,3) + qJ(6) + pkin(10);
t108 = pkin(2) * t121 + pkin(1);
t118 = sin(qJ(1));
t122 = cos(qJ(1));
t112 = sin(pkin(6));
t98 = t113 * t141 + (-pkin(8) - qJ(3)) * t112;
t138 = t118 * t108 + t122 * t98;
t135 = t117 * t122;
t134 = t118 * t112;
t133 = t118 * t117;
t132 = t118 * t121;
t131 = t121 * t122;
t130 = t122 * t112;
t99 = -t121 * t111 - t117 * t137;
t97 = t99 * t113;
t87 = t118 * t100 - t122 * t97;
t129 = t87 * pkin(3) + t138;
t115 = sin(qJ(5));
t128 = pkin(5) * t115 + pkin(9);
t126 = t113 * qJ(3) + t112 * t141 + t140;
t103 = t122 * t108;
t89 = t100 * t122 + t118 * t97;
t125 = t89 * pkin(3) - t118 * t98 + t103;
t96 = t99 * t112;
t124 = -t96 * pkin(3) + t126;
t123 = t113 * t100;
t120 = cos(qJ(4));
t119 = cos(qJ(5));
t116 = sin(qJ(4));
t107 = pkin(5) * t119 + pkin(4);
t95 = t100 * t112;
t91 = t113 * t116 - t120 * t96;
t90 = -t113 * t120 - t116 * t96;
t88 = -t118 * t123 + t122 * t99;
t86 = t118 * t99 + t122 * t123;
t83 = t116 * t134 + t120 * t89;
t82 = t116 * t89 - t120 * t134;
t81 = -t116 * t130 + t87 * t120;
t80 = t87 * t116 + t120 * t130;
t79 = -t115 * t95 + t119 * t91;
t78 = -t115 * t91 - t119 * t95;
t77 = -t115 * t88 + t119 * t83;
t76 = -t115 * t83 - t119 * t88;
t75 = -t115 * t86 + t119 * t81;
t74 = -t115 * t81 - t119 * t86;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t122 - t118 * rSges(2,2)) + g(2) * (t118 * rSges(2,1) + rSges(2,2) * t122) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (t122 * pkin(1) + (-t113 * t133 + t131) * rSges(3,1) + (-t113 * t132 - t135) * rSges(3,2)) + g(2) * (t118 * pkin(1) + (t113 * t135 + t132) * rSges(3,1) + (t113 * t131 - t133) * rSges(3,2)) + g(3) * (rSges(3,3) * t113 + t140) + (g(3) * (rSges(3,1) * t117 + rSges(3,2) * t121) + (g(1) * t118 - g(2) * t122) * (rSges(3,3) + pkin(8))) * t112) - m(4) * (g(1) * (rSges(4,1) * t89 + rSges(4,2) * t88 + t103 + (rSges(4,3) * t112 - t98) * t118) + g(2) * (t87 * rSges(4,1) + t86 * rSges(4,2) - rSges(4,3) * t130 + t138) + g(3) * (-rSges(4,1) * t96 + rSges(4,2) * t95 + rSges(4,3) * t113 + t126)) - m(5) * (g(1) * (rSges(5,1) * t83 - rSges(5,2) * t82 - t143 * t88 + t125) + g(2) * (rSges(5,1) * t81 - rSges(5,2) * t80 - t143 * t86 + t129) + g(3) * (rSges(5,1) * t91 - rSges(5,2) * t90 - t143 * t95 + t124)) - m(6) * (g(1) * (rSges(6,1) * t77 + rSges(6,2) * t76 + pkin(4) * t83 - pkin(9) * t88 + t142 * t82 + t125) + g(2) * (rSges(6,1) * t75 + rSges(6,2) * t74 + pkin(4) * t81 - pkin(9) * t86 + t142 * t80 + t129) + g(3) * (rSges(6,1) * t79 + rSges(6,2) * t78 + pkin(4) * t91 - pkin(9) * t95 + t142 * t90 + t124)) - m(7) * (g(1) * (rSges(7,1) * t77 + rSges(7,2) * t76 + t107 * t83 - t128 * t88 + t139 * t82 + t125) + g(2) * (rSges(7,1) * t75 + rSges(7,2) * t74 + t107 * t81 - t128 * t86 + t139 * t80 + t129) + g(3) * (rSges(7,1) * t79 + rSges(7,2) * t78 + t107 * t91 - t128 * t95 + t139 * t90 + t124));
U  = t1;
