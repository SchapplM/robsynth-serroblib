% Calculate potential energy for
% S6RRPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 09:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR8_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR8_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR8_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR8_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:22:46
% EndTime: 2019-03-09 09:22:47
% DurationCPUTime: 0.56s
% Computational Cost: add. (206->114), mult. (362->146), div. (0->0), fcn. (388->10), ass. (0->42)
t83 = sin(qJ(1));
t86 = cos(qJ(1));
t113 = g(1) * t86 + g(2) * t83;
t85 = cos(qJ(2));
t110 = g(3) * t85;
t109 = -rSges(6,3) - pkin(8);
t82 = sin(qJ(2));
t108 = t82 * pkin(2) + pkin(6);
t80 = cos(pkin(10));
t101 = t86 * t80;
t103 = t83 * t85;
t79 = sin(pkin(10));
t60 = t79 * t103 + t101;
t81 = sin(qJ(5));
t107 = t60 * t81;
t102 = t86 * t79;
t62 = t85 * t102 - t83 * t80;
t106 = t62 * t81;
t105 = t79 * t81;
t104 = t83 * t82;
t100 = t86 * t82;
t99 = -rSges(7,3) - pkin(9) - pkin(8);
t98 = t86 * pkin(1) + t83 * pkin(7);
t97 = qJ(3) * t82;
t96 = rSges(5,3) + qJ(4);
t95 = t108 + (pkin(3) * t80 + qJ(4) * t79) * t82;
t94 = t98 + (pkin(2) * t85 + t97) * t86;
t63 = t85 * t101 + t83 * t79;
t93 = t63 * pkin(3) + t94;
t92 = g(3) * t95;
t76 = t83 * pkin(1);
t91 = pkin(2) * t103 - t86 * pkin(7) + t83 * t97 + t76;
t61 = t80 * t103 - t102;
t90 = t61 * pkin(3) + t91;
t89 = t62 * qJ(4) + t93;
t88 = t60 * qJ(4) + t90;
t84 = cos(qJ(5));
t78 = qJ(5) + qJ(6);
t73 = cos(t78);
t72 = sin(t78);
t71 = t84 * pkin(5) + pkin(4);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t86 * rSges(2,1) - t83 * rSges(2,2)) + g(2) * (t83 * rSges(2,1) + t86 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t83 * rSges(3,3) + t98) + g(2) * (rSges(3,1) * t103 - rSges(3,2) * t104 + t76) + g(3) * (t82 * rSges(3,1) + t85 * rSges(3,2) + pkin(6)) + (g(1) * (rSges(3,1) * t85 - rSges(3,2) * t82) + g(2) * (-rSges(3,3) - pkin(7))) * t86) - m(4) * (g(1) * (t63 * rSges(4,1) - t62 * rSges(4,2) + rSges(4,3) * t100 + t94) + g(2) * (t61 * rSges(4,1) - t60 * rSges(4,2) + rSges(4,3) * t104 + t91) + g(3) * ((-rSges(4,3) - qJ(3)) * t85 + (rSges(4,1) * t80 - rSges(4,2) * t79) * t82 + t108)) - m(5) * (g(1) * (t63 * rSges(5,1) + rSges(5,2) * t100 + t96 * t62 + t93) + g(2) * (t61 * rSges(5,1) + rSges(5,2) * t104 + t96 * t60 + t90) + g(3) * ((-rSges(5,2) - qJ(3)) * t85 + (rSges(5,1) * t80 + rSges(5,3) * t79) * t82 + t95)) - m(6) * (g(1) * (t63 * pkin(4) + (t63 * t84 + t106) * rSges(6,1) + (t62 * t84 - t63 * t81) * rSges(6,2) + t89) + g(2) * (t61 * pkin(4) + (t61 * t84 + t107) * rSges(6,1) + (t60 * t84 - t61 * t81) * rSges(6,2) + t88) + t92 + (-qJ(3) - t109) * t110 + (g(3) * (t80 * pkin(4) + (t80 * t84 + t105) * rSges(6,1) + (t79 * t84 - t80 * t81) * rSges(6,2)) + t113 * t109) * t82) - m(7) * (g(1) * (t63 * t71 + pkin(5) * t106 + (t62 * t72 + t63 * t73) * rSges(7,1) + (t62 * t73 - t63 * t72) * rSges(7,2) + t89) + g(2) * (t61 * t71 + pkin(5) * t107 + (t60 * t72 + t61 * t73) * rSges(7,1) + (t60 * t73 - t61 * t72) * rSges(7,2) + t88) + t92 + (-qJ(3) - t99) * t110 + (g(3) * (t80 * t71 + pkin(5) * t105 + (t72 * t79 + t73 * t80) * rSges(7,1) + (-t72 * t80 + t73 * t79) * rSges(7,2)) + t113 * t99) * t82);
U  = t1;
