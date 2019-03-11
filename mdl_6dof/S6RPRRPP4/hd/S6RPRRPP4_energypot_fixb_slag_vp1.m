% Calculate potential energy for
% S6RPRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-03-09 04:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPP4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP4_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP4_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:38:35
% EndTime: 2019-03-09 04:38:36
% DurationCPUTime: 0.44s
% Computational Cost: add. (239->95), mult. (227->115), div. (0->0), fcn. (215->10), ass. (0->43)
t114 = rSges(7,1) + pkin(5);
t113 = rSges(5,3) + pkin(8);
t112 = rSges(7,3) + qJ(6);
t86 = sin(qJ(1));
t88 = cos(qJ(1));
t111 = g(1) * t88 + g(2) * t86;
t81 = sin(pkin(9));
t107 = t81 * pkin(2) + pkin(6);
t79 = pkin(9) + qJ(3);
t74 = sin(t79);
t106 = rSges(4,2) * t74;
t76 = cos(t79);
t105 = t86 * t76;
t80 = qJ(4) + pkin(10);
t77 = cos(t80);
t104 = t86 * t77;
t85 = sin(qJ(4));
t103 = t86 * t85;
t87 = cos(qJ(4));
t102 = t86 * t87;
t101 = t88 * t76;
t100 = t88 * t77;
t99 = t88 * t85;
t98 = t88 * t87;
t82 = cos(pkin(9));
t71 = t82 * pkin(2) + pkin(1);
t84 = -pkin(7) - qJ(2);
t95 = t86 * t71 + t88 * t84;
t94 = rSges(3,3) + qJ(2);
t67 = t88 * t71;
t93 = -t86 * t84 + t67;
t73 = t87 * pkin(4) + pkin(3);
t83 = -qJ(5) - pkin(8);
t92 = t74 * t73 + t76 * t83 + t107;
t91 = pkin(4) * t103 + t73 * t101 + t93;
t90 = rSges(3,1) * t82 - rSges(3,2) * t81 + pkin(1);
t89 = -pkin(4) * t99 + t73 * t105 + t95;
t75 = sin(t80);
t62 = t76 * t100 + t86 * t75;
t61 = t75 * t101 - t104;
t60 = t76 * t104 - t88 * t75;
t59 = t75 * t105 + t100;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t88 * rSges(2,1) - t86 * rSges(2,2)) + g(2) * (t86 * rSges(2,1) + t88 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (t81 * rSges(3,1) + t82 * rSges(3,2) + pkin(6)) + (g(1) * t90 - g(2) * t94) * t88 + (g(1) * t94 + g(2) * t90) * t86) - m(4) * (g(1) * (rSges(4,1) * t101 - t88 * t106 + t67) + g(2) * (-t88 * rSges(4,3) + t95) + g(3) * (t74 * rSges(4,1) + t76 * rSges(4,2) + t107) + (g(1) * (rSges(4,3) - t84) + g(2) * (rSges(4,1) * t76 - t106)) * t86) - m(5) * (g(1) * (pkin(3) * t101 + (t76 * t98 + t103) * rSges(5,1) + (-t76 * t99 + t102) * rSges(5,2) + t93) + g(2) * (pkin(3) * t105 + (t76 * t102 - t99) * rSges(5,1) + (-t76 * t103 - t98) * rSges(5,2) + t95) + g(3) * (-t113 * t76 + t107) + (g(3) * (rSges(5,1) * t87 - rSges(5,2) * t85 + pkin(3)) + t111 * t113) * t74) - m(6) * (g(1) * (t62 * rSges(6,1) - t61 * rSges(6,2) + t91) + g(2) * (t60 * rSges(6,1) - t59 * rSges(6,2) + t89) + g(3) * (-t76 * rSges(6,3) + t92) + (g(3) * (rSges(6,1) * t77 - rSges(6,2) * t75) + t111 * (rSges(6,3) - t83)) * t74) - m(7) * (g(1) * (t112 * t61 + t114 * t62 + t91) + g(2) * (t112 * t59 + t114 * t60 + t89) + g(3) * (-t76 * rSges(7,2) + t92) + (g(3) * (t112 * t75 + t114 * t77) + t111 * (rSges(7,2) - t83)) * t74);
U  = t1;
