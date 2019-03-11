% Calculate potential energy for
% S6RRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2019-03-09 15:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPPR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR1_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR1_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:21:07
% EndTime: 2019-03-09 15:21:08
% DurationCPUTime: 0.44s
% Computational Cost: add. (235->93), mult. (186->113), div. (0->0), fcn. (166->12), ass. (0->44)
t109 = rSges(7,3) + pkin(9) + qJ(5);
t108 = rSges(6,3) + qJ(5);
t80 = sin(qJ(1));
t82 = cos(qJ(1));
t107 = g(1) * t82 + g(2) * t80;
t83 = -pkin(8) - pkin(7);
t104 = rSges(3,3) + pkin(7);
t79 = sin(qJ(2));
t103 = t79 * pkin(2) + pkin(6);
t81 = cos(qJ(2));
t65 = t81 * pkin(2) + pkin(1);
t75 = qJ(2) + qJ(3);
t68 = pkin(10) + t75;
t61 = sin(t68);
t102 = rSges(5,2) * t61;
t62 = cos(t68);
t101 = t62 * t80;
t100 = t62 * t82;
t74 = pkin(11) + qJ(6);
t66 = sin(t74);
t99 = t66 * t80;
t98 = t66 * t82;
t67 = cos(t74);
t97 = t67 * t80;
t96 = t67 * t82;
t76 = sin(pkin(11));
t95 = t76 * t80;
t94 = t76 * t82;
t77 = cos(pkin(11));
t93 = t77 * t80;
t92 = t77 * t82;
t91 = rSges(4,3) - t83;
t70 = cos(t75);
t59 = pkin(3) * t70 + t65;
t73 = -qJ(4) + t83;
t89 = t80 * t59 + t82 * t73;
t69 = sin(t75);
t87 = pkin(3) * t69 + t103;
t58 = t82 * t59;
t86 = -t80 * t73 + t58;
t85 = rSges(3,1) * t81 - rSges(3,2) * t79 + pkin(1);
t84 = rSges(4,1) * t70 - rSges(4,2) * t69 + t65;
t63 = pkin(5) * t77 + pkin(4);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t82 - rSges(2,2) * t80) + g(2) * (rSges(2,1) * t80 + rSges(2,2) * t82) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (rSges(3,1) * t79 + rSges(3,2) * t81 + pkin(6)) + (g(1) * t85 - g(2) * t104) * t82 + (g(1) * t104 + g(2) * t85) * t80) - m(4) * (g(3) * (rSges(4,1) * t69 + rSges(4,2) * t70 + t103) + (g(1) * t84 - g(2) * t91) * t82 + (g(1) * t91 + g(2) * t84) * t80) - m(5) * (g(1) * (rSges(5,1) * t100 - t82 * t102 + t58) + g(2) * (-rSges(5,3) * t82 + t89) + g(3) * (rSges(5,1) * t61 + rSges(5,2) * t62 + t87) + (g(1) * (rSges(5,3) - t73) + g(2) * (rSges(5,1) * t62 - t102)) * t80) - m(6) * (g(1) * (pkin(4) * t100 + (t62 * t92 + t95) * rSges(6,1) + (-t62 * t94 + t93) * rSges(6,2) + t86) + g(2) * (pkin(4) * t101 + (t62 * t93 - t94) * rSges(6,1) + (-t62 * t95 - t92) * rSges(6,2) + t89) + g(3) * (-t108 * t62 + t87) + (g(3) * (rSges(6,1) * t77 - rSges(6,2) * t76 + pkin(4)) + t107 * t108) * t61) - m(7) * (g(1) * (t63 * t100 + pkin(5) * t95 + (t62 * t96 + t99) * rSges(7,1) + (-t62 * t98 + t97) * rSges(7,2) + t86) + g(2) * (t63 * t101 - pkin(5) * t94 + (t62 * t97 - t98) * rSges(7,1) + (-t62 * t99 - t96) * rSges(7,2) + t89) + g(3) * (-t109 * t62 + t87) + (g(3) * (rSges(7,1) * t67 - rSges(7,2) * t66 + t63) + t107 * t109) * t61);
U  = t1;
