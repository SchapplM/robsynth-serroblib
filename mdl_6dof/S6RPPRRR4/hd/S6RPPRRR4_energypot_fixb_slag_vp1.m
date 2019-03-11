% Calculate potential energy for
% S6RPPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 02:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRR4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR4_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR4_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:25:18
% EndTime: 2019-03-09 02:25:18
% DurationCPUTime: 0.45s
% Computational Cost: add. (183->89), mult. (273->115), div. (0->0), fcn. (305->10), ass. (0->37)
t102 = rSges(6,3) + pkin(8);
t101 = rSges(7,3) + pkin(9) + pkin(8);
t81 = sin(pkin(10));
t82 = cos(pkin(10));
t93 = sin(qJ(1));
t94 = cos(qJ(1));
t54 = -t93 * t81 - t94 * t82;
t55 = t94 * t81 - t93 * t82;
t100 = g(1) * t54 + g(2) * t55;
t71 = cos(qJ(4));
t97 = t71 * pkin(4);
t96 = rSges(5,3) + pkin(7);
t68 = sin(qJ(5));
t92 = t54 * t68;
t91 = t55 * t68;
t67 = qJ(5) + qJ(6);
t60 = sin(t67);
t90 = t60 * t71;
t61 = cos(t67);
t89 = t61 * t71;
t88 = t68 * t71;
t70 = cos(qJ(5));
t87 = t70 * t71;
t59 = pkin(5) * t70 + pkin(4);
t86 = t71 * t59;
t84 = pkin(6) - qJ(3);
t83 = t94 * pkin(1) + t93 * qJ(2);
t80 = t94 * pkin(2) + t83;
t79 = -t54 * pkin(3) + t80;
t78 = t93 * pkin(1) - t94 * qJ(2);
t69 = sin(qJ(4));
t77 = -rSges(5,1) * t71 + rSges(5,2) * t69;
t76 = t93 * pkin(2) + t78;
t75 = t55 * pkin(7) + t79;
t74 = -t55 * pkin(3) + t76;
t73 = -t54 * pkin(7) + t74;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t94 * rSges(2,1) - t93 * rSges(2,2)) + g(2) * (t93 * rSges(2,1) + t94 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t94 * rSges(3,1) + t93 * rSges(3,3) + t83) + g(2) * (t93 * rSges(3,1) - t94 * rSges(3,3) + t78) + g(3) * (pkin(6) + rSges(3,2))) - m(4) * (g(1) * (-rSges(4,1) * t54 - rSges(4,2) * t55 + t80) + g(2) * (-t55 * rSges(4,1) + t54 * rSges(4,2) + t76) + g(3) * (-rSges(4,3) + t84)) - m(5) * (g(1) * t79 + g(2) * t74 + g(3) * (-rSges(5,1) * t69 - rSges(5,2) * t71 + t84) + (g(1) * t96 + g(2) * t77) * t55 + (g(1) * t77 - g(2) * t96) * t54) - m(6) * (g(1) * (-t54 * t97 + (-t54 * t87 + t91) * rSges(6,1) + (t54 * t88 + t55 * t70) * rSges(6,2) + t75) + g(2) * (-t55 * t97 + (-t55 * t87 - t92) * rSges(6,1) + (-t54 * t70 + t55 * t88) * rSges(6,2) + t73) + g(3) * (t102 * t71 + t84) + (g(3) * (-rSges(6,1) * t70 + rSges(6,2) * t68 - pkin(4)) - t100 * t102) * t69) - m(7) * (g(1) * (-t54 * t86 + pkin(5) * t91 + (-t54 * t89 + t55 * t60) * rSges(7,1) + (t54 * t90 + t55 * t61) * rSges(7,2) + t75) + g(2) * (-t55 * t86 - pkin(5) * t92 + (-t54 * t60 - t55 * t89) * rSges(7,1) + (-t54 * t61 + t55 * t90) * rSges(7,2) + t73) + g(3) * (t101 * t71 + t84) + (g(3) * (-rSges(7,1) * t61 + rSges(7,2) * t60 - t59) - t100 * t101) * t69);
U  = t1;
