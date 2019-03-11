% Calculate potential energy for
% S6RRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRR6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR6_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR6_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:25:00
% EndTime: 2019-03-09 18:25:00
% DurationCPUTime: 0.60s
% Computational Cost: add. (245->115), mult. (240->142), div. (0->0), fcn. (228->12), ass. (0->43)
t105 = rSges(4,3) + pkin(8);
t74 = -qJ(4) - pkin(8);
t104 = rSges(5,3) - t74;
t72 = -pkin(9) + t74;
t103 = rSges(6,3) - t72;
t102 = rSges(7,3) + pkin(10) - t72;
t77 = sin(qJ(1));
t80 = cos(qJ(1));
t101 = g(1) * t80 + g(2) * t77;
t78 = cos(qJ(3));
t62 = t78 * pkin(3) + pkin(2);
t76 = sin(qJ(2));
t97 = rSges(3,2) * t76;
t75 = sin(qJ(3));
t96 = t77 * t75;
t79 = cos(qJ(2));
t95 = t77 * t79;
t94 = t79 * t80;
t73 = qJ(3) + pkin(11);
t65 = qJ(5) + t73;
t61 = qJ(6) + t65;
t55 = sin(t61);
t93 = t80 * t55;
t56 = cos(t61);
t92 = t80 * t56;
t59 = sin(t65);
t91 = t80 * t59;
t60 = cos(t65);
t90 = t80 * t60;
t63 = sin(t73);
t89 = t80 * t63;
t64 = cos(t73);
t88 = t80 * t64;
t87 = t80 * t75;
t86 = t80 * t78;
t54 = t75 * pkin(3) + pkin(4) * t63;
t82 = t80 * pkin(1) + t77 * pkin(7);
t53 = pkin(4) * t64 + t62;
t69 = t77 * pkin(1);
t81 = -t80 * pkin(7) + t69;
t52 = pkin(5) * t59 + t54;
t51 = pkin(5) * t60 + t53;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t80 * rSges(2,1) - t77 * rSges(2,2)) + g(2) * (t77 * rSges(2,1) + t80 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t77 * rSges(3,3) + t82) + g(2) * (rSges(3,1) * t95 - t77 * t97 + t69) + g(3) * (t76 * rSges(3,1) + t79 * rSges(3,2) + pkin(6)) + (g(1) * (rSges(3,1) * t79 - t97) + g(2) * (-rSges(3,3) - pkin(7))) * t80) - m(4) * (g(1) * (pkin(2) * t94 + (t79 * t86 + t96) * rSges(4,1) + (t77 * t78 - t79 * t87) * rSges(4,2) + t82) + g(2) * (pkin(2) * t95 + (t78 * t95 - t87) * rSges(4,1) + (-t75 * t95 - t86) * rSges(4,2) + t81) + g(3) * (-t105 * t79 + pkin(6)) + (g(3) * (rSges(4,1) * t78 - rSges(4,2) * t75 + pkin(2)) + t101 * t105) * t76) - m(5) * (g(1) * (t62 * t94 + pkin(3) * t96 + (t77 * t63 + t79 * t88) * rSges(5,1) + (t77 * t64 - t79 * t89) * rSges(5,2) + t82) + g(2) * (t62 * t95 - pkin(3) * t87 + (t64 * t95 - t89) * rSges(5,1) + (-t63 * t95 - t88) * rSges(5,2) + t81) + g(3) * (-t104 * t79 + pkin(6)) + (g(3) * (rSges(5,1) * t64 - rSges(5,2) * t63 + t62) + t101 * t104) * t76) - m(6) * (g(1) * (t53 * t94 + t77 * t54 + (t77 * t59 + t79 * t90) * rSges(6,1) + (t77 * t60 - t79 * t91) * rSges(6,2) + t82) + g(2) * (t53 * t95 - t80 * t54 + (t60 * t95 - t91) * rSges(6,1) + (-t59 * t95 - t90) * rSges(6,2) + t81) + g(3) * (-t103 * t79 + pkin(6)) + (g(3) * (rSges(6,1) * t60 - rSges(6,2) * t59 + t53) + t101 * t103) * t76) - m(7) * (g(1) * (t51 * t94 + t77 * t52 + (t77 * t55 + t79 * t92) * rSges(7,1) + (t77 * t56 - t79 * t93) * rSges(7,2) + t82) + g(2) * (t51 * t95 - t80 * t52 + (t56 * t95 - t93) * rSges(7,1) + (-t55 * t95 - t92) * rSges(7,2) + t81) + g(3) * (-t102 * t79 + pkin(6)) + (g(3) * (rSges(7,1) * t56 - rSges(7,2) * t55 + t51) + t101 * t102) * t76);
U  = t1;
