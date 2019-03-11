% Calculate potential energy for
% S6RRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRR4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR4_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR4_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR4_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR4_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:46:37
% EndTime: 2019-03-10 03:46:37
% DurationCPUTime: 0.60s
% Computational Cost: add. (245->115), mult. (240->142), div. (0->0), fcn. (228->12), ass. (0->43)
t105 = rSges(4,3) + pkin(8);
t80 = -pkin(9) - pkin(8);
t104 = rSges(5,3) - t80;
t72 = -pkin(10) + t80;
t103 = rSges(6,3) - t72;
t102 = rSges(7,3) + pkin(11) - t72;
t76 = sin(qJ(1));
t79 = cos(qJ(1));
t101 = g(1) * t79 + g(2) * t76;
t77 = cos(qJ(3));
t61 = t77 * pkin(3) + pkin(2);
t75 = sin(qJ(2));
t97 = rSges(3,2) * t75;
t74 = sin(qJ(3));
t96 = t76 * t74;
t78 = cos(qJ(2));
t95 = t76 * t78;
t94 = t78 * t79;
t73 = qJ(3) + qJ(4);
t65 = qJ(5) + t73;
t62 = qJ(6) + t65;
t55 = sin(t62);
t93 = t79 * t55;
t56 = cos(t62);
t92 = t79 * t56;
t59 = sin(t65);
t91 = t79 * t59;
t60 = cos(t65);
t90 = t79 * t60;
t63 = sin(t73);
t89 = t79 * t63;
t64 = cos(t73);
t88 = t79 * t64;
t87 = t79 * t74;
t86 = t79 * t77;
t54 = t74 * pkin(3) + pkin(4) * t63;
t82 = t79 * pkin(1) + t76 * pkin(7);
t53 = pkin(4) * t64 + t61;
t69 = t76 * pkin(1);
t81 = -t79 * pkin(7) + t69;
t52 = pkin(5) * t59 + t54;
t51 = pkin(5) * t60 + t53;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t79 * rSges(2,1) - t76 * rSges(2,2)) + g(2) * (t76 * rSges(2,1) + t79 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t76 * rSges(3,3) + t82) + g(2) * (rSges(3,1) * t95 - t76 * t97 + t69) + g(3) * (t75 * rSges(3,1) + t78 * rSges(3,2) + pkin(6)) + (g(1) * (rSges(3,1) * t78 - t97) + g(2) * (-rSges(3,3) - pkin(7))) * t79) - m(4) * (g(1) * (pkin(2) * t94 + (t78 * t86 + t96) * rSges(4,1) + (t76 * t77 - t78 * t87) * rSges(4,2) + t82) + g(2) * (pkin(2) * t95 + (t77 * t95 - t87) * rSges(4,1) + (-t74 * t95 - t86) * rSges(4,2) + t81) + g(3) * (-t105 * t78 + pkin(6)) + (g(3) * (rSges(4,1) * t77 - rSges(4,2) * t74 + pkin(2)) + t101 * t105) * t75) - m(5) * (g(1) * (t61 * t94 + pkin(3) * t96 + (t76 * t63 + t78 * t88) * rSges(5,1) + (t76 * t64 - t78 * t89) * rSges(5,2) + t82) + g(2) * (t61 * t95 - pkin(3) * t87 + (t64 * t95 - t89) * rSges(5,1) + (-t63 * t95 - t88) * rSges(5,2) + t81) + g(3) * (-t104 * t78 + pkin(6)) + (g(3) * (rSges(5,1) * t64 - rSges(5,2) * t63 + t61) + t101 * t104) * t75) - m(6) * (g(1) * (t53 * t94 + t76 * t54 + (t76 * t59 + t78 * t90) * rSges(6,1) + (t76 * t60 - t78 * t91) * rSges(6,2) + t82) + g(2) * (t53 * t95 - t79 * t54 + (t60 * t95 - t91) * rSges(6,1) + (-t59 * t95 - t90) * rSges(6,2) + t81) + g(3) * (-t103 * t78 + pkin(6)) + (g(3) * (rSges(6,1) * t60 - rSges(6,2) * t59 + t53) + t101 * t103) * t75) - m(7) * (g(1) * (t51 * t94 + t76 * t52 + (t76 * t55 + t78 * t92) * rSges(7,1) + (t76 * t56 - t78 * t93) * rSges(7,2) + t82) + g(2) * (t51 * t95 - t79 * t52 + (t56 * t95 - t93) * rSges(7,1) + (-t55 * t95 - t92) * rSges(7,2) + t81) + g(3) * (-t102 * t78 + pkin(6)) + (g(3) * (rSges(7,1) * t56 - rSges(7,2) * t55 + t51) + t101 * t102) * t75);
U  = t1;
