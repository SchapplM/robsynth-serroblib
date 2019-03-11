% Calculate potential energy for
% S6RRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRR5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR5_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR5_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:20:30
% EndTime: 2019-03-09 18:20:30
% DurationCPUTime: 0.46s
% Computational Cost: add. (206->98), mult. (204->120), div. (0->0), fcn. (184->10), ass. (0->40)
t107 = rSges(6,3) + pkin(9);
t106 = rSges(7,3) + pkin(10) + pkin(9);
t75 = sin(qJ(1));
t78 = cos(qJ(1));
t105 = g(1) * t78 + g(2) * t75;
t102 = rSges(3,3) + pkin(7);
t74 = sin(qJ(2));
t100 = t74 * pkin(2) + pkin(6);
t71 = qJ(5) + qJ(6);
t66 = sin(t71);
t99 = t66 * t75;
t72 = qJ(2) + qJ(3);
t67 = sin(t72);
t98 = t67 * t78;
t68 = cos(t71);
t97 = t68 * t75;
t96 = t68 * t78;
t69 = cos(t72);
t95 = t69 * t78;
t73 = sin(qJ(5));
t94 = t73 * t75;
t93 = t73 * t78;
t76 = cos(qJ(5));
t92 = t75 * t76;
t91 = t76 * t78;
t77 = cos(qJ(2));
t64 = pkin(2) * t77 + pkin(1);
t80 = -pkin(8) - pkin(7);
t89 = t75 * t64 + t78 * t80;
t88 = qJ(4) * t67;
t87 = t67 * pkin(3) + t100;
t86 = t67 * t94;
t85 = t67 * t93;
t59 = t78 * t64;
t84 = pkin(3) * t95 + t78 * t88 + t59;
t83 = t89 + (pkin(3) * t69 + t88) * t75;
t82 = -t75 * t80 + t84;
t81 = rSges(3,1) * t77 - rSges(3,2) * t74 + pkin(1);
t63 = pkin(5) * t76 + pkin(4);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t78 - rSges(2,2) * t75) + g(2) * (rSges(2,1) * t75 + rSges(2,2) * t78) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (rSges(3,1) * t74 + rSges(3,2) * t77 + pkin(6)) + (g(1) * t81 - g(2) * t102) * t78 + (g(1) * t102 + g(2) * t81) * t75) - m(4) * (g(1) * (rSges(4,1) * t95 - rSges(4,2) * t98 + t59) + g(2) * (-rSges(4,3) * t78 + t89) + g(3) * (rSges(4,1) * t67 + rSges(4,2) * t69 + t100) + (g(1) * (rSges(4,3) - t80) + g(2) * (rSges(4,1) * t69 - rSges(4,2) * t67)) * t75) - m(5) * (g(1) * (-rSges(5,2) * t95 + rSges(5,3) * t98 + t84) + g(2) * (-rSges(5,1) * t78 + t83) + g(3) * (-rSges(5,2) * t67 + (-rSges(5,3) - qJ(4)) * t69 + t87) + (g(1) * (rSges(5,1) - t80) + g(2) * (-rSges(5,2) * t69 + rSges(5,3) * t67)) * t75) - m(6) * (g(1) * (t75 * pkin(4) + (t85 + t92) * rSges(6,1) + (t67 * t91 - t94) * rSges(6,2) + t82) + g(2) * (-t78 * pkin(4) + (t86 - t91) * rSges(6,1) + (t67 * t92 + t93) * rSges(6,2) + t83) + g(3) * (t107 * t67 + t87) + (g(3) * (-rSges(6,1) * t73 - rSges(6,2) * t76 - qJ(4)) + t105 * t107) * t69) - m(7) * (g(1) * (t75 * t63 + pkin(5) * t85 + (t66 * t98 + t97) * rSges(7,1) + (t67 * t96 - t99) * rSges(7,2) + t82) + g(2) * (-t78 * t63 + pkin(5) * t86 + (t67 * t99 - t96) * rSges(7,1) + (t66 * t78 + t67 * t97) * rSges(7,2) + t83) + g(3) * (t106 * t67 + t87) + (g(3) * (-rSges(7,1) * t66 - rSges(7,2) * t68 - pkin(5) * t73 - qJ(4)) + t105 * t106) * t69);
U  = t1;
