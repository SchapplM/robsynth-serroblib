% Calculate potential energy for
% S6RPRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRR6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR6_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR6_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:12:09
% EndTime: 2019-03-09 07:12:10
% DurationCPUTime: 0.51s
% Computational Cost: add. (236->104), mult. (212->128), div. (0->0), fcn. (196->12), ass. (0->43)
t106 = rSges(5,3) + pkin(8);
t81 = -pkin(9) - pkin(8);
t105 = rSges(6,3) - t81;
t104 = rSges(7,3) + pkin(10) - t81;
t78 = sin(qJ(1));
t80 = cos(qJ(1));
t103 = g(1) * t80 + g(2) * t78;
t74 = sin(pkin(11));
t99 = t74 * pkin(2) + pkin(6);
t79 = cos(qJ(4));
t63 = t79 * pkin(4) + pkin(3);
t71 = pkin(11) + qJ(3);
t64 = sin(t71);
t98 = rSges(4,2) * t64;
t65 = cos(t71);
t97 = t65 * t78;
t96 = t65 * t80;
t73 = qJ(4) + qJ(5);
t66 = sin(t73);
t95 = t66 * t78;
t94 = t66 * t80;
t67 = cos(t73);
t93 = t67 * t78;
t92 = t67 * t80;
t77 = sin(qJ(4));
t91 = t77 * t78;
t90 = t77 * t80;
t89 = t78 * t79;
t88 = t79 * t80;
t75 = cos(pkin(11));
t59 = pkin(2) * t75 + pkin(1);
t76 = -pkin(7) - qJ(2);
t85 = t78 * t59 + t80 * t76;
t84 = rSges(3,3) + qJ(2);
t58 = t80 * t59;
t83 = -t78 * t76 + t58;
t82 = rSges(3,1) * t75 - rSges(3,2) * t74 + pkin(1);
t69 = qJ(6) + t73;
t62 = cos(t69);
t61 = sin(t69);
t56 = t77 * pkin(4) + pkin(5) * t66;
t55 = pkin(5) * t67 + t63;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t80 - rSges(2,2) * t78) + g(2) * (rSges(2,1) * t78 + rSges(2,2) * t80) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (rSges(3,1) * t74 + rSges(3,2) * t75 + pkin(6)) + (g(1) * t82 - g(2) * t84) * t80 + (g(1) * t84 + g(2) * t82) * t78) - m(4) * (g(1) * (rSges(4,1) * t96 - t80 * t98 + t58) + g(2) * (-rSges(4,3) * t80 + t85) + g(3) * (rSges(4,1) * t64 + rSges(4,2) * t65 + t99) + (g(1) * (rSges(4,3) - t76) + g(2) * (rSges(4,1) * t65 - t98)) * t78) - m(5) * (g(1) * (pkin(3) * t96 + (t65 * t88 + t91) * rSges(5,1) + (-t65 * t90 + t89) * rSges(5,2) + t83) + g(2) * (pkin(3) * t97 + (t65 * t89 - t90) * rSges(5,1) + (-t65 * t91 - t88) * rSges(5,2) + t85) + g(3) * (-t106 * t65 + t99) + (g(3) * (rSges(5,1) * t79 - rSges(5,2) * t77 + pkin(3)) + t103 * t106) * t64) - m(6) * (g(1) * (t63 * t96 + pkin(4) * t91 + (t65 * t92 + t95) * rSges(6,1) + (-t65 * t94 + t93) * rSges(6,2) + t83) + g(2) * (t63 * t97 - pkin(4) * t90 + (t65 * t93 - t94) * rSges(6,1) + (-t65 * t95 - t92) * rSges(6,2) + t85) + g(3) * (-t105 * t65 + t99) + (g(3) * (rSges(6,1) * t67 - rSges(6,2) * t66 + t63) + t103 * t105) * t64) - m(7) * (g(1) * (t55 * t96 + t78 * t56 + (t61 * t78 + t62 * t96) * rSges(7,1) + (-t61 * t96 + t62 * t78) * rSges(7,2) + t83) + g(2) * (t55 * t97 - t80 * t56 + (-t61 * t80 + t62 * t97) * rSges(7,1) + (-t61 * t97 - t62 * t80) * rSges(7,2) + t85) + g(3) * (-t104 * t65 + t99) + (g(3) * (rSges(7,1) * t62 - rSges(7,2) * t61 + t55) + t103 * t104) * t64);
U  = t1;
