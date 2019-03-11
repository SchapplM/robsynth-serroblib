% Calculate potential energy for
% S6RRRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-03-09 16:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPPR7_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR7_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR7_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR7_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:56:27
% EndTime: 2019-03-09 15:56:27
% DurationCPUTime: 0.57s
% Computational Cost: add. (206->114), mult. (362->146), div. (0->0), fcn. (388->10), ass. (0->40)
t83 = sin(qJ(1));
t86 = cos(qJ(1));
t110 = g(1) * t86 + g(2) * t83;
t85 = cos(qJ(2));
t107 = g(3) * t85;
t82 = sin(qJ(2));
t106 = t82 * pkin(2) + pkin(6);
t100 = t83 * t85;
t81 = sin(qJ(3));
t84 = cos(qJ(3));
t59 = t81 * t100 + t84 * t86;
t78 = sin(pkin(10));
t105 = t59 * t78;
t99 = t85 * t86;
t61 = t81 * t99 - t83 * t84;
t104 = t61 * t78;
t103 = t78 * t81;
t102 = t82 * t83;
t101 = t82 * t86;
t98 = -rSges(7,3) - pkin(9) - qJ(5);
t97 = t86 * pkin(1) + t83 * pkin(7);
t96 = rSges(5,3) + qJ(4);
t95 = -rSges(6,3) - qJ(5);
t94 = t106 + (pkin(3) * t84 + qJ(4) * t81) * t82;
t93 = pkin(2) * t99 + pkin(8) * t101 + t97;
t62 = t83 * t81 + t84 * t99;
t92 = t62 * pkin(3) + t93;
t91 = g(3) * t94;
t75 = t83 * pkin(1);
t90 = pkin(2) * t100 - t86 * pkin(7) + pkin(8) * t102 + t75;
t60 = t84 * t100 - t81 * t86;
t89 = t60 * pkin(3) + t90;
t88 = t61 * qJ(4) + t92;
t87 = t59 * qJ(4) + t89;
t79 = cos(pkin(10));
t77 = pkin(10) + qJ(6);
t72 = cos(t77);
t71 = sin(t77);
t69 = pkin(5) * t79 + pkin(4);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t86 - t83 * rSges(2,2)) + g(2) * (t83 * rSges(2,1) + rSges(2,2) * t86) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t83 * rSges(3,3) + t97) + g(2) * (rSges(3,1) * t100 - rSges(3,2) * t102 + t75) + g(3) * (rSges(3,1) * t82 + rSges(3,2) * t85 + pkin(6)) + (g(1) * (rSges(3,1) * t85 - rSges(3,2) * t82) + g(2) * (-rSges(3,3) - pkin(7))) * t86) - m(4) * (g(1) * (t62 * rSges(4,1) - t61 * rSges(4,2) + rSges(4,3) * t101 + t93) + g(2) * (t60 * rSges(4,1) - t59 * rSges(4,2) + rSges(4,3) * t102 + t90) + g(3) * ((-rSges(4,3) - pkin(8)) * t85 + (rSges(4,1) * t84 - rSges(4,2) * t81) * t82 + t106)) - m(5) * (g(1) * (t62 * rSges(5,1) + rSges(5,2) * t101 + t96 * t61 + t92) + g(2) * (t60 * rSges(5,1) + rSges(5,2) * t102 + t96 * t59 + t89) + g(3) * ((-rSges(5,2) - pkin(8)) * t85 + (rSges(5,1) * t84 + rSges(5,3) * t81) * t82 + t94)) - m(6) * (g(1) * (t62 * pkin(4) + (t62 * t79 + t104) * rSges(6,1) + (t61 * t79 - t62 * t78) * rSges(6,2) + t88) + g(2) * (t60 * pkin(4) + (t60 * t79 + t105) * rSges(6,1) + (t59 * t79 - t60 * t78) * rSges(6,2) + t87) + t91 + (-pkin(8) - t95) * t107 + (g(3) * (t84 * pkin(4) + (t79 * t84 + t103) * rSges(6,1) + (-t78 * t84 + t79 * t81) * rSges(6,2)) + t110 * t95) * t82) - m(7) * (g(1) * (t62 * t69 + pkin(5) * t104 + (t61 * t71 + t62 * t72) * rSges(7,1) + (t61 * t72 - t62 * t71) * rSges(7,2) + t88) + g(2) * (t60 * t69 + pkin(5) * t105 + (t59 * t71 + t60 * t72) * rSges(7,1) + (t59 * t72 - t60 * t71) * rSges(7,2) + t87) + t91 + (-pkin(8) - t98) * t107 + (g(3) * (t84 * t69 + pkin(5) * t103 + (t71 * t81 + t72 * t84) * rSges(7,1) + (-t71 * t84 + t72 * t81) * rSges(7,2)) + t110 * t98) * t82);
U  = t1;
