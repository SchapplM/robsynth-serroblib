% Calculate potential energy for
% S6RRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPPR4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR4_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR4_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:31:40
% EndTime: 2019-03-09 15:31:40
% DurationCPUTime: 0.52s
% Computational Cost: add. (244->107), mult. (306->137), div. (0->0), fcn. (314->10), ass. (0->35)
t106 = rSges(4,3) + pkin(8);
t105 = pkin(9) + rSges(7,3);
t80 = sin(qJ(1));
t84 = cos(qJ(1));
t104 = g(1) * t84 + g(2) * t80;
t79 = sin(qJ(2));
t100 = rSges(3,2) * t79;
t78 = sin(qJ(3));
t99 = t80 * t78;
t83 = cos(qJ(2));
t98 = t80 * t83;
t97 = t84 * t78;
t96 = t84 * t83;
t93 = t84 * pkin(1) + t80 * pkin(7);
t82 = cos(qJ(3));
t69 = t82 * pkin(3) + pkin(2);
t76 = -qJ(4) - pkin(8);
t91 = t79 * t69 + t83 * t76 + pkin(6);
t73 = t80 * pkin(1);
t90 = -t84 * pkin(7) + t73;
t89 = pkin(3) * t99 + t69 * t96 + t93;
t75 = qJ(3) + pkin(10);
t70 = sin(t75);
t71 = cos(t75);
t88 = t91 + (pkin(4) * t71 + qJ(5) * t70) * t79;
t59 = t70 * t96 - t80 * t71;
t60 = t80 * t70 + t71 * t96;
t87 = t60 * pkin(4) + t59 * qJ(5) + t89;
t86 = -pkin(3) * t97 + t69 * t98 + t90;
t57 = t70 * t98 + t84 * t71;
t58 = -t84 * t70 + t71 * t98;
t85 = t58 * pkin(4) + t57 * qJ(5) + t86;
t81 = cos(qJ(6));
t77 = sin(qJ(6));
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t84 * rSges(2,1) - t80 * rSges(2,2)) + g(2) * (t80 * rSges(2,1) + t84 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t80 * rSges(3,3) + t93) + g(2) * (rSges(3,1) * t98 - t80 * t100 + t73) + g(3) * (t79 * rSges(3,1) + t83 * rSges(3,2) + pkin(6)) + (g(1) * (rSges(3,1) * t83 - t100) + g(2) * (-rSges(3,3) - pkin(7))) * t84) - m(4) * (g(1) * (pkin(2) * t96 + (t82 * t96 + t99) * rSges(4,1) + (-t78 * t96 + t80 * t82) * rSges(4,2) + t93) + g(2) * (pkin(2) * t98 + (t82 * t98 - t97) * rSges(4,1) + (-t78 * t98 - t84 * t82) * rSges(4,2) + t90) + g(3) * (-t106 * t83 + pkin(6)) + (g(3) * (rSges(4,1) * t82 - rSges(4,2) * t78 + pkin(2)) + t104 * t106) * t79) - m(5) * (g(1) * (t60 * rSges(5,1) - t59 * rSges(5,2) + t89) + g(2) * (t58 * rSges(5,1) - t57 * rSges(5,2) + t86) + g(3) * (-t83 * rSges(5,3) + t91) + (g(3) * (rSges(5,1) * t71 - rSges(5,2) * t70) + t104 * (rSges(5,3) - t76)) * t79) - m(6) * (g(1) * (t60 * rSges(6,1) + t59 * rSges(6,3) + t87) + g(2) * (t58 * rSges(6,1) + t57 * rSges(6,3) + t85) + g(3) * (-t83 * rSges(6,2) + t88) + (g(3) * (rSges(6,1) * t71 + rSges(6,3) * t70) + t104 * (rSges(6,2) - t76)) * t79) - m(7) * (g(1) * (t60 * pkin(5) + (t59 * t77 + t60 * t81) * rSges(7,1) + (t59 * t81 - t60 * t77) * rSges(7,2) + t87) + g(2) * (t58 * pkin(5) + (t57 * t77 + t58 * t81) * rSges(7,1) + (t57 * t81 - t58 * t77) * rSges(7,2) + t85) + g(3) * (t105 * t83 + t88) + (g(3) * (t71 * pkin(5) + (t70 * t77 + t71 * t81) * rSges(7,1) + (t70 * t81 - t71 * t77) * rSges(7,2)) + t104 * (-t76 - t105)) * t79);
U  = t1;
