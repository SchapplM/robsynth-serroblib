% Calculate potential energy for
% S6RRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRP2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP2_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP2_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:00:36
% EndTime: 2019-03-10 01:00:36
% DurationCPUTime: 0.39s
% Computational Cost: add. (239->87), mult. (199->105), div. (0->0), fcn. (183->10), ass. (0->41)
t105 = rSges(7,1) + pkin(5);
t104 = rSges(7,3) + qJ(6);
t86 = -pkin(8) - pkin(7);
t103 = rSges(3,3) + pkin(7);
t81 = sin(qJ(2));
t102 = t81 * pkin(2) + pkin(6);
t84 = cos(qJ(2));
t72 = t84 * pkin(2) + pkin(1);
t79 = qJ(2) + qJ(3);
t75 = qJ(4) + t79;
t70 = sin(t75);
t82 = sin(qJ(1));
t101 = t70 * t82;
t85 = cos(qJ(1));
t100 = t70 * t85;
t71 = cos(t75);
t99 = t71 * t85;
t80 = sin(qJ(5));
t98 = t80 * t82;
t97 = t80 * t85;
t83 = cos(qJ(5));
t96 = t82 * t83;
t95 = t83 * t85;
t94 = rSges(4,3) - t86;
t74 = cos(t79);
t61 = pkin(3) * t74 + t72;
t78 = -pkin(9) + t86;
t93 = t82 * t61 + t85 * t78;
t73 = sin(t79);
t92 = pkin(3) * t73 + t102;
t91 = t70 * pkin(4) + t92;
t90 = t82 * t71 * pkin(4) + pkin(10) * t101 + t93;
t60 = t85 * t61;
t89 = pkin(4) * t99 + pkin(10) * t100 - t82 * t78 + t60;
t88 = rSges(3,1) * t84 - rSges(3,2) * t81 + pkin(1);
t87 = rSges(4,1) * t74 - rSges(4,2) * t73 + t72;
t58 = t71 * t95 + t98;
t57 = t71 * t97 - t96;
t56 = t71 * t96 - t97;
t55 = t71 * t98 + t95;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t85 * rSges(2,1) - t82 * rSges(2,2)) + g(2) * (t82 * rSges(2,1) + t85 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (t81 * rSges(3,1) + t84 * rSges(3,2) + pkin(6)) + (g(1) * t88 - g(2) * t103) * t85 + (g(1) * t103 + g(2) * t88) * t82) - m(4) * (g(3) * (t73 * rSges(4,1) + t74 * rSges(4,2) + t102) + (g(1) * t87 - g(2) * t94) * t85 + (g(1) * t94 + g(2) * t87) * t82) - m(5) * (g(1) * (rSges(5,1) * t99 - rSges(5,2) * t100 + t60) + g(2) * (-t85 * rSges(5,3) + t93) + g(3) * (t70 * rSges(5,1) + t71 * rSges(5,2) + t92) + (g(1) * (rSges(5,3) - t78) + g(2) * (rSges(5,1) * t71 - rSges(5,2) * t70)) * t82) - m(6) * (g(1) * (t58 * rSges(6,1) - t57 * rSges(6,2) + rSges(6,3) * t100 + t89) + g(2) * (t56 * rSges(6,1) - t55 * rSges(6,2) + rSges(6,3) * t101 + t90) + g(3) * ((-rSges(6,3) - pkin(10)) * t71 + (rSges(6,1) * t83 - rSges(6,2) * t80) * t70 + t91)) - m(7) * (g(1) * (t104 * t57 + t105 * t58 + t89) + g(2) * (t104 * t55 + t105 * t56 + t90) + g(3) * (t91 + (-rSges(7,2) - pkin(10)) * t71) + (g(3) * (t104 * t80 + t105 * t83) + (g(1) * t85 + g(2) * t82) * rSges(7,2)) * t70);
U  = t1;
