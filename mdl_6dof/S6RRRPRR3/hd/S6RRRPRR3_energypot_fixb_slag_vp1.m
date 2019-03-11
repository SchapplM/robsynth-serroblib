% Calculate potential energy for
% S6RRRPRR3
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
% Datum: 2019-03-09 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:10:47
% EndTime: 2019-03-09 18:10:48
% DurationCPUTime: 0.44s
% Computational Cost: add. (231->94), mult. (247->114), div. (0->0), fcn. (244->10), ass. (0->39)
t103 = pkin(10) + rSges(7,3);
t106 = rSges(3,3) + pkin(7);
t82 = sin(qJ(2));
t105 = t82 * pkin(2) + pkin(6);
t88 = -pkin(8) - pkin(7);
t104 = -pkin(9) - t88;
t79 = qJ(2) + qJ(3);
t75 = sin(t79);
t85 = cos(qJ(5));
t102 = t75 * t85;
t87 = cos(qJ(1));
t101 = t75 * t87;
t76 = cos(t79);
t83 = sin(qJ(1));
t100 = t76 * t83;
t99 = t76 * t87;
t86 = cos(qJ(2));
t73 = pkin(2) * t86 + pkin(1);
t98 = t83 * t73 + t87 * t88;
t97 = qJ(4) * t75;
t96 = t75 * pkin(3) + t105;
t66 = t87 * t73;
t95 = pkin(3) * t99 + t87 * t97 + t66;
t94 = pkin(3) * t100 + t83 * t97 + t98;
t93 = pkin(4) * t99 + t95;
t81 = sin(qJ(5));
t59 = t75 * t81 + t76 * t85;
t92 = pkin(4) * t100 + t87 * pkin(9) + t94;
t91 = rSges(3,1) * t86 - rSges(3,2) * t82 + pkin(1);
t80 = sin(qJ(6));
t84 = cos(qJ(6));
t90 = rSges(7,1) * t84 - rSges(7,2) * t80 + pkin(5);
t89 = t75 * pkin(4) - qJ(4) * t76 + t96;
t60 = -t76 * t81 + t102;
t58 = t59 * t87;
t57 = -t85 * t101 + t81 * t99;
t56 = t59 * t83;
t55 = t81 * t100 - t83 * t102;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t87 - rSges(2,2) * t83) + g(2) * (rSges(2,1) * t83 + rSges(2,2) * t87) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (rSges(3,1) * t82 + t86 * rSges(3,2) + pkin(6)) + (g(1) * t91 - g(2) * t106) * t87 + (g(1) * t106 + g(2) * t91) * t83) - m(4) * (g(1) * (rSges(4,1) * t99 - rSges(4,2) * t101 + t66) + g(2) * (-rSges(4,3) * t87 + t98) + g(3) * (rSges(4,1) * t75 + rSges(4,2) * t76 + t105) + (g(1) * (rSges(4,3) - t88) + g(2) * (rSges(4,1) * t76 - rSges(4,2) * t75)) * t83) - m(5) * (g(1) * (rSges(5,1) * t99 + rSges(5,3) * t101 + t95) + g(2) * (-rSges(5,2) * t87 + t94) + g(3) * (rSges(5,1) * t75 + (-rSges(5,3) - qJ(4)) * t76 + t96) + (g(1) * (rSges(5,2) - t88) + g(2) * (rSges(5,1) * t76 + rSges(5,3) * t75)) * t83) - m(6) * (g(2) * (rSges(6,1) * t56 - rSges(6,2) * t55 + rSges(6,3) * t87 + t92) + g(3) * (rSges(6,1) * t60 - rSges(6,2) * t59 + t89) + (t58 * rSges(6,1) - t57 * rSges(6,2) + t93 + (-rSges(6,3) + t104) * t83) * g(1)) - m(7) * (g(1) * (t90 * t58 + (-rSges(7,1) * t80 - rSges(7,2) * t84 + t104) * t83 + t103 * t57 + t93) + g(2) * (t56 * pkin(5) + (t56 * t84 + t80 * t87) * rSges(7,1) + (-t56 * t80 + t84 * t87) * rSges(7,2) + t103 * t55 + t92) + (t103 * t59 + t90 * t60 + t89) * g(3));
U  = t1;
