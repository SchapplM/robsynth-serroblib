% Calculate potential energy for
% S6RPRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP5_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP5_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:10:40
% EndTime: 2019-03-09 06:10:41
% DurationCPUTime: 0.37s
% Computational Cost: add. (239->87), mult. (199->105), div. (0->0), fcn. (183->10), ass. (0->41)
t105 = rSges(7,1) + pkin(5);
t104 = rSges(7,3) + qJ(6);
t80 = sin(pkin(10));
t103 = t80 * pkin(2) + pkin(6);
t81 = cos(pkin(10));
t71 = t81 * pkin(2) + pkin(1);
t79 = pkin(10) + qJ(3);
t75 = qJ(4) + t79;
t69 = sin(t75);
t84 = sin(qJ(1));
t102 = t69 * t84;
t86 = cos(qJ(1));
t101 = t69 * t86;
t70 = cos(t75);
t100 = t70 * t86;
t83 = sin(qJ(5));
t99 = t83 * t86;
t98 = t84 * t83;
t85 = cos(qJ(5));
t97 = t84 * t85;
t96 = t85 * t86;
t82 = -pkin(7) - qJ(2);
t95 = rSges(4,3) - t82;
t74 = cos(t79);
t61 = pkin(3) * t74 + t71;
t78 = -pkin(8) + t82;
t94 = t84 * t61 + t86 * t78;
t93 = rSges(3,3) + qJ(2);
t73 = sin(t79);
t92 = pkin(3) * t73 + t103;
t91 = t69 * pkin(4) + t92;
t90 = t84 * t70 * pkin(4) + pkin(9) * t102 + t94;
t60 = t86 * t61;
t89 = pkin(4) * t100 + pkin(9) * t101 - t84 * t78 + t60;
t88 = rSges(3,1) * t81 - rSges(3,2) * t80 + pkin(1);
t87 = rSges(4,1) * t74 - rSges(4,2) * t73 + t71;
t58 = t70 * t96 + t98;
t57 = t70 * t99 - t97;
t56 = t70 * t97 - t99;
t55 = t70 * t98 + t96;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t86 - t84 * rSges(2,2)) + g(2) * (t84 * rSges(2,1) + rSges(2,2) * t86) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (rSges(3,1) * t80 + rSges(3,2) * t81 + pkin(6)) + (g(1) * t88 - g(2) * t93) * t86 + (g(1) * t93 + g(2) * t88) * t84) - m(4) * (g(3) * (rSges(4,1) * t73 + rSges(4,2) * t74 + t103) + (g(1) * t87 - g(2) * t95) * t86 + (g(1) * t95 + g(2) * t87) * t84) - m(5) * (g(1) * (rSges(5,1) * t100 - rSges(5,2) * t101 + t60) + g(2) * (-rSges(5,3) * t86 + t94) + g(3) * (rSges(5,1) * t69 + rSges(5,2) * t70 + t92) + (g(1) * (rSges(5,3) - t78) + g(2) * (rSges(5,1) * t70 - rSges(5,2) * t69)) * t84) - m(6) * (g(1) * (t58 * rSges(6,1) - t57 * rSges(6,2) + rSges(6,3) * t101 + t89) + g(2) * (rSges(6,1) * t56 - rSges(6,2) * t55 + rSges(6,3) * t102 + t90) + g(3) * ((-rSges(6,3) - pkin(9)) * t70 + (rSges(6,1) * t85 - rSges(6,2) * t83) * t69 + t91)) - m(7) * (g(1) * (t104 * t57 + t105 * t58 + t89) + g(2) * (t104 * t55 + t105 * t56 + t90) + g(3) * (t91 + (-rSges(7,2) - pkin(9)) * t70) + (g(3) * (t104 * t83 + t105 * t85) + (g(1) * t86 + g(2) * t84) * rSges(7,2)) * t69);
U  = t1;
