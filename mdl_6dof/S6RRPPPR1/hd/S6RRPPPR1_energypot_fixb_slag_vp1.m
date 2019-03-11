% Calculate potential energy for
% S6RRPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
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
% Datum: 2019-03-09 08:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPPR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR1_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPPR1_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:05:59
% EndTime: 2019-03-09 08:05:59
% DurationCPUTime: 0.49s
% Computational Cost: add. (239->101), mult. (274->127), div. (0->0), fcn. (278->10), ass. (0->39)
t104 = rSges(3,3) + pkin(7);
t103 = -rSges(7,3) - pkin(8);
t80 = sin(qJ(2));
t102 = t80 * pkin(2) + pkin(6);
t75 = qJ(2) + pkin(9);
t73 = cos(t75);
t84 = cos(qJ(1));
t101 = t73 * t84;
t72 = sin(t75);
t81 = sin(qJ(1));
t100 = t81 * t72;
t76 = sin(pkin(10));
t99 = t81 * t76;
t77 = cos(pkin(10));
t98 = t81 * t77;
t97 = t84 * t72;
t96 = t84 * t76;
t95 = t84 * t77;
t83 = cos(qJ(2));
t71 = t83 * pkin(2) + pkin(1);
t78 = -qJ(3) - pkin(7);
t94 = t81 * t71 + t84 * t78;
t93 = qJ(4) * t72;
t92 = rSges(6,3) + qJ(5);
t91 = t72 * pkin(3) + t102;
t90 = t94 + (pkin(3) * t73 + t93) * t81;
t89 = t91 + (pkin(4) * t77 + qJ(5) * t76) * t72;
t57 = t73 * t98 - t96;
t88 = t57 * pkin(4) + t90;
t67 = t84 * t71;
t87 = pkin(3) * t101 - t81 * t78 + t84 * t93 + t67;
t86 = rSges(3,1) * t83 - rSges(3,2) * t80 + pkin(1);
t59 = t73 * t95 + t99;
t85 = t59 * pkin(4) + t87;
t82 = cos(qJ(6));
t79 = sin(qJ(6));
t58 = t73 * t96 - t98;
t56 = t73 * t99 + t95;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t84 * rSges(2,1) - t81 * rSges(2,2)) + g(2) * (t81 * rSges(2,1) + t84 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (t80 * rSges(3,1) + t83 * rSges(3,2) + pkin(6)) + (g(1) * t86 - g(2) * t104) * t84 + (g(1) * t104 + g(2) * t86) * t81) - m(4) * (g(1) * (rSges(4,1) * t101 - rSges(4,2) * t97 + t67) + g(2) * (-t84 * rSges(4,3) + t94) + g(3) * (t72 * rSges(4,1) + t73 * rSges(4,2) + t102) + (g(1) * (rSges(4,3) - t78) + g(2) * (rSges(4,1) * t73 - rSges(4,2) * t72)) * t81) - m(5) * (g(1) * (t59 * rSges(5,1) - t58 * rSges(5,2) + rSges(5,3) * t97 + t87) + g(2) * (t57 * rSges(5,1) - t56 * rSges(5,2) + rSges(5,3) * t100 + t90) + g(3) * ((-rSges(5,3) - qJ(4)) * t73 + (rSges(5,1) * t77 - rSges(5,2) * t76) * t72 + t91)) - m(6) * (g(1) * (t59 * rSges(6,1) + rSges(6,2) * t97 + t92 * t58 + t85) + g(2) * (t57 * rSges(6,1) + rSges(6,2) * t100 + t92 * t56 + t88) + g(3) * ((-rSges(6,2) - qJ(4)) * t73 + (rSges(6,1) * t77 + rSges(6,3) * t76) * t72 + t89)) - m(7) * (g(1) * (t59 * pkin(5) + t58 * qJ(5) + (t58 * t79 + t59 * t82) * rSges(7,1) + (t58 * t82 - t59 * t79) * rSges(7,2) + t85) + g(2) * (t57 * pkin(5) + t56 * qJ(5) + (t56 * t79 + t57 * t82) * rSges(7,1) + (t56 * t82 - t57 * t79) * rSges(7,2) + t88) + (g(1) * t84 + g(2) * t81) * t72 * t103 + (t89 + (-qJ(4) - t103) * t73 + (t77 * pkin(5) + (t76 * t79 + t77 * t82) * rSges(7,1) + (t76 * t82 - t77 * t79) * rSges(7,2)) * t72) * g(3));
U  = t1;
