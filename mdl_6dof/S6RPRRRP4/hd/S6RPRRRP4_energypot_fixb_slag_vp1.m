% Calculate potential energy for
% S6RPRRRP4
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
% Datum: 2019-03-09 06:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP4_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP4_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:06:53
% EndTime: 2019-03-09 06:06:53
% DurationCPUTime: 0.38s
% Computational Cost: add. (225->88), mult. (186->105), div. (0->0), fcn. (166->10), ass. (0->41)
t103 = rSges(6,3) + pkin(9);
t102 = rSges(7,3) + qJ(6) + pkin(9);
t79 = sin(qJ(1));
t81 = cos(qJ(1));
t101 = g(1) * t81 + g(2) * t79;
t74 = sin(pkin(10));
t97 = t74 * pkin(2) + pkin(6);
t75 = cos(pkin(10));
t65 = t75 * pkin(2) + pkin(1);
t73 = pkin(10) + qJ(3);
t69 = qJ(4) + t73;
t63 = sin(t69);
t96 = rSges(5,2) * t63;
t64 = cos(t69);
t95 = t64 * t79;
t94 = t64 * t81;
t78 = sin(qJ(5));
t93 = t78 * t81;
t92 = t79 * t78;
t80 = cos(qJ(5));
t91 = t79 * t80;
t90 = t81 * t80;
t77 = -pkin(7) - qJ(2);
t89 = rSges(4,3) - t77;
t68 = cos(t73);
t60 = pkin(3) * t68 + t65;
t72 = -pkin(8) + t77;
t87 = t79 * t60 + t81 * t72;
t86 = rSges(3,3) + qJ(2);
t67 = sin(t73);
t85 = pkin(3) * t67 + t97;
t59 = t81 * t60;
t84 = -t79 * t72 + t59;
t83 = rSges(3,1) * t75 - rSges(3,2) * t74 + pkin(1);
t82 = rSges(4,1) * t68 - rSges(4,2) * t67 + t65;
t66 = pkin(5) * t80 + pkin(4);
t57 = t64 * t90 + t92;
t56 = -t64 * t93 + t91;
t55 = t64 * t91 - t93;
t54 = -t64 * t92 - t90;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t81 - t79 * rSges(2,2)) + g(2) * (t79 * rSges(2,1) + rSges(2,2) * t81) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (t74 * rSges(3,1) + t75 * rSges(3,2) + pkin(6)) + (g(1) * t83 - g(2) * t86) * t81 + (g(1) * t86 + g(2) * t83) * t79) - m(4) * (g(3) * (rSges(4,1) * t67 + rSges(4,2) * t68 + t97) + (g(1) * t82 - g(2) * t89) * t81 + (g(1) * t89 + g(2) * t82) * t79) - m(5) * (g(1) * (rSges(5,1) * t94 - t81 * t96 + t59) + g(2) * (-rSges(5,3) * t81 + t87) + g(3) * (rSges(5,1) * t63 + rSges(5,2) * t64 + t85) + (g(1) * (rSges(5,3) - t72) + g(2) * (rSges(5,1) * t64 - t96)) * t79) - m(6) * (g(1) * (t57 * rSges(6,1) + t56 * rSges(6,2) + pkin(4) * t94 + t84) + g(2) * (t55 * rSges(6,1) + t54 * rSges(6,2) + pkin(4) * t95 + t87) + g(3) * (-t103 * t64 + t85) + (g(3) * (rSges(6,1) * t80 - rSges(6,2) * t78 + pkin(4)) + t101 * t103) * t63) - m(7) * (g(1) * (t57 * rSges(7,1) + t56 * rSges(7,2) + pkin(5) * t92 + t66 * t94 + t84) + g(2) * (t55 * rSges(7,1) + t54 * rSges(7,2) - pkin(5) * t93 + t66 * t95 + t87) + g(3) * (-t102 * t64 + t85) + (g(3) * (rSges(7,1) * t80 - rSges(7,2) * t78 + t66) + t101 * t102) * t63);
U  = t1;
