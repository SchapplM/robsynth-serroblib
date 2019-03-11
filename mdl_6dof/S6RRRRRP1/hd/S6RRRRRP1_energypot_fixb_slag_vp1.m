% Calculate potential energy for
% S6RRRRRP1
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
% Datum: 2019-03-10 00:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRP1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP1_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP1_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP1_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP1_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:55:29
% EndTime: 2019-03-10 00:55:29
% DurationCPUTime: 0.37s
% Computational Cost: add. (225->88), mult. (186->105), div. (0->0), fcn. (166->10), ass. (0->41)
t103 = rSges(6,3) + pkin(10);
t102 = rSges(7,3) + qJ(6) + pkin(10);
t77 = sin(qJ(1));
t80 = cos(qJ(1));
t101 = g(1) * t80 + g(2) * t77;
t81 = -pkin(8) - pkin(7);
t98 = rSges(3,3) + pkin(7);
t76 = sin(qJ(2));
t96 = t76 * pkin(2) + pkin(6);
t79 = cos(qJ(2));
t66 = t79 * pkin(2) + pkin(1);
t73 = qJ(2) + qJ(3);
t69 = qJ(4) + t73;
t63 = sin(t69);
t95 = rSges(5,2) * t63;
t64 = cos(t69);
t94 = t64 * t77;
t93 = t64 * t80;
t75 = sin(qJ(5));
t92 = t75 * t80;
t91 = t77 * t75;
t78 = cos(qJ(5));
t90 = t77 * t78;
t89 = t80 * t78;
t88 = rSges(4,3) - t81;
t68 = cos(t73);
t60 = pkin(3) * t68 + t66;
t72 = -pkin(9) + t81;
t86 = t77 * t60 + t80 * t72;
t67 = sin(t73);
t85 = pkin(3) * t67 + t96;
t59 = t80 * t60;
t84 = -t77 * t72 + t59;
t83 = rSges(3,1) * t79 - rSges(3,2) * t76 + pkin(1);
t82 = rSges(4,1) * t68 - rSges(4,2) * t67 + t66;
t65 = pkin(5) * t78 + pkin(4);
t57 = t64 * t89 + t91;
t56 = -t64 * t92 + t90;
t55 = t64 * t90 - t92;
t54 = -t64 * t91 - t89;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t80 - t77 * rSges(2,2)) + g(2) * (t77 * rSges(2,1) + rSges(2,2) * t80) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (rSges(3,1) * t76 + rSges(3,2) * t79 + pkin(6)) + (g(1) * t83 - g(2) * t98) * t80 + (g(1) * t98 + g(2) * t83) * t77) - m(4) * (g(3) * (rSges(4,1) * t67 + rSges(4,2) * t68 + t96) + (g(1) * t82 - g(2) * t88) * t80 + (g(1) * t88 + g(2) * t82) * t77) - m(5) * (g(1) * (rSges(5,1) * t93 - t80 * t95 + t59) + g(2) * (-rSges(5,3) * t80 + t86) + g(3) * (rSges(5,1) * t63 + rSges(5,2) * t64 + t85) + (g(1) * (rSges(5,3) - t72) + g(2) * (rSges(5,1) * t64 - t95)) * t77) - m(6) * (g(1) * (t57 * rSges(6,1) + t56 * rSges(6,2) + pkin(4) * t93 + t84) + g(2) * (rSges(6,1) * t55 + rSges(6,2) * t54 + pkin(4) * t94 + t86) + g(3) * (-t103 * t64 + t85) + (g(3) * (rSges(6,1) * t78 - rSges(6,2) * t75 + pkin(4)) + t101 * t103) * t63) - m(7) * (g(1) * (t57 * rSges(7,1) + t56 * rSges(7,2) + pkin(5) * t91 + t65 * t93 + t84) + g(2) * (t55 * rSges(7,1) + t54 * rSges(7,2) - pkin(5) * t92 + t65 * t94 + t86) + g(3) * (-t102 * t64 + t85) + (g(3) * (rSges(7,1) * t78 - rSges(7,2) * t75 + t65) + t101 * t102) * t63);
U  = t1;
