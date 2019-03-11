% Calculate potential energy for
% S6RRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPR4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR4_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR4_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:06:06
% EndTime: 2019-03-09 22:06:06
% DurationCPUTime: 0.50s
% Computational Cost: add. (236->104), mult. (212->128), div. (0->0), fcn. (196->12), ass. (0->43)
t102 = rSges(5,3) + pkin(9);
t70 = -qJ(5) - pkin(9);
t101 = rSges(6,3) - t70;
t100 = rSges(7,3) + pkin(10) - t70;
t73 = sin(qJ(1));
t76 = cos(qJ(1));
t99 = g(1) * t76 + g(2) * t73;
t96 = rSges(3,3) + pkin(7);
t72 = sin(qJ(2));
t94 = t72 * pkin(2) + pkin(6);
t74 = cos(qJ(4));
t57 = t74 * pkin(4) + pkin(3);
t69 = qJ(2) + qJ(3);
t63 = sin(t69);
t93 = rSges(4,2) * t63;
t64 = cos(t69);
t92 = t64 * t73;
t91 = t64 * t76;
t71 = sin(qJ(4));
t90 = t71 * t73;
t89 = t71 * t76;
t68 = qJ(4) + pkin(11);
t60 = sin(t68);
t88 = t73 * t60;
t61 = cos(t68);
t87 = t73 * t61;
t86 = t73 * t74;
t85 = t74 * t76;
t84 = t76 * t60;
t83 = t76 * t61;
t75 = cos(qJ(2));
t58 = pkin(2) * t75 + pkin(1);
t77 = -pkin(8) - pkin(7);
t80 = t73 * t58 + t76 * t77;
t54 = t76 * t58;
t79 = -t73 * t77 + t54;
t78 = rSges(3,1) * t75 - rSges(3,2) * t72 + pkin(1);
t62 = qJ(6) + t68;
t56 = cos(t62);
t55 = sin(t62);
t52 = t71 * pkin(4) + pkin(5) * t60;
t51 = pkin(5) * t61 + t57;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t76 - rSges(2,2) * t73) + g(2) * (rSges(2,1) * t73 + rSges(2,2) * t76) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (rSges(3,1) * t72 + t75 * rSges(3,2) + pkin(6)) + (g(1) * t78 - g(2) * t96) * t76 + (g(1) * t96 + g(2) * t78) * t73) - m(4) * (g(1) * (rSges(4,1) * t91 - t76 * t93 + t54) + g(2) * (-rSges(4,3) * t76 + t80) + g(3) * (rSges(4,1) * t63 + rSges(4,2) * t64 + t94) + (g(1) * (rSges(4,3) - t77) + g(2) * (rSges(4,1) * t64 - t93)) * t73) - m(5) * (g(1) * (pkin(3) * t91 + (t64 * t85 + t90) * rSges(5,1) + (-t64 * t89 + t86) * rSges(5,2) + t79) + g(2) * (pkin(3) * t92 + (t64 * t86 - t89) * rSges(5,1) + (-t64 * t90 - t85) * rSges(5,2) + t80) + g(3) * (-t102 * t64 + t94) + (g(3) * (rSges(5,1) * t74 - rSges(5,2) * t71 + pkin(3)) + t99 * t102) * t63) - m(6) * (g(1) * (t57 * t91 + pkin(4) * t90 + (t64 * t83 + t88) * rSges(6,1) + (-t64 * t84 + t87) * rSges(6,2) + t79) + g(2) * (t57 * t92 - pkin(4) * t89 + (t64 * t87 - t84) * rSges(6,1) + (-t64 * t88 - t83) * rSges(6,2) + t80) + g(3) * (-t101 * t64 + t94) + (g(3) * (rSges(6,1) * t61 - rSges(6,2) * t60 + t57) + t99 * t101) * t63) - m(7) * (g(1) * (t51 * t91 + t73 * t52 + (t55 * t73 + t56 * t91) * rSges(7,1) + (-t55 * t91 + t56 * t73) * rSges(7,2) + t79) + g(2) * (t51 * t92 - t76 * t52 + (-t55 * t76 + t56 * t92) * rSges(7,1) + (-t55 * t92 - t56 * t76) * rSges(7,2) + t80) + g(3) * (-t100 * t64 + t94) + (g(3) * (rSges(7,1) * t56 - rSges(7,2) * t55 + t51) + t99 * t100) * t63);
U  = t1;
