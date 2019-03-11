% Calculate potential energy for
% S6RPRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
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
% Datum: 2019-03-09 05:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPR10_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR10_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR10_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR10_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:34:17
% EndTime: 2019-03-09 05:34:18
% DurationCPUTime: 0.42s
% Computational Cost: add. (153->98), mult. (262->124), div. (0->0), fcn. (266->8), ass. (0->33)
t89 = -pkin(9) - rSges(7,3);
t92 = pkin(2) + pkin(6);
t69 = sin(qJ(1));
t91 = g(1) * t69;
t73 = cos(qJ(1));
t90 = g(2) * t73;
t72 = cos(qJ(3));
t88 = rSges(4,2) * t72;
t68 = sin(qJ(3));
t87 = t68 * t69;
t86 = t68 * t73;
t67 = sin(qJ(4));
t85 = t69 * t67;
t71 = cos(qJ(4));
t84 = t69 * t71;
t83 = t73 * t71;
t82 = t73 * pkin(1) + t69 * qJ(2);
t62 = t69 * pkin(1);
t81 = t69 * pkin(7) + t62;
t80 = t73 * pkin(7) + t82;
t79 = t72 * pkin(3) + t68 * pkin(8) + t92;
t78 = pkin(3) * t87 + t80;
t77 = t79 + (pkin(4) * t71 + qJ(5) * t67) * t72;
t50 = t68 * t85 - t83;
t51 = t67 * t73 + t68 * t84;
t76 = t51 * pkin(4) + t50 * qJ(5) + t78;
t75 = -pkin(3) * t86 + t81 + (t72 * pkin(8) - qJ(2)) * t73;
t52 = t67 * t86 + t84;
t53 = -t68 * t83 + t85;
t74 = t53 * pkin(4) - t52 * qJ(5) + t75;
t70 = cos(qJ(6));
t66 = sin(qJ(6));
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t73 - t69 * rSges(2,2)) + g(2) * (t69 * rSges(2,1) + rSges(2,2) * t73) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (-rSges(3,2) * t73 + t69 * rSges(3,3) + t82) + g(2) * (-t69 * rSges(3,2) + t62 + (-rSges(3,3) - qJ(2)) * t73) + g(3) * (pkin(6) + rSges(3,1))) - m(4) * (g(1) * (rSges(4,1) * t87 + t69 * t88 + t80) + g(2) * (t69 * rSges(4,3) + t81) + g(3) * (rSges(4,1) * t72 - rSges(4,2) * t68 + t92) + (g(1) * rSges(4,3) + g(2) * (-rSges(4,1) * t68 - qJ(2) - t88)) * t73) - m(5) * (g(1) * (rSges(5,1) * t51 - rSges(5,2) * t50 + t78) + g(2) * (t53 * rSges(5,1) + t52 * rSges(5,2) + t75) + g(3) * (rSges(5,3) * t68 + t79) + (rSges(5,3) * t90 + g(3) * (rSges(5,1) * t71 - rSges(5,2) * t67) + (-rSges(5,3) - pkin(8)) * t91) * t72) - m(6) * (g(1) * (rSges(6,1) * t51 + rSges(6,3) * t50 + t76) + g(2) * (t53 * rSges(6,1) - t52 * rSges(6,3) + t74) + g(3) * (rSges(6,2) * t68 + t77) + (rSges(6,2) * t90 + g(3) * (rSges(6,1) * t71 + rSges(6,3) * t67) + (-rSges(6,2) - pkin(8)) * t91) * t72) - m(7) * (g(1) * (t51 * pkin(5) + (t50 * t66 + t51 * t70) * rSges(7,1) + (t50 * t70 - t51 * t66) * rSges(7,2) + t76) + g(2) * (t53 * pkin(5) + (-t52 * t66 + t53 * t70) * rSges(7,1) + (-t52 * t70 - t53 * t66) * rSges(7,2) + t74) + g(3) * (t89 * t68 + t77) + (g(3) * (t71 * pkin(5) + (t66 * t67 + t70 * t71) * rSges(7,1) + (-t66 * t71 + t67 * t70) * rSges(7,2)) + t89 * t90 + (-pkin(8) - t89) * t91) * t72);
U  = t1;
