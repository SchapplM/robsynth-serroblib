% Calculate potential energy for
% S6RPRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP9_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP9_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP9_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP9_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:26:03
% EndTime: 2019-03-09 06:26:04
% DurationCPUTime: 0.42s
% Computational Cost: add. (164->97), mult. (200->116), div. (0->0), fcn. (184->8), ass. (0->39)
t96 = rSges(5,3) + pkin(8);
t71 = -pkin(9) - pkin(8);
t95 = rSges(6,3) - t71;
t94 = rSges(7,3) + qJ(6) - t71;
t67 = sin(qJ(1));
t70 = cos(qJ(1));
t93 = -g(1) * t67 + g(2) * t70;
t92 = pkin(2) + pkin(6);
t68 = cos(qJ(4));
t54 = t68 * pkin(4) + pkin(3);
t69 = cos(qJ(3));
t88 = rSges(4,2) * t69;
t66 = sin(qJ(3));
t87 = t66 * t67;
t86 = t66 * t70;
t64 = qJ(4) + qJ(5);
t55 = sin(t64);
t85 = t67 * t55;
t56 = cos(t64);
t84 = t67 * t56;
t65 = sin(qJ(4));
t83 = t67 * t65;
t82 = t67 * t68;
t81 = t70 * t55;
t80 = t70 * t56;
t79 = t70 * t65;
t78 = t70 * t68;
t75 = t70 * pkin(1) + t67 * qJ(2);
t59 = t67 * pkin(1);
t74 = t67 * pkin(7) + t59;
t73 = t70 * pkin(7) + t75;
t72 = -t70 * qJ(2) + t74;
t53 = t65 * pkin(4) + pkin(5) * t55;
t52 = pkin(5) * t56 + t54;
t51 = -t66 * t80 + t85;
t50 = t66 * t81 + t84;
t49 = t66 * t84 + t81;
t48 = -t66 * t85 + t80;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t70 * rSges(2,1) - t67 * rSges(2,2)) + g(2) * (t67 * rSges(2,1) + t70 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (-t70 * rSges(3,2) + t67 * rSges(3,3) + t75) + g(2) * (-t67 * rSges(3,2) + t59 + (-rSges(3,3) - qJ(2)) * t70) + g(3) * (pkin(6) + rSges(3,1))) - m(4) * (g(1) * (rSges(4,1) * t87 + t67 * t88 + t73) + g(2) * (t67 * rSges(4,3) + t74) + g(3) * (t69 * rSges(4,1) - t66 * rSges(4,2) + t92) + (g(1) * rSges(4,3) + g(2) * (-rSges(4,1) * t66 - qJ(2) - t88)) * t70) - m(5) * (g(1) * (pkin(3) * t87 + (t66 * t82 + t79) * rSges(5,1) + (-t66 * t83 + t78) * rSges(5,2) + t73) + g(2) * (-pkin(3) * t86 + (-t66 * t78 + t83) * rSges(5,1) + (t66 * t79 + t82) * rSges(5,2) + t72) + g(3) * (t96 * t66 + t92) + (g(3) * (rSges(5,1) * t68 - rSges(5,2) * t65 + pkin(3)) + t93 * t96) * t69) - m(6) * (g(1) * (t49 * rSges(6,1) + t48 * rSges(6,2) + pkin(4) * t79 + t54 * t87 + t73) + g(2) * (t51 * rSges(6,1) + t50 * rSges(6,2) + pkin(4) * t83 - t54 * t86 + t72) + g(3) * (t95 * t66 + t92) + (g(3) * (rSges(6,1) * t56 - rSges(6,2) * t55 + t54) + t93 * t95) * t69) - m(7) * (g(1) * (t49 * rSges(7,1) + t48 * rSges(7,2) + t52 * t87 + t70 * t53 + t73) + g(2) * (t51 * rSges(7,1) + t50 * rSges(7,2) - t52 * t86 + t67 * t53 + t72) + g(3) * (t94 * t66 + t92) + (g(3) * (rSges(7,1) * t56 - rSges(7,2) * t55 + t52) + t93 * t94) * t69);
U  = t1;
