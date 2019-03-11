% Calculate potential energy for
% S6RPRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-03-09 03:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRP9_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP9_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP9_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP9_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:27:29
% EndTime: 2019-03-09 03:27:29
% DurationCPUTime: 0.44s
% Computational Cost: add. (171->98), mult. (215->115), div. (0->0), fcn. (203->8), ass. (0->41)
t98 = rSges(7,1) + pkin(5);
t70 = -pkin(8) - qJ(4);
t97 = rSges(7,2) - t70;
t96 = rSges(5,3) + qJ(4);
t95 = rSges(7,3) + qJ(6);
t94 = pkin(2) + pkin(6);
t72 = sin(qJ(1));
t93 = g(1) * t72;
t74 = cos(qJ(1));
t92 = g(2) * t74;
t68 = sin(pkin(9));
t91 = t68 * t74;
t71 = sin(qJ(3));
t90 = t71 * t72;
t89 = t71 * t74;
t67 = pkin(9) + qJ(5);
t60 = sin(t67);
t88 = t72 * t60;
t61 = cos(t67);
t87 = t72 * t61;
t86 = t72 * t68;
t69 = cos(pkin(9));
t85 = t72 * t69;
t73 = cos(qJ(3));
t84 = t72 * t73;
t83 = t74 * t61;
t82 = rSges(6,3) - t70;
t81 = t74 * pkin(1) + t72 * qJ(2);
t64 = t72 * pkin(1);
t80 = t72 * pkin(7) + t64;
t59 = pkin(4) * t69 + pkin(3);
t78 = t73 * t59 + t94;
t77 = t74 * pkin(7) + t81;
t76 = -t74 * qJ(2) + t80;
t75 = pkin(4) * t91 + t59 * t90 + t70 * t84 + t77;
t57 = pkin(4) * t86;
t52 = -t71 * t83 + t88;
t51 = t60 * t89 + t87;
t50 = t60 * t74 + t71 * t87;
t49 = t71 * t88 - t83;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t74 - t72 * rSges(2,2)) + g(2) * (t72 * rSges(2,1) + rSges(2,2) * t74) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (-rSges(3,2) * t74 + t72 * rSges(3,3) + t81) + g(2) * (-t72 * rSges(3,2) + t64 + (-rSges(3,3) - qJ(2)) * t74) + g(3) * (pkin(6) + rSges(3,1))) - m(4) * (g(1) * (rSges(4,1) * t90 + rSges(4,2) * t84 + t77) + g(2) * (t72 * rSges(4,3) + t80) + g(3) * (rSges(4,1) * t73 - rSges(4,2) * t71 + t94) + (g(1) * rSges(4,3) + g(2) * (-rSges(4,1) * t71 - rSges(4,2) * t73 - qJ(2))) * t74) - m(5) * (g(1) * (pkin(3) * t90 + (t71 * t85 + t91) * rSges(5,1) + (t69 * t74 - t71 * t86) * rSges(5,2) + t77) + g(2) * (-pkin(3) * t89 + (-t69 * t89 + t86) * rSges(5,1) + (t68 * t89 + t85) * rSges(5,2) + t76) + g(3) * (t96 * t71 + t94) + (g(3) * (rSges(5,1) * t69 - rSges(5,2) * t68 + pkin(3)) + (-t93 + t92) * t96) * t73) - m(6) * (g(1) * (rSges(6,1) * t50 - rSges(6,2) * t49 - rSges(6,3) * t84 + t75) + g(2) * (t52 * rSges(6,1) + t51 * rSges(6,2) + t57 + t80) + g(3) * ((rSges(6,1) * t61 - rSges(6,2) * t60) * t73 + t82 * t71 + t78) + (-t59 * t71 + t82 * t73 - qJ(2)) * t92) - m(7) * (g(1) * (t95 * t49 + t98 * t50 + t75) + g(2) * (-t95 * t51 + t98 * t52 - t59 * t89 + t57 + t76) + g(3) * (t97 * t71 + t78) + (-rSges(7,2) * t93 + g(3) * (t95 * t60 + t98 * t61) + t97 * t92) * t73);
U  = t1;
