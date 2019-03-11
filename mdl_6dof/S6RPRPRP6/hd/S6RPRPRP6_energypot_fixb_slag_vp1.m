% Calculate potential energy for
% S6RPRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-03-09 03:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRP6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP6_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP6_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:17:53
% EndTime: 2019-03-09 03:17:53
% DurationCPUTime: 0.38s
% Computational Cost: add. (196->93), mult. (204->111), div. (0->0), fcn. (184->8), ass. (0->38)
t102 = rSges(6,3) + pkin(8);
t101 = rSges(7,3) + qJ(6) + pkin(8);
t76 = sin(qJ(1));
t78 = cos(qJ(1));
t100 = g(1) * t78 + g(2) * t76;
t71 = sin(pkin(9));
t96 = t71 * pkin(2) + pkin(6);
t70 = pkin(9) + qJ(3);
t67 = sin(t70);
t95 = t67 * t78;
t68 = cos(t70);
t94 = t68 * t78;
t75 = sin(qJ(5));
t93 = t76 * t75;
t77 = cos(qJ(5));
t92 = t76 * t77;
t91 = t78 * t75;
t90 = t78 * t77;
t72 = cos(pkin(9));
t64 = t72 * pkin(2) + pkin(1);
t74 = -pkin(7) - qJ(2);
t88 = t76 * t64 + t78 * t74;
t87 = qJ(4) * t67;
t86 = rSges(3,3) + qJ(2);
t85 = t67 * pkin(3) + t96;
t84 = t67 * t93;
t83 = t67 * t91;
t60 = t78 * t64;
t82 = pkin(3) * t94 + t78 * t87 + t60;
t81 = t88 + (pkin(3) * t68 + t87) * t76;
t80 = -t76 * t74 + t82;
t79 = rSges(3,1) * t72 - rSges(3,2) * t71 + pkin(1);
t66 = t77 * pkin(5) + pkin(4);
t56 = t84 - t90;
t55 = t67 * t92 + t91;
t54 = t83 + t92;
t53 = t67 * t90 - t93;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t78 * rSges(2,1) - t76 * rSges(2,2)) + g(2) * (t76 * rSges(2,1) + t78 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (t71 * rSges(3,1) + t72 * rSges(3,2) + pkin(6)) + (g(1) * t79 - g(2) * t86) * t78 + (g(1) * t86 + g(2) * t79) * t76) - m(4) * (g(1) * (rSges(4,1) * t94 - rSges(4,2) * t95 + t60) + g(2) * (-t78 * rSges(4,3) + t88) + g(3) * (t67 * rSges(4,1) + t68 * rSges(4,2) + t96) + (g(1) * (rSges(4,3) - t74) + g(2) * (rSges(4,1) * t68 - rSges(4,2) * t67)) * t76) - m(5) * (g(1) * (-rSges(5,2) * t94 + rSges(5,3) * t95 + t82) + g(2) * (-t78 * rSges(5,1) + t81) + g(3) * (-t67 * rSges(5,2) + (-rSges(5,3) - qJ(4)) * t68 + t85) + (g(1) * (rSges(5,1) - t74) + g(2) * (-rSges(5,2) * t68 + rSges(5,3) * t67)) * t76) - m(6) * (g(1) * (t54 * rSges(6,1) + t53 * rSges(6,2) + t76 * pkin(4) + t80) + g(2) * (t56 * rSges(6,1) + t55 * rSges(6,2) - t78 * pkin(4) + t81) + g(3) * (t102 * t67 + t85) + (g(3) * (-rSges(6,1) * t75 - rSges(6,2) * t77 - qJ(4)) + t100 * t102) * t68) - m(7) * (g(1) * (t54 * rSges(7,1) + t53 * rSges(7,2) + pkin(5) * t83 + t76 * t66 + t80) + g(2) * (t56 * rSges(7,1) + t55 * rSges(7,2) + pkin(5) * t84 - t78 * t66 + t81) + g(3) * (t101 * t67 + t85) + (g(3) * (-rSges(7,2) * t77 - qJ(4) + (-rSges(7,1) - pkin(5)) * t75) + t100 * t101) * t68);
U  = t1;
