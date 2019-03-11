% Calculate potential energy for
% S6RPRRRP1
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
% Datum: 2019-03-09 05:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP1_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP1_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:55:21
% EndTime: 2019-03-09 05:55:21
% DurationCPUTime: 0.31s
% Computational Cost: add. (238->87), mult. (185->105), div. (0->0), fcn. (169->10), ass. (0->37)
t95 = rSges(7,1) + pkin(5);
t94 = rSges(7,3) + qJ(6);
t93 = rSges(4,3) + pkin(7);
t71 = qJ(1) + pkin(10);
t64 = sin(t71);
t72 = qJ(3) + qJ(4);
t66 = sin(t72);
t92 = t64 * t66;
t65 = cos(t71);
t91 = t65 * t66;
t67 = cos(t72);
t90 = t65 * t67;
t73 = sin(qJ(5));
t89 = t67 * t73;
t76 = cos(qJ(5));
t88 = t67 * t76;
t87 = pkin(6) + qJ(2);
t77 = cos(qJ(3));
t63 = pkin(3) * t77 + pkin(2);
t78 = cos(qJ(1));
t70 = t78 * pkin(1);
t86 = t65 * t63 + t70;
t74 = sin(qJ(3));
t85 = t74 * pkin(3) + t87;
t75 = sin(qJ(1));
t69 = t75 * pkin(1);
t79 = -pkin(8) - pkin(7);
t84 = t64 * t63 + t65 * t79 + t69;
t83 = t66 * pkin(4) + t85;
t82 = t64 * t67 * pkin(4) + pkin(9) * t92 + t84;
t81 = rSges(4,1) * t77 - rSges(4,2) * t74 + pkin(2);
t80 = pkin(4) * t90 + pkin(9) * t91 - t64 * t79 + t86;
t53 = t64 * t73 + t65 * t88;
t52 = -t64 * t76 + t65 * t89;
t51 = t64 * t88 - t65 * t73;
t50 = t64 * t89 + t65 * t76;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t78 - rSges(2,2) * t75) + g(2) * (rSges(2,1) * t75 + rSges(2,2) * t78) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t65 - rSges(3,2) * t64 + t70) + g(2) * (rSges(3,1) * t64 + rSges(3,2) * t65 + t69) + g(3) * (rSges(3,3) + t87)) - m(4) * (g(1) * t70 + g(2) * t69 + g(3) * (rSges(4,1) * t74 + rSges(4,2) * t77 + t87) + (g(1) * t81 - g(2) * t93) * t65 + (g(1) * t93 + g(2) * t81) * t64) - m(5) * (g(1) * (rSges(5,1) * t90 - rSges(5,2) * t91 + t86) + g(2) * (-rSges(5,3) * t65 + t84) + g(3) * (rSges(5,1) * t66 + rSges(5,2) * t67 + t85) + (g(1) * (rSges(5,3) - t79) + g(2) * (rSges(5,1) * t67 - rSges(5,2) * t66)) * t64) - m(6) * (g(1) * (t53 * rSges(6,1) - t52 * rSges(6,2) + rSges(6,3) * t91 + t80) + g(2) * (rSges(6,1) * t51 - rSges(6,2) * t50 + rSges(6,3) * t92 + t82) + g(3) * ((-rSges(6,3) - pkin(9)) * t67 + (rSges(6,1) * t76 - rSges(6,2) * t73) * t66 + t83)) - m(7) * (g(1) * (t94 * t52 + t95 * t53 + t80) + g(2) * (t94 * t50 + t95 * t51 + t82) + g(3) * (t83 + (-rSges(7,2) - pkin(9)) * t67) + (g(3) * (t94 * t73 + t95 * t76) + (g(1) * t65 + g(2) * t64) * rSges(7,2)) * t66);
U  = t1;
