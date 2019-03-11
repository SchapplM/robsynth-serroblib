% Calculate potential energy for
% S6RPRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 05:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPR8_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR8_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR8_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR8_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:22:30
% EndTime: 2019-03-09 05:22:31
% DurationCPUTime: 0.48s
% Computational Cost: add. (174->102), mult. (200->124), div. (0->0), fcn. (184->10), ass. (0->34)
t90 = rSges(5,3) + pkin(8);
t63 = -qJ(5) - pkin(8);
t89 = rSges(6,3) - t63;
t88 = rSges(7,3) + pkin(9) - t63;
t66 = sin(qJ(1));
t69 = cos(qJ(1));
t87 = -g(1) * t66 + g(2) * t69;
t86 = pkin(2) + pkin(6);
t67 = cos(qJ(4));
t51 = t67 * pkin(4) + pkin(3);
t68 = cos(qJ(3));
t82 = rSges(4,2) * t68;
t64 = sin(qJ(4));
t81 = t66 * t64;
t65 = sin(qJ(3));
t80 = t66 * t65;
t79 = t66 * t67;
t78 = t69 * t64;
t77 = t69 * t65;
t76 = t69 * t67;
t73 = t69 * pkin(1) + t66 * qJ(2);
t57 = t66 * pkin(1);
t72 = t66 * pkin(7) + t57;
t62 = qJ(4) + pkin(10);
t71 = t69 * pkin(7) + t73;
t70 = -t69 * qJ(2) + t72;
t54 = qJ(6) + t62;
t53 = cos(t62);
t52 = sin(t62);
t50 = cos(t54);
t49 = sin(t54);
t48 = t64 * pkin(4) + pkin(5) * t52;
t47 = pkin(5) * t53 + t51;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t69 * rSges(2,1) - t66 * rSges(2,2)) + g(2) * (t66 * rSges(2,1) + t69 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (-t69 * rSges(3,2) + t66 * rSges(3,3) + t73) + g(2) * (-t66 * rSges(3,2) + t57 + (-rSges(3,3) - qJ(2)) * t69) + g(3) * (pkin(6) + rSges(3,1))) - m(4) * (g(1) * (rSges(4,1) * t80 + t66 * t82 + t71) + g(2) * (t66 * rSges(4,3) + t72) + g(3) * (t68 * rSges(4,1) - t65 * rSges(4,2) + t86) + (g(1) * rSges(4,3) + g(2) * (-rSges(4,1) * t65 - qJ(2) - t82)) * t69) - m(5) * (g(1) * (pkin(3) * t80 + (t65 * t79 + t78) * rSges(5,1) + (-t64 * t80 + t76) * rSges(5,2) + t71) + g(2) * (-pkin(3) * t77 + (-t65 * t76 + t81) * rSges(5,1) + (t64 * t77 + t79) * rSges(5,2) + t70) + g(3) * (t90 * t65 + t86) + (g(3) * (rSges(5,1) * t67 - rSges(5,2) * t64 + pkin(3)) + t87 * t90) * t68) - m(6) * (g(1) * (t51 * t80 + pkin(4) * t78 + (t69 * t52 + t53 * t80) * rSges(6,1) + (-t52 * t80 + t69 * t53) * rSges(6,2) + t71) + g(2) * (-t51 * t77 + pkin(4) * t81 + (t66 * t52 - t53 * t77) * rSges(6,1) + (t52 * t77 + t66 * t53) * rSges(6,2) + t70) + g(3) * (t89 * t65 + t86) + (g(3) * (rSges(6,1) * t53 - rSges(6,2) * t52 + t51) + t87 * t89) * t68) - m(7) * (g(1) * (t47 * t80 + t69 * t48 + (t69 * t49 + t50 * t80) * rSges(7,1) + (-t49 * t80 + t69 * t50) * rSges(7,2) + t71) + g(2) * (-t47 * t77 + t66 * t48 + (t66 * t49 - t50 * t77) * rSges(7,1) + (t49 * t77 + t66 * t50) * rSges(7,2) + t70) + g(3) * (t88 * t65 + t86) + (g(3) * (rSges(7,1) * t50 - rSges(7,2) * t49 + t47) + t87 * t88) * t68);
U  = t1;
