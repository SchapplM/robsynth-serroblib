% Calculate potential energy for
% S6RPRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-03-09 07:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRR9_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR9_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR9_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR9_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:22:57
% EndTime: 2019-03-09 07:22:57
% DurationCPUTime: 0.48s
% Computational Cost: add. (174->102), mult. (200->124), div. (0->0), fcn. (184->10), ass. (0->33)
t89 = rSges(5,3) + pkin(8);
t69 = -pkin(9) - pkin(8);
t88 = rSges(6,3) - t69;
t87 = rSges(7,3) + pkin(10) - t69;
t65 = sin(qJ(1));
t68 = cos(qJ(1));
t86 = -g(1) * t65 + g(2) * t68;
t85 = pkin(2) + pkin(6);
t66 = cos(qJ(4));
t51 = t66 * pkin(4) + pkin(3);
t67 = cos(qJ(3));
t81 = rSges(4,2) * t67;
t63 = sin(qJ(4));
t80 = t65 * t63;
t64 = sin(qJ(3));
t79 = t65 * t64;
t78 = t65 * t66;
t77 = t68 * t63;
t76 = t68 * t64;
t73 = t68 * pkin(1) + t65 * qJ(2);
t57 = t65 * pkin(1);
t72 = t65 * pkin(7) + t57;
t62 = qJ(4) + qJ(5);
t71 = t68 * pkin(7) + t73;
t70 = -t68 * qJ(2) + t72;
t55 = qJ(6) + t62;
t53 = cos(t62);
t52 = sin(t62);
t50 = cos(t55);
t49 = sin(t55);
t48 = t63 * pkin(4) + pkin(5) * t52;
t47 = pkin(5) * t53 + t51;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t68 - rSges(2,2) * t65) + g(2) * (rSges(2,1) * t65 + rSges(2,2) * t68) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (-rSges(3,2) * t68 + rSges(3,3) * t65 + t73) + g(2) * (-rSges(3,2) * t65 + t57 + (-rSges(3,3) - qJ(2)) * t68) + g(3) * (pkin(6) + rSges(3,1))) - m(4) * (g(1) * (rSges(4,1) * t79 + t65 * t81 + t71) + g(2) * (rSges(4,3) * t65 + t72) + g(3) * (rSges(4,1) * t67 - rSges(4,2) * t64 + t85) + (g(1) * rSges(4,3) + g(2) * (-rSges(4,1) * t64 - qJ(2) - t81)) * t68) - m(5) * (g(1) * (pkin(3) * t79 + (t64 * t78 + t77) * rSges(5,1) + (-t63 * t79 + t66 * t68) * rSges(5,2) + t71) + g(2) * (-pkin(3) * t76 + (-t66 * t76 + t80) * rSges(5,1) + (t63 * t76 + t78) * rSges(5,2) + t70) + g(3) * (t89 * t64 + t85) + (g(3) * (rSges(5,1) * t66 - rSges(5,2) * t63 + pkin(3)) + t86 * t89) * t67) - m(6) * (g(1) * (t51 * t79 + pkin(4) * t77 + (t52 * t68 + t53 * t79) * rSges(6,1) + (-t52 * t79 + t53 * t68) * rSges(6,2) + t71) + g(2) * (-t51 * t76 + pkin(4) * t80 + (t52 * t65 - t53 * t76) * rSges(6,1) + (t52 * t76 + t53 * t65) * rSges(6,2) + t70) + g(3) * (t88 * t64 + t85) + (g(3) * (rSges(6,1) * t53 - rSges(6,2) * t52 + t51) + t86 * t88) * t67) - m(7) * (g(1) * (t47 * t79 + t68 * t48 + (t49 * t68 + t50 * t79) * rSges(7,1) + (-t49 * t79 + t50 * t68) * rSges(7,2) + t71) + g(2) * (-t47 * t76 + t65 * t48 + (t49 * t65 - t50 * t76) * rSges(7,1) + (t49 * t76 + t50 * t65) * rSges(7,2) + t70) + g(3) * (t87 * t64 + t85) + (g(3) * (rSges(7,1) * t50 - rSges(7,2) * t49 + t47) + t86 * t87) * t67);
U  = t1;
