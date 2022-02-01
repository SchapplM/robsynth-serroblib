% Calculate potential energy for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% m [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR4_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:16:11
% EndTime: 2022-01-23 09:16:11
% DurationCPUTime: 0.43s
% Computational Cost: add. (152->93), mult. (170->118), div. (0->0), fcn. (158->10), ass. (0->40)
t65 = qJ(3) + pkin(6);
t88 = rSges(6,3) + pkin(7) + t65;
t66 = sin(qJ(1));
t67 = cos(qJ(1));
t68 = g(1) * t67 + g(2) * t66;
t61 = sin(pkin(9));
t85 = t61 * pkin(3);
t63 = cos(pkin(9));
t51 = t63 * pkin(3) + pkin(2);
t62 = sin(pkin(8));
t84 = rSges(3,2) * t62;
t64 = cos(pkin(8));
t83 = t64 * t66;
t82 = t64 * t67;
t60 = pkin(9) + qJ(4);
t54 = qJ(5) + t60;
t49 = sin(t54);
t81 = t66 * t49;
t50 = cos(t54);
t80 = t66 * t50;
t52 = sin(t60);
t79 = t66 * t52;
t53 = cos(t60);
t78 = t66 * t53;
t77 = t66 * t61;
t76 = t66 * t63;
t75 = t67 * t52;
t74 = t67 * t53;
t73 = t67 * t61;
t72 = t67 * t63;
t55 = t66 * qJ(2);
t70 = t67 * pkin(1) + t55;
t69 = t67 * qJ(2);
t57 = t66 * pkin(1);
t48 = qJ(2) + t85;
t47 = pkin(2) * t64 + t62 * qJ(3) + pkin(1);
t46 = pkin(4) * t52 + t85;
t45 = pkin(4) * t53 + t51;
t44 = t51 * t64 + t65 * t62 + pkin(1);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t67 * rSges(2,1) - t66 * rSges(2,2)) + g(2) * (t66 * rSges(2,1) + t67 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t66 * rSges(3,3) + t70) + g(2) * (rSges(3,1) * t83 - t66 * t84 + t57) + g(3) * (t62 * rSges(3,1) + t64 * rSges(3,2) + pkin(5)) + (g(1) * (rSges(3,1) * t64 - t84) + g(2) * (-rSges(3,3) - qJ(2))) * t67) - m(4) * (g(1) * (t47 * t67 + t55 + (t64 * t72 + t77) * rSges(4,1) + (-t64 * t73 + t76) * rSges(4,2)) + g(2) * (t47 * t66 - t69 + (t64 * t76 - t73) * rSges(4,1) + (-t64 * t77 - t72) * rSges(4,2)) + g(3) * (pkin(5) + (-rSges(4,3) - qJ(3)) * t64) + (g(3) * (rSges(4,1) * t63 - rSges(4,2) * t61 + pkin(2)) + t68 * rSges(4,3)) * t62) - m(5) * (g(1) * (t44 * t67 + t48 * t66 + (t64 * t74 + t79) * rSges(5,1) + (-t64 * t75 + t78) * rSges(5,2)) + g(2) * (t44 * t66 - t48 * t67 + (t64 * t78 - t75) * rSges(5,1) + (-t64 * t79 - t74) * rSges(5,2)) + g(3) * (pkin(5) + (-rSges(5,3) - t65) * t64) + (g(3) * (rSges(5,1) * t53 - rSges(5,2) * t52 + t51) + t68 * rSges(5,3)) * t62) - m(6) * (g(1) * (t45 * t82 + t66 * t46 + (t50 * t82 + t81) * rSges(6,1) + (-t49 * t82 + t80) * rSges(6,2) + t70) + g(2) * (t45 * t83 - t67 * t46 + t57 - t69 + (-t49 * t67 + t64 * t80) * rSges(6,1) + (-t50 * t67 - t64 * t81) * rSges(6,2)) + g(3) * (-t88 * t64 + pkin(5)) + (g(3) * (rSges(6,1) * t50 - rSges(6,2) * t49 + t45) + t68 * t88) * t62);
U = t1;
