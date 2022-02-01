% Calculate potential energy for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR5_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:24:48
% EndTime: 2022-01-23 09:24:49
% DurationCPUTime: 0.47s
% Computational Cost: add. (152->88), mult. (173->109), div. (0->0), fcn. (161->10), ass. (0->41)
t93 = rSges(4,3) + pkin(6);
t65 = qJ(4) + pkin(6);
t92 = rSges(6,3) + pkin(7) + t65;
t67 = sin(qJ(1));
t69 = cos(qJ(1));
t91 = g(1) * t69 + g(2) * t67;
t66 = sin(qJ(3));
t88 = t66 * pkin(3);
t68 = cos(qJ(3));
t87 = t68 * pkin(3);
t64 = cos(pkin(8));
t86 = pkin(2) * t64 + pkin(1);
t63 = sin(pkin(8));
t85 = t63 * pkin(2) + pkin(5);
t84 = rSges(3,2) * t63;
t83 = t64 * t67;
t62 = qJ(3) + pkin(9);
t55 = qJ(5) + t62;
t50 = sin(t55);
t82 = t67 * t50;
t51 = cos(t55);
t81 = t67 * t51;
t53 = sin(t62);
t80 = t67 * t53;
t54 = cos(t62);
t79 = t67 * t54;
t78 = t69 * t50;
t77 = t69 * t51;
t76 = t69 * t53;
t75 = t69 * t54;
t57 = t67 * qJ(2);
t73 = t69 * pkin(1) + t57;
t72 = rSges(4,1) * t68 - rSges(4,2) * t66;
t71 = rSges(4,1) * t66 + rSges(4,2) * t68;
t70 = t63 * t93 + t64 * t72 + t86;
t59 = t67 * pkin(1);
t52 = qJ(2) + t88;
t48 = pkin(4) * t53 + t88;
t47 = pkin(4) * t54 + pkin(2) + t87;
t46 = t63 * t65 + t64 * t87 + t86;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + rSges(1,3) * g(3)) - m(2) * (g(1) * (rSges(2,1) * t69 - rSges(2,2) * t67) + g(2) * (rSges(2,1) * t67 + rSges(2,2) * t69) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t67 * rSges(3,3) + t73) + g(2) * (rSges(3,1) * t83 - t67 * t84 + t59) + g(3) * (rSges(3,1) * t63 + rSges(3,2) * t64 + pkin(5)) + (g(1) * (rSges(3,1) * t64 - t84) + g(2) * (-rSges(3,3) - qJ(2))) * t69) - m(4) * (g(1) * (t67 * t71 + t69 * t70 + t57) + g(2) * ((-qJ(2) - t71) * t69 + t70 * t67) + g(3) * (t72 * t63 - t64 * t93 + t85)) - m(5) * (g(1) * (t46 * t69 + t52 * t67 + (t64 * t75 + t80) * rSges(5,1) + (-t64 * t76 + t79) * rSges(5,2)) + g(2) * (t46 * t67 - t52 * t69 + (t64 * t79 - t76) * rSges(5,1) + (-t64 * t80 - t75) * rSges(5,2)) + g(3) * (t85 + (-rSges(5,3) - t65) * t64) + (g(3) * (rSges(5,1) * t54 - rSges(5,2) * t53 + t87) + t91 * rSges(5,3)) * t63) - m(6) * (g(1) * (t69 * t64 * t47 + t67 * t48 + (t64 * t77 + t82) * rSges(6,1) + (-t64 * t78 + t81) * rSges(6,2) + t73) + g(2) * (t47 * t83 + t59 + (t64 * t81 - t78) * rSges(6,1) + (-t64 * t82 - t77) * rSges(6,2) + (-t48 - qJ(2)) * t69) + g(3) * (-t92 * t64 + pkin(5)) + (g(3) * (rSges(6,1) * t51 - rSges(6,2) * t50 + t47) + t91 * t92) * t63);
U = t1;
