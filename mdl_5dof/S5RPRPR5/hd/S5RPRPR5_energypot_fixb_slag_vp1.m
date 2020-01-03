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
% m_mdh [6x1]
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
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:41:19
% EndTime: 2020-01-03 11:41:19
% DurationCPUTime: 0.51s
% Computational Cost: add. (152->78), mult. (180->96), div. (0->0), fcn. (168->10), ass. (0->34)
t90 = rSges(4,3) + pkin(6);
t62 = -qJ(4) - pkin(6);
t89 = rSges(5,3) - t62;
t88 = rSges(6,3) + pkin(7) - t62;
t64 = sin(qJ(1));
t66 = cos(qJ(1));
t87 = g(2) * t64 - g(3) * t66;
t63 = sin(qJ(3));
t86 = -t63 * pkin(3) - qJ(2);
t85 = g(3) * pkin(1);
t56 = t64 * pkin(1);
t84 = g(2) * t56;
t61 = cos(pkin(8));
t83 = g(2) * t61;
t81 = g(3) * t61;
t65 = cos(qJ(3));
t52 = t65 * pkin(3) + pkin(2);
t77 = t64 * t63;
t76 = t64 * t65;
t75 = t66 * t61;
t72 = -rSges(3,3) - qJ(2);
t59 = qJ(3) + pkin(9);
t60 = sin(pkin(8));
t71 = rSges(3,1) * t61 - rSges(3,2) * t60;
t53 = sin(t59);
t54 = cos(t59);
t70 = rSges(5,1) * t54 - rSges(5,2) * t53 + t52;
t55 = qJ(5) + t59;
t50 = sin(t55);
t51 = cos(t55);
t69 = rSges(6,1) * t51 - rSges(6,2) * t50 + pkin(4) * t54 + t52;
t68 = -rSges(6,1) * t50 - rSges(6,2) * t51 - pkin(4) * t53 + t86;
t67 = -rSges(5,1) * t53 - rSges(5,2) * t54 + t86;
t1 = -m(1) * (g(1) * rSges(1,1) + rSges(1,2) * g(2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (pkin(5) + rSges(2,3)) + g(2) * (rSges(2,1) * t64 + rSges(2,2) * t66) + g(3) * (-rSges(2,1) * t66 + rSges(2,2) * t64)) - m(3) * (g(1) * (rSges(3,1) * t60 + rSges(3,2) * t61 + pkin(5)) + t84 + (g(2) * t71 + g(3) * t72) * t64 + (g(2) * t72 + g(3) * (-pkin(1) - t71)) * t66) - m(4) * (g(1) * (-t90 * t61 + pkin(5)) + g(2) * (t64 * t61 * pkin(2) + t56 - t66 * qJ(2) + (t61 * t76 - t63 * t66) * rSges(4,1) + (-t61 * t77 - t65 * t66) * rSges(4,2)) + g(3) * (-pkin(2) * t75 - t66 * pkin(1) - t64 * qJ(2) + (-t65 * t75 - t77) * rSges(4,1) + (t63 * t75 - t76) * rSges(4,2)) + (g(1) * (rSges(4,1) * t65 - rSges(4,2) * t63 + pkin(2)) + t87 * t90) * t60) - m(5) * (g(1) * (-t61 * t89 + pkin(5)) + t84 + (g(3) * t67 + t70 * t83) * t64 + (g(2) * t67 - t70 * t81 - t85) * t66 + (g(1) * t70 + t87 * t89) * t60) - m(6) * (g(1) * (-t61 * t88 + pkin(5)) + t84 + (g(3) * t68 + t69 * t83) * t64 + (g(2) * t68 - t69 * t81 - t85) * t66 + (g(1) * t69 + t87 * t88) * t60);
U = t1;
