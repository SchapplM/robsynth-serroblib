% Calculate potential energy for
% S5PRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRP3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP3_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP3_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:32:34
% EndTime: 2019-12-05 15:32:35
% DurationCPUTime: 0.34s
% Computational Cost: add. (141->74), mult. (154->91), div. (0->0), fcn. (138->8), ass. (0->33)
t82 = rSges(5,3) + pkin(6);
t81 = rSges(6,3) + qJ(5) + pkin(6);
t56 = sin(pkin(7));
t57 = cos(pkin(7));
t80 = g(1) * t57 + g(2) * t56;
t77 = rSges(3,3) + pkin(5);
t55 = qJ(2) + pkin(8);
t52 = sin(t55);
t75 = rSges(4,2) * t52;
t53 = cos(t55);
t74 = t53 * t56;
t73 = t53 * t57;
t60 = sin(qJ(4));
t72 = t56 * t60;
t62 = cos(qJ(4));
t71 = t56 * t62;
t70 = t57 * t60;
t69 = t57 * t62;
t63 = cos(qJ(2));
t51 = pkin(2) * t63 + pkin(1);
t59 = -qJ(3) - pkin(5);
t67 = t56 * t51 + t57 * t59;
t61 = sin(qJ(2));
t66 = t61 * pkin(2) + qJ(1);
t48 = t57 * t51;
t65 = -t56 * t59 + t48;
t64 = rSges(3,1) * t63 - rSges(3,2) * t61 + pkin(1);
t50 = pkin(4) * t62 + pkin(3);
t46 = t53 * t69 + t72;
t45 = -t53 * t70 + t71;
t44 = t53 * t71 - t70;
t43 = -t53 * t72 - t69;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t57 - rSges(2,2) * t56) + g(2) * (rSges(2,1) * t56 + rSges(2,2) * t57) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(3) * (t61 * rSges(3,1) + rSges(3,2) * t63 + qJ(1)) + (g(1) * t64 - g(2) * t77) * t57 + (g(1) * t77 + g(2) * t64) * t56) - m(4) * (g(1) * (rSges(4,1) * t73 - t57 * t75 + t48) + g(2) * (-rSges(4,3) * t57 + t67) + g(3) * (rSges(4,1) * t52 + rSges(4,2) * t53 + t66) + (g(1) * (rSges(4,3) - t59) + g(2) * (rSges(4,1) * t53 - t75)) * t56) - m(5) * (g(1) * (rSges(5,1) * t46 + rSges(5,2) * t45 + pkin(3) * t73 + t65) + g(2) * (rSges(5,1) * t44 + rSges(5,2) * t43 + pkin(3) * t74 + t67) + g(3) * (-t82 * t53 + t66) + (g(3) * (rSges(5,1) * t62 - rSges(5,2) * t60 + pkin(3)) + t80 * t82) * t52) - m(6) * (g(1) * (rSges(6,1) * t46 + rSges(6,2) * t45 + pkin(4) * t72 + t50 * t73 + t65) + g(2) * (rSges(6,1) * t44 + rSges(6,2) * t43 - pkin(4) * t70 + t50 * t74 + t67) + g(3) * (-t81 * t53 + t66) + (g(3) * (rSges(6,1) * t62 - rSges(6,2) * t60 + t50) + t80 * t81) * t52);
U = t1;
