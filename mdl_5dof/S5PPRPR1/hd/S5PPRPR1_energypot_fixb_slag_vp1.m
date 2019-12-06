% Calculate potential energy for
% S5PPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPRPR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR1_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR1_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:00:48
% EndTime: 2019-12-05 15:00:49
% DurationCPUTime: 0.39s
% Computational Cost: add. (151->79), mult. (154->99), div. (0->0), fcn. (138->10), ass. (0->31)
t79 = rSges(6,3) + pkin(6) + qJ(4);
t56 = sin(pkin(7));
t59 = cos(pkin(7));
t78 = g(1) * t59 + g(2) * t56;
t77 = rSges(5,3) + qJ(4);
t53 = pkin(8) + qJ(3);
t48 = sin(t53);
t74 = rSges(4,2) * t48;
t54 = sin(pkin(9));
t73 = t54 * t56;
t72 = t54 * t59;
t50 = cos(t53);
t71 = t56 * t50;
t57 = cos(pkin(9));
t70 = t56 * t57;
t69 = t59 * t50;
t58 = cos(pkin(8));
t46 = pkin(2) * t58 + pkin(1);
t61 = -pkin(5) - qJ(2);
t67 = t56 * t46 + t59 * t61;
t66 = rSges(3,3) + qJ(2);
t55 = sin(pkin(8));
t64 = t55 * pkin(2) + qJ(1);
t43 = t59 * t46;
t63 = -t56 * t61 + t43;
t62 = rSges(3,1) * t58 - rSges(3,2) * t55 + pkin(1);
t52 = pkin(9) + qJ(5);
t49 = cos(t52);
t47 = sin(t52);
t45 = pkin(4) * t57 + pkin(3);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t59 - rSges(2,2) * t56) + g(2) * (rSges(2,1) * t56 + rSges(2,2) * t59) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(3) * (rSges(3,1) * t55 + rSges(3,2) * t58 + qJ(1)) + (g(1) * t62 - g(2) * t66) * t59 + (g(1) * t66 + g(2) * t62) * t56) - m(4) * (g(1) * (rSges(4,1) * t69 - t59 * t74 + t43) + g(2) * (-rSges(4,3) * t59 + t67) + g(3) * (rSges(4,1) * t48 + rSges(4,2) * t50 + t64) + (g(1) * (rSges(4,3) - t61) + g(2) * (rSges(4,1) * t50 - t74)) * t56) - m(5) * (g(1) * (pkin(3) * t69 + (t57 * t69 + t73) * rSges(5,1) + (-t54 * t69 + t70) * rSges(5,2) + t63) + g(2) * (pkin(3) * t71 + (t50 * t70 - t72) * rSges(5,1) + (-t54 * t71 - t57 * t59) * rSges(5,2) + t67) + g(3) * (-t77 * t50 + t64) + (g(3) * (rSges(5,1) * t57 - rSges(5,2) * t54 + pkin(3)) + t78 * t77) * t48) - m(6) * (g(1) * (t45 * t69 + pkin(4) * t73 + (t47 * t56 + t49 * t69) * rSges(6,1) + (-t47 * t69 + t49 * t56) * rSges(6,2) + t63) + g(2) * (t45 * t71 - pkin(4) * t72 + (-t47 * t59 + t49 * t71) * rSges(6,1) + (-t47 * t71 - t49 * t59) * rSges(6,2) + t67) + g(3) * (-t79 * t50 + t64) + (g(3) * (rSges(6,1) * t49 - rSges(6,2) * t47 + t45) + t78 * t79) * t48);
U = t1;
