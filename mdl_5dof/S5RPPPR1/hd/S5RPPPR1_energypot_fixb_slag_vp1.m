% Calculate potential energy for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPPR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR1_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:28:42
% EndTime: 2019-12-05 17:28:42
% DurationCPUTime: 0.40s
% Computational Cost: add. (157->77), mult. (141->98), div. (0->0), fcn. (125->10), ass. (0->31)
t79 = rSges(6,3) + pkin(6) + qJ(4);
t52 = qJ(1) + pkin(7);
t47 = sin(t52);
t49 = cos(t52);
t78 = -g(2) * t47 + g(3) * t49;
t77 = rSges(5,3) + qJ(4);
t58 = sin(qJ(1));
t74 = t58 * pkin(1);
t54 = sin(pkin(8));
t73 = rSges(4,2) * t54;
t53 = sin(pkin(9));
t72 = t47 * t53;
t56 = cos(pkin(8));
t71 = t47 * t56;
t70 = t49 * t53;
t69 = t49 * t56;
t68 = t53 * t56;
t55 = cos(pkin(9));
t67 = t55 * t56;
t45 = pkin(4) * t55 + pkin(3);
t66 = t56 * t45;
t64 = pkin(5) + qJ(2);
t59 = cos(qJ(1));
t50 = t59 * pkin(1);
t62 = t49 * pkin(2) + t47 * qJ(3) + t50;
t61 = t49 * qJ(3) - t74;
t60 = -t47 * pkin(2) + t61;
t51 = pkin(9) + qJ(5);
t48 = cos(t51);
t46 = sin(t51);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (pkin(5) + rSges(2,3)) + g(2) * (-t58 * rSges(2,1) - rSges(2,2) * t59) + g(3) * (rSges(2,1) * t59 - t58 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) + t64) + g(2) * (-rSges(3,1) * t47 - rSges(3,2) * t49 - t74) + g(3) * (rSges(3,1) * t49 - rSges(3,2) * t47 + t50)) - m(4) * (g(1) * (rSges(4,1) * t54 + rSges(4,2) * t56 + t64) + g(2) * (rSges(4,3) * t49 + t61) + g(3) * (rSges(4,1) * t69 - t49 * t73 + t62) + (g(2) * (-rSges(4,1) * t56 - pkin(2) + t73) + g(3) * rSges(4,3)) * t47) - m(5) * (g(1) * (-t77 * t56 + t64) + g(2) * (-pkin(3) * t71 + (-t47 * t67 + t70) * rSges(5,1) + (t47 * t68 + t49 * t55) * rSges(5,2) + t60) + g(3) * (pkin(3) * t69 + (t49 * t67 + t72) * rSges(5,1) + (t47 * t55 - t49 * t68) * rSges(5,2) + t62) + (g(1) * (rSges(5,1) * t55 - rSges(5,2) * t53 + pkin(3)) + t78 * t77) * t54) - m(6) * (g(1) * (-t79 * t56 + t64) + g(2) * (-t47 * t66 + pkin(4) * t70 + (t46 * t49 - t48 * t71) * rSges(6,1) + (t46 * t71 + t48 * t49) * rSges(6,2) + t60) + g(3) * (t49 * t66 + pkin(4) * t72 + (t46 * t47 + t48 * t69) * rSges(6,1) + (-t46 * t69 + t47 * t48) * rSges(6,2) + t62) + (g(1) * (rSges(6,1) * t48 - rSges(6,2) * t46 + t45) + t78 * t79) * t54);
U = t1;
