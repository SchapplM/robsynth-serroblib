% Calculate potential energy for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR6_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:35:37
% EndTime: 2019-12-05 18:35:37
% DurationCPUTime: 0.38s
% Computational Cost: add. (157->77), mult. (141->98), div. (0->0), fcn. (125->10), ass. (0->31)
t79 = rSges(5,3) + pkin(7);
t78 = rSges(6,3) + pkin(8) + pkin(7);
t52 = qJ(1) + qJ(2);
t47 = sin(t52);
t49 = cos(t52);
t77 = -g(2) * t47 + g(3) * t49;
t76 = pkin(5) + pkin(6);
t56 = sin(qJ(1));
t73 = t56 * pkin(1);
t53 = sin(pkin(9));
t71 = rSges(4,2) * t53;
t54 = cos(pkin(9));
t70 = t47 * t54;
t55 = sin(qJ(4));
t69 = t47 * t55;
t68 = t49 * t54;
t67 = t49 * t55;
t57 = cos(qJ(4));
t45 = pkin(4) * t57 + pkin(3);
t66 = t54 * t45;
t65 = t54 * t55;
t64 = t54 * t57;
t58 = cos(qJ(1));
t50 = t58 * pkin(1);
t62 = t49 * pkin(2) + t47 * qJ(3) + t50;
t61 = t49 * qJ(3) - t73;
t60 = -t47 * pkin(2) + t61;
t51 = qJ(4) + qJ(5);
t48 = cos(t51);
t46 = sin(t51);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (pkin(5) + rSges(2,3)) + g(2) * (-rSges(2,1) * t56 - rSges(2,2) * t58) + g(3) * (rSges(2,1) * t58 - rSges(2,2) * t56)) - m(3) * (g(1) * (rSges(3,3) + t76) + g(2) * (-rSges(3,1) * t47 - rSges(3,2) * t49 - t73) + g(3) * (rSges(3,1) * t49 - rSges(3,2) * t47 + t50)) - m(4) * (g(1) * (rSges(4,1) * t53 + rSges(4,2) * t54 + t76) + g(2) * (rSges(4,3) * t49 + t61) + g(3) * (rSges(4,1) * t68 - t49 * t71 + t62) + (g(2) * (-rSges(4,1) * t54 - pkin(2) + t71) + g(3) * rSges(4,3)) * t47) - m(5) * (g(1) * (-t79 * t54 + t76) + g(2) * (-pkin(3) * t70 + (-t47 * t64 + t67) * rSges(5,1) + (t47 * t65 + t49 * t57) * rSges(5,2) + t60) + g(3) * (pkin(3) * t68 + (t49 * t64 + t69) * rSges(5,1) + (t47 * t57 - t49 * t65) * rSges(5,2) + t62) + (g(1) * (rSges(5,1) * t57 - rSges(5,2) * t55 + pkin(3)) + t77 * t79) * t53) - m(6) * (g(1) * (-t78 * t54 + t76) + g(2) * (-t47 * t66 + pkin(4) * t67 + (t46 * t49 - t48 * t70) * rSges(6,1) + (t46 * t70 + t48 * t49) * rSges(6,2) + t60) + g(3) * (t49 * t66 + pkin(4) * t69 + (t46 * t47 + t48 * t68) * rSges(6,1) + (-t46 * t68 + t47 * t48) * rSges(6,2) + t62) + (g(1) * (rSges(6,1) * t48 - rSges(6,2) * t46 + t45) + t77 * t78) * t53);
U = t1;
