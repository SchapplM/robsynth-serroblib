% Calculate potential energy for
% S5PRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRR6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR6_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR6_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:09:29
% EndTime: 2019-12-05 17:09:30
% DurationCPUTime: 0.40s
% Computational Cost: add. (151->79), mult. (154->99), div. (0->0), fcn. (138->10), ass. (0->32)
t80 = rSges(5,3) + pkin(7);
t79 = rSges(6,3) + pkin(8) + pkin(7);
t54 = sin(pkin(9));
t55 = cos(pkin(9));
t78 = g(1) * t55 + g(2) * t54;
t75 = rSges(3,3) + pkin(5);
t53 = qJ(2) + qJ(3);
t48 = sin(t53);
t73 = rSges(4,2) * t48;
t50 = cos(t53);
t72 = t54 * t50;
t56 = sin(qJ(4));
t71 = t54 * t56;
t58 = cos(qJ(4));
t70 = t54 * t58;
t69 = t55 * t50;
t68 = t55 * t56;
t67 = t55 * t58;
t59 = cos(qJ(2));
t46 = pkin(2) * t59 + pkin(1);
t61 = -pkin(6) - pkin(5);
t65 = t54 * t46 + t55 * t61;
t57 = sin(qJ(2));
t64 = t57 * pkin(2) + qJ(1);
t43 = t55 * t46;
t63 = -t54 * t61 + t43;
t62 = rSges(3,1) * t59 - rSges(3,2) * t57 + pkin(1);
t52 = qJ(4) + qJ(5);
t49 = cos(t52);
t47 = sin(t52);
t45 = pkin(4) * t58 + pkin(3);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t55 - rSges(2,2) * t54) + g(2) * (rSges(2,1) * t54 + rSges(2,2) * t55) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(3) * (rSges(3,1) * t57 + rSges(3,2) * t59 + qJ(1)) + (g(1) * t62 - g(2) * t75) * t55 + (g(1) * t75 + g(2) * t62) * t54) - m(4) * (g(1) * (rSges(4,1) * t69 - t55 * t73 + t43) + g(2) * (-rSges(4,3) * t55 + t65) + g(3) * (rSges(4,1) * t48 + rSges(4,2) * t50 + t64) + (g(1) * (rSges(4,3) - t61) + g(2) * (rSges(4,1) * t50 - t73)) * t54) - m(5) * (g(1) * (pkin(3) * t69 + (t50 * t67 + t71) * rSges(5,1) + (-t50 * t68 + t70) * rSges(5,2) + t63) + g(2) * (pkin(3) * t72 + (t50 * t70 - t68) * rSges(5,1) + (-t50 * t71 - t67) * rSges(5,2) + t65) + g(3) * (-t80 * t50 + t64) + (g(3) * (rSges(5,1) * t58 - rSges(5,2) * t56 + pkin(3)) + t78 * t80) * t48) - m(6) * (g(1) * (t45 * t69 + pkin(4) * t71 + (t47 * t54 + t49 * t69) * rSges(6,1) + (-t47 * t69 + t49 * t54) * rSges(6,2) + t63) + g(2) * (t45 * t72 - pkin(4) * t68 + (-t47 * t55 + t49 * t72) * rSges(6,1) + (-t47 * t72 - t49 * t55) * rSges(6,2) + t65) + g(3) * (-t79 * t50 + t64) + (g(3) * (rSges(6,1) * t49 - rSges(6,2) * t47 + t45) + t78 * t79) * t48);
U = t1;
