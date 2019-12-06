% Calculate potential energy for
% S5RRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
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
% Datum: 2019-12-05 18:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_energypot_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR1_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR1_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:49:42
% EndTime: 2019-12-05 18:49:42
% DurationCPUTime: 0.26s
% Computational Cost: add. (132->62), mult. (122->81), div. (0->0), fcn. (102->10), ass. (0->28)
t62 = rSges(6,3) + pkin(6);
t44 = qJ(2) + qJ(3);
t42 = qJ(4) + t44;
t38 = cos(t42);
t61 = t38 * pkin(4);
t49 = cos(qJ(2));
t39 = t49 * pkin(2) + pkin(1);
t45 = sin(qJ(5));
t47 = sin(qJ(1));
t59 = t47 * t45;
t48 = cos(qJ(5));
t58 = t47 * t48;
t50 = cos(qJ(1));
t57 = t50 * t45;
t56 = t50 * t48;
t46 = sin(qJ(2));
t55 = -t46 * pkin(2) + pkin(5);
t37 = sin(t42);
t54 = rSges(5,1) * t38 - rSges(5,2) * t37;
t40 = sin(t44);
t53 = -pkin(3) * t40 + t55;
t52 = rSges(3,1) * t49 - rSges(3,2) * t46 + pkin(1);
t41 = cos(t44);
t51 = rSges(4,1) * t41 - rSges(4,2) * t40 + t39;
t36 = pkin(3) * t41 + t39;
t35 = t50 * t36;
t34 = t47 * t36;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t50 * rSges(2,1) - t47 * rSges(2,2)) + g(2) * (t47 * rSges(2,1) + t50 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(3) * (-t46 * rSges(3,1) - t49 * rSges(3,2) + pkin(5)) + (g(2) * rSges(3,3) + g(1) * t52) * t50 + (-g(1) * rSges(3,3) + g(2) * t52) * t47) - m(4) * (g(3) * (-t40 * rSges(4,1) - t41 * rSges(4,2) + t55) + (g(2) * rSges(4,3) + g(1) * t51) * t50 + (-g(1) * rSges(4,3) + g(2) * t51) * t47) - m(5) * (g(1) * (-t47 * rSges(5,3) + t54 * t50 + t35) + g(2) * (t50 * rSges(5,3) + t54 * t47 + t34) + g(3) * (-t37 * rSges(5,1) - t38 * rSges(5,2) + t53)) - m(6) * (g(1) * (t50 * t61 + t35 + (t38 * t56 - t59) * rSges(6,1) + (-t38 * t57 - t58) * rSges(6,2)) + g(2) * (t47 * t61 + t34 + (t38 * t58 + t57) * rSges(6,1) + (-t38 * t59 + t56) * rSges(6,2)) + g(3) * (t62 * t38 + t53) + (g(3) * (-rSges(6,1) * t48 + rSges(6,2) * t45 - pkin(4)) + (g(1) * t50 + g(2) * t47) * t62) * t37);
U = t1;
