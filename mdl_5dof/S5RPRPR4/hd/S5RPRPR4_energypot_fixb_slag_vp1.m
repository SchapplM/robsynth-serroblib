% Calculate potential energy for
% S5RPRPR4
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
% Datum: 2019-12-05 17:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR4_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:53:05
% EndTime: 2019-12-05 17:53:05
% DurationCPUTime: 0.19s
% Computational Cost: add. (135->59), mult. (97->64), div. (0->0), fcn. (73->10), ass. (0->27)
t50 = sin(qJ(1));
t62 = pkin(1) * t50;
t61 = rSges(4,3) + pkin(6);
t51 = cos(qJ(3));
t36 = t51 * pkin(3) + pkin(2);
t48 = -qJ(4) - pkin(6);
t60 = rSges(5,3) - t48;
t59 = rSges(6,3) + pkin(7) - t48;
t58 = pkin(5) + qJ(2);
t46 = qJ(3) + pkin(9);
t49 = sin(qJ(3));
t57 = t49 * pkin(3) + t58;
t52 = cos(qJ(1));
t44 = t52 * pkin(1);
t56 = -g(2) * t62 + g(3) * t44;
t55 = rSges(4,1) * t51 - rSges(4,2) * t49 + pkin(2);
t37 = sin(t46);
t39 = cos(t46);
t54 = rSges(5,1) * t39 - rSges(5,2) * t37 + t36;
t41 = qJ(5) + t46;
t34 = sin(t41);
t35 = cos(t41);
t53 = rSges(6,1) * t35 - rSges(6,2) * t34 + pkin(4) * t39 + t36;
t47 = qJ(1) + pkin(8);
t40 = cos(t47);
t38 = sin(t47);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (pkin(5) + rSges(2,3)) + g(2) * (-t50 * rSges(2,1) - rSges(2,2) * t52) + g(3) * (rSges(2,1) * t52 - t50 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) + t58) + g(2) * (-rSges(3,1) * t38 - rSges(3,2) * t40 - t62) + g(3) * (rSges(3,1) * t40 - rSges(3,2) * t38 + t44)) - m(4) * (g(1) * (rSges(4,1) * t49 + rSges(4,2) * t51 + t58) + (g(2) * t61 + g(3) * t55) * t40 + (-g(2) * t55 + g(3) * t61) * t38 + t56) - m(5) * (g(1) * (rSges(5,1) * t37 + rSges(5,2) * t39 + t57) + (g(2) * t60 + g(3) * t54) * t40 + (-g(2) * t54 + g(3) * t60) * t38 + t56) - m(6) * (g(1) * (rSges(6,1) * t34 + rSges(6,2) * t35 + pkin(4) * t37 + t57) + (g(2) * t59 + g(3) * t53) * t40 + (-g(2) * t53 + g(3) * t59) * t38 + t56);
U = t1;
