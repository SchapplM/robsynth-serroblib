% Calculate potential energy for
% S5RRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRP1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP1_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP1_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:44:41
% EndTime: 2019-12-05 18:44:41
% DurationCPUTime: 0.19s
% Computational Cost: add. (132->57), mult. (110->64), div. (0->0), fcn. (86->8), ass. (0->26)
t64 = rSges(6,1) + pkin(4);
t53 = -pkin(7) - pkin(6);
t63 = rSges(3,3) + pkin(6);
t49 = sin(qJ(2));
t62 = t49 * pkin(2) + pkin(5);
t51 = cos(qJ(2));
t40 = t51 * pkin(2) + pkin(1);
t61 = rSges(4,3) - t53;
t47 = -pkin(8) + t53;
t60 = rSges(5,3) - t47;
t59 = rSges(6,3) + qJ(5) - t47;
t48 = qJ(2) + qJ(3);
t41 = sin(t48);
t58 = pkin(3) * t41 + t62;
t42 = cos(t48);
t35 = pkin(3) * t42 + t40;
t57 = rSges(3,1) * t51 - rSges(3,2) * t49 + pkin(1);
t56 = rSges(4,1) * t42 - rSges(4,2) * t41 + t40;
t44 = qJ(4) + t48;
t38 = sin(t44);
t39 = cos(t44);
t55 = rSges(5,1) * t39 - rSges(5,2) * t38 + t35;
t54 = -rSges(6,2) * t38 + t64 * t39 + t35;
t52 = cos(qJ(1));
t50 = sin(qJ(1));
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t52 * rSges(2,1) - t50 * rSges(2,2)) + g(2) * (t50 * rSges(2,1) + t52 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(3) * (t49 * rSges(3,1) + t51 * rSges(3,2) + pkin(5)) + (g(1) * t57 - g(2) * t63) * t52 + (g(1) * t63 + g(2) * t57) * t50) - m(4) * (g(3) * (t41 * rSges(4,1) + t42 * rSges(4,2) + t62) + (g(1) * t56 - g(2) * t61) * t52 + (g(1) * t61 + g(2) * t56) * t50) - m(5) * (g(3) * (t38 * rSges(5,1) + t39 * rSges(5,2) + t58) + (g(1) * t55 - g(2) * t60) * t52 + (g(1) * t60 + g(2) * t55) * t50) - m(6) * (g(3) * (t39 * rSges(6,2) + t64 * t38 + t58) + (g(1) * t54 - g(2) * t59) * t52 + (g(1) * t59 + g(2) * t54) * t50);
U = t1;
