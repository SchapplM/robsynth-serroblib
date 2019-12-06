% Calculate potential energy for
% S5RRRRP2
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
% Datum: 2019-12-05 18:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRP2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP2_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:47:29
% EndTime: 2019-12-05 18:47:29
% DurationCPUTime: 0.22s
% Computational Cost: add. (129->57), mult. (97->62), div. (0->0), fcn. (73->8), ass. (0->25)
t58 = rSges(6,1) + pkin(4);
t57 = pkin(5) + pkin(6);
t47 = -pkin(8) - pkin(7);
t44 = sin(qJ(1));
t56 = t44 * pkin(1);
t55 = rSges(4,3) + pkin(7);
t45 = cos(qJ(3));
t32 = t45 * pkin(3) + pkin(2);
t54 = rSges(5,3) - t47;
t53 = rSges(6,3) + qJ(5) - t47;
t43 = sin(qJ(3));
t52 = t43 * pkin(3) + t57;
t46 = cos(qJ(1));
t39 = t46 * pkin(1);
t51 = -g(2) * t56 + g(3) * t39;
t50 = rSges(4,1) * t45 - rSges(4,2) * t43 + pkin(2);
t41 = qJ(3) + qJ(4);
t33 = sin(t41);
t35 = cos(t41);
t49 = rSges(5,1) * t35 - rSges(5,2) * t33 + t32;
t48 = -rSges(6,2) * t33 + t58 * t35 + t32;
t42 = qJ(1) + qJ(2);
t36 = cos(t42);
t34 = sin(t42);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (pkin(5) + rSges(2,3)) + g(2) * (-t44 * rSges(2,1) - t46 * rSges(2,2)) + g(3) * (t46 * rSges(2,1) - t44 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) + t57) + g(2) * (-t34 * rSges(3,1) - t36 * rSges(3,2) - t56) + g(3) * (t36 * rSges(3,1) - t34 * rSges(3,2) + t39)) - m(4) * (g(1) * (t43 * rSges(4,1) + t45 * rSges(4,2) + t57) + (g(2) * t55 + g(3) * t50) * t36 + (-g(2) * t50 + g(3) * t55) * t34 + t51) - m(5) * (g(1) * (t33 * rSges(5,1) + t35 * rSges(5,2) + t52) + (g(2) * t54 + g(3) * t49) * t36 + (-g(2) * t49 + g(3) * t54) * t34 + t51) - m(6) * (g(1) * (t35 * rSges(6,2) + t58 * t33 + t52) + (g(2) * t53 + g(3) * t48) * t36 + (-g(2) * t48 + g(3) * t53) * t34 + t51);
U = t1;
