% Calculate potential energy for
% S5RRPRR3
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
% Datum: 2019-12-05 18:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR3_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:30:16
% EndTime: 2019-12-05 18:30:16
% DurationCPUTime: 0.14s
% Computational Cost: add. (134->54), mult. (74->56), div. (0->0), fcn. (50->10), ass. (0->25)
t54 = pkin(5) + pkin(6);
t42 = sin(qJ(1));
t53 = t42 * pkin(1);
t52 = rSges(6,3) + pkin(8);
t40 = qJ(1) + qJ(2);
t38 = cos(t40);
t44 = cos(qJ(1));
t39 = t44 * pkin(1);
t51 = pkin(2) * t38 + t39;
t50 = qJ(3) + t54;
t36 = pkin(9) + t40;
t33 = cos(t36);
t49 = pkin(3) * t33 + t51;
t48 = pkin(7) + t50;
t37 = sin(t40);
t47 = -pkin(2) * t37 - t53;
t41 = sin(qJ(5));
t43 = cos(qJ(5));
t46 = rSges(6,1) * t43 - rSges(6,2) * t41 + pkin(4);
t32 = sin(t36);
t45 = -pkin(3) * t32 + t47;
t35 = qJ(4) + t36;
t31 = cos(t35);
t30 = sin(t35);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (pkin(5) + rSges(2,3)) + g(2) * (-t42 * rSges(2,1) - t44 * rSges(2,2)) + g(3) * (t44 * rSges(2,1) - t42 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) + t54) + g(2) * (-t37 * rSges(3,1) - t38 * rSges(3,2) - t53) + g(3) * (t38 * rSges(3,1) - t37 * rSges(3,2) + t39)) - m(4) * (g(1) * (rSges(4,3) + t50) + g(2) * (-t32 * rSges(4,1) - t33 * rSges(4,2) + t47) + g(3) * (t33 * rSges(4,1) - t32 * rSges(4,2) + t51)) - m(5) * (g(1) * (rSges(5,3) + t48) + g(2) * (-t30 * rSges(5,1) - t31 * rSges(5,2) + t45) + g(3) * (t31 * rSges(5,1) - t30 * rSges(5,2) + t49)) - m(6) * (g(1) * (t41 * rSges(6,1) + t43 * rSges(6,2) + t48) + g(2) * t45 + g(3) * t49 + (g(2) * t52 + g(3) * t46) * t31 + (-g(2) * t46 + g(3) * t52) * t30);
U = t1;
