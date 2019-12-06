% Calculate potential energy for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-05 18:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRP1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP1_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:21:55
% EndTime: 2019-12-05 18:21:55
% DurationCPUTime: 0.19s
% Computational Cost: add. (130->54), mult. (85->58), div. (0->0), fcn. (61->8), ass. (0->23)
t52 = rSges(6,1) + pkin(4);
t51 = pkin(5) + pkin(6);
t39 = sin(qJ(1));
t50 = t39 * pkin(1);
t49 = rSges(5,3) + pkin(7);
t48 = rSges(6,3) + qJ(5) + pkin(7);
t36 = qJ(1) + qJ(2);
t34 = cos(t36);
t41 = cos(qJ(1));
t35 = t41 * pkin(1);
t47 = pkin(2) * t34 + t35;
t46 = qJ(3) + t51;
t33 = sin(t36);
t45 = -pkin(2) * t33 - t50;
t38 = sin(qJ(4));
t40 = cos(qJ(4));
t44 = rSges(5,1) * t40 - rSges(5,2) * t38 + pkin(3);
t43 = -rSges(6,2) * t38 + t52 * t40 + pkin(3);
t42 = g(2) * t45 + g(3) * t47;
t32 = pkin(8) + t36;
t29 = cos(t32);
t28 = sin(t32);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (pkin(5) + rSges(2,3)) + g(2) * (-t39 * rSges(2,1) - t41 * rSges(2,2)) + g(3) * (t41 * rSges(2,1) - t39 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) + t51) + g(2) * (-t33 * rSges(3,1) - t34 * rSges(3,2) - t50) + g(3) * (t34 * rSges(3,1) - t33 * rSges(3,2) + t35)) - m(4) * (g(1) * (rSges(4,3) + t46) + g(2) * (-t28 * rSges(4,1) - t29 * rSges(4,2) + t45) + g(3) * (t29 * rSges(4,1) - t28 * rSges(4,2) + t47)) - m(5) * (g(1) * (t38 * rSges(5,1) + t40 * rSges(5,2) + t46) + (g(2) * t49 + g(3) * t44) * t29 + (-g(2) * t44 + g(3) * t49) * t28 + t42) - m(6) * (g(1) * (t40 * rSges(6,2) + t52 * t38 + t46) + (g(2) * t48 + g(3) * t43) * t29 + (-g(2) * t43 + g(3) * t48) * t28 + t42);
U = t1;
