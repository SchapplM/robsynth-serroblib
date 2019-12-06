% Calculate potential energy for
% S5RPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_energypot_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR1_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR1_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:37:56
% EndTime: 2019-12-05 17:37:56
% DurationCPUTime: 0.14s
% Computational Cost: add. (80->57), mult. (91->61), div. (0->0), fcn. (67->6), ass. (0->19)
t48 = pkin(2) + pkin(5);
t47 = -rSges(5,3) - pkin(6);
t46 = -rSges(6,3) - pkin(7) - pkin(6);
t35 = sin(qJ(1));
t31 = t35 * pkin(1);
t45 = t35 * qJ(3) + t31;
t37 = cos(qJ(1));
t44 = t37 * pkin(1) + t35 * qJ(2);
t43 = pkin(3) + t48;
t42 = t37 * qJ(3) + t44;
t34 = sin(qJ(4));
t36 = cos(qJ(4));
t41 = rSges(5,1) * t34 + rSges(5,2) * t36;
t33 = qJ(4) + qJ(5);
t26 = sin(t33);
t27 = cos(t33);
t40 = rSges(6,1) * t26 + rSges(6,2) * t27 + pkin(4) * t34;
t39 = g(1) * t42 + g(2) * t45;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t37 * rSges(2,1) - t35 * rSges(2,2)) + g(2) * (t35 * rSges(2,1) + t37 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (-t37 * rSges(3,2) + t35 * rSges(3,3) + t44) + g(2) * (-t35 * rSges(3,2) + t31 + (-rSges(3,3) - qJ(2)) * t37) + g(3) * (pkin(5) + rSges(3,1))) - m(4) * (g(1) * (t35 * rSges(4,2) + t37 * rSges(4,3) + t42) + g(2) * (t35 * rSges(4,3) + (-rSges(4,2) - qJ(2)) * t37 + t45) + g(3) * (rSges(4,1) + t48)) - m(5) * (g(3) * (t36 * rSges(5,1) - t34 * rSges(5,2) + t43) + (g(1) * t47 + g(2) * t41) * t35 + (g(1) * t41 + g(2) * (-qJ(2) - t47)) * t37 + t39) - m(6) * (g(3) * (t27 * rSges(6,1) - t26 * rSges(6,2) + t36 * pkin(4) + t43) + (g(1) * t46 + g(2) * t40) * t35 + (g(1) * t40 + g(2) * (-qJ(2) - t46)) * t37 + t39);
U = t1;
