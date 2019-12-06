% Calculate Gravitation load on the joints for
% S5PRRRR4
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
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR4_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:40
% EndTime: 2019-12-05 17:07:41
% DurationCPUTime: 0.22s
% Computational Cost: add. (242->49), mult. (151->61), div. (0->0), fcn. (113->8), ass. (0->29)
t28 = cos(qJ(4));
t26 = qJ(4) + qJ(5);
t22 = sin(t26);
t23 = cos(t26);
t37 = t23 * rSges(6,1) - rSges(6,2) * t22;
t51 = t28 * pkin(4) + t37;
t50 = pkin(7) + rSges(5,3);
t49 = pkin(8) + pkin(7) + rSges(6,3);
t27 = sin(qJ(4));
t48 = rSges(5,1) * t28 - rSges(5,2) * t27;
t47 = -pkin(3) - t48;
t46 = -pkin(3) - t51;
t25 = pkin(9) + qJ(2);
t21 = qJ(3) + t25;
t16 = sin(t21);
t17 = cos(t21);
t45 = g(1) * t17 + g(2) * t16;
t19 = sin(t25);
t44 = pkin(2) * t19;
t38 = t17 * rSges(4,1) - rSges(4,2) * t16;
t36 = -rSges(4,1) * t16 - rSges(4,2) * t17;
t35 = -rSges(6,1) * t22 - rSges(6,2) * t23;
t34 = t50 * t16 - t47 * t17;
t33 = t47 * t16 + t50 * t17;
t32 = t49 * t16 - t46 * t17;
t31 = t46 * t16 + t49 * t17;
t20 = cos(t25);
t14 = pkin(2) * t20;
t1 = [(-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t19 - rSges(3,2) * t20) + g(2) * (rSges(3,1) * t20 - rSges(3,2) * t19)) - m(4) * (g(1) * (t36 - t44) + g(2) * (t14 + t38)) - m(5) * (g(1) * (t33 - t44) + g(2) * (t14 + t34)) - m(6) * (g(1) * (t31 - t44) + g(2) * (t14 + t32)), -m(4) * (g(1) * t36 + g(2) * t38) - m(5) * (g(1) * t33 + g(2) * t34) - m(6) * (g(1) * t31 + g(2) * t32), (-m(5) * t48 - m(6) * t51) * g(3) + t45 * (-m(5) * (-rSges(5,1) * t27 - rSges(5,2) * t28) - m(6) * (-pkin(4) * t27 + t35)), -m(6) * (g(3) * t37 + t45 * t35)];
taug = t1(:);
