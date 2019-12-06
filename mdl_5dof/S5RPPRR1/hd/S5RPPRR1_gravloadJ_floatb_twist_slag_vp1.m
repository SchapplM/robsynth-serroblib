% Calculate Gravitation load on the joints for
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
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR1_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:37:59
% EndTime: 2019-12-05 17:38:00
% DurationCPUTime: 0.21s
% Computational Cost: add. (101->56), mult. (145->75), div. (0->0), fcn. (111->6), ass. (0->25)
t13 = sin(qJ(1));
t15 = cos(qJ(1));
t33 = g(1) * t15 + g(2) * t13;
t11 = qJ(4) + qJ(5);
t6 = cos(t11);
t32 = rSges(6,1) * t6;
t5 = sin(t11);
t31 = rSges(6,2) * t5;
t28 = -rSges(5,3) - pkin(6);
t27 = t15 * pkin(1) + t13 * qJ(2);
t26 = -rSges(6,3) - pkin(7) - pkin(6);
t25 = -pkin(1) - qJ(3);
t24 = -m(4) - m(5) - m(6);
t23 = t15 * qJ(3) + t27;
t22 = -rSges(6,1) * t5 - rSges(6,2) * t6;
t14 = cos(qJ(4));
t21 = pkin(4) * t14 - t31;
t12 = sin(qJ(4));
t19 = rSges(5,1) * t12 + rSges(5,2) * t14;
t9 = t15 * qJ(2);
t18 = g(1) * t9 + g(2) * t23;
t17 = pkin(4) * t12 - t22;
t4 = t15 * t32;
t3 = t13 * t32;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t13 - rSges(2,2) * t15) + g(2) * (rSges(2,1) * t15 - rSges(2,2) * t13)) - m(3) * (g(1) * (rSges(3,3) * t15 + t9 + (rSges(3,2) - pkin(1)) * t13) + g(2) * (-rSges(3,2) * t15 + rSges(3,3) * t13 + t27)) - m(4) * (g(1) * (rSges(4,2) * t15 + t9) + g(2) * (rSges(4,3) * t15 + t23) + (g(1) * (-rSges(4,3) + t25) + g(2) * rSges(4,2)) * t13) - m(5) * ((g(1) * t28 + g(2) * t19) * t15 + (g(1) * (-t19 + t25) + g(2) * t28) * t13 + t18) - m(6) * ((g(1) * t26 + g(2) * t17) * t15 + (g(1) * (-t17 + t25) + g(2) * t26) * t13 + t18), (-m(3) + t24) * (g(1) * t13 - g(2) * t15), t24 * t33, -m(5) * (-g(3) * t19 + t33 * (rSges(5,1) * t14 - rSges(5,2) * t12)) - m(6) * (g(1) * (t21 * t15 + t4) + g(2) * (t21 * t13 + t3) - g(3) * t17), -m(6) * (g(1) * (-t15 * t31 + t4) + g(2) * (-t13 * t31 + t3) + g(3) * t22)];
taug = t1(:);
