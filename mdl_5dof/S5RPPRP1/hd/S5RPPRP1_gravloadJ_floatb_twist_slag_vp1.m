% Calculate Gravitation load on the joints for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:35:36
% EndTime: 2019-12-05 17:35:38
% DurationCPUTime: 0.35s
% Computational Cost: add. (182->71), mult. (194->99), div. (0->0), fcn. (179->8), ass. (0->26)
t32 = rSges(6,1) + pkin(4);
t13 = sin(pkin(8));
t14 = cos(pkin(8));
t31 = -pkin(3) * t14 - pkin(2) + (-rSges(5,3) - pkin(6)) * t13;
t17 = sin(qJ(1));
t30 = pkin(1) * t17;
t19 = cos(qJ(1));
t29 = pkin(1) * t19;
t16 = sin(qJ(4));
t28 = pkin(4) * t16;
t12 = qJ(1) + pkin(7);
t11 = cos(t12);
t27 = g(2) * t11;
t26 = t14 * t16;
t18 = cos(qJ(4));
t25 = t14 * t18;
t24 = -m(4) - m(5) - m(6);
t23 = t11 * qJ(3) - t30;
t21 = -rSges(4,1) * t14 + rSges(4,2) * t13 - pkin(2);
t10 = sin(t12);
t3 = -t10 * t18 + t11 * t26;
t1 = t10 * t26 + t11 * t18;
t20 = -t14 * (pkin(4) * t18 + pkin(3)) - pkin(2) + (-rSges(6,3) - qJ(5) - pkin(6)) * t13;
t4 = -t10 * t16 - t11 * t25;
t2 = t10 * t25 - t11 * t16;
t5 = [-m(2) * (g(2) * (-rSges(2,1) * t19 + t17 * rSges(2,2)) + g(3) * (-t17 * rSges(2,1) - rSges(2,2) * t19)) - m(3) * (g(2) * (-t11 * rSges(3,1) + t10 * rSges(3,2) - t29) + g(3) * (-rSges(3,1) * t10 - rSges(3,2) * t11 - t30)) - m(4) * (-g(2) * t29 + g(3) * t23 + (g(3) * rSges(4,3) + g(2) * t21) * t11 + (g(2) * (-rSges(4,3) - qJ(3)) + g(3) * t21) * t10) - m(5) * (g(2) * (t4 * rSges(5,1) + t3 * rSges(5,2) - t29) + g(3) * (-rSges(5,1) * t2 + rSges(5,2) * t1 + t23) + t31 * t27 + (-g(2) * qJ(3) + g(3) * t31) * t10) - m(6) * (g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) - t29) + g(3) * (-rSges(6,1) * t2 + rSges(6,2) * t1 + t23) + (g(2) * t20 + g(3) * t28) * t11 + (g(2) * (-qJ(3) - t28) + g(3) * t20) * t10), (-m(3) + t24) * g(1), t24 * (g(3) * t10 + t27), -m(5) * (g(2) * (rSges(5,1) * t1 + rSges(5,2) * t2) + g(3) * (-rSges(5,1) * t3 + rSges(5,2) * t4)) - m(6) * (g(2) * (rSges(6,2) * t2 + t32 * t1) + g(3) * (rSges(6,2) * t4 - t32 * t3)) + (-m(5) * (-rSges(5,1) * t16 - rSges(5,2) * t18) - m(6) * (-rSges(6,1) * t16 - rSges(6,2) * t18 - t28)) * g(1) * t13, -m(6) * (-g(1) * t14 + (-g(2) * t10 + g(3) * t11) * t13)];
taug = t5(:);
