% Calculate Gravitation load on the joints for
% S5RPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR6_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:57:41
% EndTime: 2019-12-31 17:57:42
% DurationCPUTime: 0.30s
% Computational Cost: add. (210->66), mult. (181->91), div. (0->0), fcn. (155->10), ass. (0->35)
t40 = rSges(6,3) + pkin(7);
t18 = sin(qJ(5));
t20 = cos(qJ(5));
t39 = m(6) * (rSges(6,1) * t20 - rSges(6,2) * t18 + pkin(4)) + m(5) * rSges(5,1);
t19 = sin(qJ(1));
t37 = pkin(1) * t19;
t14 = qJ(1) + pkin(8);
t9 = sin(t14);
t36 = t18 * t9;
t35 = t20 * t9;
t11 = cos(t14);
t21 = cos(qJ(1));
t12 = t21 * pkin(1);
t16 = cos(pkin(9));
t7 = pkin(3) * t16 + pkin(2);
t34 = t11 * t7 + t12;
t33 = t11 * t18;
t32 = t11 * t20;
t17 = -pkin(6) - qJ(3);
t31 = rSges(5,3) - t17;
t30 = rSges(4,3) + qJ(3);
t29 = -m(4) - m(5) - m(6);
t28 = g(1) * t37;
t13 = pkin(9) + qJ(4);
t10 = cos(t13);
t8 = sin(t13);
t27 = rSges(5,1) * t10 - t8 * rSges(5,2);
t26 = rSges(4,1) * t16 - rSges(4,2) * sin(pkin(9)) + pkin(2);
t24 = m(5) * rSges(5,2) - m(6) * t40;
t23 = pkin(4) * t10 + t40 * t8;
t4 = t10 * t32 + t36;
t3 = -t10 * t33 + t35;
t2 = -t10 * t35 + t33;
t1 = t10 * t36 + t32;
t5 = [-m(2) * (g(1) * (-t19 * rSges(2,1) - rSges(2,2) * t21) + g(2) * (rSges(2,1) * t21 - t19 * rSges(2,2))) - m(3) * (g(1) * (-rSges(3,1) * t9 - rSges(3,2) * t11 - t37) + g(2) * (rSges(3,1) * t11 - rSges(3,2) * t9 + t12)) - m(4) * (-t28 + g(2) * t12 + (-g(1) * t26 + g(2) * t30) * t9 + (g(1) * t30 + g(2) * t26) * t11) - m(5) * (-t28 + g(2) * t34 + (g(1) * t31 + g(2) * t27) * t11 + (g(1) * (-t27 - t7) + g(2) * t31) * t9) - m(6) * (g(1) * (rSges(6,1) * t2 + rSges(6,2) * t1 - t37) + g(2) * (rSges(6,1) * t4 + rSges(6,2) * t3 + t34) + (-g(1) * t17 + g(2) * t23) * t11 + (g(1) * (-t23 - t7) - g(2) * t17) * t9), (-m(3) + t29) * g(3), t29 * (g(1) * t9 - g(2) * t11), (-t39 * t10 + t24 * t8) * g(3) + (g(1) * t11 + g(2) * t9) * (t24 * t10 + t39 * t8), -m(6) * (g(1) * (rSges(6,1) * t3 - rSges(6,2) * t4) + g(2) * (-rSges(6,1) * t1 + rSges(6,2) * t2) + g(3) * (-rSges(6,1) * t18 - rSges(6,2) * t20) * t8)];
taug = t5(:);
