% Calculate Gravitation load on the joints for
% S5RPPRR3
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
% Datum: 2020-01-03 11:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:27:27
% EndTime: 2020-01-03 11:27:30
% DurationCPUTime: 0.27s
% Computational Cost: add. (186->55), mult. (139->73), div. (0->0), fcn. (105->10), ass. (0->33)
t19 = pkin(9) + qJ(4);
t14 = qJ(5) + t19;
t7 = sin(t14);
t8 = cos(t14);
t45 = rSges(6,1) * t7 + rSges(6,2) * t8;
t12 = cos(t19);
t33 = t8 * rSges(6,1) - rSges(6,2) * t7;
t44 = pkin(4) * t12 + t33;
t20 = qJ(1) + pkin(8);
t13 = cos(t20);
t43 = t45 * t13;
t10 = sin(t19);
t39 = pkin(4) * t10;
t11 = sin(t20);
t38 = g(2) * t11;
t22 = cos(pkin(9));
t9 = t22 * pkin(3) + pkin(2);
t23 = -pkin(6) - qJ(3);
t37 = rSges(5,3) - t23;
t36 = rSges(6,3) + pkin(7) - t23;
t35 = rSges(4,3) + qJ(3);
t34 = -m(4) - m(5) - m(6);
t24 = sin(qJ(1));
t16 = t24 * pkin(1);
t25 = cos(qJ(1));
t17 = t25 * pkin(1);
t31 = g(2) * t17 + g(3) * t16;
t30 = rSges(5,1) * t12 - rSges(5,2) * t10;
t29 = rSges(5,1) * t10 + rSges(5,2) * t12;
t28 = t9 + t44;
t27 = rSges(4,1) * t22 - rSges(4,2) * sin(pkin(9)) + pkin(2);
t26 = t30 + t9;
t1 = [-m(2) * (g(2) * (rSges(2,1) * t25 - t24 * rSges(2,2)) + g(3) * (t24 * rSges(2,1) + rSges(2,2) * t25)) - m(3) * (g(2) * (rSges(3,1) * t13 - rSges(3,2) * t11 + t17) + g(3) * (rSges(3,1) * t11 + rSges(3,2) * t13 + t16)) - m(4) * ((g(2) * t27 - g(3) * t35) * t13 + (g(2) * t35 + g(3) * t27) * t11 + t31) - m(5) * ((g(2) * t26 - g(3) * t37) * t13 + (g(2) * t37 + g(3) * t26) * t11 + t31) - m(6) * ((g(2) * t28 - g(3) * t36) * t13 + (g(2) * t36 + g(3) * t28) * t11 + t31), (-m(3) + t34) * g(1), t34 * (-g(2) * t13 - g(3) * t11), -m(5) * (g(3) * t29 * t13 + g(1) * t30) - m(6) * (g(1) * t44 + g(3) * (t13 * t39 + t43)) + (m(5) * t29 - m(6) * (-t45 - t39)) * t38, -m(6) * (g(1) * t33 + g(3) * t43 - t38 * t45)];
taug = t1(:);
