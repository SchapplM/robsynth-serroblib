% Calculate Gravitation load on the joints for
% S5RPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP10_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP10_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP10_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP10_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:50:55
% EndTime: 2019-12-31 18:50:56
% DurationCPUTime: 0.41s
% Computational Cost: add. (230->82), mult. (288->114), div. (0->0), fcn. (264->8), ass. (0->35)
t43 = rSges(6,1) + pkin(4);
t42 = rSges(5,3) + pkin(7);
t41 = rSges(6,3) + qJ(5) + pkin(7);
t17 = sin(qJ(1));
t19 = cos(qJ(1));
t40 = g(1) * t19 + g(2) * t17;
t16 = sin(qJ(4));
t18 = cos(qJ(4));
t8 = t18 * pkin(4) + pkin(3);
t39 = m(5) * (rSges(5,1) * t18 - rSges(5,2) * t16 + pkin(3)) + m(6) * (rSges(6,1) * t18 - rSges(6,2) * t16 + t8) + m(4) * rSges(4,1);
t37 = pkin(4) * t16;
t34 = t17 * t16;
t33 = t17 * t18;
t32 = t19 * t16;
t31 = t19 * t18;
t15 = -pkin(6) - qJ(2);
t30 = rSges(4,3) - t15;
t29 = rSges(3,3) + qJ(2);
t28 = -t15 + t37;
t11 = pkin(8) + qJ(3);
t10 = cos(t11);
t9 = sin(t11);
t27 = t10 * rSges(4,1) - t9 * rSges(4,2);
t13 = cos(pkin(8));
t26 = rSges(3,1) * t13 - rSges(3,2) * sin(pkin(8)) + pkin(1);
t3 = -t10 * t32 + t33;
t1 = t10 * t34 + t31;
t23 = t10 * pkin(3) + t42 * t9;
t22 = t10 * t8 + t41 * t9;
t21 = m(4) * rSges(4,2) - m(5) * t42 - m(6) * t41;
t7 = t13 * pkin(2) + pkin(1);
t5 = t19 * t7;
t4 = t10 * t31 + t34;
t2 = -t10 * t33 + t32;
t6 = [-m(2) * (g(1) * (-t17 * rSges(2,1) - t19 * rSges(2,2)) + g(2) * (t19 * rSges(2,1) - t17 * rSges(2,2))) - m(3) * ((g(1) * t29 + g(2) * t26) * t19 + (-g(1) * t26 + g(2) * t29) * t17) - m(4) * (g(2) * t5 + (g(1) * t30 + g(2) * t27) * t19 + (g(1) * (-t27 - t7) + g(2) * t30) * t17) - m(5) * (g(1) * (t2 * rSges(5,1) + t1 * rSges(5,2)) + g(2) * (t4 * rSges(5,1) + t3 * rSges(5,2) + t5) + (-g(1) * t15 + g(2) * t23) * t19 + (g(1) * (-t23 - t7) - g(2) * t15) * t17) - m(6) * (g(1) * (t2 * rSges(6,1) + t1 * rSges(6,2)) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t5) + (g(1) * t28 + g(2) * t22) * t19 + (g(1) * (-t22 - t7) + g(2) * t28) * t17), (-m(3) - m(4) - m(5) - m(6)) * (g(1) * t17 - g(2) * t19), (-t10 * t39 + t21 * t9) * g(3) + t40 * (t21 * t10 + t39 * t9), -m(5) * (g(1) * (t3 * rSges(5,1) - t4 * rSges(5,2)) + g(2) * (-t1 * rSges(5,1) + t2 * rSges(5,2))) - m(6) * (g(1) * (-t4 * rSges(6,2) + t43 * t3) + g(2) * (t2 * rSges(6,2) - t43 * t1)) + (-m(5) * (-rSges(5,1) * t16 - rSges(5,2) * t18) - m(6) * (-rSges(6,1) * t16 - rSges(6,2) * t18 - t37)) * g(3) * t9, -m(6) * (-g(3) * t10 + t40 * t9)];
taug = t6(:);
