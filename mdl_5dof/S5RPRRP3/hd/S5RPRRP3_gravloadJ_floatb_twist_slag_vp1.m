% Calculate Gravitation load on the joints for
% S5RPRRP3
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
% Datum: 2020-01-03 11:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:46:54
% EndTime: 2020-01-03 11:46:56
% DurationCPUTime: 0.30s
% Computational Cost: add. (204->71), mult. (174->91), div. (0->0), fcn. (134->8), ass. (0->37)
t20 = qJ(3) + qJ(4);
t13 = sin(t20);
t14 = cos(t20);
t35 = t14 * rSges(5,1) - rSges(5,2) * t13;
t34 = -rSges(6,2) * t13 + (rSges(6,1) + pkin(4)) * t14;
t25 = -pkin(7) - pkin(6);
t19 = qJ(1) + pkin(8);
t12 = cos(t19);
t38 = t12 * t14;
t39 = t12 * t13;
t47 = rSges(6,1) * t39 + rSges(6,2) * t38;
t46 = rSges(5,1) * t39 + rSges(5,2) * t38;
t45 = pkin(4) * t13;
t11 = sin(t19);
t44 = g(2) * t11;
t21 = sin(qJ(3));
t43 = t21 * pkin(3);
t42 = rSges(4,3) + pkin(6);
t23 = cos(qJ(3));
t16 = t23 * pkin(3);
t10 = t16 + pkin(2);
t37 = rSges(5,3) - t25;
t36 = rSges(6,3) + qJ(5) - t25;
t22 = sin(qJ(1));
t15 = t22 * pkin(1);
t24 = cos(qJ(1));
t17 = t24 * pkin(1);
t33 = g(2) * t17 + g(3) * t15;
t32 = t23 * rSges(4,1) - t21 * rSges(4,2);
t31 = rSges(4,1) * t21 + rSges(4,2) * t23;
t30 = -rSges(5,1) * t13 - rSges(5,2) * t14;
t29 = -rSges(6,1) * t13 - rSges(6,2) * t14;
t28 = pkin(2) + t32;
t27 = t10 + t34;
t26 = t10 + t35;
t6 = -t43 - t45;
t1 = [-m(2) * (g(2) * (t24 * rSges(2,1) - t22 * rSges(2,2)) + g(3) * (t22 * rSges(2,1) + t24 * rSges(2,2))) - m(3) * (g(2) * (rSges(3,1) * t12 - rSges(3,2) * t11 + t17) + g(3) * (rSges(3,1) * t11 + rSges(3,2) * t12 + t15)) - m(4) * ((g(2) * t28 - g(3) * t42) * t12 + (g(2) * t42 + g(3) * t28) * t11 + t33) - m(5) * ((g(2) * t26 - g(3) * t37) * t12 + (g(2) * t37 + g(3) * t26) * t11 + t33) - m(6) * ((g(2) * t27 - g(3) * t36) * t12 + (g(2) * t36 + g(3) * t27) * t11 + t33), (-m(3) - m(4) - m(5) - m(6)) * g(1), -m(4) * (g(3) * t31 * t12 + g(1) * t32) - m(5) * (g(1) * (t16 + t35) + g(3) * (t12 * t43 + t46)) - m(6) * (g(1) * (t16 + t34) + g(3) * (-t12 * t6 + t47)) + (m(4) * t31 - m(5) * (t30 - t43) - m(6) * (t29 + t6)) * t44, -m(5) * (g(1) * t35 + g(3) * t46) - m(6) * (g(1) * t34 + g(3) * (pkin(4) * t39 + t47)) + (-m(5) * t30 - m(6) * (t29 - t45)) * t44, -m(6) * (-g(2) * t12 - g(3) * t11)];
taug = t1(:);
