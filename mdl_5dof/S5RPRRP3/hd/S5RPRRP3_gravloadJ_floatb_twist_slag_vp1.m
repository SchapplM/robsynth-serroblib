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
% Datum: 2019-12-05 18:04
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:03:10
% EndTime: 2019-12-05 18:03:11
% DurationCPUTime: 0.31s
% Computational Cost: add. (204->71), mult. (174->92), div. (0->0), fcn. (134->8), ass. (0->35)
t18 = qJ(3) + qJ(4);
t13 = sin(t18);
t14 = cos(t18);
t33 = t14 * rSges(5,1) - t13 * rSges(5,2);
t32 = -t13 * rSges(6,2) + (rSges(6,1) + pkin(4)) * t14;
t23 = -pkin(7) - pkin(6);
t17 = qJ(1) + pkin(8);
t11 = sin(t17);
t38 = t11 * t14;
t39 = t11 * t13;
t45 = rSges(6,1) * t39 + rSges(6,2) * t38;
t44 = rSges(5,1) * t39 + rSges(5,2) * t38;
t43 = pkin(4) * t13;
t12 = cos(t17);
t42 = g(3) * t12;
t19 = sin(qJ(3));
t41 = t19 * pkin(3);
t40 = rSges(4,3) + pkin(6);
t21 = cos(qJ(3));
t15 = t21 * pkin(3);
t10 = t15 + pkin(2);
t35 = rSges(5,3) - t23;
t34 = rSges(6,3) + qJ(5) - t23;
t31 = t21 * rSges(4,1) - t19 * rSges(4,2);
t30 = rSges(4,1) * t19 + rSges(4,2) * t21;
t29 = -rSges(5,1) * t13 - rSges(5,2) * t14;
t28 = -rSges(6,1) * t13 - rSges(6,2) * t14;
t27 = -pkin(2) - t31;
t26 = -t10 - t32;
t25 = -t10 - t33;
t20 = sin(qJ(1));
t22 = cos(qJ(1));
t24 = (-g(2) * t22 - g(3) * t20) * pkin(1);
t6 = -t41 - t43;
t1 = [-m(2) * (g(2) * (-t22 * rSges(2,1) + t20 * rSges(2,2)) + g(3) * (-t20 * rSges(2,1) - t22 * rSges(2,2))) - m(3) * (g(2) * (-t12 * rSges(3,1) + t11 * rSges(3,2) - t22 * pkin(1)) + g(3) * (-t11 * rSges(3,1) - t12 * rSges(3,2) - t20 * pkin(1))) - m(4) * (t24 + (g(2) * t27 + g(3) * t40) * t12 + (-g(2) * t40 + g(3) * t27) * t11) - m(5) * (t24 + (g(2) * t25 + g(3) * t35) * t12 + (-g(2) * t35 + g(3) * t25) * t11) - m(6) * (t24 + (g(2) * t26 + g(3) * t34) * t12 + (-g(2) * t34 + g(3) * t26) * t11), (-m(3) - m(4) - m(5) - m(6)) * g(1), -m(4) * (g(2) * t30 * t11 + g(1) * t31) - m(5) * (g(1) * (t15 + t33) + g(2) * (t11 * t41 + t44)) - m(6) * (g(1) * (t15 + t32) + g(2) * (-t11 * t6 + t45)) + (m(4) * t30 - m(5) * (t29 - t41) - m(6) * (t28 + t6)) * t42, -m(5) * (g(1) * t33 + g(2) * t44) - m(6) * (g(1) * t32 + g(2) * (pkin(4) * t39 + t45)) + (-m(5) * t29 - m(6) * (t28 - t43)) * t42, -m(6) * (g(2) * t12 + g(3) * t11)];
taug = t1(:);
