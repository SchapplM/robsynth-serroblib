% Calculate Gravitation load on the joints for
% S5RPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR10_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR10_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR10_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR10_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:03:59
% EndTime: 2019-12-31 18:04:00
% DurationCPUTime: 0.35s
% Computational Cost: add. (165->77), mult. (279->116), div. (0->0), fcn. (271->8), ass. (0->36)
t49 = -rSges(5,3) - pkin(6);
t29 = sin(qJ(1));
t48 = g(2) * t29;
t26 = sin(pkin(8));
t28 = sin(qJ(4));
t47 = t26 * t28;
t31 = cos(qJ(1));
t46 = t26 * t31;
t27 = cos(pkin(8));
t45 = t27 * t31;
t44 = -rSges(6,3) - pkin(7) - pkin(6);
t43 = t31 * pkin(1) + t29 * qJ(2);
t42 = qJ(3) * t26;
t41 = -m(4) - m(5) - m(6);
t40 = pkin(2) * t45 + t31 * t42 + t43;
t25 = qJ(4) + qJ(5);
t20 = sin(t25);
t21 = cos(t25);
t39 = t27 * t20 - t26 * t21;
t38 = t26 * t20 + t27 * t21;
t30 = cos(qJ(4));
t37 = t26 * t30 - t27 * t28;
t36 = t27 * t30 + t47;
t35 = pkin(4) * t47 + (t30 * pkin(4) + pkin(3)) * t27;
t34 = -pkin(2) * t27 - pkin(1) - t42;
t14 = t37 * t31;
t5 = t39 * t29;
t6 = t38 * t29;
t7 = t39 * t31;
t8 = t38 * t31;
t33 = m(6) * (g(1) * (-t7 * rSges(6,1) - t8 * rSges(6,2)) + g(2) * (-t5 * rSges(6,1) - t6 * rSges(6,2)) + g(3) * (-t38 * rSges(6,1) + t39 * rSges(6,2)));
t23 = t31 * qJ(2);
t15 = t36 * t31;
t13 = t36 * t29;
t12 = t37 * t29;
t1 = [-m(2) * (g(1) * (-t29 * rSges(2,1) - t31 * rSges(2,2)) + g(2) * (t31 * rSges(2,1) - t29 * rSges(2,2))) - m(3) * (g(1) * (t31 * rSges(3,3) + t23) + g(2) * (rSges(3,1) * t45 - rSges(3,2) * t46 + t43) + (g(1) * (-rSges(3,1) * t27 + rSges(3,2) * t26 - pkin(1)) + g(2) * rSges(3,3)) * t29) - m(4) * (g(1) * (t31 * rSges(4,2) + t23) + g(2) * (rSges(4,1) * t45 + rSges(4,3) * t46 + t40) + (g(1) * (-rSges(4,1) * t27 - rSges(4,3) * t26 + t34) + g(2) * rSges(4,2)) * t29) - m(5) * (g(1) * (-t13 * rSges(5,1) - t12 * rSges(5,2) + t49 * t31 + t23) + g(2) * (t15 * rSges(5,1) + t14 * rSges(5,2) + pkin(3) * t45 + t40) + (g(1) * (-pkin(3) * t27 + t34) + g(2) * t49) * t29) - m(6) * (g(1) * (-t6 * rSges(6,1) + t5 * rSges(6,2) + t23) + g(2) * (t8 * rSges(6,1) - t7 * rSges(6,2) + t40) + (g(1) * t44 + g(2) * t35) * t31 + (g(1) * (t34 - t35) + g(2) * t44) * t29), (-m(3) + t41) * (g(1) * t29 - g(2) * t31), t41 * (-g(3) * t27 + (g(1) * t31 + t48) * t26), -m(5) * (g(1) * (t14 * rSges(5,1) - t15 * rSges(5,2)) + g(2) * (t12 * rSges(5,1) - t13 * rSges(5,2)) + g(3) * (-t36 * rSges(5,1) - t37 * rSges(5,2))) - t33 - m(6) * (g(1) * t14 - g(3) * t36 + t37 * t48) * pkin(4), -t33];
taug = t1(:);
