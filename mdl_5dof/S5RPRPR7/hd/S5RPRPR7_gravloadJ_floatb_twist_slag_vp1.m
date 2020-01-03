% Calculate Gravitation load on the joints for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR7_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:18:54
% EndTime: 2019-12-31 18:18:55
% DurationCPUTime: 0.32s
% Computational Cost: add. (226->75), mult. (202->98), div. (0->0), fcn. (173->10), ass. (0->40)
t44 = rSges(6,3) + pkin(7);
t17 = sin(qJ(5));
t20 = cos(qJ(5));
t43 = m(5) * rSges(5,1) + m(6) * (rSges(6,1) * t20 - rSges(6,2) * t17 + pkin(4));
t42 = -m(5) - m(6);
t18 = sin(qJ(3));
t40 = pkin(3) * t18;
t19 = sin(qJ(1));
t39 = t19 * pkin(1);
t14 = qJ(3) + pkin(9);
t8 = sin(t14);
t38 = t8 * rSges(5,2);
t15 = qJ(1) + pkin(8);
t9 = sin(t15);
t37 = t9 * t17;
t36 = t9 * t20;
t35 = rSges(4,3) + pkin(6);
t11 = cos(t15);
t22 = cos(qJ(1));
t13 = t22 * pkin(1);
t21 = cos(qJ(3));
t12 = t21 * pkin(3);
t7 = t12 + pkin(2);
t34 = t11 * t7 + t13;
t33 = t11 * t17;
t32 = t11 * t20;
t16 = -qJ(4) - pkin(6);
t31 = rSges(5,3) - t16;
t30 = g(1) * t39;
t29 = t44 * t8;
t10 = cos(t14);
t28 = t10 * rSges(5,1) - t38;
t27 = t21 * rSges(4,1) - t18 * rSges(4,2);
t26 = pkin(2) + t27;
t24 = t10 * pkin(4) + t29;
t4 = t10 * t32 + t37;
t3 = -t10 * t33 + t36;
t2 = -t10 * t36 + t33;
t1 = t10 * t37 + t32;
t5 = [-m(2) * (g(1) * (-t19 * rSges(2,1) - t22 * rSges(2,2)) + g(2) * (t22 * rSges(2,1) - t19 * rSges(2,2))) - m(3) * (g(1) * (-t9 * rSges(3,1) - t11 * rSges(3,2) - t39) + g(2) * (t11 * rSges(3,1) - t9 * rSges(3,2) + t13)) - m(4) * (-t30 + g(2) * t13 + (-g(1) * t26 + g(2) * t35) * t9 + (g(1) * t35 + g(2) * t26) * t11) - m(5) * (-t30 + g(2) * t34 + (g(1) * t31 + g(2) * t28) * t11 + (g(1) * (-t28 - t7) + g(2) * t31) * t9) - m(6) * (g(1) * (t2 * rSges(6,1) + t1 * rSges(6,2) - t39) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t34) + (-g(1) * t16 + g(2) * t24) * t11 + (g(1) * (-t24 - t7) - g(2) * t16) * t9), (-m(3) - m(4) + t42) * g(3), (-m(4) * t27 - m(5) * (t12 - t38) - m(6) * (t12 + t29) - t43 * t10) * g(3) + (g(1) * t11 + g(2) * t9) * (-m(4) * (-rSges(4,1) * t18 - rSges(4,2) * t21) - m(5) * (-rSges(5,2) * t10 - t40) - m(6) * (t44 * t10 - t40) + t43 * t8), t42 * (g(1) * t9 - g(2) * t11), -m(6) * (g(1) * (t3 * rSges(6,1) - t4 * rSges(6,2)) + g(2) * (-t1 * rSges(6,1) + t2 * rSges(6,2)) + g(3) * (-rSges(6,1) * t17 - rSges(6,2) * t20) * t8)];
taug = t5(:);
