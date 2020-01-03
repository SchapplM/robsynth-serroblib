% Calculate Gravitation load on the joints for
% S5RPPRR4
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
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:30:26
% EndTime: 2020-01-03 11:30:29
% DurationCPUTime: 0.50s
% Computational Cost: add. (217->92), mult. (259->129), div. (0->0), fcn. (249->10), ass. (0->39)
t34 = sin(qJ(1));
t35 = cos(qJ(1));
t52 = g(2) * t35 + g(3) * t34;
t28 = pkin(9) + qJ(4);
t22 = qJ(5) + t28;
t17 = sin(t22);
t18 = cos(t22);
t32 = cos(pkin(8));
t43 = t34 * t32;
t5 = -t17 * t43 - t35 * t18;
t6 = -t35 * t17 + t18 * t43;
t51 = t5 * rSges(6,1) - t6 * rSges(6,2);
t42 = t35 * t32;
t7 = t17 * t42 - t34 * t18;
t8 = t34 * t17 + t18 * t42;
t50 = t7 * rSges(6,1) + t8 * rSges(6,2);
t20 = sin(t28);
t49 = pkin(4) * t20;
t30 = sin(pkin(8));
t48 = g(1) * t30;
t29 = sin(pkin(9));
t45 = t29 * pkin(3);
t31 = cos(pkin(9));
t19 = t31 * pkin(3) + pkin(2);
t44 = rSges(3,2) * t30;
t33 = -pkin(6) - qJ(3);
t41 = t35 * pkin(1) + t34 * qJ(2);
t39 = -m(4) - m(5) - m(6);
t38 = -rSges(6,1) * t17 - rSges(6,2) * t18;
t21 = cos(t28);
t11 = t20 * t42 - t34 * t21;
t9 = -t20 * t43 - t35 * t21;
t37 = t19 * t32 + (rSges(5,3) - t33) * t30;
t36 = (pkin(4) * t21 + t19) * t32 + (rSges(6,3) + pkin(7) - t33) * t30;
t25 = t34 * pkin(1);
t15 = t45 + t49;
t12 = t34 * t20 + t21 * t42;
t10 = -t35 * t20 + t21 * t43;
t1 = [-m(2) * (g(2) * (t35 * rSges(2,1) - t34 * rSges(2,2)) + g(3) * (t34 * rSges(2,1) + t35 * rSges(2,2))) - m(3) * (g(2) * (t34 * rSges(3,3) + t41) + g(3) * (rSges(3,1) * t43 - t34 * t44 + t25) + (g(2) * (rSges(3,1) * t32 - t44) + g(3) * (-rSges(3,3) - qJ(2))) * t35) - m(4) * (g(2) * (pkin(2) * t42 + (t34 * t29 + t31 * t42) * rSges(4,1) + (-t29 * t42 + t34 * t31) * rSges(4,2) + t41) + g(3) * (t25 - t35 * qJ(2) + pkin(2) * t43 + (-t35 * t29 + t31 * t43) * rSges(4,1) + (-t29 * t43 - t35 * t31) * rSges(4,2)) + t52 * t30 * (rSges(4,3) + qJ(3))) - m(5) * (g(2) * (t12 * rSges(5,1) - t11 * rSges(5,2) + t41) + g(3) * (t10 * rSges(5,1) + t9 * rSges(5,2) + t25) + (g(2) * t45 + g(3) * t37) * t34 + (g(2) * t37 + g(3) * (-qJ(2) - t45)) * t35) - m(6) * (g(2) * (t8 * rSges(6,1) - t7 * rSges(6,2) + t41) + g(3) * (t6 * rSges(6,1) + t5 * rSges(6,2) + t25) + (g(2) * t15 + g(3) * t36) * t34 + (g(2) * t36 + g(3) * (-qJ(2) - t15)) * t35), -(-m(3) + t39) * t52, t39 * (-g(1) * t32 + (g(2) * t34 - g(3) * t35) * t30), -m(5) * (g(2) * (t9 * rSges(5,1) - t10 * rSges(5,2)) + g(3) * (t11 * rSges(5,1) + t12 * rSges(5,2))) - m(6) * (g(2) * (pkin(4) * t9 + t51) + g(3) * (pkin(4) * t11 + t50)) + (-m(5) * (-rSges(5,1) * t20 - rSges(5,2) * t21) - m(6) * (t38 - t49)) * t48, -m(6) * (g(2) * t51 + g(3) * t50 + t38 * t48)];
taug = t1(:);
