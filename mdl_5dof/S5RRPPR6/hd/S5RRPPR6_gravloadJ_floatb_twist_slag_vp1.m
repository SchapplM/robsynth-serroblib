% Calculate Gravitation load on the joints for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR6_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:31:36
% EndTime: 2019-12-31 19:31:37
% DurationCPUTime: 0.46s
% Computational Cost: add. (249->83), mult. (284->111), div. (0->0), fcn. (256->10), ass. (0->42)
t16 = qJ(2) + pkin(8);
t11 = sin(t16);
t13 = cos(t16);
t17 = sin(pkin(9));
t18 = cos(pkin(9));
t32 = rSges(5,1) * t18 - rSges(5,2) * t17 + pkin(3);
t39 = rSges(5,3) + qJ(4);
t51 = t39 * t11 + t32 * t13;
t40 = rSges(6,3) + pkin(7) + qJ(4);
t52 = t40 * t11;
t22 = sin(qJ(1));
t24 = cos(qJ(1));
t50 = g(1) * t24 + g(2) * t22;
t23 = cos(qJ(2));
t14 = t23 * pkin(2);
t9 = t14 + pkin(1);
t7 = t24 * t9;
t49 = g(2) * t7;
t48 = -m(5) - m(6);
t21 = sin(qJ(2));
t47 = pkin(2) * t21;
t44 = rSges(3,3) + pkin(6);
t43 = t22 * t13;
t42 = t24 * t13;
t19 = -qJ(3) - pkin(6);
t41 = rSges(4,3) - t19;
t38 = pkin(4) * t17 - t19;
t36 = t23 * rSges(3,1) - t21 * rSges(3,2);
t34 = t13 * rSges(4,1) - t11 * rSges(4,2);
t33 = pkin(1) + t36;
t15 = pkin(9) + qJ(5);
t10 = sin(t15);
t12 = cos(t15);
t8 = t18 * pkin(4) + pkin(3);
t31 = rSges(6,1) * t12 - rSges(6,2) * t10 + t8;
t30 = t17 * rSges(5,1) + t18 * rSges(5,2) - t19;
t28 = t13 * t8 + t52;
t5 = t22 * t10 + t12 * t42;
t4 = -t10 * t42 + t22 * t12;
t3 = t24 * t10 - t12 * t43;
t2 = t10 * t43 + t24 * t12;
t1 = [-m(2) * (g(1) * (-t22 * rSges(2,1) - t24 * rSges(2,2)) + g(2) * (t24 * rSges(2,1) - t22 * rSges(2,2))) - m(3) * ((g(1) * t44 + g(2) * t33) * t24 + (-g(1) * t33 + g(2) * t44) * t22) - m(4) * (t49 + (g(1) * t41 + g(2) * t34) * t24 + (g(1) * (-t34 - t9) + g(2) * t41) * t22) - m(5) * (t49 + (g(1) * t30 + t51 * g(2)) * t24 + (g(2) * t30 + (-t9 - t51) * g(1)) * t22) - m(6) * (g(1) * (t3 * rSges(6,1) + t2 * rSges(6,2)) + g(2) * (t5 * rSges(6,1) + t4 * rSges(6,2) + t7) + (g(1) * t38 + g(2) * t28) * t24 + (g(1) * (-t28 - t9) + g(2) * t38) * t22), -m(3) * (g(3) * t36 + t50 * (-rSges(3,1) * t21 - rSges(3,2) * t23)) - m(4) * (g(3) * (t14 + t34) + t50 * (-rSges(4,1) * t11 - rSges(4,2) * t13 - t47)) - m(5) * (g(3) * (t14 + t51) + t50 * (-t32 * t11 + t39 * t13 - t47)) - m(6) * (g(3) * (t31 * t13 + t14 + t52) + t50 * (-t31 * t11 + t40 * t13 - t47)), (-m(4) + t48) * (g(1) * t22 - g(2) * t24), t48 * (-g(3) * t13 + t50 * t11), -m(6) * (g(1) * (t4 * rSges(6,1) - t5 * rSges(6,2)) + g(2) * (-t2 * rSges(6,1) + t3 * rSges(6,2)) + g(3) * (-rSges(6,1) * t10 - rSges(6,2) * t12) * t11)];
taug = t1(:);
