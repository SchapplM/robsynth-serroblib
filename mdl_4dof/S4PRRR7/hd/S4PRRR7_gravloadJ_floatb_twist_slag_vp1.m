% Calculate Gravitation load on the joints for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRR7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(8,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR7_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:03
% EndTime: 2019-12-31 16:36:04
% DurationCPUTime: 0.30s
% Computational Cost: add. (163->70), mult. (409->118), div. (0->0), fcn. (474->10), ass. (0->37)
t41 = rSges(4,3) + pkin(6);
t40 = rSges(5,3) + pkin(7);
t19 = sin(pkin(4));
t22 = sin(qJ(2));
t39 = t19 * t22;
t24 = cos(qJ(3));
t38 = t19 * t24;
t25 = cos(qJ(2));
t37 = t19 * t25;
t36 = cos(pkin(4));
t35 = cos(pkin(8));
t34 = g(3) * (pkin(2) * t37 + pkin(6) * t39);
t18 = sin(pkin(8));
t33 = t18 * t36;
t32 = t19 * t35;
t21 = sin(qJ(3));
t31 = -rSges(4,1) * t24 + rSges(4,2) * t21;
t20 = sin(qJ(4));
t23 = cos(qJ(4));
t30 = t20 * rSges(5,1) + t23 * rSges(5,2);
t29 = t36 * t35;
t28 = rSges(5,1) * t23 - rSges(5,2) * t20 + pkin(3);
t27 = pkin(6) + t30;
t26 = -t40 * t21 - t28 * t24;
t12 = t36 * t21 + t22 * t38;
t11 = -t21 * t39 + t36 * t24;
t10 = -t22 * t33 + t35 * t25;
t9 = t35 * t22 + t25 * t33;
t8 = t18 * t25 + t22 * t29;
t7 = t18 * t22 - t25 * t29;
t6 = t9 * pkin(2);
t5 = t7 * pkin(2);
t4 = t18 * t19 * t21 + t10 * t24;
t3 = -t10 * t21 + t18 * t38;
t2 = -t21 * t32 + t8 * t24;
t1 = -t8 * t21 - t24 * t32;
t13 = [(-m(2) - m(3) - m(4) - m(5)) * g(3), -m(3) * (g(1) * (-t9 * rSges(3,1) - t10 * rSges(3,2)) + g(2) * (-t7 * rSges(3,1) - t8 * rSges(3,2))) - m(4) * (g(1) * (t41 * t10 + t31 * t9 - t6) + g(2) * (t31 * t7 + t41 * t8 - t5) + t34) - m(5) * (g(1) * (t27 * t10 + t26 * t9 - t6) + g(2) * (t26 * t7 + t27 * t8 - t5) + t34) + ((m(3) * rSges(3,2) - m(4) * rSges(4,3) - m(5) * t30) * t22 + (-m(3) * rSges(3,1) + (m(4) * rSges(4,2) - m(5) * t40) * t21 + (-m(4) * rSges(4,1) - m(5) * t28) * t24) * t25) * g(3) * t19, -m(4) * (g(1) * (t3 * rSges(4,1) - t4 * rSges(4,2)) + g(2) * (t1 * rSges(4,1) - t2 * rSges(4,2)) + g(3) * (t11 * rSges(4,1) - t12 * rSges(4,2))) - m(5) * (g(1) * (t28 * t3 + t40 * t4) + (t28 * t11 + t40 * t12) * g(3) + (t28 * t1 + t40 * t2) * g(2)), -m(5) * (g(1) * ((-t4 * t20 + t9 * t23) * rSges(5,1) + (-t9 * t20 - t4 * t23) * rSges(5,2)) + g(2) * ((-t2 * t20 + t7 * t23) * rSges(5,1) + (-t2 * t23 - t7 * t20) * rSges(5,2)) + g(3) * ((-t12 * t20 - t23 * t37) * rSges(5,1) + (-t12 * t23 + t20 * t37) * rSges(5,2)))];
taug = t13(:);
