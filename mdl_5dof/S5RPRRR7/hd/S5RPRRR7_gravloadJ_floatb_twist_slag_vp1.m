% Calculate Gravitation load on the joints for
% S5RPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR7_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:03:12
% EndTime: 2019-12-31 19:03:13
% DurationCPUTime: 0.40s
% Computational Cost: add. (283->87), mult. (290->124), div. (0->0), fcn. (267->10), ass. (0->42)
t56 = rSges(6,3) + pkin(8) + pkin(7);
t25 = sin(qJ(3));
t28 = cos(qJ(3));
t55 = t28 * rSges(4,1) - t25 * rSges(4,2);
t45 = rSges(5,3) + pkin(7);
t54 = t28 * pkin(3) + t45 * t25;
t27 = cos(qJ(4));
t16 = t27 * pkin(4) + pkin(3);
t23 = qJ(4) + qJ(5);
t19 = sin(t23);
t20 = cos(t23);
t24 = sin(qJ(4));
t53 = m(5) * (rSges(5,1) * t27 - rSges(5,2) * t24 + pkin(3)) + m(6) * (rSges(6,1) * t20 - rSges(6,2) * t19 + t16) + m(4) * rSges(4,1);
t22 = qJ(1) + pkin(9);
t17 = sin(t22);
t18 = cos(t22);
t44 = t19 * t28;
t5 = t17 * t44 + t18 * t20;
t43 = t20 * t28;
t6 = -t17 * t43 + t18 * t19;
t52 = -t5 * rSges(6,1) + t6 * rSges(6,2);
t7 = t17 * t20 - t18 * t44;
t8 = t17 * t19 + t18 * t43;
t51 = t7 * rSges(6,1) - t8 * rSges(6,2);
t49 = pkin(4) * t24;
t48 = g(3) * t25;
t26 = sin(qJ(1));
t47 = t26 * pkin(1);
t42 = t24 * t28;
t40 = t27 * t28;
t29 = cos(qJ(1));
t21 = t29 * pkin(1);
t38 = t18 * pkin(2) + t17 * pkin(6) + t21;
t37 = t18 * pkin(6) - t47;
t36 = -rSges(6,1) * t19 - rSges(6,2) * t20;
t11 = t17 * t27 - t18 * t42;
t9 = t17 * t42 + t18 * t27;
t33 = t28 * t16 + t25 * t56;
t32 = m(4) * rSges(4,2) - m(5) * t45 - m(6) * t56;
t12 = t17 * t24 + t18 * t40;
t10 = -t17 * t40 + t18 * t24;
t1 = [-m(2) * (g(1) * (-t26 * rSges(2,1) - t29 * rSges(2,2)) + g(2) * (t29 * rSges(2,1) - t26 * rSges(2,2))) - m(3) * (g(1) * (-t17 * rSges(3,1) - t18 * rSges(3,2) - t47) + g(2) * (t18 * rSges(3,1) - t17 * rSges(3,2) + t21)) - m(4) * (g(1) * (t18 * rSges(4,3) + t37) + g(2) * (t18 * t55 + t38) + (g(1) * (-pkin(2) - t55) + g(2) * rSges(4,3)) * t17) - m(5) * ((t12 * rSges(5,1) + t11 * rSges(5,2) + t18 * t54 + t38) * g(2) + (t10 * rSges(5,1) + t9 * rSges(5,2) + t37 + (-pkin(2) - t54) * t17) * g(1)) - m(6) * (g(1) * (t6 * rSges(6,1) + t5 * rSges(6,2) + t37) + g(2) * (t8 * rSges(6,1) + t7 * rSges(6,2) + t38) + (g(1) * t49 + g(2) * t33) * t18 + (g(1) * (-pkin(2) - t33) + g(2) * t49) * t17), (-m(3) - m(4) - m(5) - m(6)) * g(3), (t32 * t25 - t28 * t53) * g(3) + (g(1) * t18 + g(2) * t17) * (t25 * t53 + t32 * t28), -m(5) * (g(1) * (t11 * rSges(5,1) - t12 * rSges(5,2)) + g(2) * (-t9 * rSges(5,1) + t10 * rSges(5,2))) - m(6) * (g(1) * (t11 * pkin(4) + t51) + g(2) * (-t9 * pkin(4) + t52)) + (-m(5) * (-rSges(5,1) * t24 - rSges(5,2) * t27) - m(6) * (t36 - t49)) * t48, -m(6) * (g(1) * t51 + g(2) * t52 + t36 * t48)];
taug = t1(:);
