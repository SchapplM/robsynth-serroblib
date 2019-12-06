% Calculate Gravitation load on the joints for
% S5PRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR5_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:53:26
% EndTime: 2019-12-05 15:53:28
% DurationCPUTime: 0.31s
% Computational Cost: add. (190->59), mult. (235->89), div. (0->0), fcn. (211->10), ass. (0->30)
t17 = sin(pkin(8));
t19 = cos(pkin(8));
t44 = g(1) * t19 + g(2) * t17;
t22 = cos(qJ(2));
t38 = t17 * t22;
t15 = pkin(9) + qJ(4);
t12 = qJ(5) + t15;
t7 = sin(t12);
t8 = cos(t12);
t43 = (-t19 * t8 - t38 * t7) * rSges(6,1) + (t19 * t7 - t38 * t8) * rSges(6,2);
t37 = t19 * t22;
t42 = (t17 * t8 - t37 * t7) * rSges(6,1) + (-t17 * t7 - t37 * t8) * rSges(6,2);
t21 = sin(qJ(2));
t39 = g(3) * t21;
t18 = cos(pkin(9));
t9 = t18 * pkin(3) + pkin(2);
t20 = -pkin(6) - qJ(3);
t36 = rSges(5,3) - t20;
t35 = rSges(6,3) + pkin(7) - t20;
t34 = rSges(4,3) + qJ(3);
t33 = -m(4) - m(5) - m(6);
t32 = -rSges(6,1) * t7 - rSges(6,2) * t8;
t11 = cos(t15);
t30 = rSges(6,1) * t8 - rSges(6,2) * t7 + pkin(4) * t11 + t9;
t29 = rSges(4,1) * t18 - rSges(4,2) * sin(pkin(9)) + pkin(2);
t10 = sin(t15);
t28 = t11 * rSges(5,1) - t10 * rSges(5,2) + t9;
t27 = -t10 * t37 + t17 * t11;
t26 = -t10 * t38 - t19 * t11;
t1 = [(-m(2) - m(3) + t33) * g(3), -m(3) * (g(3) * (t22 * rSges(3,1) - t21 * rSges(3,2)) + t44 * (-rSges(3,1) * t21 - rSges(3,2) * t22)) - m(4) * (g(3) * (t21 * t34 + t22 * t29) + t44 * (-t21 * t29 + t22 * t34)) - m(5) * (g(3) * (t21 * t36 + t22 * t28) + t44 * (-t21 * t28 + t22 * t36)) - m(6) * (g(3) * (t21 * t35 + t22 * t30) + t44 * (-t21 * t30 + t22 * t35)), t33 * (-g(3) * t22 + t44 * t21), -m(5) * (g(1) * (t27 * rSges(5,1) + (-t17 * t10 - t11 * t37) * rSges(5,2)) + g(2) * (t26 * rSges(5,1) + (t19 * t10 - t11 * t38) * rSges(5,2))) - m(6) * (g(1) * (pkin(4) * t27 + t42) + g(2) * (pkin(4) * t26 + t43)) + (-m(5) * (-rSges(5,1) * t10 - rSges(5,2) * t11) - m(6) * (-pkin(4) * t10 + t32)) * t39, -m(6) * (g(1) * t42 + g(2) * t43 + t32 * t39)];
taug = t1(:);
