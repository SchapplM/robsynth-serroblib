% Calculate Gravitation load on the joints for
% S5RPRPR4
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
% m [6x1]
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
% Datum: 2022-01-23 09:23
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:22:25
% EndTime: 2022-01-23 09:22:26
% DurationCPUTime: 0.28s
% Computational Cost: add. (200->59), mult. (158->72), div. (0->0), fcn. (121->10), ass. (0->33)
t17 = qJ(3) + pkin(9);
t11 = cos(t17);
t22 = cos(qJ(3));
t14 = t22 * pkin(3);
t13 = qJ(5) + t17;
t6 = sin(t13);
t7 = cos(t13);
t32 = t7 * rSges(6,1) - t6 * rSges(6,2);
t45 = pkin(4) * t11 + t14 + t32;
t9 = sin(t17);
t44 = t11 * rSges(5,1) - t9 * rSges(5,2) + t14;
t18 = qJ(1) + pkin(8);
t10 = sin(t18);
t12 = cos(t18);
t43 = g(1) * t12 + g(2) * t10;
t42 = -m(5) - m(6);
t20 = sin(qJ(3));
t39 = t20 * pkin(3);
t21 = sin(qJ(1));
t38 = t21 * pkin(1);
t36 = rSges(4,3) + pkin(6);
t19 = -qJ(4) - pkin(6);
t34 = rSges(5,3) - t19;
t33 = rSges(6,3) + pkin(7) - t19;
t31 = -rSges(6,1) * t6 - rSges(6,2) * t7;
t29 = t22 * rSges(4,1) - t20 * rSges(4,2);
t23 = cos(qJ(1));
t15 = t23 * pkin(1);
t28 = -g(1) * t38 + g(2) * t15;
t27 = pkin(2) + t45;
t26 = pkin(2) + t44;
t25 = pkin(2) + t29;
t1 = [-m(2) * (g(1) * (-t21 * rSges(2,1) - t23 * rSges(2,2)) + g(2) * (t23 * rSges(2,1) - t21 * rSges(2,2))) - m(3) * (g(1) * (-t10 * rSges(3,1) - t12 * rSges(3,2) - t38) + g(2) * (t12 * rSges(3,1) - t10 * rSges(3,2) + t15)) - m(4) * ((g(1) * t36 + g(2) * t25) * t12 + (-g(1) * t25 + g(2) * t36) * t10 + t28) - m(5) * ((g(1) * t34 + g(2) * t26) * t12 + (-g(1) * t26 + g(2) * t34) * t10 + t28) - m(6) * ((g(1) * t33 + g(2) * t27) * t12 + (-g(1) * t27 + g(2) * t33) * t10 + t28), (-m(3) - m(4) + t42) * g(3), (-m(4) * t29 - m(5) * t44 - m(6) * t45) * g(3) + t43 * (-m(4) * (-rSges(4,1) * t20 - rSges(4,2) * t22) - m(5) * (-rSges(5,1) * t9 - rSges(5,2) * t11 - t39) - m(6) * (-pkin(4) * t9 + t31 - t39)), t42 * (g(1) * t10 - g(2) * t12), -m(6) * (g(3) * t32 + t31 * t43)];
taug = t1(:);
