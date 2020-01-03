% Calculate Gravitation load on the joints for
% S5RPRRR6
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
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR6_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:00:57
% EndTime: 2019-12-31 19:00:58
% DurationCPUTime: 0.43s
% Computational Cost: add. (274->79), mult. (245->106), div. (0->0), fcn. (210->10), ass. (0->39)
t24 = qJ(1) + pkin(9);
t18 = sin(t24);
t19 = cos(t24);
t59 = g(1) * t19 + g(2) * t18;
t63 = rSges(6,3) + pkin(8);
t25 = qJ(3) + qJ(4);
t20 = sin(t25);
t21 = cos(t25);
t60 = t21 * rSges(5,1) - t20 * rSges(5,2);
t36 = t21 * pkin(4) + t63 * t20;
t29 = cos(qJ(5));
t58 = (-rSges(6,1) * t29 - pkin(4)) * t20;
t27 = sin(qJ(3));
t57 = pkin(3) * t27;
t28 = sin(qJ(1));
t54 = t28 * pkin(1);
t53 = rSges(4,3) + pkin(6);
t30 = cos(qJ(3));
t22 = t30 * pkin(3);
t17 = t22 + pkin(2);
t31 = cos(qJ(1));
t23 = t31 * pkin(1);
t52 = t19 * t17 + t23;
t26 = sin(qJ(5));
t47 = t21 * t26;
t46 = t21 * t29;
t32 = -pkin(7) - pkin(6);
t45 = rSges(5,3) - t32;
t44 = g(1) * t54;
t40 = t30 * rSges(4,1) - t27 * rSges(4,2);
t38 = -rSges(5,1) * t20 - rSges(5,2) * t21;
t37 = pkin(2) + t40;
t35 = rSges(6,1) * t46 - rSges(6,2) * t47 + t36;
t34 = t59 * (rSges(6,2) * t20 * t26 + t21 * t63);
t4 = t18 * t26 + t19 * t46;
t3 = t18 * t29 - t19 * t47;
t2 = -t18 * t46 + t19 * t26;
t1 = t18 * t47 + t19 * t29;
t5 = [-m(2) * (g(1) * (-t28 * rSges(2,1) - t31 * rSges(2,2)) + g(2) * (t31 * rSges(2,1) - t28 * rSges(2,2))) - m(3) * (g(1) * (-t18 * rSges(3,1) - t19 * rSges(3,2) - t54) + g(2) * (t19 * rSges(3,1) - t18 * rSges(3,2) + t23)) - m(4) * (-t44 + g(2) * t23 + (g(1) * t53 + g(2) * t37) * t19 + (-g(1) * t37 + g(2) * t53) * t18) - m(5) * (-t44 + g(2) * t52 + (g(1) * t45 + g(2) * t60) * t19 + (g(1) * (-t17 - t60) + g(2) * t45) * t18) - m(6) * (g(1) * (t2 * rSges(6,1) + t1 * rSges(6,2) - t54) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t52) + (-g(1) * t32 + g(2) * t36) * t19 + (g(1) * (-t17 - t36) - g(2) * t32) * t18), (-m(3) - m(4) - m(5) - m(6)) * g(3), -m(6) * t34 + (-m(4) * t40 - m(5) * (t22 + t60) - m(6) * (t22 + t35)) * g(3) + t59 * (-m(4) * (-rSges(4,1) * t27 - rSges(4,2) * t30) - m(5) * (t38 - t57) - m(6) * (-t57 + t58)), -m(5) * (g(3) * t60 + t59 * t38) - m(6) * (g(3) * t35 + t59 * t58 + t34), -m(6) * (g(1) * (t3 * rSges(6,1) - t4 * rSges(6,2)) + g(2) * (-t1 * rSges(6,1) + t2 * rSges(6,2)) + g(3) * (-rSges(6,1) * t26 - rSges(6,2) * t29) * t20)];
taug = t5(:);
