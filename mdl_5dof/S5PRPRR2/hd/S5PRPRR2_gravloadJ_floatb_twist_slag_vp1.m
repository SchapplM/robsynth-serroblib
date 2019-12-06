% Calculate Gravitation load on the joints for
% S5PRPRR2
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
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:44:46
% EndTime: 2019-12-05 15:44:48
% DurationCPUTime: 0.28s
% Computational Cost: add. (194->54), mult. (185->76), div. (0->0), fcn. (153->10), ass. (0->37)
t23 = qJ(2) + pkin(9);
t21 = qJ(4) + t23;
t17 = sin(t21);
t18 = cos(t21);
t26 = sin(qJ(5));
t51 = rSges(6,2) * t26;
t58 = rSges(6,3) + pkin(7);
t59 = t17 * t51 + t18 * t58;
t28 = cos(qJ(5));
t57 = -t28 * rSges(6,1) - pkin(4);
t24 = sin(pkin(8));
t25 = cos(pkin(8));
t55 = g(1) * t25 + g(2) * t24;
t27 = sin(qJ(2));
t54 = pkin(2) * t27;
t48 = t24 * t26;
t47 = t24 * t28;
t46 = t25 * t26;
t45 = t25 * t28;
t20 = cos(t23);
t29 = cos(qJ(2));
t22 = t29 * pkin(2);
t43 = pkin(3) * t20 + t22;
t42 = -m(4) - m(5) - m(6);
t41 = t59 * t24;
t40 = t59 * t25;
t37 = t18 * rSges(5,1) - rSges(5,2) * t17;
t35 = -rSges(5,1) * t17 - rSges(5,2) * t18;
t34 = t35 * t24;
t33 = t35 * t25;
t32 = t58 * t17 + (-t51 - t57) * t18;
t30 = t55 * t57 * t17;
t19 = sin(t23);
t5 = -pkin(3) * t19 - t54;
t2 = t25 * t5;
t1 = t24 * t5;
t3 = [(-m(2) - m(3) + t42) * g(3), -m(3) * (g(3) * (rSges(3,1) * t29 - t27 * rSges(3,2)) + t55 * (-rSges(3,1) * t27 - rSges(3,2) * t29)) - m(4) * (g(3) * (rSges(4,1) * t20 - rSges(4,2) * t19 + t22) + t55 * (-rSges(4,1) * t19 - rSges(4,2) * t20 - t54)) - m(5) * (g(1) * (t2 + t33) + g(2) * (t1 + t34) + g(3) * (t37 + t43)) - m(6) * (g(1) * (t2 + t40) + g(2) * (t1 + t41) + g(3) * (t32 + t43) + t30), t42 * (g(1) * t24 - g(2) * t25), -m(5) * (g(1) * t33 + g(2) * t34 + g(3) * t37) - m(6) * (g(1) * t40 + g(2) * t41 + g(3) * t32 + t30), -m(6) * (g(1) * ((-t18 * t46 + t47) * rSges(6,1) + (-t18 * t45 - t48) * rSges(6,2)) + g(2) * ((-t18 * t48 - t45) * rSges(6,1) + (-t18 * t47 + t46) * rSges(6,2)) + g(3) * (-rSges(6,1) * t26 - rSges(6,2) * t28) * t17)];
taug = t3(:);
