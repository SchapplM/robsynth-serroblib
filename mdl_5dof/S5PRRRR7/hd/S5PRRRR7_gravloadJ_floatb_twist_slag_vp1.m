% Calculate Gravitation load on the joints for
% S5PRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR7_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:11:48
% EndTime: 2019-12-05 17:11:49
% DurationCPUTime: 0.38s
% Computational Cost: add. (261->84), mult. (318->122), div. (0->0), fcn. (297->10), ass. (0->41)
t22 = sin(pkin(9));
t23 = cos(pkin(9));
t58 = g(1) * t23 + g(2) * t22;
t28 = -pkin(7) - pkin(6);
t21 = qJ(3) + qJ(4);
t18 = qJ(5) + t21;
t13 = sin(t18);
t14 = cos(t18);
t27 = cos(qJ(2));
t47 = t22 * t27;
t57 = (-t13 * t47 - t23 * t14) * rSges(6,1) + (t23 * t13 - t14 * t47) * rSges(6,2);
t46 = t23 * t27;
t56 = (-t13 * t46 + t22 * t14) * rSges(6,1) + (-t22 * t13 - t14 * t46) * rSges(6,2);
t16 = sin(t21);
t17 = cos(t21);
t35 = -t16 * t47 - t23 * t17;
t55 = t35 * rSges(5,1) + (t23 * t16 - t17 * t47) * rSges(5,2);
t36 = -t16 * t46 + t22 * t17;
t54 = t36 * rSges(5,1) + (-t22 * t16 - t17 * t46) * rSges(5,2);
t53 = pkin(4) * t16;
t25 = sin(qJ(2));
t50 = g(3) * t25;
t24 = sin(qJ(3));
t49 = t24 * pkin(3);
t48 = rSges(4,3) + pkin(6);
t45 = t24 * t27;
t26 = cos(qJ(3));
t44 = t26 * t27;
t43 = rSges(5,3) - t28;
t42 = rSges(6,3) + pkin(8) - t28;
t19 = t26 * pkin(3);
t11 = pkin(4) * t17 + t19;
t40 = -rSges(5,1) * t16 - rSges(5,2) * t17;
t39 = -rSges(6,1) * t13 - rSges(6,2) * t14;
t38 = t26 * rSges(4,1) - t24 * rSges(4,2) + pkin(2);
t37 = rSges(6,1) * t14 - rSges(6,2) * t13 + pkin(2) + t11;
t34 = t22 * t26 - t23 * t45;
t33 = -t22 * t45 - t23 * t26;
t32 = rSges(5,1) * t17 - rSges(5,2) * t16 + pkin(2) + t19;
t10 = -t49 - t53;
t1 = [(-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3), -m(3) * (g(3) * (t27 * rSges(3,1) - t25 * rSges(3,2)) + t58 * (-rSges(3,1) * t25 - rSges(3,2) * t27)) - m(4) * (g(3) * (t48 * t25 + t38 * t27) + t58 * (-t38 * t25 + t48 * t27)) - m(5) * (g(3) * (t43 * t25 + t32 * t27) + t58 * (-t32 * t25 + t43 * t27)) - m(6) * (g(3) * (t42 * t25 + t37 * t27) + t58 * (-t37 * t25 + t42 * t27)), -m(4) * (g(1) * (t34 * rSges(4,1) + (-t22 * t24 - t23 * t44) * rSges(4,2)) + g(2) * (t33 * rSges(4,1) + (-t22 * t44 + t23 * t24) * rSges(4,2))) - m(5) * (g(1) * (t34 * pkin(3) + t54) + g(2) * (t33 * pkin(3) + t55)) - m(6) * (g(1) * (t10 * t46 + t22 * t11 + t56) + g(2) * (t10 * t47 - t23 * t11 + t57)) + (-m(4) * (-rSges(4,1) * t24 - rSges(4,2) * t26) - m(5) * (t40 - t49) - m(6) * (t10 + t39)) * t50, -m(5) * (g(1) * t54 + g(2) * t55) - m(6) * (g(1) * (t36 * pkin(4) + t56) + g(2) * (t35 * pkin(4) + t57)) + (-m(5) * t40 - m(6) * (t39 - t53)) * t50, -m(6) * (g(1) * t56 + g(2) * t57 + t39 * t50)];
taug = t1(:);
