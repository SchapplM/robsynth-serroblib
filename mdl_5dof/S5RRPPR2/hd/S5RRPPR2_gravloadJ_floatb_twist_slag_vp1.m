% Calculate Gravitation load on the joints for
% S5RRPPR2
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
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:57:21
% EndTime: 2020-01-03 11:57:21
% DurationCPUTime: 0.19s
% Computational Cost: add. (281->67), mult. (192->90), div. (0->0), fcn. (168->10), ass. (0->39)
t61 = rSges(6,3) + pkin(7);
t60 = -m(5) - m(6);
t37 = qJ(1) + qJ(2);
t32 = pkin(8) + t37;
t28 = sin(t32);
t38 = sin(pkin(9));
t59 = t28 * t38;
t39 = cos(pkin(9));
t58 = t28 * t39;
t29 = cos(t32);
t57 = t29 * t38;
t56 = t29 * t39;
t40 = sin(qJ(5));
t55 = t39 * t40;
t42 = cos(qJ(5));
t54 = t39 * t42;
t33 = sin(t37);
t30 = pkin(2) * t33;
t53 = t28 * pkin(3) + t30;
t34 = cos(t37);
t52 = t33 * rSges(3,1) + t34 * rSges(3,2);
t31 = pkin(2) * t34;
t51 = t29 * pkin(3) + t28 * qJ(4) + t31;
t50 = t28 * rSges(4,1) + t29 * rSges(4,2) + t30;
t49 = t34 * rSges(3,1) - rSges(3,2) * t33;
t48 = t29 * rSges(4,1) - rSges(4,2) * t28 + t31;
t7 = -t28 * t42 + t29 * t55;
t8 = t28 * t40 + t29 * t54;
t47 = t8 * rSges(6,1) - t7 * rSges(6,2) + pkin(4) * t56 + t61 * t57 + t51;
t46 = rSges(5,1) * t56 - rSges(5,2) * t57 + t28 * rSges(5,3) + t51;
t5 = -t28 * t55 - t29 * t42;
t6 = t28 * t54 - t29 * t40;
t45 = t6 * rSges(6,1) + t5 * rSges(6,2) + pkin(4) * t58 - qJ(4) * t29 + t61 * t59 + t53;
t44 = -rSges(5,2) * t59 + rSges(5,1) * t58 + (-rSges(5,3) - qJ(4)) * t29 + t53;
t43 = cos(qJ(1));
t41 = sin(qJ(1));
t36 = t43 * pkin(1);
t35 = t41 * pkin(1);
t1 = [-m(2) * (g(2) * (rSges(2,1) * t43 - t41 * rSges(2,2)) + g(3) * (t41 * rSges(2,1) + rSges(2,2) * t43)) - m(3) * (g(2) * (t36 + t49) + g(3) * (t35 + t52)) - m(4) * (g(2) * (t36 + t48) + g(3) * (t35 + t50)) - m(5) * (g(2) * (t36 + t46) + g(3) * (t35 + t44)) - m(6) * (g(2) * (t36 + t47) + g(3) * (t35 + t45)), -m(3) * (g(2) * t49 + g(3) * t52) - m(4) * (g(2) * t48 + g(3) * t50) - m(5) * (g(2) * t46 + g(3) * t44) - m(6) * (g(2) * t47 + g(3) * t45), (-m(4) + t60) * g(1), t60 * (-g(2) * t29 - g(3) * t28), -m(6) * (g(2) * (rSges(6,1) * t5 - rSges(6,2) * t6) + g(3) * (rSges(6,1) * t7 + rSges(6,2) * t8) + g(1) * (-rSges(6,1) * t40 - rSges(6,2) * t42) * t38)];
taug = t1(:);
