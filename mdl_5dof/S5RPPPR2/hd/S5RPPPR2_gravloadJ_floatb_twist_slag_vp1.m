% Calculate Gravitation load on the joints for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPPR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:22:20
% EndTime: 2020-01-03 11:22:22
% DurationCPUTime: 0.33s
% Computational Cost: add. (167->79), mult. (365->117), div. (0->0), fcn. (411->10), ass. (0->43)
t32 = cos(pkin(7));
t34 = sin(qJ(1));
t44 = cos(pkin(8));
t42 = t34 * t44;
t29 = sin(pkin(8));
t36 = cos(qJ(1));
t48 = t36 * t29;
t15 = t32 * t48 - t42;
t33 = sin(qJ(5));
t35 = cos(qJ(5));
t41 = t36 * t44;
t49 = t34 * t29;
t16 = t32 * t41 + t49;
t28 = sin(pkin(9));
t31 = cos(pkin(9));
t30 = sin(pkin(7));
t51 = t30 * t36;
t7 = t16 * t31 + t28 * t51;
t59 = -t15 * t35 + t7 * t33;
t58 = t15 * t33 + t7 * t35;
t57 = -m(5) - m(6);
t56 = rSges(6,3) + pkin(6);
t53 = t29 * t30;
t52 = t30 * t34;
t50 = t32 * t34;
t47 = t36 * pkin(1) + t34 * qJ(2);
t46 = qJ(3) * t30;
t45 = rSges(5,3) + qJ(4);
t43 = -m(4) + t57;
t40 = t47 + (pkin(2) * t32 + t46) * t36;
t39 = t16 * pkin(3) + t40;
t26 = t34 * pkin(1);
t38 = pkin(2) * t50 - t36 * qJ(2) + t34 * t46 + t26;
t14 = t32 * t42 - t48;
t37 = t14 * pkin(3) + t38;
t13 = t32 * t49 + t41;
t12 = t30 * t44 * t31 - t32 * t28;
t6 = t16 * t28 - t31 * t51;
t5 = t14 * t31 + t28 * t52;
t4 = t14 * t28 - t31 * t52;
t2 = t13 * t33 + t5 * t35;
t1 = t13 * t35 - t5 * t33;
t3 = [-m(2) * (g(2) * (t36 * rSges(2,1) - t34 * rSges(2,2)) + g(3) * (t34 * rSges(2,1) + t36 * rSges(2,2))) - m(3) * (g(2) * (t34 * rSges(3,3) + t47) + g(3) * (rSges(3,1) * t50 - rSges(3,2) * t52 + t26) + (g(2) * (rSges(3,1) * t32 - rSges(3,2) * t30) + g(3) * (-rSges(3,3) - qJ(2))) * t36) - m(4) * (g(2) * (t16 * rSges(4,1) - t15 * rSges(4,2) + rSges(4,3) * t51 + t40) + g(3) * (t14 * rSges(4,1) - t13 * rSges(4,2) + rSges(4,3) * t52 + t38)) - m(5) * (g(2) * (t7 * rSges(5,1) - t6 * rSges(5,2) + t45 * t15 + t39) + g(3) * (t5 * rSges(5,1) - t4 * rSges(5,2) + t45 * t13 + t37)) - m(6) * (g(2) * (t58 * rSges(6,1) - t59 * rSges(6,2) + t7 * pkin(4) + t15 * qJ(4) + t56 * t6 + t39) + g(3) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t5 * pkin(4) + t13 * qJ(4) + t56 * t4 + t37)), (-m(3) + t43) * (-g(2) * t36 - g(3) * t34), t43 * (-g(1) * t32 + (g(2) * t34 - g(3) * t36) * t30), t57 * (g(1) * t53 + g(2) * t13 - g(3) * t15), -m(6) * (g(1) * ((-t12 * t33 + t35 * t53) * rSges(6,1) + (-t12 * t35 - t33 * t53) * rSges(6,2)) + g(2) * (t1 * rSges(6,1) - t2 * rSges(6,2)) + g(3) * (t59 * rSges(6,1) + t58 * rSges(6,2)))];
taug = t3(:);
