% Calculate Gravitation load on the joints for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(8,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR6_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:29:22
% EndTime: 2019-12-31 17:29:24
% DurationCPUTime: 0.52s
% Computational Cost: add. (255->102), mult. (617->161), div. (0->0), fcn. (724->10), ass. (0->46)
t32 = sin(qJ(2));
t33 = sin(qJ(1));
t36 = cos(qJ(2));
t50 = cos(pkin(4));
t58 = cos(qJ(1));
t42 = t50 * t58;
t15 = t33 * t32 - t36 * t42;
t30 = sin(qJ(4));
t34 = cos(qJ(4));
t16 = t32 * t42 + t33 * t36;
t31 = sin(qJ(3));
t35 = cos(qJ(3));
t29 = sin(pkin(4));
t47 = t29 * t58;
t4 = t16 * t35 - t31 * t47;
t62 = -t15 * t34 + t4 * t30;
t61 = -t15 * t30 - t4 * t34;
t60 = rSges(4,3) + pkin(7);
t59 = rSges(5,3) + pkin(8);
t55 = t29 * t32;
t54 = t29 * t33;
t53 = t29 * t35;
t52 = t29 * t36;
t51 = t58 * pkin(1) + pkin(6) * t54;
t44 = t33 * t50;
t18 = -t32 * t44 + t58 * t36;
t49 = t18 * pkin(2) + t51;
t48 = g(3) * (pkin(2) * t52 + pkin(7) * t55);
t46 = -t33 * pkin(1) + pkin(6) * t47;
t45 = -t16 * t31 - t35 * t47;
t43 = -t16 * pkin(2) + t46;
t41 = -rSges(4,1) * t35 + rSges(4,2) * t31;
t40 = t30 * rSges(5,1) + t34 * rSges(5,2);
t39 = rSges(5,1) * t34 - rSges(5,2) * t30 + pkin(3);
t38 = pkin(7) + t40;
t37 = -t59 * t31 - t39 * t35;
t17 = t58 * t32 + t36 * t44;
t14 = t50 * t31 + t32 * t53;
t13 = -t31 * t55 + t50 * t35;
t11 = t17 * pkin(2);
t9 = t15 * pkin(2);
t8 = t18 * t35 + t31 * t54;
t7 = t18 * t31 - t33 * t53;
t2 = t17 * t30 + t8 * t34;
t1 = t17 * t34 - t8 * t30;
t3 = [-m(2) * (g(1) * (-t33 * rSges(2,1) - t58 * rSges(2,2)) + g(2) * (t58 * rSges(2,1) - t33 * rSges(2,2))) - m(3) * (g(1) * (-t16 * rSges(3,1) + t15 * rSges(3,2) + rSges(3,3) * t47 + t46) + g(2) * (t18 * rSges(3,1) - t17 * rSges(3,2) + rSges(3,3) * t54 + t51)) - m(4) * (g(1) * (-rSges(4,1) * t4 - rSges(4,2) * t45 - t60 * t15 + t43) + g(2) * (t8 * rSges(4,1) - t7 * rSges(4,2) + t60 * t17 + t49)) - m(5) * (g(1) * (t61 * rSges(5,1) + t62 * rSges(5,2) - t4 * pkin(3) - t15 * pkin(7) + t59 * t45 + t43) + g(2) * (t2 * rSges(5,1) + t1 * rSges(5,2) + t8 * pkin(3) + t17 * pkin(7) + t59 * t7 + t49)), -m(3) * (g(1) * (-t17 * rSges(3,1) - t18 * rSges(3,2)) + g(2) * (-t15 * rSges(3,1) - t16 * rSges(3,2))) - m(4) * (g(1) * (t41 * t17 + t60 * t18 - t11) + g(2) * (t41 * t15 + t60 * t16 - t9) + t48) - m(5) * (g(1) * (t17 * t37 + t18 * t38 - t11) + g(2) * (t15 * t37 + t16 * t38 - t9) + t48) + ((m(3) * rSges(3,2) - m(4) * rSges(4,3) - m(5) * t40) * t32 + (-m(3) * rSges(3,1) + (m(4) * rSges(4,2) - m(5) * t59) * t31 + (-m(4) * rSges(4,1) - m(5) * t39) * t35) * t36) * g(3) * t29, -m(4) * (g(1) * (-t7 * rSges(4,1) - t8 * rSges(4,2)) + g(2) * (rSges(4,1) * t45 - t4 * rSges(4,2)) + g(3) * (t13 * rSges(4,1) - t14 * rSges(4,2))) - m(5) * (g(1) * (-t39 * t7 + t59 * t8) + (t39 * t13 + t59 * t14) * g(3) + (t39 * t45 + t59 * t4) * g(2)), -m(5) * (g(1) * (t1 * rSges(5,1) - t2 * rSges(5,2)) + g(2) * (-t62 * rSges(5,1) + t61 * rSges(5,2)) + g(3) * ((-t14 * t30 - t34 * t52) * rSges(5,1) + (-t14 * t34 + t30 * t52) * rSges(5,2)))];
taug = t3(:);
