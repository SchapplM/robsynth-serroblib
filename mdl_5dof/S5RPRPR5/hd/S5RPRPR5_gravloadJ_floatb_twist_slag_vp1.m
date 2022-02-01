% Calculate Gravitation load on the joints for
% S5RPRPR5
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
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:24:56
% EndTime: 2022-01-23 09:24:58
% DurationCPUTime: 0.37s
% Computational Cost: add. (234->102), mult. (286->140), div. (0->0), fcn. (277->10), ass. (0->57)
t69 = -m(5) - m(6);
t39 = cos(pkin(8));
t37 = qJ(3) + pkin(9);
t30 = qJ(5) + t37;
t26 = cos(t30);
t43 = cos(qJ(1));
t54 = t43 * t26;
t25 = sin(t30);
t41 = sin(qJ(1));
t61 = t41 * t25;
t5 = t39 * t61 + t54;
t55 = t43 * t25;
t60 = t41 * t26;
t6 = -t39 * t60 + t55;
t68 = -t5 * rSges(6,1) + t6 * rSges(6,2);
t7 = -t39 * t55 + t60;
t8 = t39 * t54 + t61;
t67 = t7 * rSges(6,1) - t8 * rSges(6,2);
t38 = sin(pkin(8));
t66 = g(3) * t38;
t40 = sin(qJ(3));
t65 = t40 * pkin(3);
t42 = cos(qJ(3));
t34 = t42 * pkin(3);
t64 = pkin(2) * t39 + pkin(1);
t63 = rSges(3,2) * t38;
t62 = t39 * t43;
t28 = sin(t37);
t59 = t41 * t28;
t29 = cos(t37);
t58 = t41 * t29;
t57 = t41 * t40;
t56 = t41 * t42;
t53 = t43 * t28;
t52 = t43 * t29;
t51 = t43 * t40;
t50 = t43 * t42;
t49 = qJ(4) + pkin(6);
t21 = pkin(4) * t29 + t34;
t31 = t41 * qJ(2);
t48 = t43 * pkin(1) + t31;
t47 = t64 + (rSges(4,3) + pkin(6)) * t38;
t46 = t39 * t34 + t64 + (rSges(5,3) + t49) * t38;
t45 = -rSges(6,1) * t25 - rSges(6,2) * t26;
t17 = -t39 * t51 + t56;
t15 = t39 * t57 + t50;
t44 = (pkin(2) + t21) * t39 + (rSges(6,3) + pkin(7) + t49) * t38;
t32 = t43 * qJ(2);
t27 = qJ(2) + t65;
t20 = pkin(4) * t28 + t65;
t18 = t39 * t50 + t57;
t16 = -t39 * t56 + t51;
t12 = t39 * t52 + t59;
t11 = -t39 * t53 + t58;
t10 = -t39 * t58 + t53;
t9 = t39 * t59 + t52;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t41 - rSges(2,2) * t43) + g(2) * (rSges(2,1) * t43 - rSges(2,2) * t41)) - m(3) * (g(1) * (t43 * rSges(3,3) + t32) + g(2) * (rSges(3,1) * t62 - t43 * t63 + t48) + (g(1) * (-rSges(3,1) * t39 - pkin(1) + t63) + g(2) * rSges(3,3)) * t41) - m(4) * (g(1) * (t16 * rSges(4,1) + t15 * rSges(4,2) - t41 * t47 + t32) + g(2) * (t18 * rSges(4,1) + t17 * rSges(4,2) + t43 * t47 + t31)) - m(5) * (g(1) * (t10 * rSges(5,1) + t9 * rSges(5,2) + t27 * t43 - t41 * t46) + g(2) * (t12 * rSges(5,1) + t11 * rSges(5,2) + t27 * t41 + t43 * t46)) - m(6) * (g(1) * (t6 * rSges(6,1) + t5 * rSges(6,2) + t32) + g(2) * (t8 * rSges(6,1) + t7 * rSges(6,2) + t48) + (g(1) * t20 + g(2) * t44) * t43 + (g(1) * (-pkin(1) - t44) + g(2) * t20) * t41), (-m(3) - m(4) + t69) * (g(1) * t41 - g(2) * t43), -m(4) * (g(1) * (rSges(4,1) * t17 - rSges(4,2) * t18) + g(2) * (-rSges(4,1) * t15 + rSges(4,2) * t16)) - m(5) * (g(1) * (t11 * rSges(5,1) - t12 * rSges(5,2) + pkin(3) * t17) + g(2) * (-t9 * rSges(5,1) + t10 * rSges(5,2) - pkin(3) * t15)) - m(6) * (g(1) * (-t20 * t62 + t21 * t41 + t67) + g(2) * (-t20 * t39 * t41 - t21 * t43 + t68)) + (-m(4) * (-rSges(4,1) * t40 - rSges(4,2) * t42) - m(5) * (-rSges(5,1) * t28 - rSges(5,2) * t29 - t65) - m(6) * (-t20 + t45)) * t66, t69 * (-g(3) * t39 + (g(1) * t43 + g(2) * t41) * t38), -m(6) * (g(1) * t67 + g(2) * t68 + t45 * t66)];
taug = t1(:);
