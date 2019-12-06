% Calculate Gravitation load on the joints for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-05 19:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 19:00:03
% EndTime: 2019-12-05 19:00:04
% DurationCPUTime: 0.32s
% Computational Cost: add. (328->86), mult. (245->106), div. (0->0), fcn. (194->10), ass. (0->57)
t75 = rSges(4,3) + pkin(7);
t39 = -pkin(8) - pkin(7);
t74 = rSges(5,3) - t39;
t73 = rSges(6,3) + pkin(9) - t39;
t37 = cos(qJ(3));
t57 = t37 * rSges(4,1);
t72 = -pkin(2) - t57;
t34 = qJ(1) + qJ(2);
t27 = sin(t34);
t33 = qJ(3) + qJ(4);
t30 = qJ(5) + t33;
t23 = sin(t30);
t61 = t23 * t27;
t24 = cos(t30);
t63 = rSges(6,2) * t24;
t71 = rSges(6,1) * t61 + t27 * t63;
t26 = sin(t33);
t59 = t26 * t27;
t28 = cos(t33);
t64 = rSges(5,2) * t28;
t70 = rSges(5,1) * t59 + t27 * t64;
t69 = pkin(4) * t26;
t29 = cos(t34);
t68 = g(3) * t29;
t35 = sin(qJ(3));
t67 = t35 * pkin(3);
t36 = sin(qJ(1));
t66 = t36 * pkin(1);
t38 = cos(qJ(1));
t65 = t38 * pkin(1);
t31 = t37 * pkin(3);
t25 = t31 + pkin(2);
t62 = t23 * rSges(6,2);
t15 = t24 * rSges(6,1);
t60 = t26 * rSges(5,2);
t17 = t28 * rSges(5,1);
t58 = t35 * rSges(4,2);
t56 = t27 * t58 + t75 * t29;
t22 = pkin(4) * t28;
t55 = -t22 - t25 - t15;
t54 = -t29 * rSges(3,1) + t27 * rSges(3,2);
t53 = -t25 - t17;
t52 = t17 - t60;
t51 = t15 - t62;
t50 = t22 + t51;
t49 = -t27 * rSges(3,1) - t29 * rSges(3,2);
t48 = rSges(4,1) * t35 + rSges(4,2) * t37;
t47 = -rSges(5,1) * t26 - t64;
t46 = -rSges(6,1) * t23 - t63;
t45 = (t58 + t72) * t29;
t44 = (t55 + t62) * t29 - t73 * t27;
t43 = rSges(6,2) * t61 + t55 * t27 + t73 * t29;
t42 = rSges(5,2) * t59 + t53 * t27 + t74 * t29;
t41 = (t53 + t60) * t29 - t74 * t27;
t40 = (-g(2) * t75 + g(3) * t72) * t27;
t6 = -t67 - t69;
t1 = [-m(2) * (g(2) * (-t38 * rSges(2,1) + t36 * rSges(2,2)) + g(3) * (-t36 * rSges(2,1) - t38 * rSges(2,2))) - m(3) * (g(2) * (t54 - t65) + g(3) * (t49 - t66)) - m(4) * (g(2) * (t45 - t65) + g(3) * (t56 - t66) + t40) - m(5) * (g(2) * (t41 - t65) + g(3) * (t42 - t66)) - m(6) * (g(2) * (t44 - t65) + g(3) * (t43 - t66)), -m(3) * (g(2) * t54 + g(3) * t49) - m(4) * (g(2) * t45 + g(3) * t56 + t40) - m(5) * (g(2) * t41 + g(3) * t42) - m(6) * (g(2) * t44 + g(3) * t43), -m(4) * (g(1) * (t57 - t58) + g(2) * t48 * t27) - m(5) * (g(1) * (t31 + t52) + g(2) * (t27 * t67 + t70)) - m(6) * (g(1) * (t31 + t50) + g(2) * (-t27 * t6 + t71)) + (m(4) * t48 - m(5) * (t47 - t67) - m(6) * (t46 + t6)) * t68, -m(5) * (g(1) * t52 + g(2) * t70) - m(6) * (g(1) * t50 + g(2) * (pkin(4) * t59 + t71)) + (-m(5) * t47 - m(6) * (t46 - t69)) * t68, -m(6) * (g(1) * t51 + g(2) * t71 + t46 * t68)];
taug = t1(:);
