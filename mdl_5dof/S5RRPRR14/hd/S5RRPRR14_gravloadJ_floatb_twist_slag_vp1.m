% Calculate Gravitation load on the joints for
% S5RRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR14_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR14_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR14_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR14_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:35:47
% EndTime: 2019-12-31 20:35:50
% DurationCPUTime: 0.76s
% Computational Cost: add. (449->132), mult. (793->193), div. (0->0), fcn. (920->12), ass. (0->56)
t44 = sin(qJ(2));
t45 = sin(qJ(1));
t47 = cos(qJ(2));
t61 = cos(pkin(5));
t73 = cos(qJ(1));
t53 = t61 * t73;
t20 = t45 * t44 - t47 * t53;
t43 = sin(qJ(5));
t46 = cos(qJ(5));
t21 = t44 * t53 + t45 * t47;
t38 = pkin(10) + qJ(4);
t35 = sin(t38);
t36 = cos(t38);
t40 = sin(pkin(5));
t59 = t40 * t73;
t5 = t21 * t36 - t35 * t59;
t84 = -t20 * t46 + t43 * t5;
t83 = -t20 * t43 - t46 * t5;
t66 = t40 * t47;
t82 = g(3) * t66;
t81 = rSges(6,1) * t43 + rSges(6,2) * t46;
t56 = t45 * t61;
t22 = t73 * t44 + t47 * t56;
t77 = g(2) * t20;
t80 = g(1) * t22 + t77;
t50 = rSges(6,1) * t46 - rSges(6,2) * t43 + pkin(4);
t74 = rSges(6,3) + pkin(9);
t79 = t74 * t35 + t50 * t36;
t41 = cos(pkin(10));
t34 = t41 * pkin(3) + pkin(2);
t76 = t34 * t82;
t75 = g(3) * t40;
t68 = t40 * t44;
t67 = t40 * t45;
t42 = -pkin(8) - qJ(3);
t65 = -t20 * t34 - t21 * t42;
t23 = -t44 * t56 + t73 * t47;
t64 = -t22 * t34 - t23 * t42;
t63 = t73 * pkin(1) + pkin(7) * t67;
t62 = qJ(3) + rSges(4,3);
t39 = sin(pkin(10));
t60 = t39 * t67;
t58 = -t45 * pkin(1) + pkin(7) * t59;
t57 = -t21 * t35 - t36 * t59;
t55 = t39 * t59;
t54 = pkin(3) * t60 - t22 * t42 + t23 * t34 + t63;
t52 = rSges(5,1) * t36 - rSges(5,2) * t35;
t51 = rSges(4,1) * t41 - rSges(4,2) * t39 + pkin(2);
t49 = pkin(3) * t55 + t20 * t42 - t21 * t34 + t58;
t15 = t61 * t35 + t36 * t68;
t14 = -t35 * t68 + t61 * t36;
t9 = t23 * t36 + t35 * t67;
t8 = t23 * t35 - t36 * t67;
t2 = t22 * t43 + t9 * t46;
t1 = t22 * t46 - t9 * t43;
t3 = [-m(2) * (g(1) * (-t45 * rSges(2,1) - t73 * rSges(2,2)) + g(2) * (t73 * rSges(2,1) - t45 * rSges(2,2))) - m(3) * (g(1) * (-t21 * rSges(3,1) + t20 * rSges(3,2) + rSges(3,3) * t59 + t58) + g(2) * (rSges(3,1) * t23 - rSges(3,2) * t22 + rSges(3,3) * t67 + t63)) - m(4) * (g(1) * (-t21 * pkin(2) + (-t21 * t41 + t55) * rSges(4,1) + (t21 * t39 + t41 * t59) * rSges(4,2) - t62 * t20 + t58) + g(2) * (t23 * pkin(2) + (t23 * t41 + t60) * rSges(4,1) + (-t23 * t39 + t41 * t67) * rSges(4,2) + t62 * t22 + t63)) - m(5) * (g(1) * (-rSges(5,1) * t5 - rSges(5,2) * t57 - rSges(5,3) * t20 + t49) + g(2) * (rSges(5,1) * t9 - rSges(5,2) * t8 + rSges(5,3) * t22 + t54)) - m(6) * (g(1) * (t83 * rSges(6,1) + t84 * rSges(6,2) - t5 * pkin(4) + t74 * t57 + t49) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t9 * pkin(4) + t74 * t8 + t54)), -m(3) * (g(1) * (-rSges(3,1) * t22 - rSges(3,2) * t23) + g(2) * (-rSges(3,1) * t20 - rSges(3,2) * t21) + (rSges(3,1) * t47 - rSges(3,2) * t44) * t75) - m(4) * (g(1) * (-t51 * t22 + t62 * t23) + g(2) * t62 * t21 - t51 * t77 + (t62 * t44 + t51 * t47) * t75) - m(5) * (g(1) * (t23 * rSges(5,3) - t52 * t22 + t64) + g(2) * (t21 * rSges(5,3) - t52 * t20 + t65) + t76 + (t52 * t47 + (rSges(5,3) - t42) * t44) * t75) - m(6) * (g(1) * (t81 * t23 + t64) + g(2) * (t81 * t21 + t65) + t76 + ((-t42 + t81) * t44 + t79 * t47) * t75 - t80 * t79), (-m(4) - m(5) - m(6)) * (t80 - t82), -m(5) * (g(1) * (-rSges(5,1) * t8 - rSges(5,2) * t9) + g(2) * (rSges(5,1) * t57 - rSges(5,2) * t5) + g(3) * (rSges(5,1) * t14 - rSges(5,2) * t15)) - m(6) * (g(1) * (-t50 * t8 + t74 * t9) + (t50 * t14 + t74 * t15) * g(3) + (t74 * t5 + t50 * t57) * g(2)), -m(6) * (g(1) * (t1 * rSges(6,1) - t2 * rSges(6,2)) + g(2) * (-t84 * rSges(6,1) + t83 * rSges(6,2)) + g(3) * ((-t15 * t43 - t46 * t66) * rSges(6,1) + (-t15 * t46 + t43 * t66) * rSges(6,2)))];
taug = t3(:);
