% Calculate Gravitation load on the joints for
% S5PRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR8_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:14:49
% EndTime: 2019-12-05 17:14:52
% DurationCPUTime: 0.61s
% Computational Cost: add. (405->105), mult. (698->161), div. (0->0), fcn. (805->12), ass. (0->57)
t92 = rSges(6,3) + pkin(9);
t43 = sin(qJ(3));
t46 = cos(qJ(3));
t67 = cos(pkin(5));
t41 = sin(pkin(5));
t44 = sin(qJ(2));
t77 = t41 * t44;
t91 = -t43 * t77 + t67 * t46;
t47 = cos(qJ(2));
t40 = sin(pkin(10));
t63 = t40 * t67;
t66 = cos(pkin(10));
t29 = -t44 * t63 + t66 * t47;
t76 = t41 * t46;
t90 = -t29 * t43 + t40 * t76;
t42 = sin(qJ(5));
t45 = cos(qJ(5));
t89 = t45 * rSges(6,1) - t42 * rSges(6,2) + pkin(4);
t88 = t42 * rSges(6,1) + t45 * rSges(6,2);
t39 = qJ(3) + qJ(4);
t37 = sin(t39);
t38 = cos(t39);
t87 = t92 * t37 + t89 * t38;
t57 = t67 * t66;
t27 = t40 * t47 + t44 * t57;
t62 = t41 * t66;
t11 = -t27 * t37 - t38 * t62;
t12 = t27 * t38 - t37 * t62;
t86 = t11 * rSges(5,1) - t12 * rSges(5,2);
t78 = t40 * t41;
t13 = -t29 * t37 + t38 * t78;
t14 = t29 * t38 + t37 * t78;
t85 = t13 * rSges(5,1) - t14 * rSges(5,2);
t26 = t40 * t44 - t47 * t57;
t84 = g(2) * t26;
t36 = t46 * pkin(3) + pkin(2);
t75 = t41 * t47;
t83 = g(3) * t36 * t75;
t82 = g(3) * t41;
t81 = rSges(4,3) + pkin(7);
t48 = -pkin(8) - pkin(7);
t70 = -t26 * t36 - t27 * t48;
t28 = t66 * t44 + t47 * t63;
t69 = -t28 * t36 - t29 * t48;
t22 = -t37 * t77 + t67 * t38;
t23 = t67 * t37 + t38 * t77;
t68 = t22 * rSges(5,1) - t23 * rSges(5,2);
t60 = t90 * pkin(3);
t59 = rSges(5,1) * t38 - rSges(5,2) * t37;
t58 = t91 * pkin(3);
t56 = rSges(4,1) * t46 - rSges(4,2) * t43 + pkin(2);
t54 = -t27 * t43 - t46 * t62;
t53 = t89 * t11 + t92 * t12;
t52 = t89 * t13 + t92 * t14;
t51 = t89 * t22 + t92 * t23;
t50 = t54 * pkin(3);
t1 = [(-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3), -m(3) * (g(1) * (-t28 * rSges(3,1) - t29 * rSges(3,2)) + g(2) * (-t26 * rSges(3,1) - t27 * rSges(3,2)) + (rSges(3,1) * t47 - rSges(3,2) * t44) * t82) - m(4) * (g(1) * (-t56 * t28 + t81 * t29) + g(2) * t81 * t27 - t56 * t84 + (t81 * t44 + t56 * t47) * t82) - m(5) * (g(1) * (t29 * rSges(5,3) - t59 * t28 + t69) + g(2) * (t27 * rSges(5,3) - t59 * t26 + t70) + t83 + (t59 * t47 + (rSges(5,3) - t48) * t44) * t82) - m(6) * (g(2) * (t88 * t27 + t70) + t83 - t87 * t84 + ((-t48 + t88) * t44 + t87 * t47) * t82 + (-t28 * t87 + t88 * t29 + t69) * g(1)), -m(4) * (g(1) * (t90 * rSges(4,1) + (-t29 * t46 - t43 * t78) * rSges(4,2)) + g(2) * (t54 * rSges(4,1) + (-t27 * t46 + t43 * t62) * rSges(4,2)) + g(3) * (t91 * rSges(4,1) + (-t67 * t43 - t44 * t76) * rSges(4,2))) - m(5) * (g(1) * (t60 + t85) + g(2) * (t50 + t86) + g(3) * (t58 + t68)) - m(6) * (g(1) * (t52 + t60) + g(2) * (t50 + t53) + g(3) * (t51 + t58)), -m(5) * (g(1) * t85 + g(2) * t86 + g(3) * t68) - m(6) * (g(1) * t52 + g(2) * t53 + g(3) * t51), -m(6) * (g(1) * ((-t14 * t42 + t28 * t45) * rSges(6,1) + (-t14 * t45 - t28 * t42) * rSges(6,2)) + g(2) * ((-t12 * t42 + t26 * t45) * rSges(6,1) + (-t12 * t45 - t26 * t42) * rSges(6,2)) + g(3) * ((-t23 * t42 - t45 * t75) * rSges(6,1) + (-t23 * t45 + t42 * t75) * rSges(6,2)))];
taug = t1(:);
