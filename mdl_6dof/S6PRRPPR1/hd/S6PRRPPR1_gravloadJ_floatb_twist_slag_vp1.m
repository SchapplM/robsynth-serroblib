% Calculate Gravitation load on the joints for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPPR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPPR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:58:33
% EndTime: 2019-03-08 20:58:35
% DurationCPUTime: 0.92s
% Computational Cost: add. (541->126), mult. (888->181), div. (0->0), fcn. (1023->14), ass. (0->63)
t39 = sin(qJ(2));
t41 = cos(qJ(2));
t62 = sin(pkin(10));
t64 = cos(pkin(6));
t50 = t64 * t62;
t63 = cos(pkin(10));
t16 = -t39 * t50 + t63 * t41;
t38 = sin(qJ(3));
t40 = cos(qJ(3));
t34 = sin(pkin(6));
t58 = t34 * t62;
t92 = -t16 * t38 + t40 * t58;
t69 = t34 * t39;
t91 = -t38 * t69 + t64 * t40;
t68 = t34 * t41;
t90 = g(3) * t68;
t51 = t64 * t63;
t14 = t39 * t51 + t62 * t41;
t59 = t34 * t63;
t45 = -t14 * t38 - t40 * t59;
t44 = t45 * pkin(3);
t66 = rSges(7,3) + pkin(9) + qJ(5);
t33 = sin(pkin(12));
t35 = cos(pkin(12));
t89 = t33 * rSges(6,1) + t35 * rSges(6,2);
t15 = t63 * t39 + t41 * t50;
t13 = t62 * t39 - t41 * t51;
t79 = g(2) * t13;
t88 = g(1) * t15 + t79;
t65 = rSges(6,3) + qJ(5);
t32 = qJ(3) + pkin(11);
t28 = sin(t32);
t30 = cos(t32);
t3 = t14 * t28 + t30 * t59;
t5 = t16 * t28 - t30 * t58;
t9 = t28 * t69 - t64 * t30;
t87 = g(1) * t5 + g(2) * t3 + g(3) * t9;
t31 = pkin(12) + qJ(6);
t27 = sin(t31);
t29 = cos(t31);
t47 = rSges(7,1) * t29 - rSges(7,2) * t27 + pkin(5) * t35 + pkin(4);
t86 = t66 * t28 + t47 * t30;
t48 = rSges(6,1) * t35 - rSges(6,2) * t33 + pkin(4);
t85 = t65 * t28 + t48 * t30;
t81 = -m(6) - m(7);
t78 = g(2) * t14;
t26 = pkin(3) * t40 + pkin(2);
t77 = t26 * t90;
t76 = g(3) * t34;
t75 = rSges(4,3) + pkin(8);
t36 = -qJ(4) - pkin(8);
t74 = -t13 * t26 - t14 * t36;
t73 = -t15 * t26 - t16 * t36;
t61 = -m(5) + t81;
t56 = t92 * pkin(3);
t53 = rSges(5,1) * t30 - rSges(5,2) * t28;
t52 = t91 * pkin(3);
t49 = rSges(4,1) * t40 - rSges(4,2) * t38 + pkin(2);
t46 = t27 * rSges(7,1) + t29 * rSges(7,2) + t33 * pkin(5);
t10 = t64 * t28 + t30 * t69;
t6 = t16 * t30 + t28 * t58;
t4 = t14 * t30 - t28 * t59;
t1 = [(-m(2) - m(3) - m(4) + t61) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t15 - rSges(3,2) * t16) + g(2) * (-rSges(3,1) * t13 - rSges(3,2) * t14) + (rSges(3,1) * t41 - rSges(3,2) * t39) * t76) - m(4) * (g(1) * (-t49 * t15 + t75 * t16) + t75 * t78 - t49 * t79 + (t75 * t39 + t49 * t41) * t76) - m(5) * (g(1) * (rSges(5,3) * t16 - t53 * t15 + t73) + g(2) * (rSges(5,3) * t14 - t53 * t13 + t74) + t77 + (t53 * t41 + (rSges(5,3) - t36) * t39) * t76) - m(6) * (g(1) * (t89 * t16 + t73) + g(2) * (t89 * t14 + t74) + t77 + ((-t36 + t89) * t39 + t85 * t41) * t76 - t88 * t85) - m(7) * (g(2) * t74 + t77 + t46 * t78 + ((-t36 + t46) * t39 + t86 * t41) * t76 - t88 * t86 + (t46 * t16 + t73) * g(1)) -m(4) * (g(1) * (t92 * rSges(4,1) + (-t16 * t40 - t38 * t58) * rSges(4,2)) + g(2) * (t45 * rSges(4,1) + (-t14 * t40 + t38 * t59) * rSges(4,2)) + g(3) * (t91 * rSges(4,1) + (-t64 * t38 - t40 * t69) * rSges(4,2))) - m(5) * (g(1) * (-rSges(5,1) * t5 - rSges(5,2) * t6 + t56) + g(2) * (-t3 * rSges(5,1) - t4 * rSges(5,2) + t44) + g(3) * (-rSges(5,1) * t9 - rSges(5,2) * t10 + t52)) + (-g(1) * (t6 * t66 + t56) - g(2) * (t4 * t66 + t44) - g(3) * (t10 * t66 + t52) + t87 * t47) * m(7) + (-g(1) * (t6 * t65 + t56) - g(2) * (t4 * t65 + t44) - g(3) * (t10 * t65 + t52) + t87 * t48) * m(6), t61 * (t88 - t90) t81 * t87, -m(7) * (g(1) * ((t15 * t29 - t27 * t6) * rSges(7,1) + (-t15 * t27 - t29 * t6) * rSges(7,2)) + g(2) * ((t13 * t29 - t27 * t4) * rSges(7,1) + (-t13 * t27 - t29 * t4) * rSges(7,2)) + g(3) * ((-t10 * t27 - t29 * t68) * rSges(7,1) + (-t10 * t29 + t27 * t68) * rSges(7,2)))];
taug  = t1(:);
