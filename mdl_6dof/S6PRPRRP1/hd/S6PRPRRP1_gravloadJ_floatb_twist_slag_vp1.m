% Calculate Gravitation load on the joints for
% S6PRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:55:16
% EndTime: 2019-03-08 19:55:19
% DurationCPUTime: 0.81s
% Computational Cost: add. (599->136), mult. (1497->216), div. (0->0), fcn. (1871->12), ass. (0->79)
t45 = sin(pkin(11));
t52 = sin(qJ(2));
t55 = cos(qJ(2));
t82 = cos(pkin(11));
t34 = -t55 * t45 - t52 * t82;
t46 = sin(pkin(10));
t47 = cos(pkin(10));
t48 = cos(pkin(6));
t60 = -t52 * t45 + t55 * t82;
t57 = t60 * t48;
t23 = t34 * t47 - t46 * t57;
t89 = t48 * t55;
t62 = -t46 * t89 - t47 * t52;
t59 = t62 * pkin(2);
t108 = t23 * pkin(3) + t59;
t107 = -t46 * t52 + t47 * t89;
t81 = sin(pkin(6));
t64 = t82 * t81;
t69 = t55 * t81;
t32 = t45 * t69 + t52 * t64;
t51 = sin(qJ(4));
t54 = cos(qJ(4));
t26 = t32 * t54 + t48 * t51;
t106 = g(3) * t26;
t84 = t34 * t48;
t24 = t46 * t84 + t47 * t60;
t19 = -t46 * t60 + t47 * t84;
t70 = t54 * t81;
t11 = -t19 * t51 + t47 * t70;
t13 = t24 * t51 - t46 * t70;
t25 = t32 * t51 - t48 * t54;
t105 = g(1) * t13 + g(2) * t11 + g(3) * t25;
t104 = pkin(4) * t54;
t71 = t52 * t81;
t31 = t45 * t71 - t55 * t64;
t101 = g(3) * t31;
t100 = rSges(5,3) + pkin(8);
t99 = rSges(6,3) + pkin(9);
t50 = sin(qJ(5));
t98 = t19 * t50;
t97 = t24 * t50;
t96 = t32 * t50;
t53 = cos(qJ(5));
t44 = pkin(5) * t53 + pkin(4);
t95 = t44 * t54;
t90 = t48 * t52;
t88 = t50 * t54;
t87 = t53 * t54;
t85 = rSges(7,3) + qJ(6) + pkin(9);
t41 = pkin(2) * t69;
t83 = -t31 * pkin(3) + t41;
t78 = g(1) * t99;
t77 = g(2) * t99;
t76 = -m(4) - m(5) - m(6) - m(7);
t75 = g(1) * t85;
t74 = g(2) * t85;
t73 = t51 * t81;
t68 = t107 * pkin(2);
t67 = pkin(8) * t32 + t83;
t66 = rSges(5,1) * t54 - rSges(5,2) * t51;
t12 = -t19 * t54 - t47 * t73;
t20 = t46 * t34 + t47 * t57;
t1 = -t12 * t50 - t20 * t53;
t14 = t24 * t54 + t46 * t73;
t3 = -t14 * t50 - t23 * t53;
t9 = -t26 * t50 + t31 * t53;
t65 = t20 * pkin(3) + t68;
t58 = -pkin(8) * t19 + t65;
t56 = pkin(8) * t24 + t108;
t16 = -t31 * t87 + t96;
t15 = t31 * t88 + t32 * t53;
t10 = -t26 * t53 - t31 * t50;
t8 = t23 * t87 + t97;
t7 = -t23 * t88 + t24 * t53;
t6 = t20 * t87 - t98;
t5 = -t19 * t53 - t20 * t88;
t4 = -t14 * t53 + t23 * t50;
t2 = -t12 * t53 + t20 * t50;
t17 = [(-m(2) - m(3) + t76) * g(3), -m(3) * (g(1) * (t62 * rSges(3,1) + (t46 * t90 - t47 * t55) * rSges(3,2)) + g(2) * (t107 * rSges(3,1) + (-t46 * t55 - t47 * t90) * rSges(3,2)) + g(3) * (rSges(3,1) * t69 - rSges(3,2) * t71)) - m(4) * (g(1) * (t23 * rSges(4,1) - rSges(4,2) * t24 + t59) + g(2) * (rSges(4,1) * t20 + rSges(4,2) * t19 + t68) + g(3) * (-rSges(4,1) * t31 - rSges(4,2) * t32 + t41)) - m(5) * (g(1) * (t100 * t24 + t23 * t66 + t108) + g(2) * (-t100 * t19 + t20 * t66 + t65) + g(3) * (t100 * t32 - t31 * t66 + t83)) - m(6) * (g(1) * (t8 * rSges(6,1) + t7 * rSges(6,2) + t104 * t23 + t56) + g(2) * (rSges(6,1) * t6 + rSges(6,2) * t5 + t104 * t20 + t58) + g(3) * (rSges(6,1) * t16 + rSges(6,2) * t15 - t104 * t31 + t67) + (-t101 * t99 + t20 * t77 + t23 * t78) * t51) - m(7) * (g(1) * (t8 * rSges(7,1) + t7 * rSges(7,2) + pkin(5) * t97 + t23 * t95 + t56) + g(2) * (rSges(7,1) * t6 + rSges(7,2) * t5 - pkin(5) * t98 + t20 * t95 + t58) + g(3) * (rSges(7,1) * t16 + rSges(7,2) * t15 + pkin(5) * t96 - t31 * t95 + t67) + (-t101 * t85 + t20 * t74 + t23 * t75) * t51) t76 * (g(3) * t48 + (g(1) * t46 - g(2) * t47) * t81) -m(5) * (g(1) * (-rSges(5,1) * t13 - rSges(5,2) * t14) + g(2) * (-rSges(5,1) * t11 - rSges(5,2) * t12) + g(3) * (-rSges(5,1) * t25 - rSges(5,2) * t26)) - m(6) * (t99 * t106 + t12 * t77 + t14 * t78 + t105 * (-rSges(6,1) * t53 + rSges(6,2) * t50 - pkin(4))) - m(7) * (t85 * t106 + t12 * t74 + t14 * t75 + t105 * (-rSges(7,1) * t53 + rSges(7,2) * t50 - t44)) -m(6) * (g(1) * (rSges(6,1) * t3 + rSges(6,2) * t4) + g(2) * (rSges(6,1) * t1 + rSges(6,2) * t2) + g(3) * (rSges(6,1) * t9 + rSges(6,2) * t10)) + (-g(1) * (rSges(7,1) * t3 + rSges(7,2) * t4) - g(2) * (rSges(7,1) * t1 + rSges(7,2) * t2) - g(3) * (rSges(7,1) * t9 + rSges(7,2) * t10) - (g(1) * t3 + g(2) * t1 + g(3) * t9) * pkin(5)) * m(7), -m(7) * t105];
taug  = t17(:);
