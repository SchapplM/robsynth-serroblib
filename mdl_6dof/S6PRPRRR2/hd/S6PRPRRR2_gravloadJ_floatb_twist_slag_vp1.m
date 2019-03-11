% Calculate Gravitation load on the joints for
% S6PRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:26:46
% EndTime: 2019-03-08 20:26:49
% DurationCPUTime: 1.10s
% Computational Cost: add. (674->147), mult. (1583->239), div. (0->0), fcn. (1984->14), ass. (0->74)
t45 = sin(qJ(2));
t77 = sin(pkin(12));
t79 = cos(pkin(12));
t96 = cos(qJ(2));
t28 = -t45 * t79 - t96 * t77;
t40 = sin(pkin(11));
t41 = cos(pkin(11));
t27 = t45 * t77 - t96 * t79;
t42 = cos(pkin(6));
t51 = t42 * t27;
t17 = t28 * t41 + t40 * t51;
t73 = t42 * t96;
t54 = -t40 * t73 - t41 * t45;
t52 = t54 * pkin(2);
t109 = t17 * pkin(3) + t52;
t108 = -t40 * t45 + t41 * t73;
t78 = sin(pkin(6));
t57 = t78 * t77;
t58 = t79 * t78;
t25 = t45 * t57 - t96 * t58;
t107 = g(1) * t17 - g(3) * t25;
t81 = t28 * t42;
t18 = -t41 * t27 + t40 * t81;
t47 = cos(qJ(4));
t44 = sin(qJ(4));
t72 = t44 * t78;
t10 = t18 * t47 + t40 * t72;
t26 = t45 * t58 + t96 * t57;
t20 = t26 * t47 + t42 * t44;
t106 = g(1) * t10 + g(3) * t20;
t13 = t40 * t27 + t41 * t81;
t19 = -t26 * t44 + t42 * t47;
t70 = t47 * t78;
t7 = t13 * t44 - t41 * t70;
t9 = -t18 * t44 + t40 * t70;
t105 = g(1) * t9 + g(2) * t7 + g(3) * t19;
t99 = t47 * pkin(4);
t98 = rSges(5,3) + pkin(8);
t97 = rSges(6,3) + pkin(9);
t43 = sin(qJ(5));
t95 = t13 * t43;
t94 = t18 * t43;
t93 = t26 * t43;
t39 = qJ(5) + qJ(6);
t37 = sin(t39);
t92 = t37 * t47;
t38 = cos(t39);
t91 = t38 * t47;
t86 = t42 * t45;
t85 = t43 * t47;
t46 = cos(qJ(5));
t84 = t46 * t47;
t36 = pkin(5) * t46 + pkin(4);
t83 = t47 * t36;
t82 = rSges(7,3) + pkin(10) + pkin(9);
t64 = t96 * t78;
t34 = pkin(2) * t64;
t80 = -t25 * pkin(3) + t34;
t76 = g(2) * t97;
t75 = -m(4) - m(5) - m(6) - m(7);
t74 = g(2) * t82;
t67 = t108 * pkin(2);
t66 = t26 * pkin(8) + t80;
t14 = t40 * t28 - t41 * t51;
t8 = -t13 * t47 - t41 * t72;
t65 = -t14 * t46 - t43 * t8;
t62 = rSges(5,1) * t47 - rSges(5,2) * t44;
t61 = -t10 * t43 - t17 * t46;
t60 = -t20 * t43 + t25 * t46;
t59 = t14 * pkin(3) + t67;
t53 = -t13 * pkin(8) + t59;
t50 = m(7) * (g(1) * ((-t10 * t37 - t17 * t38) * rSges(7,1) + (-t10 * t38 + t17 * t37) * rSges(7,2)) + g(2) * ((-t14 * t38 - t37 * t8) * rSges(7,1) + (t14 * t37 - t38 * t8) * rSges(7,2)) + g(3) * ((-t20 * t37 + t25 * t38) * rSges(7,1) + (-t20 * t38 - t25 * t37) * rSges(7,2)));
t49 = pkin(8) * t18 + t109;
t1 = [(-m(2) - m(3) + t75) * g(3), -m(3) * (g(1) * (t54 * rSges(3,1) + (t40 * t86 - t41 * t96) * rSges(3,2)) + g(2) * (t108 * rSges(3,1) + (-t40 * t96 - t41 * t86) * rSges(3,2)) + g(3) * (-rSges(3,2) * t45 * t78 + rSges(3,1) * t64)) - m(4) * (g(1) * (t17 * rSges(4,1) - rSges(4,2) * t18 + t52) + g(2) * (rSges(4,1) * t14 + rSges(4,2) * t13 + t67) + g(3) * (-rSges(4,1) * t25 - rSges(4,2) * t26 + t34)) - m(5) * (g(1) * (t62 * t17 + t18 * t98 + t109) + g(2) * (-t98 * t13 + t62 * t14 + t59) + g(3) * (-t62 * t25 + t98 * t26 + t80)) - m(6) * (g(1) * (t17 * t99 + (t17 * t84 + t94) * rSges(6,1) + (-t17 * t85 + t18 * t46) * rSges(6,2) + t49) + g(2) * (t14 * t99 + (t14 * t84 - t95) * rSges(6,1) + (-t13 * t46 - t14 * t85) * rSges(6,2) + t53) + g(3) * (-t25 * t99 + (-t25 * t84 + t93) * rSges(6,1) + (t25 * t85 + t26 * t46) * rSges(6,2) + t66) + (t107 * t97 + t14 * t76) * t44) - m(7) * (g(1) * (t17 * t83 + pkin(5) * t94 + (t17 * t91 + t18 * t37) * rSges(7,1) + (-t17 * t92 + t18 * t38) * rSges(7,2) + t49) + g(2) * (t14 * t83 - pkin(5) * t95 + (-t13 * t37 + t14 * t91) * rSges(7,1) + (-t13 * t38 - t14 * t92) * rSges(7,2) + t53) + g(3) * (-t25 * t83 + pkin(5) * t93 + (-t25 * t91 + t26 * t37) * rSges(7,1) + (t25 * t92 + t26 * t38) * rSges(7,2) + t66) + (t107 * t82 + t14 * t74) * t44) t75 * (g(3) * t42 + (g(1) * t40 - g(2) * t41) * t78) -m(5) * (g(1) * (rSges(5,1) * t9 - rSges(5,2) * t10) + g(2) * (rSges(5,1) * t7 - rSges(5,2) * t8) + g(3) * (rSges(5,1) * t19 - rSges(5,2) * t20)) - m(6) * (t8 * t76 + t106 * t97 + t105 * (rSges(6,1) * t46 - rSges(6,2) * t43 + pkin(4))) - m(7) * (t8 * t74 + t106 * t82 + t105 * (rSges(7,1) * t38 - rSges(7,2) * t37 + t36)) -m(6) * (g(1) * (t61 * rSges(6,1) + (-t10 * t46 + t17 * t43) * rSges(6,2)) + g(2) * (t65 * rSges(6,1) + (t14 * t43 - t46 * t8) * rSges(6,2)) + g(3) * (t60 * rSges(6,1) + (-t20 * t46 - t25 * t43) * rSges(6,2))) - t50 - m(7) * (g(1) * t61 + g(2) * t65 + g(3) * t60) * pkin(5), -t50];
taug  = t1(:);
