% Calculate Gravitation load on the joints for
% S6PRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRPR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:15:46
% EndTime: 2019-03-08 23:15:48
% DurationCPUTime: 0.96s
% Computational Cost: add. (637->140), mult. (1198->207), div. (0->0), fcn. (1417->14), ass. (0->68)
t44 = sin(qJ(2));
t47 = cos(qJ(2));
t71 = cos(pkin(11));
t72 = cos(pkin(6));
t62 = t72 * t71;
t70 = sin(pkin(11));
t15 = t44 * t62 + t70 * t47;
t61 = t72 * t70;
t17 = -t44 * t61 + t71 * t47;
t99 = g(1) * t17 + g(2) * t15;
t14 = t70 * t44 - t47 * t62;
t16 = t71 * t44 + t47 * t61;
t98 = g(1) * t16 + g(2) * t14;
t43 = sin(qJ(3));
t46 = cos(qJ(3));
t40 = sin(pkin(6));
t67 = t40 * t70;
t10 = t17 * t43 - t46 * t67;
t76 = t40 * t44;
t18 = t43 * t76 - t72 * t46;
t68 = t40 * t71;
t8 = t15 * t43 + t46 * t68;
t97 = g(1) * t10 + g(2) * t8 + g(3) * t18;
t11 = t17 * t46 + t43 * t67;
t19 = t72 * t43 + t46 * t76;
t9 = t15 * t46 - t43 * t68;
t96 = g(1) * t11 + g(2) * t9 + g(3) * t19;
t39 = qJ(4) + pkin(12);
t35 = cos(t39);
t45 = cos(qJ(4));
t37 = t45 * pkin(4);
t22 = pkin(5) * t35 + t37;
t36 = qJ(6) + t39;
t31 = sin(t36);
t32 = cos(t36);
t56 = rSges(7,1) * t32 - rSges(7,2) * t31 + pkin(3) + t22;
t41 = -qJ(5) - pkin(9);
t73 = rSges(7,3) + pkin(10) - t41;
t95 = t73 * t43 + t56 * t46;
t34 = sin(t39);
t57 = rSges(6,1) * t35 - rSges(6,2) * t34 + pkin(3) + t37;
t74 = rSges(6,3) - t41;
t94 = t74 * t43 + t57 * t46;
t42 = sin(qJ(4));
t60 = rSges(5,1) * t45 - rSges(5,2) * t42 + pkin(3);
t77 = rSges(5,3) + pkin(9);
t93 = t77 * t43 + t60 * t46;
t92 = -m(6) - m(7);
t91 = (t14 * t32 - t31 * t9) * rSges(7,1) + (-t14 * t31 - t32 * t9) * rSges(7,2);
t90 = (-t11 * t31 + t16 * t32) * rSges(7,1) + (-t11 * t32 - t16 * t31) * rSges(7,2);
t75 = t40 * t47;
t89 = (-t19 * t31 - t32 * t75) * rSges(7,1) + (-t19 * t32 + t31 * t75) * rSges(7,2);
t80 = g(3) * t40;
t79 = t42 * pkin(4);
t78 = rSges(4,3) + pkin(8);
t69 = g(3) * (pkin(2) * t75 + pkin(8) * t76);
t66 = t14 * t45 - t42 * t9;
t65 = rSges(4,1) * t46 - rSges(4,2) * t43;
t64 = t42 * rSges(5,1) + t45 * rSges(5,2);
t63 = -t11 * t42 + t16 * t45;
t58 = -t19 * t42 - t45 * t75;
t21 = pkin(5) * t34 + t79;
t55 = t31 * rSges(7,1) + t32 * rSges(7,2) + t21;
t53 = t34 * rSges(6,1) + t35 * rSges(6,2) + t79;
t12 = t14 * pkin(2);
t13 = t16 * pkin(2);
t51 = -g(1) * t13 - g(2) * t12 + t69;
t1 = [(-m(2) - m(3) - m(4) - m(5) + t92) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t16 - rSges(3,2) * t17) + g(2) * (-rSges(3,1) * t14 - rSges(3,2) * t15) + (rSges(3,1) * t47 - rSges(3,2) * t44) * t80) - m(4) * (g(1) * (-t65 * t16 + t78 * t17 - t13) + g(2) * (-t65 * t14 + t78 * t15 - t12) + t69 + (rSges(4,3) * t44 + t65 * t47) * t80) - m(5) * ((t64 * t44 + t93 * t47) * t80 + t51 + t99 * (pkin(8) + t64) - t98 * t93) - m(6) * ((t53 * t44 + t94 * t47) * t80 + t51 + t99 * (pkin(8) + t53) - t98 * t94) - m(7) * ((t55 * t44 + t95 * t47) * t80 + t51 + t99 * (pkin(8) + t55) - t98 * t95) -m(4) * (g(1) * (-rSges(4,1) * t10 - rSges(4,2) * t11) + g(2) * (-rSges(4,1) * t8 - rSges(4,2) * t9) + g(3) * (-rSges(4,1) * t18 - rSges(4,2) * t19)) - m(5) * (-t97 * t60 + t96 * t77) - m(6) * (-t97 * t57 + t96 * t74) - m(7) * (-t97 * t56 + t96 * t73) -m(5) * (g(1) * (t63 * rSges(5,1) + (-t11 * t45 - t16 * t42) * rSges(5,2)) + g(2) * (t66 * rSges(5,1) + (-t14 * t42 - t45 * t9) * rSges(5,2)) + g(3) * (t58 * rSges(5,1) + (-t19 * t45 + t42 * t75) * rSges(5,2))) - m(7) * (g(1) * (-t11 * t21 + t16 * t22 + t90) + g(2) * (t14 * t22 - t21 * t9 + t91) + g(3) * (-t19 * t21 - t22 * t75 + t89)) + (-g(1) * ((-t11 * t34 + t16 * t35) * rSges(6,1) + (-t11 * t35 - t16 * t34) * rSges(6,2)) - g(2) * ((t14 * t35 - t34 * t9) * rSges(6,1) + (-t14 * t34 - t35 * t9) * rSges(6,2)) - g(3) * ((-t19 * t34 - t35 * t75) * rSges(6,1) + (-t19 * t35 + t34 * t75) * rSges(6,2)) - (g(1) * t63 + g(2) * t66 + g(3) * t58) * pkin(4)) * m(6), t92 * t97, -m(7) * (g(1) * t90 + g(2) * t91 + g(3) * t89)];
taug  = t1(:);
