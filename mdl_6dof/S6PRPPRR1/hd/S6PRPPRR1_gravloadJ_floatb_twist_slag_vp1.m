% Calculate Gravitation load on the joints for
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
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
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPPRR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPPRR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:13:59
% EndTime: 2019-03-08 19:14:01
% DurationCPUTime: 0.75s
% Computational Cost: add. (507->115), mult. (1061->182), div. (0->0), fcn. (1304->14), ass. (0->60)
t39 = sin(pkin(11));
t46 = sin(qJ(2));
t48 = cos(qJ(2));
t72 = cos(pkin(11));
t26 = -t48 * t39 - t46 * t72;
t40 = sin(pkin(10));
t42 = cos(pkin(10));
t43 = cos(pkin(6));
t75 = t43 * t48;
t94 = -t40 * t46 + t42 * t75;
t52 = -t46 * t39 + t48 * t72;
t49 = t52 * t43;
t14 = t42 * t26 - t40 * t49;
t73 = t26 * t43;
t15 = t40 * t73 + t42 * t52;
t41 = cos(pkin(12));
t34 = t41 * pkin(4) + pkin(3);
t44 = -pkin(8) - qJ(4);
t53 = -t40 * t75 - t42 * t46;
t92 = t53 * pkin(2);
t93 = t14 * t34 - t15 * t44 + t92;
t71 = sin(pkin(6));
t57 = t72 * t71;
t61 = t46 * t71;
t23 = t39 * t61 - t48 * t57;
t91 = g(1) * t14 - g(3) * t23;
t90 = rSges(5,3) + qJ(4);
t11 = t40 * t26 + t42 * t49;
t89 = -g(2) * t11 - t91;
t10 = -t40 * t52 + t42 * t73;
t37 = pkin(12) + qJ(5);
t36 = cos(t37);
t84 = t36 * pkin(5);
t83 = rSges(7,3) + pkin(9);
t45 = sin(qJ(6));
t82 = t36 * t45;
t47 = cos(qJ(6));
t81 = t36 * t47;
t76 = t43 * t46;
t70 = -m(5) - m(6) - m(7);
t60 = t48 * t71;
t24 = t39 * t60 + t46 * t57;
t32 = pkin(2) * t60;
t67 = -t23 * t34 - t24 * t44 + t32;
t66 = g(2) * t83;
t65 = -m(4) + t70;
t64 = t40 * t71;
t63 = t42 * t71;
t59 = t94 * pkin(2);
t35 = sin(t37);
t58 = rSges(6,1) * t36 - rSges(6,2) * t35;
t55 = rSges(7,1) * t47 - rSges(7,2) * t45 + pkin(5);
t54 = t10 * t44 + t11 * t34 + t59;
t17 = t24 * t36 + t43 * t35;
t16 = -t24 * t35 + t43 * t36;
t5 = t15 * t36 + t35 * t64;
t4 = -t15 * t35 + t36 * t64;
t3 = -t10 * t36 - t35 * t63;
t2 = t10 * t35 - t36 * t63;
t1 = [(-m(2) - m(3) + t65) * g(3), -m(3) * (g(1) * (t53 * rSges(3,1) + (t40 * t76 - t42 * t48) * rSges(3,2)) + g(2) * (t94 * rSges(3,1) + (-t40 * t48 - t42 * t76) * rSges(3,2)) + g(3) * (rSges(3,1) * t60 - rSges(3,2) * t61)) - m(4) * (g(1) * (t14 * rSges(4,1) - rSges(4,2) * t15 + t92) + g(2) * (t11 * rSges(4,1) + t10 * rSges(4,2) + t59) + g(3) * (-t23 * rSges(4,1) - t24 * rSges(4,2) + t32)) - m(6) * (g(1) * (rSges(6,3) * t15 + t58 * t14 + t93) + g(2) * (-t10 * rSges(6,3) + t58 * t11 + t54) + g(3) * (t24 * rSges(6,3) - t58 * t23 + t67)) - m(7) * (g(1) * (t14 * t84 + (t14 * t81 + t15 * t45) * rSges(7,1) + (-t14 * t82 + t15 * t47) * rSges(7,2) + t93) + g(2) * (t11 * t84 + (-t10 * t45 + t11 * t81) * rSges(7,1) + (-t10 * t47 - t11 * t82) * rSges(7,2) + t54) + g(3) * (-t23 * t84 + (-t23 * t81 + t24 * t45) * rSges(7,1) + (t23 * t82 + t24 * t47) * rSges(7,2) + t67) + (t11 * t66 + t83 * t91) * t35) + (-g(1) * (t15 * t90 + t92) - g(2) * (-t10 * t90 + t59) - g(3) * (t24 * t90 + t32) + t89 * (rSges(5,1) * t41 - rSges(5,2) * sin(pkin(12)) + pkin(3))) * m(5), t65 * (g(1) * t64 - g(2) * t63 + g(3) * t43) t70 * t89, -m(6) * (g(1) * (t4 * rSges(6,1) - t5 * rSges(6,2)) + g(2) * (t2 * rSges(6,1) - t3 * rSges(6,2)) + g(3) * (t16 * rSges(6,1) - t17 * rSges(6,2))) - m(7) * (g(1) * (t55 * t4 + t83 * t5) + t3 * t66 + g(2) * t55 * t2 + (t55 * t16 + t83 * t17) * g(3)) -m(7) * (g(1) * ((-t14 * t47 - t5 * t45) * rSges(7,1) + (t14 * t45 - t5 * t47) * rSges(7,2)) + g(2) * ((-t11 * t47 - t3 * t45) * rSges(7,1) + (t11 * t45 - t3 * t47) * rSges(7,2)) + g(3) * ((-t17 * t45 + t23 * t47) * rSges(7,1) + (-t17 * t47 - t23 * t45) * rSges(7,2)))];
taug  = t1(:);
