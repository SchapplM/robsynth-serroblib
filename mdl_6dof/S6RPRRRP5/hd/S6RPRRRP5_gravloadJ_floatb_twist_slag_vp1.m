% Calculate Gravitation load on the joints for
% S6RPRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:10:46
% EndTime: 2019-03-09 06:10:48
% DurationCPUTime: 0.75s
% Computational Cost: add. (536->134), mult. (481->174), div. (0->0), fcn. (447->10), ass. (0->66)
t45 = cos(qJ(5));
t90 = rSges(7,1) + pkin(5);
t91 = -t45 * t90 - pkin(4);
t44 = sin(qJ(1));
t39 = pkin(10) + qJ(3);
t36 = qJ(4) + t39;
t29 = sin(t36);
t43 = sin(qJ(5));
t79 = t29 * t43;
t61 = rSges(6,2) * t79;
t30 = cos(t36);
t76 = t30 * t44;
t89 = rSges(6,3) * t76 + t44 * t61;
t46 = cos(qJ(1));
t74 = t30 * t46;
t88 = rSges(6,3) * t74 + t46 * t61;
t25 = t30 * rSges(5,1);
t55 = -rSges(5,2) * t29 + t25;
t66 = t30 * pkin(4) + t29 * pkin(9);
t87 = g(1) * t46 + g(2) * t44;
t63 = rSges(7,3) + qJ(6);
t86 = t87 * t29;
t34 = sin(t39);
t85 = pkin(3) * t34;
t42 = -pkin(7) - qJ(2);
t38 = -pkin(8) + t42;
t83 = g(2) * t38;
t41 = cos(pkin(10));
t31 = t41 * pkin(2) + pkin(1);
t24 = t29 * rSges(7,2);
t23 = t29 * rSges(6,3);
t78 = t29 * t46;
t77 = t30 * t43;
t75 = t30 * t45;
t73 = t38 * t46;
t72 = t44 * t43;
t71 = t44 * t45;
t70 = t45 * t46;
t69 = t46 * t43;
t68 = rSges(4,3) - t42;
t67 = rSges(5,3) - t38;
t65 = qJ(6) * t43;
t64 = rSges(3,3) + qJ(2);
t35 = cos(t39);
t28 = pkin(3) * t35;
t8 = t28 + t31;
t5 = t46 * t8;
t62 = pkin(4) * t74 + pkin(9) * t78 + t5;
t60 = -rSges(6,1) * t45 - pkin(4);
t17 = pkin(9) * t76;
t59 = -t44 * t85 + t17;
t20 = pkin(9) * t74;
t58 = -t46 * t85 + t20;
t57 = rSges(4,1) * t35 - rSges(4,2) * t34;
t54 = -t8 - t66;
t53 = rSges(3,1) * t41 - rSges(3,2) * sin(pkin(10)) + pkin(1);
t52 = t31 + t57;
t51 = rSges(7,3) * t77 + t30 * t65 + t90 * t75 + t24 + t66;
t50 = rSges(6,1) * t75 - rSges(6,2) * t77 + t23 + t66;
t16 = rSges(7,2) * t74;
t11 = rSges(7,2) * t76;
t4 = t30 * t70 + t72;
t3 = t30 * t69 - t71;
t2 = t30 * t71 - t69;
t1 = t30 * t72 + t70;
t6 = [-m(2) * (g(1) * (-t44 * rSges(2,1) - rSges(2,2) * t46) + g(2) * (rSges(2,1) * t46 - t44 * rSges(2,2))) - m(3) * ((g(1) * t64 + g(2) * t53) * t46 + (-g(1) * t53 + g(2) * t64) * t44) - m(4) * ((g(1) * t68 + g(2) * t52) * t46 + (-g(1) * t52 + g(2) * t68) * t44) - m(5) * (g(2) * t5 + (g(1) * t67 + g(2) * t55) * t46 + (g(1) * (-t55 - t8) + g(2) * t67) * t44) - m(6) * (g(1) * (-t2 * rSges(6,1) + t1 * rSges(6,2) - t73) + g(2) * (t4 * rSges(6,1) - t3 * rSges(6,2) + rSges(6,3) * t78 + t62) + (g(1) * (t54 - t23) - t83) * t44) - m(7) * (g(1) * (-t63 * t1 - t90 * t2 - t73) + g(2) * (rSges(7,2) * t78 + t63 * t3 + t4 * t90 + t62) + (g(1) * (t54 - t24) - t83) * t44) (-m(3) - m(4) - m(5) - m(6) - m(7)) * (g(1) * t44 - g(2) * t46) -m(4) * (g(3) * t57 + t87 * (-rSges(4,1) * t34 - rSges(4,2) * t35)) - m(5) * (g(3) * (t28 + t55) + t87 * (-rSges(5,1) * t29 - rSges(5,2) * t30 - t85)) - m(6) * (g(1) * (t58 + t88) + g(2) * (t59 + t89) + g(3) * (t28 + t50) + t60 * t86) - m(7) * (g(1) * (t16 + t58) + g(2) * (t11 + t59) + g(3) * (t28 + t51) + (-t63 * t43 + t91) * t86) -m(5) * (g(3) * t25 + (-g(1) * t74 - g(2) * t76) * rSges(5,2)) - m(6) * (g(1) * (t20 + t88) + g(2) * (t17 + t89) + g(3) * t50) - m(7) * (g(1) * (t16 + t20) + g(2) * (t11 + t17) + g(3) * t51) + (m(5) * g(3) * rSges(5,2) + t87 * (m(5) * rSges(5,1) - m(6) * t60 - m(7) * (-rSges(7,3) * t43 - t65 + t91))) * t29, -m(6) * (g(1) * (-rSges(6,1) * t3 - rSges(6,2) * t4) + g(2) * (-rSges(6,1) * t1 - rSges(6,2) * t2)) - m(7) * (g(1) * (-t3 * t90 + t63 * t4) + g(2) * (-t1 * t90 + t63 * t2)) + (-m(6) * (-rSges(6,1) * t43 - rSges(6,2) * t45) - m(7) * (-t90 * t43 + t63 * t45)) * g(3) * t29, -m(7) * (g(1) * t3 + g(2) * t1 + g(3) * t79)];
taug  = t6(:);
