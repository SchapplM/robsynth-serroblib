% Calculate Gravitation load on the joints for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPPRR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPPRR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:21:29
% EndTime: 2019-03-08 19:21:30
% DurationCPUTime: 0.53s
% Computational Cost: add. (402->118), mult. (982->185), div. (0->0), fcn. (1203->12), ass. (0->56)
t78 = rSges(7,3) + pkin(9);
t81 = g(1) * t78;
t51 = cos(qJ(5));
t80 = t51 * pkin(5);
t79 = rSges(6,3) + pkin(8);
t44 = sin(pkin(6));
t48 = sin(qJ(5));
t77 = t44 * t48;
t49 = sin(qJ(2));
t76 = t44 * t49;
t75 = t44 * t51;
t52 = cos(qJ(2));
t74 = t44 * t52;
t47 = sin(qJ(6));
t73 = t47 * t51;
t50 = cos(qJ(6));
t72 = t50 * t51;
t71 = pkin(2) * t74 + qJ(3) * t76;
t70 = rSges(4,3) + qJ(3);
t69 = cos(pkin(6));
t68 = -m(5) - m(6) - m(7);
t67 = pkin(3) * t74 + t71;
t66 = g(2) * t78;
t65 = g(3) * t78;
t64 = -m(4) + t68;
t43 = sin(pkin(10));
t46 = cos(pkin(10));
t60 = t52 * t69;
t30 = t43 * t49 - t46 * t60;
t61 = t49 * t69;
t31 = t43 * t52 + t46 * t61;
t42 = sin(pkin(11));
t45 = cos(pkin(11));
t63 = -t30 * t45 + t31 * t42;
t11 = t30 * t42 + t31 * t45;
t32 = t43 * t60 + t46 * t49;
t33 = -t43 * t61 + t46 * t52;
t62 = -t32 * t45 + t33 * t42;
t15 = t32 * t42 + t33 * t45;
t24 = (t42 * t49 + t45 * t52) * t44;
t59 = t24 * pkin(4) + t67;
t27 = t30 * pkin(2);
t58 = -t30 * pkin(3) + t31 * qJ(3) - t27;
t29 = t32 * pkin(2);
t57 = -t32 * pkin(3) + t33 * qJ(3) - t29;
t56 = rSges(6,1) * t51 - rSges(6,2) * t48;
t55 = pkin(4) * t63 + t58;
t54 = pkin(4) * t62 + t57;
t25 = -t42 * t74 + t45 * t76;
t17 = t25 * t51 - t69 * t48;
t16 = -t25 * t48 - t69 * t51;
t4 = t15 * t51 - t43 * t77;
t3 = -t15 * t48 - t43 * t75;
t2 = t11 * t51 + t46 * t77;
t1 = -t11 * t48 + t46 * t75;
t5 = [(-m(2) - m(3) + t64) * g(3), -m(3) * (g(1) * (-t32 * rSges(3,1) - t33 * rSges(3,2)) + g(2) * (-t30 * rSges(3,1) - t31 * rSges(3,2)) + g(3) * (rSges(3,1) * t52 - rSges(3,2) * t49) * t44) - m(4) * (g(1) * (-t32 * rSges(4,1) + t70 * t33 - t29) + g(2) * (-t30 * rSges(4,1) + t70 * t31 - t27) + g(3) * ((rSges(4,1) * t52 + rSges(4,3) * t49) * t44 + t71)) - m(5) * (g(1) * (rSges(5,1) * t62 + rSges(5,2) * t15 + t57) + g(2) * (rSges(5,1) * t63 + rSges(5,2) * t11 + t58) + g(3) * (t24 * rSges(5,1) + t25 * rSges(5,2) + t67)) - m(6) * (g(1) * (-t15 * t79 + t56 * t62 + t54) + g(2) * (-t11 * t79 + t56 * t63 + t55) + g(3) * (t56 * t24 - t79 * t25 + t59)) - m(7) * (g(1) * (t62 * t80 - t15 * pkin(8) + (-t15 * t47 + t62 * t72) * rSges(7,1) + (-t15 * t50 - t62 * t73) * rSges(7,2) + t54) + g(2) * (t63 * t80 - t11 * pkin(8) + (-t11 * t47 + t63 * t72) * rSges(7,1) + (-t11 * t50 - t63 * t73) * rSges(7,2) + t55) + g(3) * (t24 * t80 - t25 * pkin(8) + (t24 * t72 - t25 * t47) * rSges(7,1) + (-t24 * t73 - t25 * t50) * rSges(7,2) + t59) + (t24 * t65 + t62 * t81 + t63 * t66) * t48) t64 * (g(1) * t32 + g(2) * t30 - g(3) * t74) t68 * (-g(3) * t69 + (-g(1) * t43 + g(2) * t46) * t44) -m(6) * (g(1) * (t3 * rSges(6,1) - t4 * rSges(6,2)) + g(2) * (t1 * rSges(6,1) - t2 * rSges(6,2)) + g(3) * (t16 * rSges(6,1) - t17 * rSges(6,2))) - m(7) * (t4 * t81 + t17 * t65 + t2 * t66 + (g(1) * t3 + g(2) * t1 + g(3) * t16) * (rSges(7,1) * t50 - rSges(7,2) * t47 + pkin(5))) -m(7) * (g(1) * ((-t4 * t47 + t50 * t62) * rSges(7,1) + (-t4 * t50 - t47 * t62) * rSges(7,2)) + g(2) * ((-t2 * t47 + t50 * t63) * rSges(7,1) + (-t2 * t50 - t47 * t63) * rSges(7,2)) + g(3) * ((-t17 * t47 + t24 * t50) * rSges(7,1) + (-t17 * t50 - t24 * t47) * rSges(7,2)))];
taug  = t5(:);
