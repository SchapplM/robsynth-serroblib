% Calculate Gravitation load on the joints for
% S6RRPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 11:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR9_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR9_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR9_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR9_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:55:37
% EndTime: 2019-03-09 10:55:41
% DurationCPUTime: 1.26s
% Computational Cost: add. (756->178), mult. (1214->254), div. (0->0), fcn. (1424->14), ass. (0->72)
t51 = sin(qJ(2));
t52 = sin(qJ(1));
t74 = cos(pkin(6));
t90 = cos(qJ(2));
t62 = t74 * t90;
t91 = cos(qJ(1));
t21 = t51 * t52 - t91 * t62;
t42 = pkin(12) + qJ(6);
t37 = sin(t42);
t39 = cos(t42);
t68 = t51 * t74;
t22 = t52 * t90 + t91 * t68;
t43 = pkin(11) + qJ(4);
t38 = sin(t43);
t40 = cos(t43);
t46 = sin(pkin(6));
t71 = t46 * t91;
t6 = t22 * t40 - t38 * t71;
t105 = -t21 * t39 + t37 * t6;
t104 = -t21 * t37 - t39 * t6;
t77 = pkin(10) + qJ(5) + rSges(7,3);
t75 = qJ(5) + rSges(6,3);
t23 = t91 * t51 + t52 * t62;
t102 = g(1) * t23 + g(2) * t21;
t82 = t46 * t51;
t15 = t38 * t82 - t74 * t40;
t5 = t22 * t38 + t40 * t71;
t24 = -t52 * t68 + t91 * t90;
t81 = t46 * t52;
t9 = t24 * t38 - t40 * t81;
t101 = g(1) * t9 + g(2) * t5 + g(3) * t15;
t10 = t24 * t40 + t38 * t81;
t16 = t74 * t38 + t40 * t82;
t100 = g(1) * t10 + g(2) * t6 + g(3) * t16;
t97 = -m(6) - m(7);
t44 = sin(pkin(12));
t96 = pkin(5) * t44;
t92 = g(3) * t46;
t87 = rSges(5,2) * t38;
t84 = t21 * t44;
t83 = t23 * t44;
t48 = cos(pkin(11));
t36 = pkin(3) * t48 + pkin(2);
t50 = -pkin(9) - qJ(3);
t80 = -t21 * t36 - t22 * t50;
t79 = -t23 * t36 - t24 * t50;
t78 = t91 * pkin(1) + pkin(8) * t81;
t76 = qJ(3) + rSges(4,3);
t45 = sin(pkin(11));
t73 = t45 * t81;
t72 = t40 * t90;
t70 = t46 * t90;
t69 = -t52 * pkin(1) + pkin(8) * t71;
t67 = t38 * t70;
t66 = t40 * t70;
t65 = t45 * t71;
t25 = t36 * t70;
t64 = -t50 * t82 + t25;
t63 = pkin(3) * t73 - t23 * t50 + t24 * t36 + t78;
t61 = -rSges(5,1) * t40 + t87;
t47 = cos(pkin(12));
t60 = rSges(6,1) * t44 + rSges(6,2) * t47;
t58 = -rSges(6,1) * t47 + rSges(6,2) * t44 - pkin(4);
t35 = pkin(5) * t47 + pkin(4);
t57 = -rSges(7,1) * t39 + rSges(7,2) * t37 - t35;
t56 = pkin(3) * t65 + t21 * t50 - t22 * t36 + t69;
t55 = rSges(7,1) * t37 + rSges(7,2) * t39 + t96;
t54 = -t75 * t38 + t58 * t40;
t53 = -t77 * t38 + t57 * t40;
t3 = t10 * t39 + t23 * t37;
t2 = -t10 * t37 + t23 * t39;
t1 = [-m(2) * (g(1) * (-t52 * rSges(2,1) - t91 * rSges(2,2)) + g(2) * (t91 * rSges(2,1) - t52 * rSges(2,2))) - m(3) * (g(1) * (-t22 * rSges(3,1) + t21 * rSges(3,2) + rSges(3,3) * t71 + t69) + g(2) * (t24 * rSges(3,1) - t23 * rSges(3,2) + rSges(3,3) * t81 + t78)) - m(4) * (g(1) * (-t22 * pkin(2) + (-t22 * t48 + t65) * rSges(4,1) + (t22 * t45 + t48 * t71) * rSges(4,2) - t76 * t21 + t69) + g(2) * (t24 * pkin(2) + (t24 * t48 + t73) * rSges(4,1) + (-t24 * t45 + t48 * t81) * rSges(4,2) + t76 * t23 + t78)) - m(5) * (g(1) * (-rSges(5,1) * t6 + rSges(5,2) * t5 - t21 * rSges(5,3) + t56) + g(2) * (rSges(5,1) * t10 - rSges(5,2) * t9 + rSges(5,3) * t23 + t63)) - m(6) * (g(1) * (-t6 * pkin(4) + (-t47 * t6 - t84) * rSges(6,1) + (-t21 * t47 + t44 * t6) * rSges(6,2) - t75 * t5 + t56) + g(2) * (t10 * pkin(4) + (t10 * t47 + t83) * rSges(6,1) + (-t10 * t44 + t23 * t47) * rSges(6,2) + t75 * t9 + t63)) - m(7) * (g(1) * (t104 * rSges(7,1) + t105 * rSges(7,2) - pkin(5) * t84 - t6 * t35 - t77 * t5 + t56) + g(2) * (rSges(7,1) * t3 + rSges(7,2) * t2 + pkin(5) * t83 + t10 * t35 + t77 * t9 + t63)) -m(3) * (g(1) * (-rSges(3,1) * t23 - rSges(3,2) * t24) + g(2) * (-rSges(3,1) * t21 - rSges(3,2) * t22) + (t90 * rSges(3,1) - rSges(3,2) * t51) * t92) - m(4) * ((g(1) * t24 + g(2) * t22 + t51 * t92) * t76 + (t90 * t92 - t102) * (rSges(4,1) * t48 - rSges(4,2) * t45 + pkin(2))) - m(5) * (g(1) * (rSges(5,3) * t24 + t61 * t23 + t79) + g(2) * (rSges(5,3) * t22 + t61 * t21 + t80) + g(3) * t25 + (rSges(5,1) * t72 - t90 * t87 + (rSges(5,3) - t50) * t51) * t92) - m(6) * (g(1) * (t54 * t23 + t60 * t24 + t79) + g(2) * (t54 * t21 + t60 * t22 + t80) + g(3) * (pkin(4) * t66 + t64 + t75 * t67 + ((t44 * t51 + t47 * t72) * rSges(6,1) + (-t44 * t72 + t47 * t51) * rSges(6,2)) * t46)) - m(7) * (g(1) * (t53 * t23 + t55 * t24 + t79) + g(2) * (t53 * t21 + t55 * t22 + t80) + g(3) * (t35 * t66 + t82 * t96 + t64 + t77 * t67 + ((t37 * t51 + t39 * t72) * rSges(7,1) + (-t37 * t72 + t39 * t51) * rSges(7,2)) * t46)) (-m(4) - m(5) + t97) * (-g(3) * t70 + t102) -m(5) * (g(1) * (-rSges(5,1) * t9 - rSges(5,2) * t10) + g(2) * (-rSges(5,1) * t5 - rSges(5,2) * t6) + g(3) * (-rSges(5,1) * t15 - rSges(5,2) * t16)) - m(6) * (t100 * t75 + t101 * t58) - m(7) * (t100 * t77 + t101 * t57) t97 * t101, -m(7) * (g(1) * (rSges(7,1) * t2 - rSges(7,2) * t3) + g(2) * (-t105 * rSges(7,1) + t104 * rSges(7,2)) + g(3) * ((-t16 * t37 - t39 * t70) * rSges(7,1) + (-t16 * t39 + t37 * t70) * rSges(7,2)))];
taug  = t1(:);
