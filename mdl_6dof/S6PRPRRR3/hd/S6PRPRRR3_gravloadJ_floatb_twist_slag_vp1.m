% Calculate Gravitation load on the joints for
% S6PRPRRR3
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
% Datum: 2019-03-08 20:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:31:47
% EndTime: 2019-03-08 20:31:49
% DurationCPUTime: 0.67s
% Computational Cost: add. (621->117), mult. (812->174), div. (0->0), fcn. (928->14), ass. (0->66)
t105 = rSges(7,3) + pkin(10);
t46 = pkin(12) + qJ(4);
t41 = sin(t46);
t42 = cos(t46);
t76 = cos(pkin(6));
t49 = sin(pkin(6));
t53 = sin(qJ(2));
t83 = t49 * t53;
t104 = -t41 * t83 + t76 * t42;
t55 = cos(qJ(2));
t48 = sin(pkin(11));
t71 = t48 * t76;
t75 = cos(pkin(11));
t31 = -t53 * t71 + t75 * t55;
t84 = t48 * t49;
t103 = -t31 * t41 + t42 * t84;
t52 = sin(qJ(6));
t54 = cos(qJ(6));
t102 = rSges(7,1) * t54 - rSges(7,2) * t52 + pkin(5);
t82 = t49 * t55;
t101 = g(3) * t82;
t100 = rSges(7,1) * t52 + rSges(7,2) * t54;
t30 = t75 * t53 + t55 * t71;
t65 = t76 * t75;
t28 = t48 * t53 - t55 * t65;
t94 = g(2) * t28;
t99 = g(1) * t30 + t94;
t43 = qJ(5) + t46;
t38 = sin(t43);
t39 = cos(t43);
t98 = t102 * t39 + t105 * t38;
t29 = t48 * t55 + t53 * t65;
t70 = t49 * t75;
t12 = -t29 * t38 - t39 * t70;
t13 = t29 * t39 - t38 * t70;
t97 = t12 * rSges(6,1) - t13 * rSges(6,2);
t14 = -t31 * t38 + t39 * t84;
t15 = t31 * t39 + t38 * t84;
t96 = t14 * rSges(6,1) - t15 * rSges(6,2);
t93 = g(2) * t29;
t50 = cos(pkin(12));
t40 = t50 * pkin(3) + pkin(2);
t33 = pkin(4) * t42 + t40;
t92 = t33 * t101;
t91 = g(3) * t49;
t51 = -pkin(8) - qJ(3);
t81 = rSges(5,3) - t51;
t45 = -pkin(9) + t51;
t80 = -t28 * t33 - t29 * t45;
t79 = -t30 * t33 - t31 * t45;
t23 = -t38 * t83 + t76 * t39;
t24 = t76 * t38 + t39 * t83;
t78 = t23 * rSges(6,1) - t24 * rSges(6,2);
t77 = rSges(4,3) + qJ(3);
t72 = -m(4) - m(5) - m(6) - m(7);
t68 = t103 * pkin(4);
t67 = rSges(6,1) * t39 - rSges(6,2) * t38;
t66 = t104 * pkin(4);
t64 = rSges(4,1) * t50 - rSges(4,2) * sin(pkin(12)) + pkin(2);
t62 = rSges(5,1) * t42 - rSges(5,2) * t41 + t40;
t61 = -t29 * t41 - t42 * t70;
t60 = t102 * t12 + t105 * t13;
t59 = t102 * t14 + t105 * t15;
t58 = t102 * t23 + t105 * t24;
t57 = t61 * pkin(4);
t1 = [(-m(2) - m(3) + t72) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t30 - rSges(3,2) * t31) + g(2) * (-rSges(3,1) * t28 - rSges(3,2) * t29) + (rSges(3,1) * t55 - rSges(3,2) * t53) * t91) - m(4) * (g(1) * (-t64 * t30 + t77 * t31) + t77 * t93 - t64 * t94 + (t77 * t53 + t64 * t55) * t91) - m(5) * (g(1) * (-t62 * t30 + t81 * t31) + t81 * t93 - t62 * t94 + (t81 * t53 + t62 * t55) * t91) - m(6) * (g(1) * (rSges(6,3) * t31 - t67 * t30 + t79) + g(2) * (rSges(6,3) * t29 - t67 * t28 + t80) + t92 + (t67 * t55 + (rSges(6,3) - t45) * t53) * t91) - m(7) * (g(1) * (t100 * t31 + t79) + g(2) * (t100 * t29 + t80) + t92 + ((-t45 + t100) * t53 + t98 * t55) * t91 - t99 * t98) t72 * (t99 - t101) -m(5) * (g(1) * (t103 * rSges(5,1) + (-t31 * t42 - t41 * t84) * rSges(5,2)) + g(2) * (t61 * rSges(5,1) + (-t29 * t42 + t41 * t70) * rSges(5,2)) + g(3) * (t104 * rSges(5,1) + (-t76 * t41 - t42 * t83) * rSges(5,2))) - m(6) * (g(1) * (t68 + t96) + g(2) * (t57 + t97) + g(3) * (t66 + t78)) - m(7) * (g(1) * (t59 + t68) + g(2) * (t57 + t60) + g(3) * (t58 + t66)) -m(6) * (g(1) * t96 + g(2) * t97 + g(3) * t78) - m(7) * (g(1) * t59 + g(2) * t60 + g(3) * t58) -m(7) * (g(1) * ((-t15 * t52 + t30 * t54) * rSges(7,1) + (-t15 * t54 - t30 * t52) * rSges(7,2)) + g(2) * ((-t13 * t52 + t28 * t54) * rSges(7,1) + (-t13 * t54 - t28 * t52) * rSges(7,2)) + g(3) * ((-t24 * t52 - t54 * t82) * rSges(7,1) + (-t24 * t54 + t52 * t82) * rSges(7,2)))];
taug  = t1(:);
