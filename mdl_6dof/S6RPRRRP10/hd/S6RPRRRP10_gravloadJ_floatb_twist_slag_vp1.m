% Calculate Gravitation load on the joints for
% S6RPRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP10_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP10_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP10_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP10_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:30:18
% EndTime: 2019-03-09 06:30:21
% DurationCPUTime: 0.82s
% Computational Cost: add. (370->136), mult. (527->182), div. (0->0), fcn. (513->8), ass. (0->65)
t42 = sin(qJ(3));
t45 = cos(qJ(3));
t81 = rSges(5,3) + pkin(8);
t95 = t81 * t45;
t98 = t42 * pkin(3) - t95;
t82 = rSges(7,1) + pkin(5);
t44 = cos(qJ(4));
t33 = pkin(4) * t44 + pkin(3);
t40 = qJ(4) + qJ(5);
t34 = sin(t40);
t35 = cos(t40);
t64 = rSges(7,3) + qJ(6);
t93 = t64 * t34 + t82 * t35;
t97 = -t33 - t93;
t41 = sin(qJ(4));
t43 = sin(qJ(1));
t75 = t43 * t44;
t46 = cos(qJ(1));
t76 = t42 * t46;
t17 = t41 * t76 + t75;
t73 = t44 * t46;
t77 = t42 * t43;
t15 = -t41 * t77 + t73;
t85 = g(2) * t46;
t87 = g(1) * t43;
t94 = -t85 + t87;
t91 = -pkin(1) - pkin(7);
t71 = t46 * t35;
t11 = t34 * t77 - t71;
t12 = t34 * t46 + t35 * t77;
t90 = -t11 * rSges(6,1) - t12 * rSges(6,2);
t13 = t34 * t76 + t35 * t43;
t14 = -t34 * t43 + t42 * t71;
t89 = t13 * rSges(6,1) + t14 * rSges(6,2);
t88 = pkin(4) * t41;
t86 = g(2) * t45;
t84 = g(3) * t45;
t80 = t34 * t45;
t78 = t41 * t46;
t74 = t43 * t45;
t72 = t45 * t46;
t47 = -pkin(9) - pkin(8);
t70 = t46 * t47;
t69 = rSges(7,2) - t47;
t68 = rSges(6,3) - t47;
t67 = t17 * pkin(4);
t66 = t64 * t35 * t45;
t65 = t46 * pkin(1) + t43 * qJ(2);
t37 = t46 * qJ(2);
t61 = t33 * t76 + t45 * t70 + t37;
t60 = t46 * pkin(7) + t65;
t59 = g(1) * t33 * t74 + g(2) * t42 * t70;
t57 = rSges(4,1) * t42 + rSges(4,2) * t45;
t56 = rSges(6,1) * t35 - rSges(6,2) * t34;
t55 = -rSges(6,1) * t34 - rSges(6,2) * t35;
t54 = -t82 * t11 + t64 * t12;
t53 = t15 * pkin(4);
t52 = t82 * t13 - t64 * t14;
t51 = g(1) * (-t88 + t91);
t50 = pkin(4) * t78 + t33 * t77 + t47 * t74 + t60;
t49 = rSges(5,1) * t44 - rSges(5,2) * t41 + pkin(3);
t48 = -t33 - t56;
t18 = -t41 * t43 + t42 * t73;
t16 = t42 * t75 + t78;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t43 - rSges(2,2) * t46) + g(2) * (rSges(2,1) * t46 - rSges(2,2) * t43)) - m(3) * (g(1) * (rSges(3,3) * t46 + t37 + (rSges(3,2) - pkin(1)) * t43) + g(2) * (-rSges(3,2) * t46 + rSges(3,3) * t43 + t65)) - m(4) * (g(1) * (rSges(4,1) * t76 + rSges(4,2) * t72 + t37) + g(2) * (rSges(4,3) * t46 + t60) + (g(1) * (-rSges(4,3) + t91) + g(2) * t57) * t43) - m(5) * ((rSges(5,1) * t16 + rSges(5,2) * t15 + t98 * t43 + t60) * g(2) + (rSges(5,1) * t18 - rSges(5,2) * t17 + t91 * t43 + t98 * t46 + t37) * g(1)) - m(6) * (g(1) * (rSges(6,1) * t14 - rSges(6,2) * t13 - rSges(6,3) * t72 + t61) + g(2) * (rSges(6,1) * t12 - rSges(6,2) * t11 + t50) + (-rSges(6,3) * t86 + t51) * t43) - m(7) * (g(1) * (-rSges(7,2) * t72 + t64 * t13 + t82 * t14 + t61) + g(2) * (t64 * t11 + t82 * t12 + t50) + (-rSges(7,2) * t86 + t51) * t43) (-m(3) - m(4) - m(5) - m(6) - m(7)) * t94, -m(4) * (-g(3) * t57 + t94 * (rSges(4,1) * t45 - rSges(4,2) * t42)) - m(5) * (g(3) * (-t49 * t42 + t95) + t94 * (t81 * t42 + t49 * t45)) - m(6) * ((-rSges(6,3) * t85 + g(3) * t48 + t68 * t87) * t42 + (g(3) * t68 + t48 * t85 + t56 * t87) * t45 + t59) - m(7) * ((-rSges(7,2) * t85 + g(3) * t97 + t69 * t87) * t42 + (g(3) * t69 + t97 * t85 + t93 * t87) * t45 + t59) -m(5) * (g(1) * (rSges(5,1) * t15 - rSges(5,2) * t16) + g(2) * (rSges(5,1) * t17 + rSges(5,2) * t18)) - m(6) * (g(1) * (t53 + t90) + g(2) * (t67 + t89)) - m(7) * (g(1) * (t53 + t54) + g(2) * (t52 + t67) + g(3) * t66) + (-m(5) * (-rSges(5,1) * t41 - rSges(5,2) * t44) - m(6) * (t55 - t88) - m(7) * (-t82 * t34 - t88)) * t84, -m(6) * (g(1) * t90 + g(2) * t89 + t55 * t84) - m(7) * (g(1) * t54 + g(2) * t52 + g(3) * (-t82 * t80 + t66)) -m(7) * (g(1) * t11 - g(2) * t13 + g(3) * t80)];
taug  = t1(:);
