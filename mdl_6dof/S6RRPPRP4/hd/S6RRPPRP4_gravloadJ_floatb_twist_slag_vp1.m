% Calculate Gravitation load on the joints for
% S6RRPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% Datum: 2019-03-09 08:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRP4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRP4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:37:02
% EndTime: 2019-03-09 08:37:04
% DurationCPUTime: 0.91s
% Computational Cost: add. (366->157), mult. (849->214), div. (0->0), fcn. (928->8), ass. (0->62)
t89 = rSges(7,1) + pkin(5);
t48 = sin(qJ(2));
t45 = sin(pkin(9));
t46 = cos(pkin(9));
t47 = sin(qJ(5));
t50 = cos(qJ(5));
t61 = t45 * t47 + t46 * t50;
t95 = t61 * t48;
t52 = cos(qJ(1));
t49 = sin(qJ(1));
t90 = g(2) * t49;
t94 = g(1) * t52 + t90;
t73 = rSges(7,3) + qJ(6);
t93 = t48 * (-t45 * t50 + t46 * t47);
t92 = g(1) * t49;
t51 = cos(qJ(2));
t42 = t51 * pkin(2);
t88 = -rSges(7,2) - pkin(8);
t87 = -rSges(6,3) - pkin(8);
t85 = t45 * t51;
t83 = t46 * t51;
t82 = t48 * t52;
t81 = t49 * t51;
t80 = t51 * t52;
t79 = t52 * t45;
t40 = t48 * qJ(3);
t78 = t40 + t42;
t77 = t52 * pkin(1) + t49 * pkin(7);
t76 = qJ(3) * t51;
t75 = qJ(4) * t45;
t74 = rSges(5,3) + qJ(4);
t72 = -m(5) - m(6) - m(7);
t71 = -pkin(1) - t42;
t70 = t88 * t52;
t69 = t87 * t52;
t68 = pkin(3) * t83 + t51 * t75 + t78;
t67 = pkin(2) * t80 + t52 * t40 + t77;
t24 = t45 * t81 + t46 * t52;
t25 = t46 * t81 - t79;
t43 = t52 * pkin(7);
t66 = -t25 * pkin(3) - qJ(4) * t24 + t43;
t27 = t49 * t45 + t46 * t80;
t65 = t27 * pkin(3) + t67;
t64 = pkin(4) * t83 + t68;
t63 = rSges(3,1) * t51 - rSges(3,2) * t48;
t4 = t24 * t50 - t25 * t47;
t3 = t24 * t47 + t25 * t50;
t59 = t49 * t48 * pkin(8) - t25 * pkin(4) + t66;
t26 = -t49 * t46 + t51 * t79;
t55 = t27 * pkin(4) + t26 * qJ(4) + t65;
t53 = t94 * (-t75 - pkin(2) + (-pkin(3) - pkin(4)) * t46);
t35 = t52 * t76;
t31 = t49 * t76;
t19 = t61 * t51;
t18 = t47 * t83 - t50 * t85;
t11 = t52 * t95;
t10 = t52 * t93;
t9 = t49 * t95;
t8 = t49 * t93;
t7 = t26 * t47 + t27 * t50;
t6 = -t26 * t50 + t27 * t47;
t1 = [-m(2) * (g(1) * (-t49 * rSges(2,1) - rSges(2,2) * t52) + g(2) * (rSges(2,1) * t52 - t49 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t52 + t43) + g(2) * (rSges(3,1) * t80 - rSges(3,2) * t82 + t77) + (g(1) * (-pkin(1) - t63) + g(2) * rSges(3,3)) * t49) - m(4) * (g(1) * (-rSges(4,1) * t25 + rSges(4,2) * t24 + t43) + g(2) * (t27 * rSges(4,1) - t26 * rSges(4,2) + rSges(4,3) * t82 + t67) + ((-rSges(4,3) - qJ(3)) * t48 + t71) * t92) - m(5) * (g(1) * (-rSges(5,1) * t25 - rSges(5,3) * t24 + t66) + g(2) * (t27 * rSges(5,1) + rSges(5,2) * t82 + t74 * t26 + t65) + ((-rSges(5,2) - qJ(3)) * t48 + t71) * t92) - m(6) * (g(1) * (-rSges(6,1) * t3 - rSges(6,2) * t4 + t59) + g(2) * (t7 * rSges(6,1) - t6 * rSges(6,2) + t48 * t69 + t55) + ((rSges(6,3) - qJ(3)) * t48 + t71) * t92) - m(7) * (g(1) * (-t3 * t89 + t73 * t4 + t59) + g(2) * (t48 * t70 + t73 * t6 + t89 * t7 + t55) + ((rSges(7,2) - qJ(3)) * t48 + t71) * t92) -m(3) * (g(3) * t63 + t94 * (-rSges(3,1) * t48 - rSges(3,2) * t51)) - m(4) * (g(1) * (rSges(4,3) * t80 + t35) + g(2) * (rSges(4,3) * t81 + t31) + g(3) * (rSges(4,1) * t83 - rSges(4,2) * t85 + t78) + (g(3) * rSges(4,3) + t94 * (-rSges(4,1) * t46 + rSges(4,2) * t45 - pkin(2))) * t48) - m(5) * (g(1) * (rSges(5,2) * t80 + t35) + g(2) * (rSges(5,2) * t81 + t31) + g(3) * (rSges(5,1) * t83 + rSges(5,3) * t85 + t68) + (g(3) * rSges(5,2) + t94 * (-pkin(2) + (-rSges(5,1) - pkin(3)) * t46 - t74 * t45)) * t48) - m(6) * (g(1) * (-t11 * rSges(6,1) + t10 * rSges(6,2) + t35) + g(2) * (-rSges(6,1) * t9 + rSges(6,2) * t8 + t31) + g(3) * (rSges(6,1) * t19 - rSges(6,2) * t18 + t64) + (g(1) * t69 + t87 * t90) * t51 + (g(3) * t87 + t53) * t48) - m(7) * (g(1) * (-t73 * t10 - t89 * t11 + t35) + g(2) * (-t73 * t8 - t89 * t9 + t31) + g(3) * (t73 * t18 + t89 * t19 + t64) + (g(1) * t70 + t88 * t90) * t51 + (g(3) * t88 + t53) * t48) (-m(4) + t72) * (-g(3) * t51 + t94 * t48) t72 * (g(3) * t45 * t48 + g(1) * t26 + g(2) * t24) -m(6) * (g(1) * (-rSges(6,1) * t6 - rSges(6,2) * t7) + g(2) * (rSges(6,1) * t4 - rSges(6,2) * t3) + g(3) * (-rSges(6,1) * t93 - rSges(6,2) * t95)) - m(7) * (g(1) * (-t89 * t6 + t73 * t7) + g(2) * (t73 * t3 + t4 * t89) + g(3) * (t73 * t95 - t89 * t93)) -m(7) * (g(1) * t6 - g(2) * t4 + g(3) * t93)];
taug  = t1(:);
