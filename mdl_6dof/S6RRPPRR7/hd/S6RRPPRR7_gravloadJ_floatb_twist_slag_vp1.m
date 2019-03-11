% Calculate Gravitation load on the joints for
% S6RRPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR7_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:17:16
% EndTime: 2019-03-09 09:17:18
% DurationCPUTime: 0.90s
% Computational Cost: add. (435->158), mult. (995->222), div. (0->0), fcn. (1149->10), ass. (0->61)
t43 = sin(qJ(2));
t44 = sin(qJ(1));
t47 = cos(qJ(2));
t48 = cos(qJ(1));
t71 = cos(pkin(6));
t63 = t48 * t71;
t26 = t43 * t63 + t44 * t47;
t41 = sin(qJ(6));
t45 = cos(qJ(6));
t25 = t43 * t44 - t47 * t63;
t42 = sin(qJ(5));
t46 = cos(qJ(5));
t40 = sin(pkin(6));
t80 = t40 * t48;
t5 = t25 * t46 + t42 * t80;
t96 = -t26 * t45 + t41 * t5;
t95 = -t26 * t41 - t45 * t5;
t87 = rSges(6,3) + pkin(9);
t64 = t44 * t71;
t28 = -t43 * t64 + t47 * t48;
t94 = g(1) * t28 + g(2) * t26;
t27 = t48 * t43 + t47 * t64;
t93 = g(1) * t27 + g(2) * t25;
t88 = g(3) * t40;
t86 = rSges(7,3) + pkin(10);
t83 = t40 * t43;
t82 = t40 * t44;
t81 = t40 * t47;
t13 = t25 * pkin(2);
t79 = -t25 * pkin(3) - t13;
t19 = t27 * pkin(2);
t78 = -t27 * pkin(3) - t19;
t77 = pkin(2) * t81 + qJ(3) * t83;
t76 = t48 * pkin(1) + pkin(8) * t82;
t75 = qJ(4) * t40;
t74 = rSges(5,1) + qJ(3);
t73 = rSges(4,3) + qJ(3);
t70 = -m(5) - m(6) - m(7);
t69 = t26 * pkin(4) + t79;
t68 = t28 * pkin(4) + t78;
t67 = t28 * pkin(2) + t76;
t66 = pkin(3) * t81 + t77;
t65 = -t44 * pkin(1) + pkin(8) * t80;
t62 = -t26 * pkin(2) + t65;
t61 = rSges(6,1) * t46 - rSges(6,2) * t42;
t60 = rSges(7,1) * t41 + rSges(7,2) * t45;
t59 = rSges(7,1) * t45 - rSges(7,2) * t41 + pkin(5);
t57 = -t25 * t42 + t46 * t80;
t56 = g(3) * (pkin(4) * t83 + pkin(9) * t81 + t66);
t54 = t28 * pkin(3) + qJ(3) * t27 + t67;
t53 = -t26 * pkin(3) - t25 * qJ(3) + t62;
t52 = t27 * pkin(4) - t44 * t75 + t54;
t51 = t86 * t42 + t59 * t46;
t50 = -t25 * pkin(4) - t48 * t75 + t53;
t24 = t71 * t42 + t46 * t81;
t23 = t42 * t81 - t71 * t46;
t9 = t27 * t46 - t42 * t82;
t8 = t27 * t42 + t46 * t82;
t2 = t28 * t41 + t45 * t9;
t1 = t28 * t45 - t41 * t9;
t3 = [-m(2) * (g(1) * (-t44 * rSges(2,1) - rSges(2,2) * t48) + g(2) * (rSges(2,1) * t48 - t44 * rSges(2,2))) - m(3) * (g(1) * (-t26 * rSges(3,1) + t25 * rSges(3,2) + rSges(3,3) * t80 + t65) + g(2) * (rSges(3,1) * t28 - rSges(3,2) * t27 + rSges(3,3) * t82 + t76)) - m(4) * (g(1) * (-t26 * rSges(4,1) + rSges(4,2) * t80 - t73 * t25 + t62) + g(2) * (rSges(4,1) * t28 + rSges(4,2) * t82 + t73 * t27 + t67)) - m(5) * (g(1) * (-t25 * rSges(5,1) + t26 * rSges(5,2) + t53) + g(2) * (rSges(5,1) * t27 - rSges(5,2) * t28 + t54) + (g(1) * t48 + g(2) * t44) * t40 * (-rSges(5,3) - qJ(4))) - m(6) * (g(1) * (-rSges(6,1) * t5 - rSges(6,2) * t57 - t87 * t26 + t50) + g(2) * (rSges(6,1) * t9 - rSges(6,2) * t8 + t87 * t28 + t52)) - m(7) * (g(1) * (t95 * rSges(7,1) + t96 * rSges(7,2) - t5 * pkin(5) - t26 * pkin(9) + t86 * t57 + t50) + g(2) * (rSges(7,1) * t2 + rSges(7,2) * t1 + pkin(5) * t9 + pkin(9) * t28 + t86 * t8 + t52)) -m(3) * (g(1) * (-rSges(3,1) * t27 - rSges(3,2) * t28) + g(2) * (-rSges(3,1) * t25 - rSges(3,2) * t26) + (rSges(3,1) * t47 - rSges(3,2) * t43) * t88) - m(4) * (g(1) * (-rSges(4,1) * t27 + t73 * t28 - t19) + g(2) * (-rSges(4,1) * t25 + t73 * t26 - t13) + g(3) * ((rSges(4,1) * t47 + rSges(4,3) * t43) * t40 + t77)) - m(5) * (g(1) * (rSges(5,2) * t27 + t74 * t28 + t78) + g(2) * (rSges(5,2) * t25 + t74 * t26 + t79) + g(3) * ((rSges(5,1) * t43 - rSges(5,2) * t47) * t40 + t66)) - m(6) * (g(1) * (-t87 * t27 + t68) + g(2) * (-t87 * t25 + t69) + t56 + (rSges(6,3) * t47 + t61 * t43) * t88 + t94 * (qJ(3) + t61)) - m(7) * (g(1) * t68 + g(2) * t69 + t56 + (t51 * t43 + t60 * t47) * t88 + t93 * (-pkin(9) - t60) + t94 * (qJ(3) + t51)) (-m(4) + t70) * (-g(3) * t81 + t93) t70 * (-g(3) * t71 + (-g(1) * t44 + g(2) * t48) * t40) -m(6) * (g(1) * (-rSges(6,1) * t8 - rSges(6,2) * t9) + g(2) * (rSges(6,1) * t57 - rSges(6,2) * t5) + g(3) * (rSges(6,1) * t23 + rSges(6,2) * t24)) - m(7) * (g(1) * (-t59 * t8 + t86 * t9) + (t59 * t23 - t86 * t24) * g(3) + (t86 * t5 + t57 * t59) * g(2)) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (-t96 * rSges(7,1) + t95 * rSges(7,2)) + g(3) * ((t24 * t41 + t45 * t83) * rSges(7,1) + (t24 * t45 - t41 * t83) * rSges(7,2)))];
taug  = t3(:);
