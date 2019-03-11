% Calculate Gravitation load on the joints for
% S6RRPPRR5
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
% Datum: 2019-03-09 09:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:07:29
% EndTime: 2019-03-09 09:07:30
% DurationCPUTime: 0.82s
% Computational Cost: add. (435->156), mult. (995->221), div. (0->0), fcn. (1149->10), ass. (0->63)
t42 = sin(qJ(5));
t46 = cos(qJ(5));
t93 = -rSges(6,1) * t46 + rSges(6,2) * t42;
t43 = sin(qJ(2));
t44 = sin(qJ(1));
t47 = cos(qJ(2));
t71 = cos(pkin(6));
t85 = cos(qJ(1));
t57 = t71 * t85;
t25 = t43 * t44 - t47 * t57;
t61 = t44 * t71;
t27 = t43 * t85 + t47 * t61;
t92 = g(1) * t27 + g(2) * t25;
t73 = rSges(5,2) + qJ(3);
t41 = sin(qJ(6));
t45 = cos(qJ(6));
t54 = t45 * rSges(7,1) - t41 * rSges(7,2) + pkin(5);
t86 = rSges(7,3) + pkin(10);
t91 = t42 * t86 + t46 * t54;
t28 = -t43 * t61 + t47 * t85;
t89 = g(1) * t28;
t40 = sin(pkin(6));
t87 = g(3) * t40;
t82 = t40 * t43;
t81 = t40 * t44;
t80 = t40 * t46;
t79 = t40 * t47;
t78 = pkin(9) - qJ(3);
t13 = t25 * pkin(2);
t77 = -t25 * pkin(3) - t13;
t19 = t27 * pkin(2);
t76 = -t27 * pkin(3) - t19;
t75 = pkin(2) * t79 + qJ(3) * t82;
t74 = t85 * pkin(1) + pkin(8) * t81;
t72 = rSges(4,3) + qJ(3);
t70 = -m(5) - m(6) - m(7);
t69 = rSges(6,3) + t78;
t68 = -t25 * pkin(4) + t77;
t67 = -t27 * pkin(4) + t76;
t66 = t28 * pkin(2) + t74;
t65 = pkin(3) * t79 + t75;
t64 = t40 * t85;
t63 = -t44 * pkin(1) + pkin(8) * t64;
t62 = t85 * qJ(4);
t60 = t28 * pkin(3) + t66;
t59 = g(2) * t69;
t26 = t43 * t57 + t44 * t47;
t58 = -t26 * pkin(2) + t63;
t56 = g(3) * (pkin(4) * t79 + t65);
t55 = -t26 * pkin(3) + t58;
t53 = -t41 * rSges(7,1) - t45 * rSges(7,2) - pkin(9);
t52 = -qJ(3) - t53;
t51 = -t26 * t42 + t46 * t64;
t5 = t26 * t46 + t42 * t64;
t50 = t28 * pkin(4) - qJ(4) * t81 + t60;
t49 = -t26 * pkin(4) - t40 * t62 + t55;
t24 = -t42 * t71 + t43 * t80;
t23 = -t42 * t82 - t46 * t71;
t9 = t28 * t46 - t42 * t81;
t8 = t28 * t42 + t44 * t80;
t2 = -t27 * t41 + t45 * t9;
t1 = -t27 * t45 - t41 * t9;
t3 = [-m(2) * (g(1) * (-t44 * rSges(2,1) - rSges(2,2) * t85) + g(2) * (rSges(2,1) * t85 - t44 * rSges(2,2))) - m(3) * (g(1) * (-t26 * rSges(3,1) + t25 * rSges(3,2) + rSges(3,3) * t64 + t63) + g(2) * (rSges(3,1) * t28 - rSges(3,2) * t27 + rSges(3,3) * t81 + t74)) - m(4) * (g(1) * (-t26 * rSges(4,1) + rSges(4,2) * t64 - t25 * t72 + t58) + g(2) * (rSges(4,1) * t28 + rSges(4,2) * t81 + t27 * t72 + t66)) - m(5) * (g(1) * (-t26 * rSges(5,1) - t73 * t25 + t55) + g(2) * (rSges(5,1) * t28 + t73 * t27 + t60) + (g(1) * (-rSges(5,3) * t85 - t62) + g(2) * (-rSges(5,3) - qJ(4)) * t44) * t40) - m(6) * (g(2) * (rSges(6,1) * t9 - rSges(6,2) * t8 + t50) - t27 * t59 + (-rSges(6,1) * t5 - rSges(6,2) * t51 + t69 * t25 + t49) * g(1)) - m(7) * (g(1) * (t25 * t52 - t5 * t54 + t51 * t86 + t49) + g(2) * (rSges(7,1) * t2 + rSges(7,2) * t1 + pkin(5) * t9 - t27 * t78 + t8 * t86 + t50)) -m(3) * (g(1) * (-rSges(3,1) * t27 - rSges(3,2) * t28) + g(2) * (-rSges(3,1) * t25 - rSges(3,2) * t26) + (rSges(3,1) * t47 - rSges(3,2) * t43) * t87) - m(4) * (g(1) * (-rSges(4,1) * t27 + t28 * t72 - t19) + g(2) * (-rSges(4,1) * t25 + t26 * t72 - t13) + g(3) * ((rSges(4,1) * t47 + rSges(4,3) * t43) * t40 + t75)) - m(5) * (g(1) * (-rSges(5,1) * t27 + t28 * t73 + t76) + g(2) * (-rSges(5,1) * t25 + t26 * t73 + t77) + g(3) * ((rSges(5,1) * t47 + rSges(5,2) * t43) * t40 + t65)) - m(6) * (g(1) * (t93 * t27 + t67) + g(2) * (t93 * t25 + t68) + t56 - t69 * t89 - t26 * t59 + (-t93 * t47 + (-rSges(6,3) - pkin(9)) * t43) * t87) - m(7) * (g(1) * t67 + t56 - t52 * t89 + (t53 * t43 + t91 * t47) * t87 - t92 * t91 + (-t52 * t26 + t68) * g(2)) (-m(4) + t70) * (-g(3) * t79 + t92) t70 * (-g(3) * t71 + (-g(1) * t44 + g(2) * t85) * t40) -m(6) * (g(1) * (-rSges(6,1) * t8 - rSges(6,2) * t9) + g(2) * (rSges(6,1) * t51 - rSges(6,2) * t5) + g(3) * (rSges(6,1) * t23 - rSges(6,2) * t24)) - m(7) * (g(1) * (-t54 * t8 + t86 * t9) + (t23 * t54 + t24 * t86) * g(3) + (t5 * t86 + t51 * t54) * g(2)) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * ((-t25 * t45 - t41 * t5) * rSges(7,1) + (t25 * t41 - t45 * t5) * rSges(7,2)) + g(3) * ((-t24 * t41 + t45 * t79) * rSges(7,1) + (-t24 * t45 - t41 * t79) * rSges(7,2)))];
taug  = t3(:);
