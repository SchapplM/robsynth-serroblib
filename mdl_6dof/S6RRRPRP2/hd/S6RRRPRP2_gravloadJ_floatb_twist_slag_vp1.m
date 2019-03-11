% Calculate Gravitation load on the joints for
% S6RRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:35:17
% EndTime: 2019-03-09 16:35:19
% DurationCPUTime: 0.78s
% Computational Cost: add. (580->141), mult. (532->176), div. (0->0), fcn. (492->10), ass. (0->77)
t104 = rSges(7,1) + pkin(5);
t74 = rSges(7,3) + qJ(6);
t46 = sin(qJ(1));
t43 = qJ(2) + qJ(3);
t38 = pkin(10) + t43;
t32 = sin(t38);
t44 = sin(qJ(5));
t87 = t32 * t44;
t72 = rSges(6,2) * t87;
t33 = cos(t38);
t84 = t33 * t46;
t103 = rSges(6,3) * t84 + t46 * t72;
t49 = cos(qJ(1));
t82 = t33 * t49;
t102 = rSges(6,3) * t82 + t49 * t72;
t62 = t33 * rSges(5,1) - rSges(5,2) * t32;
t101 = -t33 * pkin(4) - t32 * pkin(9);
t39 = sin(t43);
t40 = cos(t43);
t69 = t40 * rSges(4,1) - rSges(4,2) * t39;
t100 = g(1) * t49 + g(2) * t46;
t99 = t32 * t100;
t50 = -pkin(8) - pkin(7);
t45 = sin(qJ(2));
t98 = pkin(2) * t45;
t97 = pkin(3) * t39;
t42 = -qJ(4) + t50;
t95 = g(2) * t42;
t92 = rSges(3,3) + pkin(7);
t20 = pkin(9) * t84;
t11 = -t97 - t98;
t6 = t46 * t11;
t91 = t20 + t6;
t23 = pkin(9) * t82;
t7 = t49 * t11;
t90 = t23 + t7;
t48 = cos(qJ(2));
t41 = t48 * pkin(2);
t37 = t41 + pkin(1);
t27 = t32 * rSges(7,2);
t26 = t32 * rSges(6,3);
t86 = t32 * t49;
t85 = t33 * t44;
t47 = cos(qJ(5));
t83 = t33 * t47;
t81 = t42 * t49;
t80 = t44 * t46;
t79 = t46 * t47;
t78 = t47 * t49;
t77 = t49 * t44;
t76 = rSges(4,3) - t50;
t75 = rSges(5,3) - t42;
t34 = pkin(3) * t40;
t10 = t34 + t37;
t5 = t49 * t10;
t73 = pkin(4) * t82 + pkin(9) * t86 + t5;
t71 = t34 - t101;
t68 = -t46 * t97 + t20;
t67 = -t49 * t97 + t23;
t66 = t34 + t62;
t65 = rSges(3,1) * t48 - rSges(3,2) * t45;
t63 = -rSges(4,1) * t39 - rSges(4,2) * t40;
t61 = -rSges(5,1) * t32 - rSges(5,2) * t33;
t60 = -t10 + t101;
t59 = pkin(1) + t65;
t58 = t37 + t69;
t55 = t104 * t83 + t74 * t85 + t27 + t71;
t54 = rSges(6,1) * t83 - rSges(6,2) * t85 + t26 + t71;
t52 = (-rSges(6,1) * t47 - pkin(4)) * t99;
t51 = (-t104 * t47 - t44 * t74 - pkin(4)) * t99;
t19 = rSges(7,2) * t82;
t14 = rSges(7,2) * t84;
t4 = t33 * t78 + t80;
t3 = t33 * t77 - t79;
t2 = t33 * t79 - t77;
t1 = t33 * t80 + t78;
t8 = [-m(2) * (g(1) * (-rSges(2,1) * t46 - rSges(2,2) * t49) + g(2) * (rSges(2,1) * t49 - rSges(2,2) * t46)) - m(3) * ((g(1) * t92 + g(2) * t59) * t49 + (-g(1) * t59 + g(2) * t92) * t46) - m(4) * ((g(1) * t76 + g(2) * t58) * t49 + (-g(1) * t58 + g(2) * t76) * t46) - m(5) * (g(2) * t5 + (g(1) * t75 + g(2) * t62) * t49 + (g(1) * (-t10 - t62) + g(2) * t75) * t46) - m(6) * (g(1) * (-rSges(6,1) * t2 + rSges(6,2) * t1 - t81) + g(2) * (rSges(6,1) * t4 - rSges(6,2) * t3 + rSges(6,3) * t86 + t73) + (g(1) * (t60 - t26) - t95) * t46) - m(7) * (g(1) * (-t74 * t1 - t104 * t2 - t81) + g(2) * (rSges(7,2) * t86 + t104 * t4 + t74 * t3 + t73) + (g(1) * (t60 - t27) - t95) * t46) -m(3) * (g(3) * t65 + t100 * (-rSges(3,1) * t45 - rSges(3,2) * t48)) - m(4) * (g(3) * (t41 + t69) + t100 * (t63 - t98)) - m(5) * (g(1) * (t49 * t61 + t7) + g(2) * (t46 * t61 + t6) + g(3) * (t41 + t66)) - m(6) * (g(1) * (t90 + t102) + g(2) * (t91 + t103) + g(3) * (t41 + t54) + t52) - m(7) * (g(1) * (t19 + t90) + g(2) * (t14 + t91) + g(3) * (t41 + t55) + t51) -m(4) * (g(3) * t69 + t100 * t63) - m(5) * (g(3) * t66 + t100 * (t61 - t97)) - m(6) * (g(1) * (t67 + t102) + g(2) * (t68 + t103) + g(3) * t54 + t52) - m(7) * (g(1) * (t19 + t67) + g(2) * (t14 + t68) + g(3) * t55 + t51) (-m(5) - m(6) - m(7)) * (g(1) * t46 - g(2) * t49) -m(6) * (g(1) * (-rSges(6,1) * t3 - rSges(6,2) * t4) + g(2) * (-rSges(6,1) * t1 - rSges(6,2) * t2)) - m(7) * (g(1) * (-t104 * t3 + t4 * t74) + g(2) * (-t1 * t104 + t2 * t74)) + (-m(6) * (-rSges(6,1) * t44 - rSges(6,2) * t47) - m(7) * (-t104 * t44 + t74 * t47)) * g(3) * t32, -m(7) * (g(1) * t3 + g(2) * t1 + g(3) * t87)];
taug  = t8(:);
