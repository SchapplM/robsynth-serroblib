% Calculate Gravitation load on the joints for
% S6RRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2019-03-09 15:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:21:14
% EndTime: 2019-03-09 15:21:17
% DurationCPUTime: 0.95s
% Computational Cost: add. (539->134), mult. (453->167), div. (0->0), fcn. (402->12), ass. (0->75)
t44 = qJ(2) + qJ(3);
t38 = pkin(10) + t44;
t31 = sin(t38);
t32 = cos(t38);
t46 = cos(pkin(11));
t73 = -rSges(6,1) * t46 - pkin(4);
t45 = sin(pkin(11));
t89 = rSges(6,2) * t45;
t111 = (rSges(6,3) + qJ(5)) * t31 + (-t73 - t89) * t32;
t108 = qJ(5) * t32 + t31 * t89;
t33 = pkin(5) * t46 + pkin(4);
t106 = t31 * rSges(7,3) + t32 * t33;
t65 = t32 * rSges(5,1) - rSges(5,2) * t31;
t39 = sin(t44);
t40 = cos(t44);
t70 = t40 * rSges(4,1) - rSges(4,2) * t39;
t49 = sin(qJ(1));
t51 = cos(qJ(1));
t104 = g(1) * t51 + g(2) * t49;
t103 = t104 * t31;
t34 = pkin(3) * t40;
t50 = cos(qJ(2));
t41 = t50 * pkin(2);
t35 = t41 + pkin(1);
t15 = t34 + t35;
t8 = t51 * t15;
t102 = g(2) * t8;
t101 = -m(6) - m(7);
t52 = -pkin(8) - pkin(7);
t48 = sin(qJ(2));
t100 = pkin(2) * t48;
t99 = pkin(3) * t39;
t96 = rSges(3,3) + pkin(7);
t43 = pkin(11) + qJ(6);
t36 = sin(t43);
t88 = rSges(7,2) * t36;
t76 = t31 * t88;
t86 = t32 * t49;
t95 = rSges(7,3) * t86 + t49 * t76;
t85 = t32 * t51;
t94 = rSges(7,3) * t85 + t51 * t76;
t37 = cos(t43);
t92 = rSges(7,1) * t37;
t47 = -pkin(9) - qJ(5);
t87 = t31 * t47;
t84 = t36 * t49;
t83 = t36 * t51;
t82 = t37 * t49;
t81 = t37 * t51;
t80 = rSges(4,3) - t52;
t42 = -qJ(4) + t52;
t79 = rSges(5,3) - t42;
t75 = rSges(6,3) * t86 + t108 * t49;
t74 = rSges(6,3) * t85 + t108 * t51;
t72 = pkin(5) * t45 - t42;
t71 = -t33 - t92;
t69 = t34 + t65;
t68 = rSges(3,1) * t50 - rSges(3,2) * t48;
t66 = -rSges(4,1) * t39 - rSges(4,2) * t40;
t64 = -rSges(5,1) * t31 - rSges(5,2) * t32;
t63 = pkin(1) + t68;
t61 = t35 + t70;
t60 = rSges(6,1) * t45 + rSges(6,2) * t46 - t42;
t59 = t106 + t34 + (-t88 + t92) * t32;
t56 = -t87 + t106;
t55 = t34 + t111;
t53 = t73 * t103;
t16 = -t99 - t100;
t10 = t51 * t16;
t9 = t49 * t16;
t5 = t32 * t81 + t84;
t4 = -t32 * t83 + t82;
t3 = -t32 * t82 + t83;
t2 = t32 * t84 + t81;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t49 - rSges(2,2) * t51) + g(2) * (rSges(2,1) * t51 - rSges(2,2) * t49)) - m(3) * ((g(1) * t96 + g(2) * t63) * t51 + (-g(1) * t63 + g(2) * t96) * t49) - m(4) * ((g(1) * t80 + g(2) * t61) * t51 + (-g(1) * t61 + g(2) * t80) * t49) - m(5) * (t102 + (g(1) * t79 + g(2) * t65) * t51 + (g(1) * (-t15 - t65) + g(2) * t79) * t49) - m(6) * (t102 + (g(1) * t60 + g(2) * t111) * t51 + (g(2) * t60 + (-t15 - t111) * g(1)) * t49) - m(7) * (g(1) * (rSges(7,1) * t3 + rSges(7,2) * t2) + g(2) * (rSges(7,1) * t5 + rSges(7,2) * t4 + t8) + (g(1) * t72 + g(2) * t56) * t51 + (g(1) * (-t15 - t56) + g(2) * t72) * t49) -m(3) * (g(3) * t68 + t104 * (-rSges(3,1) * t48 - rSges(3,2) * t50)) - m(4) * (g(3) * (t41 + t70) + t104 * (t66 - t100)) - m(5) * (g(1) * (t64 * t51 + t10) + g(2) * (t64 * t49 + t9) + g(3) * (t41 + t69)) - m(6) * (g(1) * (t10 + t74) + g(2) * (t9 + t75) + g(3) * (t41 + t55) + t53) - m(7) * (g(1) * (-t47 * t85 + t10 + t94) + g(2) * (-t47 * t86 + t9 + t95) + g(3) * (t41 + t59) + (-g(3) * t47 + t104 * t71) * t31) -m(4) * (g(3) * t70 + t104 * t66) - m(5) * (g(3) * t69 + t104 * (t64 - t99)) - m(6) * (g(1) * (-t51 * t99 + t74) + g(2) * (-t49 * t99 + t75) + g(3) * t55 + t53) - m(7) * (g(1) * t94 + g(2) * t95 + g(3) * (t59 - t87) + t104 * (t71 * t31 - t32 * t47 - t99)) (-m(5) + t101) * (g(1) * t49 - g(2) * t51) t101 * (-g(3) * t32 + t103) -m(7) * (g(1) * (rSges(7,1) * t4 - rSges(7,2) * t5) + g(2) * (-rSges(7,1) * t2 + rSges(7,2) * t3) + g(3) * (-rSges(7,1) * t36 - rSges(7,2) * t37) * t31)];
taug  = t1(:);
