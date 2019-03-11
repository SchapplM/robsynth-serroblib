% Calculate Gravitation load on the joints for
% S6RRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:24:58
% EndTime: 2019-03-09 15:25:00
% DurationCPUTime: 0.72s
% Computational Cost: add. (491->130), mult. (425->165), div. (0->0), fcn. (370->10), ass. (0->76)
t43 = sin(qJ(6));
t46 = cos(qJ(6));
t105 = rSges(7,1) * t43 + rSges(7,2) * t46;
t42 = qJ(2) + qJ(3);
t37 = pkin(10) + t42;
t34 = cos(t37);
t94 = rSges(7,3) + pkin(9);
t104 = t34 * t94;
t103 = t105 * t34;
t33 = sin(t37);
t59 = t34 * rSges(5,1) - t33 * rSges(5,2);
t38 = sin(t42);
t39 = cos(t42);
t67 = t39 * rSges(4,1) - t38 * rSges(4,2);
t57 = -t34 * rSges(6,2) + t33 * rSges(6,3);
t48 = cos(qJ(1));
t45 = sin(qJ(1));
t98 = g(2) * t45;
t102 = g(1) * t48 + t98;
t101 = -m(6) - m(7);
t49 = -pkin(8) - pkin(7);
t100 = pkin(3) * t38;
t97 = g(3) * t34;
t31 = t34 * pkin(4);
t44 = sin(qJ(2));
t96 = t44 * pkin(2);
t95 = rSges(3,3) + pkin(7);
t41 = -qJ(4) + t49;
t93 = pkin(5) - t41;
t47 = cos(qJ(2));
t40 = t47 * pkin(2);
t36 = t40 + pkin(1);
t89 = t33 * t45;
t88 = t33 * t48;
t86 = t34 * t48;
t84 = t45 * t43;
t83 = t45 * t46;
t82 = t48 * t43;
t81 = t48 * t46;
t80 = rSges(6,1) - t41;
t79 = rSges(4,3) - t49;
t78 = rSges(5,3) - t41;
t77 = qJ(5) * t34;
t26 = t33 * qJ(5);
t76 = -pkin(4) - t94;
t15 = t45 * t77;
t75 = t103 * t45 + t15;
t35 = pkin(3) * t39;
t13 = t35 + t36;
t6 = t48 * t13;
t74 = pkin(4) * t86 + t48 * t26 + t6;
t17 = t48 * t77;
t71 = t103 * t48 + t17;
t70 = t45 * t34 * rSges(6,3) + rSges(6,2) * t89 + t15;
t69 = rSges(6,2) * t88 + rSges(6,3) * t86 + t17;
t68 = t26 + t31 + t35;
t66 = -t13 - t26;
t65 = g(1) * t76;
t64 = t35 + t59;
t63 = -pkin(4) * t33 - t100;
t62 = t47 * rSges(3,1) - t44 * rSges(3,2);
t60 = -rSges(4,1) * t38 - rSges(4,2) * t39;
t58 = -rSges(5,1) * t33 - rSges(5,2) * t34;
t56 = pkin(1) + t62;
t55 = t36 + t67;
t54 = t68 + t57;
t53 = t105 * t33 + t104 + t68;
t50 = (t48 * t65 + t76 * t98) * t33;
t14 = -t96 - t100;
t8 = t48 * t14;
t7 = t45 * t14;
t5 = -t33 * t84 + t81;
t4 = t33 * t83 + t82;
t3 = t33 * t82 + t83;
t2 = t33 * t81 - t84;
t1 = [-m(2) * (g(1) * (-t45 * rSges(2,1) - t48 * rSges(2,2)) + g(2) * (t48 * rSges(2,1) - t45 * rSges(2,2))) - m(3) * ((g(1) * t95 + g(2) * t56) * t48 + (-g(1) * t56 + g(2) * t95) * t45) - m(4) * ((g(1) * t79 + g(2) * t55) * t48 + (-g(1) * t55 + g(2) * t79) * t45) - m(5) * (g(2) * t6 + (g(1) * t78 + g(2) * t59) * t48 + (g(1) * (-t13 - t59) + g(2) * t78) * t45) - m(6) * (g(2) * t74 + (g(1) * t80 + g(2) * t57) * t48 + (g(1) * (-t57 + t66 - t31) + g(2) * t80) * t45) - m(7) * (g(1) * (t5 * rSges(7,1) - t4 * rSges(7,2)) + g(2) * (t3 * rSges(7,1) + t2 * rSges(7,2) + t74) + (g(1) * t93 + g(2) * t104) * t48 + (g(1) * t66 + g(2) * t93 + t34 * t65) * t45) -m(3) * (g(3) * t62 + t102 * (-rSges(3,1) * t44 - rSges(3,2) * t47)) - m(4) * (g(3) * (t40 + t67) + t102 * (t60 - t96)) - m(5) * (g(1) * (t58 * t48 + t8) + g(2) * (t58 * t45 + t7) + g(3) * (t40 + t64)) - m(6) * (g(1) * (-pkin(4) * t88 + t69 + t8) + g(2) * (-pkin(4) * t89 + t7 + t70) + g(3) * (t40 + t54)) - m(7) * (g(1) * (t8 + t71) + g(2) * (t7 + t75) + g(3) * (t40 + t53) + t50) -m(4) * (g(3) * t67 + t102 * t60) - m(5) * (g(3) * t64 + t102 * (t58 - t100)) - m(6) * (g(1) * (t63 * t48 + t69) + g(2) * (t63 * t45 + t70) + g(3) * t54) - m(7) * (g(1) * (-t48 * t100 + t71) + g(2) * (-t45 * t100 + t75) + g(3) * t53 + t50) (-m(5) + t101) * (g(1) * t45 - g(2) * t48) t101 * (t102 * t33 - t97) -m(7) * (g(1) * (t2 * rSges(7,1) - t3 * rSges(7,2)) + g(2) * (t4 * rSges(7,1) + t5 * rSges(7,2)) + (-rSges(7,1) * t46 + rSges(7,2) * t43) * t97)];
taug  = t1(:);
