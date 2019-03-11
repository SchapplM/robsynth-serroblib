% Calculate Gravitation load on the joints for
% S6RRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:12:05
% EndTime: 2019-03-09 10:12:07
% DurationCPUTime: 0.66s
% Computational Cost: add. (463->122), mult. (399->154), div. (0->0), fcn. (347->10), ass. (0->72)
t43 = sin(qJ(6));
t46 = cos(qJ(6));
t105 = rSges(7,1) * t43 + rSges(7,2) * t46;
t41 = qJ(2) + pkin(10);
t38 = qJ(4) + t41;
t34 = cos(t38);
t92 = rSges(7,3) + pkin(9);
t104 = t34 * t92;
t103 = t105 * t34;
t33 = sin(t38);
t102 = t33 * (rSges(6,2) - pkin(4));
t101 = t34 * rSges(5,1) - t33 * rSges(5,2);
t36 = sin(t41);
t37 = cos(t41);
t47 = cos(qJ(2));
t39 = t47 * pkin(2);
t100 = t37 * rSges(4,1) - t36 * rSges(4,2) + t39;
t59 = -t34 * rSges(6,2) + t33 * rSges(6,3);
t48 = cos(qJ(1));
t45 = sin(qJ(1));
t96 = g(2) * t45;
t99 = g(1) * t48 + t96;
t98 = -m(6) - m(7);
t95 = g(3) * t34;
t31 = t34 * pkin(4);
t44 = sin(qJ(2));
t94 = t44 * pkin(2);
t93 = rSges(3,3) + pkin(7);
t42 = -qJ(3) - pkin(7);
t40 = -pkin(8) + t42;
t91 = pkin(5) - t40;
t84 = t34 * t48;
t83 = t45 * t43;
t82 = t45 * t46;
t81 = t48 * t43;
t80 = t48 * t46;
t79 = rSges(6,1) - t40;
t78 = rSges(4,3) - t42;
t77 = rSges(5,3) - t40;
t26 = t33 * qJ(5);
t76 = t26 + t31;
t75 = pkin(3) * t37 + t39;
t74 = qJ(5) * t34;
t73 = -pkin(4) - t92;
t15 = t45 * t74;
t72 = t103 * t45 + t15;
t13 = pkin(1) + t75;
t6 = t48 * t13;
t71 = pkin(4) * t84 + t48 * t26 + t6;
t17 = t48 * t74;
t68 = t103 * t48 + t17;
t66 = -t13 - t26;
t65 = g(1) * t73;
t64 = t47 * rSges(3,1) - t44 * rSges(3,2);
t60 = -rSges(5,1) * t33 - rSges(5,2) * t34;
t58 = t76 + t59;
t57 = t105 * t33 + t104 + t76;
t56 = pkin(1) + t64;
t55 = pkin(1) + t100;
t54 = t15 + (rSges(6,3) * t34 + t102) * t45;
t53 = rSges(6,3) * t84 + t48 * t102 + t17;
t52 = t60 * t45;
t51 = t60 * t48;
t49 = (t48 * t65 + t73 * t96) * t33;
t14 = -pkin(3) * t36 - t94;
t8 = t48 * t14;
t7 = t45 * t14;
t5 = -t33 * t83 + t80;
t4 = t33 * t82 + t81;
t3 = t33 * t81 + t82;
t2 = t33 * t80 - t83;
t1 = [-m(2) * (g(1) * (-t45 * rSges(2,1) - t48 * rSges(2,2)) + g(2) * (t48 * rSges(2,1) - t45 * rSges(2,2))) - m(3) * ((g(1) * t93 + g(2) * t56) * t48 + (-g(1) * t56 + g(2) * t93) * t45) - m(4) * ((g(1) * t78 + g(2) * t55) * t48 + (-g(1) * t55 + g(2) * t78) * t45) - m(5) * (g(2) * t6 + (g(1) * t77 + g(2) * t101) * t48 + (g(1) * (-t13 - t101) + g(2) * t77) * t45) - m(6) * (g(2) * t71 + (g(1) * t79 + g(2) * t59) * t48 + (g(1) * (-t59 + t66 - t31) + g(2) * t79) * t45) - m(7) * (g(1) * (t5 * rSges(7,1) - t4 * rSges(7,2)) + g(2) * (t3 * rSges(7,1) + t2 * rSges(7,2) + t71) + (g(1) * t91 + g(2) * t104) * t48 + (g(1) * t66 + g(2) * t91 + t34 * t65) * t45) -m(3) * (g(3) * t64 + t99 * (-rSges(3,1) * t44 - rSges(3,2) * t47)) - m(4) * (g(3) * t100 + t99 * (-rSges(4,1) * t36 - rSges(4,2) * t37 - t94)) - m(5) * (g(1) * (t8 + t51) + g(2) * (t7 + t52) + g(3) * (t101 + t75)) - m(6) * (g(1) * (t53 + t8) + g(2) * (t54 + t7) + g(3) * (t58 + t75)) - m(7) * (g(1) * (t8 + t68) + g(2) * (t7 + t72) + g(3) * (t57 + t75) + t49) (-m(4) - m(5) + t98) * (g(1) * t45 - g(2) * t48) -m(5) * (g(1) * t51 + g(2) * t52 + g(3) * t101) - m(6) * (g(1) * t53 + g(2) * t54 + g(3) * t58) - m(7) * (g(1) * t68 + g(2) * t72 + g(3) * t57 + t49) t98 * (t99 * t33 - t95) -m(7) * (g(1) * (t2 * rSges(7,1) - t3 * rSges(7,2)) + g(2) * (t4 * rSges(7,1) + t5 * rSges(7,2)) + (-rSges(7,1) * t46 + rSges(7,2) * t43) * t95)];
taug  = t1(:);
