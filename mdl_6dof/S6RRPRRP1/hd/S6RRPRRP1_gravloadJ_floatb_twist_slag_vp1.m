% Calculate Gravitation load on the joints for
% S6RRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 11:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:39:00
% EndTime: 2019-03-09 11:39:02
% DurationCPUTime: 0.75s
% Computational Cost: add. (505->136), mult. (457->183), div. (0->0), fcn. (410->10), ass. (0->66)
t95 = rSges(6,3) + pkin(9);
t94 = rSges(7,1) + pkin(5);
t39 = qJ(2) + pkin(10);
t36 = qJ(4) + t39;
t30 = sin(t36);
t31 = cos(t36);
t45 = cos(qJ(5));
t32 = pkin(5) * t45 + pkin(4);
t93 = t30 * rSges(7,3) + t31 * t32;
t92 = t31 * rSges(5,1) - rSges(5,2) * t30;
t34 = sin(t39);
t35 = cos(t39);
t46 = cos(qJ(2));
t37 = t46 * pkin(2);
t91 = rSges(4,1) * t35 - rSges(4,2) * t34 + t37;
t44 = sin(qJ(1));
t47 = cos(qJ(1));
t90 = g(1) * t47 + g(2) * t44;
t52 = t31 * pkin(4) + t95 * t30;
t89 = t90 * t30;
t43 = sin(qJ(2));
t88 = pkin(2) * t43;
t42 = sin(qJ(5));
t87 = pkin(5) * t42;
t84 = rSges(3,3) + pkin(7);
t75 = t44 * t42;
t67 = t30 * t75;
t79 = t31 * t44;
t83 = rSges(7,2) * t67 + rSges(7,3) * t79;
t40 = -qJ(6) - pkin(9);
t81 = t30 * t40;
t80 = t31 * t42;
t78 = t31 * t45;
t77 = t31 * t47;
t76 = t42 * t47;
t74 = t44 * t45;
t73 = t45 * t47;
t41 = -qJ(3) - pkin(7);
t72 = rSges(4,3) - t41;
t38 = -pkin(8) + t41;
t71 = rSges(5,3) - t38;
t66 = t30 * t76;
t70 = rSges(7,2) * t66 + rSges(7,3) * t77;
t69 = pkin(3) * t35 + t37;
t68 = rSges(6,2) * t67 + t95 * t79;
t65 = rSges(6,2) * t66 + t95 * t77;
t64 = -rSges(6,1) * t45 - pkin(4);
t63 = -t38 + t87;
t62 = -rSges(7,1) * t45 - t32;
t60 = rSges(3,1) * t46 - rSges(3,2) * t43;
t56 = -rSges(5,1) * t30 - rSges(5,2) * t31;
t55 = pkin(1) + t60;
t3 = -t31 * t76 + t74;
t1 = t31 * t75 + t73;
t54 = pkin(1) + t91;
t53 = rSges(7,1) * t78 - rSges(7,2) * t80 + t93;
t51 = rSges(6,1) * t78 - rSges(6,2) * t80 + t52;
t49 = -t81 + t93;
t14 = -pkin(3) * t34 - t88;
t13 = pkin(1) + t69;
t7 = t47 * t14;
t6 = t44 * t14;
t5 = t47 * t13;
t4 = t31 * t73 + t75;
t2 = -t31 * t74 + t76;
t8 = [-m(2) * (g(1) * (-t44 * rSges(2,1) - rSges(2,2) * t47) + g(2) * (rSges(2,1) * t47 - t44 * rSges(2,2))) - m(3) * ((g(1) * t84 + g(2) * t55) * t47 + (-g(1) * t55 + g(2) * t84) * t44) - m(4) * ((g(1) * t72 + g(2) * t54) * t47 + (-g(1) * t54 + g(2) * t72) * t44) - m(5) * (g(2) * t5 + (g(1) * t71 + g(2) * t92) * t47 + (g(1) * (-t13 - t92) + g(2) * t71) * t44) - m(6) * (g(1) * (t2 * rSges(6,1) + t1 * rSges(6,2)) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t5) + (-g(1) * t38 + g(2) * t52) * t47 + (g(1) * (-t13 - t52) - g(2) * t38) * t44) - m(7) * (g(1) * (t2 * rSges(7,1) + t1 * rSges(7,2)) + g(2) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t5) + (g(1) * t63 + g(2) * t49) * t47 + (g(1) * (-t13 - t49) + g(2) * t63) * t44) -m(3) * (g(3) * t60 + t90 * (-rSges(3,1) * t43 - rSges(3,2) * t46)) - m(4) * (g(3) * t91 + t90 * (-rSges(4,1) * t34 - rSges(4,2) * t35 - t88)) - m(5) * (g(1) * (t56 * t47 + t7) + g(2) * (t56 * t44 + t6) + g(3) * (t92 + t69)) - m(6) * (g(1) * (t7 + t65) + g(2) * (t6 + t68) + g(3) * (t51 + t69) + t64 * t89) - m(7) * (g(1) * (-t40 * t77 + t7 + t70) + g(2) * (-t40 * t79 + t6 + t83) + g(3) * (t53 + t69) + (-g(3) * t40 + t62 * t90) * t30) (-m(4) - m(5) - m(6) - m(7)) * (g(1) * t44 - g(2) * t47) -m(5) * g(3) * t92 - m(6) * (g(1) * t65 + g(2) * t68 + g(3) * t51) - m(7) * (g(1) * t70 + g(2) * t83 + g(3) * (t53 - t81)) + t90 * ((m(5) * rSges(5,2) + m(7) * t40) * t31 + (m(5) * rSges(5,1) - m(6) * t64 - m(7) * t62) * t30) -m(6) * (g(1) * (rSges(6,1) * t3 - rSges(6,2) * t4) + g(2) * (-rSges(6,1) * t1 + rSges(6,2) * t2)) - m(7) * (g(1) * (-t4 * rSges(7,2) + t3 * t94) + g(2) * (t2 * rSges(7,2) - t1 * t94)) + (-m(6) * (-rSges(6,1) * t42 - rSges(6,2) * t45) - m(7) * (-rSges(7,1) * t42 - rSges(7,2) * t45 - t87)) * g(3) * t30, -m(7) * (-g(3) * t31 + t89)];
taug  = t8(:);
