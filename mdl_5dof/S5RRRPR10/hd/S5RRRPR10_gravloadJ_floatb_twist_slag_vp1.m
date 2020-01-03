% Calculate Gravitation load on the joints for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR10_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR10_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR10_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR10_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:26:39
% EndTime: 2019-12-31 21:26:42
% DurationCPUTime: 0.90s
% Computational Cost: add. (476->145), mult. (864->212), div. (0->0), fcn. (1000->12), ass. (0->65)
t47 = sin(qJ(3));
t51 = cos(qJ(3));
t73 = cos(pkin(5));
t44 = sin(pkin(5));
t48 = sin(qJ(2));
t82 = t44 * t48;
t101 = -t47 * t82 + t73 * t51;
t52 = cos(qJ(2));
t49 = sin(qJ(1));
t65 = t49 * t73;
t87 = cos(qJ(1));
t25 = -t48 * t65 + t87 * t52;
t80 = t44 * t51;
t10 = -t25 * t47 + t49 * t80;
t60 = t73 * t87;
t22 = t48 * t49 - t52 * t60;
t46 = sin(qJ(5));
t23 = t48 * t60 + t49 * t52;
t43 = qJ(3) + pkin(10);
t40 = sin(t43);
t41 = cos(t43);
t69 = t44 * t87;
t5 = t23 * t41 - t40 * t69;
t50 = cos(qJ(5));
t100 = -t22 * t50 + t46 * t5;
t99 = -t22 * t46 - t5 * t50;
t79 = t44 * t52;
t98 = g(3) * t79;
t54 = t23 * t47 + t51 * t69;
t97 = t54 * pkin(3);
t88 = pkin(9) + rSges(6,3);
t96 = t46 * rSges(6,1) + t50 * rSges(6,2);
t24 = t87 * t48 + t52 * t65;
t92 = g(2) * t22;
t95 = g(1) * t24 + t92;
t56 = t50 * rSges(6,1) - t46 * rSges(6,2) + pkin(4);
t94 = t88 * t40 + t56 * t41;
t39 = pkin(3) * t51 + pkin(2);
t91 = t39 * t98;
t90 = g(3) * t44;
t89 = -pkin(8) - rSges(4,3);
t81 = t44 * t49;
t45 = -qJ(4) - pkin(8);
t76 = -t22 * t39 - t23 * t45;
t75 = -t24 * t39 - t25 * t45;
t74 = t87 * pkin(1) + pkin(7) * t81;
t72 = t47 * t81;
t68 = -t49 * pkin(1) + pkin(7) * t69;
t67 = -t23 * t40 - t41 * t69;
t34 = t47 * t69;
t66 = -t23 * t51 + t34;
t62 = t10 * pkin(3);
t61 = pkin(3) * t72 - t24 * t45 + t25 * t39 + t74;
t59 = rSges(5,1) * t41 - rSges(5,2) * t40;
t58 = t101 * pkin(3);
t57 = t51 * rSges(4,1) - t47 * rSges(4,2) + pkin(2);
t55 = pkin(3) * t34 + t22 * t45 - t23 * t39 + t68;
t17 = t73 * t40 + t41 * t82;
t16 = -t40 * t82 + t73 * t41;
t11 = t25 * t51 + t72;
t9 = t25 * t41 + t40 * t81;
t8 = t25 * t40 - t41 * t81;
t2 = t24 * t46 + t50 * t9;
t1 = t24 * t50 - t46 * t9;
t3 = [-m(2) * (g(1) * (-t49 * rSges(2,1) - t87 * rSges(2,2)) + g(2) * (t87 * rSges(2,1) - t49 * rSges(2,2))) - m(3) * (g(1) * (-t23 * rSges(3,1) + t22 * rSges(3,2) + rSges(3,3) * t69 + t68) + g(2) * (rSges(3,1) * t25 - rSges(3,2) * t24 + rSges(3,3) * t81 + t74)) - m(4) * (g(1) * (t66 * rSges(4,1) + t54 * rSges(4,2) - t23 * pkin(2) + t89 * t22 + t68) + g(2) * (rSges(4,1) * t11 + rSges(4,2) * t10 + pkin(2) * t25 - t89 * t24 + t74)) - m(5) * (g(1) * (-rSges(5,1) * t5 - rSges(5,2) * t67 - rSges(5,3) * t22 + t55) + g(2) * (rSges(5,1) * t9 - rSges(5,2) * t8 + rSges(5,3) * t24 + t61)) - m(6) * (g(1) * (t99 * rSges(6,1) + t100 * rSges(6,2) - t5 * pkin(4) + t88 * t67 + t55) + g(2) * (rSges(6,1) * t2 + rSges(6,2) * t1 + pkin(4) * t9 + t88 * t8 + t61)), -m(3) * (g(1) * (-rSges(3,1) * t24 - rSges(3,2) * t25) + g(2) * (-rSges(3,1) * t22 - rSges(3,2) * t23) + (rSges(3,1) * t52 - rSges(3,2) * t48) * t90) - m(4) * (g(1) * (-t57 * t24 - t89 * t25) - g(2) * t89 * t23 - t57 * t92 + (-t89 * t48 + t57 * t52) * t90) - m(5) * (g(1) * (rSges(5,3) * t25 - t59 * t24 + t75) + g(2) * (rSges(5,3) * t23 - t59 * t22 + t76) + t91 + (t59 * t52 + (rSges(5,3) - t45) * t48) * t90) - m(6) * (g(1) * (t96 * t25 + t75) + g(2) * (t96 * t23 + t76) + t91 + ((-t45 + t96) * t48 + t94 * t52) * t90 - t95 * t94), -m(4) * (g(1) * (rSges(4,1) * t10 - rSges(4,2) * t11) + g(2) * (-t54 * rSges(4,1) + t66 * rSges(4,2)) + g(3) * (t101 * rSges(4,1) + (-t73 * t47 - t48 * t80) * rSges(4,2))) - m(5) * (g(1) * (-rSges(5,1) * t8 - rSges(5,2) * t9 + t62) + g(2) * (rSges(5,1) * t67 - t5 * rSges(5,2) - t97) + g(3) * (rSges(5,1) * t16 - rSges(5,2) * t17 + t58)) + (-g(1) * (t88 * t9 + t62) - g(2) * (t88 * t5 - t97) - g(3) * (t88 * t17 + t58) - (-g(1) * t8 + g(2) * t67 + g(3) * t16) * t56) * m(6), (-m(5) - m(6)) * (t95 - t98), -m(6) * (g(1) * (rSges(6,1) * t1 - rSges(6,2) * t2) + g(2) * (-t100 * rSges(6,1) + t99 * rSges(6,2)) + g(3) * ((-t17 * t46 - t50 * t79) * rSges(6,1) + (-t17 * t50 + t46 * t79) * rSges(6,2)))];
taug = t3(:);
