% Calculate Gravitation load on the joints for
% S5RRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR11_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR11_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR11_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR11_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:38:48
% EndTime: 2019-12-31 22:38:52
% DurationCPUTime: 0.94s
% Computational Cost: add. (529->142), mult. (1144->213), div. (0->0), fcn. (1364->12), ass. (0->65)
t44 = sin(qJ(2));
t45 = sin(qJ(1));
t48 = cos(qJ(2));
t71 = cos(pkin(5));
t80 = cos(qJ(1));
t63 = t71 * t80;
t24 = t44 * t63 + t45 * t48;
t43 = sin(qJ(3));
t47 = cos(qJ(3));
t41 = sin(pkin(5));
t68 = t41 * t80;
t12 = t24 * t47 - t43 * t68;
t23 = t44 * t45 - t48 * t63;
t42 = sin(qJ(4));
t46 = cos(qJ(4));
t100 = t12 * t42 - t23 * t46;
t99 = -t12 * t46 - t23 * t42;
t65 = t45 * t71;
t26 = -t44 * t65 + t80 * t48;
t98 = g(1) * t26 + g(2) * t24;
t25 = t80 * t44 + t48 * t65;
t97 = g(1) * t25 + g(2) * t23;
t75 = t41 * t47;
t15 = t26 * t43 - t45 * t75;
t77 = t41 * t44;
t21 = -t43 * t77 + t71 * t47;
t66 = -t24 * t43 - t47 * t68;
t96 = -g(1) * t15 + g(2) * t66 + g(3) * t21;
t76 = t41 * t45;
t16 = t26 * t47 + t43 * t76;
t22 = t71 * t43 + t44 * t75;
t95 = g(1) * t16 + g(2) * t12 + g(3) * t22;
t36 = pkin(4) * t46 + pkin(3);
t40 = qJ(4) + qJ(5);
t37 = sin(t40);
t38 = cos(t40);
t56 = rSges(6,1) * t38 - rSges(6,2) * t37 + t36;
t72 = pkin(10) + pkin(9) + rSges(6,3);
t94 = t72 * t43 + t56 * t47;
t59 = rSges(5,1) * t46 - rSges(5,2) * t42 + pkin(3);
t81 = pkin(9) + rSges(5,3);
t93 = t81 * t43 + t59 * t47;
t84 = g(3) * t41;
t83 = t42 * pkin(4);
t82 = rSges(4,3) + pkin(8);
t74 = t41 * t48;
t73 = t80 * pkin(1) + pkin(7) * t76;
t70 = t26 * pkin(2) + t73;
t69 = g(3) * (pkin(2) * t74 + pkin(8) * t77);
t67 = -t45 * pkin(1) + pkin(7) * t68;
t64 = -t24 * pkin(2) + t67;
t62 = rSges(4,1) * t47 - rSges(4,2) * t43;
t61 = t42 * rSges(5,1) + t46 * rSges(5,2);
t7 = -t16 * t42 + t25 * t46;
t57 = -t22 * t42 - t46 * t74;
t55 = t37 * rSges(6,1) + t38 * rSges(6,2) + t83;
t54 = pkin(8) + t55;
t17 = t23 * pkin(2);
t19 = t25 * pkin(2);
t53 = -g(1) * t19 - g(2) * t17 + t69;
t5 = -t16 * t37 + t25 * t38;
t6 = t16 * t38 + t25 * t37;
t50 = m(6) * (g(1) * (t5 * rSges(6,1) - t6 * rSges(6,2)) + g(2) * ((-t12 * t37 + t23 * t38) * rSges(6,1) + (-t12 * t38 - t23 * t37) * rSges(6,2)) + g(3) * ((-t22 * t37 - t38 * t74) * rSges(6,1) + (-t22 * t38 + t37 * t74) * rSges(6,2)));
t8 = t16 * t46 + t25 * t42;
t1 = [-m(2) * (g(1) * (-t45 * rSges(2,1) - t80 * rSges(2,2)) + g(2) * (t80 * rSges(2,1) - t45 * rSges(2,2))) - m(3) * (g(1) * (-t24 * rSges(3,1) + t23 * rSges(3,2) + rSges(3,3) * t68 + t67) + g(2) * (rSges(3,1) * t26 - rSges(3,2) * t25 + rSges(3,3) * t76 + t73)) - m(4) * (g(1) * (-rSges(4,1) * t12 - rSges(4,2) * t66 - t82 * t23 + t64) + g(2) * (rSges(4,1) * t16 - rSges(4,2) * t15 + t82 * t25 + t70)) - m(5) * (g(1) * (t99 * rSges(5,1) + t100 * rSges(5,2) - t12 * pkin(3) - t23 * pkin(8) + t81 * t66 + t64) + g(2) * (rSges(5,1) * t8 + rSges(5,2) * t7 + pkin(3) * t16 + t25 * pkin(8) + t81 * t15 + t70)) - m(6) * (g(1) * (-t12 * t56 - t54 * t23 + t66 * t72 + t64) + g(2) * (t6 * rSges(6,1) + t5 * rSges(6,2) + t16 * t36 + (pkin(8) + t83) * t25 + t72 * t15 + t70)), -m(3) * (g(1) * (-rSges(3,1) * t25 - rSges(3,2) * t26) + g(2) * (-rSges(3,1) * t23 - rSges(3,2) * t24) + (rSges(3,1) * t48 - rSges(3,2) * t44) * t84) - m(4) * (g(1) * (-t62 * t25 + t82 * t26 - t19) + g(2) * (-t62 * t23 + t82 * t24 - t17) + t69 + (rSges(4,3) * t44 + t62 * t48) * t84) - m(5) * ((t61 * t44 + t93 * t48) * t84 + t53 + t98 * (pkin(8) + t61) - t97 * t93) - m(6) * ((t55 * t44 + t94 * t48) * t84 + t53 + t98 * t54 - t97 * t94), -m(4) * (g(1) * (-rSges(4,1) * t15 - rSges(4,2) * t16) + g(2) * (rSges(4,1) * t66 - rSges(4,2) * t12) + g(3) * (rSges(4,1) * t21 - rSges(4,2) * t22)) - m(5) * (t96 * t59 + t95 * t81) - m(6) * (t96 * t56 + t95 * t72), -m(5) * (g(1) * (rSges(5,1) * t7 - rSges(5,2) * t8) + g(2) * (-rSges(5,1) * t100 + t99 * rSges(5,2)) + g(3) * (t57 * rSges(5,1) + (-t22 * t46 + t42 * t74) * rSges(5,2))) - t50 - m(6) * (g(1) * t7 - g(2) * t100 + g(3) * t57) * pkin(4), -t50];
taug = t1(:);
