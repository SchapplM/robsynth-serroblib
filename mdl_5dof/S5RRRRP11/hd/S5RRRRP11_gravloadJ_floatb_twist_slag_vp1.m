% Calculate Gravitation load on the joints for
% S5RRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP11_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP11_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP11_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP11_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:14:26
% EndTime: 2019-12-31 22:14:29
% DurationCPUTime: 0.93s
% Computational Cost: add. (510->157), mult. (1249->235), div. (0->0), fcn. (1511->10), ass. (0->71)
t62 = sin(qJ(2));
t65 = cos(qJ(2));
t84 = cos(pkin(5));
t98 = cos(qJ(1));
t75 = t84 * t98;
t97 = sin(qJ(1));
t42 = t62 * t75 + t65 * t97;
t61 = sin(qJ(3));
t64 = cos(qJ(3));
t59 = sin(pkin(5));
t80 = t59 * t98;
t18 = t42 * t64 - t61 * t80;
t41 = t62 * t97 - t65 * t75;
t60 = sin(qJ(4));
t63 = cos(qJ(4));
t1 = t18 * t60 - t41 * t63;
t2 = t18 * t63 + t41 * t60;
t101 = rSges(6,2) + pkin(9);
t102 = rSges(6,1) + pkin(4);
t85 = rSges(6,3) + qJ(5);
t67 = t102 * t63 + t60 * t85;
t104 = pkin(3) * t64;
t103 = g(3) * t59;
t100 = rSges(4,3) + pkin(8);
t99 = rSges(5,3) + pkin(9);
t95 = t41 * t61;
t74 = t84 * t97;
t43 = t62 * t98 + t65 * t74;
t93 = t43 * t61;
t92 = t59 * t62;
t91 = t59 * t65;
t90 = t60 * t64;
t89 = t63 * t64;
t88 = t64 * t65;
t79 = t59 * t97;
t87 = t98 * pkin(1) + pkin(7) * t79;
t86 = pkin(2) * t91 + pkin(8) * t92;
t83 = t61 * t91;
t82 = t60 * t91;
t44 = -t62 * t74 + t65 * t98;
t81 = t44 * pkin(2) + t87;
t78 = -t42 * t61 - t64 * t80;
t77 = t59 * pkin(3) * t88 + pkin(9) * t83 + t86;
t76 = -pkin(1) * t97 + pkin(7) * t80;
t73 = rSges(4,1) * t64 - rSges(4,2) * t61;
t72 = rSges(5,1) * t63 - rSges(5,2) * t60;
t35 = t41 * pkin(2);
t71 = t42 * pkin(8) - pkin(9) * t95 - t104 * t41 - t35;
t37 = t43 * pkin(2);
t70 = t44 * pkin(8) - pkin(9) * t93 - t104 * t43 - t37;
t69 = -pkin(2) * t42 + t76;
t22 = t44 * t64 + t61 * t79;
t68 = pkin(3) * t22 + t43 * pkin(8) + t81;
t66 = -pkin(3) * t18 - t41 * pkin(8) + t69;
t40 = t61 * t84 + t64 * t92;
t39 = -t61 * t92 + t64 * t84;
t34 = t39 * pkin(3);
t24 = (t60 * t62 + t63 * t88) * t59;
t23 = -t63 * t92 + t64 * t82;
t21 = t44 * t61 - t64 * t79;
t16 = t40 * t63 - t82;
t15 = t40 * t60 + t63 * t91;
t13 = t21 * pkin(3);
t11 = t78 * pkin(3);
t10 = -t43 * t89 + t44 * t60;
t9 = -t43 * t90 - t44 * t63;
t8 = -t41 * t89 + t42 * t60;
t7 = -t41 * t90 - t42 * t63;
t6 = t22 * t63 + t43 * t60;
t5 = t22 * t60 - t43 * t63;
t3 = [-m(2) * (g(1) * (-rSges(2,1) * t97 - rSges(2,2) * t98) + g(2) * (rSges(2,1) * t98 - rSges(2,2) * t97)) - m(3) * (g(1) * (-rSges(3,1) * t42 + rSges(3,2) * t41 + rSges(3,3) * t80 + t76) + g(2) * (rSges(3,1) * t44 - rSges(3,2) * t43 + rSges(3,3) * t79 + t87)) - m(4) * (g(1) * (-rSges(4,1) * t18 - rSges(4,2) * t78 - t100 * t41 + t69) + g(2) * (rSges(4,1) * t22 - rSges(4,2) * t21 + t100 * t43 + t81)) - m(5) * (g(1) * (-rSges(5,1) * t2 + rSges(5,2) * t1 + t78 * t99 + t66) + g(2) * (rSges(5,1) * t6 - rSges(5,2) * t5 + t21 * t99 + t68)) - m(6) * (g(1) * (-t1 * t85 + t101 * t78 - t102 * t2 + t66) + g(2) * (t101 * t21 + t102 * t6 + t5 * t85 + t68)), -m(3) * (g(1) * (-rSges(3,1) * t43 - rSges(3,2) * t44) + g(2) * (-rSges(3,1) * t41 - rSges(3,2) * t42) + (rSges(3,1) * t65 - rSges(3,2) * t62) * t103) - m(4) * (g(1) * (t100 * t44 - t43 * t73 - t37) + g(2) * (t100 * t42 - t41 * t73 - t35) + g(3) * t86 + (rSges(4,3) * t62 + t65 * t73) * t103) - m(5) * (g(1) * (rSges(5,1) * t10 - rSges(5,2) * t9 - rSges(5,3) * t93 + t70) + g(2) * (rSges(5,1) * t8 - rSges(5,2) * t7 - rSges(5,3) * t95 + t71) + g(3) * (rSges(5,1) * t24 - rSges(5,2) * t23 + rSges(5,3) * t83 + t77)) - m(6) * (g(1) * (-rSges(6,2) * t93 + t10 * t102 + t85 * t9 + t70) + g(2) * (-rSges(6,2) * t95 + t102 * t8 + t7 * t85 + t71) + g(3) * (rSges(6,2) * t83 + t102 * t24 + t23 * t85 + t77)), -m(4) * (g(1) * (-rSges(4,1) * t21 - rSges(4,2) * t22) + g(2) * (rSges(4,1) * t78 - rSges(4,2) * t18) + g(3) * (rSges(4,1) * t39 - rSges(4,2) * t40)) - m(5) * (g(1) * (-t21 * t72 + t22 * t99 - t13) + g(2) * (t18 * t99 + t72 * t78 + t11) + g(3) * (t39 * t72 + t40 * t99 + t34)) - m(6) * ((t101 * t40 + t67 * t39 + t34) * g(3) + (t101 * t18 + t67 * t78 + t11) * g(2) + (t101 * t22 - t21 * t67 - t13) * g(1)), -m(5) * (g(1) * (-rSges(5,1) * t5 - rSges(5,2) * t6) + g(2) * (-rSges(5,1) * t1 - rSges(5,2) * t2) + g(3) * (-rSges(5,1) * t15 - rSges(5,2) * t16)) - m(6) * (g(1) * (-t102 * t5 + t6 * t85) + g(2) * (-t1 * t102 + t2 * t85) + g(3) * (-t102 * t15 + t16 * t85)), -m(6) * (g(1) * t5 + g(2) * t1 + g(3) * t15)];
taug = t3(:);
