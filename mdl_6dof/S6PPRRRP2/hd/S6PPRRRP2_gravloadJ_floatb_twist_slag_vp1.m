% Calculate Gravitation load on the joints for
% S6PPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-03-08 18:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PPRRRP2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRRP2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:56:03
% EndTime: 2019-03-08 18:56:05
% DurationCPUTime: 0.55s
% Computational Cost: add. (925->130), mult. (2551->206), div. (0->0), fcn. (3280->14), ass. (0->80)
t105 = rSges(7,2) + pkin(10);
t83 = sin(pkin(12));
t84 = sin(pkin(11));
t68 = t84 * t83;
t87 = cos(pkin(12));
t88 = cos(pkin(11));
t75 = t88 * t87;
t90 = cos(pkin(6));
t62 = -t90 * t75 + t68;
t85 = sin(pkin(7));
t86 = sin(pkin(6));
t72 = t86 * t85;
t89 = cos(pkin(7));
t104 = t62 * t89 + t88 * t72;
t69 = t84 * t87;
t73 = t88 * t83;
t63 = t90 * t69 + t73;
t71 = t86 * t84;
t103 = t63 * t89 - t85 * t71;
t102 = t87 * t89 * t86 + t90 * t85;
t59 = cos(qJ(4));
t101 = pkin(4) * t59;
t100 = rSges(7,1) + pkin(5);
t99 = rSges(5,3) + pkin(9);
t98 = rSges(6,3) + pkin(10);
t97 = cos(qJ(3));
t49 = t90 * t73 + t69;
t57 = sin(qJ(3));
t30 = t104 * t97 + t49 * t57;
t56 = sin(qJ(4));
t96 = t30 * t56;
t50 = -t90 * t68 + t75;
t32 = t103 * t97 + t50 * t57;
t95 = t32 * t56;
t70 = t86 * t83;
t41 = -t102 * t97 + t57 * t70;
t94 = t41 * t56;
t55 = sin(qJ(5));
t93 = t55 * t59;
t58 = cos(qJ(5));
t92 = t58 * t59;
t91 = rSges(7,3) + qJ(6);
t82 = -m(3) - m(4) - m(5) - m(6) - m(7);
t81 = -rSges(5,1) * t59 + rSges(5,2) * t56;
t80 = rSges(6,1) * t58 - rSges(6,2) * t55;
t27 = t30 * pkin(3);
t31 = -t104 * t57 + t49 * t97;
t79 = t31 * pkin(9) - pkin(10) * t96 - t30 * t101 - t27;
t28 = t32 * pkin(3);
t33 = -t103 * t57 + t50 * t97;
t78 = t33 * pkin(9) - pkin(10) * t95 - t32 * t101 - t28;
t40 = t41 * pkin(3);
t42 = t102 * t57 + t97 * t70;
t77 = t42 * pkin(9) - pkin(10) * t94 - t41 * t101 - t40;
t74 = t88 * t86;
t48 = -t87 * t72 + t90 * t89;
t44 = t63 * t85 + t89 * t71;
t43 = t62 * t85 - t89 * t74;
t35 = t42 * t59 + t48 * t56;
t34 = -t42 * t56 + t48 * t59;
t29 = t34 * pkin(4);
t18 = -t41 * t92 + t42 * t55;
t17 = -t41 * t93 - t42 * t58;
t16 = t35 * t58 + t41 * t55;
t15 = t35 * t55 - t41 * t58;
t14 = t33 * t59 + t44 * t56;
t13 = -t33 * t56 + t44 * t59;
t12 = t31 * t59 + t43 * t56;
t11 = -t31 * t56 + t43 * t59;
t10 = t13 * pkin(4);
t9 = t11 * pkin(4);
t8 = -t32 * t92 + t33 * t55;
t7 = -t32 * t93 - t33 * t58;
t6 = -t30 * t92 + t31 * t55;
t5 = -t30 * t93 - t31 * t58;
t4 = t14 * t58 + t32 * t55;
t3 = t14 * t55 - t32 * t58;
t2 = t12 * t58 + t30 * t55;
t1 = t12 * t55 - t30 * t58;
t19 = [(-m(2) + t82) * g(3), t82 * (g(1) * t71 - g(2) * t74 + g(3) * t90) -m(4) * (g(1) * (-t32 * rSges(4,1) - t33 * rSges(4,2)) + g(2) * (-t30 * rSges(4,1) - t31 * rSges(4,2)) + g(3) * (-t41 * rSges(4,1) - t42 * rSges(4,2))) - m(5) * (g(1) * (t81 * t32 + t99 * t33 - t28) + g(2) * (t81 * t30 + t99 * t31 - t27) + g(3) * (t81 * t41 + t99 * t42 - t40)) - m(6) * (g(1) * (t8 * rSges(6,1) - t7 * rSges(6,2) - rSges(6,3) * t95 + t78) + g(2) * (t6 * rSges(6,1) - t5 * rSges(6,2) - rSges(6,3) * t96 + t79) + g(3) * (t18 * rSges(6,1) - t17 * rSges(6,2) - rSges(6,3) * t94 + t77)) - m(7) * (g(1) * (-rSges(7,2) * t95 + t100 * t8 + t91 * t7 + t78) + g(2) * (-rSges(7,2) * t96 + t100 * t6 + t91 * t5 + t79) + g(3) * (-rSges(7,2) * t94 + t100 * t18 + t91 * t17 + t77)) -m(5) * (g(1) * (t13 * rSges(5,1) - t14 * rSges(5,2)) + g(2) * (t11 * rSges(5,1) - t12 * rSges(5,2)) + g(3) * (t34 * rSges(5,1) - t35 * rSges(5,2))) - m(6) * (g(1) * (t80 * t13 + t98 * t14 + t10) + g(2) * (t80 * t11 + t98 * t12 + t9) + g(3) * (t80 * t34 + t98 * t35 + t29)) + (-g(1) * (t105 * t14 + t10) - g(2) * (t105 * t12 + t9) - g(3) * (t105 * t35 + t29) - (g(1) * t13 + g(2) * t11 + g(3) * t34) * (t100 * t58 + t91 * t55)) * m(7), -m(6) * (g(1) * (-t3 * rSges(6,1) - t4 * rSges(6,2)) + g(2) * (-t1 * rSges(6,1) - t2 * rSges(6,2)) + g(3) * (-t15 * rSges(6,1) - t16 * rSges(6,2))) - m(7) * (g(1) * (-t100 * t3 + t91 * t4) + g(2) * (-t100 * t1 + t91 * t2) + g(3) * (-t100 * t15 + t91 * t16)) -m(7) * (g(1) * t3 + g(2) * t1 + g(3) * t15)];
taug  = t19(:);
