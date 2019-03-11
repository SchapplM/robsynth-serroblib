% Calculate Gravitation load on the joints for
% S6PPRRRP1
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
% Datum: 2019-03-08 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PPRRRP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRRP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:52:21
% EndTime: 2019-03-08 18:52:23
% DurationCPUTime: 0.73s
% Computational Cost: add. (809->123), mult. (2204->197), div. (0->0), fcn. (2814->14), ass. (0->79)
t73 = sin(pkin(12));
t74 = sin(pkin(11));
t59 = t74 * t73;
t77 = cos(pkin(12));
t78 = cos(pkin(11));
t66 = t78 * t77;
t80 = cos(pkin(6));
t52 = -t80 * t66 + t59;
t75 = sin(pkin(7));
t76 = sin(pkin(6));
t63 = t76 * t75;
t79 = cos(pkin(7));
t104 = t52 * t79 + t78 * t63;
t60 = t74 * t77;
t64 = t78 * t73;
t53 = t80 * t60 + t64;
t62 = t76 * t74;
t103 = t53 * t79 - t75 * t62;
t102 = t77 * t79 * t76 + t80 * t75;
t35 = -t80 * t59 + t66;
t44 = sin(qJ(3));
t88 = cos(qJ(3));
t22 = -t103 * t44 + t35 * t88;
t43 = sin(qJ(4));
t46 = cos(qJ(4));
t48 = t53 * t75 + t79 * t62;
t11 = t22 * t43 - t48 * t46;
t61 = t76 * t73;
t29 = t102 * t44 + t88 * t61;
t51 = -t77 * t63 + t80 * t79;
t23 = t29 * t43 - t51 * t46;
t34 = t80 * t64 + t60;
t20 = -t104 * t44 + t34 * t88;
t65 = t78 * t76;
t47 = t52 * t75 - t79 * t65;
t9 = t20 * t43 - t47 * t46;
t101 = g(1) * t11 + g(2) * t9 + g(3) * t23;
t10 = t20 * t46 + t47 * t43;
t12 = t22 * t46 + t48 * t43;
t24 = t29 * t46 + t51 * t43;
t100 = g(1) * t12 + g(2) * t10 + g(3) * t24;
t19 = t104 * t88 + t34 * t44;
t21 = t103 * t88 + t35 * t44;
t28 = -t102 * t88 + t44 * t61;
t99 = (-g(1) * t21 - g(2) * t19 - g(3) * t28) * t43;
t98 = pkin(4) * t46;
t90 = rSges(5,3) + pkin(9);
t89 = rSges(6,3) + pkin(10);
t42 = sin(qJ(5));
t87 = t20 * t42;
t86 = t22 * t42;
t85 = t29 * t42;
t45 = cos(qJ(5));
t40 = pkin(5) * t45 + pkin(4);
t84 = t40 * t46;
t83 = t42 * t46;
t82 = t45 * t46;
t81 = rSges(7,3) + qJ(6) + pkin(10);
t17 = t19 * pkin(3);
t72 = t20 * pkin(9) - t17;
t18 = t21 * pkin(3);
t71 = t22 * pkin(9) - t18;
t27 = t28 * pkin(3);
t70 = t29 * pkin(9) - t27;
t69 = -m(3) - m(4) - m(5) - m(6) - m(7);
t68 = -rSges(5,1) * t46 + rSges(5,2) * t43;
t1 = -t10 * t42 + t19 * t45;
t3 = -t12 * t42 + t21 * t45;
t13 = -t24 * t42 + t28 * t45;
t16 = -t28 * t82 + t85;
t15 = t28 * t83 + t29 * t45;
t14 = -t24 * t45 - t28 * t42;
t8 = -t21 * t82 + t86;
t7 = t21 * t83 + t22 * t45;
t6 = -t19 * t82 + t87;
t5 = t19 * t83 + t20 * t45;
t4 = -t12 * t45 - t21 * t42;
t2 = -t10 * t45 - t19 * t42;
t25 = [(-m(2) + t69) * g(3), t69 * (g(1) * t62 - g(2) * t65 + g(3) * t80) -m(4) * (g(1) * (-rSges(4,1) * t21 - rSges(4,2) * t22) + g(2) * (-rSges(4,1) * t19 - rSges(4,2) * t20) + g(3) * (-rSges(4,1) * t28 - rSges(4,2) * t29)) - m(5) * (g(1) * (t68 * t21 + t90 * t22 - t18) + g(2) * (t68 * t19 + t90 * t20 - t17) + g(3) * (t68 * t28 + t90 * t29 - t27)) + (-g(1) * (t8 * rSges(7,1) + t7 * rSges(7,2) + pkin(5) * t86 - t21 * t84 + t71) - g(2) * (t6 * rSges(7,1) + t5 * rSges(7,2) + pkin(5) * t87 - t19 * t84 + t72) - g(3) * (t16 * rSges(7,1) + t15 * rSges(7,2) + pkin(5) * t85 - t28 * t84 + t70) - t81 * t99) * m(7) + (-g(1) * (t8 * rSges(6,1) + t7 * rSges(6,2) - t21 * t98 + t71) - g(2) * (t6 * rSges(6,1) + t5 * rSges(6,2) - t19 * t98 + t72) - g(3) * (t16 * rSges(6,1) + t15 * rSges(6,2) - t28 * t98 + t70) - t89 * t99) * m(6), -m(5) * (g(1) * (-rSges(5,1) * t11 - rSges(5,2) * t12) + g(2) * (-rSges(5,1) * t9 - rSges(5,2) * t10) + g(3) * (-rSges(5,1) * t23 - rSges(5,2) * t24)) - m(6) * (t100 * t89 + t101 * (-rSges(6,1) * t45 + rSges(6,2) * t42 - pkin(4))) - m(7) * (t100 * t81 + t101 * (-rSges(7,1) * t45 + rSges(7,2) * t42 - t40)) -m(6) * (g(1) * (rSges(6,1) * t3 + rSges(6,2) * t4) + g(2) * (rSges(6,1) * t1 + rSges(6,2) * t2) + g(3) * (rSges(6,1) * t13 + rSges(6,2) * t14)) + (-g(1) * (t3 * rSges(7,1) + t4 * rSges(7,2)) - g(2) * (t1 * rSges(7,1) + t2 * rSges(7,2)) - g(3) * (t13 * rSges(7,1) + t14 * rSges(7,2)) - (g(1) * t3 + g(2) * t1 + g(3) * t13) * pkin(5)) * m(7), -m(7) * t101];
taug  = t25(:);
