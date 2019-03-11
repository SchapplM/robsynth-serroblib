% Calculate Gravitation load on the joints for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRPR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:24:49
% EndTime: 2019-03-08 19:24:51
% DurationCPUTime: 0.94s
% Computational Cost: add. (549->130), mult. (1171->205), div. (0->0), fcn. (1438->14), ass. (0->63)
t40 = sin(pkin(11));
t48 = sin(qJ(2));
t51 = cos(qJ(2));
t73 = cos(pkin(11));
t26 = -t51 * t40 - t48 * t73;
t41 = sin(pkin(10));
t43 = cos(pkin(10));
t44 = cos(pkin(6));
t77 = t44 * t51;
t108 = -t41 * t48 + t43 * t77;
t55 = -t48 * t40 + t51 * t73;
t74 = t26 * t44;
t15 = t41 * t74 + t43 * t55;
t47 = sin(qJ(4));
t42 = sin(pkin(6));
t50 = cos(qJ(4));
t82 = t42 * t50;
t107 = -t15 * t47 + t41 * t82;
t24 = t26 * t42;
t106 = t24 * t47 + t44 * t50;
t56 = -t41 * t77 - t43 * t48;
t103 = t56 * pkin(2);
t52 = t44 * t55;
t14 = t43 * t26 - t41 * t52;
t36 = t50 * pkin(4) + pkin(3);
t45 = -qJ(5) - pkin(8);
t105 = t14 * t36 - t15 * t45 + t103;
t10 = -t41 * t55 + t43 * t74;
t57 = t10 * t47 - t43 * t82;
t104 = pkin(4) * t57;
t102 = rSges(5,3) + pkin(8);
t101 = rSges(7,3) + pkin(9);
t11 = t41 * t26 + t43 * t52;
t23 = t55 * t42;
t100 = -g(1) * t14 - g(2) * t11 - g(3) * t23;
t99 = -m(6) - m(7);
t39 = qJ(4) + pkin(12);
t38 = cos(t39);
t94 = t38 * pkin(5);
t46 = sin(qJ(6));
t89 = t38 * t46;
t49 = cos(qJ(6));
t88 = t38 * t49;
t86 = t41 * t42;
t84 = t42 * t43;
t83 = t42 * t47;
t79 = t44 * t48;
t33 = t42 * t51 * pkin(2);
t68 = t23 * t36 + t24 * t45 + t33;
t67 = -m(4) - m(5) + t99;
t64 = t108 * pkin(2);
t63 = t107 * pkin(4);
t62 = t106 * pkin(4);
t37 = sin(t39);
t61 = rSges(6,1) * t38 - rSges(6,2) * t37;
t58 = t10 * t45 + t11 * t36 + t64;
t17 = -t24 * t38 + t44 * t37;
t16 = t24 * t37 + t44 * t38;
t5 = t15 * t38 + t37 * t86;
t4 = -t15 * t37 + t38 * t86;
t3 = -t10 * t38 - t37 * t84;
t2 = t10 * t37 - t38 * t84;
t1 = [(-m(2) - m(3) + t67) * g(3), -m(3) * (g(1) * (t56 * rSges(3,1) + (t41 * t79 - t43 * t51) * rSges(3,2)) + g(2) * (t108 * rSges(3,1) + (-t41 * t51 - t43 * t79) * rSges(3,2)) + g(3) * (rSges(3,1) * t51 - rSges(3,2) * t48) * t42) - m(4) * (g(1) * (t14 * rSges(4,1) - rSges(4,2) * t15 + t103) + g(2) * (t11 * rSges(4,1) + t10 * rSges(4,2) + t64) + g(3) * (t23 * rSges(4,1) + t24 * rSges(4,2) + t33)) - m(5) * (g(1) * (t102 * t15 + t103) + g(2) * (-t102 * t10 + t64) + g(3) * (-t102 * t24 + t33) - t100 * (t50 * rSges(5,1) - t47 * rSges(5,2) + pkin(3))) - m(6) * (g(1) * (rSges(6,3) * t15 + t14 * t61 + t105) + g(2) * (-t10 * rSges(6,3) + t11 * t61 + t58) + g(3) * (-t24 * rSges(6,3) + t61 * t23 + t68)) - m(7) * (g(1) * (t14 * t94 + (t14 * t88 + t15 * t46) * rSges(7,1) + (-t14 * t89 + t15 * t49) * rSges(7,2) + t105) + g(2) * (t11 * t94 + (-t10 * t46 + t11 * t88) * rSges(7,1) + (-t10 * t49 - t11 * t89) * rSges(7,2) + t58) + g(3) * (t23 * t94 + (t23 * t88 - t24 * t46) * rSges(7,1) + (-t23 * t89 - t24 * t49) * rSges(7,2) + t68) - t100 * t37 * t101) t67 * (g(3) * t44 + (g(1) * t41 - g(2) * t43) * t42) -m(5) * (g(1) * (t107 * rSges(5,1) + (-t15 * t50 - t41 * t83) * rSges(5,2)) + g(2) * (t57 * rSges(5,1) + (t10 * t50 + t43 * t83) * rSges(5,2)) + g(3) * (t106 * rSges(5,1) + (t24 * t50 - t44 * t47) * rSges(5,2))) - m(6) * (g(1) * (t4 * rSges(6,1) - t5 * rSges(6,2) + t63) + g(2) * (t2 * rSges(6,1) - t3 * rSges(6,2) + t104) + g(3) * (t16 * rSges(6,1) - t17 * rSges(6,2) + t62)) + (-g(1) * (t101 * t5 + t63) - g(2) * (t101 * t3 + t104) - g(3) * (t101 * t17 + t62) - (g(1) * t4 + g(2) * t2 + g(3) * t16) * (t49 * rSges(7,1) - t46 * rSges(7,2) + pkin(5))) * m(7), t99 * t100, -m(7) * (g(1) * ((-t14 * t49 - t5 * t46) * rSges(7,1) + (t14 * t46 - t5 * t49) * rSges(7,2)) + g(2) * ((-t11 * t49 - t3 * t46) * rSges(7,1) + (t11 * t46 - t3 * t49) * rSges(7,2)) + g(3) * ((-t17 * t46 - t23 * t49) * rSges(7,1) + (-t17 * t49 + t23 * t46) * rSges(7,2)))];
taug  = t1(:);
