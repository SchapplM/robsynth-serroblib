% Calculate Gravitation load on the joints for
% S6PRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRPR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:34:30
% EndTime: 2019-03-08 19:34:32
% DurationCPUTime: 0.83s
% Computational Cost: add. (521->123), mult. (1307->186), div. (0->0), fcn. (1624->12), ass. (0->62)
t42 = sin(pkin(11));
t48 = sin(qJ(2));
t51 = cos(qJ(2));
t74 = cos(pkin(11));
t32 = -t51 * t42 - t48 * t74;
t43 = sin(pkin(10));
t44 = cos(pkin(10));
t45 = cos(pkin(6));
t56 = -t48 * t42 + t51 * t74;
t53 = t56 * t45;
t18 = t44 * t32 - t43 * t53;
t82 = t45 * t51;
t58 = -t43 * t82 - t44 * t48;
t54 = t58 * pkin(2);
t52 = t18 * pkin(3) + t54;
t47 = sin(qJ(4));
t50 = cos(qJ(4));
t97 = pkin(4) * t50 + qJ(5) * t47;
t102 = t97 * t18 + t52;
t46 = sin(qJ(6));
t49 = cos(qJ(6));
t101 = t46 * rSges(7,1) + t49 * rSges(7,2);
t98 = -t43 * t48 + t44 * t82;
t96 = rSges(7,3) + pkin(9);
t78 = t32 * t45;
t19 = t43 * t78 + t44 * t56;
t14 = -t43 * t56 + t44 * t78;
t94 = t101 * t47;
t93 = t49 * rSges(7,1) - t46 * rSges(7,2) + pkin(5) + pkin(8);
t92 = -m(6) - m(7);
t90 = rSges(6,1) + pkin(8);
t89 = rSges(5,3) + pkin(8);
t83 = t45 * t48;
t73 = sin(pkin(6));
t59 = t74 * t73;
t67 = t48 * t73;
t29 = t42 * t67 - t51 * t59;
t65 = t51 * t73;
t39 = pkin(2) * t65;
t77 = -t29 * pkin(3) + t39;
t75 = rSges(6,3) + qJ(5);
t70 = -m(4) - m(5) + t92;
t69 = t47 * t73;
t66 = t50 * t73;
t64 = -t97 * t29 + t77;
t63 = t98 * pkin(2);
t62 = rSges(5,1) * t50 - rSges(5,2) * t47;
t61 = rSges(6,2) * t50 - rSges(6,3) * t47;
t15 = t43 * t32 + t44 * t53;
t60 = t15 * pkin(3) + t63;
t55 = t97 * t15 + t60;
t30 = t42 * t65 + t48 * t59;
t22 = t30 * t50 + t45 * t47;
t21 = t30 * t47 - t45 * t50;
t20 = t21 * pkin(4);
t7 = t19 * t50 + t43 * t69;
t6 = t19 * t47 - t43 * t66;
t5 = -t14 * t50 - t44 * t69;
t4 = -t14 * t47 + t44 * t66;
t3 = t6 * pkin(4);
t2 = t4 * pkin(4);
t1 = [(-m(2) - m(3) + t70) * g(3), -m(3) * (g(1) * (t58 * rSges(3,1) + (t43 * t83 - t44 * t51) * rSges(3,2)) + g(2) * (t98 * rSges(3,1) + (-t43 * t51 - t44 * t83) * rSges(3,2)) + g(3) * (rSges(3,1) * t65 - rSges(3,2) * t67)) - m(4) * (g(1) * (t18 * rSges(4,1) - rSges(4,2) * t19 + t54) + g(2) * (t15 * rSges(4,1) + t14 * rSges(4,2) + t63) + g(3) * (-t29 * rSges(4,1) - t30 * rSges(4,2) + t39)) - m(5) * (g(1) * (t62 * t18 + t19 * t89 + t52) + g(2) * (-t89 * t14 + t62 * t15 + t60) + g(3) * (-t62 * t29 + t89 * t30 + t77)) - m(6) * (g(1) * (-t61 * t18 + t19 * t90 + t102) + g(2) * (-t90 * t14 - t61 * t15 + t55) + g(3) * (t61 * t29 + t90 * t30 + t64)) - m(7) * (g(1) * (t94 * t18 + t19 * t93 + t102) + g(2) * (-t14 * t93 + t15 * t94 + t55) + g(3) * (-t29 * t94 + t30 * t93 + t64) + (g(1) * t18 + g(2) * t15 - g(3) * t29) * t50 * t96) t70 * (g(3) * t45 + (g(1) * t43 - g(2) * t44) * t73) -m(5) * (g(1) * (-t6 * rSges(5,1) - t7 * rSges(5,2)) + g(2) * (-t4 * rSges(5,1) - t5 * rSges(5,2)) + g(3) * (-t21 * rSges(5,1) - t22 * rSges(5,2))) - m(6) * (g(1) * (t6 * rSges(6,2) + t75 * t7 - t3) + g(2) * (t4 * rSges(6,2) + t75 * t5 - t2) + g(3) * (t21 * rSges(6,2) + t75 * t22 - t20)) + (-g(1) * (-t96 * t6 - t3) - g(2) * (-t96 * t4 - t2) - g(3) * (-t96 * t21 - t20) - (g(1) * t7 + g(2) * t5 + g(3) * t22) * (qJ(5) + t101)) * m(7), t92 * (g(1) * t6 + g(2) * t4 + g(3) * t21) -m(7) * (g(1) * ((t18 * t46 + t6 * t49) * rSges(7,1) + (t18 * t49 - t6 * t46) * rSges(7,2)) + g(2) * ((t15 * t46 + t4 * t49) * rSges(7,1) + (t15 * t49 - t4 * t46) * rSges(7,2)) + g(3) * ((t21 * t49 - t29 * t46) * rSges(7,1) + (-t21 * t46 - t29 * t49) * rSges(7,2)))];
taug  = t1(:);
