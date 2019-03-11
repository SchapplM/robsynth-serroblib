% Calculate Gravitation load on the joints for
% S6RRPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 11:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR11_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR11_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR11_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR11_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:11:45
% EndTime: 2019-03-09 11:11:47
% DurationCPUTime: 0.76s
% Computational Cost: add. (374->160), mult. (530->217), div. (0->0), fcn. (498->10), ass. (0->76)
t82 = rSges(5,3) + pkin(8);
t42 = -qJ(5) - pkin(8);
t69 = rSges(6,3) - t42;
t68 = rSges(7,3) + pkin(9) - t42;
t48 = cos(qJ(1));
t45 = sin(qJ(1));
t85 = g(2) * t45;
t55 = g(1) * t48 + t85;
t90 = -m(6) - m(7);
t41 = qJ(4) + pkin(10);
t33 = qJ(6) + t41;
t29 = cos(t33);
t28 = sin(t33);
t77 = t45 * t28;
t44 = sin(qJ(2));
t78 = t44 * t48;
t5 = t29 * t78 - t77;
t76 = t45 * t29;
t6 = t28 * t78 + t76;
t89 = t5 * rSges(7,1) - t6 * rSges(7,2);
t7 = t28 * t48 + t44 * t76;
t8 = t29 * t48 - t44 * t77;
t88 = t7 * rSges(7,1) + t8 * rSges(7,2);
t43 = sin(qJ(4));
t87 = pkin(4) * t43;
t47 = cos(qJ(2));
t84 = g(3) * t47;
t46 = cos(qJ(4));
t36 = t46 * pkin(4);
t37 = t47 * pkin(2);
t81 = rSges(7,1) * t29;
t80 = rSges(4,2) * t47;
t31 = sin(t41);
t20 = pkin(5) * t31 + t87;
t79 = t20 * t44;
t75 = t45 * t31;
t32 = cos(t41);
t74 = t45 * t32;
t73 = t45 * t43;
t72 = t45 * t46;
t71 = t46 * t48;
t70 = t47 * t48;
t21 = pkin(5) * t32 + t36;
t34 = t44 * qJ(3);
t67 = t34 + t37;
t66 = t48 * pkin(1) + t45 * pkin(7);
t65 = qJ(3) * t47;
t64 = -pkin(2) - t82;
t63 = t44 * t87;
t62 = -pkin(2) - t69;
t61 = -pkin(2) - t68;
t60 = -pkin(1) - t34;
t59 = pkin(2) * t70 + t48 * t34 + t66;
t58 = g(1) * t64;
t57 = g(1) * t62;
t56 = g(1) * t61;
t54 = rSges(3,1) * t47 - rSges(3,2) * t44;
t52 = rSges(5,1) * t43 + rSges(5,2) * t46;
t15 = t44 * t71 - t73;
t17 = t43 * t48 + t44 * t72;
t51 = rSges(7,1) * t28 + rSges(7,2) * t29 + t20;
t50 = rSges(6,1) * t31 + rSges(6,2) * t32 + t87;
t23 = t45 * t65;
t25 = t48 * t65;
t49 = g(1) * t25 + g(2) * t23 + g(3) * t67;
t38 = t48 * pkin(7);
t30 = t36 + pkin(3);
t22 = t47 * t28 * rSges(7,2);
t19 = pkin(3) + t21;
t18 = -t44 * t73 + t71;
t16 = t43 * t78 + t72;
t12 = t32 * t48 - t44 * t75;
t11 = t31 * t48 + t44 * t74;
t10 = t31 * t78 + t74;
t9 = t32 * t78 - t75;
t1 = [-m(2) * (g(1) * (-t45 * rSges(2,1) - rSges(2,2) * t48) + g(2) * (rSges(2,1) * t48 - t45 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t48 + t38) + g(2) * (rSges(3,1) * t70 - rSges(3,2) * t78 + t66) + (g(1) * (-pkin(1) - t54) + g(2) * rSges(3,3)) * t45) - m(4) * (g(1) * (rSges(4,1) * t48 + t38) + g(2) * (-rSges(4,2) * t70 + rSges(4,3) * t78 + t59) + (g(1) * (-rSges(4,3) * t44 - t37 + t60 + t80) + g(2) * rSges(4,1)) * t45) - m(5) * (g(1) * (t18 * rSges(5,1) - t17 * rSges(5,2) + pkin(3) * t48 + t38) + g(2) * (t16 * rSges(5,1) + t15 * rSges(5,2) + t82 * t70 + t59) + (g(2) * pkin(3) + g(1) * t60 + t47 * t58) * t45) - m(6) * (g(1) * (t12 * rSges(6,1) - t11 * rSges(6,2) + t38) + g(2) * (t10 * rSges(6,1) + t9 * rSges(6,2) + t59) + (g(1) * t30 + g(2) * (t69 * t47 + t63)) * t48 + (g(1) * (t60 - t63) + g(2) * t30 + t47 * t57) * t45) - m(7) * (g(1) * (t8 * rSges(7,1) - t7 * rSges(7,2) + t38) + g(2) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t59) + (g(1) * t19 + g(2) * (t68 * t47 + t79)) * t48 + (g(1) * (t60 - t79) + g(2) * t19 + t47 * t56) * t45) -m(3) * (g(3) * t54 + t55 * (-rSges(3,1) * t44 - rSges(3,2) * t47)) - m(4) * (g(1) * (rSges(4,3) * t70 + t25) + g(2) * (rSges(4,3) * t45 * t47 + t23) + g(3) * (t67 - t80) + (g(3) * rSges(4,3) + t55 * (rSges(4,2) - pkin(2))) * t44) - m(5) * ((g(3) * t82 + t55 * t52) * t47 + (g(3) * t52 + t48 * t58 + t64 * t85) * t44 + t49) - m(6) * ((g(3) * t69 + t55 * t50) * t47 + (g(3) * t50 + t48 * t57 + t62 * t85) * t44 + t49) - m(7) * ((g(3) * t68 + t55 * t51) * t47 + (g(3) * t51 + t48 * t56 + t61 * t85) * t44 + t49) (-m(4) - m(5) + t90) * (t55 * t44 - t84) -m(5) * (g(1) * (rSges(5,1) * t15 - rSges(5,2) * t16) + g(2) * (rSges(5,1) * t17 + rSges(5,2) * t18)) - m(6) * (g(1) * (t9 * rSges(6,1) - t10 * rSges(6,2) + t15 * pkin(4)) + g(2) * (t11 * rSges(6,1) + t12 * rSges(6,2) + t17 * pkin(4))) - m(7) * (g(1) * (-t45 * t20 + t21 * t78 + t89) + g(2) * (t45 * t44 * t21 + t20 * t48 + t88) + g(3) * t22) + (-m(5) * (-rSges(5,1) * t46 + rSges(5,2) * t43) - m(6) * (-rSges(6,1) * t32 + rSges(6,2) * t31 - t36) - m(7) * (-t21 - t81)) * t84, t90 * (g(3) * t44 + t55 * t47) -m(7) * (g(1) * t89 + g(2) * t88 + g(3) * (-t47 * t81 + t22))];
taug  = t1(:);
