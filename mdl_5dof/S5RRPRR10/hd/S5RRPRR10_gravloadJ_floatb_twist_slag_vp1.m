% Calculate Gravitation load on the joints for
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR10_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR10_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR10_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR10_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:23:48
% EndTime: 2019-12-31 20:23:52
% DurationCPUTime: 1.03s
% Computational Cost: add. (519->146), mult. (1259->225), div. (0->0), fcn. (1563->12), ass. (0->65)
t44 = sin(pkin(10));
t49 = sin(qJ(2));
t53 = cos(qJ(2));
t70 = cos(pkin(10));
t34 = -t53 * t44 - t49 * t70;
t50 = sin(qJ(1));
t54 = cos(qJ(1));
t46 = cos(pkin(5));
t57 = -t49 * t44 + t53 * t70;
t55 = t57 * t46;
t17 = t54 * t34 - t50 * t55;
t73 = t54 * t49;
t75 = t50 * t53;
t31 = -t46 * t75 - t73;
t56 = t31 * pkin(2);
t92 = t17 * pkin(3) + t56;
t72 = t54 * t53;
t76 = t50 * t49;
t91 = t46 * t72 - t76;
t14 = t50 * t34 + t54 * t55;
t27 = t34 * t46;
t15 = -t54 * t27 + t50 * t57;
t48 = sin(qJ(4));
t52 = cos(qJ(4));
t45 = sin(pkin(5));
t80 = t45 * t54;
t4 = t15 * t52 - t48 * t80;
t47 = sin(qJ(5));
t51 = cos(qJ(5));
t90 = t14 * t51 + t4 * t47;
t89 = t14 * t47 - t4 * t51;
t88 = t52 * pkin(4);
t87 = t53 * pkin(2);
t85 = rSges(5,3) + pkin(8);
t84 = rSges(6,3) + pkin(9);
t81 = t45 * t50;
t79 = t47 * t52;
t43 = pkin(1) + t87;
t77 = t50 * t43;
t74 = t51 * t52;
t25 = t57 * t45;
t41 = t45 * t87;
t71 = t25 * pkin(3) + t41;
t67 = g(2) * t84;
t28 = t46 * t49 * pkin(2) + (-pkin(7) - qJ(3)) * t45;
t66 = rSges(4,3) * t45 - t28;
t65 = -t15 * t48 - t52 * t80;
t63 = t91 * pkin(2);
t16 = -t50 * t27 - t54 * t57;
t37 = t54 * t43;
t62 = -pkin(3) * t16 - t50 * t28 + t37;
t61 = rSges(5,1) * t52 - rSges(5,2) * t48;
t60 = t14 * pkin(3) + t63;
t59 = rSges(6,1) * t51 - rSges(6,2) * t47 + pkin(4);
t58 = -t15 * pkin(3) - t54 * t28 - t77;
t32 = -t46 * t76 + t72;
t30 = -t46 * t73 - t75;
t26 = t34 * t45;
t20 = -t26 * t52 + t46 * t48;
t19 = t26 * t48 + t46 * t52;
t8 = -t16 * t52 + t48 * t81;
t7 = -t16 * t48 - t52 * t81;
t2 = -t17 * t47 + t8 * t51;
t1 = -t17 * t51 - t8 * t47;
t3 = [-m(2) * (g(1) * (-t50 * rSges(2,1) - t54 * rSges(2,2)) + g(2) * (t54 * rSges(2,1) - t50 * rSges(2,2))) - m(3) * (g(1) * (t30 * rSges(3,1) - rSges(3,2) * t91 - t50 * pkin(1)) + g(2) * (t32 * rSges(3,1) + t31 * rSges(3,2) + t54 * pkin(1)) + (g(1) * t54 + g(2) * t50) * t45 * (rSges(3,3) + pkin(7))) - m(4) * (g(1) * (-t15 * rSges(4,1) - t14 * rSges(4,2) + t66 * t54 - t77) + g(2) * (-rSges(4,1) * t16 + t17 * rSges(4,2) + t66 * t50 + t37)) - m(5) * (g(1) * (-rSges(5,1) * t4 - rSges(5,2) * t65 + t85 * t14 + t58) + g(2) * (t8 * rSges(5,1) - t7 * rSges(5,2) - t85 * t17 + t62)) - m(6) * (g(1) * (t89 * rSges(6,1) + t90 * rSges(6,2) - t4 * pkin(4) + t14 * pkin(8) + t84 * t65 + t58) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t8 * pkin(4) - t17 * pkin(8) + t84 * t7 + t62)), -m(3) * (g(1) * (t31 * rSges(3,1) - t32 * rSges(3,2)) + g(2) * (rSges(3,1) * t91 + t30 * rSges(3,2)) + g(3) * (rSges(3,1) * t53 - rSges(3,2) * t49) * t45) - m(4) * (g(1) * (t17 * rSges(4,1) + t16 * rSges(4,2) + t56) + g(2) * (t14 * rSges(4,1) - rSges(4,2) * t15 + t63) + g(3) * (t25 * rSges(4,1) + t26 * rSges(4,2) + t41)) - m(5) * (g(1) * (-t85 * t16 + t61 * t17 + t92) + g(2) * (t61 * t14 + t15 * t85 + t60) + g(3) * (t61 * t25 - t85 * t26 + t71)) - m(6) * (g(1) * (t17 * t88 - t16 * pkin(8) + (-t16 * t47 + t17 * t74) * rSges(6,1) + (-t16 * t51 - t17 * t79) * rSges(6,2) + t92) + g(2) * (t14 * t88 + t15 * pkin(8) + (t14 * t74 + t15 * t47) * rSges(6,1) + (-t14 * t79 + t15 * t51) * rSges(6,2) + t60) + g(3) * (t25 * t88 - t26 * pkin(8) + (t25 * t74 - t26 * t47) * rSges(6,1) + (-t25 * t79 - t26 * t51) * rSges(6,2) + t71) + (t14 * t67 + (g(1) * t17 + g(3) * t25) * t84) * t48), (-m(4) - m(5) - m(6)) * (g(3) * t46 + (g(1) * t50 - g(2) * t54) * t45), -m(5) * (g(1) * (-t7 * rSges(5,1) - t8 * rSges(5,2)) + g(2) * (rSges(5,1) * t65 - t4 * rSges(5,2)) + g(3) * (t19 * rSges(5,1) - t20 * rSges(5,2))) - m(6) * (g(1) * (-t59 * t7 + t84 * t8) + t4 * t67 + g(2) * t59 * t65 + (t59 * t19 + t84 * t20) * g(3)), -m(6) * (g(1) * (t1 * rSges(6,1) - t2 * rSges(6,2)) + g(2) * (-t90 * rSges(6,1) + t89 * rSges(6,2)) + g(3) * ((-t20 * t47 - t25 * t51) * rSges(6,1) + (-t20 * t51 + t25 * t47) * rSges(6,2)))];
taug = t3(:);
