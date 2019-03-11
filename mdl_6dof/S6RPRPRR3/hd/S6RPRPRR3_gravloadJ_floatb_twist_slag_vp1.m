% Calculate Gravitation load on the joints for
% S6RPRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 03:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:40:38
% EndTime: 2019-03-09 03:40:39
% DurationCPUTime: 0.60s
% Computational Cost: add. (457->113), mult. (397->157), div. (0->0), fcn. (367->12), ass. (0->56)
t35 = sin(qJ(3));
t37 = cos(qJ(3));
t32 = sin(pkin(11));
t33 = cos(pkin(11));
t47 = rSges(5,1) * t33 - rSges(5,2) * t32 + pkin(3);
t56 = rSges(5,3) + qJ(4);
t71 = t56 * t35 + t47 * t37;
t34 = -pkin(8) - qJ(4);
t57 = rSges(7,3) + pkin(9) - t34;
t73 = t57 * t35;
t58 = rSges(6,3) - t34;
t72 = t58 * t35;
t31 = qJ(1) + pkin(10);
t23 = sin(t31);
t25 = cos(t31);
t70 = g(1) * t25 + g(2) * t23;
t30 = pkin(11) + qJ(5);
t26 = qJ(6) + t30;
t19 = sin(t26);
t20 = cos(t26);
t60 = t23 * t37;
t5 = t19 * t60 + t20 * t25;
t6 = t19 * t25 - t20 * t60;
t69 = -t5 * rSges(7,1) + t6 * rSges(7,2);
t59 = t25 * t37;
t7 = -t19 * t59 + t20 * t23;
t8 = t19 * t23 + t20 * t59;
t68 = t7 * rSges(7,1) - t8 * rSges(7,2);
t67 = pkin(4) * t32;
t22 = sin(t30);
t66 = pkin(5) * t22;
t63 = g(3) * t35;
t36 = sin(qJ(1));
t62 = t36 * pkin(1);
t21 = t33 * pkin(4) + pkin(3);
t61 = rSges(4,2) * t35;
t55 = -m(5) - m(6) - m(7);
t38 = cos(qJ(1));
t28 = t38 * pkin(1);
t54 = t25 * pkin(2) + t23 * pkin(7) + t28;
t53 = t25 * pkin(7) - t62;
t51 = rSges(4,1) * t37 - t61;
t49 = rSges(5,1) * t32 + rSges(5,2) * t33;
t48 = -rSges(7,1) * t19 - rSges(7,2) * t20;
t24 = cos(t30);
t11 = -t22 * t59 + t23 * t24;
t9 = t22 * t60 + t24 * t25;
t46 = rSges(6,1) * t24 - rSges(6,2) * t22 + t21;
t14 = pkin(5) * t24 + t21;
t45 = rSges(7,1) * t20 - rSges(7,2) * t19 + t14;
t43 = t37 * t21 + t72;
t42 = t37 * t14 + t73;
t15 = t66 + t67;
t12 = t22 * t23 + t24 * t59;
t10 = t22 * t25 - t24 * t60;
t1 = [-m(2) * (g(1) * (-t36 * rSges(2,1) - rSges(2,2) * t38) + g(2) * (rSges(2,1) * t38 - t36 * rSges(2,2))) - m(3) * (g(1) * (-rSges(3,1) * t23 - rSges(3,2) * t25 - t62) + g(2) * (rSges(3,1) * t25 - rSges(3,2) * t23 + t28)) - m(4) * (g(1) * (rSges(4,3) * t25 + t53) + g(2) * (rSges(4,1) * t59 - t25 * t61 + t54) + (g(1) * (-pkin(2) - t51) + g(2) * rSges(4,3)) * t23) - m(5) * (g(1) * t53 + g(2) * t54 + (g(1) * t49 + t71 * g(2)) * t25 + (g(2) * t49 + (-pkin(2) - t71) * g(1)) * t23) - m(6) * (g(1) * (rSges(6,1) * t10 + rSges(6,2) * t9 + t53) + g(2) * (rSges(6,1) * t12 + rSges(6,2) * t11 + t54) + (g(1) * t67 + g(2) * t43) * t25 + (g(1) * (-pkin(2) - t43) + g(2) * t67) * t23) - m(7) * (g(1) * (rSges(7,1) * t6 + rSges(7,2) * t5 + t53) + g(2) * (rSges(7,1) * t8 + rSges(7,2) * t7 + t54) + (g(1) * t15 + g(2) * t42) * t25 + (g(1) * (-pkin(2) - t42) + g(2) * t15) * t23) (-m(3) - m(4) + t55) * g(3), -m(4) * (g(3) * t51 + t70 * (-rSges(4,1) * t35 - rSges(4,2) * t37)) - m(5) * (g(3) * t71 + t70 * (-t47 * t35 + t56 * t37)) - m(6) * (g(3) * (t46 * t37 + t72) + t70 * (-t46 * t35 + t58 * t37)) - m(7) * (g(3) * (t45 * t37 + t73) + t70 * (-t45 * t35 + t57 * t37)) t55 * (-g(3) * t37 + t70 * t35) -m(6) * (g(1) * (rSges(6,1) * t11 - rSges(6,2) * t12) + g(2) * (-rSges(6,1) * t9 + rSges(6,2) * t10)) - m(7) * (g(1) * (t11 * pkin(5) + t68) + g(2) * (-t9 * pkin(5) + t69)) + (-m(6) * (-rSges(6,1) * t22 - rSges(6,2) * t24) - m(7) * (t48 - t66)) * t63, -m(7) * (g(1) * t68 + g(2) * t69 + t48 * t63)];
taug  = t1(:);
