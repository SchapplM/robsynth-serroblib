% Calculate Gravitation load on the joints for
% S6PRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
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
% Datum: 2019-03-08 19:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRPR7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR7_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:51:13
% EndTime: 2019-03-08 19:51:15
% DurationCPUTime: 0.64s
% Computational Cost: add. (319->114), mult. (770->164), div. (0->0), fcn. (890->10), ass. (0->53)
t71 = rSges(7,3) + pkin(9);
t29 = sin(pkin(10));
t33 = sin(qJ(2));
t36 = cos(qJ(2));
t54 = cos(pkin(10));
t55 = cos(pkin(6));
t45 = t55 * t54;
t17 = t29 * t36 + t33 * t45;
t49 = t29 * t55;
t19 = -t33 * t49 + t54 * t36;
t70 = g(1) * t19 + g(2) * t17;
t16 = t29 * t33 - t36 * t45;
t18 = t54 * t33 + t36 * t49;
t69 = g(1) * t18 + g(2) * t16;
t68 = -m(6) - m(7);
t32 = sin(qJ(4));
t67 = pkin(4) * t32;
t30 = sin(pkin(6));
t62 = g(3) * t30;
t61 = t29 * t30;
t60 = t30 * t33;
t59 = t30 * t36;
t58 = pkin(2) * t59 + qJ(3) * t60;
t57 = rSges(4,3) + qJ(3);
t56 = rSges(6,3) + qJ(5);
t53 = pkin(8) * t59 + t58;
t52 = -m(4) - m(5) + t68;
t13 = t16 * pkin(2);
t51 = -t16 * pkin(8) - t13;
t14 = t18 * pkin(2);
t50 = -t18 * pkin(8) - t14;
t48 = t30 * t54;
t35 = cos(qJ(4));
t47 = rSges(5,1) * t32 + rSges(5,2) * t35;
t46 = g(3) * (t60 * t67 + t53);
t31 = sin(qJ(6));
t34 = cos(qJ(6));
t44 = rSges(7,1) * t34 - rSges(7,2) * t31 + pkin(5);
t42 = rSges(7,1) * t31 + rSges(7,2) * t34 + qJ(5);
t40 = -rSges(6,2) * t32 - t56 * t35;
t38 = t71 * t32 - t42 * t35;
t21 = -t32 * t59 + t55 * t35;
t20 = t55 * t32 + t35 * t59;
t15 = t20 * pkin(4);
t10 = t19 * t67;
t9 = t17 * t67;
t8 = -t16 * t32 + t35 * t48;
t7 = t16 * t35 + t32 * t48;
t6 = t18 * t32 + t35 * t61;
t5 = -t18 * t35 + t32 * t61;
t4 = t7 * pkin(4);
t3 = t5 * pkin(4);
t1 = [(-m(2) - m(3) + t52) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t18 - rSges(3,2) * t19) + g(2) * (-rSges(3,1) * t16 - rSges(3,2) * t17) + (rSges(3,1) * t36 - rSges(3,2) * t33) * t62) - m(4) * (g(1) * (rSges(4,2) * t18 + t57 * t19 - t14) + g(2) * (rSges(4,2) * t16 + t57 * t17 - t13) + g(3) * ((-rSges(4,2) * t36 + rSges(4,3) * t33) * t30 + t58)) - m(5) * (g(1) * (-t18 * rSges(5,3) + t50) + g(2) * (-rSges(5,3) * t16 + t51) + g(3) * t53 + (rSges(5,3) * t36 + t47 * t33) * t62 + t70 * (qJ(3) + t47)) - m(6) * (g(1) * (-rSges(6,1) * t18 + t10 + t50) + g(2) * (-rSges(6,1) * t16 + t51 + t9) + t46 + (rSges(6,1) * t36 + t40 * t33) * t62 + t70 * (qJ(3) + t40)) - m(7) * (g(1) * (t10 - t14) + g(2) * (-t13 + t9) + t46 + (t38 * t33 + t44 * t36) * t62 + t69 * (-pkin(8) - t44) + t70 * (qJ(3) + t38)) t52 * (-g(3) * t59 + t69) -m(5) * (g(1) * (-rSges(5,1) * t5 - rSges(5,2) * t6) + g(2) * (rSges(5,1) * t7 + rSges(5,2) * t8) + g(3) * (-rSges(5,1) * t20 - rSges(5,2) * t21)) - m(6) * (g(1) * (rSges(6,2) * t5 + t56 * t6 - t3) + g(2) * (-rSges(6,2) * t7 - t56 * t8 + t4) + g(3) * (rSges(6,2) * t20 + t56 * t21 - t15)) + (-g(1) * (-t71 * t5 - t3) - g(2) * (t71 * t7 + t4) - g(3) * (-t71 * t20 - t15) - (g(1) * t6 - g(2) * t8 + g(3) * t21) * t42) * m(7), t68 * (g(1) * t5 - g(2) * t7 + g(3) * t20) -m(7) * (g(1) * ((-t19 * t31 + t34 * t5) * rSges(7,1) + (-t19 * t34 - t31 * t5) * rSges(7,2)) + g(2) * ((-t17 * t31 - t34 * t7) * rSges(7,1) + (-t17 * t34 + t31 * t7) * rSges(7,2)) + g(3) * ((t20 * t34 - t31 * t60) * rSges(7,1) + (-t20 * t31 - t34 * t60) * rSges(7,2)))];
taug  = t1(:);
