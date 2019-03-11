% Calculate Gravitation load on the joints for
% S6RRPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
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
% Datum: 2019-03-09 10:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPP5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPP5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:03:36
% EndTime: 2019-03-09 10:03:38
% DurationCPUTime: 0.72s
% Computational Cost: add. (281->141), mult. (593->184), div. (0->0), fcn. (584->6), ass. (0->53)
t66 = rSges(7,1) + pkin(5);
t31 = sin(qJ(1));
t34 = cos(qJ(1));
t44 = g(1) * t34 + g(2) * t31;
t73 = -m(6) - m(7);
t72 = -pkin(2) - pkin(8);
t29 = sin(qJ(4));
t71 = pkin(4) * t29;
t70 = g(1) * t31;
t33 = cos(qJ(2));
t67 = g(3) * t33;
t25 = t33 * pkin(2);
t64 = rSges(4,2) * t33;
t63 = t29 * t31;
t30 = sin(qJ(2));
t62 = t30 * t34;
t32 = cos(qJ(4));
t61 = t31 * t32;
t60 = t32 * t34;
t59 = t33 * t34;
t21 = t30 * qJ(3);
t58 = t21 + t25;
t57 = t34 * pkin(1) + t31 * pkin(7);
t26 = t34 * pkin(7);
t56 = t34 * pkin(3) + t26;
t55 = qJ(3) * t33;
t54 = rSges(7,2) + qJ(5);
t53 = rSges(6,3) + qJ(5);
t52 = -rSges(7,3) - qJ(6);
t51 = -rSges(6,2) + t72;
t50 = -rSges(5,3) + t72;
t49 = t33 * t71;
t48 = t33 * pkin(8) + t58;
t47 = -pkin(1) - t21;
t46 = -t52 + t72;
t45 = pkin(2) * t59 + t34 * t21 + t57;
t43 = rSges(3,1) * t33 - rSges(3,2) * t30;
t41 = rSges(5,1) * t29 + rSges(5,2) * t32;
t10 = -t30 * t63 + t60;
t9 = t29 * t34 + t30 * t61;
t40 = t10 * pkin(4) + qJ(5) * t9 + t56;
t39 = t31 * pkin(3) + pkin(8) * t59 + t45;
t8 = t29 * t62 + t61;
t38 = t8 * pkin(4) + t39;
t37 = rSges(6,1) * t29 - t32 * t53;
t36 = t29 * t66 - t32 * t54;
t15 = t31 * t55;
t17 = t34 * t55;
t35 = g(1) * (t34 * t49 + t17) + g(2) * (t31 * t49 + t15) + g(3) * (t30 * t71 + t48);
t7 = -t30 * t60 + t63;
t5 = t9 * pkin(4);
t3 = t7 * pkin(4);
t1 = [-m(2) * (g(1) * (-t31 * rSges(2,1) - rSges(2,2) * t34) + g(2) * (rSges(2,1) * t34 - t31 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t34 + t26) + g(2) * (rSges(3,1) * t59 - rSges(3,2) * t62 + t57) + (g(1) * (-pkin(1) - t43) + g(2) * rSges(3,3)) * t31) - m(4) * (g(1) * (rSges(4,1) * t34 + t26) + g(2) * (-rSges(4,2) * t59 + rSges(4,3) * t62 + t45) + (g(1) * (-rSges(4,3) * t30 - t25 + t47 + t64) + g(2) * rSges(4,1)) * t31) - m(5) * (g(1) * (rSges(5,1) * t10 - rSges(5,2) * t9 + t56) + g(2) * (t8 * rSges(5,1) - t7 * rSges(5,2) + rSges(5,3) * t59 + t39) + (t33 * t50 + t47) * t70) - m(6) * (g(1) * (rSges(6,1) * t10 + rSges(6,3) * t9 + t40) + g(2) * (t8 * rSges(6,1) + rSges(6,2) * t59 + t53 * t7 + t38) + (t33 * t51 + t47) * t70) - m(7) * (g(1) * (rSges(7,2) * t9 + t66 * t10 + t40) + g(2) * (t52 * t59 + t54 * t7 + t66 * t8 + t38) + (t33 * t46 + t47) * t70) -m(3) * (g(3) * t43 + t44 * (-rSges(3,1) * t30 - rSges(3,2) * t33)) - m(4) * (g(1) * (rSges(4,3) * t59 + t17) + g(2) * (rSges(4,3) * t31 * t33 + t15) + g(3) * (t58 - t64) + (g(3) * rSges(4,3) + t44 * (rSges(4,2) - pkin(2))) * t30) - m(5) * (g(1) * t17 + g(2) * t15 + g(3) * t48 + (g(3) * rSges(5,3) + t44 * t41) * t33 + (g(3) * t41 + t44 * t50) * t30) - m(6) * ((g(3) * rSges(6,2) + t44 * t37) * t33 + (g(3) * t37 + t44 * t51) * t30 + t35) - m(7) * ((g(3) * t52 + t44 * t36) * t33 + (g(3) * t36 + t44 * t46) * t30 + t35) (-m(4) - m(5) + t73) * (t30 * t44 - t67) -m(5) * (g(1) * (-rSges(5,1) * t7 - rSges(5,2) * t8) + g(2) * (rSges(5,1) * t9 + rSges(5,2) * t10)) - m(6) * (g(1) * (-rSges(6,1) * t7 + t53 * t8 - t3) + g(2) * (rSges(6,1) * t9 - t10 * t53 + t5)) - m(7) * (g(1) * (t54 * t8 - t66 * t7 - t3) + g(2) * (-t10 * t54 + t66 * t9 + t5)) + ((-m(5) * rSges(5,2) + m(6) * t53 + m(7) * t54) * t29 + (m(5) * rSges(5,1) - m(6) * (-rSges(6,1) - pkin(4)) - m(7) * (-pkin(4) - t66)) * t32) * t67, t73 * (g(1) * t7 - g(2) * t9 + t32 * t67) -m(7) * (-g(3) * t30 - t33 * t44)];
taug  = t1(:);
