% Calculate Gravitation load on the joints for
% S5RRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR15_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR15_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR15_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR15_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:41:12
% EndTime: 2019-12-31 20:41:14
% DurationCPUTime: 0.56s
% Computational Cost: add. (212->112), mult. (367->158), div. (0->0), fcn. (341->8), ass. (0->56)
t32 = cos(qJ(2));
t65 = g(3) * t32;
t63 = rSges(5,3) + pkin(7);
t51 = rSges(6,3) + pkin(8) + pkin(7);
t33 = cos(qJ(1));
t30 = sin(qJ(1));
t66 = g(2) * t30;
t69 = g(1) * t33 + t66;
t28 = sin(qJ(4));
t68 = pkin(4) * t28;
t24 = t32 * pkin(2);
t29 = sin(qJ(2));
t62 = t29 * t33;
t27 = qJ(4) + qJ(5);
t20 = sin(t27);
t61 = t30 * t20;
t21 = cos(t27);
t60 = t30 * t21;
t59 = t30 * t28;
t31 = cos(qJ(4));
t58 = t30 * t31;
t57 = t32 * rSges(4,2);
t56 = t32 * t33;
t55 = t33 * t20;
t54 = t33 * t21;
t53 = t33 * t28;
t52 = t33 * t31;
t22 = t29 * qJ(3);
t50 = t22 + t24;
t49 = t33 * pkin(1) + t30 * pkin(6);
t48 = qJ(3) * t32;
t47 = -pkin(2) - t63;
t46 = t29 * t68;
t45 = -pkin(2) - t51;
t44 = -pkin(1) - t22;
t43 = pkin(2) * t56 + t33 * t22 + t49;
t42 = g(1) * t47;
t41 = g(1) * t45;
t40 = t32 * rSges(3,1) - t29 * rSges(3,2);
t38 = rSges(5,1) * t28 + rSges(5,2) * t31;
t10 = t29 * t52 - t59;
t12 = t29 * t58 + t53;
t37 = rSges(6,1) * t20 + rSges(6,2) * t21 + t68;
t15 = t30 * t48;
t17 = t33 * t48;
t36 = g(1) * t17 + g(2) * t15 + g(3) * t50;
t5 = t29 * t54 - t61;
t6 = t29 * t55 + t60;
t7 = t29 * t60 + t55;
t8 = -t29 * t61 + t54;
t35 = g(1) * (t5 * rSges(6,1) - t6 * rSges(6,2)) + g(2) * (t7 * rSges(6,1) + t8 * rSges(6,2)) + (-rSges(6,1) * t21 + rSges(6,2) * t20) * t65;
t25 = t33 * pkin(6);
t19 = t31 * pkin(4) + pkin(3);
t13 = -t29 * t59 + t52;
t11 = t29 * t53 + t58;
t1 = [-m(2) * (g(1) * (-t30 * rSges(2,1) - t33 * rSges(2,2)) + g(2) * (t33 * rSges(2,1) - t30 * rSges(2,2))) - m(3) * (g(1) * (t33 * rSges(3,3) + t25) + g(2) * (rSges(3,1) * t56 - rSges(3,2) * t62 + t49) + (g(1) * (-pkin(1) - t40) + g(2) * rSges(3,3)) * t30) - m(4) * (g(1) * (t33 * rSges(4,1) + t25) + g(2) * (-rSges(4,2) * t56 + rSges(4,3) * t62 + t43) + (g(1) * (-t29 * rSges(4,3) - t24 + t44 + t57) + g(2) * rSges(4,1)) * t30) - m(5) * (g(1) * (t13 * rSges(5,1) - t12 * rSges(5,2) + t33 * pkin(3) + t25) + g(2) * (t11 * rSges(5,1) + t10 * rSges(5,2) + t63 * t56 + t43) + (g(2) * pkin(3) + g(1) * t44 + t32 * t42) * t30) - m(6) * (g(1) * (t8 * rSges(6,1) - t7 * rSges(6,2) + t25) + g(2) * (t6 * rSges(6,1) + t5 * rSges(6,2) + t43) + (g(1) * t19 + g(2) * (t51 * t32 + t46)) * t33 + (g(1) * (t44 - t46) + g(2) * t19 + t32 * t41) * t30), -m(3) * (g(3) * t40 + t69 * (-rSges(3,1) * t29 - rSges(3,2) * t32)) - m(4) * (g(1) * (rSges(4,3) * t56 + t17) + g(2) * (t30 * t32 * rSges(4,3) + t15) + g(3) * (t50 - t57) + (g(3) * rSges(4,3) + t69 * (rSges(4,2) - pkin(2))) * t29) - m(5) * ((g(3) * t63 + t69 * t38) * t32 + (g(3) * t38 + t33 * t42 + t47 * t66) * t29 + t36) - m(6) * ((g(3) * t51 + t69 * t37) * t32 + (g(3) * t37 + t33 * t41 + t45 * t66) * t29 + t36), (-m(4) - m(5) - m(6)) * (t69 * t29 - t65), -m(5) * (g(1) * (t10 * rSges(5,1) - t11 * rSges(5,2)) + g(2) * (t12 * rSges(5,1) + t13 * rSges(5,2)) + (-rSges(5,1) * t31 + rSges(5,2) * t28) * t65) - m(6) * ((g(1) * t10 + g(2) * t12 - t31 * t65) * pkin(4) + t35), -m(6) * t35];
taug = t1(:);
