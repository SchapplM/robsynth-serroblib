% Calculate Gravitation load on the joints for
% S6RPRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-03-09 04:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPP5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:42:43
% EndTime: 2019-03-09 04:42:44
% DurationCPUTime: 0.70s
% Computational Cost: add. (415->140), mult. (527->179), div. (0->0), fcn. (521->8), ass. (0->54)
t26 = pkin(9) + qJ(3);
t24 = sin(t26);
t65 = g(3) * t24;
t63 = rSges(7,1) + pkin(5);
t31 = sin(qJ(1));
t33 = cos(qJ(1));
t71 = g(1) * t33 + g(2) * t31;
t47 = -rSges(7,3) - qJ(6);
t70 = -m(6) - m(7);
t29 = -pkin(7) - qJ(2);
t68 = g(2) * t29;
t32 = cos(qJ(4));
t66 = t32 * qJ(5) * t65;
t19 = t24 * pkin(8);
t25 = cos(t26);
t20 = t25 * pkin(3);
t64 = -rSges(6,1) - pkin(4);
t62 = t24 * t33;
t30 = sin(qJ(4));
t61 = t25 * t30;
t60 = t25 * t31;
t59 = t25 * t32;
t58 = t25 * t33;
t57 = t29 * t33;
t56 = t31 * t30;
t55 = t31 * t32;
t54 = t32 * t33;
t53 = t33 * t30;
t52 = rSges(4,3) - t29;
t51 = t19 + t20;
t50 = rSges(7,2) + qJ(5);
t49 = rSges(3,3) + qJ(2);
t48 = rSges(6,3) + qJ(5);
t46 = -pkin(4) - t63;
t28 = cos(pkin(9));
t21 = pkin(2) * t28 + pkin(1);
t12 = t33 * t21;
t45 = pkin(3) * t58 + pkin(8) * t62 + t12;
t44 = -t21 - t20;
t43 = pkin(4) * t59 + qJ(5) * t61 + t51;
t42 = rSges(4,1) * t25 - rSges(4,2) * t24;
t40 = t44 - t19;
t39 = rSges(3,1) * t28 - rSges(3,2) * sin(pkin(9)) + pkin(1);
t6 = t25 * t56 + t54;
t7 = t25 * t55 - t53;
t37 = -t7 * pkin(4) - t6 * qJ(5) - t57;
t8 = t25 * t53 - t55;
t9 = t25 * t54 + t56;
t36 = t9 * pkin(4) + t8 * qJ(5) + t45;
t17 = pkin(8) * t58;
t14 = pkin(8) * t60;
t4 = t8 * pkin(4);
t2 = t6 * pkin(4);
t1 = [-m(2) * (g(1) * (-t31 * rSges(2,1) - rSges(2,2) * t33) + g(2) * (rSges(2,1) * t33 - t31 * rSges(2,2))) - m(3) * ((g(1) * t49 + g(2) * t39) * t33 + (-g(1) * t39 + g(2) * t49) * t31) - m(4) * (g(2) * t12 + (g(1) * t52 + g(2) * t42) * t33 + (g(1) * (-t21 - t42) + g(2) * t52) * t31) - m(5) * (g(1) * (-t7 * rSges(5,1) + t6 * rSges(5,2) - t57) + g(2) * (t9 * rSges(5,1) - t8 * rSges(5,2) + rSges(5,3) * t62 + t45) + (g(1) * (-rSges(5,3) * t24 + t40) - t68) * t31) - m(6) * (g(1) * (-t7 * rSges(6,1) - t6 * rSges(6,3) + t37) + g(2) * (t9 * rSges(6,1) + rSges(6,2) * t62 + t8 * rSges(6,3) + t36) + (g(1) * (-rSges(6,2) * t24 + t40) - t68) * t31) - m(7) * (g(1) * (-t6 * rSges(7,2) - t63 * t7 + t37) + g(2) * (t8 * rSges(7,2) + t47 * t62 + t63 * t9 + t36) + (-t68 + (t44 + (-pkin(8) - t47) * t24) * g(1)) * t31) (-m(3) - m(4) - m(5) + t70) * (g(1) * t31 - g(2) * t33) -m(4) * (g(3) * t42 + t71 * (-rSges(4,1) * t24 - rSges(4,2) * t25)) - m(5) * (g(1) * (rSges(5,3) * t58 + t17) + g(2) * (rSges(5,3) * t60 + t14) + g(3) * (rSges(5,1) * t59 - rSges(5,2) * t61 + t51) + (g(3) * rSges(5,3) + t71 * (-rSges(5,1) * t32 + rSges(5,2) * t30 - pkin(3))) * t24) - m(6) * (g(1) * (rSges(6,2) * t58 + t17) + g(2) * (rSges(6,2) * t60 + t14) + g(3) * (rSges(6,1) * t59 + rSges(6,3) * t61 + t43) + (g(3) * rSges(6,2) + t71 * (-t48 * t30 + t64 * t32 - pkin(3))) * t24) - m(7) * (g(1) * t17 + g(2) * t14 + g(3) * t43 + (g(3) * (rSges(7,2) * t30 + t63 * t32) + t71 * t47) * t25 + (g(3) * t47 + t71 * (-t50 * t30 + t46 * t32 - pkin(3))) * t24) -m(5) * (g(1) * (-rSges(5,1) * t8 - rSges(5,2) * t9) + g(2) * (-rSges(5,1) * t6 - rSges(5,2) * t7)) - m(6) * (g(1) * (-rSges(6,1) * t8 + t48 * t9 - t4) + g(2) * (-rSges(6,1) * t6 + t48 * t7 - t2) + t66) - m(7) * (g(1) * (t50 * t9 - t63 * t8 - t4) + g(2) * (t50 * t7 - t63 * t6 - t2) + t66) + ((m(5) * rSges(5,2) - m(6) * rSges(6,3) - m(7) * rSges(7,2)) * t32 + (m(5) * rSges(5,1) - m(6) * t64 - m(7) * t46) * t30) * t65, t70 * (g(1) * t8 + g(2) * t6 + t30 * t65) -m(7) * (g(3) * t25 - t71 * t24)];
taug  = t1(:);
