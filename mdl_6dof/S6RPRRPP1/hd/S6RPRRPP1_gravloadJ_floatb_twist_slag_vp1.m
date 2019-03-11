% Calculate Gravitation load on the joints for
% S6RPRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-03-09 04:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:27:54
% EndTime: 2019-03-09 04:27:56
% DurationCPUTime: 0.73s
% Computational Cost: add. (475->126), mult. (459->168), div. (0->0), fcn. (439->10), ass. (0->59)
t34 = cos(qJ(3));
t31 = sin(qJ(3));
t65 = rSges(5,3) + pkin(8);
t51 = t65 * t31;
t77 = t34 * pkin(3) + t51;
t28 = qJ(1) + pkin(9);
t23 = sin(t28);
t25 = cos(t28);
t33 = cos(qJ(4));
t30 = sin(qJ(4));
t58 = t30 * t34;
t8 = t23 * t33 - t25 * t58;
t66 = rSges(7,1) + pkin(5);
t70 = g(2) * t23;
t76 = g(1) * t25 + t70;
t54 = rSges(7,3) + qJ(6);
t27 = qJ(4) + pkin(10);
t22 = sin(t27);
t24 = cos(t27);
t75 = t54 * t22 + t66 * t24;
t74 = -m(6) - m(7);
t73 = pkin(4) * t30;
t72 = g(1) * t23;
t69 = g(3) * t31;
t32 = sin(qJ(1));
t68 = t32 * pkin(1);
t64 = rSges(4,2) * t31;
t63 = t23 * t30;
t61 = t23 * t34;
t60 = t25 * t30;
t59 = t25 * t34;
t57 = t33 * t34;
t21 = pkin(4) * t33 + pkin(3);
t17 = t34 * t21;
t29 = -qJ(5) - pkin(8);
t56 = rSges(7,2) - t29;
t55 = rSges(6,3) - t29;
t35 = cos(qJ(1));
t26 = t35 * pkin(1);
t52 = t25 * pkin(2) + t23 * pkin(7) + t26;
t50 = t25 * pkin(7) - t68;
t49 = -pkin(2) - t17;
t48 = t25 * t56;
t47 = t55 * t31;
t46 = pkin(4) * t63 + t25 * t17 + t52;
t45 = rSges(4,1) * t34 - t64;
t43 = rSges(6,1) * t24 - rSges(6,2) * t22;
t42 = t8 * pkin(4);
t41 = t23 * t31 * t29 + pkin(4) * t60 + t50;
t40 = rSges(5,1) * t33 - rSges(5,2) * t30 + pkin(3);
t6 = t23 * t58 + t25 * t33;
t39 = t6 * pkin(4);
t9 = t25 * t57 + t63;
t7 = -t23 * t57 + t60;
t4 = t22 * t23 + t24 * t59;
t3 = t22 * t59 - t23 * t24;
t2 = -t25 * t22 + t24 * t61;
t1 = t22 * t61 + t24 * t25;
t5 = [-m(2) * (g(1) * (-t32 * rSges(2,1) - t35 * rSges(2,2)) + g(2) * (t35 * rSges(2,1) - t32 * rSges(2,2))) - m(3) * (g(1) * (-rSges(3,1) * t23 - t25 * rSges(3,2) - t68) + g(2) * (t25 * rSges(3,1) - rSges(3,2) * t23 + t26)) - m(4) * (g(1) * (t25 * rSges(4,3) + t50) + g(2) * (rSges(4,1) * t59 - t25 * t64 + t52) + (g(1) * (-pkin(2) - t45) + g(2) * rSges(4,3)) * t23) - m(5) * (g(1) * (t7 * rSges(5,1) + t6 * rSges(5,2) + t50) + (-pkin(2) - t77) * t72 + (t9 * rSges(5,1) + t8 * rSges(5,2) + t77 * t25 + t52) * g(2)) - m(6) * (g(1) * (-t2 * rSges(6,1) + t1 * rSges(6,2) + t41) + g(2) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t25 * t47 + t46) + (-t31 * rSges(6,3) + t49) * t72) - m(7) * (g(1) * (-t54 * t1 - t66 * t2 + t41) + g(2) * (t54 * t3 + t31 * t48 + t66 * t4 + t46) + (-rSges(7,2) * t31 + t49) * t72) (-m(3) - m(4) - m(5) + t74) * g(3), -m(4) * (g(3) * t45 + t76 * (-rSges(4,1) * t31 - rSges(4,2) * t34)) - m(5) * (g(3) * (t40 * t34 + t51) + t76 * (-t40 * t31 + t65 * t34)) - m(6) * (g(3) * (t43 * t34 + t17 + t47) + t76 * (t55 * t34 + (-t21 - t43) * t31)) - m(7) * (g(3) * t17 + (g(1) * t48 + g(3) * t75 + t56 * t70) * t34 + (g(3) * t56 + t76 * (-t21 - t75)) * t31) -m(5) * (g(1) * (rSges(5,1) * t8 - rSges(5,2) * t9) + g(2) * (-rSges(5,1) * t6 + rSges(5,2) * t7)) - m(6) * (g(1) * (-rSges(6,1) * t3 - rSges(6,2) * t4 + t42) + g(2) * (-rSges(6,1) * t1 - rSges(6,2) * t2 - t39)) - m(7) * (g(1) * (-t66 * t3 + t54 * t4 + t42) + g(2) * (-t66 * t1 + t54 * t2 - t39)) + (-m(5) * (-rSges(5,1) * t30 - rSges(5,2) * t33) - m(6) * (-rSges(6,1) * t22 - rSges(6,2) * t24 - t73) - m(7) * (-t66 * t22 + t54 * t24 - t73)) * t69, t74 * (-g(3) * t34 + t76 * t31) -m(7) * (g(1) * t3 + g(2) * t1 + t22 * t69)];
taug  = t5(:);
