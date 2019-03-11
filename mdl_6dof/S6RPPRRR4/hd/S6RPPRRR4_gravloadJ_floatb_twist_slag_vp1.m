% Calculate Gravitation load on the joints for
% S6RPPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 02:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:25:29
% EndTime: 2019-03-09 02:25:30
% DurationCPUTime: 0.52s
% Computational Cost: add. (316->96), mult. (540->142), div. (0->0), fcn. (613->10), ass. (0->46)
t28 = sin(qJ(4));
t61 = g(3) * t28;
t26 = qJ(5) + qJ(6);
t20 = sin(t26);
t21 = cos(t26);
t66 = t20 * rSges(7,1) + t21 * rSges(7,2);
t59 = rSges(6,3) + pkin(8);
t50 = rSges(7,3) + pkin(9) + pkin(8);
t29 = cos(qJ(5));
t19 = t29 * pkin(5) + pkin(4);
t37 = rSges(7,1) * t21 - rSges(7,2) * t20 + t19;
t27 = sin(qJ(5));
t39 = t29 * rSges(6,1) - t27 * rSges(6,2) + pkin(4);
t65 = m(5) * rSges(5,1) + m(6) * t39 + m(7) * t37;
t63 = g(1) * t28;
t30 = cos(qJ(4));
t62 = g(1) * t30;
t60 = rSges(5,3) + pkin(7);
t58 = cos(qJ(1));
t57 = sin(qJ(1));
t55 = t20 * t30;
t53 = t21 * t30;
t52 = t27 * t30;
t51 = t29 * t30;
t49 = t58 * pkin(1) + t57 * qJ(2);
t48 = cos(pkin(10));
t47 = sin(pkin(10));
t46 = t58 * pkin(2) + t49;
t45 = m(4) + m(5) + m(6) + m(7);
t44 = pkin(5) * t27 + pkin(7);
t11 = -t57 * t47 - t58 * t48;
t43 = -t11 * pkin(3) + t46;
t42 = -t57 * pkin(1) + t58 * qJ(2);
t41 = t30 * rSges(5,1) - t28 * rSges(5,2);
t40 = rSges(6,1) * t27 + rSges(6,2) * t29;
t12 = t58 * t47 - t57 * t48;
t38 = -t11 * t29 + t12 * t52;
t7 = t11 * t52 + t12 * t29;
t36 = -t57 * pkin(2) + t42;
t35 = -m(5) * rSges(5,2) + m(6) * t59 + m(7) * t50;
t34 = g(1) * (t12 * pkin(3) + t36);
t5 = t11 * t55 + t12 * t21;
t6 = -t11 * t53 + t12 * t20;
t33 = g(1) * (t5 * rSges(7,1) - t6 * rSges(7,2)) + g(2) * ((-t11 * t21 + t12 * t55) * rSges(7,1) + (t11 * t20 + t12 * t53) * rSges(7,2)) + t66 * t61;
t8 = -t11 * t51 + t12 * t27;
t1 = [-m(2) * (g(1) * (-t57 * rSges(2,1) - t58 * rSges(2,2)) + g(2) * (t58 * rSges(2,1) - t57 * rSges(2,2))) - m(3) * (g(1) * (-t57 * rSges(3,1) + t58 * rSges(3,3) + t42) + g(2) * (t58 * rSges(3,1) + t57 * rSges(3,3) + t49)) - m(4) * (g(1) * (t12 * rSges(4,1) - t11 * rSges(4,2) + t36) + g(2) * (-t11 * rSges(4,1) - t12 * rSges(4,2) + t46)) - m(5) * (t34 + g(2) * t43 + (g(1) * t41 + g(2) * t60) * t12 + (g(1) * t60 - g(2) * t41) * t11) - m(6) * (t34 + g(2) * (t8 * rSges(6,1) + t7 * rSges(6,2) + t43) + (g(2) * pkin(7) + t39 * t62 + t59 * t63) * t12 + (g(1) * (pkin(7) + t40) + g(2) * (-t30 * pkin(4) - t59 * t28)) * t11) - m(7) * (t34 + g(2) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t43) + (g(2) * t44 + t37 * t62 + t50 * t63) * t12 + (g(1) * (t44 + t66) + g(2) * (-t30 * t19 - t50 * t28)) * t11) (-m(3) - t45) * (g(1) * t57 - g(2) * t58) t45 * g(3) (t35 * t28 + t65 * t30) * g(3) + (g(1) * t11 + g(2) * t12) * (-t65 * t28 + t35 * t30) -m(6) * (g(1) * (t7 * rSges(6,1) - t8 * rSges(6,2)) + g(2) * (t38 * rSges(6,1) + (t11 * t27 + t12 * t51) * rSges(6,2)) + t40 * t61) - m(7) * ((g(1) * t7 + g(2) * t38 + t27 * t61) * pkin(5) + t33) -m(7) * t33];
taug  = t1(:);
