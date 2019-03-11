% Calculate Gravitation load on the joints for
% S6RPPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 02:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:22:55
% EndTime: 2019-03-09 02:22:57
% DurationCPUTime: 0.54s
% Computational Cost: add. (329->96), mult. (322->134), div. (0->0), fcn. (291->10), ass. (0->46)
t24 = qJ(1) + pkin(10);
t19 = sin(t24);
t20 = cos(t24);
t67 = -g(1) * t19 + g(2) * t20;
t27 = sin(qJ(4));
t30 = cos(qJ(4));
t50 = rSges(6,3) + pkin(8);
t65 = t27 * pkin(4) - t30 * t50;
t43 = rSges(7,3) + pkin(9) + pkin(8);
t63 = rSges(5,1) * t27 + rSges(5,2) * t30;
t29 = cos(qJ(5));
t18 = pkin(5) * t29 + pkin(4);
t25 = qJ(5) + qJ(6);
t21 = sin(t25);
t22 = cos(t25);
t26 = sin(qJ(5));
t33 = m(6) * (rSges(6,1) * t29 - rSges(6,2) * t26 + pkin(4)) + m(7) * (rSges(7,1) * t22 - rSges(7,2) * t21 + t18) + m(5) * rSges(5,1);
t62 = -m(5) * rSges(5,2) + m(6) * t50 + m(7) * t43;
t61 = -pkin(2) - pkin(7);
t49 = t21 * t27;
t5 = -t19 * t49 + t20 * t22;
t48 = t22 * t27;
t6 = t19 * t48 + t20 * t21;
t60 = rSges(7,1) * t5 - rSges(7,2) * t6;
t7 = t19 * t22 + t20 * t49;
t8 = -t19 * t21 + t20 * t48;
t59 = rSges(7,1) * t7 + rSges(7,2) * t8;
t56 = pkin(5) * t26;
t53 = g(3) * t30;
t28 = sin(qJ(1));
t51 = t28 * pkin(1);
t47 = t26 * t27;
t45 = t27 * t29;
t31 = cos(qJ(1));
t23 = t31 * pkin(1);
t42 = pkin(2) * t20 + qJ(3) * t19 + t23;
t41 = -m(4) - m(5) - m(6) - m(7);
t40 = qJ(3) * t20 - t51;
t39 = pkin(7) * t20 + t42;
t38 = -rSges(7,1) * t21 - rSges(7,2) * t22;
t11 = t19 * t29 + t20 * t47;
t9 = -t19 * t47 + t20 * t29;
t35 = t27 * t18 - t30 * t43;
t12 = -t19 * t26 + t20 * t45;
t10 = t19 * t45 + t20 * t26;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t28 - rSges(2,2) * t31) + g(2) * (rSges(2,1) * t31 - rSges(2,2) * t28)) - m(3) * (g(1) * (-rSges(3,1) * t19 - rSges(3,2) * t20 - t51) + g(2) * (rSges(3,1) * t20 - rSges(3,2) * t19 + t23)) - m(4) * (g(1) * (t20 * rSges(4,3) + (rSges(4,2) - pkin(2)) * t19 + t40) + g(2) * (-rSges(4,2) * t20 + rSges(4,3) * t19 + t42)) - m(5) * (g(1) * (t63 * t20 + t40) + g(2) * (t20 * rSges(5,3) + t39) + (g(1) * (-rSges(5,3) + t61) + g(2) * t63) * t19) - m(6) * ((t10 * rSges(6,1) + t9 * rSges(6,2) + t19 * t65 + t39) * g(2) + (t12 * rSges(6,1) - t11 * rSges(6,2) + t61 * t19 + t20 * t65 + t40) * g(1)) - m(7) * (g(1) * (t8 * rSges(7,1) - t7 * rSges(7,2) + t40) + g(2) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t39) + (g(1) * t35 + g(2) * t56) * t20 + (g(1) * (-t56 + t61) + g(2) * t35) * t19) (-m(3) + t41) * g(3), -t41 * t67 (t27 * t33 - t30 * t62) * g(3) + t67 * (t27 * t62 + t33 * t30) -m(6) * (g(1) * (rSges(6,1) * t9 - rSges(6,2) * t10) + g(2) * (rSges(6,1) * t11 + rSges(6,2) * t12)) - m(7) * (g(1) * (pkin(5) * t9 + t60) + g(2) * (pkin(5) * t11 + t59)) + (-m(6) * (-rSges(6,1) * t26 - rSges(6,2) * t29) - m(7) * (t38 - t56)) * t53, -m(7) * (g(1) * t60 + g(2) * t59 + t38 * t53)];
taug  = t1(:);
