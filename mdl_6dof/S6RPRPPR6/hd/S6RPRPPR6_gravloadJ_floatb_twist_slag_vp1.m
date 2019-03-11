% Calculate Gravitation load on the joints for
% S6RPRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
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
% Datum: 2019-03-09 02:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPPR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR6_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:53:17
% EndTime: 2019-03-09 02:53:18
% DurationCPUTime: 0.64s
% Computational Cost: add. (274->112), mult. (324->144), div. (0->0), fcn. (288->10), ass. (0->56)
t21 = qJ(3) + pkin(9);
t16 = cos(t21);
t42 = rSges(6,3) + qJ(5);
t67 = t42 * t16;
t44 = rSges(7,3) + pkin(8) + qJ(5);
t29 = cos(qJ(1));
t56 = g(2) * t29;
t27 = sin(qJ(1));
t59 = g(1) * t27;
t38 = -t56 + t59;
t66 = g(1) * t29 + g(2) * t27;
t14 = sin(t21);
t23 = cos(pkin(10));
t11 = t23 * pkin(5) + pkin(4);
t20 = pkin(10) + qJ(6);
t13 = sin(t20);
t15 = cos(t20);
t31 = rSges(7,1) * t15 - rSges(7,2) * t13 + t11;
t65 = t44 * t14 + t31 * t16;
t22 = sin(pkin(10));
t32 = rSges(6,1) * t23 - rSges(6,2) * t22 + pkin(4);
t64 = t42 * t14 + t32 * t16;
t28 = cos(qJ(3));
t61 = pkin(3) * t28;
t9 = t27 * t61;
t63 = g(1) * t9;
t62 = -m(6) - m(7);
t60 = pkin(5) * t22;
t55 = g(3) * t14;
t26 = sin(qJ(3));
t54 = t26 * pkin(3);
t53 = rSges(4,3) + pkin(7);
t52 = t27 * t14;
t51 = t27 * t15;
t50 = t27 * t22;
t49 = t27 * t23;
t48 = t29 * t14;
t47 = t29 * t15;
t46 = t29 * t22;
t45 = t29 * t23;
t43 = t29 * pkin(1) + t27 * qJ(2);
t41 = -m(5) + t62;
t40 = t27 * t54 + t43;
t18 = t29 * qJ(2);
t24 = -qJ(4) - pkin(7);
t39 = t27 * t24 + t29 * t54 + t18;
t37 = t44 * t16;
t35 = t26 * rSges(4,1) + t28 * rSges(4,2);
t34 = rSges(5,1) * t16 - rSges(5,2) * t14;
t33 = t14 * rSges(5,1) + t16 * rSges(5,2);
t30 = t14 * t11 - t37;
t5 = -t27 * t13 + t14 * t47;
t4 = t13 * t48 + t51;
t3 = t29 * t13 + t14 * t51;
t2 = -t13 * t52 + t47;
t1 = [-m(2) * (g(1) * (-t27 * rSges(2,1) - t29 * rSges(2,2)) + g(2) * (t29 * rSges(2,1) - t27 * rSges(2,2))) - m(3) * (g(1) * (t29 * rSges(3,3) + t18 + (rSges(3,2) - pkin(1)) * t27) + g(2) * (-t29 * rSges(3,2) + t27 * rSges(3,3) + t43)) - m(4) * (g(1) * t18 + g(2) * t43 + (g(1) * t35 + g(2) * t53) * t29 + (g(1) * (-pkin(1) - t53) + g(2) * t35) * t27) - m(5) * (g(1) * t39 + g(2) * t40 + (g(1) * t33 + g(2) * (rSges(5,3) - t24)) * t29 + (g(1) * (-rSges(5,3) - pkin(1)) + g(2) * t33) * t27) - m(6) * (g(1) * (pkin(4) * t48 - t27 * pkin(1) + (t14 * t45 - t50) * rSges(6,1) + (-t14 * t46 - t49) * rSges(6,2) + t39) + g(2) * (-t29 * t24 + pkin(4) * t52 + (t14 * t49 + t46) * rSges(6,1) + (-t14 * t50 + t45) * rSges(6,2) + t40) - t66 * t67) - m(7) * (g(1) * (t5 * rSges(7,1) - t4 * rSges(7,2) + t39) + g(2) * (t3 * rSges(7,1) + t2 * rSges(7,2) + t40) + (g(1) * t30 + g(2) * (-t24 + t60)) * t29 + (g(1) * (-pkin(1) - t60) + g(2) * t30) * t27) (-m(3) - m(4) + t41) * t38, -m(4) * (-g(3) * t35 + t38 * (rSges(4,1) * t28 - rSges(4,2) * t26)) - m(5) * (g(1) * (t27 * t34 + t9) + g(3) * (-t33 - t54) + (-t34 - t61) * t56) - m(6) * (t63 + g(3) * (-t54 + t67) - t32 * t55 + t64 * t59 + (-t61 - t64) * t56) - m(7) * (t63 + g(3) * (t37 - t54) - t31 * t55 + t65 * t59 + (-t61 - t65) * t56) t41 * t66, t62 * (-t16 * t38 + t55) -m(7) * (g(1) * (t2 * rSges(7,1) - t3 * rSges(7,2)) + g(2) * (t4 * rSges(7,1) + t5 * rSges(7,2)) + g(3) * (-rSges(7,1) * t13 - rSges(7,2) * t15) * t16)];
taug  = t1(:);
