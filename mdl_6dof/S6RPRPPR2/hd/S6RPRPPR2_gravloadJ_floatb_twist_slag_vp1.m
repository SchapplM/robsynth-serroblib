% Calculate Gravitation load on the joints for
% S6RPRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-03-09 02:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPPR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:41:26
% EndTime: 2019-03-09 02:41:27
% DurationCPUTime: 0.46s
% Computational Cost: add. (343->103), mult. (287->131), div. (0->0), fcn. (247->10), ass. (0->55)
t21 = qJ(3) + pkin(10);
t15 = sin(t21);
t12 = t15 * qJ(5);
t17 = cos(t21);
t13 = t17 * pkin(4);
t63 = t12 + t13;
t22 = qJ(1) + pkin(9);
t18 = cos(t22);
t16 = sin(t22);
t58 = g(2) * t16;
t62 = g(1) * t18 + t58;
t61 = -m(6) - m(7);
t25 = sin(qJ(3));
t60 = pkin(3) * t25;
t57 = g(3) * t17;
t26 = sin(qJ(1));
t56 = t26 * pkin(1);
t55 = rSges(4,3) + pkin(7);
t54 = rSges(7,3) + pkin(8);
t23 = -qJ(4) - pkin(7);
t53 = pkin(5) - t23;
t24 = sin(qJ(6));
t52 = t16 * t24;
t27 = cos(qJ(6));
t51 = t16 * t27;
t50 = t18 * t24;
t49 = t18 * t27;
t48 = rSges(6,1) - t23;
t47 = rSges(5,3) - t23;
t28 = cos(qJ(3));
t19 = t28 * pkin(3);
t14 = t19 + pkin(2);
t29 = cos(qJ(1));
t20 = t29 * pkin(1);
t46 = t18 * t14 + t20;
t45 = qJ(5) * t17;
t44 = -m(5) + t61;
t43 = g(1) * t56;
t42 = -pkin(4) - t54;
t41 = t19 + t63;
t40 = t63 * t18 + t46;
t39 = -t14 - t12;
t38 = g(1) * t42;
t37 = t28 * rSges(4,1) - t25 * rSges(4,2);
t35 = t17 * rSges(5,1) - t15 * rSges(5,2);
t34 = rSges(7,1) * t24 + rSges(7,2) * t27;
t33 = t17 * rSges(6,2) - t15 * rSges(6,3);
t32 = pkin(2) + t37;
t9 = t18 * t45;
t7 = t16 * t45;
t5 = -t15 * t52 + t49;
t4 = t15 * t51 + t50;
t3 = t15 * t50 + t51;
t2 = t15 * t49 - t52;
t1 = [-m(2) * (g(1) * (-t26 * rSges(2,1) - t29 * rSges(2,2)) + g(2) * (t29 * rSges(2,1) - t26 * rSges(2,2))) - m(3) * (g(1) * (-t16 * rSges(3,1) - t18 * rSges(3,2) - t56) + g(2) * (t18 * rSges(3,1) - t16 * rSges(3,2) + t20)) - m(4) * (-t43 + g(2) * t20 + (g(1) * t55 + g(2) * t32) * t18 + (-g(1) * t32 + g(2) * t55) * t16) - m(5) * (-t43 + g(2) * t46 + (g(1) * t47 + g(2) * t35) * t18 + (g(1) * (-t14 - t35) + g(2) * t47) * t16) - m(6) * (-t43 + g(2) * t40 + (g(1) * t48 - g(2) * t33) * t18 + (g(1) * (t33 + t39 - t13) + g(2) * t48) * t16) - m(7) * (g(1) * (t5 * rSges(7,1) - t4 * rSges(7,2) - t56) + g(2) * (t3 * rSges(7,1) + t2 * rSges(7,2) + t40) + (g(2) * t54 * t17 + g(1) * t53) * t18 + (g(1) * t39 + g(2) * t53 + t17 * t38) * t16) (-m(3) - m(4) + t44) * g(3), -m(4) * (g(3) * t37 + t62 * (-rSges(4,1) * t25 - rSges(4,2) * t28)) - m(5) * (g(3) * (t19 + t35) + t62 * (-rSges(5,1) * t15 - rSges(5,2) * t17 - t60)) - m(6) * (g(1) * t9 + g(2) * t7 + g(3) * (-t33 + t41) + t62 * (rSges(6,3) * t17 - t60 + (rSges(6,2) - pkin(4)) * t15)) - m(7) * (g(1) * (-t18 * t60 + t9) + g(2) * (-t16 * t60 + t7) + g(3) * t41 + (g(3) * t54 + t62 * t34) * t17 + (g(3) * t34 + t18 * t38 + t42 * t58) * t15) t44 * (g(1) * t16 - g(2) * t18) t61 * (t62 * t15 - t57) -m(7) * (g(1) * (t2 * rSges(7,1) - t3 * rSges(7,2)) + g(2) * (t4 * rSges(7,1) + t5 * rSges(7,2)) + (-rSges(7,1) * t27 + rSges(7,2) * t24) * t57)];
taug  = t1(:);
