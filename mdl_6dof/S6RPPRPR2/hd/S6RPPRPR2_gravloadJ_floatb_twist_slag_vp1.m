% Calculate Gravitation load on the joints for
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
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
% Datum: 2019-03-09 01:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRPR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:41:34
% EndTime: 2019-03-09 01:41:35
% DurationCPUTime: 0.42s
% Computational Cost: add. (322->90), mult. (261->111), div. (0->0), fcn. (224->10), ass. (0->47)
t20 = pkin(10) + qJ(4);
t15 = sin(t20);
t12 = t15 * qJ(5);
t17 = cos(t20);
t13 = t17 * pkin(4);
t43 = t12 + t13;
t21 = qJ(1) + pkin(9);
t16 = sin(t21);
t18 = cos(t21);
t58 = g(1) * t18 + g(2) * t16;
t57 = -m(6) - m(7);
t54 = g(3) * t17;
t26 = sin(qJ(1));
t53 = t26 * pkin(1);
t52 = rSges(7,3) + pkin(8);
t24 = -pkin(7) - qJ(3);
t51 = pkin(5) - t24;
t25 = sin(qJ(6));
t50 = t16 * t25;
t27 = cos(qJ(6));
t49 = t16 * t27;
t48 = t18 * t25;
t47 = t18 * t27;
t46 = rSges(6,1) - t24;
t45 = rSges(5,3) - t24;
t23 = cos(pkin(10));
t14 = pkin(3) * t23 + pkin(2);
t28 = cos(qJ(1));
t19 = t28 * pkin(1);
t44 = t14 * t18 + t19;
t41 = rSges(4,3) + qJ(3);
t40 = g(1) * t53;
t39 = -pkin(4) - t52;
t38 = -m(4) - m(5) + t57;
t37 = t52 * t17;
t36 = t18 * t43 + t44;
t35 = -t14 - t12;
t34 = t58 * qJ(5) * t17;
t33 = rSges(5,1) * t17 - rSges(5,2) * t15;
t32 = rSges(7,1) * t25 + rSges(7,2) * t27;
t31 = -t17 * rSges(6,2) + t15 * rSges(6,3);
t30 = rSges(4,1) * t23 - rSges(4,2) * sin(pkin(10)) + pkin(2);
t5 = -t15 * t50 + t47;
t4 = t15 * t49 + t48;
t3 = t15 * t48 + t49;
t2 = t15 * t47 - t50;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t26 - rSges(2,2) * t28) + g(2) * (rSges(2,1) * t28 - rSges(2,2) * t26)) - m(3) * (g(1) * (-rSges(3,1) * t16 - rSges(3,2) * t18 - t53) + g(2) * (rSges(3,1) * t18 - rSges(3,2) * t16 + t19)) - m(4) * (-t40 + g(2) * t19 + (g(1) * t41 + g(2) * t30) * t18 + (-g(1) * t30 + g(2) * t41) * t16) - m(5) * (-t40 + g(2) * t44 + (g(1) * t45 + g(2) * t33) * t18 + (g(1) * (-t14 - t33) + g(2) * t45) * t16) - m(6) * (-t40 + g(2) * t36 + (g(1) * t46 + g(2) * t31) * t18 + (g(1) * (-t31 + t35 - t13) + g(2) * t46) * t16) - m(7) * ((t3 * rSges(7,1) + t2 * rSges(7,2) + t16 * t51 + t18 * t37 + t36) * g(2) + (t5 * rSges(7,1) - t4 * rSges(7,2) - t53 + t51 * t18 + (t17 * t39 + t35) * t16) * g(1)) (-m(3) + t38) * g(3), t38 * (g(1) * t16 - g(2) * t18) -m(5) * g(3) * t33 - m(6) * (g(3) * (t31 + t43) + t34) - m(7) * (g(3) * (t15 * t32 + t37 + t43) + t34) + t58 * ((m(5) * rSges(5,2) - m(6) * rSges(6,3) - m(7) * t32) * t17 + (m(5) * rSges(5,1) - m(6) * (rSges(6,2) - pkin(4)) - m(7) * t39) * t15) t57 * (t15 * t58 - t54) -m(7) * (g(1) * (rSges(7,1) * t2 - rSges(7,2) * t3) + g(2) * (rSges(7,1) * t4 + rSges(7,2) * t5) + (-rSges(7,1) * t27 + rSges(7,2) * t25) * t54)];
taug  = t1(:);
