% Calculate Gravitation load on the joints for
% S6RPPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 01:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRPR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:46:26
% EndTime: 2019-03-09 01:46:27
% DurationCPUTime: 0.42s
% Computational Cost: add. (259->85), mult. (396->117), div. (0->0), fcn. (435->10), ass. (0->42)
t21 = qJ(4) + pkin(10);
t15 = sin(t21);
t59 = rSges(7,3) + pkin(8);
t60 = t59 * t15;
t40 = sin(pkin(9));
t41 = cos(pkin(9));
t46 = sin(qJ(1));
t47 = cos(qJ(1));
t8 = -t40 * t46 - t41 * t47;
t9 = t40 * t47 - t41 * t46;
t58 = g(1) * t9 - g(2) * t8;
t23 = sin(qJ(6));
t25 = cos(qJ(6));
t57 = m(7) * (t25 * rSges(7,1) - t23 * rSges(7,2) + pkin(5)) + m(6) * rSges(6,1);
t54 = -m(6) - m(7);
t24 = sin(qJ(4));
t52 = pkin(4) * t24;
t16 = cos(t21);
t51 = t16 * pkin(5);
t26 = cos(qJ(4));
t50 = t26 * pkin(4);
t49 = rSges(5,3) + pkin(7);
t45 = t15 * rSges(6,2);
t44 = t16 * t23;
t43 = t16 * t25;
t42 = t47 * pkin(1) + t46 * qJ(2);
t39 = t47 * pkin(2) + t42;
t38 = m(4) + m(5) - t54;
t14 = pkin(3) + t50;
t22 = -qJ(5) - pkin(7);
t37 = -t8 * t14 - t9 * t22 + t39;
t36 = -pkin(1) * t46 + t47 * qJ(2);
t35 = -t26 * rSges(5,1) + t24 * rSges(5,2);
t34 = t16 * rSges(6,1) - t45;
t33 = t8 * t23 + t43 * t9;
t32 = -t8 * t25 + t44 * t9;
t31 = pkin(3) - t35;
t29 = -pkin(2) * t46 + t36;
t28 = t9 * t14 - t8 * t22 + t29;
t2 = t9 * t23 - t43 * t8;
t1 = t9 * t25 + t44 * t8;
t3 = [-m(2) * (g(1) * (-rSges(2,1) * t46 - rSges(2,2) * t47) + g(2) * (rSges(2,1) * t47 - rSges(2,2) * t46)) - m(3) * (g(1) * (-rSges(3,1) * t46 + rSges(3,3) * t47 + t36) + g(2) * (rSges(3,1) * t47 + rSges(3,3) * t46 + t42)) - m(4) * (g(1) * (t9 * rSges(4,1) - t8 * rSges(4,2) + t29) + g(2) * (-t8 * rSges(4,1) - t9 * rSges(4,2) + t39)) - m(5) * (g(1) * t29 + g(2) * t39 + (g(1) * t31 + g(2) * t49) * t9 + (g(1) * t49 - g(2) * t31) * t8) - m(6) * (g(1) * (t8 * rSges(6,3) + t34 * t9 + t28) + g(2) * (t9 * rSges(6,3) - t34 * t8 + t37)) - m(7) * (g(1) * (rSges(7,1) * t33 - rSges(7,2) * t32 + t51 * t9 + t28) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) - t51 * t8 + t37) + t58 * t60) (-m(3) - t38) * (g(1) * t46 - g(2) * t47) t38 * g(3) (-m(5) * t35 - m(6) * (t45 - t50) - m(7) * (-t50 - t60) + t57 * t16) * g(3) + (g(1) * t8 + g(2) * t9) * (-m(5) * (rSges(5,1) * t24 + rSges(5,2) * t26) - m(6) * (rSges(6,2) * t16 + t52) - m(7) * (-t59 * t16 + t52) - t57 * t15) t54 * t58, -m(7) * (g(1) * (t1 * rSges(7,1) - t2 * rSges(7,2)) + g(2) * (rSges(7,1) * t32 + rSges(7,2) * t33) + g(3) * (t23 * rSges(7,1) + t25 * rSges(7,2)) * t15)];
taug  = t3(:);
