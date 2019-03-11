% Calculate Gravitation load on the joints for
% S6RPPRRR8
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
% Datum: 2019-03-09 02:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRR8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR8_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:35:25
% EndTime: 2019-03-09 02:35:26
% DurationCPUTime: 0.58s
% Computational Cost: add. (299->103), mult. (346->140), div. (0->0), fcn. (315->10), ass. (0->52)
t32 = sin(qJ(1));
t34 = cos(qJ(1));
t75 = -g(1) * t32 + g(2) * t34;
t49 = rSges(7,3) + pkin(9) + pkin(8);
t72 = g(1) * t34 + g(2) * t32;
t33 = cos(qJ(5));
t18 = t33 * pkin(5) + pkin(4);
t27 = qJ(5) + qJ(6);
t21 = sin(t27);
t22 = cos(t27);
t31 = sin(qJ(5));
t36 = m(6) * (rSges(6,1) * t33 - rSges(6,2) * t31 + pkin(4)) + m(7) * (rSges(7,1) * t22 - rSges(7,2) * t21 + t18) + m(5) * rSges(5,1);
t58 = rSges(6,3) + pkin(8);
t71 = -m(5) * rSges(5,2) + m(6) * t58 + m(7) * t49;
t26 = pkin(10) + qJ(4);
t19 = sin(t26);
t52 = t34 * t22;
t57 = t32 * t21;
t5 = -t19 * t57 + t52;
t53 = t34 * t21;
t56 = t32 * t22;
t6 = t19 * t56 + t53;
t70 = t5 * rSges(7,1) - t6 * rSges(7,2);
t7 = t19 * t53 + t56;
t8 = t19 * t52 - t57;
t69 = t7 * rSges(7,1) + t8 * rSges(7,2);
t28 = sin(pkin(10));
t66 = pkin(3) * t28;
t65 = pkin(5) * t31;
t20 = cos(t26);
t60 = g(3) * t20;
t59 = t19 * pkin(4);
t55 = t32 * t31;
t54 = t32 * t33;
t51 = t34 * t31;
t50 = t34 * t33;
t48 = t34 * pkin(1) + t32 * qJ(2);
t47 = rSges(4,3) + qJ(3);
t46 = t32 * t66 + t48;
t24 = t34 * qJ(2);
t30 = -pkin(7) - qJ(3);
t45 = t32 * t30 + t34 * t66 + t24;
t44 = -m(4) - m(5) - m(6) - m(7);
t43 = rSges(4,1) * t28 + rSges(4,2) * cos(pkin(10));
t42 = t19 * rSges(5,1) + t20 * rSges(5,2);
t41 = -rSges(7,1) * t21 - rSges(7,2) * t22;
t11 = t19 * t51 + t54;
t9 = -t19 * t55 + t50;
t38 = t19 * t18 - t49 * t20;
t12 = t19 * t50 - t55;
t10 = t19 * t54 + t51;
t1 = [-m(2) * (g(1) * (-t32 * rSges(2,1) - t34 * rSges(2,2)) + g(2) * (t34 * rSges(2,1) - t32 * rSges(2,2))) - m(3) * (g(1) * (t34 * rSges(3,3) + t24 + (rSges(3,2) - pkin(1)) * t32) + g(2) * (-t34 * rSges(3,2) + t32 * rSges(3,3) + t48)) - m(4) * (g(1) * t24 + g(2) * t48 + (g(1) * t43 + g(2) * t47) * t34 + (g(1) * (-pkin(1) - t47) + g(2) * t43) * t32) - m(5) * (g(1) * t45 + g(2) * t46 + (g(1) * t42 + g(2) * (rSges(5,3) - t30)) * t34 + (g(1) * (-rSges(5,3) - pkin(1)) + g(2) * t42) * t32) - m(6) * (g(1) * (t12 * rSges(6,1) - t11 * rSges(6,2) - t32 * pkin(1) + t34 * t59 + t45) + g(2) * (t10 * rSges(6,1) + t9 * rSges(6,2) - t34 * t30 + t32 * t59 + t46) - t72 * t20 * t58) - m(7) * (g(1) * (t8 * rSges(7,1) - t7 * rSges(7,2) + t45) + g(2) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t46) + (g(1) * t38 + g(2) * (-t30 + t65)) * t34 + (g(1) * (-pkin(1) - t65) + g(2) * t38) * t32) -(-m(3) + t44) * t75, t44 * t72 (t19 * t36 - t20 * t71) * g(3) + t75 * (t71 * t19 + t36 * t20) -m(6) * (g(1) * (t9 * rSges(6,1) - t10 * rSges(6,2)) + g(2) * (t11 * rSges(6,1) + t12 * rSges(6,2))) - m(7) * (g(1) * (t9 * pkin(5) + t70) + g(2) * (t11 * pkin(5) + t69)) + (-m(6) * (-rSges(6,1) * t31 - rSges(6,2) * t33) - m(7) * (t41 - t65)) * t60, -m(7) * (g(1) * t70 + g(2) * t69 + t41 * t60)];
taug  = t1(:);
