% Calculate Gravitation load on the joints for
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPPRR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPPRR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:17:43
% EndTime: 2019-03-08 19:17:45
% DurationCPUTime: 0.87s
% Computational Cost: add. (430->100), mult. (1072->158), div. (0->0), fcn. (1322->12), ass. (0->51)
t36 = sin(qJ(2));
t39 = cos(qJ(2));
t60 = sin(pkin(11));
t61 = cos(pkin(11));
t23 = t36 * t60 - t39 * t61;
t31 = sin(pkin(10));
t33 = cos(pkin(10));
t62 = cos(pkin(6));
t48 = t62 * t60;
t49 = t62 * t61;
t64 = -t36 * t49 - t39 * t48;
t14 = -t33 * t23 + t31 * t64;
t24 = -t36 * t61 - t39 * t60;
t32 = sin(pkin(6));
t22 = t24 * t32;
t9 = t31 * t23 + t33 * t64;
t84 = -g(1) * t14 + g(2) * t9 + g(3) * t22;
t40 = -t36 * t48 + t39 * t49;
t13 = t24 * t33 - t31 * t40;
t55 = t39 * t62;
t44 = -t31 * t55 - t33 * t36;
t43 = t44 * pkin(2);
t42 = t13 * pkin(3) + t43;
t83 = -t31 * t36 + t33 * t55;
t81 = rSges(6,3) + pkin(8);
t10 = t31 * t24 + t33 * t40;
t21 = t23 * t32;
t80 = -g(1) * t13 - g(2) * t10 + g(3) * t21;
t72 = rSges(7,3) + pkin(9);
t35 = sin(qJ(5));
t69 = t32 * t35;
t38 = cos(qJ(5));
t68 = t32 * t38;
t29 = t32 * t39 * pkin(2);
t65 = -t21 * pkin(3) + t29;
t63 = rSges(5,3) + qJ(4);
t59 = -m(5) - m(6) - m(7);
t58 = -m(4) + t59;
t57 = t36 * t62;
t53 = t83 * pkin(2);
t50 = t10 * pkin(3) + t53;
t34 = sin(qJ(6));
t37 = cos(qJ(6));
t47 = rSges(7,1) * t37 - rSges(7,2) * t34 + pkin(5);
t16 = t21 * t35 + t62 * t38;
t15 = t21 * t38 - t62 * t35;
t5 = t10 * t35 + t33 * t68;
t4 = -t10 * t38 + t33 * t69;
t3 = -t13 * t35 + t31 * t68;
t2 = -t13 * t38 - t31 * t69;
t1 = [(-m(2) - m(3) + t58) * g(3), -m(3) * (g(1) * (t44 * rSges(3,1) + (t31 * t57 - t33 * t39) * rSges(3,2)) + g(2) * (t83 * rSges(3,1) + (-t31 * t39 - t33 * t57) * rSges(3,2)) + g(3) * (rSges(3,1) * t39 - rSges(3,2) * t36) * t32) - m(4) * (g(1) * (t13 * rSges(4,1) - rSges(4,2) * t14 + t43) + g(2) * (rSges(4,1) * t10 + rSges(4,2) * t9 + t53) + g(3) * (-rSges(4,1) * t21 + rSges(4,2) * t22 + t29)) - m(5) * (g(1) * (-t13 * rSges(5,2) + t14 * t63 + t42) + g(2) * (-rSges(5,2) * t10 - t63 * t9 + t50) + g(3) * (rSges(5,2) * t21 - t63 * t22 + t65)) + (-g(1) * t42 - g(2) * t50 - g(3) * t65 + t80 * (t34 * rSges(7,1) + t37 * rSges(7,2) + pkin(8)) + t84 * (t47 * t35 - t72 * t38 + qJ(4))) * m(7) + (-g(1) * (t81 * t13 + t42) - g(2) * (t81 * t10 + t50) - g(3) * (-t81 * t21 + t65) + t84 * (rSges(6,1) * t35 + rSges(6,2) * t38 + qJ(4))) * m(6), t58 * (g(3) * t62 + (g(1) * t31 - g(2) * t33) * t32) t59 * t80, -m(6) * (g(1) * (rSges(6,1) * t2 - rSges(6,2) * t3) + g(2) * (rSges(6,1) * t4 + rSges(6,2) * t5) + g(3) * (rSges(6,1) * t15 - rSges(6,2) * t16)) - m(7) * (g(2) * (t47 * t4 - t72 * t5) + (t47 * t15 + t72 * t16) * g(3) + (t47 * t2 + t72 * t3) * g(1)) -m(7) * (g(1) * ((t14 * t37 - t3 * t34) * rSges(7,1) + (-t14 * t34 - t3 * t37) * rSges(7,2)) + g(2) * ((t34 * t5 - t37 * t9) * rSges(7,1) + (t34 * t9 + t37 * t5) * rSges(7,2)) + g(3) * ((-t16 * t34 - t22 * t37) * rSges(7,1) + (-t16 * t37 + t22 * t34) * rSges(7,2)))];
taug  = t1(:);
