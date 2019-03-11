% Calculate Gravitation load on the joints for
% S6PRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-03-08 19:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRPR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:39:05
% EndTime: 2019-03-08 19:39:06
% DurationCPUTime: 0.60s
% Computational Cost: add. (506->102), mult. (796->149), div. (0->0), fcn. (920->14), ass. (0->57)
t33 = sin(pkin(6));
t39 = cos(qJ(2));
t58 = t33 * t39;
t80 = g(3) * t58;
t31 = sin(pkin(12));
t34 = cos(pkin(12));
t79 = rSges(6,1) * t31 + rSges(6,2) * t34;
t38 = sin(qJ(2));
t52 = sin(pkin(10));
t54 = cos(pkin(6));
t46 = t54 * t52;
t53 = cos(pkin(10));
t15 = t53 * t38 + t39 * t46;
t47 = t54 * t53;
t13 = t52 * t38 - t39 * t47;
t67 = g(2) * t13;
t78 = g(1) * t15 + t67;
t14 = t38 * t47 + t52 * t39;
t30 = pkin(11) + qJ(4);
t26 = sin(t30);
t28 = cos(t30);
t50 = t33 * t53;
t3 = t14 * t26 + t28 * t50;
t16 = -t38 * t46 + t53 * t39;
t49 = t33 * t52;
t5 = t16 * t26 - t28 * t49;
t59 = t33 * t38;
t9 = t26 * t59 - t54 * t28;
t77 = g(1) * t5 + g(2) * t3 + g(3) * t9;
t10 = t54 * t26 + t28 * t59;
t4 = t14 * t28 - t26 * t50;
t6 = t16 * t28 + t26 * t49;
t76 = g(1) * t6 + g(2) * t4 + g(3) * t10;
t29 = pkin(12) + qJ(6);
t25 = sin(t29);
t27 = cos(t29);
t43 = rSges(7,1) * t27 - rSges(7,2) * t25 + pkin(5) * t34 + pkin(4);
t57 = rSges(7,3) + pkin(9) + qJ(5);
t75 = t57 * t26 + t43 * t28;
t44 = rSges(6,1) * t34 - rSges(6,2) * t31 + pkin(4);
t55 = rSges(6,3) + qJ(5);
t74 = t55 * t26 + t44 * t28;
t69 = -m(6) - m(7);
t66 = g(2) * t14;
t35 = cos(pkin(11));
t24 = pkin(3) * t35 + pkin(2);
t65 = t24 * t80;
t64 = g(3) * t33;
t37 = -pkin(8) - qJ(3);
t63 = -t13 * t24 - t14 * t37;
t62 = -t15 * t24 - t16 * t37;
t56 = rSges(4,3) + qJ(3);
t51 = -m(4) - m(5) + t69;
t48 = rSges(5,1) * t28 - rSges(5,2) * t26;
t45 = rSges(4,1) * t35 - rSges(4,2) * sin(pkin(11)) + pkin(2);
t42 = rSges(7,1) * t25 + rSges(7,2) * t27 + pkin(5) * t31;
t1 = [(-m(2) - m(3) + t51) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t15 - rSges(3,2) * t16) + g(2) * (-rSges(3,1) * t13 - rSges(3,2) * t14) + (rSges(3,1) * t39 - rSges(3,2) * t38) * t64) - m(4) * (g(1) * (-t45 * t15 + t56 * t16) + t56 * t66 - t45 * t67 + (t56 * t38 + t45 * t39) * t64) - m(5) * (g(1) * (rSges(5,3) * t16 - t48 * t15 + t62) + g(2) * (rSges(5,3) * t14 - t48 * t13 + t63) + t65 + (t48 * t39 + (rSges(5,3) - t37) * t38) * t64) - m(6) * (g(1) * (t79 * t16 + t62) + g(2) * (t79 * t14 + t63) + t65 + ((-t37 + t79) * t38 + t74 * t39) * t64 - t78 * t74) - m(7) * (g(2) * t63 + t65 + t42 * t66 + ((-t37 + t42) * t38 + t75 * t39) * t64 - t78 * t75 + (t42 * t16 + t62) * g(1)) t51 * (t78 - t80) -m(5) * (g(1) * (-rSges(5,1) * t5 - rSges(5,2) * t6) + g(2) * (-rSges(5,1) * t3 - rSges(5,2) * t4) + g(3) * (-rSges(5,1) * t9 - rSges(5,2) * t10)) - m(6) * (-t77 * t44 + t76 * t55) - m(7) * (-t77 * t43 + t76 * t57) t69 * t77, -m(7) * (g(1) * ((t15 * t27 - t25 * t6) * rSges(7,1) + (-t15 * t25 - t27 * t6) * rSges(7,2)) + g(2) * ((t13 * t27 - t25 * t4) * rSges(7,1) + (-t13 * t25 - t27 * t4) * rSges(7,2)) + g(3) * ((-t10 * t25 - t27 * t58) * rSges(7,1) + (-t10 * t27 + t25 * t58) * rSges(7,2)))];
taug  = t1(:);
