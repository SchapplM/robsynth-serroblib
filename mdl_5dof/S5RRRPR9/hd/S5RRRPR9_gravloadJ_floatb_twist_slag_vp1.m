% Calculate Gravitation load on the joints for
% S5RRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR9_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR9_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR9_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR9_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:22:15
% EndTime: 2019-12-31 21:22:17
% DurationCPUTime: 0.57s
% Computational Cost: add. (316->120), mult. (410->171), div. (0->0), fcn. (389->10), ass. (0->58)
t34 = -qJ(4) - pkin(7);
t51 = rSges(5,3) - t34;
t50 = rSges(6,3) + pkin(8) - t34;
t37 = sin(qJ(1));
t40 = cos(qJ(1));
t70 = g(1) * t40 + g(2) * t37;
t36 = sin(qJ(2));
t39 = cos(qJ(2));
t61 = rSges(4,3) + pkin(7);
t69 = t39 * pkin(2) + t61 * t36;
t33 = qJ(3) + pkin(9);
t27 = qJ(5) + t33;
t22 = sin(t27);
t23 = cos(t27);
t56 = t40 * t23;
t59 = t37 * t39;
t5 = t22 * t59 + t56;
t57 = t40 * t22;
t6 = -t23 * t59 + t57;
t68 = -t5 * rSges(6,1) + t6 * rSges(6,2);
t7 = t37 * t23 - t39 * t57;
t8 = t37 * t22 + t39 * t56;
t67 = t7 * rSges(6,1) - t8 * rSges(6,2);
t64 = g(3) * t36;
t35 = sin(qJ(3));
t63 = t35 * pkin(3);
t60 = t36 * rSges(3,2);
t58 = t39 * t40;
t25 = sin(t33);
t55 = t40 * t25;
t26 = cos(t33);
t54 = t40 * t26;
t53 = t40 * t35;
t38 = cos(qJ(3));
t52 = t40 * t38;
t29 = t38 * pkin(3);
t20 = pkin(4) * t26 + t29;
t49 = t40 * pkin(1) + t37 * pkin(6);
t48 = t39 * rSges(3,1) - t60;
t46 = -rSges(6,1) * t22 - rSges(6,2) * t23;
t45 = rSges(4,1) * t38 - rSges(4,2) * t35 + pkin(2);
t16 = t37 * t38 - t39 * t53;
t14 = t35 * t59 + t52;
t24 = t29 + pkin(2);
t44 = rSges(5,1) * t26 - rSges(5,2) * t25 + t24;
t18 = pkin(2) + t20;
t43 = rSges(6,1) * t23 - rSges(6,2) * t22 + t18;
t42 = t39 * t24 + t51 * t36;
t41 = t39 * t18 + t50 * t36;
t30 = t40 * pkin(6);
t19 = pkin(4) * t25 + t63;
t17 = t37 * t35 + t39 * t52;
t15 = -t38 * t59 + t53;
t12 = t37 * t25 + t39 * t54;
t11 = t37 * t26 - t39 * t55;
t10 = -t26 * t59 + t55;
t9 = t25 * t59 + t54;
t1 = [-m(2) * (g(1) * (-t37 * rSges(2,1) - t40 * rSges(2,2)) + g(2) * (t40 * rSges(2,1) - t37 * rSges(2,2))) - m(3) * (g(1) * (t40 * rSges(3,3) + t30) + g(2) * (rSges(3,1) * t58 - t40 * t60 + t49) + (g(1) * (-pkin(1) - t48) + g(2) * rSges(3,3)) * t37) - m(4) * ((t17 * rSges(4,1) + t16 * rSges(4,2) + t69 * t40 + t49) * g(2) + (t15 * rSges(4,1) + t14 * rSges(4,2) + t30 + (-pkin(1) - t69) * t37) * g(1)) - m(5) * (g(1) * (t10 * rSges(5,1) + t9 * rSges(5,2) + t30) + g(2) * (t12 * rSges(5,1) + t11 * rSges(5,2) + t49) + (g(1) * t63 + g(2) * t42) * t40 + (g(1) * (-pkin(1) - t42) + g(2) * t63) * t37) - m(6) * (g(1) * (t6 * rSges(6,1) + t5 * rSges(6,2) + t30) + g(2) * (t8 * rSges(6,1) + t7 * rSges(6,2) + t49) + (g(1) * t19 + g(2) * t41) * t40 + (g(1) * (-pkin(1) - t41) + g(2) * t19) * t37), -m(3) * (g(3) * t48 + t70 * (-rSges(3,1) * t36 - rSges(3,2) * t39)) - m(4) * ((g(3) * t45 + t70 * t61) * t39 + (g(3) * t61 - t70 * t45) * t36) - m(5) * ((g(3) * t44 + t70 * t51) * t39 + (g(3) * t51 - t70 * t44) * t36) - m(6) * ((g(3) * t43 + t70 * t50) * t39 + (g(3) * t50 - t70 * t43) * t36), -m(4) * (g(1) * (t16 * rSges(4,1) - t17 * rSges(4,2)) + g(2) * (-t14 * rSges(4,1) + t15 * rSges(4,2))) - m(5) * (g(1) * (t11 * rSges(5,1) - t12 * rSges(5,2) + t16 * pkin(3)) + g(2) * (-t9 * rSges(5,1) + t10 * rSges(5,2) - t14 * pkin(3))) - m(6) * (g(1) * (-t19 * t58 + t37 * t20 + t67) + g(2) * (-t19 * t59 - t40 * t20 + t68)) + (-m(4) * (-rSges(4,1) * t35 - rSges(4,2) * t38) - m(5) * (-rSges(5,1) * t25 - rSges(5,2) * t26 - t63) - m(6) * (-t19 + t46)) * t64, (-m(5) - m(6)) * (-g(3) * t39 + t70 * t36), -m(6) * (g(1) * t67 + g(2) * t68 + t46 * t64)];
taug = t1(:);
