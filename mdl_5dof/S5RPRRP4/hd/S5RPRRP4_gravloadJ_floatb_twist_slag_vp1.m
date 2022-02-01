% Calculate Gravitation load on the joints for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% m [6x1]
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
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:31:51
% EndTime: 2022-01-23 09:31:52
% DurationCPUTime: 0.39s
% Computational Cost: add. (233->103), mult. (316->139), div. (0->0), fcn. (308->8), ass. (0->49)
t64 = pkin(7) + pkin(6);
t35 = cos(pkin(8));
t33 = qJ(3) + qJ(4);
t26 = cos(t33);
t39 = cos(qJ(1));
t48 = t39 * t26;
t25 = sin(t33);
t37 = sin(qJ(1));
t53 = t37 * t25;
t10 = t35 * t53 + t48;
t49 = t39 * t25;
t52 = t37 * t26;
t11 = -t35 * t52 + t49;
t63 = -t10 * rSges(5,1) + t11 * rSges(5,2);
t62 = -t10 * rSges(6,1) + t11 * rSges(6,2);
t12 = -t35 * t49 + t52;
t13 = t35 * t48 + t53;
t61 = t12 * rSges(6,1) - t13 * rSges(6,2);
t60 = t12 * rSges(5,1) - t13 * rSges(5,2);
t59 = pkin(4) * t25;
t34 = sin(pkin(8));
t58 = g(3) * t34;
t36 = sin(qJ(3));
t57 = t36 * pkin(3);
t38 = cos(qJ(3));
t30 = t38 * pkin(3);
t56 = pkin(2) * t35 + pkin(1);
t55 = rSges(3,2) * t34;
t54 = t35 * t39;
t51 = t37 * t36;
t50 = t37 * t38;
t47 = t39 * t36;
t46 = t39 * t38;
t20 = pkin(4) * t26 + t30;
t27 = t37 * qJ(2);
t45 = t39 * pkin(1) + t27;
t44 = t35 * t30 + t56 + (rSges(5,3) + t64) * t34;
t43 = t56 + (rSges(4,3) + pkin(6)) * t34;
t42 = -rSges(5,1) * t25 - rSges(5,2) * t26;
t41 = -rSges(6,1) * t25 - rSges(6,2) * t26;
t16 = -t35 * t47 + t50;
t14 = t35 * t51 + t46;
t40 = (pkin(2) + t20) * t35 + (rSges(6,3) + qJ(5) + t64) * t34;
t28 = t39 * qJ(2);
t23 = qJ(2) + t57;
t19 = t57 + t59;
t17 = t35 * t46 + t51;
t15 = -t35 * t50 + t47;
t1 = [-m(2) * (g(1) * (-t37 * rSges(2,1) - t39 * rSges(2,2)) + g(2) * (t39 * rSges(2,1) - t37 * rSges(2,2))) - m(3) * (g(1) * (t39 * rSges(3,3) + t28) + g(2) * (rSges(3,1) * t54 - t39 * t55 + t45) + (g(1) * (-rSges(3,1) * t35 - pkin(1) + t55) + g(2) * rSges(3,3)) * t37) - m(4) * (g(1) * (t15 * rSges(4,1) + t14 * rSges(4,2) - t43 * t37 + t28) + g(2) * (t17 * rSges(4,1) + t16 * rSges(4,2) + t43 * t39 + t27)) - m(5) * (g(1) * (t11 * rSges(5,1) + t10 * rSges(5,2) + t23 * t39 - t44 * t37) + g(2) * (t13 * rSges(5,1) + t12 * rSges(5,2) + t23 * t37 + t44 * t39)) - m(6) * (g(1) * (t11 * rSges(6,1) + t10 * rSges(6,2) + t28) + g(2) * (t13 * rSges(6,1) + t12 * rSges(6,2) + t45) + (g(1) * t19 + g(2) * t40) * t39 + (g(1) * (-pkin(1) - t40) + g(2) * t19) * t37), (-m(3) - m(4) - m(5) - m(6)) * (g(1) * t37 - g(2) * t39), -m(4) * (g(1) * (t16 * rSges(4,1) - t17 * rSges(4,2)) + g(2) * (-t14 * rSges(4,1) + t15 * rSges(4,2))) - m(5) * (g(1) * (t16 * pkin(3) + t60) + g(2) * (-t14 * pkin(3) + t63)) - m(6) * (g(1) * (-t19 * t54 + t37 * t20 + t61) + g(2) * (-t37 * t35 * t19 - t39 * t20 + t62)) + (-m(4) * (-rSges(4,1) * t36 - rSges(4,2) * t38) - m(5) * (t42 - t57) - m(6) * (-t19 + t41)) * t58, -m(5) * (g(1) * t60 + g(2) * t63) - m(6) * (g(1) * (t12 * pkin(4) + t61) + g(2) * (-t10 * pkin(4) + t62)) + (-m(5) * t42 - m(6) * (t41 - t59)) * t58, -m(6) * (-g(3) * t35 + (g(1) * t39 + g(2) * t37) * t34)];
taug = t1(:);
