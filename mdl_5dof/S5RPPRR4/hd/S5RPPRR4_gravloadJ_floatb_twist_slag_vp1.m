% Calculate Gravitation load on the joints for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:16:17
% EndTime: 2022-01-23 09:16:18
% DurationCPUTime: 0.36s
% Computational Cost: add. (217->87), mult. (249->119), div. (0->0), fcn. (239->10), ass. (0->45)
t31 = sin(pkin(9));
t32 = sin(pkin(8));
t33 = cos(pkin(9));
t34 = cos(pkin(8));
t58 = pkin(1) - (-rSges(4,3) - qJ(3)) * t32 + (rSges(4,1) * t33 - rSges(4,2) * t31 + pkin(2)) * t34;
t30 = pkin(9) + qJ(4);
t24 = qJ(5) + t30;
t21 = cos(t24);
t36 = cos(qJ(1));
t20 = sin(t24);
t35 = sin(qJ(1));
t49 = t35 * t20;
t5 = t21 * t36 + t34 * t49;
t48 = t35 * t21;
t6 = t20 * t36 - t34 * t48;
t57 = -t5 * rSges(6,1) + t6 * rSges(6,2);
t50 = t34 * t36;
t7 = -t20 * t50 + t48;
t8 = t21 * t50 + t49;
t56 = t7 * rSges(6,1) - t8 * rSges(6,2);
t55 = pkin(3) * t31;
t22 = sin(t30);
t54 = pkin(4) * t22;
t53 = g(3) * t32;
t52 = t33 * pkin(3) + pkin(2);
t51 = rSges(3,2) * t32;
t47 = t35 * t22;
t23 = cos(t30);
t46 = t35 * t23;
t45 = qJ(3) + pkin(6);
t25 = t35 * qJ(2);
t44 = t36 * pkin(1) + t25;
t43 = -m(4) - m(5) - m(6);
t41 = t52 * t34 + pkin(1) + (rSges(5,3) + t45) * t32;
t39 = rSges(4,1) * t31 + rSges(4,2) * t33;
t38 = -rSges(6,1) * t20 - rSges(6,2) * t21;
t11 = -t22 * t50 + t46;
t9 = t23 * t36 + t34 * t47;
t37 = (pkin(4) * t23 + t52) * t34 + (rSges(6,3) + pkin(7) + t45) * t32;
t27 = t36 * qJ(2);
t19 = qJ(2) + t55;
t16 = t54 + t55;
t12 = t23 * t50 + t47;
t10 = t22 * t36 - t34 * t46;
t1 = [-m(2) * (g(1) * (-t35 * rSges(2,1) - rSges(2,2) * t36) + g(2) * (rSges(2,1) * t36 - t35 * rSges(2,2))) - m(3) * (g(1) * (t36 * rSges(3,3) + t27) + g(2) * (rSges(3,1) * t50 - t36 * t51 + t44) + (g(1) * (-rSges(3,1) * t34 - pkin(1) + t51) + g(2) * rSges(3,3)) * t35) - m(4) * (g(1) * (-t58 * t35 + t39 * t36 + t27) + g(2) * (t39 * t35 + t58 * t36 + t25)) - m(5) * (g(1) * (t10 * rSges(5,1) + t9 * rSges(5,2) + t19 * t36 - t41 * t35) + g(2) * (t12 * rSges(5,1) + t11 * rSges(5,2) + t19 * t35 + t41 * t36)) - m(6) * (g(1) * (t6 * rSges(6,1) + t5 * rSges(6,2) + t27) + g(2) * (t8 * rSges(6,1) + t7 * rSges(6,2) + t44) + (g(1) * t16 + g(2) * t37) * t36 + (g(1) * (-pkin(1) - t37) + g(2) * t16) * t35), (-m(3) + t43) * (g(1) * t35 - g(2) * t36), t43 * (-g(3) * t34 + (g(1) * t36 + g(2) * t35) * t32), -m(5) * (g(1) * (rSges(5,1) * t11 - rSges(5,2) * t12) + g(2) * (-rSges(5,1) * t9 + rSges(5,2) * t10)) - m(6) * (g(1) * (t11 * pkin(4) + t56) + g(2) * (-t9 * pkin(4) + t57)) + (-m(5) * (-rSges(5,1) * t22 - rSges(5,2) * t23) - m(6) * (t38 - t54)) * t53, -m(6) * (g(1) * t56 + g(2) * t57 + t38 * t53)];
taug = t1(:);
