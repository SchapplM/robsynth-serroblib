% Calculate Gravitation load on the joints for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPPR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 08:59:04
% EndTime: 2022-01-23 08:59:06
% DurationCPUTime: 0.46s
% Computational Cost: add. (167->93), mult. (323->147), div. (0->0), fcn. (353->10), ass. (0->48)
t30 = cos(pkin(9));
t31 = cos(pkin(8));
t33 = sin(qJ(5));
t28 = sin(pkin(8));
t35 = cos(qJ(5));
t51 = t28 * t35;
t13 = t30 * t51 - t33 * t31;
t34 = sin(qJ(1));
t36 = cos(qJ(1));
t27 = sin(pkin(9));
t29 = sin(pkin(7));
t32 = cos(pkin(7));
t48 = t30 * t31;
t11 = t29 * t27 + t32 * t48;
t52 = t28 * t33;
t38 = t11 * t35 + t32 * t52;
t58 = t36 * t13 - t38 * t34;
t12 = t30 * t52 + t35 * t31;
t3 = -t11 * t33 + t32 * t51;
t57 = -t34 * t12 + t3 * t36;
t56 = -m(5) - m(6);
t54 = t29 * qJ(3) + pkin(1);
t53 = t28 * qJ(4) + pkin(2);
t50 = t29 * t34;
t49 = t29 * t36;
t46 = t34 * t28;
t45 = t34 * t31;
t44 = t36 * t28;
t43 = t36 * t31;
t42 = -m(4) + t56;
t41 = rSges(4,3) * t29 + pkin(2) * t32 + t54;
t40 = qJ(4) * t31 - qJ(2);
t39 = rSges(3,1) * t32 - rSges(3,2) * t29 + pkin(1);
t20 = t30 * pkin(4) + t27 * pkin(6) + pkin(3);
t37 = -t20 * t28 + t40;
t26 = t36 * qJ(2);
t25 = t34 * qJ(2);
t18 = -t28 * pkin(3) + t40;
t17 = t32 * t43 + t46;
t16 = t32 * t44 - t45;
t15 = -t32 * t45 + t44;
t14 = t32 * t46 + t43;
t10 = -t32 * t27 + t29 * t48;
t7 = (t31 * pkin(3) + t53) * t32 + t54;
t6 = t17 * t27 - t30 * t49;
t5 = t15 * t27 + t30 * t50;
t2 = (t20 * t31 + t53) * t32 + pkin(1) + (t27 * pkin(4) - t30 * pkin(6) + qJ(3)) * t29;
t1 = [-m(2) * (g(1) * (-t34 * rSges(2,1) - t36 * rSges(2,2)) + g(2) * (t36 * rSges(2,1) - t34 * rSges(2,2))) - m(3) * (g(1) * t26 + g(2) * t25 + (g(1) * rSges(3,3) + g(2) * t39) * t36 + (g(2) * rSges(3,3) - g(1) * t39) * t34) - m(4) * (g(1) * (t15 * rSges(4,1) + t14 * rSges(4,2) - t41 * t34 + t26) + g(2) * (t17 * rSges(4,1) - t16 * rSges(4,2) + t41 * t36 + t25)) - m(5) * (g(1) * (-t7 * t34 - t18 * t36 + (t15 * t30 - t27 * t50) * rSges(5,1) - t5 * rSges(5,2) - t14 * rSges(5,3)) + g(2) * (t7 * t36 - t18 * t34 + (t17 * t30 + t27 * t49) * rSges(5,1) - t6 * rSges(5,2) + t16 * rSges(5,3))) - m(6) * (g(1) * (-t2 * t34 - t37 * t36 + t58 * rSges(6,1) + (-t36 * t12 - t3 * t34) * rSges(6,2) + t5 * rSges(6,3)) + g(2) * (t2 * t36 - t37 * t34 + ((t11 * t36 + t30 * t46) * t35 + t16 * t33) * rSges(6,1) + t57 * rSges(6,2) + t6 * rSges(6,3))), (-m(3) + t42) * (g(1) * t34 - g(2) * t36), t42 * (-g(3) * t32 + (g(1) * t36 + g(2) * t34) * t29), t56 * (g(3) * t29 * t28 + g(1) * t16 + g(2) * t14), -m(6) * (g(1) * (t57 * rSges(6,1) + (-t34 * t13 - t36 * t38) * rSges(6,2)) + g(2) * ((-(t11 * t34 - t30 * t44) * t33 + t14 * t35) * rSges(6,1) + t58 * rSges(6,2)) + g(3) * ((-t10 * t33 + t29 * t51) * rSges(6,1) + (-t10 * t35 - t29 * t52) * rSges(6,2)))];
taug = t1(:);
