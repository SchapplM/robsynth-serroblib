% Calculate Gravitation load on the joints for
% S5RRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP9_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP9_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP9_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP9_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:05:42
% EndTime: 2019-12-31 20:05:44
% DurationCPUTime: 0.60s
% Computational Cost: add. (274->106), mult. (395->150), div. (0->0), fcn. (382->8), ass. (0->45)
t50 = rSges(6,1) + pkin(4);
t26 = cos(qJ(1));
t24 = sin(qJ(1));
t52 = g(2) * t24;
t57 = g(1) * t26 + t52;
t39 = rSges(6,3) + qJ(5);
t19 = pkin(8) + qJ(4);
t14 = sin(t19);
t15 = cos(t19);
t56 = t39 * t14 + t50 * t15;
t21 = cos(pkin(8));
t13 = pkin(3) * t21 + pkin(2);
t25 = cos(qJ(2));
t7 = t25 * t13;
t55 = g(3) * t7;
t54 = g(1) * t24;
t23 = sin(qJ(2));
t51 = g(3) * t23;
t49 = rSges(3,2) * t23;
t20 = sin(pkin(8));
t48 = t20 * t26;
t47 = t24 * t20;
t46 = t24 * t25;
t45 = t25 * t26;
t44 = t26 * t14;
t22 = -pkin(7) - qJ(3);
t43 = rSges(6,2) - t22;
t42 = rSges(5,3) - t22;
t41 = t26 * pkin(1) + t24 * pkin(6);
t40 = rSges(4,3) + qJ(3);
t17 = t26 * pkin(6);
t38 = t24 * t23 * t22 + pkin(3) * t48 + t17;
t37 = -pkin(1) - t7;
t36 = t43 * t26;
t35 = t42 * t26;
t34 = t26 * t40;
t33 = pkin(3) * t47 + t13 * t45 + t41;
t32 = rSges(3,1) * t25 - t49;
t30 = rSges(5,1) * t15 - rSges(5,2) * t14;
t29 = rSges(4,1) * t21 - rSges(4,2) * t20 + pkin(2);
t4 = t24 * t14 + t15 * t45;
t3 = -t24 * t15 + t25 * t44;
t2 = t15 * t46 - t44;
t1 = t14 * t46 + t15 * t26;
t5 = [-m(2) * (g(1) * (-t24 * rSges(2,1) - rSges(2,2) * t26) + g(2) * (rSges(2,1) * t26 - t24 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t26 + t17) + g(2) * (rSges(3,1) * t45 - t26 * t49 + t41) + (g(1) * (-pkin(1) - t32) + g(2) * rSges(3,3)) * t24) - m(4) * (g(1) * (-pkin(2) * t46 - t24 * pkin(1) + t17 + (-t21 * t46 + t48) * rSges(4,1) + (t20 * t46 + t21 * t26) * rSges(4,2)) + g(2) * (pkin(2) * t45 + (t21 * t45 + t47) * rSges(4,1) + (-t20 * t45 + t24 * t21) * rSges(4,2) + t41) + (g(2) * t34 - t40 * t54) * t23) - m(5) * (g(1) * (-rSges(5,1) * t2 + rSges(5,2) * t1 + t38) + g(2) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t23 * t35 + t33) + (-rSges(5,3) * t23 + t37) * t54) - m(6) * (g(1) * (-t39 * t1 - t50 * t2 + t38) + g(2) * (t23 * t36 + t39 * t3 + t50 * t4 + t33) + (-rSges(6,2) * t23 + t37) * t54), -m(3) * (g(3) * t32 + t57 * (-rSges(3,1) * t23 - rSges(3,2) * t25)) - m(4) * ((g(1) * t34 + g(3) * t29 + t40 * t52) * t25 + (g(3) * t40 - t57 * t29) * t23) - m(5) * (t55 + (g(1) * t35 + g(3) * t30 + t42 * t52) * t25 + (g(3) * t42 + t57 * (-t13 - t30)) * t23) - m(6) * (t55 + (g(1) * t36 + g(3) * t56 + t43 * t52) * t25 + (g(3) * t43 + t57 * (-t13 - t56)) * t23), (-m(4) - m(5) - m(6)) * (-g(3) * t25 + t57 * t23), -m(5) * (g(1) * (-rSges(5,1) * t3 - rSges(5,2) * t4) + g(2) * (-rSges(5,1) * t1 - rSges(5,2) * t2)) - m(6) * (g(1) * (-t50 * t3 + t39 * t4) + g(2) * (-t50 * t1 + t39 * t2)) + (-m(5) * (-rSges(5,1) * t14 - rSges(5,2) * t15) - m(6) * (-t50 * t14 + t39 * t15)) * t51, -m(6) * (g(1) * t3 + g(2) * t1 + t14 * t51)];
taug = t5(:);
