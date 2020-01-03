% Calculate Gravitation load on the joints for
% S5RRPRP7
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
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP7_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:59:55
% EndTime: 2019-12-31 19:59:57
% DurationCPUTime: 0.54s
% Computational Cost: add. (265->98), mult. (350->137), div. (0->0), fcn. (333->8), ass. (0->45)
t52 = rSges(6,1) + pkin(4);
t19 = qJ(2) + pkin(8);
t16 = sin(t19);
t17 = cos(t19);
t61 = t17 * pkin(3) + t16 * pkin(7);
t21 = sin(qJ(4));
t24 = cos(qJ(4));
t60 = rSges(5,1) * t24 - rSges(5,2) * t21;
t23 = sin(qJ(1));
t26 = cos(qJ(1));
t59 = g(1) * t26 + g(2) * t23;
t39 = rSges(6,3) + qJ(5);
t58 = t21 * t39 + t24 * t52;
t22 = sin(qJ(2));
t57 = pkin(2) * t22;
t20 = -qJ(3) - pkin(6);
t55 = g(2) * t20;
t53 = g(3) * t16;
t51 = rSges(3,3) + pkin(6);
t48 = t16 * t26;
t47 = t17 * t23;
t46 = t17 * t26;
t45 = t23 * t21;
t44 = t23 * t24;
t43 = t26 * t20;
t42 = t26 * t21;
t41 = t26 * t24;
t40 = rSges(4,3) - t20;
t25 = cos(qJ(2));
t18 = t25 * pkin(2);
t15 = t18 + pkin(1);
t10 = t26 * t15;
t38 = pkin(3) * t46 + pkin(7) * t48 + t10;
t37 = t18 + t61;
t36 = pkin(7) * t47 - t23 * t57;
t35 = pkin(7) * t46 - t26 * t57;
t34 = rSges(3,1) * t25 - rSges(3,2) * t22;
t32 = rSges(4,1) * t17 - rSges(4,2) * t16;
t31 = -t15 - t61;
t30 = pkin(1) + t34;
t4 = t17 * t41 + t45;
t3 = t17 * t42 - t44;
t2 = t17 * t44 - t42;
t1 = t17 * t45 + t41;
t5 = [-m(2) * (g(1) * (-rSges(2,1) * t23 - rSges(2,2) * t26) + g(2) * (rSges(2,1) * t26 - rSges(2,2) * t23)) - m(3) * ((g(1) * t51 + g(2) * t30) * t26 + (-g(1) * t30 + g(2) * t51) * t23) - m(4) * (g(2) * t10 + (g(1) * t40 + g(2) * t32) * t26 + (g(1) * (-t15 - t32) + g(2) * t40) * t23) - m(5) * (g(1) * (-t2 * rSges(5,1) + t1 * rSges(5,2) - t43) + g(2) * (rSges(5,1) * t4 - rSges(5,2) * t3 + rSges(5,3) * t48 + t38) + (g(1) * (-t16 * rSges(5,3) + t31) - t55) * t23) - m(6) * (g(1) * (-t39 * t1 - t52 * t2 - t43) + g(2) * (rSges(6,2) * t48 + t3 * t39 + t4 * t52 + t38) + (g(1) * (-t16 * rSges(6,2) + t31) - t55) * t23), -m(3) * (g(3) * t34 + t59 * (-rSges(3,1) * t22 - rSges(3,2) * t25)) - m(4) * (g(3) * (t18 + t32) + t59 * (-rSges(4,1) * t16 - rSges(4,2) * t17 - t57)) - m(5) * (g(1) * (rSges(5,3) * t46 + t35) + g(2) * (rSges(5,3) * t47 + t36) + g(3) * (t17 * t60 + t37) + (g(3) * rSges(5,3) + t59 * (-pkin(3) - t60)) * t16) - m(6) * (g(1) * t35 + g(2) * t36 + g(3) * t37 + (rSges(6,2) * t59 + g(3) * t58) * t17 + (g(3) * rSges(6,2) + t59 * (-pkin(3) - t58)) * t16), (-m(4) - m(5) - m(6)) * (g(1) * t23 - g(2) * t26), -m(5) * (g(1) * (-rSges(5,1) * t3 - rSges(5,2) * t4) + g(2) * (-rSges(5,1) * t1 - rSges(5,2) * t2)) - m(6) * (g(1) * (-t3 * t52 + t39 * t4) + g(2) * (-t1 * t52 + t2 * t39)) + (-m(5) * (-rSges(5,1) * t21 - rSges(5,2) * t24) - m(6) * (-t21 * t52 + t24 * t39)) * t53, -m(6) * (g(1) * t3 + g(2) * t1 + t21 * t53)];
taug = t5(:);
