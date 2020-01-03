% Calculate Gravitation load on the joints for
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPP6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP6_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:59:59
% EndTime: 2019-12-31 21:00:02
% DurationCPUTime: 0.75s
% Computational Cost: add. (291->115), mult. (437->159), div. (0->0), fcn. (425->8), ass. (0->52)
t26 = sin(qJ(3));
t28 = sin(qJ(1));
t29 = cos(qJ(3));
t30 = cos(qJ(2));
t31 = cos(qJ(1));
t51 = t30 * t31;
t8 = -t26 * t51 + t28 * t29;
t58 = rSges(6,1) + pkin(4);
t62 = g(2) * t28;
t68 = g(1) * t31 + t62;
t46 = rSges(6,3) + qJ(5);
t27 = sin(qJ(2));
t57 = rSges(4,3) + pkin(7);
t67 = t30 * pkin(2) + t57 * t27;
t24 = qJ(3) + pkin(8);
t19 = sin(t24);
t20 = cos(t24);
t66 = t46 * t19 + t58 * t20;
t65 = pkin(3) * t26;
t64 = g(1) * t28;
t18 = pkin(3) * t29 + pkin(2);
t11 = t30 * t18;
t61 = g(3) * t11;
t60 = g(3) * t27;
t56 = rSges(3,2) * t27;
t55 = t26 * t31;
t54 = t28 * t26;
t52 = t28 * t30;
t50 = t31 * t19;
t25 = -qJ(4) - pkin(7);
t49 = rSges(6,2) - t25;
t48 = rSges(5,3) - t25;
t47 = t31 * pkin(1) + t28 * pkin(6);
t22 = t31 * pkin(6);
t44 = t28 * t27 * t25 + pkin(3) * t55 + t22;
t43 = -pkin(1) - t11;
t42 = t49 * t31;
t41 = t48 * t31;
t40 = pkin(3) * t54 + t18 * t51 + t47;
t39 = rSges(3,1) * t30 - t56;
t37 = rSges(5,1) * t20 - rSges(5,2) * t19;
t36 = t8 * pkin(3);
t35 = rSges(4,1) * t29 - rSges(4,2) * t26 + pkin(2);
t6 = t26 * t52 + t29 * t31;
t33 = t6 * pkin(3);
t9 = t29 * t51 + t54;
t7 = -t29 * t52 + t55;
t4 = t28 * t19 + t20 * t51;
t3 = -t28 * t20 + t30 * t50;
t2 = t20 * t52 - t50;
t1 = t19 * t52 + t20 * t31;
t5 = [-m(2) * (g(1) * (-t28 * rSges(2,1) - rSges(2,2) * t31) + g(2) * (rSges(2,1) * t31 - t28 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t31 + t22) + g(2) * (rSges(3,1) * t51 - t31 * t56 + t47) + (g(1) * (-pkin(1) - t39) + g(2) * rSges(3,3)) * t28) - m(4) * (g(1) * (rSges(4,1) * t7 + rSges(4,2) * t6 + t22) + (-pkin(1) - t67) * t64 + (t9 * rSges(4,1) + t8 * rSges(4,2) + t67 * t31 + t47) * g(2)) - m(5) * (g(1) * (-rSges(5,1) * t2 + rSges(5,2) * t1 + t44) + g(2) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t27 * t41 + t40) + (-rSges(5,3) * t27 + t43) * t64) - m(6) * (g(1) * (-t46 * t1 - t58 * t2 + t44) + g(2) * (t27 * t42 + t46 * t3 + t58 * t4 + t40) + (-rSges(6,2) * t27 + t43) * t64), -m(3) * (g(3) * t39 + t68 * (-rSges(3,1) * t27 - rSges(3,2) * t30)) - m(4) * ((g(3) * t35 + t68 * t57) * t30 + (g(3) * t57 - t68 * t35) * t27) - m(5) * (t61 + (g(1) * t41 + g(3) * t37 + t48 * t62) * t30 + (g(3) * t48 + t68 * (-t18 - t37)) * t27) - m(6) * (t61 + (g(1) * t42 + g(3) * t66 + t49 * t62) * t30 + (g(3) * t49 + t68 * (-t18 - t66)) * t27), -m(4) * (g(1) * (rSges(4,1) * t8 - rSges(4,2) * t9) + g(2) * (-rSges(4,1) * t6 + rSges(4,2) * t7)) - m(5) * (g(1) * (-t3 * rSges(5,1) - t4 * rSges(5,2) + t36) + g(2) * (-t1 * rSges(5,1) - t2 * rSges(5,2) - t33)) - m(6) * (g(1) * (-t58 * t3 + t46 * t4 + t36) + g(2) * (-t58 * t1 + t46 * t2 - t33)) + (-m(4) * (-rSges(4,1) * t26 - rSges(4,2) * t29) - m(5) * (-rSges(5,1) * t19 - rSges(5,2) * t20 - t65) - m(6) * (-t58 * t19 + t46 * t20 - t65)) * t60, (-m(5) - m(6)) * (-g(3) * t30 + t68 * t27), -m(6) * (g(1) * t3 + g(2) * t1 + t19 * t60)];
taug = t5(:);
