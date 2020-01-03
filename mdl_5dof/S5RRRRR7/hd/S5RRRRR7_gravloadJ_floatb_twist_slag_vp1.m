% Calculate Gravitation load on the joints for
% S5RRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR7_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:21:22
% EndTime: 2019-12-31 22:21:24
% DurationCPUTime: 0.51s
% Computational Cost: add. (378->96), mult. (347->126), div. (0->0), fcn. (300->10), ass. (0->59)
t84 = rSges(6,3) + pkin(9);
t31 = qJ(2) + qJ(3);
t28 = qJ(4) + t31;
t23 = sin(t28);
t24 = cos(t28);
t32 = sin(qJ(5));
t72 = rSges(6,2) * t32;
t83 = t23 * t72 + t24 * t84;
t81 = t24 * rSges(5,1) - t23 * rSges(5,2);
t26 = sin(t31);
t27 = cos(t31);
t57 = t27 * rSges(4,1) - t26 * rSges(4,2);
t34 = sin(qJ(1));
t37 = cos(qJ(1));
t80 = g(1) * t37 + g(2) * t34;
t45 = t24 * pkin(4) + t84 * t23;
t35 = cos(qJ(5));
t73 = rSges(6,1) * t35;
t79 = (-pkin(4) - t73) * t23;
t38 = -pkin(7) - pkin(6);
t78 = pkin(3) * t26;
t33 = sin(qJ(2));
t75 = t33 * pkin(2);
t74 = rSges(3,3) + pkin(6);
t36 = cos(qJ(2));
t29 = t36 * pkin(2);
t25 = t29 + pkin(1);
t67 = t34 * t32;
t66 = t34 * t35;
t65 = t37 * t32;
t64 = t37 * t35;
t63 = rSges(4,3) - t38;
t30 = -pkin(8) + t38;
t62 = rSges(5,3) - t30;
t61 = t83 * t34;
t60 = t83 * t37;
t22 = pkin(3) * t27;
t55 = t22 + t81;
t54 = t36 * rSges(3,1) - t33 * rSges(3,2);
t52 = -rSges(4,1) * t26 - rSges(4,2) * t27;
t50 = -rSges(5,1) * t23 - rSges(5,2) * t24;
t49 = pkin(1) + t54;
t48 = t25 + t57;
t47 = t50 * t34;
t46 = t50 * t37;
t44 = t45 + (-t72 + t73) * t24;
t42 = t22 + t44;
t41 = g(1) * t60 + g(2) * t61;
t40 = t80 * t79;
t11 = -t75 - t78;
t10 = t22 + t25;
t7 = t37 * t11;
t6 = t34 * t11;
t5 = t37 * t10;
t4 = t24 * t64 + t67;
t3 = -t24 * t65 + t66;
t2 = -t24 * t66 + t65;
t1 = t24 * t67 + t64;
t8 = [-m(2) * (g(1) * (-t34 * rSges(2,1) - t37 * rSges(2,2)) + g(2) * (t37 * rSges(2,1) - t34 * rSges(2,2))) - m(3) * ((g(1) * t74 + g(2) * t49) * t37 + (-g(1) * t49 + g(2) * t74) * t34) - m(4) * ((g(1) * t63 + g(2) * t48) * t37 + (-g(1) * t48 + g(2) * t63) * t34) - m(5) * (g(2) * t5 + (g(1) * t62 + g(2) * t81) * t37 + (g(1) * (-t10 - t81) + g(2) * t62) * t34) - m(6) * (g(1) * (t2 * rSges(6,1) + t1 * rSges(6,2)) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t5) + (-g(1) * t30 + g(2) * t45) * t37 + (g(1) * (-t10 - t45) - g(2) * t30) * t34), -m(3) * (g(3) * t54 + t80 * (-rSges(3,1) * t33 - rSges(3,2) * t36)) - m(4) * (g(3) * (t29 + t57) + t80 * (t52 - t75)) - m(5) * (g(1) * (t7 + t46) + g(2) * (t6 + t47) + g(3) * (t29 + t55)) - m(6) * (g(1) * (t7 + t60) + g(2) * (t6 + t61) + g(3) * (t29 + t42) + t40), -m(6) * t41 + (-m(4) * t57 - m(5) * t55 - m(6) * t42) * g(3) + t80 * (-m(4) * t52 - m(5) * (t50 - t78) - m(6) * (-t78 + t79)), -m(5) * (g(1) * t46 + g(2) * t47 + g(3) * t81) - m(6) * (g(3) * t44 + t40 + t41), -m(6) * (g(1) * (t3 * rSges(6,1) - t4 * rSges(6,2)) + g(2) * (-t1 * rSges(6,1) + t2 * rSges(6,2)) + g(3) * (-rSges(6,1) * t32 - rSges(6,2) * t35) * t23)];
taug = t8(:);
