% Calculate Gravitation load on the joints for
% S5RRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP7_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:56:08
% EndTime: 2019-12-31 21:56:10
% DurationCPUTime: 0.70s
% Computational Cost: add. (343->119), mult. (441->161), div. (0->0), fcn. (415->8), ass. (0->60)
t36 = cos(qJ(4));
t80 = rSges(6,1) + pkin(4);
t82 = -t80 * t36 - pkin(3);
t35 = sin(qJ(1));
t32 = qJ(2) + qJ(3);
t29 = sin(t32);
t33 = sin(qJ(4));
t69 = t29 * t33;
t53 = rSges(5,2) * t69;
t30 = cos(t32);
t66 = t30 * t35;
t81 = rSges(5,3) * t66 + t35 * t53;
t38 = cos(qJ(1));
t64 = t30 * t38;
t79 = rSges(5,3) * t64 + t38 * t53;
t23 = t30 * rSges(4,1);
t47 = -rSges(4,2) * t29 + t23;
t57 = t30 * pkin(3) + t29 * pkin(8);
t78 = g(1) * t38 + g(2) * t35;
t55 = rSges(6,3) + qJ(5);
t77 = t78 * t29;
t34 = sin(qJ(2));
t76 = pkin(2) * t34;
t39 = -pkin(7) - pkin(6);
t73 = g(2) * t39;
t71 = rSges(3,3) + pkin(6);
t22 = t29 * rSges(6,2);
t21 = t29 * rSges(5,3);
t68 = t29 * t38;
t67 = t30 * t33;
t65 = t30 * t36;
t63 = t33 * t35;
t62 = t35 * t36;
t61 = t36 * t38;
t60 = t38 * t33;
t59 = t38 * t39;
t58 = rSges(4,3) - t39;
t56 = qJ(5) * t33;
t37 = cos(qJ(2));
t31 = t37 * pkin(2);
t28 = t31 + pkin(1);
t8 = t38 * t28;
t54 = pkin(3) * t64 + pkin(8) * t68 + t8;
t52 = -rSges(5,1) * t36 - pkin(3);
t16 = pkin(8) * t66;
t51 = -t35 * t76 + t16;
t19 = pkin(8) * t64;
t50 = -t38 * t76 + t19;
t49 = rSges(3,1) * t37 - rSges(3,2) * t34;
t46 = -t28 - t57;
t45 = pkin(1) + t49;
t44 = rSges(6,3) * t67 + t30 * t56 + t65 * t80 + t22 + t57;
t43 = rSges(5,1) * t65 - rSges(5,2) * t67 + t21 + t57;
t15 = rSges(6,2) * t64;
t10 = rSges(6,2) * t66;
t4 = t30 * t61 + t63;
t3 = t30 * t60 - t62;
t2 = t30 * t62 - t60;
t1 = t30 * t63 + t61;
t5 = [-m(2) * (g(1) * (-rSges(2,1) * t35 - rSges(2,2) * t38) + g(2) * (rSges(2,1) * t38 - rSges(2,2) * t35)) - m(3) * ((g(1) * t71 + g(2) * t45) * t38 + (-g(1) * t45 + g(2) * t71) * t35) - m(4) * (g(2) * t8 + (g(1) * t58 + g(2) * t47) * t38 + (g(1) * (-t28 - t47) + g(2) * t58) * t35) - m(5) * (g(1) * (-t2 * rSges(5,1) + t1 * rSges(5,2) - t59) + g(2) * (t4 * rSges(5,1) - t3 * rSges(5,2) + rSges(5,3) * t68 + t54) + (g(1) * (t46 - t21) - t73) * t35) - m(6) * (g(1) * (-t55 * t1 - t80 * t2 - t59) + g(2) * (rSges(6,2) * t68 + t55 * t3 + t80 * t4 + t54) + (g(1) * (t46 - t22) - t73) * t35), -m(3) * (g(3) * t49 + t78 * (-rSges(3,1) * t34 - rSges(3,2) * t37)) - m(4) * (g(3) * (t31 + t47) + t78 * (-rSges(4,1) * t29 - rSges(4,2) * t30 - t76)) - m(5) * (g(1) * (t50 + t79) + g(2) * (t51 + t81) + g(3) * (t31 + t43) + t52 * t77) - m(6) * (g(1) * (t15 + t50) + g(2) * (t10 + t51) + g(3) * (t31 + t44) + (-t55 * t33 + t82) * t77), -m(4) * (g(3) * t23 + (-g(1) * t64 - g(2) * t66) * rSges(4,2)) - m(5) * (g(1) * (t19 + t79) + g(2) * (t16 + t81) + g(3) * t43) - m(6) * (g(1) * (t15 + t19) + g(2) * (t10 + t16) + g(3) * t44) + (m(4) * g(3) * rSges(4,2) + t78 * (m(4) * rSges(4,1) - m(5) * t52 - m(6) * (-rSges(6,3) * t33 - t56 + t82))) * t29, -m(5) * (g(1) * (-rSges(5,1) * t3 - rSges(5,2) * t4) + g(2) * (-rSges(5,1) * t1 - rSges(5,2) * t2)) - m(6) * (g(1) * (-t3 * t80 + t55 * t4) + g(2) * (-t1 * t80 + t55 * t2)) + (-m(5) * (-rSges(5,1) * t33 - rSges(5,2) * t36) - m(6) * (-t80 * t33 + t55 * t36)) * g(3) * t29, -m(6) * (g(1) * t3 + g(2) * t1 + g(3) * t69)];
taug = t5(:);
