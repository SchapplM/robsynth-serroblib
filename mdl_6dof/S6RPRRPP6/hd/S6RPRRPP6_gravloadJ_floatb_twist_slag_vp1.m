% Calculate Gravitation load on the joints for
% S6RPRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
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
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPP6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP6_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:46:22
% EndTime: 2019-03-09 04:46:24
% DurationCPUTime: 0.85s
% Computational Cost: add. (316->131), mult. (471->173), div. (0->0), fcn. (451->8), ass. (0->60)
t32 = sin(qJ(3));
t35 = cos(qJ(3));
t66 = rSges(5,3) + pkin(8);
t78 = t66 * t35;
t81 = t32 * pkin(3) - t78;
t34 = cos(qJ(4));
t22 = pkin(4) * t34 + pkin(3);
t29 = qJ(4) + pkin(9);
t23 = sin(t29);
t24 = cos(t29);
t51 = rSges(7,3) + qJ(6);
t67 = rSges(7,1) + pkin(5);
t77 = t51 * t23 + t67 * t24;
t80 = -t22 - t77;
t36 = cos(qJ(1));
t58 = t34 * t36;
t31 = sin(qJ(4));
t33 = sin(qJ(1));
t61 = t33 * t31;
t6 = -t32 * t61 + t58;
t60 = t33 * t34;
t64 = t32 * t36;
t8 = t31 * t64 + t60;
t70 = g(2) * t36;
t72 = g(1) * t33;
t45 = -t70 + t72;
t75 = -m(6) - m(7);
t74 = -pkin(1) - pkin(7);
t73 = pkin(4) * t31;
t71 = g(2) * t35;
t69 = g(3) * t35;
t65 = t31 * t36;
t63 = t33 * t23;
t62 = t33 * t24;
t59 = t33 * t35;
t57 = t35 * t36;
t56 = t36 * t24;
t30 = -qJ(5) - pkin(8);
t55 = rSges(7,2) - t30;
t54 = rSges(6,3) - t30;
t53 = t8 * pkin(4);
t52 = t36 * pkin(1) + t33 * qJ(2);
t26 = t36 * qJ(2);
t48 = t22 * t64 + t30 * t57 + t26;
t47 = t36 * pkin(7) + t52;
t46 = g(1) * t22 * t59 + g(2) * t30 * t64;
t43 = rSges(4,1) * t32 + rSges(4,2) * t35;
t42 = rSges(6,1) * t24 - rSges(6,2) * t23;
t41 = t6 * pkin(4);
t40 = g(1) * (-t73 + t74);
t39 = t33 * t32 * t22 + pkin(4) * t65 + t30 * t59 + t47;
t38 = rSges(5,1) * t34 - rSges(5,2) * t31 + pkin(3);
t37 = -t22 - t42;
t9 = t32 * t58 - t61;
t7 = t32 * t60 + t65;
t4 = t32 * t56 - t63;
t3 = t23 * t64 + t62;
t2 = t23 * t36 + t32 * t62;
t1 = t32 * t63 - t56;
t5 = [-m(2) * (g(1) * (-t33 * rSges(2,1) - rSges(2,2) * t36) + g(2) * (rSges(2,1) * t36 - t33 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t36 + t26 + (rSges(3,2) - pkin(1)) * t33) + g(2) * (-rSges(3,2) * t36 + t33 * rSges(3,3) + t52)) - m(4) * (g(1) * (rSges(4,1) * t64 + rSges(4,2) * t57 + t26) + g(2) * (rSges(4,3) * t36 + t47) + (g(1) * (-rSges(4,3) + t74) + g(2) * t43) * t33) - m(5) * ((rSges(5,1) * t7 + rSges(5,2) * t6 + t33 * t81 + t47) * g(2) + (t9 * rSges(5,1) - t8 * rSges(5,2) + t74 * t33 + t36 * t81 + t26) * g(1)) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) - rSges(6,3) * t57 + t48) + g(2) * (rSges(6,1) * t2 - rSges(6,2) * t1 + t39) + (-rSges(6,3) * t71 + t40) * t33) - m(7) * (g(1) * (-rSges(7,2) * t57 + t51 * t3 + t67 * t4 + t48) + g(2) * (t51 * t1 + t67 * t2 + t39) + (-rSges(7,2) * t71 + t40) * t33) (-m(3) - m(4) - m(5) + t75) * t45, -m(4) * (-g(3) * t43 + t45 * (rSges(4,1) * t35 - rSges(4,2) * t32)) - m(5) * (g(3) * (-t38 * t32 + t78) + t45 * (t66 * t32 + t38 * t35)) - m(6) * ((-rSges(6,3) * t70 + g(3) * t37 + t54 * t72) * t32 + (g(3) * t54 + t37 * t70 + t42 * t72) * t35 + t46) - m(7) * ((-rSges(7,2) * t70 + g(3) * t80 + t55 * t72) * t32 + (g(3) * t55 + t80 * t70 + t77 * t72) * t35 + t46) -m(5) * (g(1) * (rSges(5,1) * t6 - rSges(5,2) * t7) + g(2) * (rSges(5,1) * t8 + rSges(5,2) * t9)) - m(6) * (g(1) * (-rSges(6,1) * t1 - rSges(6,2) * t2 + t41) + g(2) * (rSges(6,1) * t3 + rSges(6,2) * t4 + t53)) - m(7) * (g(1) * (-t67 * t1 + t51 * t2 + t41) + g(2) * (t67 * t3 - t51 * t4 + t53)) + (-m(5) * (-rSges(5,1) * t31 - rSges(5,2) * t34) - m(6) * (-rSges(6,1) * t23 - rSges(6,2) * t24 - t73) - m(7) * (-t67 * t23 + t51 * t24 - t73)) * t69, t75 * (g(3) * t32 - t45 * t35) -m(7) * (g(1) * t1 - g(2) * t3 + t23 * t69)];
taug  = t5(:);
