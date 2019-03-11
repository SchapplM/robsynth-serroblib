% Calculate Gravitation load on the joints for
% S6RRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 11:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:47:21
% EndTime: 2019-03-09 11:47:23
% DurationCPUTime: 0.73s
% Computational Cost: add. (487->144), mult. (511->186), div. (0->0), fcn. (480->10), ass. (0->70)
t32 = qJ(2) + pkin(10);
t25 = sin(t32);
t72 = rSges(5,3) + pkin(8);
t87 = t72 * t25;
t41 = -pkin(9) - pkin(8);
t60 = rSges(7,3) + qJ(6) - t41;
t86 = t60 * t25;
t61 = rSges(6,3) - t41;
t85 = t61 * t25;
t37 = sin(qJ(1));
t40 = cos(qJ(1));
t84 = g(1) * t40 + g(2) * t37;
t26 = cos(t32);
t33 = qJ(4) + qJ(5);
t28 = cos(t33);
t68 = t28 * t37;
t27 = sin(t33);
t69 = t27 * t40;
t10 = -t26 * t68 + t69;
t67 = t28 * t40;
t70 = t27 * t37;
t9 = t26 * t70 + t67;
t83 = -t9 * rSges(6,1) + t10 * rSges(6,2);
t82 = -t9 * rSges(7,1) + t10 * rSges(7,2);
t11 = -t26 * t69 + t68;
t12 = t26 * t67 + t70;
t81 = t11 * rSges(7,1) - t12 * rSges(7,2);
t80 = t11 * rSges(6,1) - t12 * rSges(6,2);
t36 = sin(qJ(2));
t79 = pkin(2) * t36;
t35 = sin(qJ(4));
t78 = pkin(4) * t35;
t77 = pkin(5) * t27;
t74 = g(3) * t25;
t73 = rSges(3,3) + pkin(7);
t18 = t77 + t78;
t71 = t18 * t26;
t66 = t35 * t37;
t65 = t35 * t40;
t38 = cos(qJ(4));
t64 = t37 * t38;
t63 = t38 * t40;
t34 = -qJ(3) - pkin(7);
t62 = rSges(4,3) - t34;
t59 = t18 - t34;
t29 = t38 * pkin(4);
t19 = pkin(5) * t28 + t29;
t58 = -t34 + t78;
t39 = cos(qJ(2));
t57 = rSges(3,1) * t39 - rSges(3,2) * t36;
t55 = rSges(4,1) * t26 - rSges(4,2) * t25;
t54 = -rSges(6,1) * t27 - rSges(6,2) * t28;
t53 = -rSges(7,1) * t27 - rSges(7,2) * t28;
t52 = pkin(1) + t57;
t51 = rSges(5,1) * t38 - rSges(5,2) * t35 + pkin(3);
t15 = -t26 * t65 + t64;
t13 = t26 * t66 + t63;
t23 = t29 + pkin(3);
t50 = rSges(6,1) * t28 - rSges(6,2) * t27 + t23;
t17 = pkin(3) + t19;
t49 = rSges(7,1) * t28 - rSges(7,2) * t27 + t17;
t48 = t26 * pkin(3) + t87;
t46 = t26 * t23 + t85;
t45 = t26 * t17 + t86;
t30 = t39 * pkin(2);
t24 = t30 + pkin(1);
t21 = t40 * t24;
t16 = t26 * t63 + t66;
t14 = -t26 * t64 + t65;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t37 - rSges(2,2) * t40) + g(2) * (rSges(2,1) * t40 - rSges(2,2) * t37)) - m(3) * ((g(1) * t73 + g(2) * t52) * t40 + (-g(1) * t52 + g(2) * t73) * t37) - m(4) * (g(2) * t21 + (g(1) * t62 + g(2) * t55) * t40 + (g(1) * (-t24 - t55) + g(2) * t62) * t37) - m(5) * (g(1) * (rSges(5,1) * t14 + rSges(5,2) * t13) + g(2) * (rSges(5,1) * t16 + rSges(5,2) * t15 + t21) + (-g(1) * t34 + g(2) * t48) * t40 + (g(1) * (-t24 - t48) - g(2) * t34) * t37) - m(6) * (g(1) * (t10 * rSges(6,1) + t9 * rSges(6,2)) + g(2) * (t12 * rSges(6,1) + t11 * rSges(6,2) + t21) + (g(1) * t58 + g(2) * t46) * t40 + (g(1) * (-t24 - t46) + g(2) * t58) * t37) - m(7) * (g(1) * (rSges(7,1) * t10 + rSges(7,2) * t9) + g(2) * (rSges(7,1) * t12 + rSges(7,2) * t11 + t21) + (g(1) * t59 + g(2) * t45) * t40 + (g(1) * (-t24 - t45) + g(2) * t59) * t37) -m(3) * (g(3) * t57 + t84 * (-rSges(3,1) * t36 - rSges(3,2) * t39)) - m(4) * (g(3) * (t30 + t55) + t84 * (-rSges(4,1) * t25 - rSges(4,2) * t26 - t79)) - m(5) * (g(3) * (t51 * t26 + t30 + t87) + t84 * (-t51 * t25 + t72 * t26 - t79)) - m(6) * (g(3) * (t50 * t26 + t30 + t85) + t84 * (-t50 * t25 + t61 * t26 - t79)) - m(7) * (g(3) * (t49 * t26 + t30 + t86) + t84 * (-t49 * t25 + t60 * t26 - t79)) (-m(4) - m(5) - m(6) - m(7)) * (g(1) * t37 - g(2) * t40) -m(5) * (g(1) * (rSges(5,1) * t15 - rSges(5,2) * t16) + g(2) * (-rSges(5,1) * t13 + rSges(5,2) * t14)) - m(6) * (g(1) * (t15 * pkin(4) + t80) + g(2) * (-t13 * pkin(4) + t83)) - m(7) * (g(1) * (t19 * t37 - t40 * t71 + t81) + g(2) * (-t19 * t40 - t37 * t71 + t82)) + (-m(5) * (-rSges(5,1) * t35 - rSges(5,2) * t38) - m(6) * (t54 - t78) - m(7) * (-t18 + t53)) * t74, -m(6) * (g(1) * t80 + g(2) * t83) - m(7) * (g(1) * (t11 * pkin(5) + t81) + g(2) * (-t9 * pkin(5) + t82)) + (-m(6) * t54 - m(7) * (t53 - t77)) * t74, -m(7) * (-g(3) * t26 + t84 * t25)];
taug  = t1(:);
