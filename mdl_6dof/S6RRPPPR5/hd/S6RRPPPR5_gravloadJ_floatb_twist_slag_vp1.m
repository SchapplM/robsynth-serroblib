% Calculate Gravitation load on the joints for
% S6RRPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
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
% Datum: 2019-03-09 08:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPPR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPPR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:21:21
% EndTime: 2019-03-09 08:21:23
% DurationCPUTime: 0.78s
% Computational Cost: add. (296->156), mult. (660->197), div. (0->0), fcn. (681->8), ass. (0->57)
t37 = cos(qJ(1));
t34 = sin(qJ(1));
t76 = g(2) * t34;
t80 = g(1) * t37 + t76;
t79 = -m(6) - m(7);
t78 = g(1) * t34;
t33 = sin(qJ(2));
t75 = g(3) * t33;
t36 = cos(qJ(2));
t27 = t36 * pkin(2);
t74 = pkin(8) + rSges(7,3);
t30 = sin(pkin(9));
t73 = t30 * t36;
t31 = cos(pkin(9));
t72 = t31 * t36;
t71 = t33 * t37;
t70 = t34 * t36;
t69 = t36 * t37;
t68 = t37 * t30;
t67 = -pkin(4) - qJ(3);
t66 = pkin(5) + qJ(4);
t61 = qJ(3) * t36;
t14 = t34 * t61;
t65 = pkin(4) * t70 + t14;
t17 = t37 * t61;
t64 = pkin(4) * t69 + t17;
t24 = t33 * qJ(3);
t63 = t24 + t27;
t62 = t37 * pkin(1) + t34 * pkin(7);
t60 = rSges(6,1) + qJ(4);
t59 = rSges(5,3) + qJ(4);
t58 = rSges(6,3) + qJ(5);
t57 = -m(5) + t79;
t56 = -pkin(1) - t27;
t55 = t74 * t37;
t54 = pkin(3) * t72 + qJ(4) * t73 + t63;
t53 = pkin(2) * t69 + t37 * t24 + t62;
t28 = t37 * pkin(7);
t8 = t30 * t70 + t31 * t37;
t9 = t31 * t70 - t68;
t52 = -t9 * pkin(3) - t8 * qJ(4) + t28;
t11 = t34 * t30 + t31 * t69;
t51 = t11 * pkin(3) + t53;
t32 = sin(qJ(6));
t35 = cos(qJ(6));
t50 = -t32 * t9 - t35 * t8;
t49 = t32 * t8 - t35 * t9;
t48 = rSges(3,1) * t36 - rSges(3,2) * t33;
t46 = t30 * t35 + t31 * t32;
t45 = -t30 * t32 + t31 * t35;
t44 = pkin(4) * t71 + t51;
t43 = t33 * pkin(4) + qJ(5) * t72 + t54;
t41 = -t9 * qJ(5) + t52;
t10 = -t34 * t31 + t36 * t68;
t4 = t10 * t35 + t11 * t32;
t3 = -t10 * t32 + t11 * t35;
t1 = [-m(2) * (g(1) * (-t34 * rSges(2,1) - rSges(2,2) * t37) + g(2) * (rSges(2,1) * t37 - t34 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t37 + t28) + g(2) * (rSges(3,1) * t69 - rSges(3,2) * t71 + t62) + (g(1) * (-pkin(1) - t48) + g(2) * rSges(3,3)) * t34) - m(4) * (g(1) * (-rSges(4,1) * t9 + rSges(4,2) * t8 + t28) + g(2) * (t11 * rSges(4,1) - t10 * rSges(4,2) + rSges(4,3) * t71 + t53) + ((-rSges(4,3) - qJ(3)) * t33 + t56) * t78) - m(5) * (g(1) * (rSges(5,2) * t9 - rSges(5,3) * t8 + t52) + g(2) * (rSges(5,1) * t71 - t11 * rSges(5,2) + t59 * t10 + t51) + ((-rSges(5,1) - qJ(3)) * t33 + t56) * t78) - m(6) * (g(1) * (-rSges(6,1) * t8 - rSges(6,3) * t9 + t41) + g(2) * (-rSges(6,2) * t71 + t60 * t10 + t58 * t11 + t44) + ((rSges(6,2) + t67) * t33 + t56) * t78) - m(7) * (g(1) * (t50 * rSges(7,1) + t49 * rSges(7,2) - t8 * pkin(5) + t41) + g(2) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t11 * qJ(5) + t66 * t10 + t33 * t55 + t44) + ((t67 - t74) * t33 + t56) * t78) -m(3) * (g(3) * t48 + t80 * (-rSges(3,1) * t33 - rSges(3,2) * t36)) - m(4) * (g(1) * (rSges(4,3) * t69 + t17) + g(2) * (rSges(4,3) * t70 + t14) + g(3) * (rSges(4,1) * t72 - rSges(4,2) * t73 + t63) + (g(3) * rSges(4,3) + t80 * (-rSges(4,1) * t31 + rSges(4,2) * t30 - pkin(2))) * t33) - m(5) * (g(1) * (rSges(5,1) * t69 + t17) + g(2) * (rSges(5,1) * t70 + t14) + g(3) * (-rSges(5,2) * t72 + rSges(5,3) * t73 + t54) + (g(3) * rSges(5,1) + t80 * (-pkin(2) + (rSges(5,2) - pkin(3)) * t31 - t59 * t30)) * t33) - m(6) * (g(1) * (-rSges(6,2) * t69 + t64) + g(2) * (-rSges(6,2) * t70 + t65) + g(3) * (rSges(6,1) * t73 + rSges(6,3) * t72 + t43) + (-g(3) * rSges(6,2) + t80 * (-pkin(2) - t60 * t30 + (-pkin(3) - t58) * t31)) * t33) - m(7) * (g(1) * t64 + g(2) * t65 + g(3) * t43 + (g(3) * (t46 * rSges(7,1) + t45 * rSges(7,2) + t30 * pkin(5)) + g(1) * t55 + t74 * t76) * t36 + (g(3) * t74 + t80 * (-pkin(2) + (-t32 * rSges(7,1) - t35 * rSges(7,2) - pkin(3) - qJ(5)) * t31 + (-t35 * rSges(7,1) + t32 * rSges(7,2) - t66) * t30)) * t33) (-m(4) + t57) * (-g(3) * t36 + t80 * t33) t57 * (g(1) * t10 + g(2) * t8 + t30 * t75) t79 * (g(1) * t11 + g(2) * t9 + t31 * t75) -m(7) * (g(1) * (rSges(7,1) * t3 - rSges(7,2) * t4) + g(2) * (-t49 * rSges(7,1) + t50 * rSges(7,2)) + (t45 * rSges(7,1) - t46 * rSges(7,2)) * t75)];
taug  = t1(:);
