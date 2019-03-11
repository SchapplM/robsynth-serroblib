% Calculate Gravitation load on the joints for
% S6RPRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% Datum: 2019-03-09 04:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPP7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP7_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:50:12
% EndTime: 2019-03-09 04:50:14
% DurationCPUTime: 0.72s
% Computational Cost: add. (248->140), mult. (521->183), div. (0->0), fcn. (515->6), ass. (0->56)
t25 = sin(qJ(4));
t28 = cos(qJ(4));
t48 = rSges(6,3) + qJ(5);
t64 = -rSges(6,1) - pkin(4);
t73 = -t48 * t25 + t64 * t28 - pkin(3);
t70 = -pkin(1) - pkin(7);
t63 = rSges(7,1) + pkin(5);
t30 = cos(qJ(1));
t67 = g(2) * t30;
t27 = sin(qJ(1));
t69 = g(1) * t27;
t37 = -t67 + t69;
t71 = -m(6) - m(7);
t29 = cos(qJ(3));
t68 = g(2) * t29;
t55 = t28 * t29;
t66 = g(3) * qJ(5) * t55;
t65 = g(3) * t29;
t62 = -rSges(6,2) - pkin(8);
t61 = -rSges(5,3) - pkin(8);
t26 = sin(qJ(3));
t60 = t26 * t27;
t59 = t26 * t30;
t58 = t27 * t25;
t57 = t27 * t28;
t56 = t27 * t29;
t54 = t29 * t30;
t53 = t30 * t28;
t52 = pkin(3) * t56 + pkin(8) * t60;
t21 = t30 * qJ(2);
t51 = pkin(3) * t59 + t21;
t50 = t30 * pkin(1) + t27 * qJ(2);
t49 = rSges(7,2) + qJ(5);
t47 = -rSges(7,3) - qJ(6);
t46 = -pkin(4) - t63;
t45 = pkin(8) * t54;
t44 = t25 * t56;
t43 = t27 * t55;
t42 = -pkin(8) - t47;
t41 = t30 * pkin(7) + t50;
t40 = g(1) * t70;
t39 = pkin(4) * t43 + qJ(5) * t44 + t52;
t38 = pkin(3) * t60 + t41;
t35 = rSges(4,1) * t26 + rSges(4,2) * t29;
t8 = t25 * t59 + t57;
t9 = t26 * t53 - t58;
t34 = t9 * pkin(4) + t8 * qJ(5) + t51;
t33 = -rSges(5,1) * t28 + rSges(5,2) * t25 - pkin(3);
t6 = t26 * t58 - t53;
t7 = t25 * t30 + t26 * t57;
t32 = t7 * pkin(4) + qJ(5) * t6 + t38;
t31 = -t49 * t25 + t46 * t28 - pkin(3);
t22 = t29 * pkin(8);
t4 = t8 * pkin(4);
t2 = t6 * pkin(4);
t1 = [-m(2) * (g(1) * (-t27 * rSges(2,1) - rSges(2,2) * t30) + g(2) * (rSges(2,1) * t30 - t27 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t30 + t21 + (rSges(3,2) - pkin(1)) * t27) + g(2) * (-rSges(3,2) * t30 + t27 * rSges(3,3) + t50)) - m(4) * (g(1) * (rSges(4,1) * t59 + rSges(4,2) * t54 + t21) + g(2) * (rSges(4,3) * t30 + t41) + (g(1) * (-rSges(4,3) + t70) + g(2) * t35) * t27) - m(5) * (g(1) * (t9 * rSges(5,1) - t8 * rSges(5,2) - rSges(5,3) * t54 - t45 + t51) + g(2) * (rSges(5,1) * t7 - rSges(5,2) * t6 + t38) + (t61 * t68 + t40) * t27) - m(6) * (g(1) * (t9 * rSges(6,1) - rSges(6,2) * t54 + t8 * rSges(6,3) + t34 - t45) + g(2) * (rSges(6,1) * t7 + rSges(6,3) * t6 + t32) + (t62 * t68 + t40) * t27) - m(7) * (g(1) * (t8 * rSges(7,2) + t70 * t27 + t63 * t9 + t34) + g(2) * (rSges(7,2) * t6 + t63 * t7 + t32) + (g(1) * t30 + g(2) * t27) * t29 * t42) (-m(3) - m(4) - m(5) + t71) * t37, -m(4) * (-g(3) * t35 + t37 * (rSges(4,1) * t29 - rSges(4,2) * t26)) - m(5) * (g(1) * (rSges(5,1) * t43 - rSges(5,2) * t44 + t52) + g(3) * (rSges(5,3) * t29 + t22) + (rSges(5,3) * t69 + g(3) * t33) * t26 + (t61 * t26 + t33 * t29) * t67) - m(6) * (g(1) * (rSges(6,1) * t43 + rSges(6,3) * t44 + t39) + g(3) * (rSges(6,2) * t29 + t22) + (rSges(6,2) * t69 + g(3) * t73) * t26 + (t62 * t26 + t73 * t29) * t67) - m(7) * (g(1) * t39 + g(3) * t22 + (g(3) * t47 + (rSges(7,2) * t25 + t63 * t28) * t69) * t29 + (g(3) * t31 + t47 * t69) * t26 + (t42 * t26 + t31 * t29) * t67) -m(5) * (g(1) * (-rSges(5,1) * t6 - rSges(5,2) * t7) + g(2) * (rSges(5,1) * t8 + rSges(5,2) * t9)) - m(6) * (g(1) * (-rSges(6,1) * t6 + t48 * t7 - t2) + g(2) * (rSges(6,1) * t8 - t48 * t9 + t4) + t66) - m(7) * (g(1) * (t49 * t7 - t63 * t6 - t2) + g(2) * (-t49 * t9 + t63 * t8 + t4) + t66) + ((m(5) * rSges(5,2) - m(6) * rSges(6,3) - m(7) * rSges(7,2)) * t28 + (m(5) * rSges(5,1) - m(6) * t64 - m(7) * t46) * t25) * t65, t71 * (g(1) * t6 - g(2) * t8 + t25 * t65) -m(7) * (-g(3) * t26 + t37 * t29)];
taug  = t1(:);
