% Calculate Gravitation load on the joints for
% S6RPRRPP8
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
% Datum: 2019-03-09 04:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPP8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP8_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:54:04
% EndTime: 2019-03-09 04:54:06
% DurationCPUTime: 0.76s
% Computational Cost: add. (250->141), mult. (526->184), div. (0->0), fcn. (522->6), ass. (0->57)
t25 = sin(qJ(4));
t28 = cos(qJ(4));
t47 = rSges(6,3) + qJ(5);
t61 = rSges(6,2) - pkin(4);
t74 = -t25 * t47 + t28 * t61 - pkin(3);
t69 = -pkin(1) - pkin(7);
t29 = cos(qJ(3));
t54 = t28 * t29;
t73 = g(3) * t54;
t30 = cos(qJ(1));
t66 = g(2) * t30;
t27 = sin(qJ(1));
t68 = g(1) * t27;
t71 = -t66 + t68;
t46 = rSges(7,3) + qJ(6);
t70 = -m(6) - m(7);
t67 = g(2) * t29;
t65 = qJ(5) * t73;
t64 = g(3) * t29;
t63 = -rSges(6,1) - pkin(8);
t62 = rSges(7,1) + pkin(5);
t60 = -rSges(5,3) - pkin(8);
t26 = sin(qJ(3));
t59 = t26 * t27;
t58 = t26 * t30;
t57 = t27 * t25;
t56 = t27 * t28;
t55 = t27 * t29;
t53 = t29 * t30;
t52 = t30 * t28;
t51 = pkin(3) * t55 + pkin(8) * t59;
t21 = t30 * qJ(2);
t50 = pkin(3) * t58 + t21;
t49 = t30 * pkin(1) + t27 * qJ(2);
t48 = rSges(7,2) + qJ(5);
t45 = -pkin(8) - t62;
t44 = pkin(8) * t53;
t43 = t25 * t55;
t42 = t27 * t54;
t41 = -pkin(4) - t46;
t40 = t30 * pkin(7) + t49;
t39 = g(1) * t69;
t38 = pkin(4) * t42 + qJ(5) * t43 + t51;
t37 = pkin(3) * t59 + t40;
t35 = rSges(4,1) * t26 + rSges(4,2) * t29;
t8 = t25 * t58 + t56;
t9 = t26 * t52 - t57;
t34 = t9 * pkin(4) + t8 * qJ(5) + t50;
t33 = -rSges(5,1) * t28 + rSges(5,2) * t25 - pkin(3);
t6 = t26 * t57 - t52;
t7 = t25 * t30 + t26 * t56;
t32 = t7 * pkin(4) + qJ(5) * t6 + t37;
t31 = -t25 * t48 + t28 * t41 - pkin(3);
t22 = t29 * pkin(8);
t4 = t8 * pkin(4);
t2 = t6 * pkin(4);
t1 = [-m(2) * (g(1) * (-t27 * rSges(2,1) - rSges(2,2) * t30) + g(2) * (rSges(2,1) * t30 - t27 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t30 + t21 + (rSges(3,2) - pkin(1)) * t27) + g(2) * (-rSges(3,2) * t30 + t27 * rSges(3,3) + t49)) - m(4) * (g(1) * (rSges(4,1) * t58 + rSges(4,2) * t53 + t21) + g(2) * (rSges(4,3) * t30 + t40) + (g(1) * (-rSges(4,3) + t69) + g(2) * t35) * t27) - m(5) * (g(1) * (t9 * rSges(5,1) - t8 * rSges(5,2) - rSges(5,3) * t53 - t44 + t50) + g(2) * (rSges(5,1) * t7 - rSges(5,2) * t6 + t37) + (t60 * t67 + t39) * t27) - m(6) * (g(1) * (-rSges(6,1) * t53 - t9 * rSges(6,2) + t8 * rSges(6,3) + t34 - t44) + g(2) * (-rSges(6,2) * t7 + rSges(6,3) * t6 + t32) + (t63 * t67 + t39) * t27) - m(7) * (g(1) * (t8 * rSges(7,2) + t69 * t27 + t46 * t9 + t34) + g(2) * (rSges(7,2) * t6 + t46 * t7 + t32) + (g(1) * t30 + g(2) * t27) * t29 * t45) (-m(3) - m(4) - m(5) + t70) * t71, -m(4) * (-g(3) * t35 + t71 * (rSges(4,1) * t29 - rSges(4,2) * t26)) - m(5) * (g(1) * (rSges(5,1) * t42 - rSges(5,2) * t43 + t51) + g(3) * (rSges(5,3) * t29 + t22) + (rSges(5,3) * t68 + g(3) * t33) * t26 + (t26 * t60 + t29 * t33) * t66) - m(6) * (g(1) * (-rSges(6,2) * t42 + rSges(6,3) * t43 + t38) + g(3) * (rSges(6,1) * t29 + t22) + (rSges(6,1) * t68 + g(3) * t74) * t26 + (t63 * t26 + t74 * t29) * t66) - m(7) * (g(1) * t38 + g(3) * t22 + (g(3) * t62 + (rSges(7,2) * t25 + t28 * t46) * t68) * t29 + (g(3) * t31 + t62 * t68) * t26 + (t26 * t45 + t29 * t31) * t66) -m(5) * (g(1) * (-rSges(5,1) * t6 - rSges(5,2) * t7) + g(2) * (rSges(5,1) * t8 + rSges(5,2) * t9)) - m(6) * (g(1) * (rSges(6,2) * t6 + t47 * t7 - t2) + g(2) * (-rSges(6,2) * t8 - t47 * t9 + t4) + t65) - m(7) * (g(1) * (-t46 * t6 + t48 * t7 - t2) + g(2) * (t46 * t8 - t48 * t9 + t4) + t65) + ((m(5) * rSges(5,2) - m(6) * rSges(6,3) - m(7) * rSges(7,2)) * t28 + (m(5) * rSges(5,1) - m(6) * t61 - m(7) * t41) * t25) * t64, t70 * (g(1) * t6 - g(2) * t8 + t25 * t64) -m(7) * (g(1) * t7 - g(2) * t9 + t73)];
taug  = t1(:);
