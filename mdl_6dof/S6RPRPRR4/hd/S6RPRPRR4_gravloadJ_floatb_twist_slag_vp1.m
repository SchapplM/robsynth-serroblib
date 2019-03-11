% Calculate Gravitation load on the joints for
% S6RPRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:44:27
% EndTime: 2019-03-09 03:44:29
% DurationCPUTime: 0.57s
% Computational Cost: add. (376->116), mult. (389->160), div. (0->0), fcn. (355->10), ass. (0->54)
t36 = cos(qJ(3));
t72 = g(3) * t36;
t30 = qJ(1) + pkin(10);
t24 = cos(t30);
t23 = sin(t30);
t73 = g(2) * t23;
t75 = g(1) * t24 + t73;
t70 = rSges(6,3) + pkin(8);
t60 = rSges(7,3) + pkin(9) + pkin(8);
t33 = sin(qJ(3));
t43 = -t36 * rSges(5,2) + t33 * rSges(5,3);
t46 = t36 * rSges(4,1) - t33 * rSges(4,2);
t34 = sin(qJ(1));
t71 = t34 * pkin(1);
t28 = t36 * pkin(3);
t69 = t24 * t36;
t31 = qJ(5) + qJ(6);
t25 = sin(t31);
t68 = t25 * t33;
t26 = cos(t31);
t67 = t26 * t33;
t32 = sin(qJ(5));
t66 = t32 * t33;
t35 = cos(qJ(5));
t63 = t33 * t35;
t27 = t33 * qJ(4);
t59 = t27 + t28;
t57 = -m(5) - m(6) - m(7);
t56 = -pkin(3) - t70;
t55 = pkin(5) * t66;
t54 = -pkin(3) - t60;
t37 = cos(qJ(1));
t29 = t37 * pkin(1);
t53 = t24 * pkin(2) + t23 * pkin(7) + t29;
t52 = t24 * pkin(7) - t71;
t51 = -pkin(2) - t27;
t50 = g(1) * t56;
t49 = g(1) * t54;
t48 = pkin(3) * t69 + t24 * t27 + t53;
t47 = t75 * qJ(4) * t36;
t44 = rSges(6,1) * t32 + rSges(6,2) * t35;
t10 = -t23 * t32 + t24 * t63;
t12 = t23 * t63 + t24 * t32;
t41 = rSges(7,1) * t25 + rSges(7,2) * t26 + pkin(5) * t32;
t40 = g(3) * t59 + t47;
t5 = -t23 * t25 + t24 * t67;
t6 = t23 * t26 + t24 * t68;
t7 = t23 * t67 + t24 * t25;
t8 = -t23 * t68 + t24 * t26;
t39 = g(1) * (t5 * rSges(7,1) - t6 * rSges(7,2)) + g(2) * (t7 * rSges(7,1) + t8 * rSges(7,2)) + (-rSges(7,1) * t26 + rSges(7,2) * t25) * t72;
t22 = t35 * pkin(5) + pkin(4);
t13 = -t23 * t66 + t24 * t35;
t11 = t23 * t35 + t24 * t66;
t1 = [-m(2) * (g(1) * (-t34 * rSges(2,1) - t37 * rSges(2,2)) + g(2) * (t37 * rSges(2,1) - t34 * rSges(2,2))) - m(3) * (g(1) * (-t23 * rSges(3,1) - t24 * rSges(3,2) - t71) + g(2) * (t24 * rSges(3,1) - t23 * rSges(3,2) + t29)) - m(4) * (g(1) * (t24 * rSges(4,3) + t52) + g(2) * (t46 * t24 + t53) + (g(1) * (-pkin(2) - t46) + g(2) * rSges(4,3)) * t23) - m(5) * (g(1) * (t24 * rSges(5,1) + t52) + g(2) * (t43 * t24 + t48) + (g(1) * (-t43 + t51 - t28) + g(2) * rSges(5,1)) * t23) - m(6) * (g(1) * (t13 * rSges(6,1) - t12 * rSges(6,2) + t24 * pkin(4) + t52) + g(2) * (t11 * rSges(6,1) + t10 * rSges(6,2) + t70 * t69 + t48) + (g(2) * pkin(4) + g(1) * t51 + t36 * t50) * t23) - m(7) * (g(1) * (t8 * rSges(7,1) - t7 * rSges(7,2) + t52) + g(2) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t48) + (g(1) * t22 + g(2) * (t60 * t36 + t55)) * t24 + (g(1) * (t51 - t55) + g(2) * t22 + t36 * t49) * t23) (-m(3) - m(4) + t57) * g(3), -m(4) * (g(3) * t46 + t75 * (-rSges(4,1) * t33 - rSges(4,2) * t36)) - m(5) * (g(3) * (t43 + t59) + t47 + t75 * (rSges(5,3) * t36 + (rSges(5,2) - pkin(3)) * t33)) - m(6) * ((g(3) * t70 + t75 * t44) * t36 + (g(3) * t44 + t24 * t50 + t56 * t73) * t33 + t40) - m(7) * ((g(3) * t60 + t75 * t41) * t36 + (g(3) * t41 + t24 * t49 + t54 * t73) * t33 + t40) t57 * (t75 * t33 - t72) -m(6) * (g(1) * (t10 * rSges(6,1) - t11 * rSges(6,2)) + g(2) * (t12 * rSges(6,1) + t13 * rSges(6,2)) + (-rSges(6,1) * t35 + rSges(6,2) * t32) * t72) - m(7) * ((g(1) * t10 + g(2) * t12 - t35 * t72) * pkin(5) + t39) -m(7) * t39];
taug  = t1(:);
