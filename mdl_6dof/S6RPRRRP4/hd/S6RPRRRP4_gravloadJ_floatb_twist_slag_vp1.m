% Calculate Gravitation load on the joints for
% S6RPRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:06:58
% EndTime: 2019-03-09 06:06:59
% DurationCPUTime: 0.67s
% Computational Cost: add. (489->124), mult. (432->169), div. (0->0), fcn. (388->10), ass. (0->57)
t87 = rSges(6,3) + pkin(9);
t86 = rSges(7,1) + pkin(5);
t36 = pkin(10) + qJ(3);
t33 = qJ(4) + t36;
t27 = sin(t33);
t28 = cos(t33);
t85 = t28 * rSges(5,1) - t27 * rSges(5,2);
t42 = sin(qJ(1));
t44 = cos(qJ(1));
t84 = g(1) * t44 + g(2) * t42;
t43 = cos(qJ(5));
t30 = t43 * pkin(5) + pkin(4);
t39 = -qJ(6) - pkin(9);
t48 = t28 * t30 + (rSges(7,3) - t39) * t27;
t52 = t28 * pkin(4) + t87 * t27;
t83 = t84 * t27;
t31 = sin(t36);
t82 = pkin(3) * t31;
t41 = sin(qJ(5));
t81 = pkin(5) * t41;
t38 = cos(pkin(10));
t29 = t38 * pkin(2) + pkin(1);
t76 = t28 * t41;
t75 = t28 * t42;
t74 = t28 * t43;
t73 = t28 * t44;
t72 = t42 * t41;
t71 = t42 * t43;
t70 = t44 * t41;
t69 = t44 * t43;
t40 = -pkin(7) - qJ(2);
t68 = rSges(4,3) - t40;
t35 = -pkin(8) + t40;
t67 = rSges(5,3) - t35;
t66 = rSges(3,3) + qJ(2);
t63 = t27 * t72;
t65 = rSges(6,2) * t63 + t75 * t87;
t62 = t27 * t70;
t64 = rSges(6,2) * t62 + t73 * t87;
t61 = -rSges(6,1) * t43 - pkin(4);
t60 = -t35 + t81;
t59 = -rSges(7,1) * t43 - t30;
t32 = cos(t36);
t57 = t32 * rSges(4,1) - t31 * rSges(4,2);
t54 = rSges(3,1) * t38 - rSges(3,2) * sin(pkin(10)) + pkin(1);
t3 = -t28 * t70 + t71;
t1 = t28 * t72 + t69;
t53 = t29 + t57;
t51 = rSges(6,1) * t74 - rSges(6,2) * t76 + t52;
t49 = g(1) * (rSges(7,2) * t62 + rSges(7,3) * t73) + g(2) * (rSges(7,2) * t63 + rSges(7,3) * t75);
t47 = rSges(7,1) * t74 - rSges(7,2) * t76 + t48;
t26 = pkin(3) * t32;
t10 = t26 + t29;
t5 = t44 * t10;
t4 = t28 * t69 + t72;
t2 = -t28 * t71 + t70;
t6 = [-m(2) * (g(1) * (-t42 * rSges(2,1) - t44 * rSges(2,2)) + g(2) * (t44 * rSges(2,1) - t42 * rSges(2,2))) - m(3) * ((g(1) * t66 + g(2) * t54) * t44 + (-g(1) * t54 + g(2) * t66) * t42) - m(4) * ((g(1) * t68 + g(2) * t53) * t44 + (-g(1) * t53 + g(2) * t68) * t42) - m(5) * (g(2) * t5 + (g(1) * t67 + g(2) * t85) * t44 + (g(1) * (-t10 - t85) + g(2) * t67) * t42) - m(6) * (g(1) * (t2 * rSges(6,1) + t1 * rSges(6,2)) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t5) + (-g(1) * t35 + g(2) * t52) * t44 + (g(1) * (-t10 - t52) - g(2) * t35) * t42) - m(7) * (g(1) * (t2 * rSges(7,1) + t1 * rSges(7,2)) + g(2) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t5) + (g(1) * t60 + g(2) * t48) * t44 + (g(1) * (-t10 - t48) + g(2) * t60) * t42) (-m(3) - m(4) - m(5) - m(6) - m(7)) * (g(1) * t42 - g(2) * t44) -m(4) * (g(3) * t57 + t84 * (-rSges(4,1) * t31 - rSges(4,2) * t32)) - m(5) * (g(3) * (t26 + t85) + t84 * (-rSges(5,1) * t27 - rSges(5,2) * t28 - t82)) - m(6) * (g(1) * (-t44 * t82 + t64) + g(2) * (-t42 * t82 + t65) + g(3) * (t26 + t51) + t61 * t83) - m(7) * (g(3) * (t26 + t47) + t49 + t84 * (t59 * t27 - t28 * t39 - t82)) -m(5) * g(3) * t85 - m(6) * (g(1) * t64 + g(2) * t65 + g(3) * t51) - m(7) * (g(3) * t47 + t49) + t84 * ((m(5) * rSges(5,2) + m(7) * t39) * t28 + (m(5) * rSges(5,1) - m(6) * t61 - m(7) * t59) * t27) -m(6) * (g(1) * (t3 * rSges(6,1) - t4 * rSges(6,2)) + g(2) * (-t1 * rSges(6,1) + t2 * rSges(6,2))) - m(7) * (g(1) * (-t4 * rSges(7,2) + t86 * t3) + g(2) * (t2 * rSges(7,2) - t86 * t1)) + (-m(6) * (-rSges(6,1) * t41 - rSges(6,2) * t43) - m(7) * (-rSges(7,1) * t41 - rSges(7,2) * t43 - t81)) * g(3) * t27, -m(7) * (-g(3) * t28 + t83)];
taug  = t6(:);
