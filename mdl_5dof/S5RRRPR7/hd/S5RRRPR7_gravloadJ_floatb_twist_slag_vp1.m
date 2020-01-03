% Calculate Gravitation load on the joints for
% S5RRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR7_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:16:07
% EndTime: 2019-12-31 21:16:09
% DurationCPUTime: 0.67s
% Computational Cost: add. (327->101), mult. (362->139), div. (0->0), fcn. (325->10), ass. (0->49)
t33 = qJ(2) + qJ(3);
t29 = sin(t33);
t30 = cos(t33);
t35 = cos(pkin(9));
t58 = -rSges(5,1) * t35 - pkin(3);
t34 = sin(pkin(9));
t70 = rSges(5,2) * t34;
t48 = (rSges(5,3) + qJ(4)) * t29 + (-t58 - t70) * t30;
t83 = qJ(4) * t30 + t29 * t70;
t80 = t30 * rSges(4,1) - t29 * rSges(4,2);
t38 = sin(qJ(1));
t40 = cos(qJ(1));
t79 = g(1) * t40 + g(2) * t38;
t24 = t35 * pkin(4) + pkin(3);
t36 = -pkin(8) - qJ(4);
t45 = t30 * t24 + (rSges(6,3) - t36) * t29;
t78 = t79 * t29;
t37 = sin(qJ(2));
t77 = pkin(2) * t37;
t39 = cos(qJ(2));
t31 = t39 * pkin(2);
t26 = t31 + pkin(1);
t15 = t40 * t26;
t75 = g(2) * t15;
t73 = rSges(3,3) + pkin(6);
t32 = pkin(9) + qJ(5);
t28 = cos(t32);
t71 = rSges(6,1) * t28;
t27 = sin(t32);
t69 = rSges(6,2) * t27;
t66 = t38 * t30;
t65 = t40 * t30;
t41 = -pkin(7) - pkin(6);
t64 = rSges(4,3) - t41;
t62 = rSges(5,3) * t66 + t83 * t38;
t61 = rSges(5,3) * t65 + t83 * t40;
t59 = t29 * t69;
t57 = pkin(4) * t34 - t41;
t56 = -t24 - t71;
t54 = t39 * rSges(3,1) - t37 * rSges(3,2);
t51 = pkin(1) + t54;
t49 = t34 * rSges(5,1) + t35 * rSges(5,2) - t41;
t46 = g(1) * (rSges(6,3) * t65 + t40 * t59) + g(2) * (rSges(6,3) * t66 + t38 * t59);
t44 = t45 + (-t69 + t71) * t30;
t5 = t38 * t27 + t28 * t65;
t4 = -t27 * t65 + t38 * t28;
t3 = t40 * t27 - t28 * t66;
t2 = t27 * t66 + t40 * t28;
t1 = [-m(2) * (g(1) * (-t38 * rSges(2,1) - t40 * rSges(2,2)) + g(2) * (t40 * rSges(2,1) - t38 * rSges(2,2))) - m(3) * ((g(1) * t73 + g(2) * t51) * t40 + (-g(1) * t51 + g(2) * t73) * t38) - m(4) * (t75 + (g(1) * t64 + g(2) * t80) * t40 + (g(1) * (-t26 - t80) + g(2) * t64) * t38) - m(5) * (t75 + (g(1) * t49 + g(2) * t48) * t40 + (g(2) * t49 + (-t26 - t48) * g(1)) * t38) - m(6) * (g(1) * (t3 * rSges(6,1) + t2 * rSges(6,2)) + g(2) * (t5 * rSges(6,1) + t4 * rSges(6,2) + t15) + (g(1) * t57 + g(2) * t45) * t40 + (g(1) * (-t26 - t45) + g(2) * t57) * t38), -m(3) * (g(3) * t54 + t79 * (-rSges(3,1) * t37 - rSges(3,2) * t39)) - m(4) * (g(3) * (t31 + t80) + t79 * (-rSges(4,1) * t29 - rSges(4,2) * t30 - t77)) - m(5) * (g(1) * (-t40 * t77 + t61) + g(2) * (-t38 * t77 + t62) + g(3) * (t31 + t48) + t58 * t78) - m(6) * (g(3) * (t31 + t44) + t46 + t79 * (t56 * t29 - t30 * t36 - t77)), -m(4) * g(3) * t80 - m(5) * (g(1) * t61 + g(2) * t62 + g(3) * t48) - m(6) * (g(3) * t44 + t46) + t79 * ((m(4) * rSges(4,2) + m(6) * t36) * t30 + (m(4) * rSges(4,1) - m(5) * t58 - m(6) * t56) * t29), (-m(5) - m(6)) * (-g(3) * t30 + t78), -m(6) * (g(1) * (t4 * rSges(6,1) - t5 * rSges(6,2)) + g(2) * (-t2 * rSges(6,1) + t3 * rSges(6,2)) + g(3) * (-rSges(6,1) * t27 - rSges(6,2) * t28) * t29)];
taug = t1(:);
