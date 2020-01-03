% Calculate Gravitation load on the joints for
% S5RRRPR3
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
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:08:50
% EndTime: 2020-01-03 12:08:51
% DurationCPUTime: 0.30s
% Computational Cost: add. (290->73), mult. (218->91), div. (0->0), fcn. (173->10), ass. (0->44)
t41 = qJ(3) + pkin(9);
t32 = sin(t41);
t33 = cos(t41);
t46 = cos(qJ(3));
t38 = t46 * pkin(3);
t79 = t33 * rSges(5,1) - t32 * rSges(5,2) + t38;
t78 = rSges(4,3) + pkin(7);
t43 = -qJ(4) - pkin(7);
t77 = rSges(5,3) - t43;
t76 = rSges(6,3) + pkin(8) - t43;
t34 = qJ(5) + t41;
t26 = sin(t34);
t27 = cos(t34);
t57 = t27 * rSges(6,1) - t26 * rSges(6,2);
t44 = sin(qJ(3));
t75 = t46 * rSges(4,1) - t44 * rSges(4,2);
t73 = pkin(2) + t75;
t72 = pkin(2) + t79;
t69 = t44 * pkin(3);
t71 = m(4) * (rSges(4,1) * t44 + rSges(4,2) * t46) + m(5) * (rSges(5,1) * t32 + rSges(5,2) * t33 + t69);
t42 = qJ(1) + qJ(2);
t35 = sin(t42);
t70 = g(2) * t35;
t36 = cos(t42);
t67 = t26 * t36;
t66 = t27 * t36;
t61 = pkin(4) * t33 + t38;
t60 = t35 * rSges(3,1) + t36 * rSges(3,2);
t59 = g(3) * (rSges(6,1) * t67 + rSges(6,2) * t66);
t58 = t36 * rSges(3,1) - t35 * rSges(3,2);
t55 = -rSges(6,1) * t26 - rSges(6,2) * t27;
t54 = t78 * t35 + t73 * t36;
t8 = pkin(2) + t61;
t52 = -t76 * t36 + (t8 + t57) * t35;
t51 = rSges(6,1) * t66 - rSges(6,2) * t67 + t76 * t35 + t36 * t8;
t50 = t72 * t35 - t77 * t36;
t49 = t77 * t35 + t72 * t36;
t48 = t73 * t35 - t78 * t36;
t47 = cos(qJ(1));
t45 = sin(qJ(1));
t39 = t47 * pkin(1);
t37 = t45 * pkin(1);
t9 = -pkin(4) * t32 - t69;
t1 = [-m(2) * (g(2) * (t47 * rSges(2,1) - t45 * rSges(2,2)) + g(3) * (t45 * rSges(2,1) + t47 * rSges(2,2))) - m(3) * (g(2) * (t39 + t58) + g(3) * (t37 + t60)) - m(4) * (g(2) * (t39 + t54) + g(3) * (t37 + t48)) - m(5) * (g(2) * (t39 + t49) + g(3) * (t37 + t50)) - m(6) * (g(2) * (t39 + t51) + g(3) * (t37 + t52)), -m(3) * (g(2) * t58 + g(3) * t60) - m(4) * (g(2) * t54 + g(3) * t48) - m(5) * (g(2) * t49 + g(3) * t50) - m(6) * (g(2) * t51 + g(3) * t52), -m(6) * t59 + (-m(4) * t75 - m(5) * t79 - m(6) * (t57 + t61)) * g(1) + (m(6) * t9 - t71) * g(3) * t36 + (-m(6) * (t55 + t9) + t71) * t70, (-m(5) - m(6)) * (-g(2) * t36 - g(3) * t35), -m(6) * (g(1) * t57 + t55 * t70 + t59)];
taug = t1(:);
