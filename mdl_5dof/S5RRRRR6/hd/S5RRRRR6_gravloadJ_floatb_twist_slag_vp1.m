% Calculate Gravitation load on the joints for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2020-01-03 12:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:14:46
% EndTime: 2020-01-03 12:14:47
% DurationCPUTime: 0.31s
% Computational Cost: add. (328->84), mult. (245->108), div. (0->0), fcn. (194->10), ass. (0->50)
t81 = rSges(4,3) + pkin(7);
t49 = -pkin(8) - pkin(7);
t80 = rSges(5,3) - t49;
t79 = rSges(6,3) + pkin(9) - t49;
t43 = qJ(3) + qJ(4);
t38 = qJ(5) + t43;
t31 = sin(t38);
t32 = cos(t38);
t60 = t32 * rSges(6,1) - t31 * rSges(6,2);
t34 = sin(t43);
t36 = cos(t43);
t61 = t36 * rSges(5,1) - t34 * rSges(5,2);
t45 = sin(qJ(3));
t47 = cos(qJ(3));
t78 = t47 * rSges(4,1) - t45 * rSges(4,2);
t77 = pkin(2) + t78;
t44 = qJ(1) + qJ(2);
t37 = cos(t44);
t70 = t32 * t37;
t71 = t31 * t37;
t76 = rSges(6,1) * t71 + rSges(6,2) * t70;
t75 = pkin(4) * t34;
t35 = sin(t44);
t74 = g(2) * t35;
t73 = t45 * pkin(3);
t40 = t47 * pkin(3);
t33 = t40 + pkin(2);
t68 = t34 * t37;
t67 = t36 * t37;
t64 = rSges(5,1) * t68 + rSges(5,2) * t67;
t63 = t35 * rSges(3,1) + t37 * rSges(3,2);
t62 = t37 * rSges(3,1) - t35 * rSges(3,2);
t30 = pkin(4) * t36;
t59 = t30 + t60;
t58 = rSges(4,1) * t45 + rSges(4,2) * t47;
t57 = -rSges(5,1) * t34 - rSges(5,2) * t36;
t56 = -rSges(6,1) * t31 - rSges(6,2) * t32;
t55 = t81 * t35 + t77 * t37;
t7 = t30 + t33;
t54 = -t79 * t37 + (t7 + t60) * t35;
t53 = rSges(6,1) * t70 - rSges(6,2) * t71 + t79 * t35 + t37 * t7;
t52 = -t80 * t37 + (t33 + t61) * t35;
t51 = rSges(5,1) * t67 - rSges(5,2) * t68 + t37 * t33 + t80 * t35;
t50 = t77 * t35 - t81 * t37;
t48 = cos(qJ(1));
t46 = sin(qJ(1));
t41 = t48 * pkin(1);
t39 = t46 * pkin(1);
t8 = -t73 - t75;
t1 = [-m(2) * (g(2) * (t48 * rSges(2,1) - t46 * rSges(2,2)) + g(3) * (t46 * rSges(2,1) + t48 * rSges(2,2))) - m(3) * (g(2) * (t41 + t62) + g(3) * (t39 + t63)) - m(4) * (g(2) * (t41 + t55) + g(3) * (t39 + t50)) - m(5) * (g(2) * (t41 + t51) + g(3) * (t39 + t52)) - m(6) * (g(2) * (t41 + t53) + g(3) * (t39 + t54)), -m(3) * (g(2) * t62 + g(3) * t63) - m(4) * (g(2) * t55 + g(3) * t50) - m(5) * (g(2) * t51 + g(3) * t52) - m(6) * (g(2) * t53 + g(3) * t54), -m(4) * (g(3) * t58 * t37 + g(1) * t78) - m(5) * (g(1) * (t40 + t61) + g(3) * (t37 * t73 + t64)) - m(6) * (g(1) * (t40 + t59) + g(3) * (-t37 * t8 + t76)) + (m(4) * t58 - m(5) * (t57 - t73) - m(6) * (t56 + t8)) * t74, -m(5) * (g(1) * t61 + g(3) * t64) - m(6) * (g(1) * t59 + g(3) * (pkin(4) * t68 + t76)) + (-m(5) * t57 - m(6) * (t56 - t75)) * t74, -m(6) * (g(1) * t60 + g(3) * t76 + t56 * t74)];
taug = t1(:);
