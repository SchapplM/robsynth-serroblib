% Calculate Gravitation load on the joints for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:05:30
% EndTime: 2020-01-03 12:05:31
% DurationCPUTime: 0.29s
% Computational Cost: add. (335->91), mult. (308->123), div. (0->0), fcn. (294->10), ass. (0->48)
t80 = rSges(5,3) + pkin(7);
t49 = qJ(4) + qJ(5);
t43 = sin(t49);
t45 = cos(t49);
t50 = qJ(1) + qJ(2);
t46 = cos(t50);
t44 = sin(t50);
t52 = cos(pkin(9));
t74 = t44 * t52;
t10 = -t46 * t43 + t45 * t74;
t9 = -t43 * t74 - t46 * t45;
t79 = t9 * rSges(6,1) - t10 * rSges(6,2);
t71 = t46 * t52;
t11 = t43 * t71 - t44 * t45;
t12 = t44 * t43 + t45 * t71;
t78 = t11 * rSges(6,1) + t12 * rSges(6,2);
t53 = sin(qJ(4));
t77 = pkin(4) * t53;
t51 = sin(pkin(9));
t76 = g(1) * t51;
t75 = t44 * t51;
t73 = t44 * t53;
t72 = t46 * t51;
t70 = t51 * (-pkin(8) - pkin(7));
t69 = t52 * t53;
t55 = cos(qJ(4));
t68 = t52 * t55;
t67 = t46 * pkin(2) + t44 * qJ(3);
t66 = t44 * rSges(3,1) + t46 * rSges(3,2);
t65 = t46 * rSges(3,1) - t44 * rSges(3,2);
t64 = -rSges(6,1) * t43 - rSges(6,2) * t45;
t19 = -t44 * t55 + t46 * t69;
t17 = -t44 * t69 - t46 * t55;
t20 = t46 * t68 + t73;
t63 = t20 * rSges(5,1) - t19 * rSges(5,2) + pkin(3) * t71 + t80 * t72 + t67;
t62 = rSges(4,1) * t71 - rSges(4,2) * t72 + t44 * rSges(4,3) + t67;
t18 = t44 * t68 - t46 * t53;
t40 = t44 * pkin(2);
t61 = t18 * rSges(5,1) + t17 * rSges(5,2) + pkin(3) * t74 - t46 * qJ(3) + t80 * t75 + t40;
t60 = -rSges(4,2) * t75 + rSges(4,1) * t74 + t40 + (-rSges(4,3) - qJ(3)) * t46;
t42 = t55 * pkin(4) + pkin(3);
t59 = t12 * rSges(6,1) - t11 * rSges(6,2) + rSges(6,3) * t72 + pkin(4) * t73 + t42 * t71 - t46 * t70 + t67;
t58 = -t44 * t70 + t9 * rSges(6,2) + t42 * t74 + rSges(6,3) * t75 + t10 * rSges(6,1) + t40 + (-qJ(3) - t77) * t46;
t56 = cos(qJ(1));
t54 = sin(qJ(1));
t48 = t56 * pkin(1);
t47 = t54 * pkin(1);
t1 = [-m(2) * (g(2) * (t56 * rSges(2,1) - t54 * rSges(2,2)) + g(3) * (t54 * rSges(2,1) + t56 * rSges(2,2))) - m(3) * (g(2) * (t48 + t65) + g(3) * (t47 + t66)) - m(4) * (g(2) * (t48 + t62) + g(3) * (t47 + t60)) - m(5) * (g(2) * (t48 + t63) + g(3) * (t47 + t61)) - m(6) * (g(2) * (t48 + t59) + g(3) * (t47 + t58)), -m(3) * (g(2) * t65 + g(3) * t66) - m(4) * (g(2) * t62 + g(3) * t60) - m(5) * (g(2) * t63 + g(3) * t61) - m(6) * (g(2) * t59 + g(3) * t58), (-m(4) - m(5) - m(6)) * (-g(2) * t46 - g(3) * t44), -m(5) * (g(2) * (t17 * rSges(5,1) - t18 * rSges(5,2)) + g(3) * (t19 * rSges(5,1) + t20 * rSges(5,2))) - m(6) * (g(2) * (t17 * pkin(4) + t79) + g(3) * (t19 * pkin(4) + t78)) + (-m(5) * (-rSges(5,1) * t53 - rSges(5,2) * t55) - m(6) * (t64 - t77)) * t76, -m(6) * (g(2) * t79 + g(3) * t78 + t64 * t76)];
taug = t1(:);
