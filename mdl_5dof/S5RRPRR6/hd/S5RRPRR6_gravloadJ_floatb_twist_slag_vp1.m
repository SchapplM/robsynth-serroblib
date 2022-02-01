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
% m [6x1]
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
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 11:16:50
% EndTime: 2022-01-20 11:16:51
% DurationCPUTime: 0.33s
% Computational Cost: add. (335->93), mult. (308->124), div. (0->0), fcn. (294->10), ass. (0->50)
t78 = rSges(5,3) + pkin(7);
t44 = qJ(4) + qJ(5);
t39 = sin(t44);
t41 = cos(t44);
t45 = qJ(1) + qJ(2);
t42 = cos(t45);
t40 = sin(t45);
t47 = cos(pkin(9));
t72 = t40 * t47;
t10 = t42 * t39 - t41 * t72;
t9 = t39 * t72 + t42 * t41;
t77 = -t9 * rSges(6,1) + t10 * rSges(6,2);
t69 = t42 * t47;
t11 = -t39 * t69 + t40 * t41;
t12 = t40 * t39 + t41 * t69;
t76 = t11 * rSges(6,1) - t12 * rSges(6,2);
t75 = g(1) * t40;
t46 = sin(pkin(9));
t74 = g(3) * t46;
t49 = sin(qJ(1));
t73 = t49 * pkin(1);
t48 = sin(qJ(4));
t71 = t40 * t48;
t70 = t42 * t46;
t68 = t42 * t48;
t67 = t46 * (-pkin(8) - pkin(7));
t66 = t47 * t48;
t50 = cos(qJ(4));
t65 = t47 * t50;
t64 = t42 * pkin(2) + t40 * qJ(3);
t17 = t40 * t66 + t42 * t50;
t18 = -t40 * t65 + t68;
t33 = t42 * qJ(3);
t63 = t18 * rSges(5,1) + t17 * rSges(5,2) + t33;
t62 = t42 * rSges(3,1) - t40 * rSges(3,2);
t61 = t10 * rSges(6,1) + t9 * rSges(6,2) + pkin(4) * t68 + t40 * t67 + t33;
t60 = -t40 * rSges(3,1) - t42 * rSges(3,2);
t59 = -rSges(6,1) * t39 - rSges(6,2) * t41;
t19 = t40 * t50 - t42 * t66;
t20 = t42 * t65 + t71;
t58 = t20 * rSges(5,1) + t19 * rSges(5,2) + pkin(3) * t69 + t78 * t70 + t64;
t57 = rSges(4,1) * t69 - rSges(4,2) * t70 + t40 * rSges(4,3) + t64;
t38 = t50 * pkin(4) + pkin(3);
t56 = (-rSges(6,3) * t46 - t38 * t47 - pkin(2)) * t75;
t55 = t42 * rSges(4,3) + t33 + (-rSges(4,1) * t47 + rSges(4,2) * t46 - pkin(2)) * t40;
t54 = (-pkin(3) * t47 - t78 * t46 - pkin(2)) * t75;
t53 = t12 * rSges(6,1) + t11 * rSges(6,2) + rSges(6,3) * t70 + pkin(4) * t71 + t38 * t69 - t42 * t67 + t64;
t51 = cos(qJ(1));
t43 = t51 * pkin(1);
t1 = [-m(2) * (g(1) * (-t49 * rSges(2,1) - t51 * rSges(2,2)) + g(2) * (t51 * rSges(2,1) - t49 * rSges(2,2))) - m(3) * (g(1) * (t60 - t73) + g(2) * (t43 + t62)) - m(4) * (g(1) * (t55 - t73) + g(2) * (t43 + t57)) - m(5) * (g(1) * (t63 - t73) + g(2) * (t43 + t58) + t54) - m(6) * (g(1) * (t61 - t73) + g(2) * (t43 + t53) + t56), -m(3) * (g(1) * t60 + g(2) * t62) - m(4) * (g(1) * t55 + g(2) * t57) - m(5) * (g(1) * t63 + g(2) * t58 + t54) - m(6) * (g(1) * t61 + g(2) * t53 + t56), (-m(4) - m(5) - m(6)) * (-g(2) * t42 + t75), -m(5) * (g(1) * (t19 * rSges(5,1) - t20 * rSges(5,2)) + g(2) * (-t17 * rSges(5,1) + t18 * rSges(5,2))) - m(6) * (g(1) * (t19 * pkin(4) + t76) + g(2) * (-t17 * pkin(4) + t77)) + (-m(5) * (-rSges(5,1) * t48 - rSges(5,2) * t50) - m(6) * (-pkin(4) * t48 + t59)) * t74, -m(6) * (g(1) * t76 + g(2) * t77 + t59 * t74)];
taug = t1(:);
