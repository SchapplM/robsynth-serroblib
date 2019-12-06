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
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:35:41
% EndTime: 2019-12-05 18:35:43
% DurationCPUTime: 0.35s
% Computational Cost: add. (335->93), mult. (308->124), div. (0->0), fcn. (294->10), ass. (0->51)
t36 = qJ(1) + qJ(2);
t32 = sin(t36);
t34 = cos(t36);
t41 = cos(qJ(4));
t38 = cos(pkin(9));
t39 = sin(qJ(4));
t58 = t38 * t39;
t17 = t32 * t58 + t34 * t41;
t57 = t38 * t41;
t60 = t34 * t39;
t18 = t32 * t57 - t60;
t73 = -t18 * rSges(5,1) + t17 * rSges(5,2);
t37 = sin(pkin(9));
t63 = rSges(4,2) * t37;
t72 = t34 * rSges(4,3) + t32 * t63;
t71 = -pkin(3) * t38 - pkin(2) + (-rSges(5,3) - pkin(7)) * t37;
t35 = qJ(4) + qJ(5);
t31 = sin(t35);
t33 = cos(t35);
t62 = t32 * t38;
t10 = -t31 * t34 + t33 * t62;
t9 = t31 * t62 + t33 * t34;
t70 = t9 * rSges(6,1) + t10 * rSges(6,2);
t61 = t34 * t38;
t11 = t31 * t61 - t32 * t33;
t12 = -t31 * t32 - t33 * t61;
t69 = -t11 * rSges(6,1) + t12 * rSges(6,2);
t68 = pkin(4) * t39;
t67 = g(1) * t37;
t66 = g(2) * t34;
t40 = sin(qJ(1));
t65 = t40 * pkin(1);
t42 = cos(qJ(1));
t64 = t42 * pkin(1);
t59 = t37 * (-pkin(8) - pkin(7));
t19 = -t32 * t41 + t34 * t58;
t20 = -t32 * t39 - t34 * t57;
t56 = t20 * rSges(5,1) + t19 * rSges(5,2);
t55 = t12 * rSges(6,1) + t11 * rSges(6,2) + t34 * t59;
t27 = t34 * qJ(3);
t53 = t27 - t65;
t52 = -rSges(3,1) * t34 + t32 * rSges(3,2);
t51 = -t10 * rSges(6,1) + t9 * rSges(6,2) + pkin(4) * t60 + t32 * t59 + t27;
t50 = -rSges(3,1) * t32 - rSges(3,2) * t34;
t49 = -rSges(6,1) * t31 - rSges(6,2) * t33;
t48 = -rSges(6,3) * t37 - (pkin(4) * t41 + pkin(3)) * t38 - pkin(2);
t47 = -rSges(4,1) * t61 + (-pkin(2) + t63) * t34;
t46 = (g(2) * (-rSges(4,3) - qJ(3)) + g(3) * (-rSges(4,1) * t38 - pkin(2))) * t32;
t45 = (g(2) * (-qJ(3) - t68) + g(3) * t48) * t32 + t48 * t66;
t44 = (-g(2) * qJ(3) + g(3) * t71) * t32 + t71 * t66;
t1 = [-m(2) * (g(2) * (-rSges(2,1) * t42 + rSges(2,2) * t40) + g(3) * (-rSges(2,1) * t40 - rSges(2,2) * t42)) - m(3) * (g(2) * (t52 - t64) + g(3) * (t50 - t65)) - m(4) * (g(2) * (t47 - t64) + g(3) * (t53 + t72) + t46) - m(5) * (g(2) * (t56 - t64) + g(3) * (t53 + t73) + t44) - m(6) * (g(2) * (t55 - t64) + g(3) * (t51 - t65) + t45), -m(3) * (g(2) * t52 + g(3) * t50) - m(4) * (g(2) * t47 + g(3) * (t27 + t72) + t46) - m(5) * (g(2) * t56 + g(3) * (t27 + t73) + t44) - m(6) * (g(2) * t55 + g(3) * t51 + t45), (-m(4) - m(5) - m(6)) * (g(3) * t32 + t66), -m(5) * (g(2) * (rSges(5,1) * t17 + rSges(5,2) * t18) + g(3) * (-rSges(5,1) * t19 + rSges(5,2) * t20)) - m(6) * (g(2) * (pkin(4) * t17 + t70) + g(3) * (-pkin(4) * t19 + t69)) + (-m(5) * (-rSges(5,1) * t39 - rSges(5,2) * t41) - m(6) * (t49 - t68)) * t67, -m(6) * (g(2) * t70 + g(3) * t69 + t49 * t67)];
taug = t1(:);
