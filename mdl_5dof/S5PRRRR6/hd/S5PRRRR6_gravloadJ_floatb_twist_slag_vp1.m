% Calculate Gravitation load on the joints for
% S5PRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR6_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:09:34
% EndTime: 2019-12-05 17:09:35
% DurationCPUTime: 0.38s
% Computational Cost: add. (260->79), mult. (292->116), div. (0->0), fcn. (259->10), ass. (0->46)
t76 = rSges(5,3) + pkin(7);
t35 = cos(qJ(4));
t50 = -t35 * rSges(5,1) - pkin(3);
t29 = qJ(4) + qJ(5);
t26 = cos(t29);
t49 = -rSges(6,1) * t26 - t35 * pkin(4) - pkin(3);
t31 = sin(pkin(9));
t32 = cos(pkin(9));
t75 = g(1) * t32 + g(2) * t31;
t30 = qJ(2) + qJ(3);
t27 = cos(t30);
t59 = t32 * t26;
t24 = sin(t29);
t60 = t32 * t24;
t63 = t31 * t26;
t64 = t31 * t24;
t74 = (-t27 * t64 - t59) * rSges(6,1) + (-t27 * t63 + t60) * rSges(6,2);
t73 = (-t27 * t60 + t63) * rSges(6,1) + (-t27 * t59 - t64) * rSges(6,2);
t34 = sin(qJ(2));
t72 = pkin(2) * t34;
t25 = sin(t30);
t69 = g(3) * t25;
t67 = rSges(6,2) * t24;
t66 = t27 * t31;
t65 = t27 * t32;
t33 = sin(qJ(4));
t62 = t31 * t33;
t61 = t31 * t35;
t58 = t32 * t33;
t57 = t32 * t35;
t56 = t33 * rSges(5,2);
t52 = t25 * t56;
t54 = t31 * t52 + t76 * t66;
t53 = t32 * t52 + t76 * t65;
t51 = t25 * t67;
t48 = t27 * rSges(4,1) - t25 * rSges(4,2);
t46 = -rSges(6,1) * t24 - rSges(6,2) * t26;
t45 = -t27 * t58 + t61;
t44 = -t27 * t62 - t57;
t43 = t76 * t25 + (-t56 - t50) * t27;
t41 = g(1) * (rSges(6,3) * t65 + t32 * t51) + g(2) * (rSges(6,3) * t66 + t31 * t51);
t37 = -pkin(8) - pkin(7);
t40 = (rSges(6,3) - t37) * t25 + (-t67 - t49) * t27;
t36 = cos(qJ(2));
t28 = t36 * pkin(2);
t1 = [(-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3), -m(3) * (g(3) * (t36 * rSges(3,1) - t34 * rSges(3,2)) + t75 * (-rSges(3,1) * t34 - rSges(3,2) * t36)) - m(4) * (g(3) * (t28 + t48) + t75 * (-rSges(4,1) * t25 - rSges(4,2) * t27 - t72)) - m(5) * (g(1) * (-t32 * t72 + t53) + g(2) * (-t31 * t72 + t54) + g(3) * (t28 + t43) + t75 * t25 * t50) - m(6) * (g(3) * (t28 + t40) + t41 + t75 * (t49 * t25 - t27 * t37 - t72)), -m(4) * g(3) * t48 - m(5) * (g(1) * t53 + g(2) * t54 + g(3) * t43) - m(6) * (g(3) * t40 + t41) + t75 * ((m(4) * rSges(4,2) + m(6) * t37) * t27 + (m(4) * rSges(4,1) - m(5) * t50 - m(6) * t49) * t25), -m(5) * (g(1) * (t45 * rSges(5,1) + (-t27 * t57 - t62) * rSges(5,2)) + g(2) * (t44 * rSges(5,1) + (-t27 * t61 + t58) * rSges(5,2))) - m(6) * (g(1) * (t45 * pkin(4) + t73) + g(2) * (t44 * pkin(4) + t74)) + (-m(5) * (-rSges(5,1) * t33 - rSges(5,2) * t35) - m(6) * (-pkin(4) * t33 + t46)) * t69, -m(6) * (g(1) * t73 + g(2) * t74 + t46 * t69)];
taug = t1(:);
