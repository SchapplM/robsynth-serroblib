% Calculate Gravitation load on the joints for
% S5RRRRR9
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
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR9_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR9_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR9_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR9_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:27:44
% EndTime: 2019-12-31 22:27:46
% DurationCPUTime: 0.62s
% Computational Cost: add. (370->126), mult. (458->179), div. (0->0), fcn. (439->10), ass. (0->62)
t43 = -pkin(8) - pkin(7);
t56 = rSges(5,3) - t43;
t55 = rSges(6,3) + pkin(9) - t43;
t39 = sin(qJ(1));
t42 = cos(qJ(1));
t77 = g(1) * t42 + g(2) * t39;
t38 = sin(qJ(2));
t41 = cos(qJ(2));
t67 = rSges(4,3) + pkin(7);
t76 = t41 * pkin(2) + t67 * t38;
t36 = qJ(3) + qJ(4);
t30 = qJ(5) + t36;
t25 = sin(t30);
t26 = cos(t30);
t61 = t42 * t26;
t64 = t39 * t41;
t5 = t25 * t64 + t61;
t62 = t42 * t25;
t6 = -t26 * t64 + t62;
t75 = -t5 * rSges(6,1) + t6 * rSges(6,2);
t7 = t39 * t26 - t41 * t62;
t8 = t39 * t25 + t41 * t61;
t74 = t7 * rSges(6,1) - t8 * rSges(6,2);
t28 = sin(t36);
t73 = pkin(4) * t28;
t70 = g(3) * t38;
t37 = sin(qJ(3));
t69 = t37 * pkin(3);
t29 = cos(t36);
t59 = t42 * t29;
t13 = t28 * t64 + t59;
t60 = t42 * t28;
t14 = -t29 * t64 + t60;
t66 = -t13 * rSges(5,1) + t14 * rSges(5,2);
t65 = t38 * rSges(3,2);
t63 = t41 * t42;
t58 = t42 * t37;
t40 = cos(qJ(3));
t57 = t42 * t40;
t15 = t39 * t29 - t41 * t60;
t16 = t39 * t28 + t41 * t59;
t54 = t15 * rSges(5,1) - t16 * rSges(5,2);
t32 = t40 * pkin(3);
t23 = pkin(4) * t29 + t32;
t53 = t42 * pkin(1) + t39 * pkin(6);
t52 = t41 * rSges(3,1) - t65;
t50 = -rSges(5,1) * t28 - rSges(5,2) * t29;
t49 = -rSges(6,1) * t25 - rSges(6,2) * t26;
t48 = rSges(4,1) * t40 - rSges(4,2) * t37 + pkin(2);
t19 = t39 * t40 - t41 * t58;
t17 = t37 * t64 + t57;
t27 = t32 + pkin(2);
t47 = rSges(5,1) * t29 - rSges(5,2) * t28 + t27;
t21 = pkin(2) + t23;
t46 = rSges(6,1) * t26 - rSges(6,2) * t25 + t21;
t45 = t41 * t27 + t56 * t38;
t44 = t41 * t21 + t55 * t38;
t33 = t42 * pkin(6);
t22 = t69 + t73;
t20 = t39 * t37 + t41 * t57;
t18 = -t40 * t64 + t58;
t1 = [-m(2) * (g(1) * (-t39 * rSges(2,1) - t42 * rSges(2,2)) + g(2) * (t42 * rSges(2,1) - t39 * rSges(2,2))) - m(3) * (g(1) * (t42 * rSges(3,3) + t33) + g(2) * (rSges(3,1) * t63 - t42 * t65 + t53) + (g(1) * (-pkin(1) - t52) + g(2) * rSges(3,3)) * t39) - m(4) * ((t20 * rSges(4,1) + t19 * rSges(4,2) + t76 * t42 + t53) * g(2) + (t18 * rSges(4,1) + t17 * rSges(4,2) + t33 + (-pkin(1) - t76) * t39) * g(1)) - m(5) * (g(1) * (t14 * rSges(5,1) + t13 * rSges(5,2) + t33) + g(2) * (t16 * rSges(5,1) + t15 * rSges(5,2) + t53) + (g(1) * t69 + g(2) * t45) * t42 + (g(1) * (-pkin(1) - t45) + g(2) * t69) * t39) - m(6) * (g(1) * (t6 * rSges(6,1) + t5 * rSges(6,2) + t33) + g(2) * (t8 * rSges(6,1) + t7 * rSges(6,2) + t53) + (g(1) * t22 + g(2) * t44) * t42 + (g(1) * (-pkin(1) - t44) + g(2) * t22) * t39), -m(3) * (g(3) * t52 + t77 * (-rSges(3,1) * t38 - rSges(3,2) * t41)) - m(4) * ((g(3) * t48 + t77 * t67) * t41 + (g(3) * t67 - t77 * t48) * t38) - m(5) * ((g(3) * t47 + t77 * t56) * t41 + (g(3) * t56 - t77 * t47) * t38) - m(6) * ((g(3) * t46 + t77 * t55) * t41 + (g(3) * t55 - t77 * t46) * t38), -m(4) * (g(1) * (t19 * rSges(4,1) - t20 * rSges(4,2)) + g(2) * (-t17 * rSges(4,1) + t18 * rSges(4,2))) - m(5) * (g(1) * (t19 * pkin(3) + t54) + g(2) * (-t17 * pkin(3) + t66)) - m(6) * (g(1) * (-t22 * t63 + t39 * t23 + t74) + g(2) * (-t22 * t64 - t42 * t23 + t75)) + (-m(4) * (-rSges(4,1) * t37 - rSges(4,2) * t40) - m(5) * (t50 - t69) - m(6) * (-t22 + t49)) * t70, -m(5) * (g(1) * t54 + g(2) * t66) - m(6) * (g(1) * (t15 * pkin(4) + t74) + g(2) * (-t13 * pkin(4) + t75)) + (-m(5) * t50 - m(6) * (t49 - t73)) * t70, -m(6) * (g(1) * t74 + g(2) * t75 + t49 * t70)];
taug = t1(:);
