% Calculate Gravitation load on the joints for
% S5RRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR11_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR11_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR11_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR11_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:32:26
% EndTime: 2019-12-31 21:32:29
% DurationCPUTime: 0.95s
% Computational Cost: add. (264->144), mult. (599->195), div. (0->0), fcn. (636->8), ass. (0->52)
t30 = sin(qJ(3));
t34 = cos(qJ(3));
t36 = cos(qJ(1));
t32 = sin(qJ(1));
t35 = cos(qJ(2));
t60 = t32 * t35;
t11 = t30 * t60 + t34 * t36;
t57 = t36 * t30;
t59 = t34 * t35;
t12 = t32 * t59 - t57;
t29 = sin(qJ(5));
t33 = cos(qJ(5));
t42 = t11 * t29 + t12 * t33;
t73 = t11 * t33 - t12 * t29;
t74 = rSges(6,1) * t73 - t42 * rSges(6,2);
t67 = g(2) * t32;
t71 = g(1) * t36 + t67;
t70 = -pkin(3) - pkin(4);
t69 = g(1) * t32;
t31 = sin(qJ(2));
t66 = g(3) * t31;
t26 = t35 * pkin(2);
t65 = -rSges(5,1) - pkin(3);
t64 = -pkin(8) - rSges(6,3);
t62 = t30 * t35;
t61 = t31 * t36;
t58 = t35 * t36;
t56 = t31 * pkin(7) + t26;
t55 = t36 * pkin(1) + t32 * pkin(6);
t54 = rSges(5,3) + qJ(4);
t53 = -pkin(1) - t26;
t52 = t64 * t36;
t50 = pkin(3) * t59 + qJ(4) * t62 + t56;
t49 = pkin(2) * t58 + pkin(7) * t61 + t55;
t27 = t36 * pkin(6);
t48 = -pkin(3) * t12 - t11 * qJ(4) + t27;
t14 = t30 * t32 + t34 * t58;
t47 = pkin(3) * t14 + t49;
t13 = -t32 * t34 + t35 * t57;
t2 = t13 * t33 - t14 * t29;
t3 = t13 * t29 + t14 * t33;
t46 = rSges(6,1) * t2 - rSges(6,2) * t3;
t40 = t29 * t30 + t33 * t34;
t41 = t29 * t34 - t30 * t33;
t45 = (-rSges(6,1) * t41 - rSges(6,2) * t40) * t31;
t44 = rSges(3,1) * t35 - rSges(3,2) * t31;
t20 = pkin(7) * t58;
t17 = pkin(7) * t60;
t15 = t31 * t34 * qJ(4);
t9 = t13 * pkin(3);
t7 = t11 * pkin(3);
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t32 - rSges(2,2) * t36) + g(2) * (rSges(2,1) * t36 - rSges(2,2) * t32)) - m(3) * (g(1) * (rSges(3,3) * t36 + t27) + g(2) * (rSges(3,1) * t58 - rSges(3,2) * t61 + t55) + (g(1) * (-pkin(1) - t44) + g(2) * rSges(3,3)) * t32) - m(4) * (g(1) * (-rSges(4,1) * t12 + rSges(4,2) * t11 + t27) + g(2) * (rSges(4,1) * t14 - rSges(4,2) * t13 + rSges(4,3) * t61 + t49) + ((-rSges(4,3) - pkin(7)) * t31 + t53) * t69) - m(5) * (g(1) * (-rSges(5,1) * t12 - rSges(5,3) * t11 + t48) + g(2) * (t14 * rSges(5,1) + rSges(5,2) * t61 + t13 * t54 + t47) + ((-rSges(5,2) - pkin(7)) * t31 + t53) * t69) - m(6) * (g(1) * (-rSges(6,1) * t42 - rSges(6,2) * t73 - t12 * pkin(4) + t48) + g(2) * (t3 * rSges(6,1) + t2 * rSges(6,2) + t14 * pkin(4) + t13 * qJ(4) + t31 * t52 + t47) + ((-pkin(7) - t64) * t31 + t53) * t69), -m(3) * (g(3) * t44 + t71 * (-rSges(3,1) * t31 - rSges(3,2) * t35)) - m(4) * (g(1) * (rSges(4,3) * t58 + t20) + g(2) * (rSges(4,3) * t60 + t17) + g(3) * (rSges(4,1) * t59 - rSges(4,2) * t62 + t56) + (g(3) * rSges(4,3) + t71 * (-rSges(4,1) * t34 + rSges(4,2) * t30 - pkin(2))) * t31) - m(5) * (g(1) * (rSges(5,2) * t58 + t20) + g(2) * (rSges(5,2) * t60 + t17) + g(3) * (rSges(5,1) * t59 + rSges(5,3) * t62 + t50) + (g(3) * rSges(5,2) + t71 * (-t30 * t54 + t34 * t65 - pkin(2))) * t31) - m(6) * (g(1) * t20 + g(2) * t17 + g(3) * t50 + (g(3) * (rSges(6,1) * t40 - rSges(6,2) * t41 + t34 * pkin(4)) + g(1) * t52 + t64 * t67) * t35 + (g(3) * t64 + t71 * (-pkin(2) + (-rSges(6,1) * t29 - rSges(6,2) * t33 - qJ(4)) * t30 + (-rSges(6,1) * t33 + rSges(6,2) * t29 + t70) * t34)) * t31), -m(4) * (g(1) * (-rSges(4,1) * t13 - rSges(4,2) * t14) + g(2) * (-rSges(4,1) * t11 - rSges(4,2) * t12)) - m(5) * (g(1) * (-rSges(5,1) * t13 + t14 * t54 - t9) + g(2) * (-rSges(5,1) * t11 + t12 * t54 - t7) + g(3) * t15) - m(6) * (g(1) * (-t13 * pkin(4) + t14 * qJ(4) - t46 - t9) + g(2) * (-t11 * pkin(4) + t12 * qJ(4) - t7 - t74) + g(3) * (t15 - t45)) + ((m(4) * rSges(4,2) - m(5) * rSges(5,3)) * t34 + (m(4) * rSges(4,1) - m(5) * t65 - m(6) * t70) * t30) * t66, (-m(5) - m(6)) * (g(1) * t13 + g(2) * t11 + t30 * t66), -m(6) * (g(1) * t46 + g(2) * t74 + g(3) * t45)];
taug = t1(:);
