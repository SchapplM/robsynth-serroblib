% Calculate Gravitation load on the joints for
% S5RRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR12_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR12_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR12_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR12_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:29:03
% EndTime: 2019-12-31 20:29:05
% DurationCPUTime: 0.65s
% Computational Cost: add. (227->111), mult. (507->151), div. (0->0), fcn. (531->8), ass. (0->48)
t26 = sin(qJ(5));
t29 = cos(qJ(5));
t28 = sin(qJ(1));
t27 = sin(qJ(2));
t30 = cos(qJ(2));
t57 = sin(qJ(4));
t58 = cos(qJ(4));
t9 = t27 * t58 - t30 * t57;
t4 = t9 * t28;
t31 = cos(qJ(1));
t46 = t31 * t57;
t47 = t31 * t58;
t7 = -t27 * t47 + t30 * t46;
t8 = t27 * t57 + t30 * t58;
t71 = (g(1) * t7 - g(2) * t4 + g(3) * t8) * (rSges(6,1) * t29 - rSges(6,2) * t26 + pkin(4));
t70 = g(3) * t9;
t69 = rSges(6,3) + pkin(8);
t20 = t27 * qJ(3);
t53 = t30 * pkin(2) + t20;
t68 = g(1) * t31 + g(2) * t28;
t67 = t68 * t27;
t22 = t30 * pkin(3);
t56 = rSges(4,1) * t30;
t55 = t27 * t31;
t54 = t30 * t31;
t52 = t31 * pkin(1) + t28 * pkin(6);
t51 = qJ(3) * t30;
t50 = t22 + t53;
t24 = t31 * pkin(6);
t45 = -t31 * pkin(7) + t24;
t44 = pkin(2) * t54 + t31 * t20 + t52;
t43 = pkin(3) * t54 + t44;
t5 = t8 * t28;
t42 = rSges(5,1) * t4 - rSges(5,2) * t5;
t6 = -t27 * t46 - t30 * t47;
t41 = -t7 * rSges(5,1) + t6 * rSges(5,2);
t40 = -rSges(5,1) * t8 - rSges(5,2) * t9;
t39 = t5 * t26 - t29 * t31;
t38 = -t26 * t31 - t5 * t29;
t37 = t30 * rSges(3,1) - rSges(3,2) * t27;
t34 = -pkin(1) - t53;
t33 = g(1) * (t34 - t22);
t32 = (-pkin(2) - pkin(3)) * t67;
t16 = t31 * t51;
t14 = t28 * t51;
t2 = -t26 * t28 - t29 * t6;
t1 = t26 * t6 - t28 * t29;
t3 = [-m(2) * (g(1) * (-t28 * rSges(2,1) - rSges(2,2) * t31) + g(2) * (rSges(2,1) * t31 - t28 * rSges(2,2))) - m(3) * (g(1) * (t31 * rSges(3,3) + t24) + g(2) * (rSges(3,1) * t54 - rSges(3,2) * t55 + t52) + (g(1) * (-pkin(1) - t37) + g(2) * rSges(3,3)) * t28) - m(4) * (g(1) * (rSges(4,2) * t31 + t24) + g(2) * (rSges(4,1) * t54 + rSges(4,3) * t55 + t44) + (g(1) * (-rSges(4,3) * t27 + t34 - t56) + g(2) * rSges(4,2)) * t28) - m(5) * (g(1) * (-t5 * rSges(5,1) - t4 * rSges(5,2) - rSges(5,3) * t31 + t45) + g(2) * (-rSges(5,1) * t6 - rSges(5,2) * t7 + t43) + (t33 + g(2) * (-rSges(5,3) - pkin(7))) * t28) - m(6) * (g(1) * (t38 * rSges(6,1) + t39 * rSges(6,2) - t5 * pkin(4) + t69 * t4 + t45) + g(2) * (rSges(6,1) * t2 + rSges(6,2) * t1 - pkin(4) * t6 + t69 * t7 + t43) + (-g(2) * pkin(7) + t33) * t28), -m(3) * (g(3) * t37 + t68 * (-rSges(3,1) * t27 - rSges(3,2) * t30)) - m(4) * (g(1) * (rSges(4,3) * t54 + t16) + g(2) * (t28 * t30 * rSges(4,3) + t14) + g(3) * (t53 + t56) + (g(3) * rSges(4,3) + t68 * (-rSges(4,1) - pkin(2))) * t27) - m(5) * (g(1) * (t16 - t41) + g(2) * (t14 - t42) + g(3) * (-t40 + t50) + t32) + (-g(1) * (t69 * t6 + t16) - g(2) * (-t69 * t5 + t14) - g(3) * (-t69 * t9 + t50) - t32 - t71) * m(6), (-m(4) - m(5) - m(6)) * (-g(3) * t30 + t67), -m(5) * (g(1) * t41 + g(2) * t42 + g(3) * t40) - m(6) * ((-g(1) * t6 + g(2) * t5 + t70) * t69 - t71), -m(6) * (g(1) * (rSges(6,1) * t1 - rSges(6,2) * t2) + g(2) * (-t39 * rSges(6,1) + t38 * rSges(6,2)) + (-rSges(6,1) * t26 - rSges(6,2) * t29) * t70)];
taug = t3(:);
