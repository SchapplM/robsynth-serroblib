% Calculate Gravitation load on the joints for
% S5RRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR8_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:38:08
% EndTime: 2019-12-31 19:38:09
% DurationCPUTime: 0.45s
% Computational Cost: add. (197->105), mult. (336->144), div. (0->0), fcn. (322->8), ass. (0->50)
t30 = sin(qJ(2));
t21 = t30 * qJ(3);
t32 = cos(qJ(2));
t50 = t32 * pkin(2) + t21;
t31 = sin(qJ(1));
t33 = cos(qJ(1));
t65 = g(1) * t33 + g(2) * t31;
t64 = -rSges(5,3) - qJ(4);
t63 = -m(5) - m(6);
t60 = t32 * pkin(3);
t27 = sin(pkin(8));
t57 = t30 * t27;
t56 = t30 * t33;
t55 = t32 * rSges(4,1);
t28 = cos(pkin(8));
t18 = t28 * pkin(4) + pkin(3);
t54 = t32 * t18;
t53 = t32 * t27;
t52 = t32 * t33;
t51 = -rSges(6,3) - pkin(7) - qJ(4);
t49 = t33 * pkin(1) + t31 * pkin(6);
t48 = qJ(3) * t32;
t47 = pkin(4) * t53;
t46 = t30 * (-pkin(2) - pkin(3));
t45 = pkin(2) * t52 + t33 * t21 + t49;
t26 = pkin(8) + qJ(5);
t19 = sin(t26);
t20 = cos(t26);
t39 = t32 * t19 - t30 * t20;
t1 = t39 * t31;
t38 = t30 * t19 + t32 * t20;
t2 = t38 * t31;
t44 = -t1 * rSges(6,1) - t2 * rSges(6,2);
t3 = t19 * t52 - t20 * t56;
t4 = t38 * t33;
t43 = -t3 * rSges(6,1) - t4 * rSges(6,2);
t42 = t32 * rSges(3,1) - t30 * rSges(3,2);
t40 = -rSges(6,1) * t38 + rSges(6,2) * t39;
t37 = -t30 * t28 + t53;
t36 = t32 * t28 + t57;
t35 = pkin(4) * t57 + t54;
t34 = -pkin(1) - t50;
t24 = t33 * pkin(6);
t16 = t33 * t48;
t14 = t31 * t48;
t9 = t36 * t33;
t8 = t37 * t33;
t7 = t36 * t31;
t6 = t37 * t31;
t5 = [-m(2) * (g(1) * (-t31 * rSges(2,1) - t33 * rSges(2,2)) + g(2) * (t33 * rSges(2,1) - t31 * rSges(2,2))) - m(3) * (g(1) * (t33 * rSges(3,3) + t24) + g(2) * (rSges(3,1) * t52 - rSges(3,2) * t56 + t49) + (g(1) * (-pkin(1) - t42) + g(2) * rSges(3,3)) * t31) - m(4) * (g(1) * (t33 * rSges(4,2) + t24) + g(2) * (rSges(4,1) * t52 + rSges(4,3) * t56 + t45) + (g(1) * (-t30 * rSges(4,3) + t34 - t55) + g(2) * rSges(4,2)) * t31) - m(5) * (g(1) * (-t7 * rSges(5,1) + t6 * rSges(5,2) + t64 * t33 + t24) + g(2) * (t9 * rSges(5,1) - t8 * rSges(5,2) + pkin(3) * t52 + t45) + (g(1) * (t34 - t60) + g(2) * t64) * t31) - m(6) * (g(1) * (-t2 * rSges(6,1) + t1 * rSges(6,2) + t24) + g(2) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t45) + (g(1) * t51 + g(2) * t35) * t33 + (g(1) * (t34 - t35) + g(2) * t51) * t31), -m(3) * (g(3) * t42 + t65 * (-rSges(3,1) * t30 - rSges(3,2) * t32)) - m(4) * (g(1) * (rSges(4,3) * t52 + t16) + g(2) * (t31 * t32 * rSges(4,3) + t14) + g(3) * (t50 + t55) + (g(3) * rSges(4,3) + t65 * (-rSges(4,1) - pkin(2))) * t30) - m(5) * (g(1) * (t8 * rSges(5,1) + t9 * rSges(5,2) + t33 * t46 + t16) + g(2) * (t6 * rSges(5,1) + t7 * rSges(5,2) + t31 * t46 + t14) + g(3) * (rSges(5,1) * t36 - rSges(5,2) * t37 + t50 + t60)) - m(6) * (g(1) * (t33 * t47 + t16 - t43) + g(2) * (t31 * t47 + t14 - t44) + g(3) * (-t40 + t50 + t54) + (g(3) * pkin(4) * t27 + t65 * (-pkin(2) - t18)) * t30), (-m(4) + t63) * (-g(3) * t32 + t65 * t30), t63 * (-g(1) * t31 + g(2) * t33), -m(6) * (g(1) * t43 + g(2) * t44 + g(3) * t40)];
taug = t5(:);
