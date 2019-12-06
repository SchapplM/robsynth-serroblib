% Calculate potential energy for
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRR4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRPRR4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR4_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR4_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:49:30
% EndTime: 2019-12-05 15:49:31
% DurationCPUTime: 0.52s
% Computational Cost: add. (271->109), mult. (549->142), div. (0->0), fcn. (668->12), ass. (0->45)
t31 = sin(pkin(10));
t34 = cos(pkin(10));
t39 = sin(qJ(2));
t42 = cos(qJ(2));
t21 = -t39 * t31 + t42 * t34;
t60 = rSges(5,3) + pkin(7);
t59 = pkin(8) + rSges(6,3);
t32 = sin(pkin(9));
t33 = sin(pkin(5));
t58 = t32 * t33;
t35 = cos(pkin(9));
t57 = t35 * t33;
t36 = cos(pkin(5));
t56 = t36 * t39;
t55 = t36 * t42;
t28 = t42 * pkin(2) + pkin(1);
t52 = t35 * t28 + r_base(1);
t51 = qJ(1) + r_base(3);
t19 = pkin(2) * t56 + (-pkin(6) - qJ(3)) * t33;
t50 = t35 * t19 + t32 * t28 + r_base(2);
t49 = t36 * pkin(6) + t51;
t47 = t42 * t31 + t39 * t34;
t18 = t47 * t36;
t8 = t35 * t18 + t32 * t21;
t48 = t8 * pkin(3) + t50;
t46 = t33 * t39 * pkin(2) + t36 * qJ(3) + t49;
t10 = -t32 * t18 + t35 * t21;
t45 = t10 * pkin(3) - t32 * t19 + t52;
t44 = t21 * t36;
t17 = t47 * t33;
t43 = t17 * pkin(3) + t46;
t41 = cos(qJ(4));
t40 = cos(qJ(5));
t38 = sin(qJ(4));
t37 = sin(qJ(5));
t16 = t21 * t33;
t12 = t17 * t41 + t36 * t38;
t11 = t17 * t38 - t36 * t41;
t9 = -t32 * t44 - t35 * t47;
t7 = -t32 * t47 + t35 * t44;
t4 = t10 * t41 + t38 * t58;
t3 = t10 * t38 - t41 * t58;
t2 = -t38 * t57 + t8 * t41;
t1 = t8 * t38 + t41 * t57;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t35 * rSges(2,1) - t32 * rSges(2,2) + r_base(1)) + g(2) * (t32 * rSges(2,1) + t35 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t51)) - m(3) * (g(1) * (t35 * pkin(1) + r_base(1) + (-t32 * t56 + t35 * t42) * rSges(3,1) + (-t32 * t55 - t35 * t39) * rSges(3,2)) + g(2) * (t32 * pkin(1) + r_base(2) + (t32 * t42 + t35 * t56) * rSges(3,1) + (-t32 * t39 + t35 * t55) * rSges(3,2)) + g(3) * (t36 * rSges(3,3) + t49) + (g(3) * (rSges(3,1) * t39 + rSges(3,2) * t42) + (g(1) * t32 - g(2) * t35) * (rSges(3,3) + pkin(6))) * t33) - m(4) * (g(1) * (t10 * rSges(4,1) + t9 * rSges(4,2) + (rSges(4,3) * t33 - t19) * t32 + t52) + g(2) * (t8 * rSges(4,1) + t7 * rSges(4,2) - rSges(4,3) * t57 + t50) + g(3) * (t17 * rSges(4,1) + t16 * rSges(4,2) + t36 * rSges(4,3) + t46)) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) - t60 * t9 + t45) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) - t60 * t7 + t48) + g(3) * (t12 * rSges(5,1) - t11 * rSges(5,2) - t16 * t60 + t43)) - m(6) * (g(1) * (t4 * pkin(4) - t9 * pkin(7) + (-t9 * t37 + t4 * t40) * rSges(6,1) + (-t4 * t37 - t9 * t40) * rSges(6,2) + t59 * t3 + t45) + g(2) * (t2 * pkin(4) - t7 * pkin(7) + (t2 * t40 - t7 * t37) * rSges(6,1) + (-t2 * t37 - t7 * t40) * rSges(6,2) + t59 * t1 + t48) + g(3) * (t12 * pkin(4) - t16 * pkin(7) + (t12 * t40 - t16 * t37) * rSges(6,1) + (-t12 * t37 - t16 * t40) * rSges(6,2) + t59 * t11 + t43));
U = t5;
