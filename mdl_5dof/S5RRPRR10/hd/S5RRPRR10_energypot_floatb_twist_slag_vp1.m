% Calculate potential energy for
% S5RRPRR10
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR10_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPRR10_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR10_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR10_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR10_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:23:41
% EndTime: 2019-12-31 20:23:42
% DurationCPUTime: 0.51s
% Computational Cost: add. (271->109), mult. (549->140), div. (0->0), fcn. (668->12), ass. (0->48)
t31 = sin(pkin(10));
t33 = cos(pkin(10));
t37 = sin(qJ(2));
t41 = cos(qJ(2));
t21 = -t37 * t31 + t41 * t33;
t64 = pkin(2) * t37;
t62 = rSges(5,3) + pkin(8);
t61 = pkin(9) + rSges(6,3);
t32 = sin(pkin(5));
t38 = sin(qJ(1));
t59 = t38 * t32;
t58 = t38 * t37;
t57 = t38 * t41;
t42 = cos(qJ(1));
t55 = t42 * t32;
t54 = t42 * t37;
t53 = t42 * t41;
t52 = pkin(6) + r_base(3);
t28 = pkin(2) * t41 + pkin(1);
t51 = t42 * t28 + r_base(1);
t34 = cos(pkin(5));
t50 = t34 * pkin(7) + t52;
t19 = t34 * t64 + (-pkin(7) - qJ(3)) * t32;
t49 = t42 * t19 + t38 * t28 + r_base(2);
t47 = t31 * t41 + t33 * t37;
t18 = t47 * t34;
t8 = t18 * t42 + t21 * t38;
t48 = t8 * pkin(3) + t49;
t46 = t34 * qJ(3) + t32 * t64 + t50;
t10 = -t18 * t38 + t21 * t42;
t45 = t10 * pkin(3) - t38 * t19 + t51;
t17 = t47 * t32;
t44 = t17 * pkin(3) + t46;
t43 = t21 * t34;
t40 = cos(qJ(4));
t39 = cos(qJ(5));
t36 = sin(qJ(4));
t35 = sin(qJ(5));
t16 = t21 * t32;
t12 = t17 * t40 + t34 * t36;
t11 = t17 * t36 - t34 * t40;
t9 = -t38 * t43 - t42 * t47;
t7 = -t38 * t47 + t42 * t43;
t4 = t10 * t40 + t36 * t59;
t3 = t10 * t36 - t40 * t59;
t2 = -t36 * t55 + t40 * t8;
t1 = t36 * t8 + t40 * t55;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t42 - rSges(2,2) * t38 + r_base(1)) + g(2) * (rSges(2,1) * t38 + rSges(2,2) * t42 + r_base(2)) + g(3) * (rSges(2,3) + t52)) - m(3) * (g(1) * (t42 * pkin(1) + r_base(1) + (-t34 * t58 + t53) * rSges(3,1) + (-t34 * t57 - t54) * rSges(3,2)) + g(2) * (t38 * pkin(1) + r_base(2) + (t34 * t54 + t57) * rSges(3,1) + (t34 * t53 - t58) * rSges(3,2)) + g(3) * (t34 * rSges(3,3) + t50) + (g(3) * (rSges(3,1) * t37 + rSges(3,2) * t41) + (g(1) * t38 - g(2) * t42) * (rSges(3,3) + pkin(7))) * t32) - m(4) * (g(1) * (t10 * rSges(4,1) + t9 * rSges(4,2) + (rSges(4,3) * t32 - t19) * t38 + t51) + g(2) * (rSges(4,1) * t8 + rSges(4,2) * t7 - rSges(4,3) * t55 + t49) + g(3) * (rSges(4,1) * t17 + rSges(4,2) * t16 + rSges(4,3) * t34 + t46)) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) - t62 * t9 + t45) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) - t62 * t7 + t48) + g(3) * (t12 * rSges(5,1) - t11 * rSges(5,2) - t16 * t62 + t44)) - m(6) * (g(1) * (t4 * pkin(4) - t9 * pkin(8) + (-t35 * t9 + t39 * t4) * rSges(6,1) + (-t35 * t4 - t39 * t9) * rSges(6,2) + t61 * t3 + t45) + g(2) * (t2 * pkin(4) - t7 * pkin(8) + (t2 * t39 - t35 * t7) * rSges(6,1) + (-t2 * t35 - t39 * t7) * rSges(6,2) + t61 * t1 + t48) + g(3) * (t12 * pkin(4) - t16 * pkin(8) + (t12 * t39 - t16 * t35) * rSges(6,1) + (-t12 * t35 - t16 * t39) * rSges(6,2) + t61 * t11 + t44));
U = t5;
