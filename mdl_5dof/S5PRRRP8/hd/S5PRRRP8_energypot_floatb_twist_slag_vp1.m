% Calculate potential energy for
% S5PRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRP8_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRRP8_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP8_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP8_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP8_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:58:08
% EndTime: 2019-12-05 16:58:08
% DurationCPUTime: 0.31s
% Computational Cost: add. (240->98), mult. (481->121), div. (0->0), fcn. (575->10), ass. (0->49)
t33 = sin(pkin(5));
t63 = pkin(6) * t33;
t62 = rSges(6,1) + pkin(4);
t61 = rSges(6,2) + pkin(8);
t60 = rSges(4,3) + pkin(7);
t59 = rSges(5,3) + pkin(8);
t37 = sin(qJ(3));
t58 = t33 * t37;
t38 = sin(qJ(2));
t57 = t33 * t38;
t40 = cos(qJ(3));
t56 = t33 * t40;
t41 = cos(qJ(2));
t55 = t33 * t41;
t35 = cos(pkin(5));
t54 = t35 * t38;
t53 = t35 * t41;
t52 = rSges(6,3) + qJ(5);
t32 = sin(pkin(9));
t51 = t32 * pkin(1) + r_base(2);
t50 = qJ(1) + r_base(3);
t34 = cos(pkin(9));
t49 = t34 * pkin(1) + t32 * t63 + r_base(1);
t48 = t35 * pkin(6) + t50;
t21 = -t32 * t54 + t34 * t41;
t47 = t21 * pkin(2) + t49;
t46 = pkin(2) * t57 + t48;
t19 = t32 * t41 + t34 * t54;
t45 = t19 * pkin(2) - t34 * t63 + t51;
t10 = t21 * t40 + t32 * t58;
t20 = t32 * t53 + t34 * t38;
t44 = t10 * pkin(3) + t20 * pkin(7) + t47;
t23 = t35 * t37 + t38 * t56;
t43 = t23 * pkin(3) - pkin(7) * t55 + t46;
t18 = t32 * t38 - t34 * t53;
t8 = t19 * t40 - t34 * t58;
t42 = t8 * pkin(3) + t18 * pkin(7) + t45;
t39 = cos(qJ(4));
t36 = sin(qJ(4));
t22 = -t35 * t40 + t37 * t57;
t12 = t23 * t39 - t36 * t55;
t11 = t23 * t36 + t39 * t55;
t9 = t21 * t37 - t32 * t56;
t7 = t19 * t37 + t34 * t56;
t4 = t10 * t39 + t20 * t36;
t3 = t10 * t36 - t20 * t39;
t2 = t18 * t36 + t8 * t39;
t1 = -t18 * t39 + t8 * t36;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t34 * rSges(2,1) - t32 * rSges(2,2) + r_base(1)) + g(2) * (t32 * rSges(2,1) + t34 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t50)) - m(3) * (g(1) * (t21 * rSges(3,1) - t20 * rSges(3,2) + t49) + g(2) * (t19 * rSges(3,1) - t18 * rSges(3,2) + t51) + g(3) * (t35 * rSges(3,3) + t48) + (g(1) * rSges(3,3) * t32 + g(3) * (rSges(3,1) * t38 + rSges(3,2) * t41) + g(2) * (-rSges(3,3) - pkin(6)) * t34) * t33) - m(4) * (g(1) * (t10 * rSges(4,1) - t9 * rSges(4,2) + t60 * t20 + t47) + g(2) * (t8 * rSges(4,1) - t7 * rSges(4,2) + t60 * t18 + t45) + g(3) * (t23 * rSges(4,1) - t22 * rSges(4,2) - t60 * t55 + t46)) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t59 * t9 + t44) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) + t59 * t7 + t42) + g(3) * (t12 * rSges(5,1) - t11 * rSges(5,2) + t59 * t22 + t43)) - m(6) * (g(1) * (t52 * t3 + t62 * t4 + t61 * t9 + t44) + g(2) * (t52 * t1 + t62 * t2 + t61 * t7 + t42) + g(3) * (t52 * t11 + t62 * t12 + t61 * t22 + t43));
U = t5;
