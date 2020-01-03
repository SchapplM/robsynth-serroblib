% Calculate potential energy for
% S5RRRRP11
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRP11_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRRP11_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP11_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP11_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP11_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:14:20
% EndTime: 2019-12-31 22:14:20
% DurationCPUTime: 0.30s
% Computational Cost: add. (240->98), mult. (481->119), div. (0->0), fcn. (575->10), ass. (0->51)
t65 = rSges(6,1) + pkin(4);
t64 = rSges(6,2) + pkin(9);
t63 = rSges(4,3) + pkin(8);
t62 = rSges(5,3) + pkin(9);
t32 = sin(pkin(5));
t36 = sin(qJ(2));
t61 = t32 * t36;
t37 = sin(qJ(1));
t60 = t32 * t37;
t39 = cos(qJ(3));
t59 = t32 * t39;
t40 = cos(qJ(2));
t58 = t32 * t40;
t41 = cos(qJ(1));
t57 = t32 * t41;
t56 = t37 * t36;
t55 = t37 * t40;
t54 = t41 * t36;
t53 = t41 * t40;
t52 = rSges(6,3) + qJ(5);
t51 = pkin(6) + r_base(3);
t50 = t37 * pkin(1) + r_base(2);
t33 = cos(pkin(5));
t49 = t33 * pkin(7) + t51;
t48 = t41 * pkin(1) + pkin(7) * t60 + r_base(1);
t47 = pkin(2) * t61 + t49;
t23 = -t33 * t56 + t53;
t46 = t23 * pkin(2) + t48;
t21 = t33 * t54 + t55;
t45 = t21 * pkin(2) - pkin(7) * t57 + t50;
t35 = sin(qJ(3));
t12 = t23 * t39 + t35 * t60;
t22 = t33 * t55 + t54;
t44 = t12 * pkin(3) + t22 * pkin(8) + t46;
t19 = t33 * t35 + t36 * t59;
t43 = t19 * pkin(3) - pkin(8) * t58 + t47;
t10 = t21 * t39 - t35 * t57;
t20 = -t33 * t53 + t56;
t42 = t10 * pkin(3) + t20 * pkin(8) + t45;
t38 = cos(qJ(4));
t34 = sin(qJ(4));
t18 = -t33 * t39 + t35 * t61;
t11 = t23 * t35 - t37 * t59;
t9 = t21 * t35 + t39 * t57;
t8 = t19 * t38 - t34 * t58;
t7 = t19 * t34 + t38 * t58;
t4 = t12 * t38 + t22 * t34;
t3 = t12 * t34 - t22 * t38;
t2 = t10 * t38 + t20 * t34;
t1 = t10 * t34 - t20 * t38;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t41 * rSges(2,1) - t37 * rSges(2,2) + r_base(1)) + g(2) * (t37 * rSges(2,1) + t41 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t51)) - m(3) * (g(1) * (t23 * rSges(3,1) - t22 * rSges(3,2) + t48) + g(2) * (t21 * rSges(3,1) - t20 * rSges(3,2) + t50) + g(3) * (t33 * rSges(3,3) + t49) + (g(1) * rSges(3,3) * t37 + g(3) * (rSges(3,1) * t36 + rSges(3,2) * t40) + g(2) * (-rSges(3,3) - pkin(7)) * t41) * t32) - m(4) * (g(1) * (t12 * rSges(4,1) - t11 * rSges(4,2) + t63 * t22 + t46) + g(2) * (t10 * rSges(4,1) - t9 * rSges(4,2) + t63 * t20 + t45) + g(3) * (t19 * rSges(4,1) - t18 * rSges(4,2) - t63 * t58 + t47)) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t62 * t11 + t44) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) + t62 * t9 + t42) + g(3) * (t8 * rSges(5,1) - t7 * rSges(5,2) + t62 * t18 + t43)) - m(6) * (g(1) * (t64 * t11 + t52 * t3 + t65 * t4 + t44) + g(2) * (t52 * t1 + t65 * t2 + t64 * t9 + t42) + g(3) * (t64 * t18 + t52 * t7 + t65 * t8 + t43));
U = t5;
