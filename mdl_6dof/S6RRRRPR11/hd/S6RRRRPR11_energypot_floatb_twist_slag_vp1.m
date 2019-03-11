% Calculate potential energy for
% S6RRRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPR11_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRPR11_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR11_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_energypot_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR11_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR11_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:12:14
% EndTime: 2019-03-09 23:12:14
% DurationCPUTime: 0.45s
% Computational Cost: add. (356->124), mult. (604->146), div. (0->0), fcn. (717->14), ass. (0->52)
t34 = sin(qJ(4));
t66 = t34 * pkin(4) + pkin(9);
t64 = rSges(4,3) + pkin(9);
t63 = pkin(10) + rSges(5,3);
t38 = cos(qJ(4));
t22 = t38 * pkin(4) + pkin(3);
t62 = cos(qJ(3));
t32 = sin(pkin(6));
t36 = sin(qJ(2));
t61 = t32 * t36;
t37 = sin(qJ(1));
t60 = t32 * t37;
t39 = cos(qJ(2));
t59 = t32 * t39;
t40 = cos(qJ(1));
t58 = t32 * t40;
t33 = -qJ(5) - pkin(10);
t57 = pkin(11) - t33 + rSges(7,3);
t56 = -t33 + rSges(6,3);
t55 = cos(pkin(6));
t54 = pkin(7) + r_base(3);
t31 = qJ(4) + pkin(12);
t53 = t37 * pkin(1) + r_base(2);
t52 = t32 * t62;
t51 = t55 * pkin(8) + t54;
t50 = t37 * t55;
t49 = t40 * t55;
t48 = t40 * pkin(1) + pkin(8) * t60 + r_base(1);
t47 = pkin(2) * t61 + t51;
t12 = -t36 * t50 + t40 * t39;
t46 = t12 * pkin(2) + t48;
t23 = sin(t31);
t24 = cos(t31);
t45 = t24 * rSges(6,1) - t23 * rSges(6,2) + t22;
t25 = qJ(6) + t31;
t20 = sin(t25);
t21 = cos(t25);
t44 = t21 * rSges(7,1) - t20 * rSges(7,2) + pkin(5) * t24 + t22;
t43 = t20 * rSges(7,1) + t21 * rSges(7,2) + pkin(5) * t23 + t66;
t10 = t36 * t49 + t37 * t39;
t42 = t10 * pkin(2) - pkin(8) * t58 + t53;
t41 = t23 * rSges(6,1) + t24 * rSges(6,2) + t66;
t35 = sin(qJ(3));
t11 = t40 * t36 + t39 * t50;
t9 = t37 * t36 - t39 * t49;
t8 = t55 * t35 + t36 * t52;
t7 = t35 * t61 - t55 * t62;
t4 = t12 * t62 + t35 * t60;
t3 = t12 * t35 - t37 * t52;
t2 = t10 * t62 - t35 * t58;
t1 = t10 * t35 + t40 * t52;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t40 * rSges(2,1) - t37 * rSges(2,2) + r_base(1)) + g(2) * (t37 * rSges(2,1) + t40 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t54)) - m(3) * (g(1) * (t12 * rSges(3,1) - t11 * rSges(3,2) + t48) + g(2) * (t10 * rSges(3,1) - t9 * rSges(3,2) + t53) + g(3) * (t55 * rSges(3,3) + t51) + (g(1) * rSges(3,3) * t37 + g(3) * (rSges(3,1) * t36 + rSges(3,2) * t39) + g(2) * (-rSges(3,3) - pkin(8)) * t40) * t32) - m(4) * (g(1) * (t4 * rSges(4,1) - t3 * rSges(4,2) + t64 * t11 + t46) + g(2) * (t2 * rSges(4,1) - t1 * rSges(4,2) + t64 * t9 + t42) + g(3) * (t8 * rSges(4,1) - t7 * rSges(4,2) - t64 * t59 + t47)) - m(5) * (g(1) * (t4 * pkin(3) + t11 * pkin(9) + (t11 * t34 + t4 * t38) * rSges(5,1) + (t11 * t38 - t4 * t34) * rSges(5,2) + t63 * t3 + t46) + g(2) * (t2 * pkin(3) + t9 * pkin(9) + (t2 * t38 + t9 * t34) * rSges(5,1) + (-t2 * t34 + t9 * t38) * rSges(5,2) + t63 * t1 + t42) + g(3) * (t8 * pkin(3) - pkin(9) * t59 + (-t34 * t59 + t8 * t38) * rSges(5,1) + (-t8 * t34 - t38 * t59) * rSges(5,2) + t63 * t7 + t47)) - m(6) * (g(1) * (t41 * t11 + t56 * t3 + t45 * t4 + t46) + g(2) * (t56 * t1 + t45 * t2 + t41 * t9 + t42) + g(3) * (-t41 * t59 + t45 * t8 + t56 * t7 + t47)) - m(7) * (g(1) * (t43 * t11 + t57 * t3 + t44 * t4 + t46) + g(2) * (t57 * t1 + t44 * t2 + t43 * t9 + t42) + g(3) * (-t43 * t59 + t44 * t8 + t57 * t7 + t47));
U  = t5;
