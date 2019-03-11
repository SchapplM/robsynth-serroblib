% Calculate potential energy for
% S6PRRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRPR4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRRPR4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_energypot_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR4_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR4_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:15:32
% EndTime: 2019-03-08 23:15:32
% DurationCPUTime: 0.45s
% Computational Cost: add. (356->124), mult. (604->146), div. (0->0), fcn. (717->14), ass. (0->52)
t36 = sin(qJ(4));
t66 = t36 * pkin(4) + pkin(8);
t33 = sin(pkin(6));
t65 = pkin(7) * t33;
t63 = rSges(4,3) + pkin(8);
t62 = pkin(9) + rSges(5,3);
t39 = cos(qJ(4));
t22 = t39 * pkin(4) + pkin(3);
t61 = cos(qJ(3));
t37 = sin(qJ(3));
t60 = t33 * t37;
t38 = sin(qJ(2));
t59 = t33 * t38;
t40 = cos(qJ(2));
t58 = t33 * t40;
t35 = -qJ(5) - pkin(9);
t57 = pkin(10) - t35 + rSges(7,3);
t56 = -t35 + rSges(6,3);
t55 = cos(pkin(6));
t31 = qJ(4) + pkin(12);
t32 = sin(pkin(11));
t54 = t32 * pkin(1) + r_base(2);
t53 = qJ(1) + r_base(3);
t52 = t33 * t61;
t51 = t38 * t55;
t50 = t40 * t55;
t34 = cos(pkin(11));
t49 = t34 * pkin(1) + t32 * t65 + r_base(1);
t48 = t55 * pkin(7) + t53;
t10 = -t32 * t51 + t34 * t40;
t47 = t10 * pkin(2) + t49;
t46 = pkin(2) * t59 + t48;
t23 = sin(t31);
t24 = cos(t31);
t45 = t24 * rSges(6,1) - t23 * rSges(6,2) + t22;
t25 = qJ(6) + t31;
t20 = sin(t25);
t21 = cos(t25);
t44 = t21 * rSges(7,1) - t20 * rSges(7,2) + pkin(5) * t24 + t22;
t43 = t20 * rSges(7,1) + t21 * rSges(7,2) + pkin(5) * t23 + t66;
t8 = t32 * t40 + t34 * t51;
t42 = t8 * pkin(2) - t34 * t65 + t54;
t41 = t23 * rSges(6,1) + t24 * rSges(6,2) + t66;
t12 = t55 * t37 + t38 * t52;
t11 = t37 * t59 - t55 * t61;
t9 = t32 * t50 + t34 * t38;
t7 = t32 * t38 - t34 * t50;
t4 = t10 * t61 + t32 * t60;
t3 = t10 * t37 - t32 * t52;
t2 = -t34 * t60 + t8 * t61;
t1 = t34 * t52 + t8 * t37;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t34 * rSges(2,1) - t32 * rSges(2,2) + r_base(1)) + g(2) * (t32 * rSges(2,1) + t34 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t53)) - m(3) * (g(1) * (t10 * rSges(3,1) - t9 * rSges(3,2) + t49) + g(2) * (t8 * rSges(3,1) - t7 * rSges(3,2) + t54) + g(3) * (t55 * rSges(3,3) + t48) + (g(1) * rSges(3,3) * t32 + g(3) * (rSges(3,1) * t38 + rSges(3,2) * t40) + g(2) * (-rSges(3,3) - pkin(7)) * t34) * t33) - m(4) * (g(1) * (t4 * rSges(4,1) - t3 * rSges(4,2) + t63 * t9 + t47) + g(2) * (t2 * rSges(4,1) - t1 * rSges(4,2) + t63 * t7 + t42) + g(3) * (t12 * rSges(4,1) - t11 * rSges(4,2) - t63 * t58 + t46)) - m(5) * (g(1) * (t4 * pkin(3) + t9 * pkin(8) + (t9 * t36 + t4 * t39) * rSges(5,1) + (-t4 * t36 + t9 * t39) * rSges(5,2) + t62 * t3 + t47) + g(2) * (t2 * pkin(3) + t7 * pkin(8) + (t2 * t39 + t7 * t36) * rSges(5,1) + (-t2 * t36 + t7 * t39) * rSges(5,2) + t62 * t1 + t42) + g(3) * (t12 * pkin(3) - pkin(8) * t58 + (t12 * t39 - t36 * t58) * rSges(5,1) + (-t12 * t36 - t39 * t58) * rSges(5,2) + t62 * t11 + t46)) - m(6) * (g(1) * (t56 * t3 + t45 * t4 + t41 * t9 + t47) + g(2) * (t56 * t1 + t45 * t2 + t41 * t7 + t42) + g(3) * (t56 * t11 + t45 * t12 - t41 * t58 + t46)) - m(7) * (g(1) * (t57 * t3 + t44 * t4 + t43 * t9 + t47) + g(2) * (t57 * t1 + t44 * t2 + t43 * t7 + t42) + g(3) * (t57 * t11 + t44 * t12 - t43 * t58 + t46));
U  = t5;
