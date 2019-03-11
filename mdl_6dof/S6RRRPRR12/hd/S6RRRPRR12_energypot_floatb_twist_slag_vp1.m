% Calculate potential energy for
% S6RRRPRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRR12_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPRR12_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR12_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_energypot_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR12_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR12_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:38:03
% EndTime: 2019-03-09 19:38:04
% DurationCPUTime: 0.45s
% Computational Cost: add. (356->124), mult. (604->146), div. (0->0), fcn. (717->14), ass. (0->52)
t32 = sin(pkin(12));
t66 = t32 * pkin(4) + pkin(9);
t64 = rSges(4,3) + pkin(9);
t34 = cos(pkin(12));
t22 = t34 * pkin(4) + pkin(3);
t63 = cos(qJ(3));
t33 = sin(pkin(6));
t37 = sin(qJ(2));
t62 = t33 * t37;
t38 = sin(qJ(1));
t61 = t33 * t38;
t39 = cos(qJ(2));
t60 = t33 * t39;
t40 = cos(qJ(1));
t59 = t33 * t40;
t35 = -pkin(10) - qJ(4);
t58 = pkin(11) - t35 + rSges(7,3);
t57 = -t35 + rSges(6,3);
t56 = qJ(4) + rSges(5,3);
t55 = cos(pkin(6));
t54 = pkin(7) + r_base(3);
t31 = pkin(12) + qJ(5);
t53 = t38 * pkin(1) + r_base(2);
t52 = t33 * t63;
t51 = t55 * pkin(8) + t54;
t50 = t38 * t55;
t49 = t40 * t55;
t48 = t40 * pkin(1) + pkin(8) * t61 + r_base(1);
t47 = pkin(2) * t62 + t51;
t12 = -t37 * t50 + t40 * t39;
t46 = t12 * pkin(2) + t48;
t23 = sin(t31);
t24 = cos(t31);
t45 = t24 * rSges(6,1) - t23 * rSges(6,2) + t22;
t25 = qJ(6) + t31;
t20 = sin(t25);
t21 = cos(t25);
t44 = t21 * rSges(7,1) - t20 * rSges(7,2) + pkin(5) * t24 + t22;
t43 = t20 * rSges(7,1) + t21 * rSges(7,2) + pkin(5) * t23 + t66;
t10 = t37 * t49 + t38 * t39;
t42 = t10 * pkin(2) - pkin(8) * t59 + t53;
t41 = t23 * rSges(6,1) + t24 * rSges(6,2) + t66;
t36 = sin(qJ(3));
t11 = t40 * t37 + t39 * t50;
t9 = t38 * t37 - t39 * t49;
t8 = t55 * t36 + t37 * t52;
t7 = t36 * t62 - t55 * t63;
t4 = t12 * t63 + t36 * t61;
t3 = t12 * t36 - t38 * t52;
t2 = t10 * t63 - t36 * t59;
t1 = t10 * t36 + t40 * t52;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t40 * rSges(2,1) - t38 * rSges(2,2) + r_base(1)) + g(2) * (t38 * rSges(2,1) + t40 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t54)) - m(3) * (g(1) * (t12 * rSges(3,1) - t11 * rSges(3,2) + t48) + g(2) * (t10 * rSges(3,1) - t9 * rSges(3,2) + t53) + g(3) * (t55 * rSges(3,3) + t51) + (g(1) * rSges(3,3) * t38 + g(3) * (rSges(3,1) * t37 + rSges(3,2) * t39) + g(2) * (-rSges(3,3) - pkin(8)) * t40) * t33) - m(4) * (g(1) * (t4 * rSges(4,1) - t3 * rSges(4,2) + t64 * t11 + t46) + g(2) * (t2 * rSges(4,1) - t1 * rSges(4,2) + t64 * t9 + t42) + g(3) * (t8 * rSges(4,1) - t7 * rSges(4,2) - t64 * t60 + t47)) - m(5) * (g(1) * (t4 * pkin(3) + t11 * pkin(9) + (t11 * t32 + t4 * t34) * rSges(5,1) + (t11 * t34 - t4 * t32) * rSges(5,2) + t56 * t3 + t46) + g(2) * (t2 * pkin(3) + t9 * pkin(9) + (t2 * t34 + t9 * t32) * rSges(5,1) + (-t2 * t32 + t9 * t34) * rSges(5,2) + t56 * t1 + t42) + g(3) * (t8 * pkin(3) - pkin(9) * t60 + (-t32 * t60 + t8 * t34) * rSges(5,1) + (-t8 * t32 - t34 * t60) * rSges(5,2) + t56 * t7 + t47)) - m(6) * (g(1) * (t41 * t11 + t57 * t3 + t45 * t4 + t46) + g(2) * (t57 * t1 + t45 * t2 + t41 * t9 + t42) + g(3) * (-t41 * t60 + t45 * t8 + t57 * t7 + t47)) - m(7) * (g(1) * (t43 * t11 + t58 * t3 + t44 * t4 + t46) + g(2) * (t58 * t1 + t44 * t2 + t43 * t9 + t42) + g(3) * (-t43 * t60 + t44 * t8 + t58 * t7 + t47));
U  = t5;
