% Calculate potential energy for
% S6PRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-03-08 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRPP3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRRPP3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP3_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPP3_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:54:21
% EndTime: 2019-03-08 22:54:21
% DurationCPUTime: 0.45s
% Computational Cost: add. (356->116), mult. (739->138), div. (0->0), fcn. (903->10), ass. (0->53)
t72 = rSges(7,2) + qJ(5);
t71 = rSges(7,3) + qJ(6);
t37 = sin(pkin(6));
t70 = pkin(7) * t37;
t69 = rSges(6,1) + pkin(9);
t68 = rSges(4,3) + pkin(8);
t67 = rSges(5,3) + pkin(9);
t66 = cos(qJ(3));
t40 = sin(qJ(3));
t65 = t37 * t40;
t41 = sin(qJ(2));
t64 = t37 * t41;
t43 = cos(qJ(2));
t63 = t37 * t43;
t62 = rSges(6,3) + qJ(5);
t61 = cos(pkin(6));
t36 = sin(pkin(10));
t59 = pkin(1) * t36 + r_base(2);
t58 = qJ(1) + r_base(3);
t57 = t37 * t66;
t56 = t41 * t61;
t55 = t43 * t61;
t38 = cos(pkin(10));
t54 = pkin(1) * t38 + t36 * t70 + r_base(1);
t53 = pkin(7) * t61 + t58;
t24 = -t36 * t56 + t38 * t43;
t52 = pkin(2) * t24 + t54;
t51 = pkin(2) * t64 + t53;
t22 = t36 * t43 + t38 * t56;
t50 = pkin(2) * t22 - t38 * t70 + t59;
t13 = t24 * t66 + t36 * t65;
t23 = t36 * t55 + t38 * t41;
t49 = pkin(3) * t13 + t23 * pkin(8) + t52;
t39 = sin(qJ(4));
t42 = cos(qJ(4));
t6 = t13 * t42 + t23 * t39;
t48 = pkin(4) * t6 + t49;
t26 = t40 * t61 + t41 * t57;
t47 = pkin(3) * t26 - pkin(8) * t63 + t51;
t15 = t26 * t42 - t39 * t63;
t46 = pkin(4) * t15 + t47;
t11 = t22 * t66 - t38 * t65;
t21 = t36 * t41 - t38 * t55;
t45 = pkin(3) * t11 + t21 * pkin(8) + t50;
t4 = t11 * t42 + t21 * t39;
t44 = pkin(4) * t4 + t45;
t25 = t40 * t64 - t61 * t66;
t14 = t26 * t39 + t42 * t63;
t12 = t24 * t40 - t36 * t57;
t10 = t22 * t40 + t38 * t57;
t5 = t13 * t39 - t23 * t42;
t3 = t11 * t39 - t21 * t42;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t38 - rSges(2,2) * t36 + r_base(1)) + g(2) * (rSges(2,1) * t36 + rSges(2,2) * t38 + r_base(2)) + g(3) * (rSges(2,3) + t58)) - m(3) * (g(1) * (rSges(3,1) * t24 - rSges(3,2) * t23 + t54) + g(2) * (t22 * rSges(3,1) - t21 * rSges(3,2) + t59) + g(3) * (t61 * rSges(3,3) + t53) + (g(1) * rSges(3,3) * t36 + g(3) * (rSges(3,1) * t41 + rSges(3,2) * t43) + g(2) * (-rSges(3,3) - pkin(7)) * t38) * t37) - m(4) * (g(1) * (t13 * rSges(4,1) - t12 * rSges(4,2) + t23 * t68 + t52) + g(2) * (t11 * rSges(4,1) - t10 * rSges(4,2) + t21 * t68 + t50) + g(3) * (t26 * rSges(4,1) - t25 * rSges(4,2) - t63 * t68 + t51)) - m(5) * (g(1) * (t6 * rSges(5,1) - t5 * rSges(5,2) + t12 * t67 + t49) + g(2) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t10 * t67 + t45) + g(3) * (t15 * rSges(5,1) - t14 * rSges(5,2) + t25 * t67 + t47)) - m(6) * (g(1) * (-t6 * rSges(6,2) + t12 * t69 + t5 * t62 + t48) + g(2) * (-t4 * rSges(6,2) + t10 * t69 + t3 * t62 + t44) + g(3) * (-t15 * rSges(6,2) + t14 * t62 + t25 * t69 + t46)) - m(7) * (g(1) * (t72 * t5 + t71 * t6 + t48) + g(2) * (t72 * t3 + t71 * t4 + t44) + g(3) * (t72 * t14 + t71 * t15 + t46) + (g(1) * t12 + g(2) * t10 + g(3) * t25) * (rSges(7,1) + pkin(5) + pkin(9)));
U  = t1;
