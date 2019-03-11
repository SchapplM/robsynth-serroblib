% Calculate potential energy for
% S6PRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRPR6_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRRPR6_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR6_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR6_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR6_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:31:02
% EndTime: 2019-03-08 23:31:03
% DurationCPUTime: 0.50s
% Computational Cost: add. (382->126), mult. (807->156), div. (0->0), fcn. (997->12), ass. (0->53)
t36 = sin(pkin(6));
t71 = pkin(7) * t36;
t70 = rSges(6,2) + pkin(9);
t69 = rSges(4,3) + pkin(8);
t68 = rSges(5,3) + pkin(9);
t67 = cos(qJ(3));
t40 = sin(qJ(3));
t66 = t36 * t40;
t41 = sin(qJ(2));
t65 = t36 * t41;
t44 = cos(qJ(2));
t64 = t36 * t44;
t63 = rSges(6,3) + qJ(5);
t62 = cos(pkin(6));
t35 = sin(pkin(11));
t60 = t35 * pkin(1) + r_base(2);
t59 = qJ(1) + r_base(3);
t58 = t36 * t67;
t57 = t41 * t62;
t56 = t44 * t62;
t37 = cos(pkin(11));
t55 = t37 * pkin(1) + t35 * t71 + r_base(1);
t54 = t62 * pkin(7) + t59;
t24 = -t35 * t57 + t37 * t44;
t53 = t24 * pkin(2) + t55;
t52 = pkin(2) * t65 + t54;
t22 = t35 * t44 + t37 * t57;
t51 = t22 * pkin(2) - t37 * t71 + t60;
t13 = t24 * t67 + t35 * t66;
t23 = t35 * t56 + t37 * t41;
t50 = t13 * pkin(3) + t23 * pkin(8) + t53;
t39 = sin(qJ(4));
t43 = cos(qJ(4));
t6 = t13 * t43 + t23 * t39;
t49 = t6 * pkin(4) + t50;
t26 = t62 * t40 + t41 * t58;
t48 = t26 * pkin(3) - pkin(8) * t64 + t52;
t15 = t26 * t43 - t39 * t64;
t47 = t15 * pkin(4) + t48;
t11 = t22 * t67 - t37 * t66;
t21 = t35 * t41 - t37 * t56;
t46 = t11 * pkin(3) + t21 * pkin(8) + t51;
t4 = t11 * t43 + t21 * t39;
t45 = t4 * pkin(4) + t46;
t42 = cos(qJ(6));
t38 = sin(qJ(6));
t25 = t40 * t65 - t62 * t67;
t14 = t26 * t39 + t43 * t64;
t12 = t24 * t40 - t35 * t58;
t10 = t22 * t40 + t37 * t58;
t5 = t13 * t39 - t23 * t43;
t3 = t11 * t39 - t21 * t43;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t37 * rSges(2,1) - t35 * rSges(2,2) + r_base(1)) + g(2) * (t35 * rSges(2,1) + t37 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t59)) - m(3) * (g(1) * (t24 * rSges(3,1) - t23 * rSges(3,2) + t55) + g(2) * (t22 * rSges(3,1) - t21 * rSges(3,2) + t60) + g(3) * (t62 * rSges(3,3) + t54) + (g(1) * rSges(3,3) * t35 + g(3) * (rSges(3,1) * t41 + rSges(3,2) * t44) + g(2) * (-rSges(3,3) - pkin(7)) * t37) * t36) - m(4) * (g(1) * (t13 * rSges(4,1) - t12 * rSges(4,2) + t69 * t23 + t53) + g(2) * (t11 * rSges(4,1) - t10 * rSges(4,2) + t69 * t21 + t51) + g(3) * (t26 * rSges(4,1) - t25 * rSges(4,2) - t69 * t64 + t52)) - m(5) * (g(1) * (t6 * rSges(5,1) - t5 * rSges(5,2) + t68 * t12 + t50) + g(2) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t68 * t10 + t46) + g(3) * (t15 * rSges(5,1) - t14 * rSges(5,2) + t68 * t25 + t48)) - m(6) * (g(1) * (t6 * rSges(6,1) + t70 * t12 + t63 * t5 + t49) + g(2) * (t4 * rSges(6,1) + t70 * t10 + t63 * t3 + t45) + g(3) * (t15 * rSges(6,1) + t63 * t14 + t70 * t25 + t47)) - m(7) * (g(1) * (t6 * pkin(5) + t5 * qJ(5) + (t5 * t38 + t6 * t42) * rSges(7,1) + (-t6 * t38 + t5 * t42) * rSges(7,2) + t49) + g(2) * (t4 * pkin(5) + t3 * qJ(5) + (t3 * t38 + t4 * t42) * rSges(7,1) + (t3 * t42 - t4 * t38) * rSges(7,2) + t45) + g(3) * (t15 * pkin(5) + t14 * qJ(5) + (t14 * t38 + t15 * t42) * rSges(7,1) + (t14 * t42 - t15 * t38) * rSges(7,2) + t47) + (g(1) * t12 + g(2) * t10 + g(3) * t25) * (-pkin(10) + pkin(9) - rSges(7,3)));
U  = t1;
