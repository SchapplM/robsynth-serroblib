% Calculate potential energy for
% S6PRRRPP2
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
% Datum: 2019-03-08 22:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRPP2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRRPP2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP2_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPP2_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:48:38
% EndTime: 2019-03-08 22:48:39
% DurationCPUTime: 0.46s
% Computational Cost: add. (356->116), mult. (739->138), div. (0->0), fcn. (903->10), ass. (0->53)
t71 = rSges(7,1) + pkin(5);
t70 = rSges(7,2) + qJ(5);
t36 = sin(pkin(6));
t69 = pkin(7) * t36;
t68 = rSges(6,2) + pkin(9);
t67 = rSges(4,3) + pkin(8);
t66 = rSges(5,3) + pkin(9);
t65 = cos(qJ(3));
t39 = sin(qJ(3));
t64 = t36 * t39;
t40 = sin(qJ(2));
t63 = t36 * t40;
t42 = cos(qJ(2));
t62 = t36 * t42;
t61 = rSges(6,3) + qJ(5);
t60 = cos(pkin(6));
t35 = sin(pkin(10));
t59 = t35 * pkin(1) + r_base(2);
t58 = qJ(1) + r_base(3);
t56 = t36 * t65;
t55 = t40 * t60;
t54 = t42 * t60;
t37 = cos(pkin(10));
t53 = t37 * pkin(1) + t35 * t69 + r_base(1);
t52 = t60 * pkin(7) + t58;
t24 = -t35 * t55 + t37 * t42;
t51 = t24 * pkin(2) + t53;
t50 = pkin(2) * t63 + t52;
t22 = t35 * t42 + t37 * t55;
t49 = t22 * pkin(2) - t37 * t69 + t59;
t13 = t24 * t65 + t35 * t64;
t23 = t35 * t54 + t37 * t40;
t48 = t13 * pkin(3) + t23 * pkin(8) + t51;
t38 = sin(qJ(4));
t41 = cos(qJ(4));
t6 = t13 * t41 + t23 * t38;
t47 = t6 * pkin(4) + t48;
t26 = t39 * t60 + t40 * t56;
t46 = t26 * pkin(3) - pkin(8) * t62 + t50;
t15 = t26 * t41 - t38 * t62;
t45 = t15 * pkin(4) + t46;
t11 = t22 * t65 - t37 * t64;
t21 = t35 * t40 - t37 * t54;
t44 = t11 * pkin(3) + t21 * pkin(8) + t49;
t4 = t11 * t41 + t21 * t38;
t43 = t4 * pkin(4) + t44;
t25 = t39 * t63 - t60 * t65;
t14 = t26 * t38 + t41 * t62;
t12 = t24 * t39 - t35 * t56;
t10 = t22 * t39 + t37 * t56;
t5 = t13 * t38 - t23 * t41;
t3 = t11 * t38 - t21 * t41;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t37 * rSges(2,1) - t35 * rSges(2,2) + r_base(1)) + g(2) * (t35 * rSges(2,1) + t37 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t58)) - m(3) * (g(1) * (t24 * rSges(3,1) - t23 * rSges(3,2) + t53) + g(2) * (t22 * rSges(3,1) - t21 * rSges(3,2) + t59) + g(3) * (t60 * rSges(3,3) + t52) + (g(1) * rSges(3,3) * t35 + g(3) * (rSges(3,1) * t40 + rSges(3,2) * t42) + g(2) * (-rSges(3,3) - pkin(7)) * t37) * t36) - m(4) * (g(1) * (t13 * rSges(4,1) - t12 * rSges(4,2) + t67 * t23 + t51) + g(2) * (t11 * rSges(4,1) - t10 * rSges(4,2) + t67 * t21 + t49) + g(3) * (t26 * rSges(4,1) - t25 * rSges(4,2) - t67 * t62 + t50)) - m(5) * (g(1) * (t6 * rSges(5,1) - t5 * rSges(5,2) + t66 * t12 + t48) + g(2) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t66 * t10 + t44) + g(3) * (t15 * rSges(5,1) - t14 * rSges(5,2) + t66 * t25 + t46)) - m(6) * (g(1) * (t6 * rSges(6,1) + t68 * t12 + t61 * t5 + t47) + g(2) * (t4 * rSges(6,1) + t68 * t10 + t61 * t3 + t43) + g(3) * (t15 * rSges(6,1) + t61 * t14 + t68 * t25 + t45)) - m(7) * (g(1) * (t70 * t5 + t71 * t6 + t47) + g(2) * (t70 * t3 + t71 * t4 + t43) + g(3) * (t70 * t14 + t71 * t15 + t45) + (g(1) * t12 + g(2) * t10 + g(3) * t25) * (-rSges(7,3) + pkin(9) - qJ(6)));
U  = t1;
