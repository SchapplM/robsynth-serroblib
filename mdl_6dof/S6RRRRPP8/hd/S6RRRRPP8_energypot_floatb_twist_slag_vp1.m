% Calculate potential energy for
% S6RRRRPP8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPP8_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRPP8_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP8_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP8_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP8_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:30:12
% EndTime: 2019-03-09 21:30:12
% DurationCPUTime: 0.46s
% Computational Cost: add. (356->116), mult. (739->138), div. (0->0), fcn. (903->10), ass. (0->53)
t71 = rSges(7,1) + pkin(5);
t70 = rSges(7,2) + qJ(5);
t69 = rSges(6,2) + pkin(10);
t68 = rSges(4,3) + pkin(9);
t67 = rSges(5,3) + pkin(10);
t66 = cos(qJ(3));
t35 = sin(pkin(6));
t38 = sin(qJ(2));
t65 = t35 * t38;
t39 = sin(qJ(1));
t64 = t35 * t39;
t41 = cos(qJ(2));
t63 = t35 * t41;
t42 = cos(qJ(1));
t62 = t35 * t42;
t61 = rSges(6,3) + qJ(5);
t60 = cos(pkin(6));
t59 = pkin(7) + r_base(3);
t58 = t39 * pkin(1) + r_base(2);
t56 = t35 * t66;
t55 = t60 * pkin(8) + t59;
t54 = t39 * t60;
t53 = t42 * t60;
t52 = t42 * pkin(1) + pkin(8) * t64 + r_base(1);
t51 = pkin(2) * t65 + t55;
t26 = -t38 * t54 + t42 * t41;
t50 = t26 * pkin(2) + t52;
t24 = t38 * t53 + t39 * t41;
t49 = t24 * pkin(2) - pkin(8) * t62 + t58;
t37 = sin(qJ(3));
t15 = t26 * t66 + t37 * t64;
t25 = t42 * t38 + t41 * t54;
t48 = t15 * pkin(3) + t25 * pkin(9) + t50;
t36 = sin(qJ(4));
t40 = cos(qJ(4));
t6 = t15 * t40 + t25 * t36;
t47 = t6 * pkin(4) + t48;
t22 = t37 * t60 + t38 * t56;
t46 = t22 * pkin(3) - pkin(9) * t63 + t51;
t11 = t22 * t40 - t36 * t63;
t45 = t11 * pkin(4) + t46;
t13 = t24 * t66 - t37 * t62;
t23 = t39 * t38 - t41 * t53;
t44 = t13 * pkin(3) + t23 * pkin(9) + t49;
t4 = t13 * t40 + t23 * t36;
t43 = t4 * pkin(4) + t44;
t21 = t37 * t65 - t60 * t66;
t14 = t26 * t37 - t39 * t56;
t12 = t24 * t37 + t42 * t56;
t10 = t22 * t36 + t40 * t63;
t5 = t15 * t36 - t25 * t40;
t3 = t13 * t36 - t23 * t40;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t42 * rSges(2,1) - t39 * rSges(2,2) + r_base(1)) + g(2) * (t39 * rSges(2,1) + t42 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t59)) - m(3) * (g(1) * (t26 * rSges(3,1) - t25 * rSges(3,2) + t52) + g(2) * (t24 * rSges(3,1) - t23 * rSges(3,2) + t58) + g(3) * (t60 * rSges(3,3) + t55) + (g(1) * rSges(3,3) * t39 + g(3) * (rSges(3,1) * t38 + rSges(3,2) * t41) + g(2) * (-rSges(3,3) - pkin(8)) * t42) * t35) - m(4) * (g(1) * (t15 * rSges(4,1) - t14 * rSges(4,2) + t25 * t68 + t50) + g(2) * (t13 * rSges(4,1) - t12 * rSges(4,2) + t23 * t68 + t49) + g(3) * (t22 * rSges(4,1) - t21 * rSges(4,2) - t63 * t68 + t51)) - m(5) * (g(1) * (t6 * rSges(5,1) - t5 * rSges(5,2) + t14 * t67 + t48) + g(2) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t12 * t67 + t44) + g(3) * (t11 * rSges(5,1) - t10 * rSges(5,2) + t21 * t67 + t46)) - m(6) * (g(1) * (t6 * rSges(6,1) + t14 * t69 + t61 * t5 + t47) + g(2) * (t4 * rSges(6,1) + t12 * t69 + t61 * t3 + t43) + g(3) * (t11 * rSges(6,1) + t61 * t10 + t21 * t69 + t45)) - m(7) * (g(1) * (t70 * t5 + t71 * t6 + t47) + g(2) * (t70 * t3 + t71 * t4 + t43) + g(3) * (t70 * t10 + t71 * t11 + t45) + (g(1) * t14 + g(2) * t12 + g(3) * t21) * (-rSges(7,3) + pkin(10) - qJ(6)));
U  = t1;
