% Calculate potential energy for
% S6PRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRPR8_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRRPR8_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR8_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_energypot_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR8_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR8_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:48:05
% EndTime: 2019-03-08 23:48:05
% DurationCPUTime: 0.50s
% Computational Cost: add. (553->126), mult. (1342->161), div. (0->0), fcn. (1691->14), ass. (0->63)
t41 = sin(pkin(7));
t44 = cos(pkin(7));
t45 = cos(pkin(6));
t42 = sin(pkin(6));
t51 = cos(qJ(2));
t79 = t42 * t51;
t63 = -t41 * t79 + t45 * t44;
t40 = sin(pkin(12));
t43 = cos(pkin(12));
t49 = sin(qJ(2));
t76 = t45 * t51;
t28 = -t40 * t76 - t43 * t49;
t81 = t42 * t44;
t64 = -t28 * t41 + t40 * t81;
t89 = rSges(6,1) + pkin(10);
t88 = rSges(5,3) + pkin(10);
t87 = pkin(11) + rSges(7,3);
t86 = cos(qJ(3));
t85 = cos(qJ(4));
t83 = t40 * t42;
t82 = t42 * t43;
t80 = t42 * t49;
t77 = t45 * t49;
t75 = rSges(6,3) + qJ(5);
t74 = t40 * pkin(1) + r_base(2);
t71 = qJ(1) + r_base(3);
t70 = t41 * t86;
t69 = t44 * t86;
t68 = t43 * pkin(1) + pkin(8) * t83 + r_base(1);
t67 = t45 * pkin(8) + t71;
t66 = t42 * t70;
t26 = -t40 * t49 + t43 * t76;
t65 = -t26 * t41 - t43 * t81;
t46 = sin(qJ(6));
t50 = cos(qJ(6));
t62 = rSges(7,1) * t46 + rSges(7,2) * t50 + qJ(5);
t61 = rSges(7,1) * t50 - rSges(7,2) * t46 + pkin(5) + pkin(10);
t29 = -t40 * t77 + t43 * t51;
t60 = t29 * pkin(2) + t64 * pkin(9) + t68;
t48 = sin(qJ(3));
t13 = t29 * t86 + (t28 * t44 + t41 * t83) * t48;
t59 = t13 * pkin(3) + t60;
t58 = pkin(2) * t80 + t63 * pkin(9) + t67;
t47 = sin(qJ(4));
t6 = t13 * t85 + t64 * t47;
t57 = t6 * pkin(4) + t59;
t20 = t45 * t41 * t48 + (t44 * t48 * t51 + t86 * t49) * t42;
t56 = t20 * pkin(3) + t58;
t15 = t20 * t85 + t63 * t47;
t55 = t15 * pkin(4) + t56;
t27 = t40 * t51 + t43 * t77;
t54 = t27 * pkin(2) - pkin(8) * t82 + t65 * pkin(9) + t74;
t11 = t27 * t86 + (t26 * t44 - t41 * t82) * t48;
t53 = t11 * pkin(3) + t54;
t4 = t11 * t85 + t65 * t47;
t52 = t4 * pkin(4) + t53;
t19 = -t45 * t70 + t48 * t80 - t69 * t79;
t14 = t20 * t47 - t63 * t85;
t12 = -t28 * t69 + t29 * t48 - t40 * t66;
t10 = -t26 * t69 + t27 * t48 + t43 * t66;
t5 = t13 * t47 - t64 * t85;
t3 = t11 * t47 - t65 * t85;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t43 - rSges(2,2) * t40 + r_base(1)) + g(2) * (rSges(2,1) * t40 + rSges(2,2) * t43 + r_base(2)) + g(3) * (rSges(2,3) + t71)) - m(3) * (g(1) * (rSges(3,1) * t29 + rSges(3,2) * t28 + t68) + g(2) * (t27 * rSges(3,1) + t26 * rSges(3,2) + t74) + g(3) * (t45 * rSges(3,3) + t67) + (g(1) * rSges(3,3) * t40 + g(3) * (rSges(3,1) * t49 + rSges(3,2) * t51) + g(2) * (-rSges(3,3) - pkin(8)) * t43) * t42) - m(4) * (g(1) * (t13 * rSges(4,1) - t12 * rSges(4,2) + t64 * rSges(4,3) + t60) + g(2) * (t11 * rSges(4,1) - t10 * rSges(4,2) + t65 * rSges(4,3) + t54) + g(3) * (t20 * rSges(4,1) - t19 * rSges(4,2) + t63 * rSges(4,3) + t58)) - m(5) * (g(1) * (t6 * rSges(5,1) - t5 * rSges(5,2) + t88 * t12 + t59) + g(2) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t88 * t10 + t53) + g(3) * (t15 * rSges(5,1) - t14 * rSges(5,2) + t88 * t19 + t56)) - m(6) * (g(1) * (-t6 * rSges(6,2) + t89 * t12 + t75 * t5 + t57) + g(2) * (-t4 * rSges(6,2) + t89 * t10 + t75 * t3 + t52) + g(3) * (-t15 * rSges(6,2) + t75 * t14 + t89 * t19 + t55)) - m(7) * (g(1) * (t61 * t12 + t62 * t5 + t87 * t6 + t57) + g(2) * (t61 * t10 + t62 * t3 + t87 * t4 + t52) + g(3) * (t62 * t14 + t87 * t15 + t61 * t19 + t55));
U  = t1;
