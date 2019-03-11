% Calculate potential energy for
% S6PRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRR4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRPRR4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR4_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR4_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:10:09
% EndTime: 2019-03-08 22:10:10
% DurationCPUTime: 0.44s
% Computational Cost: add. (355->124), mult. (737->148), div. (0->0), fcn. (900->12), ass. (0->57)
t38 = sin(pkin(6));
t76 = pkin(7) * t38;
t75 = rSges(5,2) + pkin(8);
t74 = rSges(4,3) + pkin(8);
t73 = rSges(6,3) - pkin(8);
t72 = pkin(10) + rSges(7,3);
t71 = cos(qJ(3));
t42 = sin(qJ(3));
t70 = t38 * t42;
t43 = sin(qJ(2));
t69 = t38 * t43;
t46 = cos(qJ(2));
t68 = t38 * t46;
t67 = rSges(5,3) + qJ(4);
t66 = cos(pkin(6));
t65 = -pkin(9) - t73;
t37 = sin(pkin(11));
t64 = t37 * pkin(1) + r_base(2);
t63 = qJ(1) + r_base(3);
t62 = t38 * t71;
t61 = t43 * t66;
t60 = t46 * t66;
t39 = cos(pkin(11));
t59 = t39 * pkin(1) + t37 * t76 + r_base(1);
t58 = t66 * pkin(7) + t63;
t25 = -t37 * t61 + t39 * t46;
t57 = t25 * pkin(2) + t59;
t56 = pkin(2) * t69 + t58;
t16 = t25 * t71 + t37 * t70;
t55 = t16 * pkin(3) + t57;
t27 = t66 * t42 + t43 * t62;
t54 = t27 * pkin(3) + t56;
t40 = sin(qJ(6));
t44 = cos(qJ(6));
t53 = t44 * rSges(7,1) - t40 * rSges(7,2) + pkin(5);
t52 = -t40 * rSges(7,1) - t44 * rSges(7,2) + pkin(8) - pkin(9);
t23 = t37 * t46 + t39 * t61;
t51 = t23 * pkin(2) - t39 * t76 + t64;
t14 = t23 * t71 - t39 * t70;
t50 = t14 * pkin(3) + t51;
t15 = t25 * t42 - t37 * t62;
t49 = t16 * pkin(4) + t15 * qJ(4) + t55;
t26 = t42 * t69 - t66 * t71;
t48 = t27 * pkin(4) + pkin(9) * t68 + t26 * qJ(4) + t54;
t13 = t23 * t42 + t39 * t62;
t47 = t14 * pkin(4) + t13 * qJ(4) + t50;
t45 = cos(qJ(5));
t41 = sin(qJ(5));
t24 = t37 * t60 + t39 * t43;
t22 = t37 * t43 - t39 * t60;
t6 = t26 * t41 + t27 * t45;
t5 = -t26 * t45 + t27 * t41;
t4 = t15 * t41 + t16 * t45;
t3 = -t15 * t45 + t16 * t41;
t2 = t13 * t41 + t14 * t45;
t1 = -t13 * t45 + t14 * t41;
t7 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t39 * rSges(2,1) - t37 * rSges(2,2) + r_base(1)) + g(2) * (t37 * rSges(2,1) + t39 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t63)) - m(3) * (g(1) * (t25 * rSges(3,1) - t24 * rSges(3,2) + t59) + g(2) * (t23 * rSges(3,1) - t22 * rSges(3,2) + t64) + g(3) * (t66 * rSges(3,3) + t58) + (g(1) * rSges(3,3) * t37 + g(3) * (rSges(3,1) * t43 + rSges(3,2) * t46) + g(2) * (-rSges(3,3) - pkin(7)) * t39) * t38) - m(4) * (g(1) * (t16 * rSges(4,1) - t15 * rSges(4,2) + t74 * t24 + t57) + g(2) * (t14 * rSges(4,1) - t13 * rSges(4,2) + t74 * t22 + t51) + g(3) * (t27 * rSges(4,1) - t26 * rSges(4,2) - t74 * t68 + t56)) - m(5) * (g(1) * (t16 * rSges(5,1) + t67 * t15 + t75 * t24 + t55) + g(2) * (t14 * rSges(5,1) + t67 * t13 + t75 * t22 + t50) + g(3) * (t27 * rSges(5,1) + t67 * t26 - t75 * t68 + t54)) - m(6) * (g(3) * (t6 * rSges(6,1) - t5 * rSges(6,2) + t73 * t68 + t48) + (t2 * rSges(6,1) - t1 * rSges(6,2) + t65 * t22 + t47) * g(2) + (t4 * rSges(6,1) - t3 * rSges(6,2) + t65 * t24 + t49) * g(1)) - m(7) * (g(1) * (t52 * t24 + t72 * t3 + t53 * t4 + t49) + g(2) * (t72 * t1 + t53 * t2 + t52 * t22 + t47) + g(3) * (t6 * pkin(5) - pkin(8) * t68 + (t40 * t68 + t6 * t44) * rSges(7,1) + (-t6 * t40 + t44 * t68) * rSges(7,2) + t72 * t5 + t48));
U  = t7;
