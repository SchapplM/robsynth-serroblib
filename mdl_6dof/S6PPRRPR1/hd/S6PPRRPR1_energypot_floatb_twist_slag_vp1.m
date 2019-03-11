% Calculate potential energy for
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
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
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PPRRPR1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PPRRPR1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_energypot_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR1_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRPR1_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:45:11
% EndTime: 2019-03-08 18:45:12
% DurationCPUTime: 0.53s
% Computational Cost: add. (584->136), mult. (1388->178), div. (0->0), fcn. (1753->16), ass. (0->62)
t44 = sin(pkin(7));
t49 = cos(pkin(7));
t50 = cos(pkin(6));
t45 = sin(pkin(6));
t47 = cos(pkin(12));
t80 = t45 * t47;
t62 = -t44 * t80 + t49 * t50;
t42 = sin(pkin(12));
t48 = cos(pkin(11));
t43 = sin(pkin(11));
t82 = t43 * t50;
t25 = -t42 * t48 - t47 * t82;
t79 = t45 * t49;
t63 = -t25 * t44 + t43 * t79;
t87 = rSges(5,3) + pkin(9);
t86 = cos(qJ(3));
t85 = cos(qJ(4));
t83 = t42 * t45;
t81 = t44 * t45;
t78 = t48 * t50;
t76 = pkin(10) + qJ(5) + rSges(7,3);
t75 = t45 * qJ(2);
t74 = qJ(5) + rSges(6,3);
t73 = t43 * pkin(1) + r_base(2);
t70 = qJ(1) + r_base(3);
t69 = t44 * t86;
t68 = t49 * t86;
t67 = t48 * pkin(1) + t43 * t75 + r_base(1);
t66 = t50 * qJ(2) + t70;
t65 = t45 * t69;
t23 = -t42 * t43 + t47 * t78;
t64 = -t23 * t44 - t48 * t79;
t40 = pkin(13) + qJ(6);
t35 = sin(t40);
t36 = cos(t40);
t46 = cos(pkin(13));
t61 = t36 * rSges(7,1) - t35 * rSges(7,2) + pkin(5) * t46 + pkin(4);
t41 = sin(pkin(13));
t60 = t35 * rSges(7,1) + t36 * rSges(7,2) + t41 * pkin(5) + pkin(9);
t26 = -t42 * t82 + t47 * t48;
t59 = t26 * pkin(2) + t63 * pkin(8) + t67;
t53 = sin(qJ(3));
t10 = t26 * t86 + (t25 * t49 + t43 * t81) * t53;
t58 = t10 * pkin(3) + t59;
t57 = pkin(2) * t83 + t62 * pkin(8) + t66;
t17 = t50 * t44 * t53 + (t47 * t49 * t53 + t86 * t42) * t45;
t56 = t17 * pkin(3) + t57;
t24 = t42 * t78 + t43 * t47;
t55 = t24 * pkin(2) + t64 * pkin(8) - t48 * t75 + t73;
t8 = t24 * t86 + (t23 * t49 - t48 * t81) * t53;
t54 = t8 * pkin(3) + t55;
t52 = sin(qJ(4));
t16 = -t50 * t69 + t53 * t83 - t68 * t80;
t12 = t17 * t85 + t62 * t52;
t11 = t17 * t52 - t62 * t85;
t9 = -t25 * t68 + t26 * t53 - t43 * t65;
t7 = -t23 * t68 + t24 * t53 + t48 * t65;
t4 = t10 * t85 + t63 * t52;
t3 = t10 * t52 - t63 * t85;
t2 = t64 * t52 + t8 * t85;
t1 = t52 * t8 - t64 * t85;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t48 - rSges(2,2) * t43 + r_base(1)) + g(2) * (rSges(2,1) * t43 + rSges(2,2) * t48 + r_base(2)) + g(3) * (rSges(2,3) + t70)) - m(3) * (g(1) * (rSges(3,1) * t26 + rSges(3,2) * t25 + t67) + g(2) * (rSges(3,1) * t24 + rSges(3,2) * t23 + t73) + g(3) * (rSges(3,3) * t50 + t66) + (g(1) * rSges(3,3) * t43 + g(3) * (rSges(3,1) * t42 + rSges(3,2) * t47) + g(2) * (-rSges(3,3) - qJ(2)) * t48) * t45) - m(4) * (g(1) * (t10 * rSges(4,1) - t9 * rSges(4,2) + t63 * rSges(4,3) + t59) + g(2) * (t8 * rSges(4,1) - t7 * rSges(4,2) + t64 * rSges(4,3) + t55) + g(3) * (t17 * rSges(4,1) - t16 * rSges(4,2) + t62 * rSges(4,3) + t57)) - m(5) * (g(1) * (rSges(5,1) * t4 - rSges(5,2) * t3 + t87 * t9 + t58) + g(2) * (rSges(5,1) * t2 - rSges(5,2) * t1 + t87 * t7 + t54) + g(3) * (rSges(5,1) * t12 - rSges(5,2) * t11 + t87 * t16 + t56)) - m(6) * (g(1) * (t4 * pkin(4) + t9 * pkin(9) + (t4 * t46 + t41 * t9) * rSges(6,1) + (-t4 * t41 + t46 * t9) * rSges(6,2) + t74 * t3 + t58) + g(2) * (t2 * pkin(4) + t7 * pkin(9) + (t2 * t46 + t41 * t7) * rSges(6,1) + (-t2 * t41 + t46 * t7) * rSges(6,2) + t74 * t1 + t54) + g(3) * (t12 * pkin(4) + t16 * pkin(9) + (t12 * t46 + t16 * t41) * rSges(6,1) + (-t12 * t41 + t16 * t46) * rSges(6,2) + t74 * t11 + t56)) - m(7) * (g(1) * (t76 * t3 + t61 * t4 + t60 * t9 + t58) + g(2) * (t76 * t1 + t61 * t2 + t60 * t7 + t54) + g(3) * (t76 * t11 + t61 * t12 + t60 * t16 + t56));
U  = t5;
