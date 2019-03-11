% Calculate potential energy for
% S6PRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
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
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPPRR1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPPRR1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_energypot_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR1_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPPRR1_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:13:49
% EndTime: 2019-03-08 19:13:49
% DurationCPUTime: 0.65s
% Computational Cost: add. (404->136), mult. (724->174), div. (0->0), fcn. (880->14), ass. (0->54)
t41 = sin(pkin(11));
t45 = cos(pkin(11));
t50 = sin(qJ(2));
t52 = cos(qJ(2));
t24 = -t41 * t50 + t45 * t52;
t74 = pkin(9) + rSges(7,3);
t42 = sin(pkin(10));
t43 = sin(pkin(6));
t73 = t42 * t43;
t46 = cos(pkin(10));
t72 = t46 * t43;
t40 = sin(pkin(12));
t47 = cos(pkin(6));
t71 = t47 * t40;
t70 = t47 * t50;
t69 = t47 * t52;
t66 = qJ(4) + rSges(5,3);
t34 = pkin(2) * t52 + pkin(1);
t65 = t34 * t46 + r_base(1);
t64 = t40 * t73;
t63 = t40 * t72;
t62 = qJ(1) + r_base(3);
t22 = pkin(2) * t70 + (-pkin(7) - qJ(3)) * t43;
t61 = t22 * t46 + t34 * t42 + r_base(2);
t60 = pkin(7) * t47 + t62;
t59 = t41 * t52 + t45 * t50;
t58 = -t42 * t22 + t65;
t57 = pkin(2) * t43 * t50 + qJ(3) * t47 + t60;
t56 = t24 * t47;
t11 = -t42 * t56 - t46 * t59;
t21 = t59 * t47;
t12 = -t21 * t42 + t24 * t46;
t44 = cos(pkin(12));
t33 = pkin(4) * t44 + pkin(3);
t48 = -pkin(8) - qJ(4);
t55 = pkin(4) * t64 + t11 * t48 + t12 * t33 + t58;
t19 = t24 * t43;
t20 = t59 * t43;
t54 = pkin(4) * t71 + t19 * t48 + t20 * t33 + t57;
t10 = t21 * t46 + t24 * t42;
t9 = -t42 * t59 + t46 * t56;
t53 = -pkin(4) * t63 + t10 * t33 + t48 * t9 + t61;
t51 = cos(qJ(6));
t49 = sin(qJ(6));
t39 = pkin(12) + qJ(5);
t36 = cos(t39);
t35 = sin(t39);
t14 = t20 * t36 + t35 * t47;
t13 = t20 * t35 - t36 * t47;
t4 = t12 * t36 + t35 * t73;
t3 = t12 * t35 - t36 * t73;
t2 = t10 * t36 - t35 * t72;
t1 = t10 * t35 + t36 * t72;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t46 - rSges(2,2) * t42 + r_base(1)) + g(2) * (rSges(2,1) * t42 + rSges(2,2) * t46 + r_base(2)) + g(3) * (rSges(2,3) + t62)) - m(3) * (g(1) * (t46 * pkin(1) + r_base(1) + (-t42 * t70 + t46 * t52) * rSges(3,1) + (-t42 * t69 - t46 * t50) * rSges(3,2)) + g(2) * (t42 * pkin(1) + r_base(2) + (t42 * t52 + t46 * t70) * rSges(3,1) + (-t42 * t50 + t46 * t69) * rSges(3,2)) + g(3) * (t47 * rSges(3,3) + t60) + (g(3) * (rSges(3,1) * t50 + rSges(3,2) * t52) + (g(1) * t42 - g(2) * t46) * (rSges(3,3) + pkin(7))) * t43) - m(4) * (g(1) * (t12 * rSges(4,1) + t11 * rSges(4,2) + (rSges(4,3) * t43 - t22) * t42 + t65) + g(2) * (rSges(4,1) * t10 + rSges(4,2) * t9 - rSges(4,3) * t72 + t61) + g(3) * (rSges(4,1) * t20 + rSges(4,2) * t19 + rSges(4,3) * t47 + t57)) - m(5) * (g(1) * (t12 * pkin(3) + (t12 * t44 + t64) * rSges(5,1) + (-t12 * t40 + t44 * t73) * rSges(5,2) - t66 * t11 + t58) + g(2) * (t10 * pkin(3) + (t10 * t44 - t63) * rSges(5,1) + (-t10 * t40 - t44 * t72) * rSges(5,2) - t66 * t9 + t61) + g(3) * (t20 * pkin(3) + (t20 * t44 + t71) * rSges(5,1) + (-t20 * t40 + t44 * t47) * rSges(5,2) - t66 * t19 + t57)) - m(6) * (g(1) * (rSges(6,1) * t4 - rSges(6,2) * t3 - rSges(6,3) * t11 + t55) + g(2) * (rSges(6,1) * t2 - rSges(6,2) * t1 - rSges(6,3) * t9 + t53) + g(3) * (rSges(6,1) * t14 - rSges(6,2) * t13 - rSges(6,3) * t19 + t54)) - m(7) * (g(1) * (t4 * pkin(5) + (-t11 * t49 + t4 * t51) * rSges(7,1) + (-t11 * t51 - t4 * t49) * rSges(7,2) + t74 * t3 + t55) + g(2) * (t2 * pkin(5) + (t2 * t51 - t49 * t9) * rSges(7,1) + (-t2 * t49 - t51 * t9) * rSges(7,2) + t74 * t1 + t53) + g(3) * (t14 * pkin(5) + (t14 * t51 - t19 * t49) * rSges(7,1) + (-t14 * t49 - t19 * t51) * rSges(7,2) + t74 * t13 + t54));
U  = t5;
