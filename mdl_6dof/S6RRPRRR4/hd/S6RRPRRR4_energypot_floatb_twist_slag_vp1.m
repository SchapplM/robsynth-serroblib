% Calculate potential energy for
% S6RRPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRR4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_energypot_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR4_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR4_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:27:39
% EndTime: 2019-03-09 13:27:39
% DurationCPUTime: 0.67s
% Computational Cost: add. (404->136), mult. (724->172), div. (0->0), fcn. (880->14), ass. (0->57)
t40 = sin(pkin(12));
t42 = cos(pkin(12));
t46 = sin(qJ(2));
t50 = cos(qJ(2));
t24 = -t46 * t40 + t50 * t42;
t78 = pkin(2) * t46;
t76 = pkin(9) + rSges(5,3);
t75 = pkin(11) + rSges(7,3);
t43 = cos(pkin(6));
t45 = sin(qJ(4));
t74 = t43 * t45;
t41 = sin(pkin(6));
t47 = sin(qJ(1));
t72 = t47 * t41;
t71 = t47 * t46;
t70 = t47 * t50;
t51 = cos(qJ(1));
t68 = t51 * t41;
t67 = t51 * t46;
t66 = t51 * t50;
t65 = pkin(7) + r_base(3);
t34 = t50 * pkin(2) + pkin(1);
t64 = t51 * t34 + r_base(1);
t63 = t45 * t72;
t62 = t45 * t68;
t61 = t43 * pkin(8) + t65;
t22 = t43 * t78 + (-pkin(8) - qJ(3)) * t41;
t60 = t51 * t22 + t47 * t34 + r_base(2);
t59 = t50 * t40 + t46 * t42;
t58 = -t47 * t22 + t64;
t57 = t43 * qJ(3) + t41 * t78 + t61;
t56 = t24 * t43;
t11 = -t47 * t56 - t51 * t59;
t21 = t59 * t43;
t12 = -t47 * t21 + t51 * t24;
t49 = cos(qJ(4));
t33 = t49 * pkin(4) + pkin(3);
t52 = -pkin(10) - pkin(9);
t55 = pkin(4) * t63 + t11 * t52 + t12 * t33 + t58;
t19 = t24 * t41;
t20 = t59 * t41;
t54 = pkin(4) * t74 + t19 * t52 + t20 * t33 + t57;
t10 = t51 * t21 + t47 * t24;
t9 = -t47 * t59 + t51 * t56;
t53 = -pkin(4) * t62 + t10 * t33 + t9 * t52 + t60;
t48 = cos(qJ(6));
t44 = sin(qJ(6));
t39 = qJ(4) + qJ(5);
t37 = cos(t39);
t36 = sin(t39);
t14 = t20 * t37 + t43 * t36;
t13 = t20 * t36 - t43 * t37;
t4 = t12 * t37 + t36 * t72;
t3 = t12 * t36 - t37 * t72;
t2 = t10 * t37 - t36 * t68;
t1 = t10 * t36 + t37 * t68;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t51 * rSges(2,1) - t47 * rSges(2,2) + r_base(1)) + g(2) * (t47 * rSges(2,1) + t51 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t65)) - m(3) * (g(1) * (t51 * pkin(1) + r_base(1) + (-t43 * t71 + t66) * rSges(3,1) + (-t43 * t70 - t67) * rSges(3,2)) + g(2) * (t47 * pkin(1) + r_base(2) + (t43 * t67 + t70) * rSges(3,1) + (t43 * t66 - t71) * rSges(3,2)) + g(3) * (t43 * rSges(3,3) + t61) + (g(3) * (rSges(3,1) * t46 + rSges(3,2) * t50) + (g(1) * t47 - g(2) * t51) * (rSges(3,3) + pkin(8))) * t41) - m(4) * (g(1) * (t12 * rSges(4,1) + t11 * rSges(4,2) + (rSges(4,3) * t41 - t22) * t47 + t64) + g(2) * (t10 * rSges(4,1) + t9 * rSges(4,2) - rSges(4,3) * t68 + t60) + g(3) * (t20 * rSges(4,1) + t19 * rSges(4,2) + t43 * rSges(4,3) + t57)) - m(5) * (g(1) * (t12 * pkin(3) + (t12 * t49 + t63) * rSges(5,1) + (-t12 * t45 + t49 * t72) * rSges(5,2) - t76 * t11 + t58) + g(2) * (t10 * pkin(3) + (t10 * t49 - t62) * rSges(5,1) + (-t10 * t45 - t49 * t68) * rSges(5,2) - t76 * t9 + t60) + g(3) * (t20 * pkin(3) + (t20 * t49 + t74) * rSges(5,1) + (-t20 * t45 + t43 * t49) * rSges(5,2) - t76 * t19 + t57)) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) - t11 * rSges(6,3) + t55) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) - t9 * rSges(6,3) + t53) + g(3) * (t14 * rSges(6,1) - t13 * rSges(6,2) - t19 * rSges(6,3) + t54)) - m(7) * (g(1) * (t4 * pkin(5) + (-t11 * t44 + t4 * t48) * rSges(7,1) + (-t11 * t48 - t4 * t44) * rSges(7,2) + t75 * t3 + t55) + g(2) * (t2 * pkin(5) + (t2 * t48 - t9 * t44) * rSges(7,1) + (-t2 * t44 - t9 * t48) * rSges(7,2) + t75 * t1 + t53) + g(3) * (t14 * pkin(5) + (t14 * t48 - t19 * t44) * rSges(7,1) + (-t14 * t44 - t19 * t48) * rSges(7,2) + t75 * t13 + t54));
U  = t5;
