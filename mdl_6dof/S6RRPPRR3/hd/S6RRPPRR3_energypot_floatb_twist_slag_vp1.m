% Calculate potential energy for
% S6RRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 09:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPPRR3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_energypot_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR3_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR3_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:54:24
% EndTime: 2019-03-09 08:54:24
% DurationCPUTime: 0.68s
% Computational Cost: add. (404->136), mult. (724->172), div. (0->0), fcn. (880->14), ass. (0->57)
t41 = sin(pkin(11));
t44 = cos(pkin(11));
t48 = sin(qJ(2));
t51 = cos(qJ(2));
t24 = -t48 * t41 + t51 * t44;
t78 = pkin(2) * t48;
t76 = pkin(10) + rSges(7,3);
t40 = sin(pkin(12));
t45 = cos(pkin(6));
t75 = t45 * t40;
t42 = sin(pkin(6));
t49 = sin(qJ(1));
t73 = t49 * t42;
t72 = t49 * t48;
t71 = t49 * t51;
t52 = cos(qJ(1));
t69 = t52 * t42;
t68 = t52 * t48;
t67 = t52 * t51;
t66 = qJ(4) + rSges(5,3);
t65 = pkin(7) + r_base(3);
t34 = t51 * pkin(2) + pkin(1);
t64 = t52 * t34 + r_base(1);
t63 = t40 * t73;
t62 = t40 * t69;
t61 = t45 * pkin(8) + t65;
t22 = t45 * t78 + (-pkin(8) - qJ(3)) * t42;
t60 = t52 * t22 + t49 * t34 + r_base(2);
t59 = t51 * t41 + t48 * t44;
t58 = -t49 * t22 + t64;
t57 = t45 * qJ(3) + t42 * t78 + t61;
t56 = t24 * t45;
t11 = -t49 * t56 - t52 * t59;
t21 = t59 * t45;
t12 = -t49 * t21 + t52 * t24;
t43 = cos(pkin(12));
t33 = t43 * pkin(4) + pkin(3);
t46 = -pkin(9) - qJ(4);
t55 = pkin(4) * t63 + t11 * t46 + t12 * t33 + t58;
t19 = t24 * t42;
t20 = t59 * t42;
t54 = pkin(4) * t75 + t19 * t46 + t20 * t33 + t57;
t10 = t52 * t21 + t49 * t24;
t9 = -t49 * t59 + t52 * t56;
t53 = -pkin(4) * t62 + t10 * t33 + t9 * t46 + t60;
t50 = cos(qJ(6));
t47 = sin(qJ(6));
t39 = pkin(12) + qJ(5);
t36 = cos(t39);
t35 = sin(t39);
t14 = t20 * t36 + t45 * t35;
t13 = t20 * t35 - t45 * t36;
t4 = t12 * t36 + t35 * t73;
t3 = t12 * t35 - t36 * t73;
t2 = t10 * t36 - t35 * t69;
t1 = t10 * t35 + t36 * t69;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t52 * rSges(2,1) - t49 * rSges(2,2) + r_base(1)) + g(2) * (t49 * rSges(2,1) + t52 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t65)) - m(3) * (g(1) * (t52 * pkin(1) + r_base(1) + (-t45 * t72 + t67) * rSges(3,1) + (-t45 * t71 - t68) * rSges(3,2)) + g(2) * (t49 * pkin(1) + r_base(2) + (t45 * t68 + t71) * rSges(3,1) + (t45 * t67 - t72) * rSges(3,2)) + g(3) * (t45 * rSges(3,3) + t61) + (g(3) * (rSges(3,1) * t48 + rSges(3,2) * t51) + (g(1) * t49 - g(2) * t52) * (rSges(3,3) + pkin(8))) * t42) - m(4) * (g(1) * (t12 * rSges(4,1) + t11 * rSges(4,2) + (rSges(4,3) * t42 - t22) * t49 + t64) + g(2) * (t10 * rSges(4,1) + t9 * rSges(4,2) - rSges(4,3) * t69 + t60) + g(3) * (t20 * rSges(4,1) + t19 * rSges(4,2) + t45 * rSges(4,3) + t57)) - m(5) * (g(1) * (t12 * pkin(3) + (t12 * t43 + t63) * rSges(5,1) + (-t12 * t40 + t43 * t73) * rSges(5,2) - t66 * t11 + t58) + g(2) * (t10 * pkin(3) + (t10 * t43 - t62) * rSges(5,1) + (-t10 * t40 - t43 * t69) * rSges(5,2) - t66 * t9 + t60) + g(3) * (t20 * pkin(3) + (t20 * t43 + t75) * rSges(5,1) + (-t20 * t40 + t45 * t43) * rSges(5,2) - t66 * t19 + t57)) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) - t11 * rSges(6,3) + t55) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) - t9 * rSges(6,3) + t53) + g(3) * (t14 * rSges(6,1) - t13 * rSges(6,2) - t19 * rSges(6,3) + t54)) - m(7) * (g(1) * (t4 * pkin(5) + (-t11 * t47 + t4 * t50) * rSges(7,1) + (-t11 * t50 - t4 * t47) * rSges(7,2) + t76 * t3 + t55) + g(2) * (t2 * pkin(5) + (t2 * t50 - t9 * t47) * rSges(7,1) + (-t2 * t47 - t9 * t50) * rSges(7,2) + t76 * t1 + t53) + g(3) * (t14 * pkin(5) + (t14 * t50 - t19 * t47) * rSges(7,1) + (-t14 * t47 - t19 * t50) * rSges(7,2) + t76 * t13 + t54));
U  = t5;
