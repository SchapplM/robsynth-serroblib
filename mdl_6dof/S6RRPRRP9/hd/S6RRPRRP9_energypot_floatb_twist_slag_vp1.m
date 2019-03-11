% Calculate potential energy for
% S6RRPRRP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP9_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRP9_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP9_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP9_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP9_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:27:24
% EndTime: 2019-03-09 12:27:24
% DurationCPUTime: 0.40s
% Computational Cost: add. (364->131), mult. (545->158), div. (0->0), fcn. (633->12), ass. (0->56)
t37 = sin(pkin(11));
t71 = pkin(3) * t37;
t70 = rSges(6,3) + pkin(10);
t40 = cos(pkin(6));
t47 = cos(qJ(2));
t48 = cos(qJ(1));
t60 = t48 * t47;
t44 = sin(qJ(2));
t45 = sin(qJ(1));
t63 = t45 * t44;
t17 = -t40 * t60 + t63;
t43 = sin(qJ(5));
t69 = t17 * t43;
t61 = t48 * t44;
t62 = t45 * t47;
t19 = t40 * t62 + t61;
t68 = t19 * t43;
t38 = sin(pkin(6));
t67 = t38 * t44;
t66 = t38 * t45;
t65 = t38 * t47;
t64 = t38 * t48;
t59 = rSges(7,3) + qJ(6) + pkin(10);
t58 = qJ(3) + rSges(4,3);
t57 = pkin(7) + r_base(3);
t56 = t45 * pkin(1) + r_base(2);
t55 = t43 * t65;
t54 = t37 * t66;
t53 = t40 * pkin(8) + t57;
t52 = t48 * pkin(1) + pkin(8) * t66 + r_base(1);
t39 = cos(pkin(11));
t29 = t39 * pkin(3) + pkin(2);
t42 = -pkin(9) - qJ(3);
t51 = t29 * t67 + t40 * t71 + t42 * t65 + t53;
t20 = -t40 * t63 + t60;
t50 = pkin(3) * t54 - t19 * t42 + t20 * t29 + t52;
t18 = t40 * t61 + t62;
t49 = t18 * t29 + (-pkin(8) - t71) * t64 - t17 * t42 + t56;
t46 = cos(qJ(5));
t36 = pkin(11) + qJ(4);
t32 = cos(t36);
t31 = sin(t36);
t30 = t46 * pkin(5) + pkin(4);
t14 = t40 * t31 + t32 * t67;
t13 = t31 * t67 - t40 * t32;
t10 = t20 * t32 + t31 * t66;
t9 = t20 * t31 - t32 * t66;
t8 = t18 * t32 - t31 * t64;
t7 = t18 * t31 + t32 * t64;
t6 = t14 * t46 - t55;
t5 = -t14 * t43 - t46 * t65;
t4 = t10 * t46 + t68;
t3 = -t10 * t43 + t19 * t46;
t2 = t8 * t46 + t69;
t1 = t17 * t46 - t8 * t43;
t11 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t48 * rSges(2,1) - t45 * rSges(2,2) + r_base(1)) + g(2) * (t45 * rSges(2,1) + t48 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t57)) - m(3) * (g(1) * (t20 * rSges(3,1) - t19 * rSges(3,2) + t52) + g(2) * (t18 * rSges(3,1) - t17 * rSges(3,2) + t56) + g(3) * (t40 * rSges(3,3) + t53) + (g(1) * rSges(3,3) * t45 + g(3) * (rSges(3,1) * t44 + rSges(3,2) * t47) + g(2) * (-rSges(3,3) - pkin(8)) * t48) * t38) - m(4) * (g(1) * (t20 * pkin(2) + (t20 * t39 + t54) * rSges(4,1) + (-t20 * t37 + t39 * t66) * rSges(4,2) + t58 * t19 + t52) + g(2) * (t18 * pkin(2) - pkin(8) * t64 + (t18 * t39 - t37 * t64) * rSges(4,1) + (-t18 * t37 - t39 * t64) * rSges(4,2) + t58 * t17 + t56) + g(3) * ((t37 * rSges(4,1) + t39 * rSges(4,2)) * t40 + (-t58 * t47 + (t39 * rSges(4,1) - t37 * rSges(4,2) + pkin(2)) * t44) * t38 + t53)) - m(5) * (g(1) * (t10 * rSges(5,1) - t9 * rSges(5,2) + t19 * rSges(5,3) + t50) + g(2) * (t8 * rSges(5,1) - t7 * rSges(5,2) + t17 * rSges(5,3) + t49) + g(3) * (t14 * rSges(5,1) - t13 * rSges(5,2) - rSges(5,3) * t65 + t51)) - m(6) * (g(1) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t10 * pkin(4) + t70 * t9 + t50) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t8 * pkin(4) + t70 * t7 + t49) + g(3) * (t6 * rSges(6,1) + t5 * rSges(6,2) + t14 * pkin(4) + t70 * t13 + t51)) - m(7) * (g(1) * (t4 * rSges(7,1) + t3 * rSges(7,2) + pkin(5) * t68 + t10 * t30 + t59 * t9 + t50) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + pkin(5) * t69 + t8 * t30 + t59 * t7 + t49) + g(3) * (t6 * rSges(7,1) + t5 * rSges(7,2) - pkin(5) * t55 + t59 * t13 + t14 * t30 + t51));
U  = t11;
